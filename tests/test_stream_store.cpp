#include "../StreamStore.h"

#include <assert.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include <string>
#include <vector>

struct TestRec {
    std::string name;
    std::string seq;
    std::string qual;
    bool has_qual;
};

static std::string make_complex_seq(size_t len, uint64_t seed)
{
    static const char alphabet[] = "ACGTACGTACGTACGTNacgtnRYKMSWBDHVU";
    std::string s(len, 'A');
    uint64_t x = seed ? seed : 0x9e3779b97f4a7c15ULL;
    for (size_t i = 0; i < len; ++i) {
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        uint64_t r = x * 2685821657736338717ULL;
        if ((r & 63ULL) == 0) {
            s[i] = alphabet[r % (sizeof(alphabet) - 1)];
        } else {
            s[i] = "ACGT"[r & 3ULL];
        }
    }
    if (len > 0) s[0] = 'N';
    if (len > 1) s[len - 1] = 'n';
    if (len > 3) s[len / 2] = 'R';
    return s;
}

static std::string make_qual(size_t len, uint64_t seed)
{
    std::string q(len, 'I');
    uint64_t x = seed ^ 0xdeadbeefcafebabeULL;
    for (size_t i = 0; i < len; ++i) {
        x = x * 6364136223846793005ULL + 1442695040888963407ULL;
        q[i] = (char)(33 + ((x >> 32) % 60));
    }
    return q;
}

static int count_visible_files(const char *dir)
{
    DIR *d = opendir(dir);
    assert(d != NULL);
    int count = 0;
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) continue;
        count++;
    }
    closedir(d);
    return count;
}

static void run_roundtrip_case(double max_ram_gb, int n_threads, const char *tmp_subdir)
{
    std::vector<std::vector<TestRec>> batches(4);
    size_t lengths[] = {1, 2, 3, 4, 5, 7, 15, 16, 17, 63, 64, 65, 511, 4096, 120000, 650000, 600000};

    for (size_t b = 0; b < batches.size(); ++b) {
        size_t n_in_batch = (b == 2) ? 3 : (sizeof(lengths) / sizeof(lengths[0]));
        for (size_t i = 0; i < n_in_batch; ++i) {
            size_t l = lengths[(i + b * 3) % (sizeof(lengths) / sizeof(lengths[0]))];
            TestRec rec;
            rec.name = "read_b" + std::to_string(b) + "_i" + std::to_string(i) + "_len" + std::to_string(l);
            rec.seq = make_complex_seq(l, (b + 1) * 10007ULL + i * 313ULL);
            if (b == 1) {
                rec.has_qual = false;
            } else if (b == 2) {
                rec.has_qual = (i % 2 == 0);
            } else {
                rec.has_qual = true;
            }
            if (rec.has_qual) {
                rec.qual = make_qual(l, (b + 1) * 99991ULL + i * 733ULL);
            }
            batches[b].push_back(rec);
        }
    }

    setenv("TMPDIR", tmp_subdir, 1);
    ha_stream_store_t *st = ha_stream_store_create(2, max_ram_gb, NULL);
    assert(st != NULL);
    assert(!ha_stream_store_is_captured(st, 0));
    assert(!ha_stream_store_is_captured(st, 1));

    for (size_t b = 0; b < batches.size(); ++b) {
        int n_seq = (int)batches[b].size();
        std::vector<char *> seq_ptrs(n_seq);
        std::vector<int> seq_lens(n_seq);
        std::vector<ha_stream_raw_rec_t> raw_recs(n_seq);
        for (int i = 0; i < n_seq; ++i) {
            seq_ptrs[i] = const_cast<char *>(batches[b][i].seq.data());
            seq_lens[i] = (int)batches[b][i].seq.size();
            raw_recs[i].name = const_cast<char *>(batches[b][i].name.data());
            raw_recs[i].name_len = (uint32_t)batches[b][i].name.size();
            raw_recs[i].has_qual = batches[b][i].has_qual ? 1 : 0;
            raw_recs[i].qual_len = batches[b][i].has_qual ? (uint32_t)batches[b][i].qual.size() : 0;
            raw_recs[i].qual = batches[b][i].has_qual ? const_cast<char *>(batches[b][i].qual.data()) : NULL;
        }
        ha_stream_store_append_batch(st, 0, n_seq, seq_ptrs.data(), seq_lens.data(), raw_recs.data(), n_threads);
    }
    ha_stream_store_mark_captured(st, 0);
    ha_stream_store_mark_captured(st, 1);
    assert(ha_stream_store_is_captured(st, 0));

    if (max_ram_gb == 0.0) {
        assert(ha_stream_store_ram_bytes(st) == 0);
        assert(ha_stream_store_disk_bytes(st) > 0);
        assert(count_visible_files(tmp_subdir) == 0);
    } else if (max_ram_gb >= 1.0) {
        assert(ha_stream_store_disk_bytes(st) == 0);
        assert(ha_stream_store_ram_bytes(st) > 0);
    }

    /* Replay Pass 2 (non-consuming) and Pass 3 (consuming) */
    for (int pass = 2; pass <= 3; ++pass) {
        int consume = (pass == 3) ? 1 : 0;
        ha_stream_store_rewind_file(st, 0);
        for (size_t b = 0; b < batches.size(); ++b) {
            ha_stream_decoded_batch_t out;
            int ok = ha_stream_store_replay_batch(st, 0, consume, n_threads, &out);
            assert(ok == 1);
            assert(out.n_seq == batches[b].size());
            for (uint32_t i = 0; i < out.n_seq; ++i) {
                const TestRec &exp = batches[b][i];
                assert((size_t)out.seq_lens[i] == exp.seq.size());
                assert(memcmp(out.seqs[i], exp.seq.data(), exp.seq.size()) == 0);
                assert(out.name_lens[i] == exp.name.size());
                assert(memcmp(out.names_buf + out.name_offs[i], exp.name.data(), exp.name.size()) == 0);
                assert((out.has_qual[i] != 0) == exp.has_qual);
                if (exp.has_qual) {
                    assert(out.quals_buf != NULL);
                    assert(memcmp(out.quals_buf + out.qual_offs[i], exp.qual.data(), exp.qual.size()) == 0);
                }
                free(out.seqs[i]);
            }
            free(out.seqs);
            free(out.seq_lens);
            ha_stream_decoded_batch_free_aux(&out);
        }
        ha_stream_decoded_batch_t eof_out;
        assert(ha_stream_store_replay_batch(st, 0, consume, n_threads, &eof_out) == 0);
        assert(ha_stream_store_replay_batch(st, 1, consume, n_threads, &eof_out) == 0);
    }

    assert(ha_stream_store_ram_bytes(st) == 0);
    ha_stream_store_destroy(&st);
    assert(st == NULL);
    assert(count_visible_files(tmp_subdir) == 0);
}

int main(void)
{
    const char *base_tmp = getenv("TMPDIR");
    assert(base_tmp != NULL && base_tmp[0] != '\0');
    std::string sub = std::string(base_tmp) + "/ut_stream_store_" + std::to_string(getpid());
    mkdir(sub.c_str(), 0700);

    /* 1. Pure RAM mode at 1 thread and 8 threads */
    run_roundtrip_case(8.0, 1, sub.c_str());
    run_roundtrip_case(8.0, 8, sub.c_str());

    /* 2. Forced 100% O_TMPFILE spillover mode (max_ram_gb = 0) at 1 and 8 threads */
    run_roundtrip_case(0.0, 1, sub.c_str());
    run_roundtrip_case(0.0, 8, sub.c_str());

    /* 3. Hybrid RAM + O_TMPFILE spillover boundary (~512 KB cap) */
    run_roundtrip_case(0.0005, 4, sub.c_str());

    rmdir(sub.c_str());
    fprintf(stderr, "[PASS] test_stream_store passed all RAM and O_TMPFILE roundtrip invariants\n");
    return 0;
}
