#!/usr/bin/env python3
"""Generate a tiny, deterministic, HiFi-like read set for the CI build smoke test.

Deterministic across platforms and architectures (Python's Mersenne Twister with
a fixed seed is platform independent), so the resulting assembly graph can be
compared byte for byte between x86_64 and aarch64.

Usage: python3 test/gen_hifi_reads.py [out.fq.gz]
"""
import sys
import gzip
import random


def main(out="reads.fq.gz"):
    random.seed(20240709)
    genome_len, read_len, coverage, err_rate = 60000, 15000, 40, 0.001
    bases = "ACGT"
    genome = ''.join(random.choice(bases) for _ in range(genome_len))
    nreads = (genome_len * coverage) // read_len
    recs = []
    for i in range(nreads):
        start = random.randint(0, genome_len - read_len)
        seq = list(genome[start:start + read_len])
        for j in range(len(seq)):
            if random.random() < err_rate:
                seq[j] = random.choice(bases)
        seq = ''.join(seq)
        recs.append("@read%d\n%s\n+\n%s\n" % (i, seq, 'I' * len(seq)))
    data = ''.join(recs).encode()
    opener = gzip.open if out.endswith('.gz') else open
    with opener(out, 'wb') as f:
        f.write(data)
    sys.stderr.write("wrote %d reads (%d bp genome, ~%dx) to %s\n"
                     % (nreads, genome_len, coverage, out))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else "reads.fq.gz")
