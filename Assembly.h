#ifndef __ASSEMBLY__
#define __ASSEMBLY__
#include "CommandLines.h"
#include "Overlaps.h"
#include "Process_Read.h"
#include "Hash_Table.h"
#include "Correct.h"

#define FORWARD 0
#define REVERSE_COMPLEMENT (0x8000000000000000)

#define Get_Cigar_Type(RECORD) (RECORD&3)
#define Get_Cigar_Length(RECORD) (RECORD>>2)

#define RESEED_DP 4
#define RESEED_PEAK_RATE 0.15
#define RESEED_LEN 2000
#define RESEED_HP_RATE 0.9

typedef struct {
	int is_final, save_ov;
	// chaining and overlapping related buffers
	UC_Read self_read, ovlp_read;
	Candidates_list clist;
	overlap_region_alloc olist;
    overlap_region_alloc olist_hp;
	ha_abuf_t *ab;
    ha_abufl_t *abl;
	// error correction related buffers
	int64_t num_read_base, num_correct_base, num_recorrect_base;
	Cigar_record cigar1;
	Graph POA_Graph;
	Graph DAGCon;
	Correct_dumy correct;
	haplotype_evdience_alloc hap;
	Round2_alignment round2;
    kvec_t_u32_warp b_buf;
    kvec_t_u64_warp r_buf;
    kvec_t_u8_warp k_flag;
    overlap_region tmp_region;
    ma_utg_v *ua;
    st_mt_t sp;
	bit_extz_t exz;
} ha_ovec_buf_t;

int ha_assemble(void);
int ha_assemble_pair(void);
void ug_idx_build(ma_ug_t *ug, int hap_n);
ha_ovec_buf_t *ha_ovec_init(int is_final, int save_ov, int is_ug);
ha_ovec_buf_t *ha_ovec_buf_init(void *km, int is_final, int save_ov, int is_ug);
void ha_ovec_destroy(ha_ovec_buf_t *b);
int64_t ha_ovec_mem(const ha_ovec_buf_t *b, int64_t *mem_a);
int ha_assemble_ovec(void);

#endif
