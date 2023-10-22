#include <stdio.h>
#include <stdlib.h>
#include "CommandLines.h"
#include "Process_Read.h"
#include "Assembly.h"
#include "Levenshtein_distance.h"
#include "htab.h"

int main(int argc, char *argv[])
{
	int i, ret;
	yak_reset_realtime();
    init_opt(&asm_opt);
    if (!CommandLine_process(argc, argv, &asm_opt)) return 0;
	
	// bit_extz_t exz, exz64; init_bit_extz_t(&exz, 2); init_bit_extz_t(&exz64, 2);
	
	// char *pstr = "GACCCAG", *tsrt = "GTTGTTAATTCCAT"; int32_t thre = 14; clear_align(exz); clear_align(exz64);
	// ed_band_cal_extension_64_0_w_trace((char*)pstr, strlen(pstr), (char*)tsrt, strlen(tsrt), thre, &exz);
	// // // ed_band_cal_semi_64_w_absent_diag((char*)pstr, strlen(pstr), (char*)tsrt, strlen(tsrt), thre, 0, &exz);
	// fprintf(stderr, "\n[M::%s::] exz.err::%d, exz.ps::%d, exz.pe::%d, exz.ts::%d, exz.te::%d\n", __func__, 
    //         exz.err, exz.ps, exz.pe, exz.ts, exz.te);
	// cigar_check((char*)pstr, (char*)tsrt, &(exz));
	// ed_band_cal_extension_64_0_w((char*)pstr, strlen(pstr), (char*)tsrt, strlen(tsrt), thre, &exz);
	// fprintf(stderr, "\n[M::%s::] exz.err::%d, exz.ps::%d, exz.pe::%d, exz.ts::%d, exz.te::%d\n", __func__, 
    //         exz.err, exz.ps, exz.pe, exz.ts, exz.te);
	// ed_band_cal_semi_infi_w((char*)pstr, strlen(pstr), (char*)tsrt, strlen(tsrt), thre, NULL, &exz);
	// ed_band_cal_semi_64_w((char*)pstr, strlen(pstr), (char*)tsrt, strlen(tsrt), thre, &exz64);
	// fprintf(stderr, "\n[M::%s::] exz.err::%d, exz64.err::%d, exz.ps::%d, exz64.ps::%d, exz.pe::%d, exz64.pe::%d, exz.ts::%d, exz64.ts::%d, exz.te::%d, exz64.te::%d\n", __func__, 
    //         exz.err, exz64.err, exz.ps, exz64.ps, exz.pe, exz64.pe, exz.ts, exz64.ts, exz.te, exz64.te);
	// ed_band_cal_semi_64_w_trace((char*)pstr, strlen(pstr), (char*)tsrt, strlen(tsrt), thre, &exz64);
	// cigar_check((char*)pstr, (char*)tsrt, &(exz64));
	
	
	// char *pstr = "TGT", *tsrt = "CTGT"; int32_t thre = 1;
	// ed_band_cal_global_infi_w((char*)pstr, strlen(pstr), (char*)tsrt, strlen(tsrt), thre, NULL, &exz);
	// ed_band_cal_global_64_w((char*)pstr, strlen(pstr), (char*)tsrt, strlen(tsrt), thre, &exz64);
	// ed_band_cal_global_64_w_trace((char*)pstr, strlen(pstr), (char*)tsrt, strlen(tsrt), thre, &exz64);
	// fprintf(stderr, "\n[M::%s::] exz.err::%d, exz64.err::%d, exz.ps::%d, exz64.ps::%d, exz.pe::%d, exz64.pe::%d, exz.ts::%d, exz64.ts::%d, exz.te::%d, exz64.te::%d\n", __func__, 
    //         exz.err, exz64.err, exz.ps, exz64.ps, exz.pe, exz64.pe, exz.ts, exz64.ts, exz.te, exz64.te);
	// cigar_check((char*)pstr, (char*)tsrt, &(exz64));
	
	// ed_band_cal_extension_infi0_w((char *)"AAT", 3, (char *)"ACTTTTTT", 8, 2, NULL, &exz);
    // ed_band_cal_extension_64_w((char *)"AAT", 3, (char *)"ACTTTTTT", 8, 2, &exz64);
	// fprintf(stderr, "\n[M::%s::] exz.err::%d, exz64.err::%d, exz.ps::%d, exz64.ps::%d, exz.pe::%d, exz64.pe::%d, exz.ts::%d, exz64.ts::%d, exz.te::%d, exz64.te::%d\n", __func__, 
    //         exz.err, exz64.err, exz.ps, exz64.ps, exz.pe, exz64.pe, exz.ts, exz64.ts, exz.te, exz64.te);

	//bit_extz_t exz; ///ed_band_cal_global_128bit(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
	// ed_band_cal_extension_128bit((char *)"AAGTTTA", 7, (char *)"CCTTTTTT", 8, 4, &exz);
	// ed_band_cal_extension_128bit((char *)"AA", 2, (char *)"ACTTTTTT", 8, 1, &exz);
	// fprintf(stderr, "ed_extension::%d, pe::%d, te::%d\n", exz.err, exz.pe, exz.te);
	// exit(1);
	
	

	// fprintf(stderr, "[M::%s::] ed_global::%d, ed_global_128bit::%d\n", __func__, 
	// ed_band_cal_global((char *)"ACT", 3, (char *)"AAT", 3, 1), 
	// ed_band_cal_global_128bit((char *)"ACT", 3, (char *)"AAT", 3, 1));

	// fprintf(stderr, "[M::%s::] ed_global::%d, ed_global_128bit::%d\n", __func__, 
	// ed_band_cal_global((char*)"ACTTTTTT", 8, (char*)"AATTTT", 6, 3),
	// ed_band_cal_global_128bit((char*)"ACTTTTTT", 8, (char*)"AATTTT", 6, 3));
	// exit(1);
	if(asm_opt.sec_in) ret = ha_assemble_pair();
	else if(asm_opt.dbg_ovec_cal) ret = ha_assemble_ovec();
	else ret = ha_assemble();
	
    destory_opt(&asm_opt);
	fprintf(stderr, "[M::%s] Version: %s\n", __func__, HA_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, yak_realtime(), yak_cputime(), yak_peakrss_in_gb());
    return ret;
}
