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
	
	ret = ha_assemble();
    destory_opt(&asm_opt);
	fprintf(stderr, "[M::%s] Version: %s\n", __func__, HA_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, yak_realtime(), yak_cputime(), yak_peakrss_in_gb());
    return ret;
}
