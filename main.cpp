#include <stdio.h>
#include <stdlib.h>
#include "CommandLines.h"
#include "Process_Read.h"
#include "Assembly.h"
#include "Levenshtein_distance.h"
#include "yak.h"

int main(int argc, char *argv[])
{
	int i;
	void *flt_tab;

    init_opt(&asm_opt);

    if (!CommandLine_process(argc, argv, &asm_opt)) return 1;

	yak_reset_realtime();
	flt_tab = ha_gen_flt_tab(&asm_opt);
	ha_hf_destroy(flt_tab);
	if (0) {
	    Correct_Reads(asm_opt.number_of_round);
	}
    destory_opt(&asm_opt);
	fprintf(stderr, "[M::%s] Version: %s\n", __func__, "dummy");
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, yak_realtime(), yak_cputime(), yak_peakrss() / 1024.0 / 1024.0 / 1024.0);
    return 0;
}
