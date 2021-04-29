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
	ret = ha_assemble();
    destory_opt(&asm_opt);
	fprintf(stderr, "[M::%s] Version: %s\n", __func__, HA_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, yak_realtime(), yak_cputime(), yak_peakrss_in_gb());
    return ret;
}
