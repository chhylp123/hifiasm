#include <sys/resource.h>
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include "htab.h"

int yak_verbose = 3;

static double yak_realtime0;

double yak_cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

static inline double yak_realtime_core(void)
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

void yak_reset_realtime(void)
{
	yak_realtime0 = yak_realtime_core();
}

double yak_realtime(void)
{
	return yak_realtime_core() - yak_realtime0;
}

double yak_realtime_0(void)
{
	return yak_realtime_core();
}


long yak_peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

double yak_peakrss_in_gb(void)
{
	return yak_peakrss() / 1073741824.0;
}

double yak_cpu_usage(void)
{
	return (yak_cputime() + 1e-9) / (yak_realtime() + 1e-9);
}

double yak_available_mem_gb(void)
{
#ifdef __linux__
	FILE *fp;
	char buf[256];
	double avail = -1.0;
	fp = fopen("/proc/meminfo", "r");
	if (!fp) return -1.0;
	while (fgets(buf, sizeof(buf), fp)) {
		if (strncmp(buf, "MemAvailable:", 13) == 0) {
			avail = atof(buf + 13) / 1048576.0;
			break;
		}
	}
	fclose(fp);
	return avail;
#else
	return -1.0;
#endif
}
