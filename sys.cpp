#include <sys/resource.h>
#ifdef __GLIBC__
#include <malloc.h>
#endif
#include <sys/time.h>
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

#ifdef HA_MEMSTAT
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
static unsigned long long ha_memstat_val(const char *p, const char *key)
{
	unsigned long long v = 0;
	const char *q = strstr(p, key);
	if (q && (q = strstr(q, "size=\"")) != NULL) v = strtoull(q + 6, NULL, 10);
	return v;
}

///a diagnostic build (-DHA_MEMSTAT): the RSS now, and the bytes of the malloc arenas from the system, free,
///and in mmap blocks (malloc_info(): the totals after the last heap)
void ha_memstat(const char *tag)
{
	char *buf = NULL; size_t len = 0;
	const char *p, *q;
	long pages = 0, rss = 0;
	double g = 1073741824.0, sys, fr, mm;
	FILE *fp = open_memstream(&buf, &len);
	if (fp == NULL) return;
	malloc_info(0, fp); fclose(fp);
	for (p = buf; (q = strstr(p, "</heap>")) != NULL; p = q + 7);
	sys = ha_memstat_val(p, "<system type=\"current\"");
	fr = ha_memstat_val(p, "<total type=\"fast\"") + ha_memstat_val(p, "<total type=\"rest\"");
	mm = ha_memstat_val(p, "<total type=\"mmap\"");
	free(buf);
	fp = fopen("/proc/self/statm", "r");
	if (fp) {
		if (fscanf(fp, "%ld %ld", &pages, &rss) != 2) rss = 0;
		fclose(fp);
	}
	fprintf(stderr, "[M::ha_memstat::%.3f] %s: rss %.2f GB; malloc: arenas %.2f GB (free %.2f GB), mmap blocks %.2f GB\n",
			yak_realtime(), tag, rss * (double)sysconf(_SC_PAGESIZE) / g, sys / g, fr / g, mm / g);
}
#endif

#ifdef __GLIBC__
#include <pthread.h>
#include <stdio.h>
static pthread_t ha_trim_t;
static int ha_trim_on = 0;

static void *ha_trim_run(void *data)
{
#ifdef HA_MEMSTAT
	double t0 = yak_realtime();
	malloc_trim(0);
	fprintf(stderr, "[M::ha_memstat::%.3f] malloc_trim: %.3f s\n", yak_realtime(), yak_realtime() - t0);
#else
	malloc_trim(0);
#endif
	return NULL;
}

static void ha_trim_join(void)
{
	if (ha_trim_on) pthread_join(ha_trim_t, NULL);
	ha_trim_on = 0;
}
#endif

///give the free memory of the malloc arenas back to the system, in a thread (T2: 1.5-2.4 s; it locks one
///arena at a time). glibc keeps the freed small blocks in the arena of each thread, and a large mmap()
///of a later stage adds new pages on top of them.
void ha_mem_trim(void)
{
#ifdef __GLIBC__
	static int exit_set = 0;
	ha_trim_join();
	if (!exit_set) {
		atexit(ha_trim_join); exit_set = 1;
	}
	if (pthread_create(&ha_trim_t, NULL, ha_trim_run, NULL) == 0) ha_trim_on = 1;
	else malloc_trim(0);
#endif
}
