#define __STDC_LIMIT_MACROS
#include <pthread.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include "kthread.h"

#if (defined(WIN32) || defined(_WIN32)) && defined(_MSC_VER)
#define __sync_fetch_and_add(ptr, addend)     _InterlockedExchangeAdd((void*)ptr, addend)
#endif

/************
 * kt_for() *
 ************/

struct kt_for_t;

typedef struct {
	struct kt_for_t *t;
	long i;
} ktf_worker_t;

typedef struct kt_for_t {
	int n_threads;
	long n;
	ktf_worker_t *w;
	void (*func)(void*,long,int);
	void *data;
} kt_for_t;

static inline long steal_work(kt_for_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	return k >= t->n? -1 : k;
}

static void ktf_run(ktf_worker_t *w)
{
	long i;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		if (i >= w->t->n) break;
		w->t->func(w->t->data, i, w - w->t->w);
	}
	while ((i = steal_work(w->t)) >= 0)
		w->t->func(w->t->data, i, w - w->t->w);
}

static void *ktf_worker(void *data)
{
	ktf_run((ktf_worker_t*)data);
	pthread_exit(0);
}

/* A pool of n_threads threads that stay alive between kt_for() calls: the creation and the join of
 * the threads cost about 2.8 ms for each call with 64 threads. Each job runs exactly as with new
 * threads (worker i is thread i). A pool runs one job at a time: a call that finds no free pool
 * (concurrent or nested calls) makes a new pool, up to KTF_MAX_POOLS, else it creates threads. */
#define KTF_MAX_POOLS 8

#ifdef HA_KTF_STAT
#include <stdio.h>
#include <sys/time.h>
/* the check build HA_KTF_STAT: for each caller of kt_for() (the return address, as an offset in the
 * executable for addr2line), the calls, the wall time, the tail (from the first worker without items
 * to the end) and the idle thread-seconds (the sum over the workers of the time from their last item
 * to the end), printed at exit */
extern char __executable_start;
static double ktf_now(void) { struct timeval tv; gettimeofday(&tv, 0); return tv.tv_sec + tv.tv_usec * 1e-6; }
typedef struct { void *a; long n; double wall, tail, idle; } ktf_stat_t;
static ktf_stat_t ktf_st[512];
static pthread_mutex_t ktf_st_mutex = PTHREAD_MUTEX_INITIALIZER;
static int ktf_st_cmp(const void *x, const void *y)
{
	double a = ((const ktf_stat_t*)x)->idle, b = ((const ktf_stat_t*)y)->idle;
	return a < b? 1 : a > b? -1 : 0;
}
static void ktf_st_report(void)
{
	int i;
	qsort(ktf_st, 512, sizeof(ktf_stat_t), ktf_st_cmp);
	for (i = 0; i < 40 && ktf_st[i].a; ++i)
		fprintf(stderr, "[M::ha_ktf_stat] 0x%lx calls %ld wall %.3f tail %.3f idle %.3f\n", (unsigned long)((char*)ktf_st[i].a - &__executable_start),
				ktf_st[i].n, ktf_st[i].wall, ktf_st[i].tail, ktf_st[i].idle);
}
static int ktf_st_set = atexit(ktf_st_report);
static void ktf_st_add(void *a, double wall, double tail, double idle)
{
	int i;
	pthread_mutex_lock(&ktf_st_mutex);
	for (i = (int)(((unsigned long)a>>4)&511); ktf_st[i].a && ktf_st[i].a != a; i = (i+1)&511);
	ktf_st[i].a = a, ++ktf_st[i].n, ktf_st[i].wall += wall, ktf_st[i].tail += tail, ktf_st[i].idle += idle;
	pthread_mutex_unlock(&ktf_st_mutex);
}
#endif

typedef struct ktf_pool_t {
	int n_threads, busy;
#ifdef HA_KTF_STAT
	double *td;
#endif
	kt_for_t t;
	pthread_mutex_t mutex;
	pthread_cond_t cv_go, cv_done;
	long gen;
	int n_done;
	struct ktf_pool_t *next;
} ktf_pool_t;

typedef struct {
	ktf_pool_t *p;
	int i;
} ktf_pool_arg_t;

static pthread_mutex_t ktf_pools_mutex = PTHREAD_MUTEX_INITIALIZER;
static ktf_pool_t *ktf_pools = 0;
static int ktf_n_pools = 0;

static void *ktf_pool_worker(void *data)
{
	ktf_pool_arg_t *a = (ktf_pool_arg_t*)data;
	ktf_pool_t *p = a->p;
	int i = a->i;
	long gen = 0;
	free(a);
	for (;;) {
		pthread_mutex_lock(&p->mutex);
		while (p->gen == gen) pthread_cond_wait(&p->cv_go, &p->mutex);
		gen = p->gen;
		pthread_mutex_unlock(&p->mutex);
		ktf_run(&p->t.w[i]);
#ifdef HA_KTF_STAT
		p->td[i] = ktf_now();
#endif
		pthread_mutex_lock(&p->mutex);
		if (++p->n_done == p->n_threads) pthread_cond_signal(&p->cv_done);
		pthread_mutex_unlock(&p->mutex);
	}
	return 0;
}

static ktf_pool_t *ktf_pool_get(int n_threads)
{
	ktf_pool_t *p;
	int i;
	pthread_mutex_lock(&ktf_pools_mutex);
	for (p = ktf_pools; p; p = p->next)
		if (!p->busy && p->n_threads == n_threads) break;
	if (p == 0 && ktf_n_pools < KTF_MAX_POOLS) {
		p = (ktf_pool_t*)calloc(1, sizeof(ktf_pool_t));
		p->n_threads = n_threads;
		p->t.n_threads = n_threads;
		p->t.w = (ktf_worker_t*)calloc(n_threads, sizeof(ktf_worker_t));
#ifdef HA_KTF_STAT
		p->td = (double*)calloc(n_threads, sizeof(double));
#endif
		pthread_mutex_init(&p->mutex, 0);
		pthread_cond_init(&p->cv_go, 0);
		pthread_cond_init(&p->cv_done, 0);
		for (i = 0; i < n_threads; ++i) {
			pthread_t tid;
			ktf_pool_arg_t *a = (ktf_pool_arg_t*)malloc(sizeof(ktf_pool_arg_t));
			a->p = p, a->i = i;
			p->t.w[i].t = &p->t;
			if (pthread_create(&tid, 0, ktf_pool_worker, a) != 0) {
				free(a);
				break;
			}
			pthread_detach(tid);
		}
		p->next = ktf_pools, ktf_pools = p;
		++ktf_n_pools;
		if (i < n_threads) { /* a thread was not created: the pool stays busy and is never used; its workers wait */
			p->busy = 1;
			p = 0;
		}
	}
	if (p) p->busy = 1;
	pthread_mutex_unlock(&ktf_pools_mutex);
	return p;
}

static void ktf_pool_put(ktf_pool_t *p)
{
	pthread_mutex_lock(&ktf_pools_mutex);
	p->busy = 0;
	pthread_mutex_unlock(&ktf_pools_mutex);
}

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n)
{
	if (n_threads > 1) {
		int i;
		kt_for_t t;
		ktf_pool_t *p = ktf_pool_get(n_threads);
		if (p) {
#ifdef HA_KTF_STAT
			double t0 = ktf_now(), t1, mn, idle = 0;
#endif
			p->t.func = func, p->t.data = data, p->t.n = n;
			for (i = 0; i < n_threads; ++i) p->t.w[i].i = i;
			pthread_mutex_lock(&p->mutex);
			p->n_done = 0, ++p->gen;
			pthread_cond_broadcast(&p->cv_go);
			while (p->n_done < n_threads) pthread_cond_wait(&p->cv_done, &p->mutex);
			pthread_mutex_unlock(&p->mutex);
#ifdef HA_KTF_STAT
			t1 = ktf_now(); mn = t1;
			for (i = 0; i < n_threads; ++i) { if (p->td[i] < mn) mn = p->td[i]; idle += t1 - p->td[i]; }
			ktf_st_add((void*)func, t1 - t0, t1 - mn, idle);
#endif
			ktf_pool_put(p);
			return;
		}
		pthread_t *tid;
		char *ok;
		t.func = func, t.data = data, t.n_threads = n_threads, t.n = n;
		t.w = (ktf_worker_t*)calloc(n_threads, sizeof(ktf_worker_t));
		tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
		for (i = 0; i < n_threads; ++i)
			t.w[i].t = &t, t.w[i].i = i;
		ok = (char*)calloc(n_threads, 1);
		for (i = 0; i < n_threads; ++i) ok[i] = (pthread_create(&tid[i], 0, ktf_worker, &t.w[i]) == 0);
		for (i = 0; i < n_threads; ++i) /* the caller runs the workers that were not created */
			if (!ok[i]) ktf_run(&t.w[i]);
		for (i = 0; i < n_threads; ++i)
			if (ok[i]) pthread_join(tid[i], 0);
		free(ok); free(tid); free(t.w);
	} else {
		long j;
		for (j = 0; j < n; ++j) func(data, j, 0);
	}
}

/*****************
 * kt_pipeline() *
 *****************/

struct ktp_t;

typedef struct {
	struct ktp_t *pl;
	int64_t index;
	int step;
	void *data;
} ktp_worker_t;

typedef struct ktp_t {
	void *shared;
	void *(*func)(void*, int, void*);
	int64_t index;
	int n_workers, n_steps;
	ktp_worker_t *workers;
	pthread_mutex_t mutex;
	pthread_cond_t cv;
} ktp_t;

static void *ktp_worker(void *data)
{
	ktp_worker_t *w = (ktp_worker_t*)data;
	ktp_t *p = w->pl;
	while (w->step < p->n_steps) {
		// test whether we can kick off the job with this worker
		pthread_mutex_lock(&p->mutex);
		for (;;) {
			int i;
			// test whether another worker is doing the same step
			for (i = 0; i < p->n_workers; ++i) {
				if (w == &p->workers[i]) continue; // ignore itself
				if (p->workers[i].step <= w->step && p->workers[i].index < w->index)
					break;
			}
			if (i == p->n_workers) break; // no workers with smaller indices are doing w->step or the previous steps
			pthread_cond_wait(&p->cv, &p->mutex);
		}
		pthread_mutex_unlock(&p->mutex);

		// working on w->step
		w->data = p->func(p->shared, w->step, w->step? w->data : 0); // for the first step, input is NULL

		// update step and let other workers know
		pthread_mutex_lock(&p->mutex);
		w->step = w->step == p->n_steps - 1 || w->data? (w->step + 1) % p->n_steps : p->n_steps;
		if (w->step == 0) w->index = p->index++;
		pthread_cond_broadcast(&p->cv);
		pthread_mutex_unlock(&p->mutex);
	}
	pthread_exit(0);
}

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps)
{
	ktp_t aux;
	pthread_t *tid;
	int i;

	if (n_threads < 1) n_threads = 1;
	aux.n_workers = n_threads;
	aux.n_steps = n_steps;
	aux.func = func;
	aux.shared = shared_data;
	aux.index = 0;
	pthread_mutex_init(&aux.mutex, 0);
	pthread_cond_init(&aux.cv, 0);

	aux.workers = (ktp_worker_t*)calloc(n_threads, sizeof(ktp_worker_t));
	for (i = 0; i < n_threads; ++i) {
		ktp_worker_t *w = &aux.workers[i];
		w->step = 0; w->pl = &aux; w->data = 0;
		w->index = aux.index++;
	}

	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktp_worker, &aux.workers[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
	free(tid); free(aux.workers);

	pthread_mutex_destroy(&aux.mutex);
	pthread_cond_destroy(&aux.cv);
}
