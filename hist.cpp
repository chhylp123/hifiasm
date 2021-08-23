#include <stdio.h>
#include <assert.h>
#include "htab.h"

static void ha_hist_line(int c, int x, int exceed, int64_t cnt)
{
	int j;
	if (c >= 0) fprintf(stderr, "[M::%s] %5d: ", __func__, c);
	else fprintf(stderr, "[M::%s] %5s: ", __func__, "rest");
	for (j = 0; j < x; ++j) fputc('*', stderr);
	if (exceed) fputc('>', stderr);
	fprintf(stderr, " %lld\n", (long long)cnt);
}

void print_hist_lines(int n_cnt, int start_cnt, const int64_t *cnt)
{
	const int hist_max = 100;
	int i, start, low_i, max_i, max;
	// determine the start point
	assert(n_cnt > start_cnt);
	start = cnt[1] > 0? 1 : 2;

	// find the low point from the left
	low_i = start > start_cnt? start : start_cnt;
	for (i = low_i; i < n_cnt; ++i)
		if (cnt[i] > cnt[i-1]) break;
	low_i = i - 1;
	fprintf(stderr, "[M::%s] lowest: count[%d] = %ld\n", __func__, low_i, (long)cnt[low_i]);

	// find the highest peak
	max_i = start > start_cnt? start : start_cnt, max = cnt[max_i];
	for (i = max_i; i < n_cnt; ++i)
		if (cnt[i] > max)
			max = cnt[i], max_i = i;
	fprintf(stderr, "[M::%s] highest: count[%d] = %ld\n", __func__, max_i, (long)cnt[max_i]);

	for (i = start; i < n_cnt; ++i) {
		int x, exceed = 0;
		x = (int)((double)hist_max * cnt[i] / cnt[max_i] + .499);
		if (x > hist_max) exceed = 1, x = hist_max; // may happen if cnt[2] is higher
		if (i > max_i && x == 0) break;
		ha_hist_line(i, x, exceed, cnt[i]);
	}
}

int adj_m_peak_hom(int m_peak_hom, int max_i, int max2_i, int max3_i, int *peak_het)
{
	int64_t mm[3], d, min_i, min_d, i;
	mm[0] = max2_i; mm[1] = max_i; mm[2] = max3_i;
	for (i = 0, min_i = -1, min_d = -1; i < 3; i++){
		if(mm[i] <= 0) continue;
		d = (mm[i] >= m_peak_hom?mm[i]-m_peak_hom:m_peak_hom-mm[i]);
		if(min_d == -1 || min_d > d || (min_d == d && i == 1)){
			min_d = d; min_i = i;
		}
	}
	if(min_i < 0) return m_peak_hom;
	if(mm[min_i] < m_peak_hom){
		d = m_peak_hom - mm[min_i];
		if(d >= mm[min_i]*0.51) {
			*peak_het = mm[min_i];
			return m_peak_hom;
		}
	}

	for (i = min_i-1; i >= 0; i--){
		if(mm[i] <= 0) continue;
		*peak_het = mm[i];
		break;
	}
	return mm[min_i];
}

int ha_analyze_count(int n_cnt, int start_cnt, int m_peak_hom, const int64_t *cnt, int *peak_het)
{
	const int hist_max = 100;
	int i, start, low_i, max_i, max2_i, max3_i;
	int64_t max, max2, max3, min;

	// determine the start point
	assert(n_cnt > start_cnt);
	*peak_het = -1;
	start = cnt[1] > 0? 1 : 2;

	// find the low point from the left
	low_i = start > start_cnt? start : start_cnt;
	for (i = low_i + 1; i < n_cnt; ++i)
		if (cnt[i] > cnt[i-1]) break;
	low_i = i - 1;
	fprintf(stderr, "[M::%s] lowest: count[%d] = %ld\n", __func__, low_i, (long)cnt[low_i]);
	if (low_i == n_cnt - 1) return -1; // low coverage

	// find the highest peak
	max_i = low_i + 1, max = cnt[max_i];
	for (i = low_i + 1; i < n_cnt; ++i)
		if (cnt[i] > max)
			max = cnt[i], max_i = i;
	fprintf(stderr, "[M::%s] highest: count[%d] = %ld\n", __func__, max_i, (long)cnt[max_i]);

	// print histogram
	for (i = start; i < n_cnt; ++i) {
		int x, exceed = 0;
		x = (int)((double)hist_max * cnt[i] / cnt[max_i] + .499);
		if (x > hist_max) exceed = 1, x = hist_max; // may happen if cnt[2] is higher
		if (i > max_i && x == 0) break;
		ha_hist_line(i, x, exceed, cnt[i]);
	}
	{
		int x, exceed = 0;
		int64_t rest = 0;
		for (; i < n_cnt; ++i) rest += cnt[i];
		x = (int)((double)hist_max * rest / cnt[max_i] + .499);
		if (x > hist_max) exceed = 1, x = hist_max;
		ha_hist_line(-1, x, exceed, rest);
	}

	// look for smaller peak on the low end
	max2 = -1; max2_i = -1;
	for (i = max_i - 1; i > low_i; --i) {
		if (cnt[i] >= cnt[i-1] && cnt[i] >= cnt[i+1]) {
			if (cnt[i] > max2) max2 = cnt[i], max2_i = i;
		}
	}
	if (max2_i > low_i && max2_i < max_i) {
		for (i = max2_i + 1, min = max; i < max_i; ++i)
			if (cnt[i] < min) min = cnt[i];
		if (max2 < max * 0.05 || min > max2 * 0.95)
			max2 = -1, max2_i = -1;
	}
	if (max2 > 0) fprintf(stderr, "[M::%s] left: count[%d] = %ld\n", __func__, max2_i, (long)cnt[max2_i]);
	else fprintf(stderr, "[M::%s] left: none\n", __func__);

	// look for smaller peak on the high end
	max3 = -1; max3_i = -1;
	for (i = max_i + 1; i < n_cnt - 1; ++i) {
		if (cnt[i] >= cnt[i-1] && cnt[i] >= cnt[i+1]) {
			if (cnt[i] > max3) max3 = cnt[i], max3_i = i;
		}
	}
	if (max3_i > max_i) {
		for (i = max_i + 1, min = max; i < max3_i; ++i)
			if (cnt[i] < min) min = cnt[i];
		if (max3 < max * 0.05 || min > max3 * 0.95 || max3_i > max_i * 2.5)
			max3 = -1, max3_i = -1;
	}
	if (max3 > 0) fprintf(stderr, "[M::%s] right: count[%d] = %ld\n", __func__, max3_i, (long)cnt[max3_i]);
	else fprintf(stderr, "[M::%s] right: none\n", __func__);

	if(m_peak_hom > 0) return adj_m_peak_hom(m_peak_hom, max_i, max2_i, max3_i, peak_het);
	if (max3_i > 0) {
		*peak_het = max_i;
		return max3_i;
	} else {
		if (max2_i > 0) *peak_het = max2_i;
		return max_i;
	}
}
