// =============================================================================
// CD-HI/CD-HIT
//
// Cluster Database at High Identity Threshold
//
// CD-HIT clusters protein sequence database at high sequence identity threshold.
// This program can remove the high sequence redundance efficiently.
//
// program written by
//                    Weizhong Li
//                    UCSD, San Diego Supercomputer Center
//                    La Jolla, CA, 92093
//                    Email liwz@sdsc.edu
//
//                 at
//                    Adam Godzik's lab
//                    The Burnham Institute
//                    La Jolla, CA, 92037
//                    Email adam@burnham-inst.org
//
// modified by:
//                    Limin Fu
//                    Center for Research in Biological Systems (CRBS), UCSD
//                    La Jolla, CA, 92093
//                    Email: l2fu@ucsd.edu, fu@daovm.net
// =============================================================================

// #ifdef _WIN32
// typedef int pid_t;
// typedef _off_t off_t;
// #endif

import std;

//#include <stdio.h>
//#include <ctype.h>
//#include <stdint.h>
//#include <time.h>

#ifndef NO_OPENMP

import <omp.h>;

inline constexpr auto WITH_OPENMP = "(+OpenMP)";

#else

constexpr string WITH_OPENMP = "";
#define omp_set_num_threads(T) (T = T)
#define omp_get_thread_num() 0

#endif

#ifndef NO_ZLIB
// #define Z_LARGE64 1
// #define Z_WANT64 1
import <zlib.h>;
#endif

import ScoreMatrix;

#define MAX_NA 6
#define MAX_UAA 21
#define MAX_DIAG (MAX_SEQ << 1) // MAX_DIAG be twice of MAX_SEQ
#define MAX_GAP MAX_SEQ			// MAX_GAP <= MAX_SEQ
#define MAX_DES 300000
#define MAX_LINE_SIZE 300000
#define MAX_FILE_NAME 1280
#define MAX_SEG 50
#define MAX_BIN_SWAP 2E9
#define MAX_TABLE_SIZE 50000000
#define CLOCK_TICKS 100
#define FAILED_FUNC 1
#define OK_FUNC 0

#define IS_REP 1
#define IS_REDUNDANT 2
#define IS_PROCESSED 16
#define IS_MINUS_STRAND 32

// if the longset sequence is longer than 65535, I use INT4
#if defined(LONG_SEQ) && defined(SHORT_SEQ)
#error "LONG_SEQ and SHORT_SEQ cannot be defined at the same time"
#endif

#ifdef LONG_SEQ
#ifdef SHORT_SEQ
#error "LONG_SEQ and SHORT_SEQ cannot be defined at the same time"
#endif
typedef short UINT4;
#else
#ifdef SHORT_SEQ
typedef byte INTs;
#else
typedef unsigned short INTs;
#endif
#endif

using namespace std;

// class function definition
constexpr char aa[] = { "ARNDCQEGHILKMFPSTWYVBZX" };
//{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,2,6,20};
extern int aa2idx[];

constexpr auto MEGA_MiBi = 1'000'000;

// for primitive types only
template <class TYPE>
class NVector
{
public:
	TYPE* items;
	int size;
	int capacity;

	NVector()
	{
		size = capacity = 0;
		items = NULL;
	}
	NVector(int n, const TYPE& v = TYPE())
	{
		size = capacity = 0;
		items = NULL;
		Resize(n, v);
	}
	NVector(const NVector& other)
	{
		size = capacity = 0;
		items = NULL;
		if (other.items)
		{
			Resize(other.size);
			memcpy(items, other.items, other.size * sizeof(TYPE));
		}
	}

	~NVector()
	{
		if (items)
			free(items);
	}

	int Size() const { return size; }
	void Clear()
	{
		if (items)
			free(items);
		size = capacity = 0;
		items = NULL;
	}

	void Resize(int n, const TYPE& value = TYPE())
	{
		if (n == size && capacity > 0)
			return;
		int i;
		// When resize() is called, probably this is the intended size,
		// and will not be changed frequently.
		if (n != capacity)
		{
			capacity = n;
			items = (TYPE*)realloc(items, capacity * sizeof(TYPE));
		}
		for (i = size; i < n; i++)
			items[i] = value;
		size = n;
	}
	void Append(const TYPE& item)
	{
		if (size + 1 >= capacity)
		{
			capacity = size + size / 5 + 1;
			items = (TYPE*)realloc(items, capacity * sizeof(TYPE));
		}
		items[size] = item;
		size++;
	}

	TYPE& operator[](const int i)
	{
		// if( i <0 or i >= size ) printf( "out of range\n" );
		return items[i];
	}
	TYPE& operator[](const int i) const
	{
		// if( i <0 or i >= size ) printf( "out of range\n" );
		return items[i];
	}
};
typedef NVector<int> VectorInt;
typedef std::vector<VectorInt> MatrixInt;

typedef NVector<int64_t> VectorInt64;
typedef std::vector<VectorInt64> MatrixInt64;

////////// Class definition //////////

typedef NVector<INTs> VectorIntX;
typedef std::vector<VectorIntX> MatrixIntX;

constexpr size_t NAA0 = 1;
constexpr size_t NAA1 = MAX_UAA;
constexpr size_t NAA2 = NAA1 * NAA1;
constexpr size_t NAA3 = NAA1 * NAA2;
constexpr size_t NAA4 = NAA2 * NAA2;
constexpr size_t NAA5 = NAA2 * NAA3;
constexpr size_t NAA6 = NAA3 * NAA3;
constexpr size_t NAA7 = NAA3 * NAA4;
constexpr size_t NAA8 = NAA4 * NAA4;
constexpr size_t NAA9 = NAA4 * NAA5;
constexpr size_t NAA10 = NAA5 * NAA5;
constexpr size_t NAA11 = NAA5 * NAA6;
constexpr size_t NAA12 = NAA6 * NAA6;
constexpr size_t NAAN_array[13]{
	NAA0,
	NAA1,
	NAA2,
	NAA3,
	NAA4,
	NAA5,
	NAA6,
	NAA7,
	NAA8,
	NAA9,
	NAA10,
	NAA11,
	NAA12,
};

extern int naa_stat_start_percent;
extern int naa_stat[5][61][4];

struct IndexCount
{
	int index;
	int count;

	IndexCount(int i = 0, int c = 0) { index = i, count = c; }
};

struct Sequence;

void bomb_error(const char* message);

struct Sequence
{
	// real sequence, if it is not stored swap file:
	char* data;
	// length of the sequence:
	size_t size;
	size_t bufsize;
	size_t size_R2; // size = size.R1 + size.R2 for back-to-back merged seq

	// uint32_t stats;

	// if swap != NULL, the sequence is stored in file.
	// swap is opened as temporary file, which will be deleted automatically
	// after the program is finished:
	FILE* swap;
	// stream offset of the sequence:
	size_t offset;

	// stream offset of the description string in the database:
	size_t des_begin, des_begin2;
	// total record length
	size_t tot_length;
	size_t tot_length2;

	char* identifier;

	// index of the sequence in the original database:
	size_t index;
	short state;
	size_t cluster_id;
	float identity;
	float distance;
	size_t coverage[4];

	Sequence();
	Sequence(const Sequence& other);
	Sequence(const Sequence& other, const Sequence& other2, size_t mode);
	~Sequence();

	void Clear();

	void operator=(const char* s);
	void operator+=(const char* s);

	void Resize(size_t n);
	void Reserve(size_t n);

	void Swap(Sequence& other);
	int Format();

	void ConvertBases();
	void trim(size_t trim_len);

	void SwapIn();
	void SwapOut();
	void PrintInfo(size_t id, FILE* fout, char* buf);
};

struct WorkingParam
{
	double aa1_cutoff = 0;
	double aas_cutoff = 0;
	double aan_cutoff = 0;
	int len_upper_bound = 0;
	int len_lower_bound = 0;

	WorkingParam(double a1 = 0, double a2 = 0, double an = 0)
	{
		aa1_cutoff = a1;
		aas_cutoff = a2;
		aan_cutoff = an;
	}

	int len_eff;
	int aln_cover_flag;
	unsigned int min_aln_lenS;
	unsigned int min_aln_lenL;
	int required_aa1;
	int required_aas; /* or aa2 */
	int required_aan;

	void ControlShortCoverage(int lenion);
	void ControlLongCoverage(int lenion);
	void ComputeRequiredBases(int NAA, int ssion);
};

// #define MAX_TABLE_SEQ (1<<22)
#define MAX_TABLE_SEQ 4000000

enum
{
	DP_BACK_NONE = 0,
	DP_BACK_LEFT_TOP = 1,
	DP_BACK_LEFT = 2,
	DP_BACK_TOP = 3
};

void bomb_error(const char* message);
void bomb_error(const char* message, const char* message2);
void bomb_warning(const char* message);
void bomb_warning(const char* message, const char* message2);
void format_seq(char* seq);

void strrev(const char* p);
int print_usage(const char* arg);
int print_usage_2d(const char* arg);
int print_usage_est(const char* arg);
int print_usage_div(const char* arg);
int print_usage_est_2d(const char* arg);
int print_usage_454(const char* arg);

void cal_aax_cutoff(double& aa1_cutoff, double& aa2_cutoff, double& aan_cutoff,	double cluster_thd, int tolerance, int naa_stat_start_percent,	int naa_stat[5][61][4], int NAA)
{
	aa1_cutoff = cluster_thd;
	aa2_cutoff = 1 - (1 - cluster_thd) * 2;
	aan_cutoff = 1 - (1 - cluster_thd) * NAA;
	if (tolerance == 0) return;

	int clstr_idx = (int)(cluster_thd * 100) - naa_stat_start_percent;
	if (clstr_idx < 0) clstr_idx = 0;
	double d2 = ((double)(naa_stat[tolerance - 1][clstr_idx][3])) / 100;
	double dn = ((double)(naa_stat[tolerance - 1][clstr_idx][5 - NAA])) / 100;
	aa2_cutoff = d2 > aa2_cutoff ? d2 : aa2_cutoff;
	aan_cutoff = dn > aan_cutoff ? dn : aan_cutoff;
	return;
}

void update_aax_cutoff(double& aa1_cutoff, double& aa2_cutoff, double& aan_cutoff, int tolerance, int naa_stat_start_percent, int naa_stat[5][61][4], int NAA, double cluster_thd)
{
	if (cluster_thd > 1.0) cluster_thd = 1.00;

	double aa1_t, aa2_t, aan_t;
	cal_aax_cutoff(aa1_t, aa2_t, aan_t, cluster_thd, tolerance, naa_stat_start_percent,
		naa_stat, NAA);
	if (aa1_t > aa1_cutoff) aa1_cutoff = aa1_t;
	if (aa2_t > aa2_cutoff) aa2_cutoff = aa2_t;
	if (aan_t > aan_cutoff) aan_cutoff = aan_t;
	return;
}

float current_time();

// some functions from very old cd-hit

void setaa_to_na();
void make_comp_iseq(int len, char* iseq_comp, char* iseq);


int calc_ann_list(int len, char* seqi, int NAA, int& aan_no, vector<int>& aan_list, vector<INTs>& aan_list_no, bool est)
{
	// check_aan_list
	aan_no = len - NAA + 1;
	for (auto j = 0; j < aan_no; j++) {
		aan_list[j] = 0;
		for (auto k = 0, k1 = NAA - 1; k < NAA; k++, k1--)
			aan_list[j] += seqi[j + k] * NAAN_array[k1];
	}
	if (est) {
		// for the short word containing 'N', mask it to '-1'
		for (auto j = 0; j < len; j++) {
			if (seqi[j] >= 4) {                      // here N is 4
				auto i0 = (j - NAA + 1 > 0) ? j - NAA + 1 : 0;
				auto i1 = j < aan_no ? j : aan_no - 1;
				for (auto i = i0; i <= i1; i++)
					aan_list[i] = -1;
			}
		}
	}

	std::sort(aan_list.begin(), aan_list.begin() + aan_no);
	for (auto j = 0; j < aan_no; j++)
		aan_list_no[j] = 1;
	for (auto j = aan_no - 1; j; j--) {
		if (aan_list[j] == aan_list[j - 1]) {
			aan_list_no[j - 1] += aan_list_no[j];
			aan_list_no[j] = 0;
		}
	}
	return OK_FUNC;
}


void make_comp_short_word_index(int NAA, int* NAAN_array, vector<int>& Comp_AAN_idx) {
	int i, j, k, icomp, k1;
	int c[4] = { 3,2,1,0 };
	unsigned char short_word[32]; //short_word[12] is enough

	int NAA1 = NAAN_array[1];
	int NAAN = NAAN_array[NAA];

	for (i = 0; i < NAAN; i++) {
		// decompose i back to short_word
		for (k = i, j = 0; j < NAA; j++) {
			short_word[j] = (unsigned char)(k % NAA1);
			k = k / NAA1;
		}

		// calc_comp_aan_list
		icomp = 0;
		for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) icomp += c[short_word[k1]] * NAAN_array[k];

		Comp_AAN_idx[i] = icomp;
	}
} // make_comp_short_word_index

//quick_sort_idx calling (a, idx, 0, no-1)
//sort a with another array idx
//so that idx rearranged
int quick_sort_idx(int* a, int* idx, int lo0, int hi0) {
	int lo = lo0;
	int hi = hi0;
	int mid;
	int tmp;

	if (hi0 > lo0) {
		mid = a[(lo0 + hi0) / 2];

		while (lo <= hi) {
			while ((lo < hi0) && (a[lo] < mid)) lo++;
			while ((hi > lo0) && (a[hi] > mid)) hi--;
			if (lo <= hi) {
				tmp = a[lo];   a[lo] = a[hi];     a[hi] = tmp;
				tmp = idx[lo]; idx[lo] = idx[hi]; idx[hi] = tmp;
				lo++; hi--;
			}
		} // while

		if (lo0 < hi) quick_sort_idx(a, idx, lo0, hi);
		if (lo < hi0) quick_sort_idx(a, idx, lo, hi0);
	} // if ( hi0 > lo0)
	return 0;
} // quick_sort_idx


//decreasing can not use reverse of quick_sort_idx due to tie
//quick_sort_idxr calling (a, idx, 0, no-1)
//sort a with another array idx
//so that idx rearranged
int quick_sort_idxr(int* a, int* idx, int lo0, int hi0) {
	int lo = lo0;
	int hi = hi0;
	int mid;
	int tmp;

	if (hi0 > lo0) {
		mid = a[(lo0 + hi0) / 2];

		while (lo <= hi) {
			while ((lo < hi0) && (a[lo] > mid)) lo++;
			while ((hi > lo0) && (a[hi] < mid)) hi--;
			if (lo <= hi) {
				tmp = a[lo];   a[lo] = a[hi];     a[hi] = tmp;
				tmp = idx[lo]; idx[lo] = idx[hi]; idx[hi] = tmp;
				lo++; hi--;
			}
		} // while

		if (lo0 < hi) quick_sort_idxr(a, idx, lo0, hi);
		if (lo < hi0) quick_sort_idxr(a, idx, lo, hi0);
	} // if ( hi0 > lo0)
	return 0;
} // quick_sort_idxr
