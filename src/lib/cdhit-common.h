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

#ifdef _WIN32
typedef int pid_t;
#include <sys/types.h>
typedef _off_t off_t;
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <string_view>
#include <vector>
#include <valarray>
#include <map>

#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <time.h>

#ifndef NO_OPENMP

#include <omp.h>
constexpr std::string_view WITH_OPENMP = "(+OpenMP)";

#else

constexpr string WITH_OPENMP = "";
#define omp_set_num_threads(T) (T = T)
#define omp_get_thread_num() 0

#endif

#ifdef WITH_ZLIB
// #define Z_LARGE64 1
// #define Z_WANT64 1
#include <zlib.h>
#endif

#include "ScoreMatrix.h"

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
constexpr char aa[] = {"ARNDCQEGHILKMFPSTWYVBZX"};
//{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,2,6,20};
extern int aa2idx[];

constexpr auto MEGA_MiBi = 1'000'000;

// the parent containter must guarantee continuous memory allocation.
// std::valarray could be used instead of std::vector.
template <class TYPE>
class Vector : public vector<TYPE>
{
public:
	Vector() : vector<TYPE>() {}
	Vector(size_t size) : vector<TYPE>(size) {}
	Vector(size_t size, const TYPE &deft) : vector<TYPE>(size, deft) {}

	void Append(const TYPE &item)
	{
		size_t n = this->size();
		if (n + 1 >= this->capacity())
			this->reserve(n + n / 5 + 1);
		this->push_back(item);
	}
};

// for primitive types only
template <class TYPE>
class NVector
{
public:
	TYPE *items;
	int size;
	int capacity;

	NVector()
	{
		size = capacity = 0;
		items = NULL;
	}
	NVector(int n, const TYPE &v = TYPE())
	{
		size = capacity = 0;
		items = NULL;
		Resize(n, v);
	}
	NVector(const NVector &other)
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

	void Resize(int n, const TYPE &value = TYPE())
	{
		if (n == size && capacity > 0)
			return;
		int i;
		// When resize() is called, probably this is the intended size,
		// and will not be changed frequently.
		if (n != capacity)
		{
			capacity = n;
			items = (TYPE *)realloc(items, capacity * sizeof(TYPE));
		}
		for (i = size; i < n; i++)
			items[i] = value;
		size = n;
	}
	void Append(const TYPE &item)
	{
		if (size + 1 >= capacity)
		{
			capacity = size + size / 5 + 1;
			items = (TYPE *)realloc(items, capacity * sizeof(TYPE));
		}
		items[size] = item;
		size++;
	}

	TYPE &operator[](const int i)
	{
		// if( i <0 or i >= size ) printf( "out of range\n" );
		return items[i];
	}
	TYPE &operator[](const int i) const
	{
		// if( i <0 or i >= size ) printf( "out of range\n" );
		return items[i];
	}
};
typedef NVector<int> VectorInt;
typedef Vector<VectorInt> MatrixInt;

typedef NVector<int64_t> VectorInt64;
typedef Vector<VectorInt64> MatrixInt64;

////////// Class definition //////////

typedef NVector<INTs> VectorIntX;
typedef Vector<VectorIntX> MatrixIntX;

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

class WordTable
{
private:
public:
	Vector<NVector<IndexCount>> indexCounts; // hold index and word counts of seqs
	Vector<Sequence *> sequences;
	int NAA;	// length of word
	int NAAN;	// rows of table
	char is_aa; // aa is for prot
	size_t size;
	int frag_count;

public:
	WordTable(int naa = 0, int naan = 0);
	void Init(int, int);
	void Clear();
	void SetDNA();
	int AddWordCounts(NVector<IndexCount> &counts, Sequence *seq, bool skipN = false);
	int AddWordCountsFrag(NVector<IndexCount> &counts, int frag, int frag_size, int repfrag);

	int AddWordCounts(int aan_no, Vector<int> &word_encodes,
					  Vector<INTs> &word_encodes_no, int idx, bool skipN = false);
	int AddWordCountsFrag(int aan_no, Vector<int> &word_encodes,
						  Vector<INTs> &word_encodes_no, int frag, int frag_size);
	int CountWords(int aan_no, Vector<int> &aan_list, Vector<INTs> &aan_list_no,
				   NVector<IndexCount> &lookCounts, NVector<uint32_t> &indexMapping,
				   bool est = false, int min = 0);
	void PrintAll();
}; // END class INDEX_TBL

void bomb_error(const char *message);

struct Sequence
{
	// real sequence, if it is not stored swap file:
	char *data;
	// length of the sequence:
	size_t size;
	size_t bufsize;
	size_t size_R2; // size = size.R1 + size.R2 for back-to-back merged seq

	// uint32_t stats;

	// if swap != NULL, the sequence is stored in file.
	// swap is opened as temporary file, which will be deleted automatically
	// after the program is finished:
	FILE *swap;
	// stream offset of the sequence:
	size_t offset;

	// stream offset of the description string in the database:
	size_t des_begin, des_begin2;
	// total record length
	size_t tot_length;
	size_t tot_length2;

	char *identifier;

	// index of the sequence in the original database:
	size_t index;
	short state;
	size_t cluster_id;
	float identity;
	float distance;
	size_t coverage[4];

	Sequence();
	Sequence(const Sequence &other);
	Sequence(const Sequence &other, const Sequence &other2, size_t mode);
	~Sequence();

	void Clear();

	void operator=(const char *s);
	void operator+=(const char *s);

	void Resize(size_t n);
	void Reserve(size_t n);

	void Swap(Sequence &other);
	int Format();

	void ConvertBases();
	void trim(size_t trim_len);

	void SwapIn();
	void SwapOut();
	void PrintInfo(size_t id, FILE *fout, char *buf);
};

struct WorkingParam
{
	double aa1_cutoff;
	double aas_cutoff; /* or aa2 */
	double aan_cutoff;
	int len_upper_bound;
	int len_lower_bound;

	WorkingParam(double a1 = 0, double a2 = 0, double an = 0)
	{
		Set(a1, a2, an);
	}
	void Set(double a1 = 0, double a2 = 0, double an = 0)
	{
		aa1_cutoff = a1;
		aas_cutoff = a2;
		aan_cutoff = an;
		len_upper_bound = 0;
		len_lower_bound = 0;
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

struct WorkingBuffer
{
	Vector<int> taap;
	Vector<int> word_encodes;
	Vector<int> word_encodes_backup;
	Vector<INTs> word_encodes_no;
	Vector<INTs> aap_list;
	Vector<INTs> aap_begin;
	// Vector<IndexCount>  indexCounts;
	NVector<IndexCount> lookCounts;
	NVector<uint32_t> indexMapping;
	MatrixInt64 score_mat;
	MatrixInt back_mat;
	Vector<int> diag_score;
	Vector<int> diag_score2;
	Vector<int> aan_list_comp;
	Vector<char> seqi_comp;
	int total_bytes;

	WorkingBuffer(size_t frag = 0, size_t maxlen = 0)
	{
		Set(frag, maxlen);
		seqi_comp.resize(MAX_SEQ);
	}
	void Set(size_t frag, size_t maxlen, const bool isEST, const size_t band_width)
	{
		bool est = isEST;
		size_t m = MAX_UAA * MAX_UAA;
		size_t max_len = maxlen;
		size_t band = max_len * max_len;
		if (est)
			m = m * m;
		if (band > band_width)
			band = band_width;
		taap.resize(m);
		aap_list.resize(max_len);
		aap_begin.resize(m);
		// indexCounts.resize( max_len );
		word_encodes.resize(max_len);
		word_encodes_no.resize(max_len);
		word_encodes_backup.resize(max_len);
		/* each table can not contain more than MAX_TABLE_SEQ representatives or fragments! */
		if (frag > MAX_TABLE_SEQ)
			frag = MAX_TABLE_SEQ;
		lookCounts.Resize(frag + 2);
		indexMapping.Resize(frag + 2);
		diag_score.resize(MAX_DIAG);
		diag_score2.resize(MAX_DIAG);
		aan_list_comp.resize(max_len);
		total_bytes = max_len;
		total_bytes += taap.size() * sizeof(int);
		total_bytes += word_encodes.size() * sizeof(int);
		total_bytes += word_encodes_backup.size() * sizeof(int);
		total_bytes += diag_score.size() * sizeof(int);
		total_bytes += diag_score2.size() * sizeof(int);
		total_bytes += aan_list_comp.size() * sizeof(int);
		total_bytes += word_encodes_no.size() * sizeof(INTs);
		total_bytes += aap_list.size() * sizeof(INTs);
		total_bytes += aap_begin.size() * sizeof(INTs);
		total_bytes += indexMapping.Size() * sizeof(uint32_t);
		// total_bytes += indexCounts.size()*sizeof(IndexCount);
		total_bytes += lookCounts.Size() * sizeof(IndexCount);
		total_bytes += max_len * (band * sizeof(int) + sizeof(VectorInt));
		total_bytes += max_len * (band * sizeof(int) + sizeof(VectorInt64));
	}

	int EncodeWords(Sequence *seq, int NA, bool est = false);
	void ComputeAAP(const char *seqi, int size);
	void ComputeAAP2(const char *seqi, int size);
};
extern Vector<size_t> Comp_AAN_idx;
extern ScoreMatrix mat;

class SequenceDB
{
public:
	size_t NAAN;
	Vector<Sequence *> sequences;
	Vector<int> rep_seqs;

	long long total_letter;
	long long total_desc;
	size_t max_len;
	size_t min_len;
	size_t len_n50;

	void Clear()
	{
		for (int i = 0; i < sequences.size(); i++)
			delete sequences[i];
		sequences.clear();
		rep_seqs.clear();
	}

	SequenceDB()
	{
		total_letter = 0;
		total_desc = 0;
		min_len = 0;
		max_len = 0;
		len_n50 = 0;
	}
	~SequenceDB() { Clear(); }

	void Read(const char *file);
	void Readgz(string file);

	void Read(const char *file, const char *file2);
	void Readgz(string file, string file2);

	void WriteClusters(const char *db, const char *newdb);
	void WriteClustersgz(const char *db, const char *newdb);

	void WriteClusters(const char *db, const char *db_pe, const char *newdb, const char *newdb_pe);
	void WriteClustersgz(const char *db, const char *db_pe, const char *newdb, const char *newdb_pe);

	void WriteExtra1D();
	void WriteExtra2D(SequenceDB &other);
	void DivideSave(const char *db, const char *newdb, int n);

	void SwapIn(int seg, bool reponly = false);
	void SwapOut(int seg);

	// void Sort( int first, int last );
	void SortDivide(bool sort = true);
	void MakeWordTable();

	size_t MinimalMemory(int frag_no, int bsize, int T, size_t extra = 0);

	void ClusterOne(Sequence *seq, int id, WordTable &table,
					WorkingParam &param, WorkingBuffer &buf);

	// void SelfComparing( int start, int end, WordTable & table,
	//     WorkingParam & param, WorkingBuffer & buf );

	void ComputeDistance();
	void DoClustering();
	void DoClustering(int T);
	void ClusterTo(SequenceDB &other);
	int CheckOne(Sequence *seq, WordTable &tab, WorkingParam &par, WorkingBuffer &buf);
	int CheckOneEST(Sequence *seq, WordTable &tab, WorkingParam &par, WorkingBuffer &buf);
	int CheckOneAA(Sequence *seq, WordTable &tab, WorkingParam &par, WorkingBuffer &buf);
};

void bomb_error(const char *message);
void bomb_error(const char *message, const char *message2);
void bomb_warning(const char *message);
void bomb_warning(const char *message, const char *message2);
void format_seq(char *seq);
int diag_test_aapn(int NAA1, char iseq2[], int len1, int len2,
				   WorkingBuffer &buffer, int &best_sum,
				   int band_width, int &band_left, int &band_center, int &band_right, int required_aa1);
int diag_test_aapn_est(int NAA1, char iseq2[], int len1, int len2,
					   WorkingBuffer &buffer, int &best_sum,
					   int band_width, int &band_left, int &band_center, int &band_right, int required_aa1);
int local_band_align(char query[], char ref[], int qlen, int rlen, ScoreMatrix &mat,
					 int &best_score, int &iden_no, int &alnln, float &dist, unsigned int *alninfo,
					 int band_left, int band_center, int band_right, WorkingBuffer &buffer);

void strrev(const char *p);
int print_usage(const char *arg);
int print_usage_2d(const char *arg);
int print_usage_est(const char *arg);
int print_usage_div(const char *arg);
int print_usage_est_2d(const char *arg);
int print_usage_454(const char *arg);

void cal_aax_cutoff(double &aa1_cutoff, double &aa2_cutoff, double &aan_cutoff,
					double NR_clstr, int tolerance, int naa_stat_start_percent,
					int naa_stat[5][61][4], int NAA);
void update_aax_cutoff(double &aa1_cutoff, double &aa2_cutoff, double &aan_cutoff,
					   int tolerance, int naa_stat_start_percent,
					   int naa_stat[5][61][4], int NAA, int iden);

int calc_ann_list(int len, char *seqi, int NAA, int &aan_no, Vector<int> &aan_list, Vector<INTs> &aan_list_no, bool est = false);

float current_time();

// some functions from very old cd-hit
int quick_sort_idx(int *a, int *idx, int lo0, int hi0);
int quick_sort_idxr(int *a, int *idx, int lo0, int hi0);

void setaa_to_na();
void make_comp_short_word_index(size_t NAA);
void make_comp_iseq(int len, char *iseq_comp, char *iseq);
