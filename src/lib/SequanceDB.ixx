export module SequenceDB;

#define assert(condition) ((void)0)

import std;

import common;
import Options;
import ScoreMatrix;
import Sequence;
import WordTable;
import WorkingParam;
import gzfstream;

#define Z_LARGE64 1
#define Z_WANT64 1


#define omp_set_num_threads(T) void(T = T)
#define omp_get_thread_num() 0


using HitFileStream::HitFStream;
using HitFileStream::LINE_BUF_SIZE;
using HitFileStream::MAX_LINE_SIZE;
using std::cout;
using std::endl;
using std::fclose;
using std::feof;
using std::fgets;
using std::fopen;
using std::fprintf;
using std::fread;
using std::fseek;
using std::ftell;
using std::fwrite;
using std::printf;
using std::sprintf;
using std::strcmp;
using std::string;
using std::strlen;
using std::uint32_t;
using std::uint64_t;
using std::uint8_t;

constexpr uint8_t naa_stat_start_percent = 40;
export std::vector<int> Comp_AAN_idx;

template <class T>
struct tsv_view {
    const T& ref;
};

template <class T>
struct std::formatter<tsv_view<T>> : std::range_formatter<typename T::value_type> {
    constexpr formatter() {
        this->set_brackets("", "");
        this->set_separator("\t");
    }
};

constexpr uint8_t get_naa_stat(size_t i, size_t j, size_t k);

constexpr int calc_ann_list(const std::string_view seqi, size_t& aan_no, VectorInt& aan_list, VectorIntX& aan_list_no, bool est) {
    // check_aan_list
    aan_no = seqi.size() - options.NAA + 1;
    for (auto j = 0; j < aan_no; j++) {
        aan_list[j] = 0;
        for (size_t k = 0, k1 = options.NAA - 1; k < options.NAA; k++, k1--)
            aan_list[j] += seqi[j + k] * NAAN_array[k1];
    }
    if (est) {
        // for the short word containing 'N', mask it to '-1'
        for (auto j = 0; j < seqi.size(); j++) {
            if (seqi[j] >= 4) { // here N is 4
                auto i0 = (j - options.NAA + 1 > 0) ? j - options.NAA + 1 : 0;
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


constexpr int upper_bound_length_rep(int len, double opt_s, int opt_S, double opt_aL, int opt_AL) {
    size_t len_upper_bound = 99999999;
    double r1 = (opt_s > opt_aL) ? opt_s : opt_aL;
    int a2 = (opt_S < opt_AL) ? opt_S : opt_AL;
    if (r1 > 0.0)
        len_upper_bound = static_cast<size_t>(len / r1);
    if (len + a2 < len_upper_bound)
        len_upper_bound = len + a2;

    return len_upper_bound;
}


constexpr int upper_bound_length_rep(int len) {
    double opt_s = options.diff_cutoff;
    int opt_S = options.diff_cutoff_aa;
    double opt_aL = options.long_coverage;
    int opt_AL = options.long_control;
    return upper_bound_length_rep(len, opt_s, opt_S, opt_aL, opt_AL);
}


// when alignment coverage such as -aL is specified
// if a existing rep is too long, it won't be qulified

constexpr size_t MemoryLimit(const size_t mem_need) {
    size_t mem_limit = (options.max_memory - mem_need) / sizeof(IndexCount);

    // printf( "Table limit with the given memory limit:\n" );
    if (options.max_memory == 0) {
        mem_limit = options.max_entries;
        if (mem_limit > MAX_TABLE_SIZE)
            mem_limit = MAX_TABLE_SIZE;
    }
    // printf( "Max number of representatives: %zu\n", mem_limit );
    // printf( "Max number of word counting entries: %zu\n\n", mem_limit );
    return mem_limit;
}


////For smiple len1 <= len2, len2 is for existing representative
////walk along all diag path of two sequences,
////find the diags with most aap
////return top n diags
////added on 2006 11 13
////band 0                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                            XXXXXXXXXXXXXXX                  seq1
////band 1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                             XXXXXXXXXXXXXXX                 seq1
////extreme right (+)           XXXXXXXXXXXXXXXXXX               seq2, rep seq
////    band = len2-1                            XXXXXXXXXXXXXXX seq1
////band-1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                           XXXXXXXXXXXXXXX                   seq1
////extreme left (-)            XXXXXXXXXXXXXXXXXX               seq2, rep seq
////              XXXXXXXXXXXXXXX   band = -(len1-1)             seq1
////index of diag_score = band+len1-1;
int diag_test_aapn(int NAA1, char iseq2[], int len1, int len2, WorkingBuffer& buffer, int& best_sum, int band_width, int& band_left, int& band_center,
                   int& band_right, int required_aa1) {
    int i1, j, k;
    size_t nall = len1 + len2 - 1;

    if (nall > StaticOptions::MAX_DIAG)
        throw std::invalid_argument("in diag_test_aapn, MAX_DIAG reached");
    for (auto* pp = &buffer.diag_score[0], i = nall; i; i--, pp++)
        *pp = 0;
    for (auto* pp = &buffer.diag_score2[0], i = nall; i; i--, pp++)
        *pp = 0;

    int c22, cpx;
    // INTs *bip;
    int len11 = len1 - 1;
    int len22 = len2 - 1;
    i1 = len11;
    for (auto i = 0; i < len22; i++, i1++) {
        c22 = iseq2[i] * NAA1 + iseq2[i + 1];
        cpx = 1 + (iseq2[i] != iseq2[i + 1]);
        if ((j = buffer.taap[c22]) == 0)
            continue;
        int m = buffer.aap_begin[c22];
        for (int k = 0; k < j; k++) {
            buffer.diag_score[i1 - buffer.aap_list[m + k]]++;
            buffer.diag_score2[i1 - buffer.aap_list[m + k]] += cpx;
        }
    }

    // find the best band range
    //   int band_b = required_aa1;
    int band_b = required_aa1 - 1 >= 0 ? required_aa1 - 1 : 0; // on dec 21 2001
    int band_e = nall - band_b;

    int band_m = (band_b + band_width - 1 < band_e) ? band_b + band_width - 1 : band_e;
    int best_score = 0;
    int best_score2 = 0;
    int max_diag2 = 0;
    int imax_diag = 0;

#ifdef USE_AVX2
    __m256i best_score_vec = _mm256_setzero_si256();
    __m256i best_score2_vec = _mm256_setzero_si256();
    __m256i max_diag2_vec = _mm256_setzero_si256();
    int i = band_b;
    for (; i <= band_m - 7; i += 8) {
        __m256i diag_score_vec = _mm256_loadu_si256((__m256i const*)&buffer.diag_score[i]);
        __m256i diag_score2_vec = _mm256_loadu_si256((__m256i const*)&buffer.diag_score2[i]);
        best_score_vec = _mm256_add_epi32(best_score_vec, diag_score_vec);
        best_score2_vec = _mm256_add_epi32(best_score2_vec, diag_score2_vec);
        max_diag2_vec = _mm256_max_epi32(max_diag2_vec, diag_score2_vec);
    }
    int temp_best_score[8], temp_best_score2[8], temp_max_diag2[8];
    _mm256_storeu_si256((__m256i*)temp_best_score, best_score_vec);
    _mm256_storeu_si256((__m256i*)temp_best_score2, best_score2_vec);
    _mm256_storeu_si256((__m256i*)temp_max_diag2, max_diag2_vec);
    for(int k=0; k<8; ++k) {
        best_score += temp_best_score[k];
        best_score2 += temp_best_score2[k];
        if(temp_max_diag2[k] > max_diag2) {
            max_diag2 = temp_max_diag2[k];
        }
    }

    // Find imax_diag
    for (int k=0; k<band_m - band_b + 1; ++k) {
        if(buffer.diag_score2[band_b+k] == max_diag2) {
            imax_diag = band_b+k;
            break;
        }
    }

    for (; i <= band_m; i++) {
        best_score += buffer.diag_score[i];
        best_score2 += buffer.diag_score2[i];
        if (buffer.diag_score2[i] > max_diag2) {
            max_diag2 = buffer.diag_score2[i];
            imax_diag = i;
        }
    }
#else
    for (auto i = band_b; i <= band_m; i++) {
        best_score += buffer.diag_score[i];
        best_score2 += buffer.diag_score2[i];
        if (buffer.diag_score2[i] > max_diag2) {
            max_diag2 = buffer.diag_score2[i];
            imax_diag = i;
        }
    }
#endif
    int from = band_b;
    int end = band_m;
    int score = best_score;
    int score2 = best_score2;
    for (k = from, j = band_m + 1; j < band_e; j++, k++) {
        score -= buffer.diag_score[k];
        score += buffer.diag_score[j];
        score2 -= buffer.diag_score2[k];
        score2 += buffer.diag_score2[j];
        if (score2 > best_score2) {
            from = k + 1;
            end = j;
            best_score = score;
            best_score2 = score2;
            if (buffer.diag_score2[j] > max_diag2) {
                max_diag2 = buffer.diag_score2[j];
                // max_diag = diag_score[j];
                imax_diag = j;
            }
        }
    }
    int mlen = imax_diag;
    if (imax_diag > len1)
        mlen = nall - imax_diag;
    int emax = int((1.0 - options.cluster_thd) * mlen) + 1;
    for (j = from; j < imax_diag; j++) { // if aap pairs fail to open gap
        if ((imax_diag - j) > emax || buffer.diag_score[j] < 1) {
            best_score -= buffer.diag_score[j];
            from++;
        }
        else
            break;
    }
    for (j = end; j > imax_diag; j--) { // if aap pairs fail to open gap
        if ((j - imax_diag) > emax || buffer.diag_score[j] < 1) {
            best_score -= buffer.diag_score[j];
            end--;
        }
        else
            break;
    }

    //  delete [] diag_score;
    band_left = from - len1 + 1;
    band_right = end - len1 + 1;
    band_center = imax_diag - len1 + 1;
    best_sum = best_score;
    return OK_FUNC;
}
// END diag_test_aapn


int diag_test_aapn_est(int NAA1, const char iseq2[], int len1, int len2, WorkingBuffer& buffer, int& best_sum, int band_width, int& band_left, int& band_center,
                       int& band_right, int required_aa1) {
    size_t nall = len1 + len2 - 1;
    int NAA2 = NAA1 * NAA1;
    int NAA3 = NAA2 * NAA1;

    if (nall > StaticOptions::MAX_DIAG)
        throw std::invalid_argument("in diag_test_aapn_est, MAX_DIAG reached");
    auto* pp = &buffer.diag_score[0];
    auto* pp2 = &buffer.diag_score2[0];
    for (auto i = nall; i; i--, pp++, pp2++)
        *pp = *pp2 = 0;

    INTs* bip;
    int c22, cpx;
    int len22 = len2 - 3;
    auto i1 = len1 - 1;
    for (auto i = 0; i < len22; i++, i1++, iseq2++) {
        unsigned char c0 = iseq2[0];
        unsigned char c1 = iseq2[1];
        unsigned char c2 = iseq2[2];
        unsigned char c3 = iseq2[3];
        if ((c0 >= 4) || (c1 >= 4) || (c2 >= 4) || (c3 >= 4))
            continue; // skip N

        c22 = c0 * NAA3 + c1 * NAA2 + c2 * NAA1 + c3;
        auto j = buffer.taap[c22];
        if (j == 0)
            continue;
        cpx = 1 + int(c0 != c1) + int(c1 != c2) + (c2 != c3);
        bip = &buffer.aap_list[buffer.aap_begin[c22]]; //    bi = aap_begin[c22];
        for (; j; j--, bip++) {
            buffer.diag_score[i1 - *bip]++;
            buffer.diag_score2[i1 - *bip] += cpx;
        }
    }
#if 0
	int mmax = 0;
	int immax = 0;
	for (i = 0; i <= nall; i++) {
		if (i % len2 == 0 or i == nall) printf("\n");
		printf("%3i ", diag_score[i]);
		if (diag_score[i] > mmax) {
			mmax = diag_score[i];
			immax = i;
		}
	}
#endif

    // find the best band range
    //   int band_b = required_aa1;
    size_t band_b = required_aa1 >= 1 ? required_aa1 - 1 : 0; // on dec 21 2001
    assert(nall > band_b);
    size_t band_e = nall - band_b;

    if (options.is454) {
        band_b = len1 - band_width;
        band_e = len1 + band_width;
        if (band_b < 0)
            band_b = 0;
        if (band_e > nall)
            band_e = nall;
    }

    auto band_m = (band_b + band_width - 1 < band_e) ? band_b + band_width - 1 : band_e;
    int best_score = 0;
    int best_score2 = 0;
    int max_diag2 = 0;
    int imax_diag = 0;

#ifdef USE_AVX2
    __m256i best_score_vec = _mm256_setzero_si256();
    __m256i best_score2_vec = _mm256_setzero_si256();
    __m256i max_diag2_vec = _mm256_setzero_si256();
    size_t i = band_b;
    for (; i <= band_m - 7; i += 8) {
        __m256i diag_score_vec = _mm256_loadu_si256((__m256i const*)&buffer.diag_score[i]);
        __m256i diag_score2_vec = _mm256_loadu_si256((__m256i const*)&buffer.diag_score2[i]);
        best_score_vec = _mm256_add_epi32(best_score_vec, diag_score_vec);
        best_score2_vec = _mm256_add_epi32(best_score2_vec, diag_score2_vec);
        max_diag2_vec = _mm256_max_epi32(max_diag2_vec, diag_score2_vec);
    }
    int temp_best_score[8], temp_best_score2[8], temp_max_diag2[8];
    _mm256_storeu_si256((__m256i*)temp_best_score, best_score_vec);
    _mm256_storeu_si256((__m256i*)temp_best_score2, best_score2_vec);
    _mm256_storeu_si256((__m256i*)temp_max_diag2, max_diag2_vec);
    for(int k=0; k<8; ++k) {
        best_score += temp_best_score[k];
        best_score2 += temp_best_score2[k];
        if(temp_max_diag2[k] > max_diag2) {
            max_diag2 = temp_max_diag2[k];
        }
    }
    for (; i <= band_m; i++) {
        best_score += buffer.diag_score[i];
        best_score2 += buffer.diag_score2[i];
        if (buffer.diag_score2[i] > max_diag2) {
            max_diag2 = buffer.diag_score2[i];
            imax_diag = i;
        }
    }
    // Find imax_diag
    for (size_t k=0; k<band_m - band_b + 1; ++k) {
        if(buffer.diag_score2[band_b+k] == max_diag2) {
            imax_diag = band_b+k;
            break;
        }
    }
#else
    for (auto i = band_b; i <= band_m; i++) {
        best_score += buffer.diag_score[i];
        best_score2 += buffer.diag_score2[i];
        if (buffer.diag_score2[i] > max_diag2) {
            max_diag2 = buffer.diag_score2[i];
            imax_diag = i;
        }
    }
#endif
    int from = band_b;
    int end = band_m;
    int score = best_score;
    int score2 = best_score2;

    for (size_t k = from, j = band_m + 1; j < band_e; j++, k++) {
        score -= buffer.diag_score[k];
        score += buffer.diag_score[j];
        score2 -= buffer.diag_score2[k];
        score2 += buffer.diag_score2[j];
        if (score2 > best_score2) {
            from = k + 1;
            end = j;
            best_score = score;
            best_score2 = score2;
            if (buffer.diag_score2[j] > max_diag2) {
                max_diag2 = buffer.diag_score2[j];
                // max_diag = diag_score[j];
                imax_diag = j;
            }
        }
    }
#if 0
	printf("%i %i\n", required_aa1, from);
	printf("max=%3i  imax=%3i; band:  %3i  %3i  %i\n", max_diag, imax_diag, band_b, band_e, band_m);
	printf("best: %i\n", best_score);
	printf("from: %i, end: %i,  best: %i\n", from, end, best_score);
#endif
    int mlen = imax_diag;
    if (imax_diag > len1)
        mlen = nall - imax_diag;
    int emax = int((1.0 - options.cluster_thd) * mlen) + 1;
    for (auto j = from; j < imax_diag; j++) { // if aap pairs fail to open gap
        if ((imax_diag - j) > emax || buffer.diag_score[j] < 1) {
            best_score -= buffer.diag_score[j];
            from++;
        }
        else
            break;
    }
    for (auto j = end; j > imax_diag; j--) { // if aap pairs fail to open gap
        if ((j - imax_diag) > emax || buffer.diag_score[j] < 1) {
            best_score -= buffer.diag_score[j];
            end--;
        }
        else
            break;
    }

    band_left = from - len1 + 1;
    band_right = end - len1 + 1;
    band_center = imax_diag - len1 + 1;
    best_sum = best_score;
    if (options.is454) {
        if (band_left > 0)
            best_sum = 0;
        if (band_right < 0)
            best_sum = 0;
    }
#if 0
	printf("%3i:  best: %i,  %i  %i  %i\n", required_aa1, best_score, band_left, band_right, band_width);
	printf("max=%3i  imax=%3i; band:  %3i  %3i  %i\n", mmax, immax, band_b, band_e, band_m);
#endif
    return OK_FUNC;
}
// END diag_test_aapn_est

/*
local alignment of two sequence within a diag band
for band 0 means direction (0,0) -> (1,1)
         1 means direction (0,1) -> (1,2)
        -1 means direction (1,0) -> (2,1)
added on 2006 11 13
band 0                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                            XXXXXXXXXXXXXXX                  seq1
band 1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                             XXXXXXXXXXXXXXX                 seq1
extreme right (+)           XXXXXXXXXXXXXXXXXX               seq2, rep seq
    band = len2-1                            XXXXXXXXXXXXXXX seq1
band-1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                           XXXXXXXXXXXXXXX                   seq1
extreme left (-)            XXXXXXXXXXXXXXXXXX               seq2, rep seq
              XXXXXXXXXXXXXXX   band = -(len1-1)             seq1
iseq len are integer sequence and its length,
mat is matrix, return ALN_PAIR class

       band:  -101   seq2 len2 = 17
                \\\1234567890123456
              0  \xxxxxxxxxxxxxxxxx
              1   xxxxxxxxxxxxxxxxx\ most right band = len2-1
              2   xxxxxxxxxxxxxxxxx
    seq1      3   xxxxxxxxxxxxxxxxx
    len1 = 11 4   xxxxxxxxxxxxxxxxx
              5   xxxxxxxxxxxxxxxxx
              6   xxxxxxxxxxxxxxxxx
              7   xxxxxxxxxxxxxxxxx
              8   xxxxxxxxxxxxxxxxx
              9   xxxxxxxxxxxxxxxxx
              0   xxxxxxxxxxxxxxxxx
                  \
                   most left band = -(len1-1)

*/

int local_band_align(const char iseq1[], const char iseq2[], int ilen1, int ilen2, int& best_score, int& iden_no, int& alnln, float& dist, unsigned int* alninfo, int band_left, int band_center, int band_right, WorkingBuffer& buffer) {
    size_t best_score1;
    iden_no = 0;

    assert(ilen1 > 0 && ilen2 > 0);

    if ((band_right >= ilen2) || (band_left <= -ilen1) || (band_left > band_right))
        return FAILED_FUNC;

    auto len1 = static_cast<size_t>(ilen1);
    auto len2 = static_cast<size_t>(ilen2);

    // allocate mem for score_mat[len1][len2] etc
    int band_width = band_right - band_left + 1;
    int band_width1 = band_width + 1;

    // score_mat, back_mat [i][j]: i index of seqi (0 to len(seqi)-1), j index of band (0 to band_width-1)
    MatrixInt64& score_mat = buffer.score_mat;
    MatrixInt& back_mat = buffer.back_mat;

    // printf( "%i  %i\n", band_right, band_left );

    if (score_mat.size() <= len1) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= len1) {
            score_mat.push_back(row2);
            back_mat.push_back(row);
        }
    }
    for (size_t i = 0; i <= len1; i++) {
        if (score_mat[i].size() < band_width1)
            score_mat[i].resize(band_width1);
        if (back_mat[i].size() < band_width1)
            back_mat[i].resize(band_width1);
    }

    best_score = 0;
    /* seq1 is query, seq2 is rep
                  seq2    len2 = 17       seq2    len2 = 17    seq2    len2 = 17
                  01234567890123456       01234567890123456    01234567890123456
       0          xxxxxxxxxxxxxxxxx \\\\\\XXXxxxxxxxxxxxxxx    xXXXXXXXxxxxxxxxx
       1     \\\\\Xxxxxxxxxxxxxxxxx  \\\\\Xxx\xxxxxxxxxxxxx    xx\xxxxx\xxxxxxxx
       2      \\\\X\xxxxxxxxxxxxxxx   \\\\Xxxx\xxxxxxxxxxxx    xxx\xxxxx\xxxxxxx
  seq1 3       \\\Xx\xxxxxxxxxxxxxx    \\\Xxxxx\xxxxxxxxxxx    xxxx\xxxxx\xxxxxx
  len1 4        \\Xxx\xxxxxxxxxxxxx     \\Xxxxxx\xxxxxxxxxx    xxxxx\xxxxx\xxxxx
  = 11 5         \Xxxx\xxxxxxxxxxxx      \Xxxxxxx\xxxxxxxxx    xxxxxx\xxxxx\xxxx
       6          Xxxxx\xxxxxxxxxxx       Xxxxxxxx\xxxxxxxx    xxxxxxx\xxxxx\xxx
       7          x\xxxx\xxxxxxxxxx       x\xxxxxxx\xxxxxxx    xxxxxxxx\xxxxx\xx
       8          xx\xxxx\xxxxxxxxx       xx\xxxxxxx\xxxxxx    xxxxxxxxx\xxxxx\x
       9          xxx\xxxx\xxxxxxxx       xxx\xxxxxxx\xxxxx    xxxxxxxxxx\xxxxx\
       0          xxxx\xxxx\xxxxxxx       xxxx\xxxxxxx\xxxx    xxxxxxxxxxx\xxxxx
                  band_left < 0 (-6)      band_left < 0 (-6)   band_left >=0 (1)
                  band_right < 0 (-1)     band_right >=0 (2)   band_right >=0(7)
                  band_width 6            band_width 9         band_width 7
       init score_mat, and iden_mat (place with upper 'X')
     */

    if (band_left < 0) { // set score to left border of the matrix within band
        auto tband = (band_right < 0) ? band_right : 0;
        // for (k=band_left; k<tband; k++)
        for (auto k = band_left; k <= tband; k++) { // fixed on 2006 11 14
            auto i = -k;
            auto j1 = k - band_left;
            // penalty for leading gap opening = penalty for gap extension
            // each of the left side query hunging residues give ext_gap (-1)
            score_mat[i][j1] = mat.get_extra_gap() * i;
            back_mat[i][j1] = DP_BACK_TOP;
        }
        back_mat[-tband][tband - band_left] = DP_BACK_NONE;
    }

    if (band_right >= 0) { // set score to top border of the matrix within band
        size_t tband = (band_left > 0) ? band_left : 0;
        for (size_t j = tband; j <= band_right; j++) {
            size_t j1 = j - band_left;
            score_mat[0][j1] = mat.get_extra_gap() * j;
            back_mat[0][j1] = DP_BACK_LEFT;
        }
        back_mat[0][tband - band_left] = DP_BACK_NONE;
    }

    int gap_open[2] = {mat.get_gap(), mat.get_extra_gap()};
    int max_diag = band_center - band_left;
    int extra_score[4] = {4, 3, 2, 1};

#ifdef USE_AVX2
    for (size_t i = 1; i <= len1; i++) {
        int J0 = 1 - band_left - i;
        int J1 = len2 - band_left - i;
        if (J0 < 0) J0 = 0;
        if (J1 >= band_width) J1 = band_width;

        int j1 = J0;
        for (; j1 <= J1 - 8; j1 += 8) {
            int j = j1 + i + band_left;

            // 1. Load substitution scores
            int s_scores[8];
            for(int k=0; k<8; ++k) {
                int ci = iseq1[i - 1];
                int cj = iseq2[j + k - 1];
                s_scores[k] = mat.matrix[ci][cj];
                size_t raw_dist = std::abs((j1+k) - max_diag);
                int extra = extra_score[raw_dist & 3];
                s_scores[k] += extra * (s_scores[k] > 0);
            }
            __m256i sij_vec = _mm256_loadu_si256((__m256i const*)s_scores);

            // 2. Calculate diagonal scores
            __m256i diag_scores = _mm256_loadu_si256((__m256i const*)&score_mat[i - 1][j1]);
            diag_scores = _mm256_add_epi32(diag_scores, sij_vec);

            // 3. Calculate scores from top
            int gap0 = gap_open[(i == len1) | (j == len2)];
            __m256i gap_vec = _mm256_set1_epi32(gap0);
            __m256i top_scores = _mm256_loadu_si256((__m256i const*)&score_mat[i - 1][j1 + 1]);
            top_scores = _mm256_add_epi32(top_scores, gap_vec);

            // 4. Calculate scores from left
            __m256i left_scores = _mm256_loadu_si256((__m256i const*)&score_mat[i][j1 - 1]);
            left_scores = _mm256_add_epi32(left_scores, gap_vec);

            // 5. Find max
            __m256i best_scores = _mm256_max_epi32(diag_scores, _mm256_max_epi32(top_scores, left_scores));

            // 6. Store scores
            _mm256_storeu_si256((__m256i*)&score_mat[i][j1], best_scores);

            // 7. Update back pointers (simplified)
            __m256i back_vec = _mm256_setzero_si256();
            __m256i cmp_diag = _mm256_cmpeq_epi32(best_scores, diag_scores);
            __m256i cmp_left = _mm256_cmpeq_epi32(best_scores, left_scores);
            __m256i cmp_top = _mm256_cmpeq_epi32(best_scores, top_scores);

            back_vec = _mm256_or_si256(back_vec, _mm256_and_si256(cmp_diag, _mm256_set1_epi32(DP_BACK_LEFT_TOP)));
            back_vec = _mm256_or_si256(back_vec, _mm256_and_si256(cmp_left, _mm256_set1_epi32(DP_BACK_LEFT)));
            back_vec = _mm256_or_si256(back_vec, _mm256_and_si256(cmp_top, _mm256_set1_epi32(DP_BACK_TOP)));
            _mm256_storeu_si256((__m256i*)&back_mat[i][j1], back_vec);
        }

        // Process remaining elements serially
        for (; j1 <= J1; j1++) {
            int j = j1 + i + band_left;

            int ci = iseq1[i - 1];
            int cj = iseq2[j - 1];
            int sij = mat.matrix[ci][cj];

            /* extra score according to the distance to the best diagonal */
            size_t raw_dist = std::abs(j1 - max_diag);
            int extra = extra_score[raw_dist & 3]; // max distance 3
            sij += extra * (sij > 0);

            int back = DP_BACK_LEFT_TOP;
            best_score1 = score_mat[i - 1][j1] + sij;
            int gap0 = gap_open[(i == len1) | (j == len2)];

            if (j1 > 0) {
                int gap = gap0;
                if (back_mat[i][j1 - 1] == DP_BACK_LEFT)
                    gap = mat.get_extra_gap();
                size_t score = score_mat[i][j1 - 1] + gap;
                if (score > best_score1) {
                    back = DP_BACK_LEFT;
                    best_score1 = score;
                }
            }
            if (j1 + 1 < band_width) {
                int gap = gap0;
                if (back_mat[i - 1][j1 + 1] == DP_BACK_TOP)
                    gap = mat.get_extra_gap();
                size_t score = score_mat[i - 1][j1 + 1] + gap;
                if (score > best_score1) {
                    back = DP_BACK_TOP;
                    best_score1 = score;
                }
            }
            score_mat[i][j1] = best_score1;
            back_mat[i][j1] = back;
        }
    }
#else
    for (size_t i = 1; i <= len1; i++) {
        int J0 = 1 - band_left - i;
        int J1 = len2 - band_left - i;
        if (J0 < 0) J0 = 0;
        if (J1 >= band_width) J1 = band_width;

        for (int j1 = J0; j1 <= J1; j1++) {
            int j = j1 + i + band_left;
            int ci = iseq1[i - 1];
            int cj = iseq2[j - 1];
            int sij = mat.matrix[ci][cj];

            size_t raw_dist = std::abs(j1 - max_diag);
            int extra = extra_score[raw_dist & 3];
            sij += extra * (sij > 0);

            int back = DP_BACK_LEFT_TOP;
            best_score1 = score_mat[i - 1][j1] + sij;
            int gap0 = gap_open[(i == len1) | (j == len2)];

            if (j1 > 0) {
                int gap = gap0;
                if (back_mat[i][j1 - 1] == DP_BACK_LEFT) gap = mat.get_extra_gap();
                size_t score = score_mat[i][j1 - 1] + gap;
                if (score > best_score1) {
                    back = DP_BACK_LEFT;
                    best_score1 = score;
                }
            }
            if (j1 + 1 < band_width) {
                int gap = gap0;
                if (back_mat[i - 1][j1 + 1] == DP_BACK_TOP) gap = mat.get_extra_gap();
                size_t score = score_mat[i - 1][j1 + 1] + gap;
                if (score > best_score1) {
                    back = DP_BACK_TOP;
                    best_score1 = score;
                }
            }
            score_mat[i][j1] = best_score1;
            back_mat[i][j1] = back;
        }
    }
#endif
#else
    for (size_t i = 1; i <= len1; i++) {
        int J0 = 1 - band_left - i;
        int J1 = len2 - band_left - i;
        if (J0 < 0)
            J0 = 0;
        if (J1 >= band_width)
            J1 = band_width;
        for (int j1 = J0; j1 <= J1; j1++) {
            int j = j1 + i + band_left;

            int ci = iseq1[i - 1];
            int cj = iseq2[j - 1];
            int sij = mat.matrix[ci][cj];

            /* extra score according to the distance to the best diagonal */
            size_t raw_dist = std::abs(j1 - max_diag);
            int extra = extra_score[raw_dist & 3]; // max distance 3
            sij += extra * (sij > 0);

            int back = DP_BACK_LEFT_TOP;
            best_score1 = score_mat[i - 1][j1] + sij;
            int gap0 = gap_open[(i == len1) | (j == len2)];

            if (j1 > 0) {
                int gap = gap0;
                if (back_mat[i][j1 - 1] == DP_BACK_LEFT)
                    gap = mat.get_extra_gap();
                size_t score = score_mat[i][j1 - 1] + gap;
                if (score > best_score1) {
                    back = DP_BACK_LEFT;
                    best_score1 = score;
                }
            }
            if (j1 + 1 < band_width) {
                int gap = gap0;
                if (back_mat[i - 1][j1 + 1] == DP_BACK_TOP)
                    gap = mat.get_extra_gap();
                size_t score = score_mat[i - 1][j1 + 1] + gap;
                if (score > best_score1) {
                    back = DP_BACK_TOP;
                    best_score1 = score;
                }
            }
            score_mat[i][j1] = best_score1;
            back_mat[i][j1] = back;
        }
    }
#endif
    size_t i = 0;
    size_t j = 0;
    if (len2 - band_left < len1) {
        i = len2 - band_left;
        j = len2;
    }
    else if (len1 + band_right < len2) {
        i = len1;
        j = len1 + band_right;
    }
    else {
        i = len1;
        j = len2;
    }
    size_t j1 = j - i - band_left;
    best_score = score_mat[i][j1];
    best_score1 = score_mat[i][j1];

    // #if 1
    // 	const char *letters = "acgtnx";
    // 	const char *letters2 = "ACGTNX";
    // #else
    // 	const char *letters = "arndcqeghilkmfpstwyvbzx";
    // 	const char *letters2 = "ARNDCQEGHILKMFPSTWYVBZX";
    // #endif

    int back = back_mat[i][j1];
    int last = back;
    int count = 0, count2 = 0, count3 = 0;
    int match, begin1, begin2, end1, end2;
    int gbegin1 = 0, gbegin2 = 0, gend1 = 0, gend2 = 0;
    std::int64_t score, smin = best_score1, smax = best_score1 - 1;
    int posmin, posmax, pos = 0;
    int bl, dlen = 0, dcount = 0;
    posmin = posmax = 0;
    begin1 = begin2 = end1 = end2 = 0;

#ifdef PRINT
#define PRINT
    printf("%i %i\n", best_score, score_mat[i][j1]);
    printf("%i %i %i\n", band_left, band_center, band_right);
    printf("%i %i %i %i\n", i, j, j1, len2);
#endif
#ifdef MAKEALIGN
#define MAKEALIGN
    char AA[MAX_SEQ], BB[MAX_SEQ];
    int NN = 0;
    int IA, IB;
    for (IA = len1; IA > i; IA--) {
        AA[NN] = letters[iseq1[IA - 1]];
        BB[NN++] = '-';
    }
    for (IB = len2; IB > j; IB--) {
        AA[NN] = '-';
        BB[NN++] = letters[iseq2[IB - 1]];
    }
#endif

    int masked = 0;
    int indels = 0;
    int max_indels = 0;
    while (back != DP_BACK_NONE) {
        switch (back) {
            case DP_BACK_TOP:
#ifdef PRINT
                printf("%5i: %c %c %9i\n", pos, letters[iseq1[i - 1]], '|', score_mat[i][j1]);
#endif
#ifdef MAKEALIGN
                AA[NN] = letters[iseq1[i - 1]];
                BB[NN++] = '-';
#endif
                bl = (last != back) & (j != 1) & (j != len2);
                dlen += bl;
                dcount += bl;
                score = score_mat[i][j1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                i -= 1;
                j1 += 1;
                break;
            case DP_BACK_LEFT:
#ifdef PRINT
                printf("%5i: %c %c %9i\n", pos, '|', letters[iseq2[j - 1]], score_mat[i][j1]);
#endif
#ifdef MAKEALIGN
                AA[NN] = '-';
                BB[NN++] = letters[iseq2[j - 1]];
#endif
                bl = (last != back) & (i != 1) & (i != len1);
                dlen += bl;
                dcount += bl;
                score = score_mat[i][j1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                j1 -= 1;
                j -= 1;
                break;
            case DP_BACK_LEFT_TOP:
#ifdef PRINT
                if (iseq1[i - 1] == iseq2[j - 1]) {
                    printf("%5i: %c %c %9i\n", pos, letters2[iseq1[i - 1]], letters2[iseq2[j - 1]], score_mat[i][j1]);
                }
                else {
                    printf("%5i: %c %c %9i\n", pos, letters[iseq1[i - 1]], letters[iseq2[j - 1]], score_mat[i][j1]);
                }
#endif
#ifdef MAKEALIGN
                if (iseq1[i - 1] == iseq2[j - 1]) {
                    AA[NN] = letters2[iseq1[i - 1]];
                    BB[NN++] = letters2[iseq2[j - 1]];
                }
                else {
                    AA[NN] = letters[iseq1[i - 1]];
                    BB[NN++] = letters[iseq2[j - 1]];
                }
#endif
                if (alninfo && options.global_identity) {
                    if (i == 1 || j == 1) {
                        gbegin1 = i - 1;
                        gbegin2 = j - 1;
                    }
                    else if (i == len1 || j == len2) {
                        gend1 = i - 1;
                        gend2 = j - 1;
                    }
                }
                score = score_mat[i][j1];
                i -= 1;
                j -= 1;
                match = iseq1[i] == iseq2[j];
                if (score > smax) {
                    count = 0;
                    smax = score;
                    posmax = pos;
                    end1 = i;
                    end2 = j;
                }
                if (options.isEST && (iseq1[i] > 4 || iseq2[j] > 4)) {
                    masked += 1;
                }
                else {
                    dlen += 1;
                    dcount += (match != 0);
                    count += match;
                    count2 += match;
                    count3 += match;
                }
                if (score < smin) {
                    int mm = match == 0;
                    count2 = 0;
                    smin = score;
                    posmin = pos - mm;
                    begin1 = i + mm;
                    begin2 = j + mm;
                }
                break;
            default:
                printf("%i\n", back);
                break;
        }
        if (options.is454) {
            if (back == DP_BACK_LEFT_TOP) {
                if (indels > max_indels)
                    max_indels = indels;
                indels = 0;
            }
            else {
                if (last == DP_BACK_LEFT_TOP) {
                    indels = 1;
                }
                else if (indels) {
                    indels += 1;
                }
            }
        }
        pos += 1;
        last = back;
        back = back_mat[i][j1];
    }
    if (options.is454 and max_indels > options.max_indel)
        return FAILED_FUNC;
    iden_no = options.global_identity ? count3 : count - count2;
    alnln = posmin - posmax + 1 - masked;
    dist = dcount / (float)dlen;
    // dist = - 0.75 * log( 1.0 - dist * 4.0 / 3.0 );
    size_t umtail1 = len1 - 1 - end1;
    size_t umtail2 = len2 - 1 - end2;
    size_t umhead = begin1 < begin2 ? begin1 : begin2;
    size_t umtail = umtail1 < umtail2 ? umtail1 : umtail2;
    size_t umlen = umhead + umtail;
    if (umlen > options.unmatch_len)
        return FAILED_FUNC;
    if (umlen > len1 * options.short_unmatch_per)
        return FAILED_FUNC;
    if (umlen > len2 * options.long_unmatch_per)
        return FAILED_FUNC;
    if (alninfo) {
        alninfo[0] = begin1;
        alninfo[1] = end1;
        alninfo[2] = begin2;
        alninfo[3] = end2;
        alninfo[4] = masked;
        if (options.global_identity) {
            alninfo[0] = gbegin1;
            alninfo[1] = gend1;
            alninfo[2] = gbegin2;
            alninfo[3] = gend2;
        }
    }
#ifdef PRINT
    printf("%6i %6i:  %4i %4i %4i %4i\n", alnln, iden_no, begin1, end1, begin2, end2);
    printf("%6i %6i:  %4i %4i\n", posmin, posmax, posmin - posmax, count - count2);
    printf("smin = %9i, smax = %9i\n", smin, smax);
    printf("dlen = %5i, dcount = %5i, dist = %.3f\n", dlen, dcount, dcount / (float)dlen);
#endif
#ifdef MAKEALIGN
    float identity = iden_no / (float)(options.global_identity ? (len1 - masked) : alnln);
    if (identity < options.cluster_thd)
        return OK_FUNC;
    while (i--) {
        AA[NN] = letters[iseq1[i - 1]];
        BB[NN++] = '-';
    }
    while (j--) {
        AA[NN] = '-';
        BB[NN++] = letters[iseq2[j - 1]];
    }
    AA[NN] = '\0';
    BB[NN] = '\0';
    for (i = 0; i < NN / 2; i++) {
        char aa = AA[i], bb = BB[i];
        AA[i] = AA[NN - i - 1];
        BB[i] = BB[NN - i - 1];
        AA[NN - i - 1] = aa;
        BB[NN - i - 1] = bb;
    }
    static int fcount = 0;
    fcount += 1;
    FILE* fout = fopen("alignments.txt", "a");
    if (fout == nullptr) {
        if (fcount <= 1)
            printf("alignment files open failed\n");
        return OK_FUNC;
    }
    fprintf(fout, "\n\n######################################################\n");
    fprintf(fout, "# length X = %i\n", len2);
    fprintf(fout, "# length Y = %i\n", len1);
    fprintf(fout, "# best align X: %i-%i\n", begin2 + 1, end2 + 1);
    fprintf(fout, "# best align Y: %i-%i\n", begin1 + 1, end1 + 1);
    if (alninfo) {
        fprintf(fout, "# align X: %i-%i\n", alninfo[2] + 1, alninfo[3] + 1);
        fprintf(fout, "# align Y: %i-%i\n", alninfo[0] + 1, alninfo[1] + 1);
    }
    fprintf(fout, "# alignment length: %i\n", alnln);
    fprintf(fout, "# identity count: %i\n", iden_no);
    fprintf(fout, "# identity: %g\n", identity);
    fprintf(fout, "# distance: %g\n", dist);
    if (options.is454)
        fprintf(fout, "# max indel: %i\n", max_indels);
#if 0
	fprintf(fout, "%i %s\n", seq1->index, AA);
	fprintf(fout, "%i %s\n", seq2->index, BB);
#else
    bool printaa = true;
    IB = IA = 0;
    fprintf(fout, "\n\nX ");
    while (IA < NN) {
        if (printaa) {
            fprintf(fout, "%c", BB[IB]);
            IB += 1;
            if (IB % 75 == 0 or IB == NN)
                printaa = false, fprintf(fout, "\nY ");
        }
        else {
            fprintf(fout, "%c", AA[IA]);
            IA += 1;
            if (IA % 75 == 0)
                printaa = true, fprintf(fout, "\n\nX ");
        }
    }
#endif
    fclose(fout);
#endif

    return OK_FUNC;
} // END int local_band_align


void cal_aax_cutoff(double& aa1_cutoff, double& aa2_cutoff, double& aan_cutoff, double cluster_thd, int tolerance) {
    aa1_cutoff = cluster_thd;
    aa2_cutoff = 1 - (1 - cluster_thd) * 2;
    aan_cutoff = 1 - (1 - cluster_thd) * options.NAA;
    if (tolerance == 0)
        return;

    int clstr_idx = (int)(cluster_thd * 100) - naa_stat_start_percent;
    if (clstr_idx < 0)
        clstr_idx = 0;
    double d2 = get_naa_stat(tolerance - 1, clstr_idx, 3) / 100.0f;
    double dn = get_naa_stat(tolerance - 1, clstr_idx, 5 - options.NAA) / 100.0f;
    aa2_cutoff = d2 > aa2_cutoff ? d2 : aa2_cutoff;
    aan_cutoff = dn > aan_cutoff ? dn : aan_cutoff;
}

void update_aax_cutoff(double& aa1_cutoff, double& aa2_cutoff, double& aan_cutoff, int tolerance, double cluster_thd) {
    cluster_thd = std::min(1.0, cluster_thd);

    double aa1_t, aa2_t, aan_t;
    cal_aax_cutoff(aa1_t, aa2_t, aan_t, cluster_thd, tolerance);
    if (aa1_t > aa1_cutoff)
        aa1_cutoff = aa1_t;
    if (aa2_t > aa2_cutoff)
        aa2_cutoff = aa2_t;
    if (aan_t > aan_cutoff)
        aan_cutoff = aan_t;
}


constexpr void make_comp_iseq(int len, char* iseq_comp, const char* iseq) {
    std::array<char, 8> c{3, 2, 1, 0, 4, 5};
    for (auto i = 0; i < len; i++)
        iseq_comp[i] = c[iseq[len - i - 1]];
}


// liwz 2016 01, seq is too short for the shortest (longer) seq in word_table to satisfy -aL option
// longer seqeunce * -aL -band_width
constexpr bool check_covrage(const WordTable& table, size_t len) {
    if (!table.sequences.empty()) {
        auto cov_size = static_cast<unsigned int>(options.long_coverage * table.sequences.back()->size);
        unsigned int min_red = std::min(cov_size - options.band_width, 0u);
        if (len < min_red)
            return false;
    }
    return true;
}


constexpr auto copy_sequence_chunked(std::istream& fin, std::ostream& fout, size_t to_copy) -> void {
    std::array<char, MAX_LINE_SIZE> buf;
    while (to_copy > 0) {
        size_t to_read = std::min(to_copy, buf.size());
        if (!fin.read(buf.data(), to_read))
            throw std::runtime_error("Failed to read input");
        if (!fout.write(buf.data(), to_read))
            throw std::runtime_error("Failed to write output");
        to_copy -= to_read;
    }
}

typedef std::vector<Sequence> seqvec;
static_assert(std::movable<Sequence>);

export class SequenceDB {
public:
    size_t NAAN = 0;
    std::vector<Sequence> sequences;
    std::vector<size_t> rep_seqs;
    std::map<size_t, std::vector<Sequence>> clusters;

    size_t total_letter = 0;
    size_t total_desc = 0;
    size_t max_len = 0;
    size_t min_len = 0;
    size_t len_n50 = 0;

    void Clear() {
        sequences.clear();
        rep_seqs.clear();
    }

    SequenceDB() = default; // default constructor
    ~SequenceDB() { Clear(); }


    void DoClustering(int T) {
        int i, k;
        double aa1_cutoff = options.cluster_thd;
        double aas_cutoff = 1 - (1 - options.cluster_thd) * 4;
        double aan_cutoff = 1 - (1 - options.cluster_thd) * options.NAA;
        int seq_no = sequences.size();
        int frag_no = seq_no;
        int frag_size = options.frag_size;
        // int len, len_bound;
        // int flag;
        std::valarray<size_t> letters(T);

        // printf( "%li\n".mem_limit );

        if (frag_size) {
            frag_no = 0;
            for (i = 0; i < seq_no; i++)
                frag_no += (sequences[i].size - options.NAA) / frag_size + 1;
        }

        if (not options.isEST)
            cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd, options.tolerance);

        std::vector params{T, WorkingParam(aa1_cutoff, aas_cutoff, aan_cutoff)};
        std::vector buffers(T, WorkingBuffer(frag_no, max_len));

        // word_table as self comparing table and table buffer:
        WordTable word_table(options.NAA, NAAN);

        WordTable last_table(options.NAA, NAAN);

        int N = sequences.size();
        // int K = N - 100 * T;
        size_t mem_need = MinimalMemory(frag_no, buffers[0].total_bytes, T);
        // TODO: size_t mem_limit = MemoryLimit( mem_need );
        size_t mem;
        size_t tabsize = 0;
        int remaining = 0;

        options.ComputeTableLimits(min_len, max_len, len_n50, mem_need);

        // ##########
        // PARALELIZE
        // ##########
        for (i = 0; i < N;) {
            int start = i;
            int m = i;
            size_t sum = remaining;
            float redundancy = (rep_seqs.size() + 1.0) / (i + 1.0);
            size_t max_items = options.max_entries;
            size_t max_seqs = options.max_sequences;
            size_t items = 0;
            if (i == 0 && max_seqs > 1000) { // first SCB with small size
                max_items /= 8;
                max_seqs /= 8;
            }
            while (m < N && (sum * redundancy) < max_seqs && items < max_items) {
                Sequence& seq = sequences[m];
                if (!(seq.state & seq_state::IS_REDUNDANT)) {
                    // items += seq.size;
                    items += static_cast<size_t>(seq.size * redundancy);
                    sum += 1;
                }
                m++;
            }
            if ((m > i + 1E4) && (m > i + (N - i) / (2 + T)))
                m = i + (N - i) / (2 + T);
            if (m == i || m >= N) {
                m = N;
                if (m > i + 1E3)
                    m = i + (N - i) / (2 + T);
            }
            // printf( "m = %i  %i,  %i\n", i, m, m-i );
            printf("\r# comparing sequences from  %9i  to  %9i\n", i, m);
            if (last_table.size) {
                int print = (m - i) / 20 + 1;
                // #pragma omp parallel for schedule( dynamic, 1 )
                for (int j = i; j < m; j++) {
                    Sequence& seq = sequences[j];
                    if (seq.state & seq_state::IS_REDUNDANT)
                        continue;
                    int tid = omp_get_thread_num();
                    CheckOne(seq, last_table, params[tid], buffers[tid]);
                    if (j % print == 0) {
                        std::cout << ".";
                        std::cout.flush();
                    }
                }
                int may_stop = 0;
                int self_stop = 0;
                float p0 = 0;
                int min = last_table.sequences[last_table.sequences.size() - 1]->size;
                int m0 = m;
                bool stop = false;
                // #pragma omp parallel for schedule( dynamic, 1 )
                for (int j = m - 1; j < N; j++) {
                    // #pragma omp flush (stop)
                    if (!stop) {
                        if (j + 1 == N)
                            may_stop = 1;
                        if (j == (m0 - 1)) { // use m0 to avoid other iterations satisfying the condition:
                            int tid = omp_get_thread_num();
                            for (int ks = i; ks < m; ks++) {
                                Sequence& seq = sequences[ks];
                                i = ks + 1;
                                if (seq.state & seq_state::IS_REDUNDANT)
                                    continue;
                                ClusterOne(seq, ks, word_table, params[tid], buffers[tid]);
                                if (may_stop and word_table.sequences.size() >= 100)
                                    break;
                                if (word_table.size >= max_items)
                                    break;
                                int tmax = max_seqs - (frag_size ? seq.size / frag_size + 1 : 0);
                                if (word_table.sequences.size() >= tmax)
                                    break;
                            }
                            self_stop = 1;
                        }
                        else {
                            Sequence& seq = sequences[j];
                            if (seq.state & seq_state::IS_REDUNDANT)
                                continue;
                            int tid = omp_get_thread_num();
                            CheckOne(seq, last_table, params[tid], buffers[tid]);
                            if (min > params[tid].len_upper_bound) {
                                may_stop = 1;
                                stop = true;
                            }
                            if (self_stop && tid == 1) {
                                float p = (100.0 * j) / N;
                                if (p > p0 + 1E-1) { // print only if the percentage changed
                                    printf("\r%4.1f%%", p);
                                    p0 = p;
                                }
                            }
                        }
                    }
                }
            }
            if (i == start || m == N) {
                // printf( "comparing the first or last or very small group ...\n" ); fflush( stdout );
                for (k = i; k < m;) {
                    int mm = k;
                    auto sum = 0;
                    while (mm < m && sum < 1E5) {
                        if (!(sequences[mm].state & seq_state::IS_REDUNDANT))
                            sum += sequences[mm].size;
                        mm += 1;
                    }
                    if (mm < k + 1000)
                        mm = k + 1000;
                    if (mm > m)
                        mm = m;
                    // #pragma omp parallel for schedule( dynamic, 1 )
                    for (auto kk = k; kk < mm; kk++) {
                        Sequence& seq = sequences[kk];
                        if (seq.state & seq_state::IS_REDUNDANT)
                            continue;
                        int tid = omp_get_thread_num();
                        CheckOne(seq, word_table, params[tid], buffers[tid]);
                    }
                    bool bk = false;
                    for (int ks = k; ks < mm; ks++) {
                        Sequence& seq = sequences[ks];
                        i = k = ks + 1;
                        if (seq.state & seq_state::IS_REDUNDANT)
                            continue;
                        ClusterOne(seq, ks, word_table, params[0], buffers[0]);
                        bk = true;
                        if (word_table.size >= max_items)
                            break;
                        int tmax = max_seqs - (frag_size ? seq.size / frag_size + 1 : 0);
                        if (word_table.sequences.size() >= tmax)
                            break;
                        bk = false;
                    }
                    if (bk)
                        break;
                }
            }
            else if (i < m) {
                remaining = remaining / 2 + (m - i);
                printf("\r---------- %6i remaining sequences to the next cycle\n", m - i);
            }
            std::println("---------- new table with {:8} representatives", word_table.sequences.size());
            if ((last_table.size + word_table.size) > tabsize)
                tabsize = last_table.size + word_table.size;
            last_table.Clear();
            last_table.sequences.swap(word_table.sequences);
            last_table.indexCounts.swap(word_table.indexCounts);
            last_table.size = word_table.size;
            word_table.size = 0;
        }
        std::println("\n{:9}  finished  {:9}  clusters", sequences.size(), rep_seqs.size());
        mem = (mem_need + tabsize * sizeof(IndexCount)) / MEGA_MiBi;
        std::println("\nApproximated maximum memory consumption: {}M", mem);
        last_table.Clear();
        word_table.Clear();
    }


    static int CheckOne(Sequence& seq, WordTable& table, WorkingParam& param, WorkingBuffer& buf) {
        int len = seq.size;
        param.len_upper_bound = upper_bound_length_rep(len);
        auto ret = (options.isEST) ? CheckOneEST(seq, table, param, buf) : CheckOneAA(seq, table, param, buf);
        return ret;
    }


    static int CheckOneAA(Sequence& seq, WordTable& table, WorkingParam& param, WorkingBuffer& buf) {
        auto aa1_cutoff = param.aa1_cutoff;
        auto aa2_cutoff = param.aas_cutoff;
        auto aan_cutoff = param.aan_cutoff;

        const auto& seqi = seq.get_data();
        auto flag = false;
        auto frag_size = options.frag_size;

        auto S = table.sequences.size();
        auto len_eff = seqi.size();

        if (S) {
            unsigned int min = table.sequences[S - 1]->size;
            if (min < seqi.size()) {
                if (seqi.size() * options.diff_cutoff2 > min)
                    min = (int)(seqi.size() * options.diff_cutoff2);
                if ((seqi.size() - options.diff_cutoff_aa2) > min)
                    min = seqi.size() - options.diff_cutoff_aa2;
                len_eff = min;
            }
        }

        // liwz 2016 01, seq is too short for the shortest (longer) seq in word_table to satisfy -aL option
        // longer seqeunce * -aL -band_width
        if (S) {
            auto cov_size = static_cast<unsigned int>(options.long_coverage * table.sequences[S - 1]->size);
            unsigned int min_red = std::min(cov_size - options.band_width, 0u);
            if (seqi.size() < min_red)
                return 0;
        }

        param.ControlShortCoverage(len_eff);
        param.ComputeRequiredBases(options.NAA, 2);

        buf.EncodeWords(seq, options.NAA, false);

        // if minimal alignment length > len, return
        // I can not return earlier, because I need to calc the word_encodes etc
        if (options.min_control > seqi.size())
            return 0;

        // lookup_aan
        int aan_no = seqi.size() - options.NAA + 1;
        // int M = frag_size ? table.frag_count : S;
        table.CountWords(aan_no, buf, false, param.required_aan);

        // contained_in_old_lib()
        unsigned int len_upper_bound = param.len_upper_bound;
        unsigned int len_lower_bound = param.len_lower_bound;
        int band_left, band_right, best_score, best_sum, alnln;
        int tiden_no, band_center;
        float distance = 0;
        int frg2 = frag_size ? (seqi.size() - options.NAA + options.band_width) / frag_size + 1 + 1 : 0;
        int has_aa2 = 0;

        for (auto& ic : buf.lookCounts) {
            unsigned int talign_info[5];
            if (!frag_size) {
                buf.indexMapping[ic.index] = 0;
                if (ic.count < param.required_aan)
                    continue;
            }

            Sequence* rep = table.sequences[ic.index];
            auto len2 = rep->size;
            if (len2 > len_upper_bound)
                continue;
            if (options.has2D && len2 < len_lower_bound)
                continue;
            if (frag_size) {
                auto* ims = &buf.indexMapping[ic.index];
                auto k = (len2 - options.NAA) / frag_size + 1;
                auto sum = 0;
                for (auto j1 = 0; j1 < frg2; j1++) {
                    auto im = ims[j1];
                    if (im)
                        sum += buf.lookCounts[im - 1].count;
                }
                auto count = sum;
                for (auto j1 = frg2; j1 < k; j1++) {
                    auto im1 = ims[j1];
                    auto im2 = ims[j1 - frg2];
                    if (im1)
                        sum += buf.lookCounts[im1 - 1].count;
                    if (im2)
                        sum -= buf.lookCounts[im2 - 1].count;
                    if (sum > count)
                        count = sum;
                }
                if (count < param.required_aan)
                    continue;
            }

            param.ControlLongCoverage(len2);

            if (has_aa2 == 0) { // calculate AAP array
                buf.ComputeAAP(seqi);
                has_aa2 = 1;
            }
            auto* seqj = rep->data; // NR_seq[NR90_idx[j]];
            const auto len = seq.size;
            int band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
            diag_test_aapn(NAA1, seqj, len, len2, buf, best_sum, band_width1, band_left, band_center, band_right, param.required_aa1);
            if (best_sum < param.required_aas)
                continue;

            int rc = FAILED_FUNC;
            if (options.print || param.aln_cover_flag) // return overlap region
                rc = local_band_align(seqi.data(), seqj, len, len2, best_score, tiden_no, alnln, distance, talign_info, band_left, band_center, band_right, buf);
            else
                rc = local_band_align(seqi.data(), seqj, len, len2, best_score, tiden_no, alnln, distance, talign_info, band_left, band_center, band_right, buf);
            if (rc == FAILED_FUNC)
                continue;
            if (tiden_no < param.required_aa1)
                continue;
            auto lens = len;
            if (options.has2D && len > len2)
                lens = len2;
            int len_eff1 = (options.global_identity == 0) ? alnln : (lens - talign_info[4]);
            float tiden_pc = tiden_no / static_cast<float>(len_eff1);
            if (options.useDistance()) {
                if (distance > options.distance_thd)
                    continue;
                if (distance >= seq.distance)
                    continue; // existing distance
            }
            else {
                if (tiden_pc < options.cluster_thd)
                    continue;
                if (tiden_pc <= seq.identity)
                    continue; // existing iden_no
            }
            if (param.aln_cover_flag) {
                if (talign_info[3] - talign_info[2] + 1 < param.min_aln_lenL)
                    continue;
                if (talign_info[1] - talign_info[0] + 1 < param.min_aln_lenS)
                    continue;
            }
            if (options.has2D)
                seq.state |= seq_state::IS_REDUNDANT;
            // Todo: why do we flag?
            flag = true;
            seq.identity = tiden_pc;
            seq.cluster_id = rep->cluster_id;
            seq.distance = distance;
            seq.coverage[0] = talign_info[0] + 1;
            seq.coverage[1] = talign_info[1] + 1;
            seq.coverage[2] = talign_info[2] + 1;
            seq.coverage[3] = talign_info[3] + 1;
            if (not options.cluster_best)
                break;
            update_aax_cutoff(aa1_cutoff, aa2_cutoff, aan_cutoff, options.tolerance, tiden_pc);
            param.ComputeRequiredBases(options.NAA, 2);

        }
        auto ic = buf.lookCounts.end();
        if (frag_size)
            ic = buf.lookCounts.begin();
        for (;ic != buf.lookCounts.end(); ++ic)
            buf.indexMapping[ic->index] = 0;
        buf.lookCounts.clear();
        if (flag) { // if similar to old one delete it
            if (!options.cluster_best) {
                seq.Clear();
                seq.state |= seq_state::IS_REDUNDANT;
            }
        }
        return flag;
    }


    static int CheckOneEST(Sequence& seq, WordTable& table, WorkingParam& param, WorkingBuffer& buf) {
        char* seqi_comp = &buf.seqi_comp[0];

        char* seqi = seq.data;
        size_t len = seq.size;
        int flag = 0;
        int S = table.sequences.size();
        int len_eff = len;
        if (S) {
            size_t min = table.sequences[S - 1]->size;
            if (min < len) {
                if (len * options.diff_cutoff2 > min)
                    min = static_cast<size_t>(len * options.diff_cutoff2);
                if ((len - options.diff_cutoff_aa2) > min)
                    min = len - options.diff_cutoff_aa2;
                len_eff = min;
            }
        }


        if (!check_covrage(table, len))
            return 0;


        param.ControlShortCoverage(len_eff);
        param.ComputeRequiredBases(options.NAA, 4);
        int skip = buf.EncodeWords(seq, options.NAA, true);
        param.required_aan -= skip;
        param.required_aas -= skip;
        param.required_aa1 -= skip;
        if (param.required_aan <= 0)
            param.required_aan = 1;
        if (param.required_aas <= 0)
            param.required_aas = 1;
        if (param.required_aa1 <= 0)
            param.required_aa1 = 1;

        // if minimal alignment length > len, return
        // I can not return earlier, because I need to calc the word_encodes etc
        if (options.min_control > len)
            return 0; // return flag=0

        int aan_no = len - options.NAA + 1;

        // contained_in_old_lib()
        unsigned int len_upper_bound = param.len_upper_bound;
        unsigned int len_lower_bound = param.len_lower_bound;
        int band_left, band_right, best_score, best_sum, alnln;
        int tiden_no, band_center;
        float distance = 0;
        unsigned int talign_info[5];

        for (auto comp = 0; comp < 2; comp++) {
            if (comp) {
                for (auto j0 = 0; j0 < aan_no; j0++) {
                    auto j = buf.word_encodes[j0];
                    if (j < 0)
                        buf.aan_list_comp[j0] = j;
                    else
                        buf.aan_list_comp[j0] = Comp_AAN_idx[j];
                }

                make_comp_iseq(len, seqi_comp, seqi);
                seqi = seqi_comp;
            }
            int has_aas = 0;

            if (comp) {
                table.CountWords(aan_no, buf, true, param.required_aan);
            }
            else {
                table.CountWords(aan_no, buf, true, param.required_aan);
            }

            auto ic = buf.lookCounts.begin();
            for (; ic->count; ic++) {
                buf.indexMapping[ic->index] = 0;
                if (ic->count < param.required_aan)
                    continue;
                Sequence* rep = table.sequences[ic->index];

                auto len2 = rep->size;
                if (len2 > len_upper_bound)
                    continue;
                if (options.has2D && len2 < len_lower_bound)
                    continue;

                auto seqj = rep->get_data();

                param.ControlLongCoverage(len2);

                if (has_aas == 0) { // calculate AAP array
                    buf.ComputeAAP2(seqj);
                    has_aas = 1;
                }

                auto band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
                diag_test_aapn_est(NAA1, seqj.data(), len, len2, buf, best_sum, band_width1, band_left, band_center, band_right, param.required_aa1);
                if (best_sum < param.required_aas)
                    continue;
                // if( comp and flag and (not options.cluster_best) and j > rep->cluster_id ) goto Break;

                int rc = FAILED_FUNC;
                if (options.print || param.aln_cover_flag) { // return overlap region
                    rc = local_band_align(
                        seqi, seqj.data(), len, len2, best_score, tiden_no, alnln, distance, talign_info, band_left, band_center, band_right, buf);
                    if (comp) {
                        talign_info[0] = len - talign_info[0] - 1;
                        talign_info[1] = len - talign_info[1] - 1;
                    }
                }
                else {
                    // printf( "%5i %5i %5i %5i\n", band_width1, band_right-band_left, band_left, band_right );
                    rc = local_band_align(
                        seqi, seqj.data(), len, len2, best_score, tiden_no, alnln, distance, talign_info, band_left, band_center, band_right, buf);
                }
                if (rc == FAILED_FUNC)
                    continue;
                // printf( "%i  %i  %i\n", best_score, tiden_no, required_aa1 );
                if (tiden_no < param.required_aa1)
                    continue;
                if (options.is454) {
                    if (talign_info[2] != talign_info[0])
                        continue; // same start
                    if (talign_info[0] > 1)
                        continue; // one mismatch allowed at beginning
                    if ((len - talign_info[1]) > 2)
                        continue; // one mismatch allowed at end
                }

                auto lens = len;
                if (options.has2D && len > len2)
                    lens = len2;
                int len_eff1 = (options.global_identity == 0) ? alnln : (lens - talign_info[4]);
                auto tiden_pc = tiden_no / static_cast<float>(len_eff1);
                // printf( "%i %f\n", tiden_no, tiden_pc );
                if (options.useDistance()) {
                    if (distance > options.distance_thd)
                        continue;
                    if (options.cluster_best and distance >= seq.distance)
                        continue; // existing distance
                }
                else {
                    if (tiden_pc < options.cluster_thd)
                        continue;
                    if (options.cluster_best and tiden_pc < seq.identity)
                        continue; // existing iden_no
                }
                if (param.aln_cover_flag) {
                    if (talign_info[3] - talign_info[2] + 1 < param.min_aln_lenL)
                        continue;
                    if (comp) {
                        if (talign_info[0] - talign_info[1] + 1 < param.min_aln_lenS)
                            continue;
                    }
                    else {
                        if (talign_info[1] - talign_info[0] + 1 < param.min_aln_lenS)
                            continue;
                    }
                }
                if (options.cluster_best and std::fabs(tiden_pc - seq.identity) < 1E-9 and rep->cluster_id >= seq.cluster_id)
                    continue;
                if ((not options.cluster_best) and flag != 0 and rep->cluster_id >= seq.cluster_id)
                    continue;
                flag = comp ? -1 : 1;
                seq.identity = tiden_pc;
                seq.distance = distance;
                seq.cluster_id = rep->cluster_id;
                seq.coverage[0] = talign_info[0] + 1;
                seq.coverage[1] = talign_info[1] + 1;
                seq.coverage[2] = talign_info[2] + 1;
                seq.coverage[3] = talign_info[3] + 1;
                if (not options.cluster_best)
                    break;
            }
            while (ic->count) {
                buf.indexMapping[ic->index] = 0;
                ic += 1;
            }
            buf.lookCounts.clear();
            if (not options.option_r)
                break;
        }
        if ((flag == 1) || (flag == -1)) { // if similar to old one delete it
            if (!options.cluster_best) {
                seq.Clear();
                seq.state |= seq_state::IS_REDUNDANT;
            }
            if (flag == -1)
                seq.state |= seq_state::IS_MINUS_STRAND;
            else
                seq.state -= seq_state::IS_MINUS_STRAND;
        }
        return flag;
    }


    void ComputeDistance() const {
        int N = sequences.size();
        int best_score, best_sum;
        int band_width1, band_left, band_center, band_right;
        int tiden_no, alnln;
        float distance;
        WorkingBuffer buf(N, max_len);

        std::vector dists(N, std::vector<float>(N));

        Sequence comseq(sequences[0]);

        for (unsigned int i = 0; i < N; i++) {
            unsigned int talign_info[5];
            const Sequence& seq = sequences[i];
            const auto seqi = seq.get_data();
            unsigned int len = seq.size;
            buf.EncodeWords(seq, options.NAA, false);
            buf.ComputeAAP2(seqi);
            dists[i][i] = 0.0;
            if ((i + 1) % 1000 == 0)
                printf("%9i\n", (i + 1));
            for (unsigned int j = 0; j < i; j++) {
                const Sequence& rep = sequences[j];
                const char* rep_data = rep.get_data().data();
                unsigned int len2 = rep.size;
                band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
                diag_test_aapn_est(NAA1, rep_data, len, len2, buf, best_sum, band_width1, band_left, band_center, band_right, 0);
                local_band_align(seqi.data(), rep_data, len, len2, best_score, tiden_no, alnln, distance, talign_info, band_left, band_center, band_right, buf);
                dists[seq.index][rep.index] = dists[rep.index][seq.index] = distance;
            }
            if (not options.option_r)
                break;
            comseq.index = seq.index;
            comseq.size = len;
            for (unsigned int j = 0; j < len; j++)
                comseq.data[i] = seq.get_data().data()[len - i - 1];
            const auto seqi1 = comseq.data;
            buf.EncodeWords(comseq, options.NAA, false);
            buf.ComputeAAP2({seqi1, seq.size});
            for (unsigned int j = 0; j < i; j++) {
                const Sequence& rep = sequences[j];
                const char* seqj = rep.get_data().data();
                unsigned int len2 = rep.size;
                band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
                diag_test_aapn_est(NAA1, seqj, len, len2, buf, best_sum, band_width1, band_left, band_center, band_right, 0);
                local_band_align(seqi1, seqj, len, len2, best_score, tiden_no, alnln, distance, talign_info, band_left, band_center, band_right, buf);
                if (distance < dists[seq.index][rep.index])
                    dists[seq.index][rep.index] = dists[rep.index][seq.index] = distance;
            }
        }

        auto path = std::filesystem::path{options.output}.replace_extension("dist");
        std::ofstream sffout(path, std::ios::out | std::ios::trunc);

        auto iota_view = std::views::iota(0, N);
        std::println(sffout, "{}", iota_view);

        for (const auto& d : dists)
            std::println(sffout, "{}", d);
        sffout << std::flush;
    }

    void DoClustering() {
        double aa1_cutoff = options.cluster_thd;
        double aas_cutoff = 1 - (1 - options.cluster_thd) * 4;
        double aan_cutoff = 1 - (1 - options.cluster_thd) * options.NAA;
        size_t seq_no = sequences.size();
        size_t frag_no = seq_no;
        size_t frag_size = options.frag_size;
        // int len, len_bound;
        // int flag;

#if 0
		ComputeDistance();
		return;
#endif

        if (options.threads > 1) {
            DoClustering(options.threads);
            return;
        }

        if (frag_size) {
            frag_no = 0;
            for (size_t i = 0; i < seq_no; i++)
                frag_no += (sequences[i].size - options.NAA) / frag_size + 1;
        }

        if (not options.isEST)
            cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd, options.tolerance);

        WorkingParam param(aa1_cutoff, aas_cutoff, aan_cutoff);
        WorkingBuffer buffer(frag_no, max_len);

        WordTable word_table(options.NAA, NAAN);

        size_t mem_need = MinimalMemory(frag_no, buffer.total_bytes, 1);
        // TODO: size_t mem_limit = MemoryLimit( mem_need );
        size_t N = sequences.size();

        size_t total_letters = total_letter;
        size_t tabsize = 0;

        options.ComputeTableLimits(min_len, max_len, len_n50, mem_need);

        for (size_t i = 0; i < N;) {
            float redundancy = (rep_seqs.size() + 1.0) / (i + 1.0);
            size_t m = i;
            size_t sum = 0;
            size_t max_items = options.max_entries;
            size_t max_seqs = options.max_sequences;
            size_t items = 0;

            // find a block from i to m, so that this block can fit into a word table
            //     ...
            //  i  ++++++++++++++++++++++++++
            //     ++++++++++++++++++++
            //     ++++++++++++++++
            //  m  +++++++++++++
            //     ...
            while (m < N && (sum * redundancy) < max_seqs && items < max_items) {
                Sequence& seq = sequences[m];
                if (!(seq.state & seq_state::IS_REDUNDANT)) {
                    items += (size_t)(seq.size * redundancy);
                    sum += 1;
                }
                m++;
            }
            m = std::min(m, N);
            std::println("{}comparing sequences from  {:9}  to  {:9}", '\r', i, m);
            for (int ks = i; ks < m; ks++) { // clustering this block
                Sequence& seq = sequences[ks];
                i = ks + 1;
                if (seq.state & seq_state::IS_REDUNDANT)
                    continue;
                ClusterOne(seq, ks, word_table, param, buffer);
                total_letters -= seq.size;
                if (word_table.size >= max_items)
                    break;
                int tmax = max_seqs - (frag_size ? seq.size / frag_size + 1 : 0);
                if (word_table.sequences.size() >= tmax)
                    break;
            } // finishing word table from this block
            m = i;
            if (word_table.size == 0)
                continue;
            float p0 = 0;
            for (size_t j = m; j < N; j++) { // use this word table to screen rest sequences m->N
                Sequence& seq = sequences[j];
                if (seq.state & seq_state::IS_REDUNDANT)
                    continue;
                CheckOne(seq, word_table, param, buffer);
                total_letters -= seq.size;
                size_t len_bound = param.len_upper_bound;
                if (word_table.sequences[word_table.sequences.size() - 1]->size > len_bound) {
                    break;
                }
                float p = (100.0 * j) / N;
                if (p > p0 + 1E-1) { // print only if the percentage changed
                    printf("\r%4.1f%%", p);
                    p0 = p;
                }
            }
            if (word_table.size > tabsize)
                tabsize = word_table.size;
            // if( i && i < m ) printf( "\r---------- %6i remaining sequences to the next cycle\n", m-i );
            word_table.Clear();
        }
        auto mem = (mem_need + tabsize * sizeof(IndexCount)) / MEGA_MiBi;
        std::cout << std::endl;
        std::println("{}  finished  {}  clusters", sequences.size(), rep_seqs.size());
        std::cout << "Approximated maximum memory consumption: " << mem << "M";
        std::cout << std::endl;
        word_table.Clear();

#if 0
		int zeros = 0;
		for (i = 0; i < word_table.indexCounts.size(); i++) zeros += word_table.indexCounts[i].Size() == 0;
		printf("%9i  empty entries out of  %9i\n", zeros, word_table.indexCounts.size());
#endif
    }

    void ClusterTo(SequenceDB& other) {
        double aa1_cutoff = options.cluster_thd;
        double aas_cutoff = 1 - (1 - options.cluster_thd) * 4;
        double aan_cutoff = 1 - (1 - options.cluster_thd) * options.NAA;
        VectorInt word_encodes(MAX_SEQ);
        VectorIntX word_encodes_no(MAX_SEQ);

        if (not options.isEST) {
            cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd, options.tolerance);
        }

        auto N = other.sequences.size();
        auto M = sequences.size();
        auto T = options.threads;

        std::valarray<size_t> counts(T);
        std::vector<WorkingParam> params(T);
        std::vector<WorkingBuffer> buffers(T);
        WorkingBuffer& buffer = buffers[0];
        for (auto i = 0; i < T; i++) {
            params[i] = WorkingParam(aa1_cutoff, aas_cutoff, aan_cutoff);
            buffers[i].Set(N, max_len);
        }
        if (T > 1)
            omp_set_num_threads(T);

        size_t mem_need = MinimalMemory(N, buffer.total_bytes, T, other.total_letter + other.total_desc);
        // TODO: size_t mem_limit = MemoryLimit( mem_need );

        options.ComputeTableLimits(min_len, max_len, len_n50, mem_need);

        WordTable word_table(options.NAA, NAAN);

        size_t max_items = options.max_entries;
        size_t max_seqs = options.max_sequences;
        size_t NR2_red_no = 0;
        for (size_t i = 0; i < N;) {
            size_t items = 0;
            size_t sum = 0;
            size_t m = i;
            while (m < N && sum < max_seqs && items < max_items) {
                Sequence& seq = other.sequences[m];
                if (!(seq.state & seq_state::IS_REDUNDANT)) {
                    items += seq.size;
                    sum += 1;
                }
                m++;
            }
            m = std::min(m, N);
            for (size_t ks = i; ks < m; ks++) {
                Sequence& seq = other.sequences[ks];
                size_t out_ann_no = 0;
                calc_ann_list(seq.get_data(), out_ann_no, word_encodes, word_encodes_no, options.isEST);
                word_table.AddWordCounts(out_ann_no, word_encodes, word_encodes_no, ks - i);
                word_table.sequences.push_back(&seq);

                set_representative(ks, seq);

                auto ks1 = ks + 1;
                if (ks1 % 1000 == 0) {
                    std::cout << "." << std::flush;
                    if (ks1 % 10000 == 0)
                        std::println("{:9}  finished", ks1);
                }
            }
            double p0 = 0.0;
            auto JM = M;
            counts = 0;
            // #pragma omp parallel for schedule( dynamic, 1 )
            for (int j = 0; j < JM; j++) {
                Sequence& seq = sequences[j];
                if (seq.state & seq_state::IS_REDUNDANT)
                    continue;

                int tid = omp_get_thread_num();
                params[tid].len_upper_bound = upper_bound_length_rep(seq.size);
                params[tid].len_lower_bound = std::min<size_t>(seq.cutoff2_len(), seq.size - options.diff_cutoff_aa2);

                if (word_table.sequences[word_table.sequences.size() - 1]->size > params[tid].len_upper_bound) {
                    JM = 0;
                    continue;
                }

                auto flag = other.CheckOne(seq, word_table, params[tid], buffers[tid]);
                if ((flag == 1) || (flag == -1)) { // if similar to old one delete it
                    if (!options.cluster_best) {
                        seq.Clear();
                        seq.state |= seq_state::IS_REDUNDANT;
                        counts[tid]++;
                    }
                    if (flag == -1)
                        seq.state |= seq_state::IS_MINUS_STRAND; // for EST only
                }
                auto p = (100.0 * j) / M;
                if (p > p0 + 1E-1) { // print only if the percentage changed
                    std::cout << '\r' << std::fixed << std::setprecision(1) << std::setw(4) << p << '%' << std::flush;
                    p0 = p;
                }
            }
            for (int j = 0; j < T; j++)
                NR2_red_no += counts[j];
            std::println("\r..........{:9}  compared  {:9}  clusters", i, NR2_red_no);
            word_table.Clear();
            word_table.size = 0;
            i = m;
        }

        if (options.cluster_best) { // delete redundant sequences in options.cluster_best mode
            for (auto& seq : sequences) {
                if (seq.identity > 0) {
                    seq.state |= seq_state::IS_REDUNDANT;
                    NR2_red_no++;
                }
            }
        }
        for (auto i = 0; i < sequences.size(); i++) {
            Sequence& seq = sequences[i];
            if (seq.identity < 0)
                seq.identity *= -1;
            if (not(seq.state & seq_state::IS_REDUNDANT))
                rep_seqs.push_back(i);
        }

        std::cout << std::endl;
        std::cout << sequences.size() << " compared\t" << NR2_red_no << " clustered" << std::endl;
    }


    void set_representative(size_t id, Sequence& seq, bool is_first = false) {
        if (is_first) {
            size_t new_cluster_id = rep_seqs.size();
            rep_seqs.push_back(id);
            id = new_cluster_id;
            clusters.insert({id, {}});
            seq.identity = 0;
        }
        seq.cluster_id = id;
        clusters[id].push_back(seq);
        seq.state |= seq_state::IS_REP;
    }


    void ClusterOne(Sequence& seq, int id, WordTable& table, WorkingParam& param, WorkingBuffer& buffer) {
        if (seq.state & seq_state::IS_REDUNDANT)
            return;

        param.len_upper_bound = upper_bound_length_rep(seq.size);
        ;
        auto flag = CheckOne(seq, table, param, buffer);

        if (flag == 0) {
            if ((seq.identity > 0) && options.cluster_best) {
                // because of the -g option, this seq is similar to seqs in old SEGs
                seq.state |= seq_state::IS_REDUNDANT;
                seq.Clear();
            }
            else { // else add to NR90 db
                set_representative(id, seq, true);

                size_t len_nna = seq.size - options.NAA;
                size_t aan_no = len_nna + 1;

                /* not used for EST */
                if (options.frag_size) {
                    size_t num_frags = len_nna / options.frag_size + 1;
                    table.AddWordCountsFrag(aan_no, buffer.word_encodes_backup, buffer.word_encodes_no, num_frags);
                }
                else {
                    table.AddWordCounts(aan_no, buffer.word_encodes, buffer.word_encodes_no, table.sequences.size());
                }

                table.sequences.emplace_back(&seq);

                if (options.frag_size) {
                    while (table.sequences.size() < table.frag_count)
                        table.sequences.push_back(&seq);
                }
            }
        }
        auto id1 = id + 1;
        if (id1 % 1000 == 0) {
            std::cout << "." << std::flush;
            if (id1 % 10000 == 0)
                std::println("\r..........{:9}  finished  {:9}clusters", id1, rep_seqs.size());
        }
    }


    size_t MinimalMemory(int frag_no, int bsize, int T, size_t extra = 0) {
        int N = sequences.size();
        int F = frag_no < MAX_TABLE_SEQ ? frag_no : MAX_TABLE_SEQ;
        size_t mem_need = 0;
        int table = T > 1 ? 2 : 1;

        std::cout << std::endl << "Approximated minimal memory consumption:" << std::endl;
        auto mem = N * sizeof(Sequence) + total_desc + N + extra;
        if (options.store_disk == false)
            mem += total_letter + N;
        size_t mem_mb = mem / MEGA_MiBi;
        std::print("Sequence        : {}M\n", mem_mb);
        mem_need += mem;

        mem = bsize;
        std::print("Buffer          : {} X {}M = {}uM\n", T, mem_mb, T * mem_mb);
        mem_need += T * mem;

        mem = F * (sizeof(Sequence*) + sizeof(IndexCount)) + NAAN * sizeof(std::vector<IndexCount>);
        std::print("Table           : {} X {}M = {}M\n", table, mem_mb, table * mem_mb);
        mem_need += table * mem;

        mem = sequences.capacity() * sizeof(Sequence*) + N * sizeof(int);
        mem += Comp_AAN_idx.size() * sizeof(int);
        std::print("Miscellaneous   : {}M\n", mem_mb);
        mem_need += mem;

        size_t i = mem_need / MEGA_MiBi;
        std::print("Total           : {}M\n\n", i);

        auto limit = mem_need + 50 * table;
        if (options.max_memory and options.max_memory < limit) {
            size_t limit_mb = limit / MEGA_MiBi;
            auto msg = std::format("not enough memory, please set -M option greater than {}\n", limit_mb);
            throw std::invalid_argument(msg.c_str());
        }
        return mem_need;
    }

    // by liwz gzip version 2019-02
    // by liwz
    // disable swap option
    // change des_begin, des_length, des_length2, dat_length => des_begin, tot_length
    // where des_begin is the FILE pointer of sequence record start
    //       tot_length is the total bytes of sequence record
    void Read(const std::string_view filename) {
        Sequence id_builder;


        Clear();

        std::string buffer(2 << 10, '\0');
        auto fin = HitFStream(filename);
        while (!fin.eof() && !fin.fail()) {
            std::getline(fin, buffer);
            if (buffer.empty())
                continue;

            if (!(buffer[0] == '>' || buffer[0] == '@'))
                throw std::invalid_argument("Invalid sequence file format! found neither '>' nor '@'\n" + buffer);

            Sequence& read_one = sequences.emplace_back();
            read_one.set_id(std::move(buffer));
            buffer.clear();
            // we assume standard compliance we read the rest of the sequence
            std::getline(fin, buffer);
            while (!buffer.empty() && std::isalpha(buffer[0])) {
                read_one += buffer;
                std::getline(fin, buffer);
            }
            read_one.index = sequences.size();

            // If this is FASTQ
            if (buffer[0] == '+') {
                // Ignore this line `same sequence identifier (and any description) again`
                buffer.clear();
                std::getline(fin, buffer);
                // Ignore quality values
                buffer.clear();
            }

            // Validate
            if (read_one.FormatAndValidate()) {
                std::println("Warning: from file \"{}\" Discarding invalid sequence or sequence without identifier and description!", filename);
                sequences.pop_back();
                if (read_one.id_empty())
                    std::println("Seq:\t {}", read_one.get_identifier());
                else
                    std::println("ID:\t {}", read_one.get_data());
                std::cout.flush();
            }
        }
    }


    // by liwz gzip version 2019-02
    // PE reads liwz, disable swap option
    void Read(const std::string_view file, const std::string_view file2) {
        // Sequence one, two;
        // Sequence des;
        // const auto fin = HitFStream(file);
        // const auto fin2 = HitFStream(file2);
        // char* buffer = nullptr;
        // char* buffer2 = nullptr;
        // char* res = nullptr;
        // char* res2 = nullptr;
        // size_t option_l = options.min_length;

        // Clear();

        // buffer = new char[MAX_LINE_SIZE + 1];
        // buffer2 = new char[MAX_LINE_SIZE + 1];
        //
        // while (((not gzeof(fin)) && (not gzeof(fin2))) || (one.size && two.size)) { /* do not break when the last sequence is not handled */
        // 	buffer[0] = '>'; res = gzgets(fin, buffer, MAX_LINE_SIZE);
        // 	buffer2[0] = '>'; res2 = gzgets(fin2, buffer2, MAX_LINE_SIZE);
        //
        // 	if ((res == nullptr) && (res2 != nullptr)) throw std::invalid_argument("Paired input files have different number sequences");
        // 	if ((res != nullptr) && (res2 == nullptr)) throw std::invalid_argument("Paired input files have different number sequences");
        // 	if ((one.size == 0) && (two.size > 0)) throw std::invalid_argument("Paired input files have different number sequences");
        // 	if ((one.size > 0) && (two.size == 0)) throw std::invalid_argument("Paired input files have different number sequences");
        // 	if ((res == nullptr) && (one.size == 0)) break;
        //
        // 	if (buffer[0] == '+') { // fastq 3rd line
        // 		// file 1
        // 		int len = strlen(buffer);
        // 		int len2 = len;
        // 		while (len2 && buffer[len2 - 1] != '\n') { // read until the end of the line
        // 			if ((res = gzgets(fin, buffer, MAX_LINE_SIZE)) == nullptr) break;
        // 			len2 = strlen(buffer);
        // 			len += len2;
        // 		}
        // 		one.tot_length += len;
        //
        // 		// read next line quality score
        // 		if ((res = gzgets(fin, buffer, MAX_LINE_SIZE)) == nullptr) throw std::invalid_argument("can not read quality score after");
        // 		len = strlen(buffer);
        // 		len2 = len;
        // 		while (len2 && buffer[len2 - 1] != '\n') {
        // 			if ((res = gzgets(fin, buffer, MAX_LINE_SIZE)) == nullptr) break;
        // 			len2 = strlen(buffer);
        // 			len += len2;
        // 		}
        // 		one.tot_length += len;
        //
        // 		// file 2
        // 		len = strlen(buffer2);
        // 		len2 = len;
        // 		while (len2 && buffer2[len2 - 1] != '\n') { // read until the end of the line
        // 			if ((res2 = gzgets(fin2, buffer2, MAX_LINE_SIZE)) == nullptr) break;
        // 			len2 = strlen(buffer2);
        // 			len += len2;
        // 		}
        // 		two.tot_length += len;
        //
        // 		// read next line quality score
        // 		if ((res2 = gzgets(fin2, buffer2, MAX_LINE_SIZE)) == nullptr) throw std::invalid_argument("can not read quality score after");
        // 		len = strlen(buffer2);
        // 		len2 = len;
        // 		while (len2 && buffer2[len2 - 1] != '\n') {
        // 			if ((res2 = gzgets(fin2, buffer2, MAX_LINE_SIZE)) == nullptr) break;
        // 			len2 = strlen(buffer2);
        // 			len += len2;
        // 		}
        // 		two.tot_length += len;
        //
        // 	}
        // 	else if (buffer[0] == '>' || buffer[0] == '@' || (res == nullptr && one.size)) {
        // 		if (one.size && two.size) { // write previous record
        // 			if (one.id_empty() || one.Format()) {
        // 				std::println("Warning: from file \"{}\",", file);
        // 				std::println("Discarding invalid sequence or sequence without identifier and description!\n");
        // 				if (!one.id_empty()) std::println(one.get_identifier());
        // 				std::println( one.data);
        // 				one.size = 0; two.size = 0;
        // 			}
        // 			if (two.id_empty() || two.Format()) {
        // 				std::println("Warning: from file \"{}\",", file2);
        // 				std::println("Discarding invalid sequence or sequence without identifier and description!\n");
        // 				if (!two.id_empty()) std::println(two.get_identifier());
        // 				std::println(std::string_view{two.data});
        // 				one.size = 0; two.size = 0;
        // 			}
        // 			one.index = sequences.size();
        // 			if ((one.size + two.size) > option_l) {
        // 				if (options.trim_len > 0) one.trim(options.trim_len);
        // 				if (options.trim_len_R2 > 0) two.trim(options.trim_len_R2);
        // 				sequences.push_back(new Sequence(one, two, 1));
        // 			}
        // 		}
        // 		// R1
        // 		one.size = 0;
        // 		one.tot_length = 0;
        //
        // 		int len = strlen(buffer);
        // 		int len2 = len;
        // 		des.size = 0;
        // 		des += buffer;
        // 		while (len2 && buffer[len2 - 1] != '\n') {
        // 			if ((res = gzgets(fin, buffer, MAX_LINE_SIZE)) == nullptr) break;
        // 			des += buffer;
        // 			len2 = strlen(buffer);
        // 			len += len2;
        // 		}
        // 		size_t offset = gztell(fin);
        // 		one.des_begin = offset - len; // offset of ">" or "@"
        // 		one.tot_length += len;              // count first line
        //
        // 		size_t i = 0;
        // 		if (des.data[i] == '>' || des.data[i] == '@' || des.data[i] == '+') i += 1;
        // 		if (des.data[i] == ' ' or des.data[i] == '\t') i += 1;
        // 		if (options.des_len and options.des_len < des.size) des.size = options.des_len;
        // 		while (i < des.size and !std::isspace(des.data[i])) i += 1;
        // 		des.data[i] = 0;                   // find first non-space letter
        // 		one.identifier = std::move(des._data);
        //
        // 		// R2
        // 		two.size = 0;
        // 		two.tot_length = 0;
        //
        // 		len = strlen(buffer2);
        // 		len2 = len;
        // 		while (len2 && buffer2[len2 - 1] != '\n') {
        // 			if ((res = gzgets(fin2, buffer2, MAX_LINE_SIZE)) == nullptr) break;
        // 			len2 = strlen(buffer2);
        // 			len += len2;
        // 		}
        // 		offset = gztell(fin2);
        // 		two.des_begin = offset - len;
        // 		two.tot_length += len;              // count first line
        // 		two.identifier = std::move(des._data);
        // 	}
        // 	else {
        // 		one.tot_length += strlen(buffer);  one += buffer;
        // 		two.tot_length += strlen(buffer2); two += buffer2;
        // 	}
        // }
        // one.identifier = nullptr;
        // two.identifier = nullptr;
        // delete[] buffer;
        // gzclose(fin);
        // delete[] buffer2;
        // gzclose(fin2);
    }

    // PE reads liwz, disable swap option
    void _Read_(const std::string_view sfile, const std::string_view sfile2) {
        const auto file = sfile.data();
        const auto file2 = sfile2.data();
        int f_len = strlen(file);
        int f_len2 = strlen(file2);
        if (strcmp(file + f_len - 3, ".gz") == 0) {
            if (strcmp(file2 + f_len2 - 3, ".gz"))
                throw std::invalid_argument("Both input files need to be in .gz format");
            // Readgz(file, file2);
            return;
        }

        Sequence one, two;
        Sequence des;
        auto fin = std::ifstream(file);
        if (fin.fail())
            throw std::invalid_argument("Failed to open the database file");
        auto fin2 = std::ifstream(file2);
        if (fin2.fail())
            throw std::invalid_argument("Failed to open the database file");
        Clear();

        // Todo: merge with ReadOne
        /*
        char* buffer = nullptr;
        char* buffer2 = nullptr;
        char* res = nullptr;
        char* res2 = nullptr;
        size_t option_l = options.min_length;
        buffer = new char[MAX_LINE_SIZE + 1];
        buffer2 = new char[MAX_LINE_SIZE + 1];

        while ((not fin.eof() && not fin2.eof()) || (one.size && two.size)) { // do not break when the last sequence is not handled
            buffer[0] = '>';
            res = fgets(buffer, MAX_LINE_SIZE, fin);
            buffer2[0] = '>';
            res2 = fgets(buffer2, MAX_LINE_SIZE, fin2);

            if ((res == nullptr) && (res2 != nullptr))
                throw std::invalid_argument("Paired input files have different number sequences");
            if ((res != nullptr) && (res2 == nullptr))
                throw std::invalid_argument("Paired input files have different number sequences");
            if ((one.size == 0) && (two.size > 0))
                throw std::invalid_argument("Paired input files have different number sequences");
            if ((one.size > 0) && (two.size == 0))
                throw std::invalid_argument("Paired input files have different number sequences");
            if ((res == nullptr) && (one.size == 0))
                break;

            if (buffer[0] == '+') { // fastq 3rd line
                // file 1
                int len = strlen(buffer);
                int len2 = len;
                while (len2 && buffer[len2 - 1] != '\n') { // read until the end of the line
                    if ((res = fgets(buffer, MAX_LINE_SIZE, fin)) == nullptr)
                        break;
                    len2 = strlen(buffer);
                    len += len2;
                }
                one.tot_length += len;

                // read next line quality score
                if ((res = fgets(buffer, MAX_LINE_SIZE, fin)) == nullptr)
                    throw std::invalid_argument("can not read quality score after");
                len = strlen(buffer);
                len2 = len;
                while (len2 && buffer[len2 - 1] != '\n') {
                    if ((res = fgets(buffer, MAX_LINE_SIZE, fin)) == nullptr)
                        break;
                    len2 = strlen(buffer);
                    len += len2;
                }
                one.tot_length += len;

                // file 2
                len = strlen(buffer2);
                len2 = len;
                while (len2 && buffer2[len2 - 1] != '\n') { // read until the end of the line
                    if ((res2 = fgets(buffer2, MAX_LINE_SIZE, fin2)) == nullptr)
                        break;
                    len2 = strlen(buffer2);
                    len += len2;
                }
                two.tot_length += len;

                // read next line quality score
                if ((res2 = fgets(buffer2, MAX_LINE_SIZE, fin2)) == nullptr)
                    throw std::invalid_argument("can not read quality score after");
                len = strlen(buffer2);
                len2 = len;
                while (len2 && buffer2[len2 - 1] != '\n') {
                    if ((res2 = fgets(buffer2, MAX_LINE_SIZE, fin2)) == nullptr)
                        break;
                    len2 = strlen(buffer2);
                    len += len2;
                }
                two.tot_length += len;
            }
            else if (buffer[0] == '>' || buffer[0] == '@' || (res == nullptr && one.size)) {
                if (one.size && two.size) { // write previous record
                    if (one.id_empty() || one.FormatAndValidate()) {
                        std::println("Warning: from file \"{}\",", file);
                        std::println("Discarding invalid sequence or sequence without identifier and description!\n");
                        if (one.identifier)
                            std::cout << one.get_identifier();
                        std::cout << one.get_data();
                        one.size = 0;
                        two.size = 0;
                    }
                    if (two.id_empty() || two.FormatAndValidate()) {
                        std::println("Warning: from file \"{}\",", file2);
                        std::println("Discarding invalid sequence or sequence without identifier and description!\n");
                        if (two.identifier)
                            std::cout << two.get_identifier();
                        std::cout << two.get_data();
                        one.size = 0;
                        two.size = 0;
                    }
                    one.index = sequences.size();
                    if ((one.size + two.size) > option_l) {
                        if (options.trim_len > 0)
                            one.trim(options.trim_len);
                        if (options.trim_len_R2 > 0)
                            two.trim(options.trim_len_R2);
                        sequences.emplace_back(std::move(one), std::move(two));
                    }
                }
                // R1
                one.size = 0;
                one.tot_length = 0;

                size_t len = strlen(buffer);
                size_t len2 = len;
                des.size = 0;
                des += buffer;
                while (len2 && buffer[len2 - 1] != '\n') {
                    if ((res = fgets(buffer, MAX_LINE_SIZE, fin)) == nullptr)
                        break;
                    des += buffer;
                    len2 = strlen(buffer);
                    len += len2;
                }
                size_t offset = ftell(fin);
                one.des_begin = offset - len; // offset of ">" or "@"
                one.tot_length += len; // count first line

                size_t i = 0;
                if (des.data[i] == '>' || des.data[i] == '@' || des.data[i] == '+')
                    i += 1;
                if (des.data[i] == ' ' or des.data[i] == '\t')
                    i += 1;
                if (options.des_len and options.des_len < des.size)
                    des.size = options.des_len;
                while (i < des.size and !std::isspace(des.data[i]))
                    i += 1;
                des.data[i] = 0; // find first non-space letter
                one.set_id(des);

                // R2
                two.size = 0;
                two.tot_length = 0;

                len = strlen(buffer2);
                len2 = len;
                while (len2 && buffer2[len2 - 1] != '\n') {
                    if ((res = fgets(buffer2, MAX_LINE_SIZE, fin2)) == nullptr)
                        break;
                    len2 = strlen(buffer2);
                    len += len2;
                }
                offset = ftell(fin2);
                two.des_begin = offset - len;
                two.tot_length += len; // count first line
                two.set_id(des);
            }
            else {
                one.tot_length += strlen(buffer);
                one += buffer;
                two.tot_length += strlen(buffer2);
                two += buffer2;
            }
        }
            */
    }

#if 0
	void Sort(int first, int last)
	{
		int lower = first + 1, upper = last;
		Sequence* pivot;
		Sequence* val;
		if (first >= last) return;
		val = sequences[first];
		sequences[first] = sequences[(first + last) / 2];
		sequences[(first + last) / 2] = val;
		pivot = sequences[first];

		while (lower <= upper) {
			while (lower <= last && sequences[lower]->stats < pivot->stats) lower++;
			while (pivot->stats < sequences[upper]->stats) upper--;
			if (lower < upper) {
				val = sequences[lower];
				sequences[lower] = sequences[upper];
				sequences[upper] = val;
				upper--;
			}
			lower++;
		}
		val = sequences[first];
		sequences[first] = sequences[upper];
		sequences[upper] = val;
		if (first < upper - 1) Sort(first, upper - 1);
		if (upper + 1 < last) Sort(upper + 1, last);
	}
#endif


    void SortDivide(bool sort = true) {
        total_letter = 0;
        total_desc = 0;
        max_len = 0;
        min_len = 0;
        for (auto& seq : sequences) {
            size_t len = seq.size;
            total_letter += len;
            if (len > max_len)
                max_len = len;
            if (len < min_len)
                min_len = len;
            seq.ConvertBases();
            total_desc += seq.get_identifier().size();
        }
        options.max_entries = max_len * MAX_TABLE_SEQ;

        if (max_len >= std::numeric_limits<INTs>::max() || max_len > MAX_SEQ)
            throw std::logic_error("Some sequences are too long, please rebuild the program with make parameter "
                                   "LONG_SEQ=1 MAX_SEQ=new-maximum-length (e.g. make MAX_SEQ=10000000)");

        cout << "longest and shortest : " << max_len << " and " << min_len << endl;
        cout << "Total letters: " << total_letter << endl;
        // END change all the NR_seq to iseq

        len_n50 = (max_len + min_len) / 2; // will be properly set, if sort is true;
        if (sort) {
            // **************************** Form NR_idx[], Sort them from Long to short
            size_t sum = 0;
            std::map<size_t, size_t> counter; // count for each size = max_len - i

            size_t sequences_number = sequences.size();

            for (const auto& seq : sequences)
                counter[seq.size]++;
            for (auto [sizer, number]: counter) {
                sum += sizer * number;
                if (sum >= (total_letter / 2)) {
                    len_n50 = sizer;
                    break;
                }
            }

            // template <>
            // struct std::less<Sequence> : std::ranges::less<Sequence, Sequence, bool> {
            //     [[nodiscard]] bool operator()(Sequence& _Left, Sequence& _Right) const {
            //         return _Left.get_data_size() > _Right.get_data_size();
            //     }
            //
            // };

            struct SeqBigger {
                bool operator()(const Sequence& a, const Sequence& b) const {
                    return a.get_data_size() > b.get_data_size();
                }
            };

            std::ranges::sort(sequences, SeqBigger{});

            options.max_entries = 0;
            for (size_t i = 0; i < sequences.size(); i++) {
                if (i < MAX_TABLE_SEQ)
                    options.max_entries += sequences[i].size;
            }
            cout << "Sequences have been sorted" << endl;
        }
    }


    void DivideSave(const std::string_view db, const std::string_view newdb, const int _n) {
        if (_n == 0 or sequences.empty())
            return;

        size_t max_seg = total_letter / _n + sequences[0].size;
        if (max_seg >= MAX_BIN_SWAP)
            max_seg = MAX_BIN_SWAP;

        auto fin = HitFStream(db);
        size_t seg_size = 0;
        int seg = 0;
        auto outfile = std::format("{}-{}", newdb, seg);
        auto fout = std::ofstream(outfile);
        auto n = sequences.size();
        for (int i = 0; i < n; i++) {
            Sequence& seq = sequences[i];
            fin.seek(seq.des_begin);

            if (seg_size + seq.size >= max_seg) {
                seg += 1;
                fout.close();
                auto outfile_new = std::format("{}-{}", newdb, seg);
                fout = std::ofstream(outfile_new);
                seg_size = 0;
            }
            seg_size += seq.size;

            string sbuf;
            while (!fin.eof()) {
                fin >> sbuf;
                if (fin.fail())
                    throw std::invalid_argument("Can not swap in sequence");
                fout << sbuf;
                if (fout.fail())
                    throw std::invalid_argument("Can not swap in sequence");
            }
        }
    }


    void WriteClusters(const std::string_view db, const std::string_view newdb, const std::string_view db_pe = {}, const std::string_view newdb_pe = {}) {
        std::ranges::sort(sequences, {}, [](const Sequence& s) { return s.index; });
        std::ranges::sort(sequences, {}, [](const Sequence& s) { return s.cluster_id; });

        auto fin = HitFStream(db);
        std::ofstream fout(newdb.data());

        for (const Sequence& seq : sequences) {
            // R1
            fin.seek(seq.des_begin);
            copy_sequence_chunked(fin, fout, seq.tot_length);
        }

        // R2
        if (newdb_pe.empty()) {
            if (db_pe.empty())
                return;
            else
                throw std::invalid_argument("PE output file is missing");
        }
        auto fin_pe = HitFStream(db_pe);
        std::ofstream fout_pe(newdb_pe.data());

        for (const Sequence& seq : sequences) {
            fin_pe.seek(seq.des_begin2);
            copy_sequence_chunked(fin_pe, fout_pe, seq.tot_length);
        }
    }

    void WriteExtra1D() const {
        // WriteExtra2D(this->sequences, this->sequences);
    }
    //
    //	void WriteExtra2D(const SequenceDB& other) const {
    //		WriteExtra2D(this->sequences, other.sequences);
    //	}
    //
    //	typedef std::vector<Sequence*> seqvec;
    //	static void WriteExtra2D(const seqvec& data, const seqvec& index) {
    //		string db_clstr = options.output + ".clstr";
    //		string db_clstr_bak = options.output + ".bak.clstr";
    //
    //		if (options.backupFile) {
    //			std::ofstream foutbk{db_clstr_bak};
    //			for (const auto seq: index) {
    //				foutbk << *seq;
    //			}
    //			if (data != index) {
    //				for (const auto seq : data)
    //					foutbk << *seq;
    //			}
    //		}
    //
    //		seqvec scoper_indx{index};
    //		std::ranges::sort(scoper_indx, {}, [](const Sequence* s) { return s->index; });
    //		auto sorted_seq_indx = &scoper_indx;
    //
    //		seqvec scoper_data;
    //        if (data != index) {
    //            scoper_data = data;
    //			std::ranges::sort(scoper_data, {}, [](const Sequence* s) { return s->index; });
    //			auto sorted_seq_data = &scoper_data;
    //		}
    //        else {
    //			auto sorted_seq_data = &scoper_indx;
    //	    }
    //
    //		cout << "writing clustering information" << endl;
    //		std::map<size_t, const Sequence*> clusters;
    //		for (auto s: data) {
    //			clusters.insert({s->cluster_id, s});
    //		}
    //
    //		size_t i = 0;
    //		ofstream fout(db_clstr);
    //		for (const auto& c : sorted_seq_data) {
    //			std::println(fout, ">Cluster {}{}", i++, c.first);
    //			for (const Sequence* s : c.second)
    //				fout << *s <<endl;
    //		}
    //	}


    //	static void WriteExtra2D4(const SequenceDB& thisone, const SequenceDB& other)
    //	{
    //		string db_clstr = options.output + ".clstr";
    //		string db_clstr_bak = options.output + ".bak.clstr";
    //		int i, k, N = other.sequences.size();
    //		int N2 = thisone.sequences.size();
    //		vector<long long> sorting(N);
    //		for (i = 0; i < N; i++) sorting[i] = ((long long)other.sequences[i]->index << 32) | i;
    //		std::ranges::sort(sorting);
    //
    //		FILE* fout;
    //		array<char, MAX_LINE_SIZE + 1> buf;
    //		if (options.backupFile) {
    //			fout = fopen(db_clstr_bak.c_str(), "w+");
    //			for (i = 0; i < N; i++) {
    //				Sequence* seq = other.sequences[sorting[i] & 0xffffffff];
    //				seq.PrintInfo(seq.cluster_id, fout, buf.data());
    //			}
    //			for (i = 0; i < N2; i++) {
    //				Sequence* seq = thisone.sequences[i];
    //				if (seq.state & seq_state::IS_REDUNDANT) seq.PrintInfo(seq.cluster_id, fout, buf.data());
    //			}
    //			fclose(fout);
    //		}
    //
    //		cout << "writing clustering information" << endl;
    //		std::vector<std::vector<int> > clusters(N);
    //		for (i = 0; i < N2; i++) {
    //			int id = thisone.sequences[i]->cluster_id;
    //			if (thisone.sequences[i]->state & seq_state::IS_REDUNDANT) clusters[id].push_back(i);
    //		}
    //
    //		fout = fopen(db_clstr.c_str(), "w+");
    //		for (i = 0; i < N; i++) {
    //			Sequence* seq = other.sequences[i];
    //			fprintf(fout, ">Cluster %i\n", i);
    //			seq.PrintInfo(0, fout, buf.data());
    //			for (k = 0; k < (int)clusters[i].size(); k++)
    //				thisone.sequences[clusters[i][k]]->PrintInfo(k + 1, fout, buf.data());
    //		}
    //	}
};


export void make_comp_short_word_index(size_t NAA) {
    int i, j, k, icomp, k1;
    int c[4] = {3, 2, 1, 0};
    unsigned char short_word[32]; // short_word[12] is enough

    int NAA1 = NAAN_array[1];
    int NAAN = NAAN_array[NAA];

    for (i = 0; i < NAAN; i++) {
        // decompose i back to short_word
        for (k = i, j = 0; j < NAA; j++) {
            short_word[j] = static_cast<unsigned char>(k % NAA1);
            k = k / NAA1;
        }

        // calc_comp_aan_list
        icomp = 0;
        for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--)
            icomp += c[short_word[k1]] * NAAN_array[k];

        Comp_AAN_idx[i] = icomp;
    }
}

inline constexpr std::array packed_blob{
    0x00000007, 0x00000008, 0x00000009, 0x00000009, 0x0000010a, 0x0000010b, 0x0000010c, 0x0000020d, 0x0000020e, 0x00000410, 0x00000410, 0x00000511, 0x00000512,
    0x00000714, 0x00010715, 0x00010715, 0x00020817, 0x00020819, 0x00020a19, 0x00030a1a, 0x00040d1c, 0x00050d1e, 0x00050e1e, 0x01060f21, 0x02071122, 0x02071123,
    0x02091425, 0x040a1425, 0x040b1628, 0x050c1829, 0x050c192a, 0x06101b2b, 0x08101b2d, 0x09111d2f, 0x0a121f2f, 0x0a142032, 0x0c142033, 0x0e162436, 0x0f182537,
    0x111a293a, 0x121d293b, 0x141e2d3c, 0x1823303e, 0x1a243040, 0x1b263341, 0x1f2b3644, 0x232b3746, 0x24303c47, 0x24323d49, 0x28323d4b, 0x2d36414b, 0x343c464f,
    0x353e4751, 0x39424b54, 0x39424c55, 0x40474e55, 0x464b5259, 0x4d51565c, 0x52565a5e, 0x53575b5f, 0x5b5d5f61, 0x00000109, 0x0000020a, 0x0000020b, 0x0000030c,
    0x0000030c, 0x0000040e, 0x0000040e, 0x00010510, 0x00010611, 0x00020713, 0x00020813, 0x00020814, 0x00020915, 0x00040a17, 0x01040b18, 0x01040b18, 0x01050d1a,
    0x02050d1b, 0x02060f1d, 0x02070f1e, 0x0308101f, 0x04081220, 0x04091221, 0x050b1424, 0x060c1625, 0x060c1626, 0x080e1828, 0x080f1929, 0x0a101b2a, 0x0a121c2d,
    0x0b121d2d, 0x0e151f2f, 0x0e162030, 0x0e162132, 0x11182434, 0x11192434, 0x121b2736, 0x141d2938, 0x151f2a3a, 0x151f2e3c, 0x1b232e3c, 0x1c25323f, 0x1f263240,
    0x222b3542, 0x242d3643, 0x29323c46, 0x2b333c47, 0x2d363f4a, 0x3037404b, 0x363c444e, 0x373e4750, 0x383f4750, 0x40464c54, 0x454a5056, 0x494e5358, 0x4a4e5459,
    0x5054575b, 0x53565a5d, 0x56595c5f, 0x5b5d5f61, 0x5c5d5f61, 0x0000020b, 0x0000030c, 0x0000030c, 0x0001040d, 0x0001050e, 0x0001050f, 0x00010610, 0x00020712,
    0x00020712, 0x00030914, 0x01040914, 0x01040a15, 0x01040b17, 0x02050c18, 0x02050c19, 0x02060d1a, 0x03070e1c, 0x03070f1c, 0x0408101e, 0x0509111f, 0x050a1220,
    0x060b1423, 0x060b1423, 0x070d1626, 0x080e1727, 0x080f1827, 0x0a101a2a, 0x0a111b2a, 0x0c131d2c, 0x0d141e2e, 0x0d151f2f, 0x10172130, 0x12192232, 0x121a2433,
    0x131c2635, 0x141d2635, 0x171e2938, 0x18212b39, 0x1a222d3b, 0x1c25303d, 0x1e25303e, 0x212a3440, 0x232b3541, 0x262f3844, 0x282f3844, 0x2c353d47, 0x2d353e49,
    0x323a424b, 0x333a424c, 0x393f474f, 0x3c424851, 0x3e444b53, 0x464a5055, 0x4a4e5258, 0x55575a5c, 0x56585a5c, 0x57595b5d, 0x57595c5e, 0x595b5d60, 0x5d5e6061,
    0x5e5f6162, 0x0001040d, 0x0001050d, 0x0001050e, 0x0002060f, 0x00020610, 0x00020711, 0x01030812, 0x01040914, 0x01040914, 0x02050b16, 0x02050b16, 0x02060c18,
    0x03060d19, 0x03070e1a, 0x04080e1b, 0x04080f1c, 0x0509111e, 0x0509111e, 0x060b1320, 0x070c1422, 0x080c1422, 0x090e1625, 0x090e1725, 0x0a101927, 0x0b111a29,
    0x0c121b29, 0x0d141c2b, 0x0e151e2d, 0x0f161f2e, 0x11182130, 0x11182230, 0x131a2432, 0x141b2533, 0x151d2735, 0x171f2937, 0x171f2937, 0x1a222c3a, 0x1c242e3b,
    0x1d252f3c, 0x2229323e, 0x222a333f, 0x262d3742, 0x272e3743, 0x2c333c46, 0x2c333c46, 0x31384049, 0x3239404a, 0x393f454d, 0x3a40464e, 0x44474c52, 0x44484d53,
    0x4b4f5155, 0x5657595a, 0x58595a5c, 0x5a5b5c5d, 0x5b5c5d5e, 0x5c5e5e5f, 0x5d5e5f60, 0x5e5f5f60, 0x5e5f6062, 0x5f606162, 0x0102060f, 0x01030710, 0x01030811,
    0x02040912, 0x02040913, 0x02050a14, 0x03050a15, 0x03060c16, 0x03060c17, 0x04080e19, 0x04080e19, 0x05080f1a, 0x0509101b, 0x060a111d, 0x060b121e, 0x070b121f,
    0x080c1420, 0x080d1421, 0x0a0e1623, 0x0a0f1725, 0x0b101825, 0x0c121a27, 0x0d121a28, 0x0e141c2a, 0x10161e2b, 0x10161f2c, 0x1117202d, 0x1219212f, 0x131a2330,
    0x151b2432, 0x161d2533, 0x181e2734, 0x19202935, 0x1a212a37, 0x1d232c39, 0x1d242d39, 0x2027303c, 0x2229323d, 0x242b333e, 0x282e3641, 0x282e3641, 0x2e343b44,
    0x2e343c45, 0x353b4149, 0x363c4249, 0x3f43494e, 0x44474b4f, 0x4e505255, 0x4f515355, 0x53555657, 0x55565759, 0x5658595a, 0x58595a5b, 0x5a5a5b5c, 0x5b5c5c5d,
    0x5c5d5e5e, 0x5e5e5f5f, 0x5f5f6060, 0x5f606161, 0x60606162, 0x61626263};


consteval static std::array<std::array<std::array<uint8_t, 4>, 61>, 5> make_naa_stat() {
    std::array<std::array<std::array<uint8_t, 4>, 61>, 5> naa_stat{};

    size_t n = 0; // index in packed_blob
    for (auto& coverlevel : naa_stat) {
        for (auto& N2_5 : coverlevel) {
            uint32_t packed = packed_blob[n++];
            N2_5[0] = (packed >> 0) & 0xFF;
            N2_5[1] = (packed >> 8) & 0xFF;
            N2_5[2] = (packed >> 16) & 0xFF;
            N2_5[3] = (packed >> 24) & 0xFF;
        }
    }

    return naa_stat;
}

inline constexpr auto naa_stat = make_naa_stat();

constexpr uint8_t get_naa_stat(size_t i, size_t j, size_t k) { return naa_stat[i][j][k]; }
