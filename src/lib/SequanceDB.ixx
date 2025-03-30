export module SequenceDB;

import std;
import Options;
import WorkingBuffer;
import WordTable;
import ScoreMatrix;

import "cdhit-common.h";
import <cassert>;

using std::string;

std::vector<int> Comp_AAN_idx;

int upper_bound_length_rep(int len, double opt_s, int opt_S, double opt_aL, int opt_AL)
{
	int len_upper_bound = 99999999;
	double r1 = (opt_s > opt_aL) ? opt_s : opt_aL;
	int    a2 = (opt_S < opt_AL) ? opt_S : opt_AL;
	if (r1 > 0.0) len_upper_bound = (int)(((float)len) / r1);
	if ((len + a2) < len_upper_bound)  len_upper_bound = len + a2;

	return len_upper_bound;
} // END upper_bound_length_rep
int upper_bound_length_rep(int len)
{
	double opt_s = options.diff_cutoff;
	int    opt_S = options.diff_cutoff_aa;
	double opt_aL = options.long_coverage;
	int    opt_AL = options.long_control;
	return upper_bound_length_rep(len, opt_s, opt_S, opt_aL, opt_AL);
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
int diag_test_aapn(int NAA1, char iseq2[], int len1, int len2, WorkingBuffer& buffer, int& best_sum, int band_width, int& band_left, int& band_center, int& band_right, int required_aa1)
{
	int i, i1, j, k;
	int* pp;
	size_t nall = len1 + len2 - 1;
	std::vector<int>& taap = buffer.taap;
	std::vector<INTs>& aap_begin = buffer.aap_begin;
	std::vector<INTs>& aap_list = buffer.aap_list;
	std::vector<int>& diag_score = buffer.diag_score;
	std::vector<int>& diag_score2 = buffer.diag_score2;

	if (nall > MAX_DIAG) bomb_error("in diag_test_aapn, MAX_DIAG reached");
	for (pp = &diag_score[0], i = nall; i; i--, pp++) *pp = 0;
	for (pp = &diag_score2[0], i = nall; i; i--, pp++) *pp = 0;

	int c22, cpx;
	// INTs *bip;
	int len11 = len1 - 1;
	int len22 = len2 - 1;
	i1 = len11;
	for (i = 0; i < len22; i++, i1++) {
		c22 = iseq2[i] * NAA1 + iseq2[i + 1];
		cpx = 1 + (iseq2[i] != iseq2[i + 1]);
		if ((j = taap[c22]) == 0) continue;
		int m = aap_begin[c22];
		for (int k = 0; k < j; k++) {
			diag_score[i1 - aap_list[m + k]]++;
			diag_score2[i1 - aap_list[m + k]] += cpx;
		}
	}

	//find the best band range
	//  int band_b = required_aa1;
	int band_b = required_aa1 - 1 >= 0 ? required_aa1 - 1 : 0;  // on dec 21 2001
	int band_e = nall - band_b;

	int band_m = (band_b + band_width - 1 < band_e) ? band_b + band_width - 1 : band_e;
	int best_score = 0;
	int best_score2 = 0;
	// int max_diag(0);
	int max_diag2 = 0;
	int imax_diag = 0;
	for (i = band_b; i <= band_m; i++) {
		best_score += diag_score[i];
		best_score2 += diag_score2[i];
		if (diag_score2[i] > max_diag2) {
			max_diag2 = diag_score2[i];
			// max_diag = diag_score[i];
			imax_diag = i;
		}
	}
	int from = band_b;
	int end = band_m;
	int score = best_score;
	int score2 = best_score2;
	for (k = from, j = band_m + 1; j < band_e; j++, k++) {
		score -= diag_score[k];
		score += diag_score[j];
		score2 -= diag_score2[k];
		score2 += diag_score2[j];
		if (score2 > best_score2) {
			from = k + 1;
			end = j;
			best_score = score;
			best_score2 = score2;
			if (diag_score2[j] > max_diag2) {
				max_diag2 = diag_score2[j];
				// max_diag = diag_score[j];
				imax_diag = j;
			}
		}
	}
	int mlen = imax_diag;
	if (imax_diag > len1) mlen = nall - imax_diag;
	int emax = int((1.0 - options.cluster_thd) * mlen) + 1;
	for (j = from; j < imax_diag; j++) { // if aap pairs fail to open gap
		if ((imax_diag - j) > emax || diag_score[j] < 1) {
			best_score -= diag_score[j]; from++;
		}
		else break;
	}
	for (j = end; j > imax_diag; j--) { // if aap pairs fail to open gap
		if ((j - imax_diag) > emax || diag_score[j] < 1) {
			best_score -= diag_score[j]; end--;
		}
		else break;
	}

	//  delete [] diag_score;
	band_left = from - len1 + 1;
	band_right = end - len1 + 1;
	band_center = imax_diag - len1 + 1;
	best_sum = best_score;
	return OK_FUNC;
}
// END diag_test_aapn


int diag_test_aapn_est(int NAA1, char iseq2[], int len1, int len2, WorkingBuffer& buffer,
	int& best_sum, int band_width, int& band_left, int& band_center, int& band_right, int required_aa1)
{
	int* pp, * pp2;
	size_t nall = len1 + len2 - 1;
	int NAA2 = NAA1 * NAA1;
	int NAA3 = NAA2 * NAA1;
	const std::vector<int>& taap = buffer.taap;
	std::vector<INTs>& aap_begin = buffer.aap_begin;
	std::vector<INTs>& aap_list = buffer.aap_list;
	std::vector<int>& diag_score = buffer.diag_score;
	std::vector<int>& diag_score2 = buffer.diag_score2;

	if (nall > MAX_DIAG) bomb_error("in diag_test_aapn_est, MAX_DIAG reached");
	pp = &diag_score[0];
	pp2 = &diag_score2[0];
	for (auto i = nall; i; i--, pp++, pp2++) *pp = *pp2 = 0;

	INTs* bip;
	int c22, cpx;
	int len22 = len2 - 3;
	auto i1 = len1 - 1;
	for (auto i = 0; i < len22; i++, i1++, iseq2++) {
		unsigned char c0 = iseq2[0];
		unsigned char c1 = iseq2[1];
		unsigned char c2 = iseq2[2];
		unsigned char c3 = iseq2[3];
		if ((c0 >= 4) || (c1 >= 4) || (c2 >= 4) || (c3 >= 4)) continue; //skip N

		c22 = c0 * NAA3 + c1 * NAA2 + c2 * NAA1 + c3;
		auto j = taap[c22];
		if (j == 0)
			continue;
		cpx = 1 + int(c0 != c1) + int(c1 != c2) + (c2 != c3);
		bip = &aap_list[aap_begin[c22]];     //    bi = aap_begin[c22];
		for (; j; j--, bip++) {
			diag_score[i1 - *bip]++;
			diag_score2[i1 - *bip] += cpx;
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

	//find the best band range
	//  int band_b = required_aa1;
	size_t band_b = required_aa1 >= 1 ? required_aa1 - 1 : 0;  // on dec 21 2001
	assert(nall > band_b);
	size_t band_e = nall - band_b;

	if (options.is454) {
		band_b = len1 - band_width;
		band_e = len1 + band_width;
		if (band_b < 0) band_b = 0;
		if (band_e > nall) band_e = nall;
	}

	auto  band_m = (band_b + band_width - 1 < band_e) ? band_b + band_width - 1 : band_e;
	int best_score = 0;
	int best_score2 = 0;
	// int max_diag = 0;
	int max_diag2 = 0;
	int imax_diag = 0;
	for (auto i = band_b; i <= band_m; i++) {
		best_score += diag_score[i];
		best_score2 += diag_score2[i];
		if (diag_score2[i] > max_diag2) {
			max_diag2 = diag_score2[i];
			// max_diag = diag_score[i];
			imax_diag = i;
		}
	}
	int from = band_b;
	int end = band_m;
	int score = best_score;
	int score2 = best_score2;

	for (size_t k = from, j = band_m + 1; j < band_e; j++, k++)
	{
		score -= diag_score[k];
		score += diag_score[j];
		score2 -= diag_score2[k];
		score2 += diag_score2[j];
		if (score2 > best_score2) {
			from = k + 1;
			end = j;
			best_score = score;
			best_score2 = score2;
			if (diag_score2[j] > max_diag2) {
				max_diag2 = diag_score2[j];
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
	if (imax_diag > len1) mlen = nall - imax_diag;
	int emax = int((1.0 - options.cluster_thd) * mlen) + 1;
	for (auto j = from; j < imax_diag; j++) { // if aap pairs fail to open gap
		if ((imax_diag - j) > emax || diag_score[j] < 1) {
			best_score -= diag_score[j]; from++;
		}
		else break;
	}
	for (auto j = end; j > imax_diag; j--) { // if aap pairs fail to open gap
		if ((j - imax_diag) > emax || diag_score[j] < 1) {
			best_score -= diag_score[j]; end--;
		}
		else break;
	}

	band_left = from - len1 + 1;
	band_right = end - len1 + 1;
	band_center = imax_diag - len1 + 1;
	best_sum = best_score;
	if (options.is454) {
		if (band_left > 0) best_sum = 0;
		if (band_right < 0) best_sum = 0;
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

int local_band_align(char iseq1[], char iseq2[], int ilen1, int ilen2, ScoreMatrix& mat,
	int& best_score, int& iden_no, int& alnln, float& dist, unsigned int* alninfo,
	int band_left, int band_center, int band_right, WorkingBuffer& buffer)
{
	size_t best_score1;
	iden_no = 0;

	assert(ilen1 > 0 && ilen2 > 0);

	if ((band_right >= ilen2) ||
		(band_left <= -ilen1) ||
		(band_left > band_right))
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

	if (score_mat.size() <= len1)
	{
		VectorInt row(band_width1, 0);
		VectorInt64 row2(band_width1, 0);
		while (score_mat.size() <= len1)
		{
			score_mat.push_back(row2);
			back_mat.push_back(row);
		}
	}
	for (size_t i = 0; i <= len1; i++)
	{
		if (score_mat[i].Size() < band_width1)
			score_mat[i].Resize(band_width1);
		if (back_mat[i].Size() < band_width1)
			back_mat[i].Resize(band_width1);
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

	if (band_left < 0)
	{ // set score to left border of the matrix within band
		size_t tband = (band_right < 0) ? band_right : 0;
		// for (k=band_left; k<tband; k++)
		for (size_t k = band_left; k <= tband; k++)
		{ // fixed on 2006 11 14
			size_t i = -k;
			size_t j1 = k - band_left;
			// penalty for leading gap opening = penalty for gap extension
			// each of the left side query hunging residues give ext_gap (-1)
			score_mat[i][j1] = mat.ext_gap * i;
			back_mat[i][j1] = DP_BACK_TOP;
		}
		back_mat[-tband][tband - band_left] = DP_BACK_NONE;
	}

	if (band_right >= 0)
	{ // set score to top border of the matrix within band
		size_t tband = (band_left > 0) ? band_left : 0;
		for (size_t j = tband; j <= band_right; j++)
		{
			size_t j1 = j - band_left;
			score_mat[0][j1] = mat.ext_gap * j;
			back_mat[0][j1] = DP_BACK_LEFT;
		}
		back_mat[0][tband - band_left] = DP_BACK_NONE;
	}

	int gap_open[2] = { mat.gap, mat.ext_gap };
	int max_diag = band_center - band_left;
	int extra_score[4] = { 4, 3, 2, 1 };
	for (size_t i = 1; i <= len1; i++)
	{
		int J0 = 1 - band_left - i;
		int J1 = len2 - band_left - i;
		if (J0 < 0)
			J0 = 0;
		if (J1 >= band_width)
			J1 = band_width;
		for (int j1 = J0; j1 <= J1; j1++)
		{
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

			if (j1 > 0)
			{
				int gap = gap0;
				if (back_mat[i][j1 - 1] == DP_BACK_LEFT)
					gap = mat.ext_gap;
				size_t score = score_mat[i][j1 - 1] + gap;
				if (score > best_score1)
				{
					back = DP_BACK_LEFT;
					best_score1 = score;
				}
			}
			if (j1 + 1 < band_width)
			{
				int gap = gap0;
				if (back_mat[i - 1][j1 + 1] == DP_BACK_TOP)
					gap = mat.ext_gap;
				size_t score = score_mat[i - 1][j1 + 1] + gap;
				if (score > best_score1)
				{
					back = DP_BACK_TOP;
					best_score1 = score;
				}
			}
			score_mat[i][j1] = best_score1;
			back_mat[i][j1] = back;
			// printf( "%2i(%2i) ", best_score1, iden_no1 );
		}
		// printf( "\n" );
	}
	size_t i = 0;
	size_t j = 0;
	if (len2 - band_left < len1)
	{
		i = len2 - band_left;
		j = len2;
	}
	else if (len1 + band_right < len2)
	{
		i = len1;
		j = len1 + band_right;
	}
	else
	{
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
	int64_t score, smin = best_score1, smax = best_score1 - 1;
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
	for (IA = len1; IA > i; IA--)
	{
		AA[NN] = letters[iseq1[IA - 1]];
		BB[NN++] = '-';
	}
	for (IB = len2; IB > j; IB--)
	{
		AA[NN] = '-';
		BB[NN++] = letters[iseq2[IB - 1]];
	}
#endif

	int masked = 0;
	int indels = 0;
	int max_indels = 0;
	while (back != DP_BACK_NONE)
	{
		switch (back)
		{
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
			if (score < smin)
			{
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
			if (score < smin)
			{
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
			if (iseq1[i - 1] == iseq2[j - 1])
			{
				printf("%5i: %c %c %9i\n", pos, letters2[iseq1[i - 1]], letters2[iseq2[j - 1]], score_mat[i][j1]);
			}
			else
			{
				printf("%5i: %c %c %9i\n", pos, letters[iseq1[i - 1]], letters[iseq2[j - 1]], score_mat[i][j1]);
			}
#endif
#ifdef MAKEALIGN
			if (iseq1[i - 1] == iseq2[j - 1])
			{
				AA[NN] = letters2[iseq1[i - 1]];
				BB[NN++] = letters2[iseq2[j - 1]];
			}
			else
			{
				AA[NN] = letters[iseq1[i - 1]];
				BB[NN++] = letters[iseq2[j - 1]];
			}
#endif
			if (alninfo && options.global_identity)
			{
				if (i == 1 || j == 1)
				{
					gbegin1 = i - 1;
					gbegin2 = j - 1;
				}
				else if (i == len1 || j == len2)
				{
					gend1 = i - 1;
					gend2 = j - 1;
				}
			}
			score = score_mat[i][j1];
			i -= 1;
			j -= 1;
			match = iseq1[i] == iseq2[j];
			if (score > smax)
			{
				count = 0;
				smax = score;
				posmax = pos;
				end1 = i;
				end2 = j;
			}
			if (options.isEST && (iseq1[i] > 4 || iseq2[j] > 4))
			{
				masked += 1;
			}
			else
			{
				dlen += 1;
				dcount += (match != 0);
				count += match;
				count2 += match;
				count3 += match;
			}
			if (score < smin)
			{
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
		if (options.is454)
		{
			if (back == DP_BACK_LEFT_TOP)
			{
				if (indels > max_indels)
					max_indels = indels;
				indels = 0;
			}
			else
			{
				if (last == DP_BACK_LEFT_TOP)
				{
					indels = 1;
				}
				else if (indels)
				{
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
	if (alninfo)
	{
		alninfo[0] = begin1;
		alninfo[1] = end1;
		alninfo[2] = begin2;
		alninfo[3] = end2;
		alninfo[4] = masked;
		if (options.global_identity)
		{
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
	while (i--)
	{
		AA[NN] = letters[iseq1[i - 1]];
		BB[NN++] = '-';
	}
	while (j--)
	{
		AA[NN] = '-';
		BB[NN++] = letters[iseq2[j - 1]];
	}
	AA[NN] = '\0';
	BB[NN] = '\0';
	for (i = 0; i < NN / 2; i++)
	{
		char aa = AA[i], bb = BB[i];
		AA[i] = AA[NN - i - 1];
		BB[i] = BB[NN - i - 1];
		AA[NN - i - 1] = aa;
		BB[NN - i - 1] = bb;
	}
	static int fcount = 0;
	fcount += 1;
	FILE* fout = fopen("alignments.txt", "a");
	if (fout == NULL)
	{
		if (fcount <= 1)
			printf("alignment files open failed\n");
		return OK_FUNC;
	}
	fprintf(fout, "\n\n######################################################\n");
	fprintf(fout, "# length X = %i\n", len2);
	fprintf(fout, "# length Y = %i\n", len1);
	fprintf(fout, "# best align X: %i-%i\n", begin2 + 1, end2 + 1);
	fprintf(fout, "# best align Y: %i-%i\n", begin1 + 1, end1 + 1);
	if (alninfo)
	{
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
	while (IA < NN)
	{
		if (printaa)
		{
			fprintf(fout, "%c", BB[IB]);
			IB += 1;
			if (IB % 75 == 0 or IB == NN)
				printaa = false, fprintf(fout, "\nY ");
		}
		else
		{
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



export class SequenceDB
{
public:
	size_t NAAN;
	std::vector<Sequence *> sequences;
	std::vector<int> rep_seqs;

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

	void DoClustering(int T)
	{
		int i, k;
		int NAA = options.NAA;
		double aa1_cutoff = options.cluster_thd;
		double aas_cutoff = 1 - (1 - options.cluster_thd) * 4;
		double aan_cutoff = 1 - (1 - options.cluster_thd) * options.NAA;
		int seq_no = sequences.size();
		int frag_no = seq_no;
		int frag_size = options.frag_size;
		// int len, len_bound;
		// int flag;
		std::valarray<size_t>  letters(T);

		//printf( "%li\n".mem_limit );

		if (frag_size) {
			frag_no = 0;
			for (i = 0; i < seq_no; i++) frag_no += (sequences[i]->size - NAA) / frag_size + 1;
		}

		if (not options.isEST)
			cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd,	options.tolerance, naa_stat_start_percent, naa_stat, NAA);

		std::vector<WorkingParam> params(T);
		std::vector<WorkingBuffer> buffers(T);
		for (i = 0; i < T; i++) {
			params[i] = WorkingParam(aa1_cutoff, aas_cutoff, aan_cutoff);
			buffers[i].Set(frag_no, max_len);
		}

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

		omp_set_num_threads(T);
		for (i = 0; i < N; ) {
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
				Sequence* seq = sequences[m];
				if (!(seq->state & IS_REDUNDANT)) {
					if (options.store_disk) seq->SwapIn();
					//items += seq->size;
					items += (size_t)(seq->size * redundancy);
					sum += 1;
				}
				m++;
			}
			if ((m > i + 1E4) && (m > i + (N - i) / (2 + T))) m = i + (N - i) / (2 + T);
			if (m == i || m >= N) {
				m = N;
				if (m > i + 1E3) m = i + (N - i) / (2 + T);
			}
			//printf( "m = %i  %i,  %i\n", i, m, m-i );
			printf("\r# comparing sequences from  %9i  to  %9i\n", i, m);
			if (last_table.size) {
				int print = (m - i) / 20 + 1;
#pragma omp parallel for schedule( dynamic, 1 )
				for (int j = i; j < m; j++) {
					Sequence* seq = sequences[j];
					if (seq->state & IS_REDUNDANT) continue;
					int tid = omp_get_thread_num();
					CheckOne(seq, last_table, params[tid], buffers[tid]);
					if (options.store_disk && (seq->state & IS_REDUNDANT)) seq->SwapOut();
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
#pragma omp parallel for schedule( dynamic, 1 )
				for (int j = m - 1; j < N; j++) {
#pragma omp flush (stop)
					if (!stop) {
						if (j + 1 == N) may_stop = 1;
						if (j == (m0 - 1)) { // use m0 to avoid other iterations satisfying the condition:
							int tid = omp_get_thread_num();
							for (int ks = i; ks < m; ks++) {
								Sequence* seq = sequences[ks];
								i = ks + 1;
								if (seq->state & IS_REDUNDANT) continue;
								ClusterOne(seq, ks, word_table, params[tid], buffers[tid]);
								if (options.store_disk && (seq->state & IS_REDUNDANT)) seq->SwapOut();
								if (may_stop and word_table.sequences.size() >= 100) break;
								if (word_table.size >= max_items) break;
								int tmax = max_seqs - (frag_size ? seq->size / frag_size + 1 : 0);
								if (word_table.sequences.size() >= tmax) break;
							}
							self_stop = 1;
						}
						else {
							Sequence* seq = sequences[j];
							if (seq->state & IS_REDUNDANT) continue;
							if (options.store_disk) {
#pragma omp critical
								seq->SwapIn();
							}
							int tid = omp_get_thread_num();
							CheckOne(seq, last_table, params[tid], buffers[tid]);
							if (options.store_disk && (seq->state & IS_REDUNDANT)) seq->SwapOut();
							if (min > params[tid].len_upper_bound) {
								may_stop = 1;
								stop = true;
#pragma omp flush (stop)
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
				//printf( "comparing the first or last or very small group ...\n" ); fflush( stdout );
				for (k = i; k < m; ) {
					int kk, mm = k, sum = 0;
					while (mm < m && sum < 1E5) {
						if (!(sequences[mm]->state & IS_REDUNDANT)) sum += sequences[mm]->size;
						mm += 1;
					}
					if (mm < k + 1000) mm = k + 1000;
					if (mm > m) mm = m;
#pragma omp parallel for schedule( dynamic, 1 )
					for (kk = k; kk < mm; kk++) {
						Sequence* seq = sequences[kk];
						if (seq->state & IS_REDUNDANT) continue;
						int tid = omp_get_thread_num();
						CheckOne(seq, word_table, params[tid], buffers[tid]);
						if (options.store_disk && (seq->state & IS_REDUNDANT)) seq->SwapOut();
					}
					bool bk = false;
					for (int ks = k; ks < mm; ks++) {
						Sequence* seq = sequences[ks];
						i = k = ks + 1;
						if (seq->state & IS_REDUNDANT) continue;
						ClusterOne(seq, ks, word_table, params[0], buffers[0]);
						bk = true;
						if (options.store_disk && (seq->state & IS_REDUNDANT)) seq->SwapOut();
						if (word_table.size >= max_items) break;
						int tmax = max_seqs - (frag_size ? seq->size / frag_size + 1 : 0);
						if (word_table.sequences.size() >= tmax) break;
						bk = false;
					}
					if (bk) break;
				}
			}
			else if (i < m) {
				remaining = remaining / 2 + (m - i);
				printf("\r---------- %6i remaining sequences to the next cycle\n", m - i);
			}
			printf("---------- new table with %8i representatives\n", word_table.sequences.size());
			if ((last_table.size + word_table.size) > tabsize)
				tabsize = last_table.size + word_table.size;
			last_table.Clear();
			last_table.sequences.swap(word_table.sequences);
			last_table.indexCounts.swap(word_table.indexCounts);
			last_table.size = word_table.size;
			word_table.size = 0;
		}
		printf("\n%9i  finished  %9i  clusters\n", sequences.size(), rep_seqs.size());
		mem = (mem_need + tabsize * sizeof(IndexCount)) / MEGA_MiBi;
		printf("\nApproximated maximum memory consumption: %zuM\n", mem);
		last_table.Clear();
		word_table.Clear();
	}

	int CheckOne(Sequence* seq, WordTable& table, WorkingParam& param, WorkingBuffer& buf)
	{
		int len = seq->size;
		param.len_upper_bound = upper_bound_length_rep(len);
		if (options.isEST) return CheckOneEST(seq, table, param, buf);
		return CheckOneAA(seq, table, param, buf);
	}
	int CheckOneAA(Sequence* seq, WordTable& table, WorkingParam& param, WorkingBuffer& buf)
	{
		NVector<IndexCount>& lookCounts = buf.lookCounts;
		NVector<uint32_t>& indexMapping = buf.indexMapping;
		std::vector<INTs>& word_encodes_no = buf.word_encodes_no;
		// std::vector<INTs> & aap_list = buf.aap_list;
		// std::vector<INTs> & aap_begin = buf.aap_begin;
		std::vector<int>& word_encodes = buf.word_encodes;
		// std::vector<int>  & taap = buf.taap;
		double aa1_cutoff = param.aa1_cutoff;
		double aa2_cutoff = param.aas_cutoff;
		double aan_cutoff = param.aan_cutoff;

		char* seqi = seq->data;
		int k, j1;
		unsigned int len = seq->size;
		int flag = 0;
		int frag_size = options.frag_size;
		int& aln_cover_flag = param.aln_cover_flag;
		int& required_aa1 = param.required_aa1;
		int& required_aa2 = param.required_aas;
		int& required_aan = param.required_aan;
		unsigned int& min_aln_lenS = param.min_aln_lenS;
		unsigned int& min_aln_lenL = param.min_aln_lenL;

		int NAA = options.NAA;
		int S = table.sequences.size();
		int len_eff = len;

		if (S) {
			unsigned int min = table.sequences[S - 1]->size;
			if (min < len) {
				if (len * options.diff_cutoff2 > min) min = (int)(len * options.diff_cutoff2);
				if ((len - options.diff_cutoff_aa2) > min) min = len - options.diff_cutoff_aa2;
				len_eff = min;
			}
		}

		//liwz 2016 01, seq is too short for the shortest (longer) seq in word_table to satisfy -aL option
		//longer seqeunce * -aL -band_width
		if (S) {
			unsigned int min = table.sequences[S - 1]->size;
			unsigned int min_red = min * options.long_coverage - options.band_width;
			if (len < min_red) return 0; // return flag=0
		}

		param.ControlShortCoverage(len_eff);
		param.ComputeRequiredBases(options.NAA, 2);

		buf.EncodeWords(seq, options.NAA, false);

		// if minimal alignment length > len, return
		// I can not return earlier, because I need to calc the word_encodes etc
		if (options.min_control > len) return 0; // return flag=0

		// lookup_aan
		int aan_no = len - options.NAA + 1;
		// int M = frag_size ? table.frag_count : S;
		table.CountWords(aan_no, word_encodes, word_encodes_no, lookCounts, indexMapping, false, required_aan);

		// contained_in_old_lib()
		unsigned int len_upper_bound = param.len_upper_bound;
		unsigned int len_lower_bound = param.len_lower_bound;
		int band_left, band_right, best_score, band_width1, best_sum, alnln, len_eff1;
		unsigned int len2;
		int tiden_no, band_center;
		float tiden_pc, distance = 0;
		unsigned int talign_info[5];
		int sum;
		// INTs *lookptr;
		char* seqj;
		int frg2 = frag_size ? (len - NAA + options.band_width) / frag_size + 1 + 1 : 0;
		int lens;
		int has_aa2 = 0;

		IndexCount* ic = lookCounts.items;
		ic = lookCounts.items;
		for (; ic->count; ic++) {
			if (!frag_size) {
				indexMapping[ic->index] = 0;
				if (ic->count < required_aan) continue;
			}

			Sequence* rep = table.sequences[ic->index];
			len2 = rep->size;
			if (len2 > len_upper_bound) continue;
			if (options.has2D && len2 < len_lower_bound) continue;
			if (frag_size) {
				uint32_t* ims = &indexMapping[ic->index];
				int count = ic->count;
				k = (len2 - NAA) / frag_size + 1;
				sum = 0;
				for (j1 = 0; j1 < frg2; j1++) {
					uint32_t im = ims[j1];
					if (im) sum += lookCounts[im - 1].count;
				}
				count = sum;
				for (j1 = frg2; j1 < k; j1++) {
					uint32_t im1 = ims[j1];
					uint32_t im2 = ims[j1 - frg2];
					if (im1) sum += lookCounts[im1 - 1].count;
					if (im2) sum -= lookCounts[im2 - 1].count;
					if (sum > count) count = sum;
				}
				if (count < required_aan) continue;
			}

			param.ControlLongCoverage(len2);

			if (has_aa2 == 0) { // calculate AAP array
				buf.ComputeAAP(seqi, seq->size);
				has_aa2 = 1;
			}
			seqj = rep->data; //NR_seq[NR90_idx[j]];

			band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
			diag_test_aapn(NAA1, seqj, len, len2, buf, best_sum,
				band_width1, band_left, band_center, band_right, required_aa1);
			if (best_sum < required_aa2) continue;

			int rc = FAILED_FUNC;
			if (options.print || aln_cover_flag) //return overlap region
				rc = local_band_align(seqi, seqj, len, len2, mat,
					best_score, tiden_no, alnln, distance, talign_info,
					band_left, band_center, band_right, buf);
			else
				rc = local_band_align(seqi, seqj, len, len2, mat,
					best_score, tiden_no, alnln, distance, talign_info,
					band_left, band_center, band_right, buf);
			if (rc == FAILED_FUNC) continue;
			if (tiden_no < required_aa1) continue;
			lens = len;
			if (options.has2D && len > len2) lens = len2;
			len_eff1 = (options.global_identity == 0) ? alnln : (lens - talign_info[4]);
			tiden_pc = tiden_no / (float)len_eff1;
			if (options.useDistance) {
				if (distance > options.distance_thd) continue;
				if (distance >= seq->distance) continue; // existing distance
			}
			else {
				if (tiden_pc < options.cluster_thd) continue;
				if (tiden_pc <= seq->identity) continue; // existing iden_no
			}
			if (aln_cover_flag) {
				if (talign_info[3] - talign_info[2] + 1 < min_aln_lenL) continue;
				if (talign_info[1] - talign_info[0] + 1 < min_aln_lenS) continue;
			}
			if (options.has2D) seq->state |= IS_REDUNDANT;
			flag = 1; seq->identity = tiden_pc; seq->cluster_id = rep->cluster_id;
			seq->distance = distance;
			seq->coverage[0] = talign_info[0] + 1;
			seq->coverage[1] = talign_info[1] + 1;
			seq->coverage[2] = talign_info[2] + 1;
			seq->coverage[3] = talign_info[3] + 1;
			if (not options.cluster_best) break;
			update_aax_cutoff(aa1_cutoff, aa2_cutoff, aan_cutoff,
				options.tolerance, naa_stat_start_percent, naa_stat, NAA, tiden_pc);
			param.ComputeRequiredBases(options.NAA, 2);
		}
		if (frag_size) ic = lookCounts.items;
		while (ic->count) {
			indexMapping[ic->index] = 0;
			ic += 1;
		}
		lookCounts.size = 0;
		if (flag == 1) { // if similar to old one delete it
			if (!options.cluster_best) {
				seq->Clear();
				seq->state |= IS_REDUNDANT;
			}
		}
		return flag;

	}
	int CheckOneEST(Sequence* seq, WordTable& table, WorkingParam& param, WorkingBuffer& buf)
	{
		NVector<IndexCount>& lookCounts = buf.lookCounts;
		NVector<uint32_t>& indexMapping = buf.indexMapping;
		std::vector<INTs>& word_encodes_no = buf.word_encodes_no;
		// std::vector<INTs> & aap_list = buf.aap_list;
		// std::vector<INTs> & aap_begin = buf.aap_begin;
		std::vector<int>& word_encodes = buf.word_encodes;
		// std::vector<int>  & taap = buf.taap;
		std::vector<int>& aan_list_comp = buf.aan_list_comp;
		char* seqi_comp = &buf.seqi_comp[0];

		int& aln_cover_flag = param.aln_cover_flag;
		int& required_aa1 = param.required_aa1;
		int& required_aas = param.required_aas;
		int& required_aan = param.required_aan;
		unsigned int& min_aln_lenS = param.min_aln_lenS;
		unsigned int& min_aln_lenL = param.min_aln_lenL;

		char* seqi = seq->data;
		int j;
		unsigned int len = seq->size;
		int flag = 0;
		int S = table.sequences.size();
		int len_eff = len;
		if (S) {
			unsigned int min = table.sequences[S - 1]->size;
			if (min < len) {
				if (len * options.diff_cutoff2 > min) min = (int)(len * options.diff_cutoff2);
				if ((len - options.diff_cutoff_aa2) > min) min = len - options.diff_cutoff_aa2;
				len_eff = min;
			}
		}


		//liwz 2016 01, seq is too short for the shortest (longer) seq in word_table to satisfy -aL option
		//longer seqeunce * -aL -band_width
		if (S) {
			unsigned int min = table.sequences[S - 1]->size;
			unsigned int min_red = min * options.long_coverage - options.band_width;
			if (len < min_red) return 0; // return flag=0
		}


		param.ControlShortCoverage(len_eff);
		param.ComputeRequiredBases(options.NAA, 4);
		int skip = buf.EncodeWords(seq, options.NAA, true);
		required_aan -= skip;
		required_aas -= skip;
		required_aa1 -= skip;
		if (required_aan <= 0) required_aan = 1;
		if (required_aas <= 0) required_aas = 1;
		if (required_aa1 <= 0) required_aa1 = 1;

		// if minimal alignment length > len, return
		// I can not return earlier, because I need to calc the word_encodes etc
		if (options.min_control > len) return 0; // return flag=0

		int aan_no = len - options.NAA + 1;

		// contained_in_old_lib()
		unsigned int len_upper_bound = param.len_upper_bound;
		unsigned int len_lower_bound = param.len_lower_bound;
		int band_left, band_right, best_score, band_width1, best_sum, alnln, len_eff1;
		unsigned int len2;
		int tiden_no, band_center;
		float tiden_pc, distance = 0;
		unsigned int talign_info[5];
		int j0, comp, lens;
		char* seqj;

		for (comp = 0; comp < 2; comp++) {
			if (comp) {
				for (j0 = 0; j0 < aan_no; j0++) {
					j = word_encodes[j0];
					if (j < 0) aan_list_comp[j0] = j;
					else       aan_list_comp[j0] = Comp_AAN_idx[j];
				}
				make_comp_iseq(len, seqi_comp, seqi);
				seqi = seqi_comp;
			}
			int has_aas = 0;

			if (comp) {
				table.CountWords(aan_no, aan_list_comp, word_encodes_no, lookCounts, indexMapping, true, required_aan);
			}
			else {
				table.CountWords(aan_no, word_encodes, word_encodes_no, lookCounts, indexMapping, true, required_aan);
			}

			IndexCount* ic = lookCounts.items;
			ic = lookCounts.items;
			for (; ic->count; ic++) {
				indexMapping[ic->index] = 0;
				if (ic->count < required_aan) continue;
				Sequence* rep = table.sequences[ic->index];

				len2 = rep->size;
				if (len2 > len_upper_bound) continue;
				if (options.has2D && len2 < len_lower_bound) continue;

				seqj = rep->data;

				param.ControlLongCoverage(len2);

				if (has_aas == 0) { // calculate AAP array
					buf.ComputeAAP2(seqi, seq->size);
					has_aas = 1;
				}

				band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
				diag_test_aapn_est(NAA1, seqj, len, len2, buf, best_sum,
					band_width1, band_left, band_center, band_right, required_aa1);
				if (best_sum < required_aas) continue;
				//if( comp and flag and (not options.cluster_best) and j > rep->cluster_id ) goto Break;

				int rc = FAILED_FUNC;
				if (options.print || aln_cover_flag) { //return overlap region
					rc = local_band_align(seqi, seqj, len, len2, mat,
						best_score, tiden_no, alnln, distance, talign_info,
						band_left, band_center, band_right, buf);
					if (comp) {
						talign_info[0] = len - talign_info[0] - 1;
						talign_info[1] = len - talign_info[1] - 1;
					}
				}
				else {
					//printf( "%5i %5i %5i %5i\n", band_width1, band_right-band_left, band_left, band_right );
					rc = local_band_align(seqi, seqj, len, len2, mat,
						best_score, tiden_no, alnln, distance, talign_info,
						band_left, band_center, band_right, buf);
				}
				if (rc == FAILED_FUNC) continue;
				//printf( "%i  %i  %i\n", best_score, tiden_no, required_aa1 );
				if (tiden_no < required_aa1) continue;
				if (options.is454) {
					if (talign_info[2] != talign_info[0]) continue; // same start
					if (talign_info[0] > 1) continue; // one mismatch allowed at beginning
					if ((len - talign_info[1]) > 2) continue; // one mismatch allowed at end
				}

				lens = len;
				if (options.has2D && len > len2) lens = len2;
				len_eff1 = (options.global_identity == 0) ? alnln : (lens - talign_info[4]);
				tiden_pc = tiden_no / (float)len_eff1;
				//printf( "%i %f\n", tiden_no, tiden_pc );
				if (options.useDistance) {
					if (distance > options.distance_thd) continue;
					if (options.cluster_best and distance >= seq->distance) continue; // existing distance
				}
				else {
					if (tiden_pc < options.cluster_thd) continue;
					if (options.cluster_best and tiden_pc < seq->identity) continue; // existing iden_no
				}
				if (aln_cover_flag) {
					if (talign_info[3] - talign_info[2] + 1 < min_aln_lenL) continue;
					if (comp) {
						if (talign_info[0] - talign_info[1] + 1 < min_aln_lenS) continue;
					}
					else {
						if (talign_info[1] - talign_info[0] + 1 < min_aln_lenS) continue;
					}
				}
				if (options.cluster_best and fabs(tiden_pc - seq->identity) < 1E-9 and rep->cluster_id >= seq->cluster_id) continue;
				if ((not options.cluster_best) and flag != 0 and rep->cluster_id >= seq->cluster_id) continue;
				flag = comp ? -1 : 1;
				seq->identity = tiden_pc;
				seq->distance = distance;
				seq->cluster_id = rep->cluster_id;
				seq->coverage[0] = talign_info[0] + 1;
				seq->coverage[1] = talign_info[1] + 1;
				seq->coverage[2] = talign_info[2] + 1;
				seq->coverage[3] = talign_info[3] + 1;
				if (not options.cluster_best) break;
			}
			while (ic->count) {
				indexMapping[ic->index] = 0;
				ic += 1;
			}
			lookCounts.size = 0;
			if (not options.option_r) break;
		}
		if ((flag == 1) || (flag == -1)) { // if similar to old one delete it
			if (!options.cluster_best) {
				seq->Clear();
				seq->state |= IS_REDUNDANT;
			}
			if (flag == -1)
				seq->state |= IS_MINUS_STRAND;
			else
				seq->state &= ~IS_MINUS_STRAND;
		}
		return flag;
	}
	void ComputeDistance()
	{
		unsigned int N = sequences.size();
		int best_score, best_sum;
		int band_width1, band_left, band_center, band_right;
		int tiden_no, alnln;
		unsigned int talign_info[5];
		float distance;
		WorkingBuffer buf(N, max_len);

		std::vector<NVector<float> > dists(N, NVector<float>(N));

		Sequence comseq(*sequences[0]);

		for (unsigned int i = 0; i < N; i++) {
			Sequence* seq = sequences[i];
			char* seqi = seq->data;
			unsigned int len = seq->size;
			buf.EncodeWords(seq, options.NAA, false);
			buf.ComputeAAP2(seqi, seq->size);
			dists[i][i] = 0.0;
			if ((i + 1) % 1000 == 0) printf("%9i\n", (i + 1));
			for (unsigned int j = 0; j < i; j++) {
				Sequence* rep = sequences[j];
				char* seqj = rep->data;
				unsigned int len2 = rep->size;
				band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
				diag_test_aapn_est(NAA1, seqj, len, len2, buf, best_sum,
					band_width1, band_left, band_center, band_right, 0);
				local_band_align(seqi, seqj, len, len2, mat,
					best_score, tiden_no, alnln, distance, talign_info,
					band_left, band_center, band_right, buf);
				dists[seq->index][rep->index] = dists[rep->index][seq->index] = distance;
			}
			if (not options.option_r) break;
			comseq.index = seq->index;
			comseq.size = len;
			for (unsigned int j = 0; j < len; j++) comseq.data[i] = seq->data[len - i - 1];
			seqi = comseq.data;
			buf.EncodeWords(&comseq, options.NAA, false);
			buf.ComputeAAP2(seqi, seq->size);
			for (unsigned int j = 0; j < i; j++) {
				Sequence* rep = sequences[j];
				char* seqj = rep->data;
				unsigned int len2 = rep->size;
				band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
				diag_test_aapn_est(NAA1, seqj, len, len2, buf, best_sum,
					band_width1, band_left, band_center, band_right, 0);
				local_band_align(seqi, seqj, len, len2, mat,
					best_score, tiden_no, alnln, distance, talign_info,
					band_left, band_center, band_right, buf);
				if (distance < dists[seq->index][rep->index])
					dists[seq->index][rep->index] = dists[rep->index][seq->index] = distance;
			}
		}
		std::string output = options.output + ".dist";
		FILE* fout = fopen(output.c_str(), "w+");
		fprintf(fout, "1");
		for (unsigned int i = 1; i < N; i++) fprintf(fout, "\t%i", i + 1);
		fprintf(fout, "\n");
		for (unsigned int i = 0; i < N; i++) {
			fprintf(fout, "%g", dists[i][0]);
			for (unsigned int j = 1; j < N; j++) fprintf(fout, "\t%g", dists[i][j]);
			fprintf(fout, "\n");
		}
		fclose(fout);
	}
	void DoClustering()
	{
		size_t NAA = options.NAA;
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
			for (size_t i = 0; i < seq_no; i++) frag_no += (sequences[i]->size - NAA) / frag_size + 1;
		}

		if (not options.isEST)
			cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd, options.tolerance, naa_stat_start_percent, naa_stat, NAA);

		WorkingParam param(aa1_cutoff, aas_cutoff, aan_cutoff);
		WorkingBuffer buffer(frag_no, max_len);

		WordTable word_table(options.NAA, NAAN);

		size_t mem_need = MinimalMemory(frag_no, buffer.total_bytes, 1);
		// TODO: size_t mem_limit = MemoryLimit( mem_need );
		size_t N = sequences.size();

		size_t total_letters = total_letter;
		size_t tabsize = 0;

		options.ComputeTableLimits(min_len, max_len, len_n50, mem_need);

		for (size_t i = 0; i < N; ) {
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
				Sequence* seq = sequences[m];
				if (!(seq->state & IS_REDUNDANT)) {
					if (options.store_disk) seq->SwapIn();
					items += (size_t)(seq->size * redundancy);
					sum += 1;
				}
				m++;
			}
			if (m > N) m = N;
			printf("\rcomparing sequences from  %9i  to  %9i\n", i, m);
			for (int ks = i; ks < m; ks++) { // clustering this block
				Sequence* seq = sequences[ks];
				i = ks + 1;
				if (seq->state & IS_REDUNDANT) continue;
				ClusterOne(seq, ks, word_table, param, buffer);
				total_letters -= seq->size;
				if (options.store_disk && (seq->state & IS_REDUNDANT)) seq->SwapOut();
				if (word_table.size >= max_items) break;
				int tmax = max_seqs - (frag_size ? seq->size / frag_size + 1 : 0);
				if (word_table.sequences.size() >= tmax) break;
			} // finishing word table from this block
			m = i;
			if (word_table.size == 0) continue;
			float p0 = 0;
			for (size_t j = m; j < N; j++) { // use this word table to screen rest sequences m->N
				Sequence* seq = sequences[j];
				if (seq->state & IS_REDUNDANT) continue;
				if (options.store_disk) seq->SwapIn();
				CheckOne(seq, word_table, param, buffer);
				total_letters -= seq->size;
				if (options.store_disk && (seq->state & IS_REDUNDANT)) seq->SwapOut();
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
			if (word_table.size > tabsize) tabsize = word_table.size;
			//if( i && i < m ) printf( "\r---------- %6i remaining sequences to the next cycle\n", m-i );
			word_table.Clear();
		}
		auto mem = (mem_need + tabsize * sizeof(IndexCount)) / MEGA_MiBi;
		std::cout << std::endl;
		std::cout << sequences.size() << "  finished  " << rep_seqs.size() << "  clusters";
		std::cout << std::endl;
		std::cout << "Approximated maximum memory consumption: " << mem << "M";
		std::cout << std::endl;
		word_table.Clear();

#if 0
		int zeros = 0;
		for (i = 0; i < word_table.indexCounts.size(); i++) zeros += word_table.indexCounts[i].Size() == 0;
		printf("%9i  empty entries out of  %9i\n", zeros, word_table.indexCounts.size());
#endif

	}

	void ClusterTo(SequenceDB& other)
	{
		int i, flag;
		int len, len_tmp, len_lower_bound, len_upper_bound;
		int NR2_red_no = 0;
		int aan_no = 0;
		char* seqi;
		int NAA = options.NAA;
		double aa1_cutoff = options.cluster_thd;
		double aas_cutoff = 1 - (1 - options.cluster_thd) * 4;
		double aan_cutoff = 1 - (1 - options.cluster_thd) * options.NAA;
		std::vector<int>  word_encodes(MAX_SEQ);
		std::vector<INTs> word_encodes_no(MAX_SEQ);

		if (not options.isEST) {
			cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd, options.tolerance, naa_stat_start_percent, naa_stat, NAA);
		}


		int N = other.sequences.size();
		int M = sequences.size();
		int T = options.threads;

		std::valarray<size_t>  counts(T);
		std::vector<WorkingParam> params(T);
		std::vector<WorkingBuffer> buffers(T);
		WorkingParam& param = params[0];
		WorkingBuffer& buffer = buffers[0];
		for (i = 0; i < T; i++) {
			params[i] = WorkingParam(aa1_cutoff, aas_cutoff, aan_cutoff);
			buffers[i].Set(N, max_len);
		}
		if (T > 1) omp_set_num_threads(T);

		size_t mem_need = MinimalMemory(N, buffer.total_bytes, T, other.total_letter + other.total_desc);
		// TODO: size_t mem_limit = MemoryLimit( mem_need );

		options.ComputeTableLimits(min_len, max_len, len_n50, mem_need);

		WordTable word_table(options.NAA, NAAN);

		size_t max_items = options.max_entries;
		size_t max_seqs = options.max_sequences;
		for (i = 0; i < N; ) {
			size_t items = 0;
			size_t sum = 0;
			int m = i;
			while (m < N && sum < max_seqs && items < max_items) {
				Sequence* seq = other.sequences[m];
				if (!(seq->state & IS_REDUNDANT)) {
					if (options.store_disk) seq->SwapIn();
					items += seq->size;
					sum += 1;
				}
				m++;
			}
			if (m > N) m = N;
			//printf( "m = %i  %i,  %i\n", i, m, m-i );
			for (int ks = i; ks < m; ks++) {
				Sequence* seq = other.sequences[ks];
				len = seq->size;
				seqi = seq->data;
				calc_ann_list(len, seqi, NAA, aan_no, word_encodes, word_encodes_no, options.isEST);
				word_table.AddWordCounts(aan_no, word_encodes, word_encodes_no, ks - i, options.isEST);
				word_table.sequences.push_back(seq);
				seq->cluster_id = ks;
				seq->state |= IS_REP;
				auto ks1 = ks + 1;
				if (ks1 % 1000 == 0) {
					std::cout << "." << std::flush;
					if (ks1 % 10000 == 0)
						std::cout << std::setw(9) << ks1 << "  finished\n" << std::flush;
				}
			}
			float p0 = 0;
			if (T > 1) {
				int JM = M;
				counts = 0;
#pragma omp parallel for schedule( dynamic, 1 )
				for (int j = 0; j < JM; j++) {
					Sequence* seq = sequences[j];
					if (seq->state & IS_REDUNDANT) continue;
					int len = seq->size;
					// char *seqi = seq->data;
					int len_upper_bound = upper_bound_length_rep(len);
					int len_lower_bound = len - options.diff_cutoff_aa2;

					int len_tmp = (int)(((double)len) * options.diff_cutoff2);
					if (len_tmp < len_lower_bound) len_lower_bound = len_tmp;

					int tid = omp_get_thread_num();
					params[tid].len_upper_bound = len_upper_bound;
					params[tid].len_lower_bound = len_lower_bound;

					if (word_table.sequences[word_table.sequences.size() - 1]->size > len_upper_bound) {
						JM = 0;
						continue;
					}

					int flag = other.CheckOne(seq, word_table, params[tid], buffers[tid]);
					if ((flag == 1) || (flag == -1)) { // if similar to old one delete it
						if (!options.cluster_best) {
							seq->Clear();
							seq->state |= IS_REDUNDANT;
							counts[tid]++;
						}
						if (flag == -1) seq->state |= IS_MINUS_STRAND; // for EST only
					}
					float p = (100.0 * j) / M;
					if (p > p0 + 1E-1) { // print only if the percentage changed
						std::cout << '\r' << std::fixed << std::setprecision(1) << std::setw(4) << p << '%' << std::flush;
						p0 = p;
					}
				}
				for (int j = 0; j < T; j++) NR2_red_no += counts[j];
			}
			else {
				for (int j = 0; j < M; j++) {
					Sequence* seq = sequences[j];
					if (seq->state & IS_REDUNDANT) continue;
					len = seq->size;
					seqi = seq->data;
					len_upper_bound = upper_bound_length_rep(len);
					len_lower_bound = len - options.diff_cutoff_aa2;

					len_tmp = (int)(((double)len) * options.diff_cutoff2);
					if (len_tmp < len_lower_bound) len_lower_bound = len_tmp;
					param.len_upper_bound = len_upper_bound;
					param.len_lower_bound = len_lower_bound;

					if (word_table.sequences[word_table.sequences.size() - 1]->size > len_upper_bound) {
						break;
					}

					flag = other.CheckOne(seq, word_table, param, buffer);
					if ((flag == 1) || (flag == -1)) { // if similar to old one delete it
						if (!options.cluster_best) {
							seq->Clear();
							seq->state |= IS_REDUNDANT;
							NR2_red_no++;
						}
						if (flag == -1) seq->state |= IS_MINUS_STRAND; // for EST only
					}
					float p = (100.0 * j) / M;
					if (p > p0 + 1E-1) { // print only if the percentage changed
						std::cout << '\r' << std::fixed << std::setprecision(1) << std::setw(4) << p << '%' << std::flush;
						p0 = p;
					}
				}
			}
			printf("\r..........%9i  compared  %9i  clusters\n", i, NR2_red_no);
			word_table.Clear();
			word_table.size = 0;
			i = m;
		}

		if (options.cluster_best) {//delete redundant sequences in options.cluster_best mode
			for (i = 0; i < (int)sequences.size(); i++) {
				Sequence* seq = sequences[i];
				if (seq->identity > 0) {
					seq->state |= IS_REDUNDANT;
					NR2_red_no++;
				}
			}
		}
		for (i = 0; i < (int)sequences.size(); i++) {
			Sequence* seq = sequences[i];
			if (seq->identity < 0) seq->identity *= -1;
			if (not(seq->state & IS_REDUNDANT)) rep_seqs.push_back(i);
		}

		std::cout << std::endl;
		std::cout << sequences.size() << " compared\t" << NR2_red_no << " clustered" << std::endl;
	}

	void ClusterOne(Sequence* seq, int id, WordTable& table, WorkingParam& param, WorkingBuffer& buffer)
	{
		if (seq->state & IS_REDUNDANT) return;
		int frag_size = options.frag_size;
		int NAA = options.NAA;
		int len = seq->size;
		int len_bound = upper_bound_length_rep(len);
		param.len_upper_bound = len_bound;
		int flag = CheckOne(seq, table, param, buffer);

		if (flag == 0) {
			if ((seq->identity > 0) && (options.cluster_best)) {
				// because of the -g option, this seq is similar to seqs in old SEGs
				seq->state |= IS_REDUNDANT;
				seq->Clear();
			}
			else {                  // else add to NR90 db
				int aan_no = len - NAA + 1;
				int size = rep_seqs.size();
				rep_seqs.push_back(id);
				seq->cluster_id = size;
				seq->identity = 0;
				seq->state |= IS_REP;
				if (frag_size) { /* not used for EST */
					int frg1 = (len - NAA) / frag_size + 1;
					table.AddWordCountsFrag(aan_no, buffer.word_encodes_backup,
						buffer.word_encodes_no, frg1, frag_size);
				}
				else {
					table.AddWordCounts(aan_no, buffer.word_encodes, buffer.word_encodes_no, table.sequences.size(), options.isEST);
				}
				table.sequences.push_back(seq);
				if (frag_size) {
					while (table.sequences.size() < table.frag_count)
						table.sequences.push_back(seq);
				}
			}
		}
		auto id1 = id + 1;
		if (id1 % 1000 == 0) {
			int size = rep_seqs.size();
			std::cout << "." << std::flush;
			if (id1 % 10000 == 0)
				std::cout << "\r.........." << std::setw(9) << id1 << "  finished  " << std::setw(9) << size << "clusters\n" << std::flush;

		}
	}


	size_t MinimalMemory(int frag_no, int bsize, int T, size_t extra = 0)
	{
		int N = sequences.size();
		int F = frag_no < MAX_TABLE_SEQ ? frag_no : MAX_TABLE_SEQ;
		size_t mem_need = 0;
		size_t mem;
		int table = T > 1 ? 2 : 1;

		printf("\nApproximated minimal memory consumption:\n");
		mem = N * sizeof(Sequence) + total_desc + N + extra;
		if (options.store_disk == false) mem += total_letter + N;
		printf("%-16s: %zuM\n", "Sequence", mem / MEGA_MiBi);
		mem_need += mem;

		mem = bsize;
		printf("%-16s: %i X %zuM = %zuM\n", "Buffer", T, mem / MEGA_MiBi, T * mem / MEGA_MiBi);
		mem_need += T * mem;

		mem = F * (sizeof(Sequence*) + sizeof(IndexCount)) + NAAN * sizeof(NVector<IndexCount>);
		printf("%-16s: %i X %zuM = %zuM\n", "Table", table, mem / MEGA_MiBi, table * mem / MEGA_MiBi);
		mem_need += table * mem;

		mem = sequences.capacity() * sizeof(Sequence*) + N * sizeof(int);
		mem += Comp_AAN_idx.size() * sizeof(int);
		printf("%-16s: %zuM\n", "Miscellaneous", mem / MEGA_MiBi);
		mem_need += mem;

		printf("%-16s: %zuM\n\n", "Total", mem_need / MEGA_MiBi);

		if (options.max_memory and options.max_memory < mem_need + 50 * table) {
			char msg[200];
			sprintf(msg, "not enough memory, please set -M option greater than %zu\n", 50 * table + mem_need / MEGA_MiBi);
			bomb_error(msg);
		}
		return mem_need;
	}

};
