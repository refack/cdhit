#include "cdhit-common.h"
import Options;
#include "ScoreMatrix.h"
#include <assert.h>
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

int local_band_align(char iseq1[], char iseq2[], int ilen1, int ilen2, ScoreMatrix &mat,
                     int &best_score, int &iden_no, int &alnln, float &dist, unsigned int *alninfo,
                     int band_left, int band_center, int band_right, WorkingBuffer &buffer)
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
    MatrixInt64 &score_mat = buffer.score_mat;
    MatrixInt &back_mat = buffer.back_mat;

    // printf( "%i  %i\n", band_right, band_left );

    if (score_mat.size() <= len1)
    {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= len1)
        {
            score_mat.Append(row2);
            back_mat.Append(row);
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

    int gap_open[2] = {mat.gap, mat.ext_gap};
    int max_diag = band_center - band_left;
    int extra_score[4] = {4, 3, 2, 1};
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
                dcount += !match;
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
    FILE *fout = fopen("alignments.txt", "a");
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
fprintf( fout, "%i %s\n", seq1->index, AA );
fprintf( fout, "%i %s\n", seq2->index, BB );
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
