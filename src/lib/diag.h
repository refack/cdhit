import WorkingBuffer;

int diag_test_aapn(int NAA1, char iseq2[], int len1, int len2,
    WorkingBuffer &buffer, int &best_sum,
    int band_width, int &band_left, int &band_center, int &band_right, int required_aa1);
int diag_test_aapn_est(int NAA1, char iseq2[], int len1, int len2,
        WorkingBuffer &buffer, int &best_sum,
        int band_width, int &band_left, int &band_center, int &band_right, int required_aa1);
int local_band_align(char query[], char ref[], int qlen, int rlen, ScoreMatrix &mat,
      int &best_score, int &iden_no, int &alnln, float &dist, unsigned int *alninfo,
      int band_left, int band_center, int band_right, WorkingBuffer &buffer);