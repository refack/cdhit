#pragma once

constexpr unsigned int MAX_AA{23};
constexpr unsigned int MAX_SEQ{655360};

class ScoreMatrix
{
private:
public:
    int matrix[MAX_AA][MAX_AA];
    int gap, ext_gap;

    ScoreMatrix();
    void init();
    void set_gap(int gap1, int ext_gap1);
    void set_matrix(const int *mat1);
    void set_to_na();
    void set_match(int score);
    void set_mismatch(int score);
};
