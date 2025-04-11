export module ScoreMatrix;

import std;

import common;

constexpr std::size_t MAX_AA = 23;

constexpr std::array BLOSUM62_PROTEIN{
	4,                                                                  // A
   -1, 5,                                                               // R
   -2, 0, 6,                                                            // N
   -2,-2, 1, 6,                                                         // D
	0,-3,-3,-3, 9,                                                      // C
   -1, 1, 0, 0,-3, 5,                                                   // Q
   -1, 0, 0, 2,-4, 2, 5,                                                // E
	0,-2, 0,-1,-3,-2,-2, 6,                                             // G
   -2, 0, 1,-1,-3, 0, 0,-2, 8,                                          // H
   -1,-3,-3,-3,-1,-3,-3,-4,-3, 4,                                       // I
   -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,                                    // L
   -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,                                 // K
   -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5,                              // M
   -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,                           // F
   -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,                        // P
	1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4,                     // S
	0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,                  // T
   -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11,               // W
   -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,            // Y
	0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,         // V
   -2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4,      // B
   -1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,   // Z
	0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1 // X
  //A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
  //0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  2  6 20
};

constexpr size_t MAX_NA  = 6;

constexpr std::array BLOSUM62_DNA = {
	2,                  // A
   -2, 2,               // C
   -2,-2, 2,            // G
   -2,-2,-2, 2,         // T
   -2,-2,-2, 1, 2,      // U
   -2,-2,-2,-2,-2, 1,   // N
	0, 0, 0, 0, 0, 0, 1 // X
  //A  C  G  T  U  N  X
  //0  1  2  3  3  4  5
};

export template <bool is_dna>
class ScoreMatrix {
    static consteval std::size_t get_max() { return is_dna ? MAX_NA : MAX_AA; }
    static consteval std::size_t get_blosum_dim() { return is_dna ? BLOSUM62_DNA.size() : BLOSUM62_PROTEIN.size(); }
    using matrix_t = std::array<std::array<int, get_blosum_dim()>, get_blosum_dim()>;
    static consteval int get_static_gap() { return static_cast<int>(MAX_SEQ) * (is_dna ? -6 : -11); }
    static consteval matrix_t get_first() {
        matrix_t matrix{};
        unsigned k = 0;
        for (unsigned i = 0; i < MAX; i++)
            // matrix is lower triangular
            for (unsigned j = 0; j <= i; j++)
                matrix[j][i] = matrix[i][j] = MAX_SEQ * (is_dna ? BLOSUM62_DNA[k++] : BLOSUM62_PROTEIN[k++]);
        return matrix;
    }

    static constexpr std::size_t MAX = get_max();

    // constexpr int ext_gap = -MAX_SEQ;
    int gap = get_static_gap();
    int extra_gap = 0;

public:
    int get_gap() const { return gap; }
    int get_extra_gap() const { return extra_gap; }
    static matrix_t matrix;

    constexpr ScoreMatrix() {}
    constexpr ScoreMatrix(bool is_est) {}

    // Only for est
    void set_match(int score) {
        for (auto i = 0; i < 5; i++)
            matrix[i][i] = MAX_SEQ * score;
        // matrix[3][4] = matrix[4][3] = MAX_SEQ * score;
    }
    // Only for est

    void set_mismatch(int score) {
        for (unsigned i = 0; i < MAX_AA; i++)
            for (unsigned j = 0; j < i; j++)
                matrix[j][i] = matrix[i][j] = MAX_SEQ * score;
        matrix[3][4] = matrix[4][3] = MAX_SEQ;
    }

    void set_gaps(int _gap, int _extra_gap) {
        gap = MAX_SEQ * _gap;
        extra_gap = MAX_SEQ * _extra_gap;
    }
};

auto ScoreMatrix<true>::matrix = get_first();
auto ScoreMatrix<false>::matrix = get_first();

export ScoreMatrix<false> mat_aa{};
export ScoreMatrix<true> mat_na{};

// TODO: find usages and verify DNA/Protein usage
export auto& mat = mat_na;
