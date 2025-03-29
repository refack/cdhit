export module Options;

import std;

#include "cdhit-common.h"
#include "ScoreMatrix.h"

#include <cstddef>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <charconv>
#include <unordered_map>

#include <omp.h>
#include <stdexcept>
#include <iostream>

export extern Options options;
export const char *temp_dir = "";
export ScoreMatrix mat;

// idx for      A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
int na2idx[] = {0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4};

// idx for      A  B  C  D  E  F   G  H  I  J   K   L   M   N  O   P   Q  R  S   T   U   V   W   X   Y   Z
int aa2idx[] = {0, 2, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 6};

// so  aa2idx[ X  - 'A'] => idx_of_X, eg aa2idx['A' - 'A'] => 0,
// and aa2idx['M' - 'A'] => 12

export float try_parse_float(string_view value, string_view name)
{
    float fval;
    auto result = std::from_chars(value.data(), value.data() + value.size(), fval);
    if (result.ec != std::errc())
        throw std::invalid_argument(string{"Invalid float value: "} + value.data() + "for " + name.data());
    return fval;
}

export size_t try_parse_ull(string_view value, string_view name)
{
    size_t llval;
    auto result = std::from_chars(value.data(), value.data() + value.size(), llval);
    if (result.ec != std::errc())
        throw std::invalid_argument(string{"Invalid float value: "} + value.data() + "for " + name.data());
    return llval;
}

export size_t get_threads(int intval)
{
#ifndef NO_OPENMP
    auto cpus = omp_get_num_procs();
    auto threads = intval;
    if (threads > cpus)
    {
        threads = cpus;
        printf("Warning: total number of CPUs in the system is %i asked for %i\n", cpus, intval);
    }
    else if (threads < 0)
    {
        threads += cpus;
        if (threads < 0)
            threads = 0;
    }
    if (threads == 0)
    {
        threads = cpus;
        printf("total number of CPUs in the system is %i asked for %i\n", cpus, intval);
    }
    if (threads != intval)
        printf("Actual number of CPUs to be used: %i\n\n", threads);
    if (threads <= 0)
    {
        throw invalid_argument("Invalid number of threads: " + std::to_string(intval));
    }
    return (size_t)threads;
#else
    printf("Option -T is ignored: multi-threading with OpenMP is NOT enabled!\n");
    return 0;
#endif
}

export struct Options
{
    size_t NAA = 5;
    size_t NAAN;
    size_t NAA_top_limit = 5;

    size_t max_memory = 800000000;            // -M: 400,000,000 in bytes
    size_t min_length = 10;                  // -l: 10 bases
    bool cluster_best = false;               // -g: 0, the first; 1, the best
    bool global_identity = true;             // -G:
    bool store_disk = false;                 // -B:
    unsigned int band_width = 20;            // -b: 20
    double cluster_thd = 0.9;                // -c
    double distance_thd = 0.0;               // -D
    double diff_cutoff = 0.0;                // -s: 0.0
    double diff_cutoff2 = 1.0;               // -s2: 1.0
    size_t diff_cutoff_aa = 99999999;        // -S: 999999
    unsigned int diff_cutoff_aa2 = 0;        // -S2: 0
    size_t tolerance = 2;                    // -t: 2
    double long_coverage = 0.0;              // -aL:
    unsigned int long_control = 99999999;    // -AL:
    double short_coverage = 0.0;             // -aS:
    unsigned int short_control = 99999999;   // -AS:
    unsigned int min_control = 0;            // -A:
    double long_unmatch_per = 1.0;           // -uL
    double short_unmatch_per = 1.0;          // -uS
    size_t unmatch_len = 99999999;           // -U
    int max_indel = 1;                       // -D
    int print = 0;
    size_t des_len = 20;
    size_t frag_size = 0;
    int option_r = 1;
    size_t threads = 1;
    int PE_mode = 0;                         // -P
    int trim_len = 0;                        // -cx
    int trim_len_R2 = 0;                     // -cy
    size_t align_pos = 0;                    // -ap for alignment position

    size_t max_entries;
    size_t max_sequences = 1 << 20;
    size_t mem_limit = 100000000;

    bool has2D = false;
    bool isEST = false;
    bool is454 = false;
    bool useIdentity = false;
    bool useDistance = false;
    bool backupFile = false;

    string input;
    string input_pe;
    string input2;
    string input2_pe;
    string output;
    string output_pe;

    int sort_output = 0;  // -sc
    int sort_outputf = 0; // -sf

    Options() = default;
};

export void setaa_to_na()
{
    int i;
    for (i = 0; i < 26; i++)
        aa2idx[i] = na2idx[i];
}
