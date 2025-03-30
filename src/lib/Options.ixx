export module Options;

import std;

import "cdhit-common.h";
import ScoreMatrix;

import <cstddef>;
import <cstdlib>;
import <cstdio>;
import <charconv>;
import <unordered_map>;

import <omp.h>;

using std::string;
using std::string_view;

export string temp_dir = "";

// idx for  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
auto na2idx{0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4};

// idx for  A  B  C  D  E  F   G  H  I  J   K   L   M   N  O   P   Q  R  S   T   U   V   W   X   Y   Z
auto aa2idx{0, 2, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 6};

// so  aa2idx[ X  - 'A'] => idx_of_X, eg aa2idx['A' - 'A'] => 0,
// and aa2idx['M' - 'A'] => 12

float try_parse_float(std::string_view value, std::string_view name)
{
    float fval;
    auto result = std::from_chars(value.data(), value.data() + value.size(), fval);
    if (result.ec != std::errc())
        throw std::invalid_argument(string{"Invalid float value: "} + value.data() + "for " + name.data());
    return fval;
}

size_t try_parse_ull(std::string_view value, std::string_view name)
{
    size_t llval;
    auto result = std::from_chars(value.data(), value.data() + value.size(), llval);
    if (result.ec != std::errc())
        throw std::invalid_argument(string{"Invalid float value: "} + value.data() + "for " + name.data());
    return llval;
}

size_t get_threads(int intval)
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
        throw std::invalid_argument("Invalid number of threads: " + std::to_string(intval));
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
    bool SetOptionCommon(string_view sflag, string_view value)
    {
        string flag{sflag};
        int intval;
        [[maybe_unused]] auto result_i = std::from_chars(input.data(), input.data() + input.size(), intval);

        cout << "flag: " << flag << " value: " << value << endl;
        if (flag == "-i")
            input = value;
        else if (flag == "-j")
            input_pe = value;
        else if (flag == "-o")
            output = value, cout << "matched output" << endl;
        else if (flag == "-op")
            output_pe = value;
        else if (flag == "-M")
            max_memory = try_parse_ull(value, flag) * 1'000'000;
        else if (flag == "-l")
            min_length = intval;
        else if (flag == "-c")
        {
            cluster_thd = try_parse_float(value, flag);
            useIdentity = true;
        }
        else if (flag == "-D")
        {
            distance_thd = try_parse_float(value, flag);
            useIdentity = true;
        }
        else if (flag == "-b")
            band_width = intval;
        else if (flag == "-n")
            NAA = intval;
        else if (flag == "-d")
            des_len = intval;
        else if (flag == "-s")
            diff_cutoff = try_parse_float(value, flag);
        else if (flag == "-S")
            diff_cutoff_aa = intval;
        else if (flag == "-B")
            store_disk = intval;
        else if (flag == "-P")
            PE_mode = intval;
        else if (flag == "-cx")
            trim_len = intval;
        else if (flag == "-cy")
            trim_len_R2 = intval;
        else if (flag == "-ap")
            align_pos = intval;
        else if (flag == "-sc")
            sort_output = intval;
        else if (flag == "-sf")
            sort_outputf = intval;
        else if (flag == "-p")
            print = intval;
        else if (flag == "-g")
            cluster_best = intval;
        else if (flag == "-G")
            global_identity = intval;
        else if (flag == "-aL")
            long_coverage = try_parse_float(value, flag);
        else if (flag == "-AL")
            long_control = intval;
        else if (flag == "-aS")
            short_coverage = try_parse_float(value, flag);
        else if (flag == "-AS")
            short_control = intval;
        else if (flag == "-A")
            min_control = intval;
        else if (flag == "-uL")
            long_unmatch_per = try_parse_float(value, flag);
        else if (flag == "-uS")
            short_unmatch_per = try_parse_float(value, flag);
        else if (flag == "-U")
            unmatch_len = intval;
        else if (flag == "-tmp")
            temp_dir = value.data();
        else if (flag == "-bak")
            backupFile = intval;
        else if (flag == "-T")
            threads = get_threads(intval);
        else
            return false;
        return true;
    }

    bool SetOption(string_view sflag, string_view value)
    {
        if (is454)
        {
            if (sflag == "-s")
                return false;
            else if (sflag == "-S")
                return false;
            else if (sflag == "-G")
                return false;
            else if (sflag == "-A")
                return false;
            else if (sflag == "-r")
                return false;
            else if (sflag == "-D")
            {
                max_indel = atoi(value.data());
                return true;
            }
        }
        if (SetOptionCommon(sflag, value))
            return true;
        if (sflag == "-t")
            tolerance = atoi(value.data());
        else if (sflag == "-F")
            frag_size = atoi(value.data());
        else if (has2D && SetOption2D(sflag.data(), value.data()))
            return true;
        else if (isEST && SetOptionEST(sflag.data(), value.data()))
            return true;
        else
            return false;
        return true;
    }

    bool SetOption2D(const string_view sflag, const string_view value)
    {
        if (SetOptionCommon(sflag, value))
            return true;
        if (sflag == "-i2")
            input2 = value;
        else if (sflag == "-j2")
            input2_pe = value;
        else if (sflag == "-s2")
            diff_cutoff2 = atof(value.data());
        else if (sflag == "-S2")
            diff_cutoff_aa2 = atoi(value.data());
        else
            return false;
        return true;
    }

    bool SetOptionEST(const string_view sflag, const string_view value)
    {
        NAA_top_limit = 12;
        if (SetOptionCommon(sflag, value))
            return true;
        if (sflag == "-r")
            option_r = atoi(value.data());
        else if (sflag == "-gap")
            mat.gap = MAX_SEQ * atoi(value.data());
        else if (sflag == "-gap-ext")
            mat.ext_gap = MAX_SEQ * atoi(value.data());
        else if (sflag == "-match")
            mat.set_match(atoi(value.data()));
        else if (sflag == "-mismatch")
            mat.set_mismatch(atoi(value.data()));
        else if (sflag == "-mask")
        {
            string letters{value};
            for (const auto &ch_raw : letters)
            {
                char ch = toupper(ch_raw);
                if (ch < 'A' || ch > 'Z')
                    continue;
                na2idx[ch - 'A'] = 5;
            }
            setaa_to_na();
        }
        else
            return false;
        return true;
    }

    bool SetOptions(size_t argc, const char **argv, bool twod, bool est)
    {
        const std::span<const char *> args{argv, argc};
        return Options::SetOptions(args, twod, est);
    }

    bool SetOptions(const std::span<const char *> args, bool twod, bool est)
    {
        printf("================================================================\n");
        cout << "Program: CD-HIT, V" CDHIT_VERSION " " << WITH_OPENMP << ", " << __DATE__ << ", " << __TIME__ << endl;
        printf("Command:");
        auto n = 9;
        for (const char *opt : args)
        {
            n += strlen(opt) + 1;
            if (n >= 64)
            {
                printf("\n         %s", opt);
                n = strlen(opt) + 9;
            }
            else
            {
                printf(" %s", opt);
            }
        }
        printf("\n\n");
        time_t tm = time(NULL);
        printf("Started: %s", ctime(&tm));
        printf("================================================================\n");
        printf("                            Output                              \n");
        printf("----------------------------------------------------------------\n");
        has2D = twod;
        isEST = est;

        std::unordered_map<std::string_view, std::string_view> arg_map;
        for (size_t i = 1; i < args.size(); i += 2)
        {
            auto key = args[i];
            auto value = (i + 1 < args.size()) ? args[i + 1] : "";
            arg_map[key] = value;
        }

        for (const auto &[key, value] : arg_map)
        {
            auto erred = SetOption(key, value);
            if (erred)
            {
                cerr << "Error: Invalid option: " << key << " " << value << endl;
                return false;
            }
        }

        return true;
    }

    void Validate()
    {
        if (useIdentity and useDistance)
            bomb_error("can not use both identity cutoff and distance cutoff");
        if (useDistance)
        {
            if ((distance_thd > 1.0) || (distance_thd < 0.0))
                bomb_error("invalid distance threshold");
        }
        else if (isEST)
        {
            if ((cluster_thd > 1.0) || (cluster_thd < 0.8))
                bomb_error("invalid clstr threshold, should >=0.8");
        }
        else
        {
            if ((cluster_thd > 1.0) || (cluster_thd < 0.4))
                bomb_error("invalid clstr");
        }

        if (input.size() == 0)
            bomb_error("no input file");
        if (output.size() == 0)
            bomb_error("no output file");
        if (PE_mode)
        {
            if (input_pe.size() == 0)
                bomb_error("no input file for R2 sequences in PE mode");
            if (output_pe.size() == 0)
                bomb_error("no output file for R2 sequences in PE mode");
        }
        if (isEST && (align_pos == 1))
            option_r = 0;

        if (band_width < 1)
            bomb_error("invalid band width");
        if (NAA < 2 || NAA > NAA_top_limit)
            bomb_error("invalid word length");
        if (des_len < 0)
            bomb_error("too short description, not enough to identify sequences");
        if (not isEST && (tolerance < 0 || tolerance > 5))
            bomb_error("invalid tolerance");
        if ((diff_cutoff < 0) || (diff_cutoff > 1))
            bomb_error("invalid value for -s");
        if (diff_cutoff_aa < 0)
            bomb_error("invalid value for -S");
        if (has2D)
        {
            if ((diff_cutoff2 < 0) || (diff_cutoff2 > 1))
                bomb_error("invalid value for -s2");
            if (diff_cutoff_aa2 < 0)
                bomb_error("invalid value for -S2");
            if (PE_mode)
            {
                if (input2_pe.size() == 0)
                    bomb_error("no input file for R2 sequences for 2nd db in PE mode");
            }
        }
        if (global_identity == 0)
            print = 1;
        if (short_coverage < long_coverage)
            short_coverage = long_coverage;
        if (short_control > long_control)
            short_control = long_control;
        if ((global_identity == 0) && (short_coverage == 0.0) && (min_control == 0))
            bomb_error("You are using local identity, but no -aS -aL -A option");
        if (frag_size < 0)
            bomb_error("invalid fragment size");

#if 0
	if( useDistance ){
		/* when required_aan becomes zero */
		if( distance_thd * NAA >= 1 )
			bomb_warning( "word length is too long for the distance cutoff" );
	}else{
		/* when required_aan becomes zero */
		if( cluster_thd <= 1.0 - 1.0 / NAA )
			bomb_warning( "word length is too long for the identity cutoff" );
	}
#endif

        const char *message = "Your word length is %i, using %i may be faster!\n";
        if (not isEST && tolerance)
        {
            int clstr_idx = (int)(cluster_thd * 100) - naa_stat_start_percent;
            int tcutoff = naa_stat[tolerance - 1][clstr_idx][5 - NAA];

            if (tcutoff < 5)
                bomb_error("Too low cluster threshold for the word length.\n"
                           "Increase the threshold or the tolerance, or decrease the word length.");
            for (size_t i = 5; i > NAA; i--)
            {
                if (naa_stat[tolerance - 1][clstr_idx][5 - i] > 10)
                {
                    printf(message, NAA, i);
                    break;
                }
            }
        }
        else if (isEST)
        {
            if (cluster_thd > 0.9 && NAA < 8)
                printf(message, NAA, 8);
            else if (cluster_thd > 0.87 && NAA < 5)
                printf(message, NAA, 5);
            else if (cluster_thd > 0.80 && NAA < 4)
                printf(message, NAA, 4);
            else if (cluster_thd > 0.75 && NAA < 3)
                printf(message, NAA, 3);
        }
        else
        {
            if (cluster_thd > 0.85 && NAA < 5)
                printf(message, NAA, 5);
            else if (cluster_thd > 0.80 && NAA < 4)
                printf(message, NAA, 4);
            else if (cluster_thd > 0.75 && NAA < 3)
                printf(message, NAA, 3);
        }

        if ((min_length + 1) < NAA)
            bomb_error("Too short -l, redefine it");
    }

    void ComputeTableLimits(size_t min_len, size_t max_len, size_t typical_len, size_t mem_need)
    {
        // liwz Fri Jan 15 15:44:47 PST 2016
        // T=1 scale=1
        // T=2 scale=0.6035
        // T=4 scale=0.375
        // T=8 scale=0.2392
        // T=16 scale=0.1562
        // T=32 scale=0.104
        // T=64 scale=0.0703

        double scale = 0.5 / threads + 0.5 / sqrt(threads);
        max_sequences = (size_t)(scale * MAX_TABLE_SEQ);
        max_entries = (size_t)(scale * (500 * max_len + 500000 * typical_len + 50000000));
        if (max_memory)
        {
            double frac = max_sequences / (double)max_entries;
            max_entries = (options.max_memory - mem_need) / sizeof(IndexCount);
            max_sequences = (size_t)(max_entries * frac);
            if (max_sequences < MAX_TABLE_SEQ / 100)
                max_sequences = MAX_TABLE_SEQ / 100;
            if (max_sequences > MAX_TABLE_SEQ)
                max_sequences = MAX_TABLE_SEQ;
        }
        printf("Table limit with the given memory limit:\n");
        printf("Max number of representatives: %zu\n", max_sequences);
        printf("Max number of word counting entries: %zu\n\n", max_entries);
    }

    // void Options::Print()
    // {
    //     printf("isEST = %i\n", isEST);
    //     printf("has2D = %i\n", has2D);
    //     printf("NAA = %i\n", NAA);
    //     printf("NAA_top_limit = %i\n", NAA_top_limit);
    //     printf("min_length = %i\n", min_length);
    //     printf("cluster_best = %i\n", cluster_best);
    //     printf("global_identity = %i\n", global_identity);
    //     printf("cluster_thd = %g\n", cluster_thd);
    //     printf("diff_cutoff = %g\n", diff_cutoff);
    //     printf("diff_cutoff_aa = %i\n", diff_cutoff_aa);
    //     printf("tolerance = %i\n", tolerance);
    //     printf("long_coverage = %g\n", long_coverage);
    //     printf("long_control = %i\n", long_control);
    //     printf("short_coverage = %g\n", short_coverage);
    //     printf("short_control = %i\n", short_control);
    //     printf("frag_size = %i\n", frag_size);
    //     printf("option_r = %i\n", option_r);
    //     printf("print = %i\n", print);
    // }
};

void setaa_to_na()
{
    int i;
    for (i = 0; i < 26; i++)
        aa2idx[i] = na2idx[i];
}
