export module Options;

import std;

import common;
import ScoreMatrix;

using std::string;
using std::string_view;
using std::printf;
using std::atoi;
using std::cout;
using std::endl;

export string temp_dir = "";

export typedef std::map<string, string> ArgsMap;

enum class match_mask {
	none = 0,
	A = 1,
	C = 2,
	G = 4,
	T = 8,
	U = 16,
	N = 5
};

// idx for                                    A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
export std::array<unsigned int,26> na2idx = {0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4};

// idx for                                    A  B  C  D  E  F   G  H  I  J   K   L   M   N  O   P   Q  R  S   T   U   V   W   X   Y   Z
export std::array<unsigned int,26> aa2idx = {0,2,4,3,6,13,7,8,9,20,11,10,12,2,20,14,5,1,15,16,20,19,17,20,18,6};

// so  aa2idx[ X  - 'A'] => idx_of_X, eg aa2idx['A' - 'A'] => 0,
// and aa2idx['M' - 'A'] => 12

export constexpr void setaa_to_na()
{
	aa2idx = na2idx;
}


size_t get_threads(const std::string& sval) {
	auto intval = std::stoi(sval);
	if(intval < 0) {
		throw std::invalid_argument("Invalid number of threads: " + sval);
	}
	auto cpus = std::thread::hardware_concurrency();
	auto threads = intval;
	if(threads > cpus)
	{
		threads = cpus;
		std::println("Warning: total number of CPUs in the system is {} asked for {}",cpus,intval);
	}
	if(threads == 0)
	{
		threads = cpus;
		std::println("total number of CPUs in the system is %i asked for {}",cpus,intval);
	}
	if(threads != intval)
		std::println("Actual number of CPUs to be used: {}",threads);
	if(threads <= 0) {
		throw std::invalid_argument("Invalid number of threads: " + std::to_string(intval));
	}
	return threads;
}

bool parse_bool(const string& sv) {
    // 1) Empty string => false
    if (sv.empty())
        return false;

    // 2) Make a lowercase copy
    std::string s{sv};
	for (char& c : s)
		c = static_cast<char>(std::tolower(c));

	// std::transform(s.begin(), s.end(), s.begin(), 
    //                [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    // 3) Check known "false" values
    if (s == "0" || s == "false")
        return false;

	// 4) Check known "true" values
	if (s == "1" || s == "true")
		return true;

	// 4) Otherwise explore
	throw std::invalid_argument("Invalid boolean value: " + sv);
}


// non templated statics
export string_view executable_name{};
export const auto start_time = std::chrono::system_clock::now();

export struct StaticOptions {
    // static constinit size_t MAX_UAA;
    static constexpr size_t MAX_GAP = MAX_SEQ;			                     // MAX_GAP <= MAX_SEQ
    static constexpr size_t MAX_DIAG = (MAX_GAP * 2);                        // MAX_DIAG be twice of MAX_SEQ
};
// end non templated statics


export template <bool TisEST>
struct CodeConfig {
    static constexpr bool PE_mode = true;       // -P
    static constexpr bool isEST = TisEST;
};

export template<typename CodeConfig = CodeConfig<false>>
struct Options
{
    static constinit bool has2D;
    static constinit bool is454;
    static constinit bool backupFile;

    bool isEST = CodeConfig::isEST;
    bool PE_mode = CodeConfig::PE_mode;
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
	size_t frag_size = 1;
	int option_r = 1;
	size_t threads = 1;
	int trim_len = 0;                        // -cx
	int trim_len_R2 = 0;                     // -cy
	size_t align_pos = 0;                    // -ap for alignment position
	
	size_t max_entries;
	size_t max_sequences = 1 << 20;
	size_t mem_limit = 100000000;

	string input;
	string input_pe;
	string input2;
	string input2_pe;
	string output;
	string output_pe;

	int sort_output = 0;  // -sc
	int sort_outputf = 0; // -sf

private:
    int est_gap = 0;
    int est_ext_gap = 0;
    int est_match = 0;
    int est_mismatch = 0;

    // No use case for this on the API
    Options(Options&&) = default;

    Options(Options&& defaults, const ArgsMap& args): Options(defaults)
    {
        for(const auto& opt : args) {
            SetOption(opt.first, opt.second);
        }
        mat.set_gaps(est_gap, est_ext_gap);
        mat.set_match(est_match);
        mat.set_mismatch(est_mismatch);

        Validate();
    }

public:
    // container for defaults
    Options() = default;

    Options& operator=(Options&&) = default;
    Options& operator=(const Options&) = default;
    Options(const Options&) = default;

	void Validate()
	{
		std::vector<std::string> errs;
		if(distance_thd > 0.0) {
			if(cluster_thd > 0) {
				errs.push_back("can not use both identity cutoff and distance cutoff");
			}
			if((distance_thd > 1.0) || (distance_thd < 0.0)) {
				errs.push_back("invalid distance threshold");
			}
		} else if constexpr(CodeConfig::isEST) {
			if(!(0.8 <= cluster_thd && cluster_thd <= 1.0)) {
				errs.push_back(std::format("invalid clstr threshold, should by 0.8 <= {} <= 1.0",cluster_thd));
			}
		} else {
			if((cluster_thd > 1.0) || (cluster_thd < 0.4)) {
				errs.push_back("invalid clstr");
			}
		}

		if(input.empty())
			errs.push_back("no input file");
		if(output.empty())
			errs.push_back("no output file");
		if constexpr (CodeConfig::PE_mode) {
			if(input_pe.empty())
				errs.push_back("no input file for R2 sequences in PE mode");
			if(output_pe.empty())
				errs.push_back("no output file for R2 sequences in PE mode");
		}
		if constexpr(CodeConfig::isEST) if (align_pos == 1)	option_r = 0;

		if(band_width < 1)
			errs.push_back("invalid band width");
		if(NAA < 2 || NAA > NAA_top_limit)
			errs.push_back("invalid word length");
		if(des_len < 0)
			errs.push_back("too short description, not enough to identify sequences");
		if constexpr (not CodeConfig::isEST) if (tolerance < 0 || tolerance > 5)
		    errs.push_back("invalid tolerance");
		if((diff_cutoff < 0) || (diff_cutoff > 1))
			errs.push_back("invalid value for -s");
		if(diff_cutoff_aa < 0)
			errs.push_back("invalid value for -S");
		if(has2D) {
			if((diff_cutoff2 < 0) || (diff_cutoff2 > 1))
				errs.push_back("invalid value for -s2");
			if(diff_cutoff_aa2 < 0)
				errs.push_back("invalid value for -S2");
			if constexpr (CodeConfig::PE_mode)
			    if (input2_pe.empty()) errs.push_back("no input file for R2 sequences for 2nd db in PE mode");
		}
		if(global_identity == 0)
			print = 1;
		if(short_coverage < long_coverage)
			short_coverage = long_coverage;
		if(short_control > long_control)
			short_control = long_control;
		if((global_identity == 0) && (short_coverage == 0.0) && (min_control == 0))
			errs.push_back("You are using local identity, but no -aS -aL -A option");
		if(frag_size < 0)
			errs.push_back("invalid fragment size");
		constexpr auto message = "Your word length is {}, using {} may be faster!";

		#if 0
		if(useDistance){
			/* when required_aan becomes zero */
			if(distance_thd * NAA >= 1)
				std::println(std::cerr, "\nWarning:\nword length is too long for the distance cutoff\nNot fatal, but may affect results !!\n");
		} else{
			/* when required_aan becomes zero */
			if(cluster_thd <= 1.0 - 1.0 / NAA)
				std::println(std::cerr, "\nWarning:\nword length is too long for the identity cutoff\nNot fatal, but may affect results !!\n");
		}

		if(not isEST && tolerance)
		{
			int clstr_idx = (int)(cluster_thd * 100) - naa_stat_start_percent;
			int tcutoff = naa_stat[tolerance - 1][clstr_idx][5 - NAA];

			if(tcutoff < 5)
				errs.emplace("Too low cluster threshold for the word length.\n"
						   "Increase the threshold or the tolerance, or decrease the word length.");
			for(size_t i = 5; i > NAA; i--)
			{
				if(naa_stat[tolerance - 1][clstr_idx][5 - i] > 10)
				{
					printf(message,NAA,i);
					break;
				}
			}
		}
		#endif

		if constexpr(CodeConfig::isEST)
		{
			if(cluster_thd > 0.9 && NAA < 8)
				std::println(message,NAA,8);
			else if(cluster_thd > 0.87 && NAA < 5)
				std::println(message,NAA,5);
			else if(cluster_thd > 0.80 && NAA < 4)
				std::println(message,NAA,4);
			else if(cluster_thd > 0.75 && NAA < 3)
				std::println(message,NAA,3);
		} else
		{
			if(cluster_thd > 0.85 && NAA < 5)
				std::println(message,NAA,5);
			else if(cluster_thd > 0.80 && NAA < 4)
				std::println(message,NAA,4);
			else if(cluster_thd > 0.75 && NAA < 3)
				std::println(message,NAA,3);
		}

		if((min_length + 1) < NAA)
			errs.push_back("Too short -l, redefine it");

		if(!errs.empty()) {
			std::cerr << "Errors:\n";
			for(const auto &err : errs)
				std::cerr << "  " << err << "\n";
			throw std::invalid_argument("Invalid options");
		}
	}

        constexpr bool useDistance() const {
        return (distance_thd > 0.0);
    }
	bool SetOptionCommon(const string_view flag, const string& value)
	{
		std:println("flag: {} value: {}", flag, value);

		if(flag == "-i")
			input = value;
		else if(flag == "-j")
			input_pe = value;
		else if(flag == "-o")
			output = value;
		else if(flag == "-op")
			output_pe = value;
		else if(flag == "-M")
			max_memory = std::stoll(value) * MEGA_MiBi;
		else if(flag == "-l")
			min_length = std::stoi(value);
		else if(flag == "-c")
			cluster_thd = std::stod(value);
		else if(flag == "-D")
			distance_thd = std::stod(value);
		else if(flag == "-b")
			band_width = std::stoi(value);
		else if(flag == "-n")
			NAA = std::stoi(value);
		else if(flag == "-d")
			des_len = std::stoi(value);
		else if(flag == "-s")
			diff_cutoff = std::stod(value);
		else if(flag == "-S")
			diff_cutoff_aa = std::stoi(value);
		else if(flag == "-B")
			store_disk = std::stoi(value);
		else if(flag == "-P")
			PE_mode = parse_bool(value);
		else if(flag == "-cx")
			trim_len = std::stoi(value);
		else if(flag == "-cy")
			trim_len_R2 = std::stoi(value);
		else if(flag == "-ap")
			align_pos = std::stoi(value);
		else if(flag == "-sc")
			sort_output = std::stoi(value);
		else if(flag == "-sf")
			sort_outputf = std::stoi(value);
		else if(flag == "-p")
			print = std::stoi(value);
		else if(flag == "-g")
			cluster_best = std::stoi(value);
		else if(flag == "-G")
			global_identity = std::stoi(value);
		else if(flag == "-aL")
			long_coverage = std::stod(value);
		else if(flag == "-AL")
			long_control = std::stoi(value);
		else if(flag == "-aS")
			short_coverage = std::stod(value);
		else if(flag == "-AS")
			short_control = std::stoi(value);
		else if(flag == "-A")
			min_control = std::stoi(value);
		else if(flag == "-uL")
			long_unmatch_per = std::stod(value);
		else if(flag == "-uS")
			short_unmatch_per = std::stod(value);
		else if(flag == "-U")
			unmatch_len = std::stoi(value);
		else if(flag == "-tmp")
			temp_dir = value;
		else if(flag == "-bak")
			backupFile = parse_bool(value);
		else if(flag == "-T")
			threads = get_threads(value);
		else
			return false;
		return true;
	}

	bool SetOption(const string& sflag, const string& value)
	{
		if(is454)
		{
			if(sflag == "-s")
				return false;
			if(sflag == "-S")
				return false;
			if(sflag == "-G")
				return false;
			if(sflag == "-A")
				return false;
			if(sflag == "-r")
				return false;
			if(sflag == "-D")
			{
				max_indel = std::stoi(value);
				return true;
			}
		}
		if(SetOptionCommon(sflag,value))
			return true;
		if(sflag == "-t")
			tolerance = std::stoi(value);
		else if(sflag == "-F")
			frag_size = std::stoi(value);
		else if(has2D && SetOption2D(sflag,value))
			return true;
		else if constexpr(CodeConfig::isEST) if (SetOptionEST(sflag,value))
			return true;
		else
			return false;
		return true;
	}

	bool SetOption2D(const string_view sflag, const string& value)
	{
		if(SetOptionCommon(sflag,value))
			return true;
		if(sflag == "-i2")
			input2 = value;
		else if(sflag == "-j2")
			input2_pe = value;
		else if(sflag == "-s2")
			diff_cutoff2 = std::atof(value.data());
		else if(sflag == "-S2")
			diff_cutoff_aa2 = atoi(value.data());
		else
			return false;
		return true;
	}

	bool SetOptionEST(const string_view sflag, const string& value)
	{
		NAA_top_limit = 12;
		if(SetOptionCommon(sflag,value))
			return true;
		if(sflag == "-r")
			option_r = atoi(value.data());
		else if(sflag == "-mask") {
		    string letters{value};
		    for(const auto ch_raw : letters)
		    {
		        char ch = std::toupper(ch_raw);
		        size_t idx = ch - 'A';
		        if(ch < 'A' || ch > 'Z') {
		            std::cerr << "Error: Invalid mask letter: " << ch_raw << endl;
		            continue;
		        }
		        na2idx[idx] = static_cast<size_t>(match_mask::N);
		    }
		    setaa_to_na();
		} else if(sflag == "-gap")
			est_gap = MAX_SEQ * atoi(value.data());
		else if(sflag == "-gap-ext")
			est_ext_gap = MAX_SEQ * atoi(value.data());
		else if(sflag == "-match")
			est_match = (atoi(value.data()));
		else if(sflag == "-mismatch")
			est_mismatch = (atoi(value.data()));
		else
			return false;
		return true;
	}


	void ComputeTableLimits(size_t min_len,size_t max_len,size_t typical_len,size_t mem_need);

    static ArgsMap&& ParseArgs(std::span<const char*> all_args)
    {
        executable_name = all_args[0];
        auto args = all_args.subspan(1, all_args.size() - 1);
        if (args.size() < 4)
            throw std::invalid_argument(std::format("Invalid number of arguments {}", args.size()));

        ArgsMap arg_map;
        for(size_t i = 0; i < args.size(); i++)
        {
            string key = args[i];
            if(!key.starts_with('-'))
                throw std::invalid_argument(std::format("Error: weird arg, should start with '-': {} = {}", i, key));
            if(i + 1 >= args.size())
                throw std::invalid_argument(std::format("Error: Missing value for last option: {} = {}", i, key));

            string value{args[++i]};
            if(value.starts_with('-'))
                throw std::invalid_argument(std::format("Error: all parameters need a value {} = {}, {} = {}", i - 1, key, i, value));

            arg_map[std::move(key)] = std::move(value);
        }
        return std::move(arg_map);
    }

    static void RollStart(std::span<const char*> args);
};

constinit bool Options<CodeConfig<true>>::has2D = false;
constinit bool Options<CodeConfig<true>>::is454 = false;
constinit bool Options<CodeConfig<true>>::backupFile = false;

constinit bool Options<CodeConfig<false>>::has2D = false;
constinit bool Options<CodeConfig<false>>::is454 = false;
constinit bool Options<CodeConfig<false>>::backupFile = false;


export Options options;

template<typename CodeConfig = CodeConfig<false>>
void Options<CodeConfig>::RollStart(std::span<const char*> all_args) {
    std::println("================================================================");
    std::println("Program: CD-HIT, V{} {}, {}", CDHIT_VERSION, WITH_OPENMP, __DATE__);
    std::println("Started: {0:%c} {0:%Z}", start_time);

    const ArgsMap& argsmaps = Options::ParseArgs(all_args);

    std::println("Command line arguments:");
    std::println("Executable: {}", executable_name);
    for (const auto& arg : argsmaps) {
        std::println("  {} = {}", arg.first, arg.second);
    }

    // double spacer
    std::println("\n");

    std::println("================================================================");
    std::println("                            Output                              ");
    std::println("----------------------------------------------------------------");

    options = std::move(Options<CodeConfig>(std::move(options), argsmaps));
}


template<typename CodeConfig>
void Options<CodeConfig>::ComputeTableLimits(size_t min_len, const size_t max_len, const size_t typical_len, const size_t mem_need) {
    // liwz Fri Jan 15 15:44:47 PST 2016
    // T=1 scale=1
    // T=2 scale=0.6035
    // T=4 scale=0.375
    // T=8 scale=0.2392
    // T=16 scale=0.1562
    // T=32 scale=0.104
    // T=64 scale=0.0703

    double scale = 0.5 / threads + 0.5 / std::sqrt(threads);
    max_sequences = static_cast<size_t>(scale * MAX_TABLE_SEQ);
    max_entries = static_cast<size_t>(scale * (500 * max_len + 500000 * typical_len + 50000000));
    if (max_memory) {
        double frac = max_sequences / static_cast<double>(max_entries);
        max_entries = (options.max_memory - mem_need) / sizeof(IndexCount);
        max_sequences = static_cast<size_t>(max_entries * frac);
        if (max_sequences < MAX_TABLE_SEQ / 100)
            max_sequences = MAX_TABLE_SEQ / 100;
        if (max_sequences > MAX_TABLE_SEQ)
            max_sequences = MAX_TABLE_SEQ;
    }
    printf("Table limit with the given memory limit:\n");
    printf("Max number of representatives: %zu\n", max_sequences);
    printf("Max number of word counting entries: %zu\n\n", max_entries);
}


// information

consteval std::string cd_hit_ver() {
	std::string s = std::string{"\t\t====== CD-HIT version "};
	s += CDHIT_VERSION;
	s += "built on ";
	s += BUILD_DATE;
	s += ") ======";
	return s;
}


constexpr std::string_view cd_hit_ref1 = R"a("CD-HIT: a fast program for clustering and comparing large sets of protein or nucleotide sequences", Weizhong Li & Adam Godzik. Bioinformatics, (2006) 22:1658-1659
"CD-HIT: accelerated for clustering the next generation sequencing data", Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu & Weizhong Li. Bioinformatics, (2012) 28:3150-3152;)a";
constexpr std::string_view cd_hit_ref3 = R"a("Beifang Niu, Limin Fu, Shulei Sun and Weizhong Li. Artificial and natural duplicates in pyrosequencing reads of metagenomic data. BMC Bioinformatics (2010) 11:187")a";
//

char contacts[] =
R"a(
	Questions, bugs, contact Weizhong Li at liwz@sdsc.edu
	For updated versions and information, please visit: http://cd-hit.org
	                                                 or https://github.com/weizhongli/cdhit

	cd-hit web server is also available from http://cd-hit.org
	If you find cd-hit useful, please kindly cite:

)a";
char txt_option_i[] = "\tinput filename in fasta format, required, can be in .gz format\n";
char txt_option_j[] =
"\tinput filename in fasta/fastq format for R2 reads if input are paired end (PE) files\n \
\t -i R1.fq -j R2.fq -o output_R1 -op output_R2 or\n \
\t -i R1.fa -j R2.fa -o output_R1 -op output_R2 \n";
char txt_option_i_2d[] = "\tinput filename for db1 in fasta format, required, can be in .gz format\n";
char txt_option_i2[] = "\tinput filename for db2 in fasta format, required, can be in .gz format\n";
char txt_option_j2[] =
"\tinput filename in fasta/fastq format for R2 reads if input are paired end (PE) files\n \
\t -i db1-R1.fq -j db1-R2.fq -i2 db2-R1.fq -j2 db2-R2.fq -o output_R1 -op output_R2 or\n \
\t -i db1-R1.fa -j db1-R2.fa -i2 db2-R1.fq -j2 db2-R2.fq -o output_R1 -op output_R2 \n";
char txt_option_o[] = "\toutput filename, required\n";
char txt_option_op[] = "\toutput filename for R2 reads if input are paired end (PE) files\n";
char txt_option_c[] =
"\tsequence identity threshold, default 0.9\n \
\tthis is the default cd-hit's \"global sequence identity\" calculated as:\n \
\tnumber of identical amino acids or bases in alignment\n \
\tdivided by the full length of the shorter sequence\n";
char txt_option_G[] =
"\tuse global sequence identity, default 1\n \
\tif set to 0, then use local sequence identity, calculated as :\n \
\tnumber of identical amino acids or bases in alignment\n \
\tdivided by the length of the alignment\n \
\tNOTE!!! don't use -G 0 unless you use alignment coverage controls\n \
\tsee options -aL, -AL, -aS, -AS\n";
char txt_option_g[] =
"\t1 or 0, default 0\n \
\tby cd-hit's default algorithm, a sequence is clustered to the first \n \
\tcluster that meet the threshold (fast cluster). If set to 1, the program\n \
\twill cluster it into the most similar cluster that meet the threshold\n \
\t(accurate but slow mode)\n \
\tbut either 1 or 0 won't change the representatives of final clusters\n";
char txt_option_b[] = "\tband_width of alignment, default 20\n";
char txt_option_M[] = "\tmemory limit (in MB) for the program, default 800; 0 for unlimitted;\n";
char txt_option_n[] = "\tword_length, default 5, see user's guide for choosing it\n";
char txt_option_n_est[] = "\tword_length, default 10, see user's guide for choosing it\n";
char txt_option_l[] = "\tlength of throw_away_sequences, default 10\n";
char txt_option_t[] = "\ttolerance for redundance, default 2\n";
char txt_option_T[] = "\tnumber of threads, default 1; with 0, all CPUs will be used\n";
char txt_option_d[] =
"\tlength of description in .clstr file, default 20\n \
\tif set to 0, it takes the fasta defline and stops at first space\n";
char txt_option_s[] =
"\tlength difference cutoff, default 0.0\n \
\tif set to 0.9, the shorter sequences need to be\n \
\tat least 90% length of the representative of the cluster\n";
char txt_option_S[] =
"\tlength difference cutoff in amino acid, default 999999\n \
\tif set to 60, the length difference between the shorter sequences\n \
\tand the representative of the cluster can not be bigger than 60\n";
char txt_option_s2[] =
"\tlength difference cutoff for db1, default 1.0\n \
\tby default, seqs in db1 >= seqs in db2 in a same cluster\n \
\tif set to 0.9, seqs in db1 may just >= 90% seqs in db2\n";
char txt_option_S2[] =
"\tlength difference cutoff, default 0\n \
\tby default, seqs in db1 >= seqs in db2 in a same cluster\n \
\tif set to 60, seqs in db2 may 60aa longer than seqs in db1\n";
char txt_option_aL[] =
"\talignment coverage for the longer sequence, default 0.0\n \
\tif set to 0.9, the alignment must covers 90% of the sequence\n";
char txt_option_AL[] =
"\talignment coverage control for the longer sequence, default 99999999\n \
\tif set to 60, and the length of the sequence is 400,\n \
\tthen the alignment must be >= 340 (400-60) residues\n";
char txt_option_aS[] =
"\talignment coverage for the shorter sequence, default 0.0\n \
\tif set to 0.9, the alignment must covers 90% of the sequence\n";
char txt_option_AS[] =
"\talignment coverage control for the shorter sequence, default 99999999\n \
\tif set to 60, and the length of the sequence is 400,\n \
\tthen the alignment must be >= 340 (400-60) residues\n";
char txt_option_A[] =
"\tminimal alignment coverage control for the both sequences, default 0\n \
\talignment must cover >= this value for both sequences \n";
char txt_option_B[] =
"\t1 or 0, default 0, by default, sequences are stored in RAM\n \
\tif set to 1, sequence are stored on hard drive\n \
\t!! No longer supported !!\n";
char txt_option_P[] =
"\tinput paired end (PE) reads, default 0, single file\n \
\tif set to 1, please use -i R1 -j R2 to input both PE files\n";
char txt_option_cx[] =
"\tlength to keep after trimming the tail of sequence, default 0, not trimming\n \
\tif set to 50, the program only uses the first 50 letters of input sequence\n";
char txt_option_cy[] =
R"a(	length to keep after trimming the tail of R2 sequence, default 0, not trimming
	if set to 50, the program only uses the first 50 letters of input R2 sequence
	e.g. -cx 100 -cy 80 for paired end reads
)a";
char txt_option_ap[] =
R"a(	alignment position constrains,  default 0, no constrain
	if set to 1, the program will force sequences to align at beginnings
	when set to 1, the program only does +/+ alignment
)a";
char txt_option_uL[] =
R"a(	maximum unmatched percentage for the longer sequence, default 1.0
	if set to 0.1, the unmatched region (excluding leading and tailing gaps)
	must not be more than 10% of the sequence
)a";
char txt_option_uS[] =
R"a(	maximum unmatched percentage for the shorter sequence, default 1.0
	if set to 0.1, the unmatched region (excluding leading and tailing gaps)
	must not be more than 10% of the sequence
)a";
char txt_option_U[] =
R"a(	maximum unmatched length, default 99999999
	if set to 10, the unmatched region (excluding leading and tailing gaps)
	must not be more than 10 bases
)a";
char txt_option_p[] = "\t1 or 0, default 0\n \tif set to 1, print alignment overlap in .clstr file\n";
char txt_option_r[] =
R"a(	1 or 0, default 1, by default do both +/+ & +/- alignments
	if set to 0, only +/+ strand alignment
)a";
char txt_option_bak[] = "	write backup cluster file (1 or 0, default 0)\n";
char txt_option_sc[] =
R"a("	sort clusters by size (number of sequences), default 0, output clusters by decreasing length
	if set to 1, output clusters by decreasing size
)a";
char txt_option_sf[] =
R"a(	sort fasta/fastq by cluster size (number of sequences), default 0, no sorting
	if set to 1, output sequences by decreasing cluster size
	this can be very slow if the input is in .gz format
)a";

char txt_option_mask[] = "\tmasking letters (e.g. -mask NX, to mask out both 'N' and 'X')\n";
char txt_option_match[] = "\tmatching score, default 2 (1 for T-U and N-N)\n";
char txt_option_match2[] = "\tmatching score, default 2\n";
char txt_option_mismatch[] = "\tmismatching score, default -2\n";
char txt_option_mismatch2[] = "\tmismatching score, default -1\n";
char txt_option_gap[] = "\tgap opening score, default -6\n";
char txt_option_gap2[] = "\tgap opening score, default -3\n";
char txt_option_gap_ext[] = "\tgap extension score, default -1\n";
char mytxt_option_c[] =
R"a(	sequence identity threshold, default 0.98
	this is a "global sequence identity" calculated as :
	number of identical amino acids or bases in alignment
	divided by the full length of the shorter sequence + gaps
)a";
char mytxt_option_b[] = "\tband_width of alignment, default 10\n";
char mytxt_option_n_est[] = "\tword_length, default 10, see user's guide for choosing it\n";
char mytxt_option_D[] = "\tmax size per indel, default 1\n";

export inline int print_usage(const char* arg) {
	cout << cd_hit_ver() << "\n\n";
	cout << "Usage: " << arg << " [Options] \n\nOptions\n\n";
	cout << "   -i" << txt_option_i;
	cout << "   -o" << txt_option_o;
	cout << "   -c" << txt_option_c;
	cout << "   -G" << txt_option_G;
	cout << "   -b" << txt_option_b;
	cout << "   -M" << txt_option_M;
	cout << "   -T" << txt_option_T;
	cout << "   -n" << txt_option_n;
	cout << "   -l" << txt_option_l;
	cout << "   -t" << txt_option_t;
	cout << "   -d" << txt_option_d;
	cout << "   -s" << txt_option_s;
	cout << "   -S" << txt_option_S;
	cout << "   -aL" << txt_option_aL;
	cout << "   -AL" << txt_option_AL;
	cout << "   -aS" << txt_option_aS;
	cout << "   -AS" << txt_option_AS;
	cout << "   -A" << txt_option_A;
	cout << "   -uL" << txt_option_uL;
	cout << "   -uS" << txt_option_uS;
	cout << "   -U" << txt_option_U;
	cout << "   -B" << txt_option_B;
	cout << "   -p" << txt_option_p;
	cout << "   -g" << txt_option_g;
	cout << "   -sc" << txt_option_sc;
	cout << "   -sf" << txt_option_sf;
	cout << "   -bak" << txt_option_bak;
	cout << "   -h\tprint this help\n\n";
	cout << contacts;
	cout << "   " << cd_hit_ref1 << endl;
	std::println("\n\n");
	std::exit(1);
}


export inline int print_usage_2d(const char* arg) {
	cout << cd_hit_ver() << "\n\n";
	cout << "Usage: " << arg << " [Options] \n\nOptions\n\n";
	cout << "   -i" << txt_option_i_2d;
	cout << "   -i2" << txt_option_i2;
	cout << "   -o" << txt_option_o;
	cout << "   -c" << txt_option_c;
	cout << "   -G" << txt_option_G;
	cout << "   -b" << txt_option_b;
	cout << "   -M" << txt_option_M;
	cout << "   -T" << txt_option_T;
	cout << "   -n" << txt_option_n;
	cout << "   -l" << txt_option_l;
	cout << "   -t" << txt_option_t;
	cout << "   -d" << txt_option_d;
	cout << "   -s" << txt_option_s;
	cout << "   -S" << txt_option_S;
	cout << "   -s2" << txt_option_s2;
	cout << "   -S2" << txt_option_S2;
	cout << "   -aL" << txt_option_aL;
	cout << "   -AL" << txt_option_AL;
	cout << "   -aS" << txt_option_aS;
	cout << "   -AS" << txt_option_AS;
	cout << "   -A" << txt_option_A;
	cout << "   -uL" << txt_option_uL;
	cout << "   -uS" << txt_option_uS;
	cout << "   -U" << txt_option_U;
	cout << "   -B" << txt_option_B;
	cout << "   -p" << txt_option_p;
	cout << "   -g" << txt_option_g;
	cout << "   -bak" << txt_option_bak;
	cout << "   -h\tprint this help\n\n";
	cout << "   Questions, bugs, contact Weizhong Li at liwz@sdsc.edu\n\n";
	cout << "   If you find cd-hit useful, please kindly cite:\n\n";
	cout << "   " << cd_hit_ref1 << endl;
	std::println("\n\n");
	std::exit(1);
}


export inline int print_usage_est(const std::string_view arg) {
	cout << cd_hit_ver() << "\n\n";
	cout << "Usage: " << arg << " [Options] \n\nOptions\n\n";
	cout << "   -i" << txt_option_i;
	cout << "   -j" << txt_option_j;
	cout << "   -o" << txt_option_o;
	cout << "   -op" << txt_option_op;
	cout << "   -c" << txt_option_c;
	cout << "   -G" << txt_option_G;
	cout << "   -b" << txt_option_b;
	cout << "   -M" << txt_option_M;
	cout << "   -T" << txt_option_T;
	cout << "   -n" << txt_option_n_est;
	cout << "   -l" << txt_option_l;
	cout << "   -d" << txt_option_d;
	cout << "   -s" << txt_option_s;
	cout << "   -S" << txt_option_S;
	cout << "   -aL" << txt_option_aL;
	cout << "   -AL" << txt_option_AL;
	cout << "   -aS" << txt_option_aS;
	cout << "   -AS" << txt_option_AS;
	cout << "   -A" << txt_option_A;
	cout << "   -uL" << txt_option_uL;
	cout << "   -uS" << txt_option_uS;
	cout << "   -U" << txt_option_U;
	cout << "   -B" << txt_option_B;
	cout << "   -P" << txt_option_P;
	cout << "   -cx" << txt_option_cx;
	cout << "   -cy" << txt_option_cy;
	cout << "   -ap" << txt_option_ap;
	cout << "   -p" << txt_option_p;
	cout << "   -g" << txt_option_g;
	cout << "   -r" << txt_option_r;
	cout << "   -mask" << txt_option_mask;
	cout << "   -match" << txt_option_match;
	cout << "   -mismatch" << txt_option_mismatch;
	cout << "   -gap" << txt_option_gap;
	cout << "   -gap-ext" << txt_option_gap_ext;
	cout << "   -bak" << txt_option_bak;
	cout << "   -sc" << txt_option_sc;
	cout << "   -sf" << txt_option_sf;
	cout << "   -h\tprint this help\n\n";
	cout << contacts;
	cout << "   " << cd_hit_ref1 << endl;
	std::println("\n\n");
	std::exit(1);
}


export inline int print_usage_est_2d(const char* arg) {
	cout << cd_hit_ver() << "\n\n";
	cout << "Usage: " << arg << " [Options] \n\nOptions\n\n";
	cout << "   -i" << txt_option_i_2d;
	cout << "   -i2" << txt_option_i2;
	cout << "   -j, -j2" << txt_option_j2;
	cout << "   -o" << txt_option_o;
	cout << "   -op" << txt_option_op;
	cout << "   -c" << txt_option_c;
	cout << "   -G" << txt_option_G;
	cout << "   -b" << txt_option_b;
	cout << "   -M" << txt_option_M;
	cout << "   -T" << txt_option_T;
	cout << "   -n" << txt_option_n_est;
	cout << "   -l" << txt_option_l;
	cout << "   -d" << txt_option_d;
	cout << "   -s" << txt_option_s;
	cout << "   -S" << txt_option_S;
	cout << "   -s2" << txt_option_s2;
	cout << "   -S2" << txt_option_S2;
	cout << "   -aL" << txt_option_aL;
	cout << "   -AL" << txt_option_AL;
	cout << "   -aS" << txt_option_aS;
	cout << "   -AS" << txt_option_AS;
	cout << "   -A" << txt_option_A;
	cout << "   -uL" << txt_option_uL;
	cout << "   -uS" << txt_option_uS;
	cout << "   -U" << txt_option_U;
	cout << "   -B" << txt_option_B;
	cout << "   -P" << txt_option_P;
	cout << "   -cx" << txt_option_cx;
	cout << "   -cy" << txt_option_cy;
	cout << "   -p" << txt_option_p;
	cout << "   -g" << txt_option_g;
	cout << "   -r" << txt_option_r;
	cout << "   -mask" << txt_option_mask;
	cout << "   -match" << txt_option_match;
	cout << "   -mismatch" << txt_option_mismatch;
	cout << "   -gap" << txt_option_gap;
	cout << "   -gap-ext" << txt_option_gap_ext;
	cout << "   -bak" << txt_option_bak;
	cout << "   -h\tprint this help\n\n";
	cout << contacts;
	cout << "   " << cd_hit_ref1 << endl;
	std::println("\n\n");
	std::exit(1);
}


export inline int print_usage_div(const char* arg) {
	cout << cd_hit_ver() << "\n\n";
	cout << "Usage: " << arg << " [Options] \n\nOptions\n\n";
	cout << "Options " << endl << endl;
	cout << "   -i in_dbname, required" << endl;
	cout << "   -o out_dbname, required" << endl;
	cout << "   -div number of divide, required " << endl;
	//  cout << "   -dbmax max size of your db\n\n\n";
	std::exit(1);
}


export inline int print_usage_454(const char* arg)
{
	cout << cd_hit_ver() << "\n\n";
	cout << "Usage: " << arg << " [Options] \n\nOptions\n\n";
	cout << "   -i" << txt_option_i;
	cout << "   -o" << txt_option_o;
	cout << "   -c" << mytxt_option_c;
	cout << "   -b" << mytxt_option_b;
	cout << "   -M" << txt_option_M;
	cout << "   -T" << txt_option_T;
	cout << "   -n" << mytxt_option_n_est;
	cout << "   -aL" << txt_option_aL;
	cout << "   -AL" << txt_option_AL;
	cout << "   -aS" << txt_option_aS;
	cout << "   -AS" << txt_option_AS;
	cout << "   -B" << txt_option_B;
	cout << "   -g" << txt_option_g;
	cout << "   -D" << mytxt_option_D;
	cout << "   -match" << txt_option_match2;
	cout << "   -mismatch" << txt_option_mismatch2;
	cout << "   -gap" << txt_option_gap2;
	cout << "   -gap-ext" << txt_option_gap_ext;
	cout << "   -bak" << txt_option_bak;
	cout << "   -h\tprint this help\n\n";
	cout << "   Questions, bugs, contact Weizhong Li at liwz@sdsc.edu\n\n";
	cout << "   If you find cd-hit useful, please kindly cite:\n\n";
	cout << "   " << cd_hit_ref1 << endl;
	std::println("\n\n");
	cout << "   " << cd_hit_ref3 << "\n\n\n";
	std::exit(1);
}

#ifdef _MAX_SEQ_
auto Options::MAX_SEQ = _MAX_SEQ_;
#endif
#undef _MAX_SEQ_

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
