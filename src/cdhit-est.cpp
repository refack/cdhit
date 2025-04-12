// =============================================================================
// CD-HI-EST
// http://cd-hit.org/
// Cluster Database at High Identity (EST version)
// modified from CD-HI
//
// program written by
//                    Weizhong Li
//                    UCSD, San Diego Supercomputer Center
//                    La Jolla, CA, 92093
//                    Email liwz@sdsc.edu
//                 at
//                    Adam Godzik's lab
//                    The Burnham Institute
//                    La Jolla, CA, 92037
//                    Email adam@burnham-inst.org
//
// Modified by:
//                    Limin Fu
//                    Center for Research in Biological Systems (CRBS), UCSD
//                    La Jolla, CA, 92093
//                    Email: l2fu@ucsd.edu, fu@daovm.net
// =============================================================================

#define SHORT_SEQ

import std;

import common;
import Options;
import SequenceDB;
import ScoreMatrix;

static SequenceDB seq_db;

// Figure out how to do this statically
// static ScoreMatrix<true> mat{};
// mat.set_to_na(); // mat.set_gap(-6,-1);

// Compile time consts (better #define MACROS)
// constinit auto ScoreMatrixOptions::MAX_SEQ = 200;
using SetOptions = Options<CodeConfig<false>>;

////////////////////////////////////  MAIN /////////////////////////////////////
// template<typename SetOptions>
int main(std::size_t argc, const char* argv[])
{
    SetGlobals(4);

    options.cluster_thd = 0.95;
	options.NAA = 10;
	options.NAAN = NAA8;
	seq_db.NAAN = NAA8;
	options.NAA_top_limit = 12;
	setaa_to_na();

    try {
        SetOptions::RollStart({argv, argc});
    }catch (std::invalid_argument& e) {
        std::cerr << std::endl << e.what() << std::endl;
        std::print("run {} -h for help\n", executable_name);
        std::exit(1);
    }

 //    const ArgsMap& args = options.ParseArgs(args_span);
	// if (!options.SetOptions(args, false, true))
	// 	print_usage_est(prog_name);
	// options.Validate();

	seq_db.NAAN = NAAN_array[options.NAA];

	if (options.option_r) {
		Comp_AAN_idx.resize(seq_db.NAAN);
		make_comp_short_word_index(options.NAA);
	}

	if (options.PE_mode) {
		seq_db.Read(options.input, options.input_pe);
	} else {
		seq_db.Read(options.input);
	}

	std::cout << "total seq: " << seq_db.sequences.size() << std::endl;
	seq_db.SortDivide();
	seq_db.DoClustering();

	std::println("writing new database");
	if (options.PE_mode) {
		seq_db.WriteClusters(options.input, options.input_pe, options.output, options.output_pe);
	} else {
		seq_db.WriteClusters(options.input, options.output);
	}

	// write a backup clstr file in case next step crashes
	seq_db.WriteExtra1D();
	
	const auto end_time = std::chrono::system_clock::now();
	const auto elapsed_time = end_time - start_time;
	std::println("program completed in {} seconds", elapsed_time);
	return 0;
}
