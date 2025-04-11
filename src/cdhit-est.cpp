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

constexpr size_t MAX_UAA = 4;
constexpr size_t MAX_SEQ = 200;

import common;
import Options;
import SequenceDB;
import ScoreMatrix;

static SequenceDB seq_db;

// Figure out how to do this statically
// static ScoreMatrix<true> mat{};
// mat.set_to_na(); // mat.set_gap(-6,-1);

////////////////////////////////////  MAIN /////////////////////////////////////
int main(std::size_t _argc, const char* _argv[])
{
	const auto begin_time = std::time(nullptr);

    const std::span all_args{_argv, _argc};
	const std::string_view prog_name{all_args[0]};
    auto args_span = all_args.subspan(1, all_args.size() - 1);
    if (args_span.size() < 4)
		print_usage_est(prog_name);

	options.cluster_thd = 0.95;
	options.NAA = 10;
	options.NAAN = NAA8;
	seq_db.NAAN = NAA8;
	options.NAA_top_limit = 12;
	setaa_to_na();

    const std::vector<std::string_view> args(args_span.begin(), args_span.end());
	if (!options.SetOptions(args, false, true))
		print_usage_est(prog_name);
	options.Validate();

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
	
	const auto end_time = std::time(nullptr);
	const auto elapsed_time = end_time - begin_time;
	std::println("program completed in {} seconds", elapsed_time);
	return 0;
}
