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

#include "lib/cdhit-common.h"
import std;
import Options;

// over-write some defs in cd-hi.h for est version
#undef MAX_UAA
constexpr size_t MAX_UAA = 4;
// over-write some defs in cd-hi-init.h for est version

SequenceDB seq_db;

////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, const char **argv)
{
	if (argc < 5)
		print_usage_est(argv[0]);
	const std::span<const char*> args(argv, argc);

	options.cluster_thd = 0.95;
	options.NAA = 10;
	options.NAAN = NAA8;
	seq_db.NAAN = NAA8;
	options.NAA_top_limit = 12;
	setaa_to_na();
	mat.set_to_na(); // mat.set_gap(-6,-1);

	const float begin_time = current_time();

	// ***********************************    parse command line and open file
	if (options.SetOptions(argc, argv, false, true))
		print_usage_est(argv[0]);
	options.Validate();

	const auto db_in = options.input;
	const auto db_in_pe = options.input_pe;
	const auto db_out = options.output;
	const auto db_out_pe = options.output_pe;

	seq_db.NAAN = NAAN_array[options.NAA];

	if (options.option_r)
	{
		Comp_AAN_idx.resize(seq_db.NAAN);
		make_comp_short_word_index(options.NAA);
	}

	if (options.PE_mode)
	{
		seq_db.Read(db_in.c_str(), db_in_pe.c_str());
	}
	else
	{
		seq_db.Read(db_in.c_str());
	}

	cout << "total seq: " << seq_db.sequences.size() << endl;
	seq_db.SortDivide();
	seq_db.DoClustering();

	printf("writing new database\n");
	if (options.PE_mode)
	{
		seq_db.WriteClusters(db_in.c_str(), db_in_pe.c_str(), db_out.c_str(), db_out_pe.c_str());
	}
	else
	{
		seq_db.WriteClusters(db_in.c_str(), db_out.c_str());
	}

	// write a backup clstr file in case next step crashes
	seq_db.WriteExtra1D();
	cout << "program completed !" << endl
		 << endl;
	const auto end_time = current_time();
	printf("Total CPU time %.2f\n", end_time - begin_time);
	return 0;
}
