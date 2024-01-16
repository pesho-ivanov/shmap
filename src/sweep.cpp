#include "sweep.h"

int main(int argc, char **argv) {
	initialize_LUT();

	Counters C;
	Timers T;

	T.start("total");
	C.inc("seeds_limit_reached", 0);
	C.inc("matches_limit_reached", 0);
	C.inc("unmapped_reads", 0);

	params_t params;
	if(!prsArgs(argc, argv, &params)) {
		dsHlp();
		return 1;
	}

	unordered_map<hash_t, char> bLstmers;
	if (!params.bLstFl.empty()) {
		bLstmers = readBlstKmers(params.bLstFl);
	}

	T.start("indexing");
	T.start("index_reading");
	cerr << "Reading index " << params.tFile << "..." << endl;
	string ref_name;
	pos_t T_sz;
	Sketch t;
	read_fasta_klib(params.tFile, [&t, &params, &bLstmers, &ref_name, &T_sz, &T](kseq_t *seq) {
		assert(ref_name.empty());  // TODO: support multi-sequence files
		T.start("index_sketching");
		t = buildFMHSketch(seq->seq.s, params.k, params.hFrac, bLstmers);
		T.stop("index_sketching");
		T_sz = (pos_t)seq->seq.l;
		ref_name = seq->name.s;
	});
	T.stop("index_reading");

	T.start("index_initializing");
	SketchIndex tidx(t, T_sz, ref_name);
	T.stop("index_initializing");

	T.stop("indexing");
	C.inc("T_sz", tidx.T_sz);

	if (!params.paramsFile.empty()) {
		cerr << "Writing parameters to " << params.paramsFile << "..." << endl;
		auto fout = ofstream(params.paramsFile);
		params.print(fout, false);
	} else {
		params.print(cerr, true);
	}

	SweepMap mapper(tidx, params, &T, &C);
	mapper.map(params.pFile, bLstmers);

	T.stop("total");
	mapper.print_report(C, T);

	return 0;
}
