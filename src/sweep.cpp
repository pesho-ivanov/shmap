#include "sweep.h"

int main(int argc, char **argv) {
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
	SketchIndex tidx(params.tFile, params.k, params.w, bLstmers);
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
