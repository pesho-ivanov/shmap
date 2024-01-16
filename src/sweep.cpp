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
	string ref;
	ref.reserve(35000000000);
	cerr << "Reading index " << params.tFile << "..." << endl;
	if (readFASTA(params.tFile, &ref)) {
		cerr << "ERROR: Reference file could not be read" << endl;
	}
	ifstream fStr(params.tFile);
	string line;
	getline(fStr, line);
	string name = line.substr(1);
	T.stop("index_reading");

	T.start("index_sketching");
	SketchIndex tidx(name, ref, params, bLstmers);
	T.stop("index_sketching");
	//SketchIndex tidx(params.tFile, params, bLstmers);
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
