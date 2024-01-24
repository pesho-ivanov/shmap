#include "sweepmap.h"
using namespace sweepmap;

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
	params.print_display(std::cerr);

	T.start("indexing");
	cerr << "Indexing " << params.tFile << "..." << endl;
	SketchIndex tidx(params);

	T.start("index_reading");
	read_fasta_klib(params.tFile, [&tidx, &params, &T](kseq_t *seq) {
		T.stop("index_reading");
		T.start("index_sketching");
		Sketch t = buildFMHSketch(seq->seq.s, params.k, params.hFrac);
		T.stop("index_sketching");

		T.start("index_initializing");
		tidx.add_segment(t, seq->name.s, seq->seq.s); // TODO: don't pass seq.s
		T.stop("index_initializing");

		T.start("index_reading");
	});
	T.stop("index_reading");
	T.stop("indexing");

	tidx.print_kmer_hist();

	C.inc("T_sz", tidx.total_size);
	//C.inc("index_memory_MB", get_current_memory_MB());

	if (!params.paramsFile.empty()) {
		cerr << "Writing parameters to " << params.paramsFile << "..." << endl;
		auto fout = std::ofstream(params.paramsFile);
		params.print(fout, false);
	} else {
		params.print(cerr, true);
	}

	cerr << "Mapping reads " << params.pFile << "..." << endl;
	SweepMap sweepmap(tidx, params, &T, &C);
	sweepmap.map(params.pFile);

	T.stop("total");
	//C.inc("total_memory_MB", get_current_memory_MB());
	sweepmap.print_report();

	return 0;
}
