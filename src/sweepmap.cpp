#include "index.h"
#include "sweepmap.h"

using namespace sweepmap;

void print_time_stats(Timers *T, Counters *C) {
	cerr << std::fixed << std::setprecision(1);
	cerr << "Time [sec]:           "             << setw(5) << right << T->secs("total")             << " (" << setw(4) << right << C->count("reads") / T->secs("total")      << " reads per sec)" << endl;
	cerr << " | Index:                 "         << setw(5) << right << T->secs("indexing")          << " (" << setw(4) << right << T->perc("indexing", "total")              << "\%)" << endl;
	cerr << " |  | loading:                "     << setw(5) << right << T->secs("index_reading")     << " (" << setw(4) << right << T->perc("index_reading", "indexing")      << "\%)" << endl;
	cerr << " |  | sketch:                 "     << setw(5) << right << T->secs("index_sketching")   << " (" << setw(4) << right << T->perc("index_sketching", "indexing")    << "\%)" << endl;
	cerr << " |  | initialize:             "     << setw(5) << right << T->secs("index_initializing")<< " (" << setw(4) << right << T->perc("index_initializing", "indexing") << "\%)" << endl;
	cerr << " | Map:                   "         << setw(5) << right << T->secs("mapping")           << " (" << setw(4) << right << T->perc("mapping", "total")               << "\%, " << setw(5) << right << T->range_ratio("query_mapping") << "x)" << endl;
	cerr << " |  | load queries:           "     << setw(5) << right << T->secs("query_reading")     << " (" << setw(4) << right << T->perc("query_reading", "mapping")       << "\%, " << setw(5) << right << T->range_ratio("query_reading") << "x)" << endl;
	cerr << " |  | sketch reads:           "     << setw(5) << right << T->secs("sketching")         << " (" << setw(4) << right << T->perc("sketching", "mapping")           << "\%, " << setw(5) << right << T->range_ratio("sketching") << "x)" << endl;
	cerr << " |  | seeding:                "     << setw(5) << right << T->secs("seeding")           << " (" << setw(4) << right << T->perc("seeding", "mapping")             << "\%, " << setw(5) << right << T->range_ratio("seeding") << "x)" << endl;
	cerr << " |  |  | collect seed info:       " << setw(5) << right << T->secs("collect_seed_info") << " (" << setw(4) << right << T->perc("collect_seed_info", "seeding")   << "\%, " << setw(5) << right << T->range_ratio("collect_seed_info") << "x)" << endl;
	cerr << " |  |  | thin sketch:             " << setw(5) << right << T->secs("thin_sketch")       << " (" << setw(4) << right << T->perc("thin_sketch", "seeding")         << "\%, " << setw(5) << right << T->range_ratio("thin_sketch") << "x)" << endl;
	cerr << " |  |  | sort seeds:              " << setw(5) << right << T->secs("sort_seeds")        << " (" << setw(4) << right << T->perc("sort_seeds", "seeding")          << "\%, " << setw(5) << right << T->range_ratio("sort_seeds") << "x)" << endl;
	cerr << " |  |  | unique seeds:            " << setw(5) << right << T->secs("unique_seeds")      << " (" << setw(4) << right << T->perc("unique_seeds", "seeding")        << "\%, " << setw(5) << right << T->range_ratio("unique_seeds") << "x)" << endl;
	cerr << " |  | matching seeds:         "     << setw(5) << right << T->secs("matching")          << " (" << setw(4) << right << T->perc("matching", "mapping")            << "\%, " << setw(5) << right << T->range_ratio("matching") << "x)" << endl;
	cerr << " |  |  | collect matches:         " << setw(5) << right << T->secs("collect_matches")   << " (" << setw(4) << right << T->perc("collect_matches", "matching")    << "\%, " << setw(5) << right << T->range_ratio("collect_matches") << "x)" << endl;
	cerr << " |  |  | sort matches:            " << setw(5) << right << T->secs("sort_matches")      << " (" << setw(4) << right << T->perc("sort_matches", "matching")       << "\%, " << setw(5) << right << T->range_ratio("sort_matches") << "x)" << endl;
	cerr << " |  | sweep:                  "     << setw(5) << right << T->secs("sweep")             << " (" << setw(4) << right << T->perc("sweep", "mapping")               << "\%, " << setw(5) << right << T->range_ratio("sweep") << "x)" << endl;
	cerr << " |  | post proc:              "     << setw(5) << right << T->secs("postproc")          << " (" << setw(4) << right << T->perc("postproc", "mapping")            << "\%, " << setw(5) << right << T->range_ratio("postproc") << "x)" << endl;
//		cerr << "Virtual memory [MB]:  "             << setw(5) << right << C->count("total_memory_MB")  << endl;
//		cerr << " | Index:                 "         << setw(5) << right << C->count("index_memory_MB") << " (" << setw(4) << right << C->perc("index_memory_MB", "total_memory_MB") << "\%)" << endl;
	printMemoryUsage();
}

int main(int argc, char **argv) {
	Counters C;
	Timers T;
	params_t params;

	Sketch::initialize_LUT();
	Sketch::params = &params;
	Sketch::T = &T;
	Sketch::C = &C;

	T.start("total");

	if(!prsArgs(argc, argv, &params)) {
		dsHlp();
		return 1;
	}
	params.print_display(std::cerr);

	SketchIndex tidx(params, &T, &C);
	tidx.index(params.tFile);

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
	Sketch::print_stats();
	print_time_stats(&T, &C);

	return 0;
}
