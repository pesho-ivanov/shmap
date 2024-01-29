#include "index.h"
#include "sweepmap.h"
using namespace sweepmap;

int main(int argc, char **argv) {
	initialize_LUT();
	
	Counters C;
	Timers T;

	T.start("total");

	params_t params;
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
	//C.inc("total_memory_MB", get_current_memory_MB());
	sweepmap.print_report();

	return 0;
}
