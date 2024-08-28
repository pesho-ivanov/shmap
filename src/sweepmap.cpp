#include "index.h"
#include "sweepmap.h"
#include "handler.h"

using namespace sweepmap;

int main(int argc, char **argv) {
	params_t params;
	if(!prsArgs(argc, argv, &params)) {
		dsHlp();
		exit(1);
	}

	Handler H(params);

	SketchIndex tidx(H.params, &H.C, &H.T, H.sketcher);
	tidx.build_index(H.params.tFile);

	if (!H.params.paramsFile.empty()) {
		cerr << "Writing parameters to " << H.params.paramsFile << "..." << endl;
		auto fout = std::ofstream(H.params.paramsFile);
		H.params.print(fout, false);
	} else {
		H.params.print(cerr, true);
	}

	cerr << "Mapping reads " << H.params.pFile << "..." << endl;
	SweepMap sweepmap(tidx, H.params, &H.C, &H.T, H.sketcher);
	sweepmap.map(H.params.pFile);

	return 0;
}
