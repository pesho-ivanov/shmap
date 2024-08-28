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
	SketchIndex tidx(&H);
	tidx.build_index(H.params.tFile);
	SweepMap sweepmap(tidx, &H);
	sweepmap.map(H.params.pFile);

	return 0;
}
