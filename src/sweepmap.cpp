#include "index.h"
#include "mapper.h"
#include "handler.h"

using namespace sweepmap;

int main(int argc, char **argv) {
	params_t params;
	if(!params.prsArgs(argc, argv)) {
		params.dsHlp();
		exit(1);
	}

	Handler H(params);
	SketchIndex tidx(&H);
	tidx.build_index(H.params.tFile);
	Mapper* mapper = MapperFactory::createMapper("sweepmap", tidx, &H);
	mapper->map(H.params.pFile);

	return 0;
}
