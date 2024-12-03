#include "index.h"
#include "mapper.h"
#include "handler.h"

using namespace sweepmap;

int main(int argc, char **argv) {
	params_t params;
	if(!params.prsArgs(argc, argv)) {
		params.dsHlp();
		return 1;
	}

	Handler H(params);
	SketchIndex tidx(&H);
	tidx.build_index(H.params.tFile);
	Mapper* mapper = MapperFactory::createMapper(H.params.mapper, tidx, &H);
	if (mapper == nullptr) return 1;
	mapper->map_all_reads(H.params.pFile);
	H.print_sketching_stats();

	return 0;
}
