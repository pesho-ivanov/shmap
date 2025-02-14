#include "index.h"
#include "mapper.h"
#include "handler.h"

using namespace sweepmap;

int main(int argc, char **argv) {
	params_t params;
	if(!params.prsArgs(argc, argv))
		return 1;

	Handler H(params);
	SketchIndex tidx(&H);
	tidx.build_index(H.params.tFile);
	Mapper* mapper = MapperFactory::createMapper(tidx, &H);
	if (mapper == nullptr) return 1;
	mapper->map_reads(H.params.pFile);

	return 0;
}
