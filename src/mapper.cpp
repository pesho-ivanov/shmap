#include <iostream>

#include "mapper.h"
#include "sweepmap.h"
#include "bucketmap.h"

namespace sweepmap {

Mapper* MapperFactory::createMapper(const std::string& type, const SketchIndex& tidx, Handler* H) {
    if (type == "sweep") {
        return new SweepMapper(tidx, H);
    } else if (type == "bucket") {
        return new BucketMapper(tidx, H);
    } else {
		cerr << "Mapper '" << type << "' is not supported. Choose 'sweep' or 'bucket'." << endl;
        return nullptr;
    }
}

}