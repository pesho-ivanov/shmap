#include "mapper.h"
#include "sweepmap.h"
#include "bucketmap.h"

namespace sweepmap {

Mapper* MapperFactory::createMapper(const std::string& type, const SketchIndex& tidx, Handler* H) {
    if (type == "sweepmap") {
        return new SweepMapper(tidx, H);
    } else if (type == "bucketmap") {
        return new BucketMapper(tidx, H);
    } else {
        return nullptr;
    }
}

}