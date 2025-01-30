#include <iostream>

#include "mapper.h"
//#include "sweepmap.h"
//#include "bucketmap.h"
//#include "rmqmap.h"
#include "shmap.h"

namespace sweepmap {

Mapper* MapperFactory::createMapper(const SketchIndex& tidx, Handler* H) {
    bool nbp = H->params.no_bucket_pruning;
    bool os  = H->params.one_sweep;
    bool ap  = H->params.abs_pos;

    int index = (nbp ? 1 : 0)
                + (os  ? 2 : 0)
                + (ap  ? 4 : 0);

    switch (index) {
        case 0: return new SHMapper<false, false, false>(tidx, H);
//            case 1: return new SHMapper< true, false, false>(tidx, H);
//            case 2: return new SHMapper<false,  true, false>(tidx, H);
//            case 4: return new SHMapper<false, false,  true>(tidx, H);
//            case 5: return new SHMapper< true, false,  true>(tidx, H);
//            case 6: return new SHMapper<false,  true,  true>(tidx, H);
//            case 7: return new SHMapper< true,  true,  true>(tidx, H);
    }
    throw std::runtime_error("Invalid templated/precompiled index for SHMapper");
    return nullptr;  // should never happen
}

}