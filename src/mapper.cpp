#include <iostream>

#include "mapper.h"
//#include "sweepmap.h"
//#include "bucketmap.h"
//#include "rmqmap.h"
#include "shmap.h"

namespace sweepmap {

Mapper* MapperFactory::createMapper(const std::string& type, const SketchIndex& tidx, Handler* H) {
//    if (type == "sweep") {
//        return new SweepMapper(tidx, H);
//    } else if (type == "bucket") {
//        //if (H->params.max_matches != -1) {
//        //    cerr << "Bucket mapper does not need limiting the number of matches. Do not set max_matches or set to -1." << endl;
//        //    return nullptr;
//        //}
//        if (H->params.max_seeds != -1) {
//            cerr << "Bucket mapper does not need limiting the number of seeds. Do not set max_seeds or set to -1." << endl;
//            return nullptr;
//        }
//        return new BucketMapper(tidx, H);
//    } else if (type == "rmq") {
//        //if (H->params.max_matches != -1) {
//        //    cerr << "Bucket mapper does not need limiting the number of matches. Do not set max_matches or set to -1." << endl;
//        //    return nullptr;
//        //}
//        if (H->params.max_seeds != -1) {
//            cerr << "Bucket mapper does not need limiting the number of seeds. Do not set max_seeds or set to -1." << endl;
//            return nullptr;
//        }
//        return new RMQMapper(tidx, H);
//    } else
    if (type == "shmap") {
        //if (H->params.max_matches != -1) {
        //    cerr << "Bucket mapper does not need limiting the number of matches. Do not set max_matches or set to -1." << endl;
        //    return nullptr;
        //}
        //if (H->params.max_seeds != -1) {
        //    cerr << "SHMapper does not need limiting the number of seeds. Do not set max_seeds or set to -1." << endl;
        //    return nullptr;
        //}
        //return new SHMapper<H->params.tmpl>(tidx, H);
        bool nbp = H->params.no_bucket_pruning;
        bool os  = H->params.one_sweep;
        bool sh  = H->params.only_sh;
        bool ap  = H->params.abs_pos;

        int index = (nbp ? 1 : 0)
                  + (os  ? 2 : 0)
                  + (sh  ? 4 : 0)
                  + (ap  ? 8 : 0);

        switch (index) {
            case 0: return new SHMapper<false, false, false, false>(tidx, H);
//            case 1: return new SHMapper< true, false, false, false>(tidx, H);
//            case 3: return new SHMapper< true,  true, false, false>(tidx, H);
            case 4: return new SHMapper<false, false,  true, false>(tidx, H);
//            case 5: return new SHMapper< true, false,  true, false>(tidx, H);
//            case 6: return new SHMapper<false,  true,  true, false>(tidx, H);
//            case 7: return new SHMapper< true,  true,  true, false>(tidx, H);
//            case 8: return new SHMapper<false, false,  true,  true>(tidx, H);
//            case 9: return new SHMapper< true, false,  true,  true>(tidx, H);
//            case 10: return new SHMapper<false,  true,  true,  true>(tidx, H);
//            case 11: return new SHMapper< true,  true,  true,  true>(tidx, H);
//            case 12: return new SHMapper<false, false,  true,  true>(tidx, H);
//            case 13: return new SHMapper< true, false,  true,  true>(tidx, H);
//            case 14: return new SHMapper<false,  true,  true,  true>(tidx, H);
//            case 15: return new SHMapper< true,  true,  true,  true>(tidx, H);
        }
        throw std::runtime_error("Invalid templated/precompiled index for SHMapper");
        return nullptr;  // should never happen
    } else {
		cerr << "Mapper '" << type << "' is not supported. Choose 'sweep' or 'bucket'." << endl;
        return nullptr;
    }
}

}