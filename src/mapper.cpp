#include <iostream>

#include "mapper.h"
#include "sweepmap.h"
//#include "bucketmap.h"
//#include "rmqmap.h"
#include "jaccmap.h"

namespace sweepmap {

Mapper* MapperFactory::createMapper(const std::string& type, const SketchIndex& tidx, Handler* H) {
    if (type == "sweep") {
        return new SweepMapper(tidx, H);
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
    } else if (type == "jacc") {
        //if (H->params.max_matches != -1) {
        //    cerr << "Bucket mapper does not need limiting the number of matches. Do not set max_matches or set to -1." << endl;
        //    return nullptr;
        //}
        if (H->params.max_seeds != -1) {
            cerr << "JaccMapper does not need limiting the number of seeds. Do not set max_seeds or set to -1." << endl;
            return nullptr;
        }
        return new JaccMapper(tidx, H);
    } else {
		cerr << "Mapper '" << type << "' is not supported. Choose 'sweep' or 'bucket'." << endl;
        return nullptr;
    }
}

}