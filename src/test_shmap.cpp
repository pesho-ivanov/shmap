#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../ext/doctest.h"

//#include "shmap.h"
//#include "handler.h"

//using namespace sweepmap;

//    CHECK(factorial(2) == 2);
//    CHECK(factorial(3) == 6);
//    CHECK(factorial(10) == 3628800);

//#include <gtest/gtest.h>
//#include "shmap.h"
//#include "handler.h"
//
//namespace sweepmap {
//
//class SHMapperTest : public ::testing::Test {
//protected:
//    Handler* H;
//    SketchIndex* tidx;
//    SHMapper* mapper;
//
//    void SetUp() override {
//        params_t params;
//        params.k = 25;
//        params.hFrac = 0.05;
//        params.theta = 0.7;
//
//        H = new Handler(params);
//        tidx = new SketchIndex(H);
//        mapper = new SHMapper(*tidx, H);
//    }
//
//    void TearDown() override {
//        delete mapper;
//        delete tidx;
//        delete H;
//    }
//};
//
//} // namespace sweepmap