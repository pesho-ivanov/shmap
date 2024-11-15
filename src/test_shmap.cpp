#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../ext/doctest.h"

#include "shmap.h"
#include "handler.h"

using namespace sweepmap;

	TEST_CASE("testing the seeding") {
		Handler* H;
		SketchIndex* tidx;
		SHMapper* mapper;

		params_t params;
		params.k = 25;
		params.hFrac = 0.05;
		params.theta = 0.7;

		H = new Handler(params);
		tidx = new SketchIndex(H);
		mapper = new SHMapper(*tidx, H);

		CHECK(1 == 1);

		sketch_t p;
		p.emplace_back(1, 0, false);
		p.emplace_back(1, 1, false);
		p.emplace_back(2, 2, false);
		p.emplace_back(3, 3, false);
		p.emplace_back(3, 4, false);

		int nonzero = 0;
		SHMapper::Seeds kmers = mapper->select_kmers(p, nonzero);

		// Verify results
		CHECK(kmers.size() == 3); // Should have 3 unique kmers
		
		// First kmer (h=1) should have 2 occurrences
		CHECK(kmers[0].kmer.h == 1);
		CHECK(kmers[0].occs_in_p == 2);
		
		// Second kmer (h=2) should have 1 occurrence
		CHECK(kmers[1].kmer.h == 2);
		CHECK(kmers[1].occs_in_p == 1);
		
		// Third kmer (h=3) should have 2 occurrences
		CHECK(kmers[2].kmer.h == 3);
		CHECK(kmers[2].occs_in_p == 2);

		// Verify kmers are sorted by hits_in_T
		for (size_t i = 1; i < kmers.size(); i++) {
			CHECK(kmers[i-1].hits_in_T == kmers[i].hits_in_T);
		}
	}

int factorial(int number) { return number <= 1 ? number : factorial(number - 1) * number; }

TEST_CASE("testing the factorial function") {
    CHECK(factorial(1) == 1);
    CHECK(factorial(2) == 2);
    CHECK(factorial(3) == 6);
    CHECK(factorial(10) == 3628800);
}

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

TEST_F(SHMapperTest, SelectKmers) {
}
//
//} // namespace sweepmap