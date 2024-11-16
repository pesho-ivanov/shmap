#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../ext/doctest.h"

#include "shmap.h"

using namespace sweepmap;

TEST_CASE("Counters") {
    Counters C;
    CHECK_THROWS_AS(C.count("c1"), std::out_of_range);
    C.inc("c1");
    CHECK(C.count("c1") == 1);
    C.inc("c2", 2);
    CHECK(C.count("c2") == 2);
    CHECK(C.frac("c1", "c2") == doctest::Approx(0.5));
    CHECK(C.perc("c1", "c2") == doctest::Approx(50.0));
}

TEST_CASE("Timers") {
}

TEST_CASE("FracMinHash sketching") {
    int k = 25;
    double hFrac = 0.05;
    Counters C;
    Timers T;
    FracMinHash sketcher(k, hFrac, &C, &T);

    SUBCASE("Sketching a simple sequence") {
        std::string s    = "ACGGT";
        std::string s_rc = "ACGT";
        sketch_t sk_s    = sketcher.sketch(s);
        sketch_t sk_s_rc = sketcher.sketch(s_rc);
        for (size_t i = 0; i < sk_s.size(); i++) {
            CHECK(sk_s[i].r == i);
            CHECK(sk_s[i].h == sk_s_rc[i].h);
            CHECK(sk_s[i].strand == false);
        }
        CHECK(sk_s == sk_s_rc);
//        int C_sketch_seqs = C.count("sketched_seqs");
//        CHECK(C_sketch_seqs == 1);
    }
}

TEST_CASE("Parameters in handler") {
	Handler* H;

	params_t params;
	params.k = 25;
	params.hFrac = 0.05;
	params.theta = 0.7;

	H = new Handler(params);
}

TEST_CASE("Indexing a toy sequence") {
	Handler* H;
	SketchIndex* tidx;

	params_t params;
	params.k = 25;
	params.hFrac = 0.05;
	params.theta = 0.7;

	H = new Handler(params);
	tidx = new SketchIndex(H);
}

TEST_CASE("Mapping a toy read") {
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

	sketch_t p;
	p.emplace_back(10, 111111, false);
	p.emplace_back(20, 222222, false);
	p.emplace_back(30, 111111, false);
	p.emplace_back(40, 444444, false);
	p.emplace_back(50, 555555, false);

	int nonzero = 0;
	SHMapper::Seeds kmers = mapper->select_kmers(p, nonzero);

	INFO("Seeding");
	CHECK_MESSAGE(kmers.size() == 4, "wrong number of read kmers");
	for (int i=0; i<(int)kmers.size(); i++) {
		CHECK_MESSAGE(kmers[i-1].hits_in_T <= kmers[i].hits_in_T,
			"read kmers[", i-1, "].hits_in_T=", kmers[i-1].hits_in_T,
			" and read kmers[", i, "].hits_in_T=", kmers[i].hits_in_T,
			" not ordered by increasing hits");
	}
	return;
	
	CHECK(kmers[0].kmer.h == 1);
	CHECK(kmers[0].occs_in_p == 2);
	
	CHECK(kmers[1].kmer.h == 2);
	CHECK(kmers[1].occs_in_p == 1);
	
	CHECK(kmers[2].kmer.h == 3);
	CHECK(kmers[2].occs_in_p == 2);

	for (size_t i = 1; i < kmers.size(); i++)
		CHECK(kmers[i-1].hits_in_T == kmers[i].hits_in_T);
}