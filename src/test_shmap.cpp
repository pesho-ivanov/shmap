#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "../ext/doctest.h"
#include <thread>
#include <chrono>

#include "shmap.h"

using namespace sweepmap;
using namespace doctest;

TEST_SUITE("Counting") {
    TEST_CASE("Counter") {
        Counter c;
        CHECK(c.count() == 0);
        c.inc(5);
        CHECK(c.count() == 5);
    }
    TEST_CASE("Counters") {
        Counters C;
        CHECK_THROWS_AS(C.count("c1"), std::out_of_range);
        C.inc("c1");
        CHECK(C.count("c1") == 1);
        C.inc("c2", 2);
        CHECK(C.count("c2") == 2);
        CHECK(C.frac("c1", "c2") == Approx(0.5));
        CHECK(C.perc("c1", "c2") == Approx(50.0));
        CHECK_THROWS(C.frac("c1", "c3"));
        CHECK_THROWS(C.perc("c1", "c3"));
        C.inc("c4", 0);
        CHECK(C.frac("c1", "c4") == std::numeric_limits<double>::infinity());
        CHECK(C.perc("c1", "c4") == std::numeric_limits<double>::infinity());
    }
}

TEST_SUITE("Timing") {
    TEST_CASE("Timer") {
        sweepmap::Timer t;
        CHECK(t.secs() == Approx(0.00));
        CHECK_THROWS(t.range_ratio());
        t.start();
        CHECK_THROWS(t.secs());
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        t.stop();
        CHECK(t.secs() == Approx(0.01).epsilon(0.001));
        CHECK(t.range_ratio() == Approx(1.0).epsilon(0.01));

        t.start();
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        t.stop();
        CHECK(t.range_ratio() == Approx(2.0).epsilon(0.1));
    }
    TEST_CASE("Timers") {
        Timers T;
        CHECK_THROWS(T.range_ratio("t1"));

        T.start("t1");
        T.start("t2");
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        T.stop("t1");
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        T.stop("t2");

        CHECK(T.secs("t1") == Approx(0.01).epsilon(0.001));
        CHECK(T.secs("t2") == Approx(0.02).epsilon(0.001));
        CHECK(T.perc("t1", "t2") == Approx(50.0).epsilon(1.0));
        CHECK(T.perc("t2", "t1") == Approx(200.0).epsilon(4.0));

        CHECK(T.range_ratio("t2") == Approx(1.0).epsilon(0.1));
        T.start("t2");
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        T.stop("t2");
        CHECK(T.range_ratio("t2") == Approx(2.0).epsilon(0.1));
    }
}

TEST_CASE("FracMinHash sketching") {
    int k = 4;
    double hFrac = 1.0;
    Counters C;
    Timers T;
    FracMinHash sketcher(k, hFrac, &C, &T);

    SUBCASE("Sketching a sequence shorter than k") {
        CHECK(sketcher.sketch("ACC").size() == 0);
    }

    SUBCASE("Sketching is the same for a reverse-complement of a sequence (except for coordinates)") {
        std::string s    = "ACGGT";
        std::string s_rc = "ACCGT";
        sketch_t sk_s    = sketcher.sketch(s);
        sketch_t sk_s_rc = sketcher.sketch(s_rc);
        std::reverse(sk_s_rc.begin(), sk_s_rc.end());

        REQUIRE_MESSAGE(sk_s.size() == sk_s_rc.size(), "the sketches of reverse-complement strings should have the same size");
        for (int i = 0; i < sk_s.size(); i++) {
            CHECK(sk_s[i].r == i+k-1);
            if (i < (int)sk_s_rc.size()) {
                CHECK(sk_s[i].r == (int)sk_s.size() - sk_s_rc[i].r + k + 1);
                CHECK(sk_s[i].h == sk_s_rc[i].h);
            }
        }
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
	for (int i=1; i<(int)kmers.size(); i++) {
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