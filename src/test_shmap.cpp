#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "../ext/doctest.h"
#include <thread>
#include <chrono>
#include <functional>

#include "types.h"
#include "shmap.h"
#include "io.h"  // For PaulFile and pos_t definitions

using namespace sweepmap;
using namespace doctest;

TEST_SUITE("Counting") {
    TEST_CASE("Counter") {
        Counter c;
        CHECK(c.count() == 0);
        c.inc(5);
        CHECK(c.count() == 5);
    }

    TEST_CASE("Counter operator[]") {
        Counter c;
        CHECK(c[0] == 0);  // Initial value
        
        c[0] = 5;  // Write through operator[]
        CHECK(c[0] == 5);  // Read through operator[]
        CHECK(c.count() == 5);  // Verify count() sees the change
        
        c.inc(3);  // Modify through inc()
        CHECK(c[0] == 8);  // Verify operator[] sees the change
    }

    TEST_CASE("Counter operator+=") {
        Counter c;
        CHECK(c.count() == 0);

        // Basic addition
        c += 5;
        CHECK(c.count() == 5);

        // Chaining
        (c += 3) += 2;
        CHECK(c.count() == 10);

        // Interaction with inc()
        c.inc(5);
        CHECK(c.count() == 15);
        c += 5;
        CHECK(c.count() == 20);

        // Return value correctness
        Counter& ref = (c += 1);
        CHECK(&ref == &c);  // Should return reference to self
        CHECK(c.count() == 21);

        // Counter += Counter
        Counter c2;
        c2.inc(5);
        c += c2;
        CHECK(c.count() == 26);  // 21 + 5

        // Chaining with Counter
        Counter c3;
        c3.inc(3);
        (c += c2) += c3;
        CHECK(c.count() == 34);  // 26 + 5 + 3

        // Original counters unchanged
        CHECK(c2.count() == 5);
        CHECK(c3.count() == 3);
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

    TEST_CASE("Counters operator[]") {
        Counters C;
        
        // Auto-creation and initialization
        CHECK(C["new_counter"][0] == 0);
        
        // Write through operator[]
        C["counter1"][0] = 5;
        CHECK(C["counter1"][0] == 5);
        CHECK(C.count("counter1") == 5);
        
        // Modify through inc() and verify through operator[]
        C.inc("counter1", 3);
        CHECK(C["counter1"][0] == 8);
        
        // Multiple counters
        C["counter2"][0] = 10;
        CHECK(C["counter2"][0] == 10);
        CHECK(C.frac("counter1", "counter2") == Approx(0.8));
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
        for (int i = 0; i < (int)sk_s.size(); i++) {
            CHECK(sk_s[i].r == i+k-1);
            if (i < (int)sk_s_rc.size()) {
                CHECK(sk_s[i].r == (int)sk_s.size() - sk_s_rc[i].r + k + 1);
                CHECK(sk_s[i].h == sk_s_rc[i].h);
            }
        }
    }
}

TEST_CASE("Indexing") {
    params_t params;
    params.k = 4;
    params.hFrac = 1.0;
    Handler H(params);
    SketchIndex tidx(&H);

    SUBCASE("Single segment") {
        string ref_file = "ref.fa";
        string T = "ACCAGTACCA";
        vector<int> cnt = {2, 1, 1, 1, 1, 1, 2};
    
        std::ofstream refFile(ref_file);
        refFile << ">ref\n" << T << std::endl;
        refFile.close();
        tidx.build_index(ref_file);
        std::remove(ref_file.c_str());

        sketch_t t = H.sketcher.sketch(T);
        REQUIRE(t.size() == 7);
        for (int i=0; i<(int)t.size(); i++)
            CHECK(tidx.count(t[i].h) == cnt[i]);
        CHECK(tidx.H->C.count("indexed_hits") == 7);
        CHECK(tidx.H->C.count("indexed_kmers") == 6);
    }

    SUBCASE("Two segments") {
        string ref_file = "ref2.fa";
        string segm1 = "ACCAGTACCA";
        string segm2 = "GGACCA";
        vector<int> cnt = {3, 1, 1, 1, 1, 1, 3};
    
        std::ofstream refFile(ref_file);
        refFile << ">segm1\n" << segm1 << std::endl;
        refFile << ">segm2\n" << segm2 << std::endl;
        refFile.close();
        tidx.build_index(ref_file);
        std::remove(ref_file.c_str());

        sketch_t t1 = H.sketcher.sketch(segm1);
        REQUIRE(t1.size() == 7);
        for (int i=0; i<(int)t1.size(); i++)
            CHECK(tidx.count(t1[i].h) == cnt[i]);
        CHECK(tidx.H->C.count("indexed_hits") == 10);
        CHECK(tidx.H->C.count("indexed_kmers") == 8);
    }
}

TEST_CASE("SHSingleReadMapper::select_kmers selects and sorts kmers correctly") {
    // Setup test dependencies
    params_t params;
    params.k = 4;
    params.hFrac = 1.0;
    params.theta = 0.7;
    Handler H(params);
    SketchIndex tidx(&H);
    
    string query_id = "test_read";
    string sequence = "ACGTACGTACGT";
    ofstream paulout;
    Matcher matcher(tidx);
    
    SHSingleReadMapper mapper(tidx, &H, matcher, query_id, sequence, paulout);
    
    // Create test sketch
    sketch_t p;
    p.emplace_back(0, 0x1234, false);  // pos, hash, rc
    p.emplace_back(4, 0x1234, false);  // duplicate hash
    p.emplace_back(2, 0x5678, false);
    p.emplace_back(6, 0x9ABC, false);
    
    // Call select_kmers
    auto [kmers, m] = mapper.select_elements(p);
    
    // Verify results
    REQUIRE(m == 4); // Total kmers including duplicates
    REQUIRE(kmers.size() == 3); // Unique kmers
    
    // Verify kmers are sorted by increasing number of hits in tidx
    for (size_t i = 1; i < kmers.size(); i++) {
        REQUIRE(kmers[i-1].hits_in_T <= kmers[i].hits_in_T);
    }
    
    // Verify duplicate kmer handling
    auto first_kmer = std::find_if(kmers.begin(), kmers.end(),
        [](const Seed& s) { return s.el.h == 0x1234; });
    REQUIRE(first_kmer != kmers.end());
    REQUIRE(first_kmer->occs_in_p == 2); // Should count duplicates
}

TEST_CASE("Bucketing") {
    int l = 10;
    Buckets B(l);
    CHECK(B.get_bucket_len() == 10);
    segm_t segm_id = 0;
    SUBCASE("Bucket") {
        Buckets::Bucket b0(segm_id, 0, &B);
        CHECK(b0.begin() == 0);
        CHECK(b0.end() == 2*l);
        Buckets::Bucket b1(segm_id, 1, &B);
        CHECK(b1.begin() == l);
        CHECK(b1.end() == l+2*l);
        CHECK(!(b0 == b1));
    }
    SUBCASE("Add matches to bucket") {
        B.add_to_bucket(Buckets::Bucket(segm_id, 2, &B), 10);
        CHECK(B.get_matches(Buckets::Bucket(segm_id, 2, &B)) == 10);
    }
    SUBCASE("Add matches to position") {
        B.add_to_pos(Buckets::Pos(segm_id, 5), 1);
        CHECK(B.size() == 1);
        CHECK(B.get_matches(Buckets::Bucket(segm_id, 0, &B)) == 1);
        B.add_to_pos(Buckets::Pos(0, 15), 1);
        CHECK(B.size() == 2);
        CHECK(B.get_matches(Buckets::Bucket(segm_id, 0, &B)) == 2);
        CHECK(B.get_matches(Buckets::Bucket(segm_id, 1, &B)) == 1);
    }
    SUBCASE("Iteration") {
        vector<int> matches = {1, 2, 3, 4, 5, 6, 5, 4, 2, 2};
        for (int i = 0; i < (int)matches.size(); i++)
            B.add_to_bucket(Buckets::Bucket(segm_id, i, &B), matches[i]);
        int cnt = 0;
        SUBCASE("Unordered") {
            for (auto it = B.unordered_begin(); it != B.unordered_end(); ++it) {
                REQUIRE(it->first.b >= 0);
                REQUIRE(it->first.b < matches.size());
                CHECK(it->second > 0);
                CHECK(matches[it->first.b] == it->second);
                matches[it->first.b] = -1;
                ++cnt;
            }
        }
        SUBCASE("Ordered") {
            int prev_matches = std::numeric_limits<int>::max();
            for (auto it = B.ordered_begin(); it != B.ordered_end(); ++it) {
                REQUIRE(it->first.b >= 0);
                REQUIRE(it->first.b < matches.size());
                CHECK(it->second > 0);
                CHECK(matches[it->first.b] == it->second);
                matches[it->first.b] = -1;
                ++cnt;
                // specific to ordered iteration
                CHECK(prev_matches >= it->second);
                prev_matches = it->second;
            }
        }
        CHECK(cnt == matches.size());
    }
}

//TEST_CASE("choose_seeds") {
//    Seeds kmers;
//    kmers.emplace_back(Seed(El(10, 111111, true), 1, 1, 0));
//    kmers.emplace_back(Seed(El{20, 222222, true}, 1, 2, 0)); 
//    kmers.emplace_back(Seed(El{30, 333333, true}, 1, 3, 0));
//    qpos_t m = 3;
//    double theta = 0.5;
//    auto [required, remaining, seeds_with_repeats] = SHSingleReadMapper::choose_seeds(kmers, m, theta);
//    CHECK(required.size() == 2);
//    CHECK(remaining.size() == 1);
//    CHECK(seeds_with_repeats == 0);
//}

TEST_CASE("Mapping a toy read") {
	params_t params;
	params.k = 25;
	params.hFrac = 0.05;
	params.theta = 0.7;

	Handler H(params);
	SketchIndex tidx(&H);
	
	// Define the missing variables
	string query_id = "NONE";
	string P = "ACCATCG";
	std::ofstream paulout;

	Matcher matcher(tidx);

	SHSingleReadMapper mapper(tidx, &H, matcher, query_id, P, paulout);

	sketch_t p;
	p.emplace_back(10, 111111, false);
	p.emplace_back(20, 222222, false);
	p.emplace_back(30, 111111, false);
	p.emplace_back(40, 444444, false);
	p.emplace_back(50, 555555, false);

	auto [kmers, m] = mapper.select_elements(p);

    CHECK(m == 5);

	INFO("Seeding");
	CHECK_MESSAGE(kmers.size() == 4, "wrong number of read kmers");
	for (int i=1; i<(int)kmers.size(); i++) {
		CHECK_MESSAGE(kmers[i-1].hits_in_T <= kmers[i].hits_in_T,
			"read kmers[", i-1, "].hits_in_T=", kmers[i-1].hits_in_T,
			" and read kmers[", i, "].hits_in_T=", kmers[i].hits_in_T,
			" not ordered by increasing hits");
    }
}
