#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
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
        //C.init("c1", "c2", "c3", "c4");
        CHECK_THROWS(C.count("c1"));
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
        //CHECK_THROWS(t.secs());
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

//TEST_CASE("Parameters in handler") {
//	Handler* H;
//
//	params_t params;
//	params.k = 25;
//	params.hFrac = 0.05;
//	params.theta = 0.7;
//
//	H = new Handler(params);
//}

TEST_CASE("Indexing") {
	Handler* H;
	SketchIndex* tidx;

	params_t params;
	params.k = 4;
	params.hFrac = 1.0; 
	H = new Handler(params);
	tidx = new SketchIndex(H);

    SUBCASE("Single segment") {
        string ref_file = "ref.fa";
        string T = "ACCAGTACCA";
        vector<int> cnt = {2, 1, 1, 1, 1, 1, 2};
    
        std::ofstream refFile(ref_file);
        refFile << ">ref\n" << T << std::endl;
        refFile.close();
        tidx->build_index(ref_file);
        std::remove(ref_file.c_str());

        sketch_t t = H->sketcher.sketch(T);
        REQUIRE(t.size() == 7);
        for (int i=0; i<(int)t.size(); i++)
            CHECK(tidx->count(t[i].h) == cnt[i]);
        CHECK(tidx->H->C.count("indexed_hits") == 7);
        CHECK(tidx->H->C.count("indexed_kmers") == 6);
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
        tidx->build_index(ref_file);
        std::remove(ref_file.c_str());

        sketch_t t1 = H->sketcher.sketch(segm1);
        REQUIRE(t1.size() == 7);
        for (int i=0; i<(int)t1.size(); i++)
            CHECK(tidx->count(t1[i].h) == cnt[i]);
        CHECK(tidx->H->C.count("indexed_hits") == 10);
        CHECK(tidx->H->C.count("indexed_kmers") == 8);
    }

    SUBCASE("Get matches") {
//    Matches matches;
//    Kmer kmer(0, 0, false);
//    rpos_t hits_in_t = 1; 
//    Seed el(kmer, hits_in_t, kmers.size());
//    tidx->get_matches(&matches, seed);
    }
}

//TEST_CASE("Bucketing") {
//    typedef Buckets<false> Buckets;
//    int l = 10;
//	Buckets B(l);
//    CHECK(B.get_bucket_len() == 10);
//    segm_t segm_id = 0;
//    SUBCASE("Bucket") {
//        BucketLoc b0(segm_id, 0);
//        CHECK(B.begin(b0) == 0);
//        CHECK(B.end(b0) == 2*l);
//        BucketLoc b1(segm_id, 1);
//        CHECK(B.begin(b1) == l);
//        CHECK(B.end(b1) == l+2*l);
//        CHECK(!(b0 == b1));
//    }
//    SUBCASE("Add matches to bucket") {
//        B.add_to_bucket(BucketLoc(segm_id, 2), 10);
//        CHECK(B.get_matches(BucketLoc(segm_id, 2)) == 10);
//    }
//    SUBCASE("Add matches to position") {
//        B.add_to_pos(Pos(segm_id, 5), 1);
//        CHECK(B.size() == 1);
//        CHECK(B.get_matches(BucketLoc(segm_id, 0, &B)) == 1);
//        B.add_to_pos(Pos(0, 15), 1);
//        CHECK(B.size() == 2);
//        CHECK(B.get_matches(BucketLoc(segm_id, 0, &B)) == 2);
//        CHECK(B.get_matches(BucketLoc(segm_id, 1, &B)) == 1);
//    }
//    SUBCASE("Iteration") {
//        vector<int> matches = {1, 2, 3, 4, 5, 6, 5, 4, 2, 2};
//        for (int i = 0; i < (int)matches.size(); i++)
//            B.add_to_bucket(BucketLoc(segm_id, i, &B), matches[i]);
//        int cnt = 0;
//        SUBCASE("Unordered") {
//            for (auto it = B.unordered_begin(); it != B.unordered_end(); ++it) {
//                REQUIRE(it->first.b >= 0);
//                REQUIRE(it->first.b < matches.size());
//                CHECK(it->second > 0);
//                CHECK(matches[it->first.b] == it->second);
//                matches[it->first.b] = -1;
//                ++cnt;
//            }
//        }
//        SUBCASE("Ordered") {
//            int prev_matches = std::numeric_limits<int>::max();
//            for (auto it = B.ordered_begin(); it != B.ordered_end(); ++it) {
//                REQUIRE(it->first.b >= 0);
//                REQUIRE(it->first.b < matches.size());
//                CHECK(it->second > 0);
//                CHECK(matches[it->first.b] == it->second);
//                matches[it->first.b] = -1;
//                ++cnt;
//                // specific to ordered iteration
//                CHECK(prev_matches >= it->second);
//                prev_matches = it->second;
//            }
//        }
//        CHECK(cnt == matches.size());
//    }
//}

//TEST_SUITE("Metrics") {
//    qpos_t intersection = 5;
//    qpos_t m = 10;
//    qpos_t s_sz = 15;
//    TEST_CASE("Jaccard") {
//        double ans = 1.0*intersection / (m + s_sz - intersection);
//	    CHECK(SHMapper::Jaccard(intersection, m, s_sz) == Approx(ans));
//    }
//    TEST_CASE("Containment index") {
//        double ans = 1.0*intersection / m;
//	    CHECK(SHMapper::ContainmentIndex(intersection, m) == Approx(ans));
//    }
//}

//TEST_CASE("Mapping a toy read") {
//	Handler* H;
//	SketchIndex* tidx;
//	SHMapper* mapper;
//
//	params_t params;
//	params.k = 25;
//	params.hFrac = 0.05;
//	params.theta = 0.7;
//
//	H = new Handler(params);
//	tidx = new SketchIndex(H);
//	mapper = new SHMapper(*tidx, H);
//
//	sketch_t p;
//	p.emplace_back(10, 111111, false);
//	p.emplace_back(20, 222222, false);
//	p.emplace_back(30, 111111, false);
//	p.emplace_back(40, 444444, false);
//	p.emplace_back(50, 555555, false);
//
//	int nonzero = 0;
//	Seeds kmers = mapper->select_kmers(p, nonzero);
//
//	INFO("Seeding");
//	CHECK_MESSAGE(kmers.size() == 4, "wrong number of read kmers");
//	for (int i=1; i<(int)kmers.size(); i++) {
//		CHECK_MESSAGE(kmers[i-1].hits_in_T <= kmers[i].hits_in_T,
//			"read kmers[", i-1, "].hits_in_T=", kmers[i-1].hits_in_T,
//			" and read kmers[", i, "].hits_in_T=", kmers[i].hits_in_T,
//			" not ordered by increasing hits");
//	}
//	return;
//	
//	CHECK(kmers[0].kmer.h == 1);
//	CHECK(kmers[0].occs_in_p == 2);
//	
//	CHECK(kmers[1].kmer.h == 2);
//	CHECK(kmers[1].occs_in_p == 1);
//	
//	CHECK(kmers[2].kmer.h == 3);
//	CHECK(kmers[2].occs_in_p == 2);
//
//	for (size_t i = 1; i < kmers.size(); i++)
//		CHECK(kmers[i-1].hits_in_T == kmers[i].hits_in_T);
//}

TEST_CASE("test_select_kmers_pmatches") {
	params_t params;
	params.k = 25;
	params.hFrac = 0.05;
	params.theta = 0.7;

	auto H = new Handler(params);
	auto tidx = new SketchIndex(H);
	auto mapper = new SHMapper<false, false, false>(*tidx, H);

    sketch_t p;
    p.emplace_back(60, 0x111111, false);
    p.emplace_back(70, 0x222222, false);
    p.emplace_back(10, 0x111111, true);
    p.emplace_back(20, 0x222222, true);
    p.emplace_back(30, 0x111111, true);
    p.emplace_back(40, 0x444444, false);
    p.emplace_back(50, 0x555555, false);

    set< vector<qpos_t> > pmatches_gt;
    pmatches_gt.insert({60, 30, 10});
    pmatches_gt.insert({70, 20});
    pmatches_gt.insert({40});
    pmatches_gt.insert({50});

    set< vector<qpos_t> > pmatches_res;
    Seeds S = mapper->unique_elements_with_info(p);
    for (auto s : S) {
        //for (auto m : s.pmatches)
        //    std::cout << m << " ";
        //std::cout << std::endl;
        pmatches_res.insert(s.pmatches);
    }

    CHECK(pmatches_res.size() == pmatches_gt.size());
    CHECK(pmatches_res == pmatches_gt);
}

TEST_CASE("lcs") {
    sketch_t t;
    t.emplace_back(10, 0x111111, true);
    t.emplace_back(20, 0x222222, true);
    t.emplace_back(30, 0x111111, true);
    t.emplace_back(40, 0x444444, false);
    t.emplace_back(50, 0x555555, false);
    t.emplace_back(60, 0x111111, false);
    t.emplace_back(70, 0x222222, false);

    params_t params; // Create dummy params for Handler
    Handler H(params);
    SketchIndex tidx(&H);
    // Add a minimal segment to make buckets work
    tidx.add_segment("test", "AAAAAAAAAA", t);
    Buckets<false> B(tidx);
    B.set_halflen(2);
    segm_t segm_id = 0;
    rpos_t block = 1;  // t[2..6) including kmers (30, 0x111111), (40, 0x444444), (50, 0x555555), (60, 0x111111)
                       // result:                      1,              3,              4,              5
    const BucketLoc bucket(segm_id, block);
    h2seed_t p_ht;
    // Seed(kmer, hits_in_T, occs_in_p, seed_num, pmatches)
    p_ht.insert(std::make_pair(hash_t(0x111111), Seed(Kmer(1, 0x111111, false), 99,  3, 0, vector<qpos_t>{5, 1})));
    p_ht.insert(std::make_pair(hash_t(0x222222), Seed(Kmer(2, 0x222222, false), 999, 2, 1, vector<qpos_t>{4, 2})));
    p_ht.insert(std::make_pair(hash_t(0x444444), Seed(Kmer(3, 0x444444, false), 9,   1, 2, vector<qpos_t>{3})));
    // => seq = {5,1,3,5,1}  -- kmers from p in bucket t[2..6)
    // => lis = {1,3,5}

    params.k = 25;
    params.hFrac = 0.05;
    params.theta = 0.7;
    auto mapper = new SHMapper<false, false, false>(tidx, &H);
    qpos_t lcs_cnt = mapper->lcs(t, B, bucket, p_ht);
    CHECK(lcs_cnt == 3);
}

std::string exec(const std::string& cmd) {
    std::array<char, 128> buffer;
    std::string result;
    #ifdef _WIN32
      FILE* pipe = _popen(cmd.c_str(), "r");
    #else
      FILE* pipe = popen(cmd.c_str(), "r");
    #endif
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
        result += buffer.data();
    }
    #ifdef _WIN32
      _pclose(pipe);
    #else
      pclose(pipe);
    #endif
    return result;
}

TEST_CASE("test_execute_sweepmap_from_terminal") {
    // Skip integration test - requires release binary and test files
    return;
    
    std::string make_res = exec("make");
//    std::cout << make_res << std::endl;
//    CHECK(make_res.empty());
    std::string tool = "./release/shmap";
    std::string cmd = tool + " -p test.fa -t test.fa -k 25 -hFrac 0.05 -theta 0.7 -metric bucket_LCS -o test.sam";
    std::string result = exec(cmd);
    std::cout << result << std::endl;
//    CHECK(result.empty());
}

