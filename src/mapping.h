#pragma once

#include <climits>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <memory>
#include "cmath"

#include "utils.h"
#include "buckets.h"
#include "types.h"

namespace sweepmap {

// the minimal mapping info for PAF output
struct MappingPAF {
    // local information: depending on the specific mapping
    qpos_t P_start;
    qpos_t P_end;
    char strand;         // '+' or '-'
    const char *segm_name;   // name of the reference segment
    rpos_t segm_sz;     // number of nucleotides in the segment
    int segm_id;
    rpos_t T_l;         // the position of the leftmost nucleotide of the mapping
    rpos_t T_r;         // the position of the rightmost nucleotide of the mapping
    int mapq;           // from 0 to 60

    // global info: depending on the read
    const char* query_id;   // the read name
    qpos_t P_sz;          // pattern size |P| bp 

    MappingPAF() : P_start(-1), P_end(-1), strand('?'), segm_name(""), segm_sz(-1), segm_id(-1), T_l(-1), T_r(-1), mapq(255), query_id(""), P_sz(-1) {}

    MappingPAF(qpos_t P_start, qpos_t P_end, char strand, const char *segm_name, rpos_t segm_sz, int segm_id, rpos_t T_l, rpos_t T_r)
        : P_start(P_start), P_end(P_end), strand(strand), segm_name(segm_name), segm_sz(segm_sz), segm_id(segm_id), T_l(T_l), T_r(T_r), mapq(255), query_id(""), P_sz(-1) {}

    friend std::ostream& operator<<(std::ostream& os, const MappingPAF& paf) {
        if (paf.P_sz == 0) throw std::runtime_error("MappingPAF is not set at time of output.");
        os << std::fixed << std::setprecision(3);
        os  << paf.query_id             // query sequence name
            << "\t" << paf.P_sz         // query sequence length
            << "\t" << paf.P_start      // query start (0-based; closed)
            << "\t" << paf.P_end        // query end (0-based; open)
            << "\t" << paf.strand       // '+' or '-'
            << "\t" << paf.segm_name    // reference segment name
            << "\t" << paf.segm_sz      // T_sz -- target sequence length
            << "\t" << paf.T_l          // target start on original strand (0-based)
            << "\t" << paf.T_r          // target end on original strand (0-based)
            << "\t" << paf.P_sz         // TODO: fix; Number of residue matches (number of nucleotide matches)
            << "\t" << paf.P_sz         // TODO: fix; Alignment block length: total number of sequence matches, mismatches, insertions and deletions in the alignment
            << "\t" << paf.mapq;        // mapping quality (0-255; 255 for missing)
        return os;
    }
};

struct LocalMappingStats {
    rpos_t s_sz;      // the position of the rightmost nucleotide of the mapping
    qpos_t intersection;     // the number of kmers in the intersection between the pattern and its mapping in `t'
    double map_time;
    int same_strand_seeds;  // diff between positive or negative matches
    bool unreasonable;  // reserved for filtering matches
    std::vector<Match>::const_iterator l, r;
    double J, J2;     // Similarity in [0;1] for the best and for the second best mapping
    qpos_t p_sz;     // number of seeds (subset of the sketch kmers)
    Buckets::Bucket bucket;          // the bucket where the mapping is found
    Buckets::Bucket bucket2;  // the bucket of the second best mapping
    qpos_t intersection2; // number of matches in the second best mapping
    double sigmas_diff;  // how many sigmas is the diff between intersection1 and intersection2

    LocalMappingStats() {
        s_sz = -1;
        intersection = -1;
        map_time = -1.0;
        same_strand_seeds = -1;
        unreasonable = false;
        J = -1.0;
        J2 = -1.0;
        p_sz = -1;
        intersection2 = -1;
        sigmas_diff = -1.0;
    }

    friend std::ostream& operator<<(std::ostream& os, const LocalMappingStats& stats) {
        // per mapping stats
        os  << "\t" << "p:i:"           << stats.p_sz  // sketches
            << "\t" << "s:i:"           << stats.s_sz
            << "\t" << "I:i:"           << stats.intersection  // intersection of `p` and `s` [kmers]
            << "\t" << "I2:i:"          << stats.intersection2
            << "\t" << "Idiff:i:"       << stats.intersection - stats.intersection2
            << "\t" << "Isdiff:f:"      << stats.sigmas_diff
            << "\t" << "J:f:"           << stats.J   // Similarity [0; 1]
            << "\t" << "J2:f:"          << stats.J2   // second best mapping similarity [0; 1]
            << "\t" << "Jdiff:f:"       << stats.J - stats.J2   // second best mapping similarity [0; 1]
            << "\t" << "strand:i:"      << stats.same_strand_seeds
            << "\t" << "t:f:"           << stats.map_time
            << "\t" << "b:s:"           << stats.bucket.segm_id << "," << stats.bucket.b
            << "\t" << "b2:s:"          << stats.bucket2.segm_id << "," << stats.bucket2.b;
        return os;
    }
};

struct GlobalMappingStats {
    // internal stats
    int k;     // kmer size
    qpos_t seeds;               // number of seeds before pruning
    rpos_t max_seed_matches;// number of matches of the most frequent seed
    rpos_t seed_matches;    // number of matches of seeds
    rpos_t total_matches;   // number of matches of all kmers
    double match_inefficiency;   // intersection / total_matches
    rpos_t seeded_buckets;     // the initial number of buckets
    rpos_t final_buckets;   // the number of buckets after pruning
    double FPTP;         // false discovery rate: FP / PP
    double gt_J, gt_C, gt_C_bucket;  // Jaccard and Containment index of the ground truth mapping

    GlobalMappingStats() {
        k = -1;
        seeds = -1;
        total_matches = -1;
        max_seed_matches = -1;
        seed_matches = -1;
        seeded_buckets = -1;
        final_buckets = -1;
        FPTP = -1.0;
        match_inefficiency = -1.0;
        gt_J = -1.0;
        gt_C = -1.0;
        gt_C_bucket = -1.0;
    }

    friend std::ostream& operator<<(std::ostream& os, const GlobalMappingStats& stats) {
        os << std::fixed << std::setprecision(3);
        os  << "\t" << "k:i:"           << stats.k
            << "\t" << "gt_J:f:"        << stats.gt_J
            << "\t" << "gt_C:f:"        << stats.gt_C
            << "\t" << "gt_C_bucket:f:" << stats.gt_C_bucket
            << "\t" << "Seeds:i:"       << stats.seeds
            << "\t" << "MSeedmax:i:"    << stats.max_seed_matches
            << "\t" << "MSeed:i:"       << stats.seed_matches
            << "\t" << "M:i:"           << stats.total_matches
            << "\t" << "Mineff:f:"      << stats.match_inefficiency
            << "\t" << "Bmax:i:"        << stats.seeded_buckets
            << "\t" << "Bfinal:i:"      << stats.final_buckets
            << "\t" << "FPTP:f:"        << stats.FPTP;
        return os;
    }
};

class Mapping {
public:
    MappingPAF paf;
    LocalMappingStats local_stats;
    std::unique_ptr<GlobalMappingStats> global_stats;

    Mapping() {}
    Mapping(const Mapping &m) : paf(m.paf), local_stats(m.local_stats), global_stats(m.global_stats ? std::make_unique<GlobalMappingStats>(*m.global_stats) : nullptr) {}

    void update(qpos_t P_start, qpos_t P_end, rpos_t T_l, rpos_t T_r, const RefSegment &segm, qpos_t intersection, double new_score,
                int same_strand_seeds, std::vector<Match>::const_iterator l, std::vector<Match>::const_iterator r) {
        if (new_score > score()) {
            paf = MappingPAF(P_start, P_end, same_strand_seeds > 0 ? '+' : '-', segm.name.c_str(), segm.sz, segm.id, T_l, T_r);
            local_stats.s_sz = r->hit.tpos - l->hit.tpos + 1;
            local_stats.same_strand_seeds = same_strand_seeds;
            local_stats.intersection = intersection;
            local_stats.J = new_score;
            local_stats.l = l;
            local_stats.r = r;
        }
    }

    double score() const { return local_stats.J; }
    double score2() const { return local_stats.J2; }
    int segm_id() const { return paf.segm_id; }

    void set_score2(double score2) {
        local_stats.J2 = score2;
    }

    qpos_t intersection() const {
        return local_stats.intersection;
    }

    int calc_mapq(double theta, double best_score_delta) {
        assert(score() >= 0.0);
        if (score2() < theta)
            set_score2(theta);
        if (score() - score2() > best_score_delta) {
            if (abs(local_stats.same_strand_seeds) < local_stats.intersection/2)
                return 5;
            return 60;
        } else {
            return 0;
        }
    }

    double sigmas_diff(qpos_t X, qpos_t Y) {
        if (X==-1) X=0;
        if (Y==-1) Y=0;
        return std::abs(X - Y) / std::sqrt(X + Y);
    }

    void set_mapq(double theta, double best_score_delta) {
        paf.mapq = calc_mapq(theta, best_score_delta);
    }

    double mapq() const {
        return paf.mapq;
    }

    Buckets::Bucket bucket() const {
        return local_stats.bucket;
    }

    void set_bucket(const Buckets::Bucket &bucket) {
        local_stats.bucket = bucket;
    }

    void set_second_best(const Mapping &m2) {
        local_stats.bucket2 = m2.local_stats.bucket;
        local_stats.intersection2 = m2.local_stats.intersection;

        set_score2(m2.score());
        local_stats.bucket2 = m2.local_stats.bucket;
        local_stats.intersection2 = m2.local_stats.intersection;
        local_stats.sigmas_diff = sigmas_diff(local_stats.intersection2, local_stats.intersection);
    }

    void set_global_stats(const Counters &C, const char* query_id, qpos_t P_sz, int k, qpos_t seeds, double FPTP, const std::string &segm_name, rpos_t segm_sz, double map_time) {
        paf.query_id = query_id;
        paf.P_sz = P_sz;
        local_stats.map_time = map_time;

        global_stats = std::make_unique<GlobalMappingStats>();
        global_stats->k = k;
        global_stats->seeds = seeds;
        global_stats->total_matches = C["total_matches"];
        global_stats->max_seed_matches = C["max_seed_matches"];
        global_stats->seed_matches = C["seed_matches"];
        global_stats->seeded_buckets = C["seeded_buckets"];
        global_stats->final_buckets = C["final_buckets"];
        global_stats->FPTP = FPTP;
        global_stats->gt_J = -1.0;
        global_stats->gt_C = -1.0;
        global_stats->gt_C_bucket = -1.0;
        global_stats->match_inefficiency = 1.0 * C["total_matches"] / local_stats.intersection;
    }

    static bool do_overlap(const Mapping &a, const Mapping &b) {
        if (a.paf.segm_id != b.paf.segm_id)
            return false;
        return a.paf.T_r >= b.paf.T_l && a.paf.T_l <= b.paf.T_r;
    }

    static double overlap(const Mapping &a, const Mapping &b) {
        if (a.paf.segm_id != b.paf.segm_id)
            return 0.0;
        int cap = std::max(0, std::min(a.paf.T_r, b.paf.T_r) - std::max(a.paf.T_l, b.paf.T_l));
        int cup = std::max(a.paf.T_r, b.paf.T_r) - std::min(a.paf.T_l, b.paf.T_l);
        assert(cup >= 0 && cap >=0 && cup >= cap);
        return 1.0 * cap / cup;
    }

    void print_paf() const {
        std::cout << *this;
    }

    static std::string reverseComplement(const std::string& seq) {
        std::string revComp;
        revComp.reserve(seq.size());
        for (auto i = rpos_t(seq.size()) - 1; i >= 0; --i) {
            switch (seq[i]) {
                case 'A':
                    revComp.push_back('T');
                    break;
                case 'C':
                    revComp.push_back('G');
                    break;
                case 'G':
                    revComp.push_back('C');
                    break;
                case 'T':
                    revComp.push_back('A');
                    break;
                default:
                    revComp.push_back('N');
                    break;
            }
        }
        return revComp;
    }

    friend std::ostream& operator<<(std::ostream& os, const Mapping& mapping) {
        os << mapping.paf;
        os << mapping.local_stats;
        if (mapping.global_stats != nullptr) os << *mapping.global_stats;
        os << std::endl;
        return os;
    }
};

} // namespace sweepmap
