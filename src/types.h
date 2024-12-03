#pragma once

#include <climits>
#include <iostream>
#include <string>
#include <vector>

#include "utils.h"

namespace sweepmap {

// Kmer -- a kmer with metadata (a position in the sequence, hash, strand)
struct Kmer {
    rpos_t r;      // kmer resides [l, r), where l+k=r
    hash_t h;
    bool strand;  // false: forward, true: reverse
    Kmer(rpos_t r, hash_t h, bool strand) : r(r), h(h), strand(strand) {}
    bool operator==(const Kmer &other) const { return h == other.h; }
    friend std::ostream& operator<<(std::ostream& os, const Kmer& kmer) {
        os << "Kmer(r=" << kmer.r << ", h=" << kmer.h << ", strand=" << kmer.strand << ")";
        return os;
    }
};

using sketch_t = std::vector<Kmer>;

// Reference segment structure
struct RefSegment {
    sketch_t kmers;
    std::string name;
    std::string seq;   // empty if only mapping and no alignment
    rpos_t sz;
    int id;
    RefSegment(const sketch_t &sk, const std::string &name, const std::string &seq, rpos_t sz, int id)
        : kmers(sk), name(name), seq(seq), sz(sz), id(id) {}
};

// Seed -- a kmer with metadata (a position in the query P and number of hits in the reference T)
struct Seed {
    Kmer kmer;
    rpos_t hits_in_T;   // number of hits in the reference sketch
    qpos_t occs_in_p;   // number of occurrences in the query sketch
    qpos_t seed_num;    // seed number in the query sketch (unique)
    Seed(const Kmer &kmer, rpos_t hits_in_T, qpos_t occs_in_p, qpos_t seed_num) :
        kmer(kmer), hits_in_T(hits_in_T), occs_in_p(occs_in_p), seed_num(seed_num) {}    
    friend std::ostream& operator<<(std::ostream& os, const Seed& seed) {
        os << "Seed(" << seed.kmer << ", hits_in_T=" << seed.hits_in_T << ", occs_in_p=" << seed.occs_in_p << ", seed_num=" << seed.seed_num << ")";
        return os;
    }
};
using Seeds = std::vector<Seed>;

// Hit -- a kmer hit in the reference T
struct Hit {  // TODO: compress into a 32bit field
    rpos_t r;            // right end of the kmer [l, r), where l+k=r
    rpos_t tpos;         // position in the reference sketch
    bool strand;
    segm_t segm_id;
    Hit() {}
    Hit(const Kmer &kmer, rpos_t tpos, segm_t segm_id)
        : r(kmer.r), tpos(tpos), strand(kmer.strand), segm_id(segm_id) {}
    friend std::ostream& operator<<(std::ostream& os, const Hit& hit) {
        os << "Hit(r=" << hit.r << ", tpos=" << hit.tpos << ", strand=" << hit.strand << ", segm_id=" << hit.segm_id << ")";
        return os;
    }
};

// Match -- a pair of a seed and a hit
struct Match {
    Seed seed;
    Hit hit;
    Match(const Seed &seed, const Hit &hit)
        : seed(seed), hit(hit) {}
    
    inline bool is_same_strand() const {
        return seed.kmer.strand == hit.strand;
    }
    friend std::ostream& operator<<(std::ostream& os, const Match& match) {
        os << "Match(" << match.seed << ", " << match.hit << ")";
        return os;
    }
};
using Matches = std::vector<Match>;

} // namespace sweepmap
