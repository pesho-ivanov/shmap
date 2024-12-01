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

} // namespace sweepmap
