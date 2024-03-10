# distutils: language=c++

cimport cython
from cython.operator cimport dereference as deref

from libc cimport stdint
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set
from libcpp.bit cimport rotl, rotr, countr_zero
from libcpp.algorithm cimport nth_element
from libcpp.utility cimport pair

ctypedef stdint.uint64_t hash_t
ctypedef stdint.uint32_t pos_t
ctypedef stdint.int8_t segm_t

cdef struct Kmer:
    pos_t pos
    hash_t hash
    bool strand

ctypedef vector[Kmer] sketch_t

# TODO: Meaningful names
cdef hash_t[256] LUT_fw, LUT_rc

# TODO: What the hell is LUT
cdef initialize_lut():
    # https://gist.github.com/Daniel-Liu-c0deb0t/7078ebca04569068f15507aa856be6e8
    LUT_fw[<unsigned char>b'a'] = LUT_fw[<unsigned char>b'A'] = 0x3c8bfbb395c60474
    LUT_fw[<unsigned char>b'c'] = LUT_fw[<unsigned char>b'C'] = 0x3193c18562a02b4c
    LUT_fw[<unsigned char>b'g'] = LUT_fw[<unsigned char>b'G'] = 0x20323ed082572324
    LUT_fw[<unsigned char>b't'] = LUT_fw[<unsigned char>b'T'] = 0x295549f54be24456
    LUT_rc[<unsigned char>b'a'] = LUT_rc[<unsigned char>b'A'] = LUT_fw[<unsigned char>b'T']
    LUT_rc[<unsigned char>b'c'] = LUT_rc[<unsigned char>b'C'] = LUT_fw[<unsigned char>b'G']
    LUT_rc[<unsigned char>b'g'] = LUT_rc[<unsigned char>b'G'] = LUT_fw[<unsigned char>b'C']
    LUT_rc[<unsigned char>b't'] = LUT_rc[<unsigned char>b'T'] = LUT_fw[<unsigned char>b'A']

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef (sketch_t, pos_t) build_fmh_sketch(const unsigned char[:] s, pos_t k, double h_frac):
    initialize_lut()

    cdef size_t length = s.shape[0]
    cdef sketch_t kmers
    kmers.reserve(<size_t>(1.1 * <float>length * <float>h_frac))

    if length < k:
        return kmers, 0

    cdef hash_t h, h_fw = 0, h_rc = 0
    cdef hash_t h_thresh = <hash_t>(h_frac * <double>stdint.UINT64_MAX)

    cdef pos_t r
    for r in range(k):
        h_fw ^= rotl(LUT_fw[s[r]], k-r-1)
        h_rc ^= rotl(LUT_rc[s[r]], r)

    r = k
    while True:
        # HACK! the lowest differing bit is not expected to correlate much with (h < hThres)
        # TODO: tie break

        first_diff_bit = 1 << countr_zero(h_fw ^ h_rc)
        strand = h_fw & first_diff_bit
        h = h_rc if strand else h_fw

        if h < h_thresh:
            kmers.push_back(
                Kmer(pos=r-k+1, hash=h, strand=strand)
            )

        if r >= length:
            break

        h_fw = rotl(h_fw, 1) ^ rotl(LUT_fw[s[r-k]], k) ^ LUT_fw[s[r]]
        h_rc = rotr(h_rc, 1) ^ rotr(LUT_rc[s[r-k]], 1) ^ rotl(LUT_rc[s[r]], k-1)

        r += 1

    return kmers, r-k+1

cdef class Sketch:
    cdef sketch_t _sketch
    cdef readonly pos_t size

    def __init__(self, const unsigned char[:] s, pos_t k, double h_frac):
        self._sketch, self.size = build_fmh_sketch(s, k, h_frac)

    def state(self):
        return self._sketch


cdef struct RefPos:
    segm_t segm
    pos_t pos


cdef class Index:
    cdef vector[pos_t] segment_size
    cdef unordered_map[hash_t, RefPos] hit_first
    cdef unordered_map[hash_t, vector[RefPos]] hit_rest

    cdef _add(self, sketch_t sketch, pos_t size, pos_t max_matches=0):
        cdef size_t segm_id = self.segment_size.size()
        self.segment_size.push_back(size)
        for i in range(sketch.size()):
            kmer = sketch[i]
            hit = RefPos(segm_id, kmer.pos)
            if not self.hit_first.contains(kmer.hash):
                self.hit_first[kmer.hash] = hit
            else:
                self.hit_rest[kmer.hash].push_back(hit)

        # Apply max_matches
        cdef vector[hash_t] to_remove
        cdef hash_t h
        if max_matches > 0:
            for (h, _) in self.hit_rest:
                if self.hit_rest[h].size() > max_matches:
                    to_remove.push_back(h)
            for i in range(to_remove.size()):
                h = to_remove[i]
                self.hit_first.erase(h)
                self.hit_rest.erase(h)

    cdef pos_t hits_count(self, hash_t hash):
        if self.hit_first.contains(hash):
            return 1 + self.hit_rest[hash].size()
        else:
            return 1 if self.hit_first.contains(hash) else 0

    def add(self, Sketch sketch, max_matches=None):
        self._add(sketch._sketch, sketch.size, max_matches=max_matches if isinstance(max_matches, int) else 0)

    def state(self):
        return (self.hit_first, self.hit_rest)


cdef struct Mapping:
    segm_t segm
    pos_t pos

cdef bool compare_hits((hash_t, pos_t) a, (hash_t, pos_t) b):
    if a[1] == b[1]:
        # TODO: Add random order for equal hits
        # Using kmer might be a bad idea
        return a[0] < b[0]
    else:
        return a[1] < b[1]

cdef bool compare_matches((hash_t, RefPos) a, (hash_t, RefPos) b):
    if a[1].segm == b[1].segm:
        return a[1].pos < b[1].pos
    else:
        return a[1].segm < b[1].segm

cdef extern from "pdqsort.h" nogil:
    void pdqsort[Iter, Compare](Iter begin, Iter end, Compare comp) except +
    void pdqsort_branchless[Iter, Compare](Iter begin, Iter end, Compare comp) except +


cdef class DiffHist:
    """
    Intersection of two histograms. 

    H1
    value: +1       I += min(max(-D, 0), value)
        D>=0 => 0
        D<0  => I+1
    value: -1       I += min(0, value + max(0, D))
        D>0  => 0     I += min(0, value + D)
        D<=0 => I-1   I += value
    H2
    value: +1       I += min(max(D, 0), value)
        D>0  => I+1   
        D<=0 => 0
    value: -1       I += min(0, value + max(0, -D))
        D>=0  => I-1  I += value
        D<0   => 0    I += min(0, value - D)
    """
    cdef unordered_map[hash_t, int] diff  # Difference between corresponding bins H1[i]-H2[i]
    cdef readonly int intersection

    def state(self):
        return (self.intersection, self.diff)

    def __cinit__(self):
        self.intersection = 0

    cdef update(self, hash_t key, int value, int hist):
        cdef unordered_map[hash_t, int].iterator iter = self.diff.insert(pair[hash_t, int](key, 0)).first
        cdef int d = deref(iter).second

        if value > 0:
            self.intersection += min(max((-1)*hist*d, 0), value)
        else:
            self.intersection += min(0, value + max(0, hist*d))
        
        deref(iter).second = d + hist*value

    cdef update_first(self, hash_t key, int value):
        self.update(key, value, 1)

    cdef update_second(self, hash_t key, int value):
        self.update(key, value, -1)

    cdef reserve(self, size_t size):
        self.diff.reserve(size)


cpdef sweep_map(Index reference_index, Sketch pattern_sketch_, pos_t max_hashes, pos_t width):
    cdef sketch_t pattern_sketch = pattern_sketch_._sketch
    cdef vector[Mapping] mapping

    # Build a limited set of Kmer hashes from the pattern.
    # We pick a random subset of hashes with the lowest hits in the reference index.
    cdef vector[(hash_t, pos_t)] hash_with_hits          # (hash, hits in reference)
    cdef unordered_map[hash_t, pos_t] hash_pattern_hits  # whether a hash has been added to hash_with_hits
    hash_with_hits.reserve(pattern_sketch.size())
    hash_pattern_hits.reserve(pattern_sketch.size())

    cdef unordered_map[hash_t, pos_t].iterator iter
    cdef pos_t hits

    for pos in range(pattern_sketch.size()):
        hash = pattern_sketch[pos].hash
        hits = reference_index.hits_count(hash)
        if hits > 0:
            iter = hash_pattern_hits.insert(pair[hash_t, pos_t](hash, 0)).first
            if deref(iter).second == 0:
                hash_with_hits.push_back((hash, hits))
            deref(iter).second = deref(iter).second + 1

    cdef size_t total_hashes = hash_with_hits.size()
    if total_hashes > max_hashes:
        total_hashes = max_hashes
    
    # This makes sure that matched_hashes[:total_hashes] has elements with the lowest hits
    nth_element(
        hash_with_hits.begin(),
        hash_with_hits.begin() + total_hashes,
        hash_with_hits.end(),
        compare_hits,
    )
    
    cdef size_t total_matches = 0
    for i in range(total_hashes):
        total_matches += hash_with_hits[i][1]

    # Estimate the total number of mathes to reserve space for them
    cdef unordered_set[hash_t] selected_hashes
    for i in range(total_hashes):
        selected_hashes.insert(hash_with_hits[i][0])

    # Build the list of matched kmers in the reference sorted by position
    cdef vector[(hash_t, RefPos)] reference_matches  # (hash, ref_pos)
    reference_matches.reserve(total_matches)
    for hash in selected_hashes:
        # TODO: remove unnecessary lookups
        if reference_index.hit_first.contains(hash):
            reference_matches.push_back((hash, reference_index.hit_first[hash]))
            for match in reference_index.hit_rest[hash]:
                reference_matches.push_back((hash, match))

    pdqsort_branchless(reference_matches.begin(), reference_matches.end(), compare_matches)

    # Initialize diff-histogram with Pattern
    cdef DiffHist hist = DiffHist()
    hist.reserve(total_hashes)
    for i in range(total_hashes):
        h = hash_with_hits[i][0]
        hist.update_first(h, hash_pattern_hits[h])

    # index in reference_matches
    cdef pos_t left = 0
    cdef pos_t right = 0 
    
    cdef pos_t best_intersection = 0
    cdef pos_t best_left = 0
    #all_mappings = []

    for left in range(reference_matches.size()):
        
        # Extend the right end of the window
        while (
            right < reference_matches.size() and
            reference_matches[left][1].segm == reference_matches[right][1].segm and
            reference_matches[left][1].pos + width > reference_matches[right][1].pos
        ):
            # Add [right] to the window
            hist.update_second(reference_matches[right][0], +1)

            right += 1

        if reference_matches[left][1].pos + width < reference_index.segment_size[reference_matches[left][1].segm]:
            # Process the window
            if hist.intersection > best_intersection:
                best_intersection = hist.intersection
                best_left = left
            #all_mappings.append((reference_matches[left][1].pos, hist.intersection))

        # Remove [left] from the window
        hist.update_second(reference_matches[left][0], -1)

    return reference_matches[best_left][1].pos, best_intersection
    #return (reference_matches[best_left][1].pos, best_intersection), all_mappings
    # return mapping
