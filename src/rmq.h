#pragma once

#include <bits/stdc++.h>
#include <vector>

#include "sketch.h"

class SegmentTree {
    int n;
    std::vector<int> T;

    int query(int qs, int qe, int i, int s, int e) const {
        if(qs>e || s>qe) return -1;
        if(qs<=s && e<=qe) return T[i];
        int m = (s+e)/2;
        int l = query(qs, qe, 2*i,   s,   m);
        int r = query(qs, qe, 2*i+1, m+1, e);
        assert(l != -1 || r != -1);
        return std::max(l, r);
    }

    int incRange(int rs, int re, int i, int s, int e, int val) {
        if(rs>e || s>re) return -1;
        if(rs<=s && e<=re) return T[i]=T[i]+val;
        int m = (s+e)/2;
        incRange(rs, re, 2*i,   s,   m, val);
        incRange(rs, re, 2*i+1, m+1, e, val);
        assert(T[2*i] != -1 || T[2*i+1] != -1);
        return T[i] = std::max(T[2*i], T[2*i+1]);
    }

    void clear(int i) {
        if (i < (int)T.size() && T[i]) {
            T[i] = 0;
            clear(2*i);
            clear(2*i+1);
        }
    }

public:
    int P_sz;
    SegmentTree(int n) : n(n), T(4*n+1) {}

    inline int size() const { return n; }
    inline int from(const sweepmap::Hit &hit) const { return std::max(hit.r-P_sz, 0); }
    inline int to(const sweepmap::Hit &hit) const  { return std::min(hit.r, size()-1); }

    inline void clear() {
        clear(1);
    }

    inline int query(int qs, int qe) const {
        return query(qs, qe, 1, 0, n-1);
    }

    inline int query(const sweepmap::Hit &hit) const {
        return query(from(hit), to(hit));
    }

    inline int incRange(int rs, int re) {
        return incRange(rs, re, 1, 0, n-1, 1);
    }

    inline int incRange(const sweepmap::Hit &hit) {
        return incRange(from(hit), to(hit));
    }
};