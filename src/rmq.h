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
        return std::max(l,r);
    }

    int incRange(int rs, int re, int i, int s, int e) {
        if(rs>e || s>re) return 0;
        if(rs<=s && e<=re) return ++T[i];
        int m = (s+e)/2;
        incRange(rs, re, 2*i,   s,   m);
        incRange(rs, re, 2*i+1, m+1, e);
        return T[i] = std::max(T[2*i], T[2*i+1]);
    }

public:
    int P_sz;
    SegmentTree(int n) : n(n), T(4*n+1) {}

    inline int size() const { return n; }
    inline int from(const sweepmap::Hit &hit) const { return std::max(hit.r-P_sz, 0); }
    inline int to(const sweepmap::Hit &hit) const  { return std::min(hit.r+P_sz, size()-1); }

    inline int query(int qs, int qe) const {
        return query(qs, qe, 1, 0, n-1);
    }

    inline int query(const sweepmap::Hit &hit) const {
        return query(from(hit), to(hit));
    }

    inline int incRange(int rs, int re) {
        return incRange(rs, re, 1, 0, n-1);
    }

    inline int incRange(const sweepmap::Hit &hit) {
        return incRange(from(hit), to(hit));
    }
};