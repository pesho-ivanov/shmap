#pragma once

#include <bits/stdc++.h>
#include <vector>

// This file is a derivative of a work under Attribution-ShareAlike 4.0 International
// https://cp-algorithms.com/data_structures/segment_tree.html#adding-on-segments-querying-for-maximum

#include "sketch.h"

// TODO: make sure to use [l,r) intervals instead of [l,r]
class SegmentTree {
    int n;
    std::vector<int> t, lazy;
    const int MINUS_INF = std::numeric_limits<int>::min();

    void push(int v) {
        t[v*2] += lazy[v];
        lazy[v*2] += lazy[v];
        t[v*2+1] += lazy[v];
        lazy[v*2+1] += lazy[v];
        lazy[v] = 0;
    }

    void update(int v, int tl, int tr, int l, int r, int addend) {
        if (l > r) 
            return;
        if (l == tl && tr == r) {
            t[v] += addend;
            lazy[v] += addend;
        } else {
            push(v);
            int tm = (tl + tr) / 2;
            update(v*2, tl, tm, l, min(r, tm), addend);
            update(v*2+1, tm+1, tr, max(l, tm+1), r, addend);
            t[v] = max(t[v*2], t[v*2+1]);
        }
    }

    int query(int v, int tl, int tr, int l, int r) {
        if (l > r)
            return MINUS_INF;
        if (l == tl && tr == r)
            return t[v];
        push(v);
        int tm = (tl + tr) / 2;
        return max(query(v*2, tl, tm, l, min(r, tm)), 
                query(v*2+1, tm+1, tr, max(l, tm+1), r));
    }

    void clear(int i) {
        if (i >= (int)t.size()) return;
        
        if (t[i] != 0 || lazy[i] != 0) {
            t[i] = 0;
            lazy[i] = 0;
            clear(2*i);
            clear(2*i+1);
        }
    }

public:
    int P_sz;
    SegmentTree(int n) : n(n), t(4*n+1), lazy(4*n+1) {}

    int size() const { return n; }
    int from(const sweepmap::Hit &hit) const { return std::max(hit.r-P_sz, 0); }
    int to(const sweepmap::Hit &hit) const  { assert(hit.r < size()); return std::min(hit.r, size()-1); }

    void clear() {
        clear(1);
    }

    int query(int l, int r) {  // [l,r)
        return query(1, 0, n-1, l, r-1);  // [l,r-1]
    }

    void incRange(int l, int r) {  // [l,r)
#ifdef DEBUG
        cerr << "incRange(" << l << "," << r << ")" << endl;
#endif
        update(1, 0, n-1, l, r-1, 1);  // [l,r-1]
    }
};