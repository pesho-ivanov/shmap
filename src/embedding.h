//#include <immintrin.h>

#include "sketch.h"

namespace sweepmap {

using std::vector;

int D = 8*8;

class FastRand {
public:
    FastRand(unsigned long long seed = 1) : val(seed) {}

    float nextFloat() {
        val = 6364136223846793005ULL * val + 1;
        return float((val >> 40) * 2.3283064365386963e-10 * 2.0 - 1.0);
    }

private:
    unsigned long long val;
};


class Embedding {
	vector<float> v;

public:
    Embedding() {
        v.resize(D, 0.0);
    }

	Embedding(const Kmer &kmer) {
        v.resize(D);
        FastRand randGen(kmer.h);

        for (int i = 0; i < D; i += 8)
            v[i] = randGen.nextFloat();
	}

    Embedding operator+=(const Embedding &other) {
        for (int i = 0; i < D; ++i)
            v[i] += other.v[i];
        return *this;
    }

    Embedding operator-=(const Embedding &other) {
        for (int i = 0; i < D; ++i)
            v[i] -= other.v[i];
        return *this;
    }

	float dot(const Embedding &other) const {
		float sum = 0.0;
		for (int i = 0; i < D; ++i)
			sum += v[i] * other.v[i];
		return sum;
	}

	float l2(const Embedding &other) const {
		float sum = 0.0;
		for (int i = 0; i < D; ++i) {
            float d = v[i] - other.v[i];
			sum += d*d;
        }
		return 1.0/sum;
	}
};

}