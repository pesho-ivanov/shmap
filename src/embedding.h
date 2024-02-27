#include <random>

#include "sketch.h"

namespace sweepmap {

using std::vector;

int D = 50;
std::uniform_real_distribution<float> distr(-1.0, 1.0);

class Embedding {
	vector<float> v;

public:
    Embedding() {
        v.resize(D, 0.0);
    }

	Embedding(const Kmer &kmer) {
        std::minstd_rand gen(kmer.h);

		v.reserve(D);
        for (int i = 0; i < D; ++i)
            v.push_back(distr(gen));
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
};

}