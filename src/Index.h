#ifndef INDEX_HPP
#define INDEX_HPP

#include <algorithm>
#include "Sketch.h"

#define INDEX_DEFAULT_DUMP_FILE "indexDump.idx"

//A compare function to sort elements in a hash position vector primarily by hash and secondarily by position in ascending order
inline const bool smHshnPos(const pair<uint64_t, uint32_t>& hpa, const pair<uint64_t, uint32_t>& hpb){
	return hpa.first != hpb.first ? hpa.first < hpb.first : hpa.second < hpb.second;
}

//A compare function to sort elements in a hash position vector by ascending positions only
inline const bool smPos(const pair<uint64_t, uint32_t>& hpa, const pair<uint64_t, uint32_t>& hpb){
	return hpa.second < hpb.second;
}

//This function generates the L array need for alpha-homology detection
const vector<pair<uint64_t, uint32_t>> genL(const unordered_map<uint64_t, uint32_t>& phshs, const mm_idx_t *tidx);

class Index {

public:
	//Index constructor
	Index(const Sketch& sk){
		uint32_t i = 0, currPos = 0, start, end;
		uint64_t lastHash;
		vector<pair<uint64_t, uint32_t>>::const_iterator e;

		//Reserve enough space for vector
		this->pos.reserve(sk.size());

		//Store hash position pairs in pos
		for(Sketch::const_iterator h = sk.begin(); h != sk.end(); ++h, ++i) this->pos.push_back(make_pair(*h, i));

		//Resize pos after filling
		this->pos.shrink_to_fit();
		//Sort elements
		sort(this->pos.begin(), this->pos.end(), smHshnPos);
		//Fill hash table
		e = this->pos.begin();
		lastHash = e->first;
		start = 0;
		end = 0;

		while(e != this->pos.end()){
			if(lastHash != e->first){
				this->pints[lastHash] = make_pair(start, end);
				start = currPos;
				end = currPos;
				lastHash = e->first;
			}

			++currPos;
			++end;
			++e;
		}

		this->pints[lastHash] = make_pair(start, end);
	}

protected:
		//A hash table storing occurrence intervals of k-mers in pos (hash -> [start, end))
		unordered_map<uint64_t, pair<uint32_t, uint32_t>> pints;
		//Vector of k-mer hash position pairs
		vector<pair<uint64_t, uint32_t>> pos;
};

//This function generates the L array need for alpha-homology detection
const vector<pair<uint32_t, uint32_t>> genL(const unordered_map<uint64_t, int32_t>& phshs, const mm_idx_t *tidx){
	assert(false);
	uint32_t i;

	unordered_map<uint64_t, uint32_t> hash2ord;
	uint32_t ord = 0;
	
	//Number of times a hash can be found inside the text sketch
	int nHits;
	//An iterator to iterate over the sketch
	unordered_map<uint64_t, int32_t>::const_iterator sI;
	//The L array
	vector<pair<uint32_t, uint32_t>> L;

	//Reserve some space for L
	L.reserve(P_MULTIPLICITY * phshs.size());

	//Query index for all hashes occurring in the read sketch
	for(sI = phshs.begin(); sI != phshs.end(); ++sI){
		//Query current hash in index
		const uint64_t *idx_p = mm_idx_get(tidx, sI->first, &nHits);

		if (hash2ord.find(sI->first) == hash2ord.end())
			hash2ord[sI->first] = ord++;

		const uint32_t kmer_num = hash2ord[sI->first];

		//Testing
		// if(sI->first == 49938615){
		// 	cout << "genL: Current hash is " << sI->first << endl;
		// 	cout << "genL: Number of occurrences in text according to index: " << nHits << endl;
		// }

		//Check if hash could be found
		if(nHits > 0){
			//Iterate over all occurrences
			for(i = 0; i < nHits; ++i){
				//Push (hash value, k-mer position in sketch) pair
				L.push_back(make_pair(kmer_num, ((uint32_t)(*idx_p))>>1));
				//Move to next occurrence
				idx_p++;
			}
		}
	}

	//Sort L by ascending positions in text
	sort(L.begin(), L.end(), smPos);

	return L;
}

#endif