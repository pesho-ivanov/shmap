#pragma once

#include <string>
#include "index.h"
#include "handler.h"

namespace sweepmap {

class Mapper {
public:
	virtual void map_reads(const string &) = 0;
    virtual void print_stats() = 0;
};

class MapperFactory {
public:
    static Mapper* createMapper(const SketchIndex& tidx, Handler* H);
};

} // namespace sweepmap