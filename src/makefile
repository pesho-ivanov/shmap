CC = g++
CXX_STANDARD = -std=c++2a
CFLAGS = -g -O3 -march=native -lm -lpthread -Wno-unused-result
LIBS = ../minimap2/libminimap2.a -lz
DEPFLAGS = -MMD -MP

# Define the source files
SRC_FILES := sweep.cpp

# Generate corresponding dependency file names
DEP_FILES := $(SRC_FILES:.cpp=.d)

# Include the dependency files if they exist
-include $(DEP_FILES)

all: sweep

sweep: sweep.cpp
	$(CC) $(CXX_STANDARD) $(CFLAGS) $(DEPFLAGS) -c $<
	$(CC) $(CXX_STANDARD) $< -o $@ $(LIBS)
#	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

run100:
	mkdir -p ../out
	./sweep -s ../simulations/genomes/t2thumanChrY.fasta -p ../simulations/reads/100.fasta -k 15 >../out/sweep100.out 2>../out/sweep100.er

clean:
	rm -f sweep $(DEP_FILES)