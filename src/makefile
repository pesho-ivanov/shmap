CC = g++
CXX_STANDARD = -std=c++2a
CFLAGS = -g -O3 -march=native -lm -lpthread -Wno-unused-result
LIBS = ../minimap2/libminimap2.a -lz
DEPFLAGS = -MMD -MP

SRC_FILES := sweep.cpp

all: sweep

sweep: sweep.cpp IO.h Index.h Measures.h Sketch.h
	$(CC) $(CXX_STANDARD) $(CFLAGS) $< -o $@ $(LIBS)

run100: sweep
	mkdir -p ../out
	./sweep -s ../simulations/genomes/t2thumanChrY.fasta -p ../simulations/reads/100.fasta -k 15 >../out/sweep100.out

clean:
	rm -f sweep