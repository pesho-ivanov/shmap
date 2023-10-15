CC = g++
CXX_STANDARD = -std=c++2a
CFLAGS = -g -O3 -march=native -lm -lpthread -Wno-unused-result
LIBS = minimap2/libminimap2.a -lz
DEPFLAGS = -MMD -MP

SRCS = src/sweep.cpp src/IO.h src/Index.h src/Measures.h src/Sketch.h
SWEEP_BIN = ./sweep

all: sweep

$(SWEEP_BIN): $(SRCS)
	$(CC) $(CXX_STANDARD) $(CFLAGS) $< -o $@ $(LIBS)

run100: sweep
	mkdir -p ../out
	$(SWEEP_BIN) -s simulations/genomes/t2thumanChrY.fasta -p simulations/reads/100.fasta -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt >out/sweep100.out
	$(SWEEP_BIN) -s simulations/genomes/t2thumanChrY.fasta -p simulations/reads/100.fasta -k 15 >out/sweep100.out
	$(SWEEP_BIN) -s simulations/genomes/t2thumanChrY.fasta -p simulations/reads/100.fasta -k 15 -a extend_equally >out/sweep100.out

clean:
	rm -f sweep