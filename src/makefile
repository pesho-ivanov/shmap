CC = g++
#CFLAGS = -std=c++2a -O3 -lz -lm -lpthread -Wno-unused-result
CFLAGS = -g -std=c++2a -O3 -march=native -lz -lm -lpthread -Wno-unused-result

astar:
	$(CC) astar.cpp ../minimap2/libminimap2.a $(CFLAGS) -o astar

sweep: sweep.cpp
	$(CC) sweep.cpp ../minimap2/libminimap2.a $(CFLAGS) -o sweep

run100:
	mkdir -p ../out
	./sweep -s ../simulations/genomes/t2thumanChrY.fasta -p ../simulations/reads/100.fasta -k 15 >../out/sweep100.out 2>../out/sweep100.er
