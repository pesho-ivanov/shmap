CC = g++
#CFLAGS = -std=c++2a -O3 -lz -lm -lpthread -Wno-unused-result
CFLAGS = -pg -std=c++2a -O3 -DNDEBUG -march=native -lz -lm -lpthread -Wno-unused-result

eskemap:
	$(CC) Eskemap.cpp ../minimap2/libminimap2.a $(CFLAGS) -o eskemap

astar:
	$(CC) astar.cpp ../minimap2/libminimap2.a $(CFLAGS) -o astar

sweep: sweep.cpp
	$(CC) sweep.cpp ../minimap2/libminimap2.a $(CFLAGS) -o sweep
