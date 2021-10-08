
all:
	g++ -ggdb -std=c++11 -c -o dna_map.o dna_map.cpp
	g++ -ggdb -std=c++11 -c -o bins.o bins.cpp
	g++ -ggdb -std=c++11 -c -o main.o main.cpp
	g++ -ggdb -std=c++11 -c -o chromosome.o chromosome.cpp
	g++ -ggdb -std=c++11 -c -o dna.o dna.cpp
	g++ -ggdb -std=c++11 -o drv dna_map.o main.o bins.o chromosome.o dna.o

clean:
	rm -rf drv dna_map.o main.o bins.o chromosome.o
