clusters: main.o clustering_algorithm.o
	g++ -o $@ $^

main.o: main.cpp

clustering_algorithm.o: clustering_algorithm.cpp
