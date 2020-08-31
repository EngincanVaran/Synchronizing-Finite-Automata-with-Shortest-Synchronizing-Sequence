all: ams_lex.o data_structures.o
	g++ shortest_ams_lex.cpp ams_lex.o data_structures.o -O3 -fopenmp -std=c++11 -o shortest_ams_lex

ams_lex.o: ams_lex.cpp ams_lex.h
	g++ ams_lex.cpp -O3 -fopenmp -std=c++11 -c -o ams_lex.o

data_structures.o: data_structures.cpp data_structures.h
	g++ data_structures.cpp -O3 -fopenmp -std=c++11 -c -o data_structures.o

