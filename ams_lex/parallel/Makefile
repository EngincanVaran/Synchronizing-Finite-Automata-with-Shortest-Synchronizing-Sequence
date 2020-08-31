all: ams_lex.o
	g++ shortest_ams_lex.cpp ams_lex.o -O3 -fopenmp -std=c++11 -o shortest_ams_lex

test: ams_lex.o
	g++ shortest_ams_lex_second.cpp ams_lex.o -O3 -fopenmp -std=c++11 -o shortest_ams_lex_second

ams_lex.o: ams_lex.cpp ams_lex.h
	g++ ams_lex.cpp -O3 -fopenmp -std=c++11 -c -o ams_lex.o
