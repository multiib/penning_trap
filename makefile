all: compile link run

compile:
	g++ -c src/penningtrap.cpp -std=c++11 -I Include
	g++ -c src/particle.cpp -std=c++11 -I Include
	g++ -c src/main.cpp -std=c++11 -I Include

link:
	g++ -o main.o $(wildcard *.o) -larmadillo

run:
	./main

quick:
	g++ src/main.cpp src/penningtrap.cpp src/particle.cpp -std=c++11 -I Include/ -larmadillo -o main.exe
	./main.exe