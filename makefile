all: compile link run
compile:
  g++ -c main.cpp

link:
  g++ -o main.out main.o

run:
  ./main.out