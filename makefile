compile:
	g++ src/main.cpp src/penningtrap.cpp src/particle.cpp src/analytical.cpp -std=c++11 -I Include/ -larmadillo -O3 -o main.exe
	./main.exe

plot:
	python res/*.py

compile plot:
	g++ src/main.cpp src/penningtrap.cpp src/particle.cpp src/analytical.cpp -std=c++11 -I Include/ -larmadillo -O3 -o main.exe
	./main.exe
	python res/*.py

