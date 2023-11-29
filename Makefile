main: main.o CTRNN.o TSearch.o Sniffer.o random.o Fluid.o
	g++ -std=c++11 -pthread -o main main.o CTRNN.o TSearch.o Sniffer.o random.o Fluid.o
Fluid.o: Fluid.cpp Fluid.h 
	g++ -std=c++11 -pthread -c -O3 Fluid.cpp
random.o: random.cpp random.h VectorMatrix.h
	g++ -std=c++11 -pthread -c -O3 random.cpp
CTRNN.o: CTRNN.cpp random.h CTRNN.h
	g++ -std=c++11 -pthread -c -O3 CTRNN.cpp
TSearch.o: TSearch.cpp TSearch.h
	g++ -std=c++11 -pthread -c -O3 TSearch.cpp
Sniffer.o: Sniffer.cpp Sniffer.h TSearch.h CTRNN.h random.h VectorMatrix.h
	g++ -std=c++11 -pthread -c -O3 Sniffer.cpp
main.o: main.cpp CTRNN.h Sniffer.h TSearch.h Fluid.h
	g++ -std=c++11 -pthread -c -O3 main.cpp
clean:
	rm *.o main
