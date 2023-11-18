main: main.o CTRNN.o TSearch.o Sniffer.o random.o OdorPuff.o
	g++ -pthread -o main main.o CTRNN.o TSearch.o Sniffer.o random.o OdorPuff.o
OdorPuff.o: OdorPuff.cpp OdorPuff.h 
	g++ -pthread -c -O3 OdorPuff.cpp
random.o: random.cpp random.h VectorMatrix.h
	g++ -pthread -c -O3 random.cpp
CTRNN.o: CTRNN.cpp random.h CTRNN.h
	g++ -pthread -c -O3 CTRNN.cpp
TSearch.o: TSearch.cpp TSearch.h
	g++ -pthread -c -O3 TSearch.cpp
Sniffer.o: Sniffer.cpp Sniffer.h TSearch.h CTRNN.h random.h VectorMatrix.h
	g++ -pthread -c -O3 Sniffer.cpp
main.o: main.cpp CTRNN.h Sniffer.h TSearch.h OdorPuff.h
	g++ -pthread -c -O3 main.cpp
clean:
	rm *.o main
