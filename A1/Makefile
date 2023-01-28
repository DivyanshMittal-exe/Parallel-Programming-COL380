#g++ -std=c++11 -O2 -mavx2 -c matrify.cpp
CC=g++
CFLAGS=-std=c++11 -O2 -g -fopenmp

sources=main.cpp matrify.cpp

objects=$(sources:.cpp=.o)

matrify:$(objects)
	$(CC) $(CFLAGS) $^ -o $@

run: matrify
	./matrify


#matrify.o: matrify.h
#	touch matrify.o
#	$(CC) $(CFLAGS) -c $<
#main.o: matrify.h
#	$(CC) $(CFLAGS) -c $<

%.o: %.cpp matrify.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o matrify
