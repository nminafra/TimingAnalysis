INS=-Iinclude

CC=g++
CFLAGS=-g -std=gnu++0x -Wall -DNUM_DOUBLE
SOFLAGS=-fPIC

ROOT=`root-config --cflags --libs`

all:Oscilloscope_analyzeData.a

src/Oscilloscope_analyzeData.o: src/Oscilloscope_analyzeData.cxx include/Oscilloscope_analyzeData.h
	$(CC) $(SOFLAGS) $(CFLAGS) -c -o $@ $< $(ROOT) $(INS)
	
src/pulse.o: src/pulse.cxx include/pulse.h
	$(CC) $(SOFLAGS) $(CFLAGS) -c -o $@ $< $(ROOT) $(INS)

objects=Oscilloscope_analyzeData.o pulse.o
f_0_objects=$(objects:%=src/%)

Oscilloscope_analyzeData.a:$(f_0_objects)
	ar -r $@ $(f_0_objects)

.cxx.o:
	$(CC) $(CFLAGS) $(SOFLAGS) $(INS) -c $< -o $@ $< $(ROOT) $(INS)

clean:
	rm -f *.a src/*.o





