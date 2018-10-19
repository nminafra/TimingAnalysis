INS=-Iinclude

CC=g++
CFLAGS=-g -std=gnu++0x -Wall -DNUM_DOUBLE
SOFLAGS=-fPIC

ROOT=`root-config --cflags --libs`

all:TimingAnalysis.a

src/TimingAnalysis.o: src/TimingAnalysis.cxx include/TimingAnalysis.h
	$(CC) $(SOFLAGS) $(CFLAGS) -c -o $@ $< $(ROOT) $(INS)
	
src/pulse.o: src/pulse.cxx include/pulse.h
	$(CC) $(SOFLAGS) $(CFLAGS) -c -o $@ $< $(ROOT) $(INS)

objects=TimingAnalysis.o pulse.o
f_0_objects=$(objects:%=src/%)

TimingAnalysis.a:$(f_0_objects) include/TimingAnalysis.h include/pulse.h
	ar -r $@ $(f_0_objects)

.cxx.o:
	$(CC) $(CFLAGS) $(SOFLAGS) $(INS) -c $< -o $@ $< $(ROOT) $(INS)

clean:
	rm -f *.a src/*.o
