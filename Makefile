ROOT=`root-config --cflags --libs`

LIB_DIRECTORY=TimingAnalysis
LIB_NAME=TimingAnalysis
CHECKFORUPDATE=$(LIB_DIRECTORY)/$(LIB_NAME).a $(LIB_DIRECTORY)/include/$(LIB_NAME).h $(LIB_DIRECTORY)/include/timingAlgorithm.h

progs=$(LIB_DIRECTORY)/$(LIB_NAME).a example_analyzeData analyzeDataVsMCP

all:$(progs)

$(LIB_DIRECTORY)/$(LIB_NAME).a: $(objects:%=$(LIB_DIRECTORY)/src/%) $(LIB_DIRECTORY)/include/timingAlgorithm.h
	$(MAKE) -C $(LIB_DIRECTORY)

example_analyzeData:%:%.cxx $(CHECKFORUPDATE)
	g++ -o $@ $< -I$(LIB_DIRECTORY)/include  $(LIB_DIRECTORY)/$(LIB_NAME).a $(ROOT) -z muldefs -O3 -std=c++11

analyzeDataVsMCP:%:%.cxx $(CHECKFORUPDATE)
	g++ -o $@ $< -I$(LIB_DIRECTORY)/include  $(LIB_DIRECTORY)/$(LIB_NAME).a $(ROOT) -z muldefs -O3 -std=c++11

clean:
	rm -r example_analyzeData $(LIB_DIRECTORY)/*.a $(LIB_DIRECTORY)/src/*.o
