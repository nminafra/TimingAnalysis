#ifndef Oscilloscope_tree_cxx
#define Oscilloscope_tree_cxx

#include "Oscilloscope_analyzeData.h"

//----------constructor-----------//

Oscilloscope_analyzeData::Oscilloscope_analyzeData(TChain * tree): pulse(tree) {}

//----------analysis methods------------//

void Oscilloscope_analyzeData::initialize() {
  //Starting the watch
  std::cout << "*****************************************" << std::endl;
  m_timer.Start();
}

void Oscilloscope_analyzeData::finalize() {
}

#endif
