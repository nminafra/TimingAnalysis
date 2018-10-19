#ifndef TimingAnalysis_cxx
#define TimingAnalysis_cxx

#include "TimingAnalysis.h"

//----------constructor-----------//

TimingAnalysis::TimingAnalysis(TChain * tree): pulse(tree) {}

//----------analysis methods------------//

void TimingAnalysis::initialize() {
  //Starting the watch
  std::cout << "*****************************************" << std::endl;
  m_timer.Start();
}

void TimingAnalysis::finalize() {
}

#endif
