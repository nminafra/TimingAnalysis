#include <Oscilloscope_analyzeData.h>
#include <timingAlgorithm.h>

#include <TFileCollection.h>
#include <boost/lexical_cast.hpp>
#include "boost/program_options.hpp"

namespace po = boost::program_options;

int main (int argc, char** argv)
{
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
  ("help,h", "produce help message")
  ("firstchannel,f", po::value<int>()->default_value(0), "First channel to analyze")
  ("secondchannel,s", po::value<int>()->default_value(1), "Second channel to analyze")
  ("cfd_threshold,c", po::value<double>()->default_value(0.4), "CFD fraction")
  ("threshold_ch1,t", po::value<double>()->default_value(-0.1), "Threshold for ch 1, negative for negative signals (V)")
  ("threshold_ch2,w", po::value<double>()->default_value(-0.1), "Threshold for ch 2, negative for negative signals (V)")
  ("lowpass,p", po::value<double>()->default_value(0), "Lowpass filter frequency (Hz)")
  ("hysteresis,h", po::value<double>()->default_value(1e-3), "Hysteresis for the discriminator (V)")
  ("treename,n", po::value<std::string>()->default_value("pulse"), "Name of the TTree")
  ("outputdir,o", po::value<std::string>()->default_value("./Results"), "output directory")
  ("filename,i", po::value<std::string>(), "input file");
  
  po::positional_options_description positionalOptions; 
  positionalOptions.add("filename", 1); 
  
  po::variables_map vm; 
  po::store(po::command_line_parser(argc, argv).options(desc).positional(positionalOptions).run(), vm);
  
  if ( vm.count("help") ) { 
    std::cout << desc << std::endl; 
    return 0; 
  } 
  
  // Creating the analysis object from data TTree
  std::string filename(vm["filename"].as<std::string>());
  TFileCollection * fc_example = new TFileCollection();
  TChain * c_example = new TChain(vm["treename"].as<std::string>().c_str());
  fc_example->Add(filename.c_str());
  c_example->AddFileInfoList((TCollection*)fc_example->GetList());

  Oscilloscope_analyzeData example_analyzeData(c_example);

  // Output file
  TString filenameTail("_result_");
  filenameTail+=vm["firstchannel"].as<int>();
  filenameTail+="_";
  filenameTail+=vm["secondchannel"].as<int>();
  filename.insert(filename.size()-5,filenameTail.Data());
  filename.erase(0,filename.find_last_of('/',filename.size()));
  filename.insert(0,vm["outputdir"].as<std::string>());
  std::cout<<filename<<std::endl;
  TFile * f_root = new TFile (filename.c_str(),"RECREATE");  
    
  std::cout<<"Filling "<<filename<<" with "<<vm["filename"].as<std::string>()<<std::endl;    
   

  AlgorithmParameters par(vm["cfd_threshold"].as<double>(),vm["threshold_ch1"].as<double>(),vm["threshold_ch2"].as<double>(),vm["lowpass"].as<double>(),vm["hysteresis"].as<double>(),-0.3,0.3,-0.4,0.4);

   double timeres_ps = example_analyzeData.executeTimeDifference<AlgorithmParameters>(f_root, ComputeExactTimeCFD, par, vm["firstchannel"].as<int>(), vm["secondchannel"].as<int>())*1e12;
   std::cout << "\t\t\tTime difference: " << timeres_ps << " ps" << std::endl;
  
  f_root->Close();
  return 0;
}
