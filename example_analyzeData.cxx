#include <TimingAnalysis.h>
#include <timingAlgorithm.h>

#include <TFileCollection.h>

int main (int argc, char** argv)
{
  std::string filename("");
  std::string outputdir("./Results");
  int firstchannel=0;
  int secondchannel=1;
  float cfd_threshold=0.4;
  float threshold_ch1=-0.1;
  float threshold_ch2=-0.1;
  float lowpass=0;
  float hysteresis=1e-3;

  for (int i=1; i<argc; ++i)
  {
    if (argv[i][0] == '-') {
      // Option found
      std::string option(argv[i]);
      if ( option == "-h" || option == "--help" ) {
        std::cout << "List of options: " << std::endl;
        std::cout << "-h [ --help ]                         produce help message" << std::endl;
        std::cout << "-f [ --firstchannel ] arg (=0)        First channel to analyze" << std::endl;
        std::cout << "-s [ --secondchannel ] arg (=1)       Second channel to analyze" << std::endl;
        std::cout << "-c [ --cfd_threshold ] arg (=0.4)     CFD fraction" << std::endl;
        std::cout << "-t [ --threshold_ch1 ] arg (=-0.1)    Threshold for ch 1, negative for" << std::endl;
        std::cout << "                                      negative signals (V)" << std::endl;
        std::cout << "-w [ --threshold_ch2 ] arg (=-0.1)    Threshold for ch 2, negative for" << std::endl;
        std::cout << "                                      negative signals (V)" << std::endl;
        std::cout << "-p [ --lowpass ] arg (=0)             Lowpass filter frequency (Hz)" << std::endl;
        std::cout << "-o [ --outputdir ] arg (=./Results)   output directory" << std::endl;
        std::cout << "-i [ --filename ] arg                 input file" << std::endl;
        return 0;
      }
      std::string value(argv[++i]);

      if ( option == "-f" || option == "--firstchannel" )
        firstchannel = std::stoi(value);
      if ( option == "-s" || option == "--secondchannel" )
        secondchannel = std::stoi(value);
      if ( option == "-i" || option == "--filename" )
        filename = value;
      if ( option == "-c" || option == "--cfd_threshold" )
        cfd_threshold = std::stof(value);
      if ( option == "-t" || option == "--threshold_ch1" )
        threshold_ch1 = std::stof(value);
      if ( option == "-w" || option == "--threshold_ch2" )
        threshold_ch2 = std::stof(value);
      if ( option == "-p" || option == "--lowpass" )
        lowpass = std::stof(value);
      if ( option == "-o" || option == "--outputdir" )
        outputdir = value;
      if ( option == "-i" || option == "--filename" )
        filename = value;
    }
  }

  if (filename == "") {
    std::cout << "Input file required! For help use:" << std::endl;
    std::cout << argv[0] << " --help" << std::endl;
    return 0;
  }

  // Creating the analysis object from data TTree
  TFileCollection * fc_example = new TFileCollection();
  TChain * c_example = new TChain("pulse");
  fc_example->Add(filename.c_str());
  c_example->AddFileInfoList((TCollection*)fc_example->GetList());

  TimingAnalysis example_analyzeData(c_example);

  // Output file
  TString filenameTail("_result_");
  filenameTail+=firstchannel;
  filenameTail+="_";
  filenameTail+=secondchannel;
  filename.insert(filename.size()-5,filenameTail.Data());
  filename.erase(0,filename.find_last_of('/',filename.size()));
  filename.insert(0,outputdir);
  std::cout<<filename<<std::endl;
  TFile * f_root = new TFile (filename.c_str(),"RECREATE");

  std::cout<<"Filling "<<filename<<" with "<<filename<<std::endl;

//AlgorithmParameters(const double cfdRatio, const double threshold_ch0=.0, const double threshold_ch1=.0, const double sigma=0, const double hysteresis=1e-3, const double minCh0=-10., const double maxCh0=10., const double minCh1=-10., const double maxCh1=10., const double baseline_n=0.15)
  AlgorithmParameters par(cfd_threshold,threshold_ch1,threshold_ch2,lowpass,hysteresis,-0.3,0.3,-0.4,0.4,0.1);

   double timeres_ps = example_analyzeData.executeTimeDifference<AlgorithmParameters>(f_root, ComputeExactTimeCFD, par, firstchannel, secondchannel)*1e12;
   std::cout << "\t\t\tTime difference: " << timeres_ps << " ps" << std::endl;

  f_root->Close();
  return 0;
}
