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
  float threshold_MCP=0.01;
  float threshold=-0.1;
  float lowpass=0;

  // Additional parameters
  float hysteresis=1e-3;
  float min_amplitude_MCP=0;
  float max_amplitude_MCP=0.034;
  float min_amplitude_ch2=0;
  float max_amplitude_ch2=0.2;
  float baseline_p=0.1; //Use first 10% of the samples to compute baseline

  for (int i=1; i<argc; ++i)
  {
    if (argv[i][0] == '-') {
      // Option found
      std::string option(argv[i]);
      if ( option == "-h" || option == "--help" ) {
        std::cout << "List of options: " << std::endl;
        std::cout << "-h [ --help ]                         produce help message" << std::endl;
        std::cout << "-f [ --channel ] arg (=0)             channel to analyze" << std::endl;
        std::cout << "-c [ --cfd_threshold ] arg (=0.4)     CFD fraction" << std::endl;
        std::cout << "                                      a negative value will start a scan with" << std::endl;
        std::cout << "                                      a step equal to |cfd_threshold|" << std::endl;
        std::cout << "-t [ --threshold ] arg (=-0.1)        Threshold, negative for" << std::endl;
        std::cout << "                                      negative signals (V)" << std::endl;
        std::cout << "-p [ --lowpass ] arg (=0)             Lowpass filter frequency (Hz)" << std::endl;
        std::cout << "-o [ --outputdir ] arg (=./Results)   output directory" << std::endl;
        std::cout << "-i [ --filename ] arg                 input file" << std::endl;
        return 0;
      }
      std::string value(argv[++i]);

      if ( option == "-f" || option == "--channel" )
        secondchannel = std::stoi(value);
      if ( option == "-i" || option == "--filename" )
        filename = value;
      if ( option == "-c" || option == "--cfd_threshold" )
        cfd_threshold = std::stof(value);
      if ( option == "-t" || option == "--threshold" )
        threshold = std::stof(value);
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
  TString filenameTail("_result_MCP_");
  filenameTail+=secondchannel;
  filename.insert(filename.size()-5,filenameTail.Data());
  filename.erase(0,filename.find_last_of('/',filename.size()));
  filename.insert(0,outputdir);
  TFile * f_root = new TFile (filename.c_str(),"RECREATE");

  std::cout<<"Filling "<<filename<<" with "<<filename<<std::endl;

  if (cfd_threshold > 0) {
    AlgorithmParameters par(0.5, cfd_threshold, threshold_MCP,threshold,lowpass,hysteresis,min_amplitude_MCP, max_amplitude_MCP, min_amplitude_ch2, max_amplitude_ch2,baseline_p);
    double timeres_ps = example_analyzeData.executeTimeDifference<AlgorithmParameters>(f_root, ComputeExactTimeCFD, par, firstchannel, secondchannel)*1e12;
    std::cout << "\t\t\tTime difference: " << timeres_ps << " ps" << std::endl;
  }
  else {
    TDirectory* cfd_dir = f_root->mkdir("cfd_scan");
    TGraph cfd_graph;
    cfd_graph.SetName("cfd_scan");
    int cfd_counter=0;
    for (float cfd_th=std::abs(cfd_threshold); cfd_th<1; cfd_th+=std::abs(cfd_threshold)) {
      std::cout<< "############ CFD fraction " << cfd_th << " ############" <<std::endl;
      TString cfd_tmpdir_name("cfd_");
      cfd_tmpdir_name += (int) (100*cfd_th);
      cfd_tmpdir_name += "_percent";
      TDirectory* cfd_tmpdir = f_root->mkdir(cfd_tmpdir_name);
      cfd_tmpdir->cd();
      AlgorithmParameters par( 0.5, cfd_th, threshold_MCP,threshold,lowpass,hysteresis,min_amplitude_MCP, max_amplitude_MCP, min_amplitude_ch2, max_amplitude_ch2,baseline_p);  //Fixed cfd_threshold for ch0
      double timeres_ps = example_analyzeData.executeTimeDifference<AlgorithmParameters>(f_root, ComputeExactTimeCFD, par, firstchannel, secondchannel)*1e12;
      std::cout << "\t\t\tTime difference: " << timeres_ps << " ps" << std::endl;
      cfd_graph.SetPoint(cfd_counter++, cfd_th, timeres_ps);
    }
    cfd_dir->cd();
    cfd_graph.Write();
  }


  f_root->Close();
  return 0;
}
