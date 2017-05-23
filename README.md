# TimingAnalysis

Code for timing analysis. Author: Nicola Minafra

Usage:
&> make
&> ./example_analyzeData --help
  Allowed options:
  -h [ --help ]                         produce help message
  -f [ --firstchannel ] arg (=0)        First channel to analyze
  -s [ --secondchannel ] arg (=1)       Second channel to analyze
  -c [ --cfd_threshold ] arg (=0.40000000000000002)
                                        CFD fraction
  -t [ --threshold_ch1 ] arg (=-0.10000000000000001)
                                        Threshold for ch 1, negative for 
                                        negative signals (V)
  -w [ --threshold_ch2 ] arg (=-0.10000000000000001)
                                        Threshold for ch 2, negative for 
                                        negative signals (V)
  -f [ --lowpass ] arg (=0)             Lowpass filter frequency (Hz)
  -h [ --hysteresis ] arg (=0.0030000000000000001)
                                        Hysteresis for the discriminator (V)
  -t [ --treename ] arg (=pulse)        Name of the TTree
  -o [ --outputdir ] arg (=~/Work/public/Fermilab/Results)
                                        output directory
  -i [ --filename ] arg                 input file
