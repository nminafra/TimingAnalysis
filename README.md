# TimingAnalysis

Code for timing analysis. Author: Nicola Minafra

Usage:
$> make
$> ./example_analyzeData --help
  Allowed options:
  -h [ --help ]                         produce help message
  -f [ --firstchannel ] arg (=0)        First channel to analyze
  -s [ --secondchannel ] arg (=1)       Second channel to analyze
  -c [ --cfd_threshold ] arg (=0.4)
                                        CFD fraction
  -t [ --threshold_ch1 ] arg (=-0.1)
                                        Threshold for ch 1, negative for
                                        negative signals (V)
  -w [ --threshold_ch2 ] arg (=-0.1)
                                        Threshold for ch 2, negative for
                                        negative signals (V)
  -p [ --lowpass ] arg (=0)             Lowpass filter frequency (Hz)
  -j [ --hysteresis ] arg (=0.001)
                                        Hysteresis for the discriminator (V)
  -t [ --treename ] arg (=pulse)        Name of the TTree
  -o [ --outputdir ] arg (=~/Work/public/Fermilab/Results)
                                        output directory
  -i [ --filename ] arg                 input file


For Fermilab analysis:
$> ./example_analyzeData -i /afs/cern.ch/work/n/nminafra/public/Fermilab/ -t -0.1 -w -0.01 -f 0 -s 1 -c 0.5 --lowpass 700e6

The output root file is saved in Results/ with the same name of the input file, plus _result_ch0_ch1.root
Important plots:
evN: graph of waveforms
h_deltat_Smart: time difference between the two channels with the method ComputeExactTimeCFD in /TimingAnalysis/include/timingAlgorithm.h
h_max_selected_DetN: amplitude of only events selected for timing analysis
h_SNR_DetN: SNR computed event by event


NOTE: Root required, on lxplus:
$> . /afs/cern.ch/sw/lcg/app/releases/ROOT/6.00.02/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh
$> . /afs/cern.ch/sw/lcg/external/gcc/4.9.3o3/x86_64-slc6-gcc49-opt/setup.sh
$> make
