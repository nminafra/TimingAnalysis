//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov  9 11:22:10 2018 by ROOT version 6.08/07
// from TTree pulse/Digitized waveforms
// found on file: DataVMETiming_Run249_3.root
//////////////////////////////////////////////////////////

#ifndef pulse_h
#define pulse_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class pulse {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          i_evt;
   Float_t         channel[36][1024];
   Float_t         time[4][1024];
   Float_t         baseline[36];
   Float_t         baseline_RMS[36];
   Float_t         noise[36];
   Float_t         amp[36];
   Float_t         t_peak[36];
   Float_t         integral[36];
   Float_t         intfull[36];
   Float_t         risetime[36];
   Float_t         decaytime[36];
   Float_t         LP2_20[36];
   Float_t         LP2_30mV[36];
   Float_t         InterpolatedAmp[36];
   Float_t         t0_30[36];
   Float_t         t1_30[36];
   Float_t         tot_30[36];
   Float_t         t0CFD_20[36];
   Float_t         t1CFD_20[36];
   Float_t         totCFD_20[36];
   Int_t           triggerNumber;
   UShort_t        corruption;
   Float_t         xIntercept;
   Float_t         yIntercept;
   Float_t         xSlope;
   Float_t         ySlope;
   Float_t         x_dut[4];
   Float_t         y_dut[4];
   Float_t         chi2;
   Int_t           ntracks;

   // List of branches
   TBranch        *b_i_evt;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_time;   //!
   TBranch        *b_baseline;   //!
   TBranch        *b_baseline_RMS;   //!
   TBranch        *b_noise;   //!
   TBranch        *b_amp;   //!
   TBranch        *b_t_peak;   //!
   TBranch        *b_integral;   //!
   TBranch        *b_intfull;   //!
   TBranch        *b_risetime;   //!
   TBranch        *b_decaytime;   //!
   TBranch        *b_LP2_20;   //!
   TBranch        *b_LP2_30mV;   //!
   TBranch        *b_InterpolatedAmp;   //!
   TBranch        *b_t0_30;   //!
   TBranch        *b_t1_30;   //!
   TBranch        *b_tot_30;   //!
   TBranch        *b_t0CFD_20;   //!
   TBranch        *b_t1CFD_20;   //!
   TBranch        *b_totCFD_20;   //!
   TBranch        *b_triggerNumber;   //!
   TBranch        *b_corruption;   //!
   TBranch        *b_xIntercept;   //!
   TBranch        *b_yIntercept;   //!
   TBranch        *b_xSlope;   //!
   TBranch        *b_ySlope;   //!
   TBranch        *b_x_dut;   //!
   TBranch        *b_y_dut;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_ntracks;   //!

   pulse(TTree *tree=0);
   virtual ~pulse();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef pulse_cxx
pulse::pulse(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DataVMETiming_Run249_3.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("DataVMETiming_Run249_3.root");
      }
      f->GetObject("pulse",tree);

   }
   Init(tree);
}

pulse::~pulse()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pulse::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pulse::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void pulse::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("i_evt", &i_evt, &b_i_evt);
   fChain->SetBranchAddress("channel", channel, &b_channel);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("baseline", baseline, &b_baseline);
   fChain->SetBranchAddress("baseline_RMS", baseline_RMS, &b_baseline_RMS);
   fChain->SetBranchAddress("noise", noise, &b_noise);
   fChain->SetBranchAddress("amp", amp, &b_amp);
   fChain->SetBranchAddress("t_peak", t_peak, &b_t_peak);
   fChain->SetBranchAddress("integral", integral, &b_integral);
   fChain->SetBranchAddress("intfull", intfull, &b_intfull);
   fChain->SetBranchAddress("risetime", risetime, &b_risetime);
   fChain->SetBranchAddress("decaytime", decaytime, &b_decaytime);
   fChain->SetBranchAddress("LP2_20", LP2_20, &b_LP2_20);
   fChain->SetBranchAddress("LP2_30mV", LP2_30mV, &b_LP2_30mV);
   fChain->SetBranchAddress("InterpolatedAmp", InterpolatedAmp, &b_InterpolatedAmp);
   fChain->SetBranchAddress("t0_30", t0_30, &b_t0_30);
   fChain->SetBranchAddress("t1_30", t1_30, &b_t1_30);
   fChain->SetBranchAddress("tot_30", tot_30, &b_tot_30);
   fChain->SetBranchAddress("t0CFD_20", t0CFD_20, &b_t0CFD_20);
   fChain->SetBranchAddress("t1CFD_20", t1CFD_20, &b_t1CFD_20);
   fChain->SetBranchAddress("totCFD_20", totCFD_20, &b_totCFD_20);
   fChain->SetBranchAddress("triggerNumber", &triggerNumber, &b_triggerNumber);
   fChain->SetBranchAddress("corruption", &corruption, &b_corruption);
   fChain->SetBranchAddress("xIntercept", &xIntercept, &b_xIntercept);
   fChain->SetBranchAddress("yIntercept", &yIntercept, &b_yIntercept);
   fChain->SetBranchAddress("xSlope", &xSlope, &b_xSlope);
   fChain->SetBranchAddress("ySlope", &ySlope, &b_ySlope);
   fChain->SetBranchAddress("x_dut", x_dut, &b_x_dut);
   fChain->SetBranchAddress("y_dut", y_dut, &b_y_dut);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   Notify();
}

Bool_t pulse::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pulse::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pulse::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef pulse_cxx
