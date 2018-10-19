//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 18 21:37:10 2018 by ROOT version 6.06/08
// from TTree pulse/Digitized waveforms
// found on file: DataNetScope_Run1276.root
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
   Float_t         channel[4][1000];
   Float_t         time[1][1000];
   Float_t         baseline[4];
   Float_t         baseline_RMS[4];
   Float_t         amp[4];
   Float_t         t_peak[4];
   Float_t         integral[4];
   Float_t         intfull[4];
   Float_t         risetime[4];
   Float_t         decaytime[4];
   Float_t         gaus_mean[4];
   Float_t         gaus_sigma[4];
   Float_t         gaus_chi2[4];
   Float_t         LP1_10[4];
   Float_t         LP1_20[4];
   Float_t         LP1_30[4];
   Float_t         LP2_10[4];
   Float_t         LP2_20[4];
   Float_t         LP2_30[4];

   // List of branches
   TBranch        *b_i_evt;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_time;   //!
   TBranch        *b_baseline;   //!
   TBranch        *b_baseline_RMS;   //!
   TBranch        *b_amp;   //!
   TBranch        *b_t_peak;   //!
   TBranch        *b_integral;   //!
   TBranch        *b_intfull;   //!
   TBranch        *b_risetime;   //!
   TBranch        *b_decaytime;   //!
   TBranch        *b_gaus_mean;   //!
   TBranch        *b_gaus_sigma;   //!
   TBranch        *b_gaus_chi2;   //!
   TBranch        *b_LP1_10;   //!
   TBranch        *b_LP1_20;   //!
   TBranch        *b_LP1_30;   //!
   TBranch        *b_LP2_10;   //!
   TBranch        *b_LP2_20;   //!
   TBranch        *b_LP2_30;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DataNetScope_Run1276.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("DataNetScope_Run1276.root");
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
   fChain->SetBranchAddress("amp", amp, &b_amp);
   fChain->SetBranchAddress("t_peak", t_peak, &b_t_peak);
   fChain->SetBranchAddress("integral", integral, &b_integral);
   fChain->SetBranchAddress("intfull", intfull, &b_intfull);
   fChain->SetBranchAddress("risetime", risetime, &b_risetime);
   fChain->SetBranchAddress("decaytime", decaytime, &b_decaytime);
   fChain->SetBranchAddress("gaus_mean", gaus_mean, &b_gaus_mean);
   fChain->SetBranchAddress("gaus_sigma", gaus_sigma, &b_gaus_sigma);
   fChain->SetBranchAddress("gaus_chi2", gaus_chi2, &b_gaus_chi2);
   fChain->SetBranchAddress("LP1_10", LP1_10, &b_LP1_10);
   fChain->SetBranchAddress("LP1_20", LP1_20, &b_LP1_20);
   fChain->SetBranchAddress("LP1_30", LP1_30, &b_LP1_30);
   fChain->SetBranchAddress("LP2_10", LP2_10, &b_LP2_10);
   fChain->SetBranchAddress("LP2_20", LP2_20, &b_LP2_20);
   fChain->SetBranchAddress("LP2_30", LP2_30, &b_LP2_30);
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
