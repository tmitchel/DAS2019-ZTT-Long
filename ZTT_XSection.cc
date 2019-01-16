//////////////////////////////////////////////////////////////////////////
//   Compiling the code:   ./Make.sh ZTT_XSection.cc                    //
//   Running the code:     ./ZTT_XSection.exe OutPut.root   Input.root  //
//////////////////////////////////////////////////////////////////////////
#include <ostream>
#include <string>
#include "TreeReader.h"
#include "WeightCalculator.h"

int main(int argc, char **argv) {
  using std::cout;
  using std::endl;

  std::string out = *(argv + 1);

  cout << "\n\n\n OUTPUT NAME IS:    " << out << endl;  // PRINTING THE OUTPUT FILE NAME
  TFile *fout = TFile::Open(out.c_str(), "RECREATE");

  std::string input = *(argv + 2);
  cout << "\n\n\n INPUT NAME IS:    " << input << endl;  // PRINTING THE INPUT FILE NAME
  TFile *myFile = TFile::Open(input.c_str());
  TH1F *HistoTot = reinterpret_cast<TFile*>(myFile->Get("hcount"));

  // add the histrograms of muon and tau visible mass (both for opposite sign and same sign pair )
  TH1F *visibleMassOS = new TH1F("visibleMassOS", "visibleMassOS", 30, 0, 300);
  TH1F *visibleMassSS = new TH1F("visibleMassSS", "visibleMassSS", 30, 0, 300);

  TTree *Run_Tree = reinterpret_cast<TTree*>(myFile->Get("EventTree"));
  cout.setf(ios::fixed, ios::floatfield);

  // end of analysis code, close and write histograms/file
  fout->cd();
  visibleMassOS->Write();
  visibleMassSS->Write();
  visibleMassOSRelaxedTauIso->Write();
  visibleMassSSRelaxedTauIso->Write();
  fout->Close();
}
