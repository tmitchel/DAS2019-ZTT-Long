////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Compiling the code:   ./Make.sh ZTT_XSection.cc
//   Running the code:     ./ZTT_XSection.exe OutPut.root   Input.root
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <ostream>
#include <string>
#include "TreeReader.h"
#include "WeightCalculator.h"

// Do some definitions up front
void setBranches(TTree *);
static float MuMass(0.10565837);
static float eleMass(0.000511);

int main(int argc, char **argv) {
  using std::cout;
  using std::endl;

  std::string out = *(argv + 1);
  cout << "OUTPUT NAME IS:    " << out << endl;  // PRINTING THE OUTPUT FILE NAME
  TFile *fout = TFile::Open(out.c_str(), "RECREATE");

  std::string input = *(argv + 2);
  cout << "INPUT NAME IS:    " << input << endl;  // PRINTING THE INPUT FILE NAME
  TFile *myFile = TFile::Open(input.c_str());

  TH1F *HistoTot = reinterpret_cast<TH1F *>(myFile->Get("hcount"));

  // add the histrograms of muon and tau visible mass (both for opposite sign and same sign pair )
  TH1F *visibleMassOS = new TH1F("visibleMassOS", "visibleMassOS", 30, 0, 300);
  TH1F *visibleMassSS = new TH1F("visibleMassSS", "visibleMassSS", 30, 0, 300);
  TH1F *visibleMassOSRelaxedTauIso = new TH1F("visibleMassOSRelaxedTauIso", "visibleMassOSRelaxedTauIso", 30, 0, 300);
  TH1F *visibleMassSSRelaxedTauIso = new TH1F("visibleMassSSRelaxedTauIso", "visibleMassSSRelaxedTauIso", 30, 0, 300);

  TTree *Run_Tree = reinterpret_cast<TTree *>(myFile->Get("EventTree"));
  cout.setf(ios::fixed, ios::floatfield);

  setBranches(Run_Tree);

  auto nentries_wtn = Run_Tree->GetEntries();
  cout << "nentries_wtn====" << nentries_wtn << "\n";
  for (Int_t i = 0; i < nentries_wtn; i++) {
    Run_Tree->GetEntry(i);

    if (i % 1000 == 0) {
      fprintf(stdout, "\r  Processed events: %8d of %8lld ", i, nentries_wtn);
      fflush(stdout);
    }

    TLorentzVector Mu4Momentum, Tau4Momentum;

    /////////////////////////////////////////////////
    // Important Analysis Loop Will Happen Here!!! //
    /////////////////////////////////////////////////
  }  // End Processing all entries

  // end of analysis code, close and write histograms/file
  fout->cd();
  visibleMassOS->Write();
  visibleMassSS->Write();
  visibleMassOSRelaxedTauIso->Write();
  visibleMassSSRelaxedTauIso->Write();
  fout->Close();
}

void setBranches(TTree *Run_Tree) {
  //########################################   General Info
  Run_Tree->SetBranchAddress("isData", &isData);
  Run_Tree->SetBranchAddress("run", &run);
  Run_Tree->SetBranchAddress("lumis", &lumis);
  Run_Tree->SetBranchAddress("event", &event);
  Run_Tree->SetBranchAddress("genWeight", &genWeight);
  Run_Tree->SetBranchAddress("HLTEleMuX", &HLTEleMuX);
  Run_Tree->SetBranchAddress("puTrue", &puTrue);
  Run_Tree->SetBranchAddress("nVtx", &nVtx);

  //########################################   MC Info
  Run_Tree->SetBranchAddress("nMC", &nMC);
  Run_Tree->SetBranchAddress("mcPID", &mcPID);
  Run_Tree->SetBranchAddress("mcStatus", &mcStatus);
  Run_Tree->SetBranchAddress("mcPt", &mcPt);
  Run_Tree->SetBranchAddress("mcEta", &mcEta);
  Run_Tree->SetBranchAddress("mcPhi", &mcPhi);
  Run_Tree->SetBranchAddress("mcE", &mcE);
  Run_Tree->SetBranchAddress("mcMass", &mcMass);
  Run_Tree->SetBranchAddress("mcMomPID", &mcMomPID);
  Run_Tree->SetBranchAddress("mcGMomPID", &mcGMomPID);

  //########################################   Tau Info
  Run_Tree->SetBranchAddress("nTau", &nTau);
  Run_Tree->SetBranchAddress("tauPt", &tauPt);
  Run_Tree->SetBranchAddress("tauEta", &tauEta);
  Run_Tree->SetBranchAddress("tauPhi", &tauPhi);
  Run_Tree->SetBranchAddress("tauMass", &tauMass);
  Run_Tree->SetBranchAddress("tauCharge", &tauCharge);

  Run_Tree->SetBranchAddress("taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding);

  Run_Tree->SetBranchAddress("tauByTightMuonRejection3", &tauByTightMuonRejection3);
  Run_Tree->SetBranchAddress("tauByLooseMuonRejection3", &tauByLooseMuonRejection3);

  Run_Tree->SetBranchAddress("tauByMVA6TightElectronRejection", &tauByMVA6TightElectronRejection);
  Run_Tree->SetBranchAddress("tauByMVA6MediumElectronRejection", &tauByMVA6MediumElectronRejection);
  Run_Tree->SetBranchAddress("tauByMVA6LooseElectronRejection", &tauByMVA6LooseElectronRejection);

  Run_Tree->SetBranchAddress("tauDxy", &tauDxy);
  Run_Tree->SetBranchAddress("tauDecayMode", &tauDecayMode);

  Run_Tree->SetBranchAddress("tauByLooseIsolationMVArun2v1DBoldDMwLT", &tauByLooseIsolationMVArun2v1DBoldDMwLT);
  Run_Tree->SetBranchAddress("tauByVLooseIsolationMVArun2v1DBoldDMwLT", &tauByVLooseIsolationMVArun2v1DBoldDMwLT);
  Run_Tree->SetBranchAddress("tauByTightIsolationMVArun2v1DBoldDMwLT", &tauByTightIsolationMVArun2v1DBoldDMwLT);

  //########################################   Mu Info
  Run_Tree->SetBranchAddress("nMu", &nMu);
  Run_Tree->SetBranchAddress("muPt", &muPt);
  Run_Tree->SetBranchAddress("muEta", &muEta);
  Run_Tree->SetBranchAddress("muPhi", &muPhi);
  Run_Tree->SetBranchAddress("muIsoTrk", &muIsoTrk);
  Run_Tree->SetBranchAddress("muCharge", &muCharge);
  Run_Tree->SetBranchAddress("muIDbit", &muIDbit);  // NEW
  Run_Tree->SetBranchAddress("muPFChIso", &muPFChIso);
  Run_Tree->SetBranchAddress("muPFPhoIso", &muPFPhoIso);
  Run_Tree->SetBranchAddress("muPFNeuIso", &muPFNeuIso);
  Run_Tree->SetBranchAddress("muPFPUIso", &muPFPUIso);
  Run_Tree->SetBranchAddress("muD0", &muD0);
  Run_Tree->SetBranchAddress("muDz", &muDz);

  //########################################   Ele Info
  Run_Tree->SetBranchAddress("nEle", &nEle);
  Run_Tree->SetBranchAddress("elePt", &elePt);
  Run_Tree->SetBranchAddress("eleEta", &eleEta);
  Run_Tree->SetBranchAddress("elePhi", &elePhi);
  Run_Tree->SetBranchAddress("elePFChIso", &elePFChIso);
  Run_Tree->SetBranchAddress("eleIDMVA", &eleIDMVA);  // NEW
  Run_Tree->SetBranchAddress("eleCharge", &eleCharge);
  Run_Tree->SetBranchAddress("eleSCEta", &eleSCEta);
  Run_Tree->SetBranchAddress("elePFChIso", &elePFChIso);
  Run_Tree->SetBranchAddress("elePFPhoIso", &elePFPhoIso);
  Run_Tree->SetBranchAddress("elePFNeuIso", &elePFNeuIso);
  Run_Tree->SetBranchAddress("elePFPUIso", &elePFPUIso);
  Run_Tree->SetBranchAddress("eleD0", &eleD0);
  Run_Tree->SetBranchAddress("eleDz", &eleDz);
  Run_Tree->SetBranchAddress("eleMissHits", &eleMissHits);
  Run_Tree->SetBranchAddress("eleConvVeto", &eleConvVeto);
  Run_Tree->SetBranchAddress("eleSCEta", &eleSCEta);

  //########################################   Jet Info
  Run_Tree->SetBranchAddress("nJet", &nJet);
  Run_Tree->SetBranchAddress("jetPt", &jetPt);
  Run_Tree->SetBranchAddress("jetEta", &jetEta);
  Run_Tree->SetBranchAddress("jetPhi", &jetPhi);
  Run_Tree->SetBranchAddress("jetEn", &jetEn);
  Run_Tree->SetBranchAddress("jetCSV2BJetTags", &jetCSV2BJetTags);
  Run_Tree->SetBranchAddress("jetPFLooseId", &jetPFLooseId);
  Run_Tree->SetBranchAddress("jetPUID", &jetPUID);
  Run_Tree->SetBranchAddress("jetRawPt", &jetRawPt);
  Run_Tree->SetBranchAddress("jetJECUnc", &jetJECUnc);
  Run_Tree->SetBranchAddress("jetRawEn", &jetRawEn);
  Run_Tree->SetBranchAddress("jetHadFlvr", &jetHadFlvr);

  //########################################   MET Info
  Run_Tree->SetBranchAddress("pfMET", &pfMET);
  Run_Tree->SetBranchAddress("pfMETPhi", &pfMETPhi);
  Run_Tree->SetBranchAddress("metFilters", &metFilters);
  Run_Tree->SetBranchAddress("genHT", &genHT);
}
