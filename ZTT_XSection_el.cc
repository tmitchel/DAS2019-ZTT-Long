#include <ostream>
#include <string>
#include "TreeReader.h"
#include "WeightCalculator.h"

//////////////////////////////////////////////////////////////////////////
//   Compiling the code:   ./Make.sh ZTT_XSection.cc                    //
//   Running the code:     ./ZTT_XSection.exe OutPut.root   Input.root  //
//////////////////////////////////////////////////////////////////////////

// Do some definitions up front
void setBranches(TTree *);
static float MuMass(0.10565837);
static float eleMass(0.000511);

int main(int argc, char **argv) {
  using std::cout;
  using std::endl;

  std::string out = *(argv + 1);
  cout << "OUTPUT NAME IS:    " << out << endl;  // PRINTING THE OUTPUT FILE NAME
  TFile *fout = TFile::Open((out+".root").c_str(), "RECREATE");

  std::string input = *(argv + 2);
  cout << "INPUT NAME IS:    " << input << endl;  // PRINTING THE INPUT FILE NAME
  TFile *myFile = TFile::Open(input.c_str());

  // Number of events processed.
  TH1F *HistoTot = reinterpret_cast<TH1F *>(myFile->Get("hEvents"));

  // Open pileup input files and read reweighting histograms.
  TFile *PUData = new TFile("data/MyDataPileupHistogram2016.root");
  TFile *PUMC = new TFile("data/mcMoriondPU.root");
  TH1F *HistoPUData = reinterpret_cast<TH1F *>(PUData->Get("pileup"));
  TH1F *HistoPUMC = reinterpret_cast<TH1F *>(PUMC->Get("pileup"));

  // Normalize histograms.
  HistoPUData->Scale(1.0 / HistoPUData->Integral());
  HistoPUMC->Scale(1.0 / HistoPUMC->Integral());

  // add the histrograms of muon and tau visible mass (both for opposite sign and same sign pair )
  fout->cd();
  TH1F *visibleMassOS = new TH1F("visibleMassOS", "visibleMassOS", 30, 0, 300);
  TH1F *visibleMassSS = new TH1F("visibleMassSS", "visibleMassSS", 30, 0, 300);
  TH1F *visibleMassOSAntiIso = new TH1F("visibleMassOSAntiIso", "visibleMassOSAntiIso", 30, 0, 300);
  TH1F *visibleMassSSAntiIso = new TH1F("visibleMassSSAntiIso", "visibleMassSSAntiIso", 30, 0, 300);
  TH1F *visibleMassOSRelaxedTauIso = new TH1F("visibleMassOSRelaxedTauIso", "visibleMassOSRelaxedTauIso", 30, 0, 300);
  TH1F *visibleMassSSRelaxedTauIso = new TH1F("visibleMassSSRelaxedTauIso", "visibleMassSSRelaxedTauIso", 30, 0, 300);

  TTree *Run_Tree = reinterpret_cast<TTree *>(myFile->Get("EventTree"));
  cout.setf(ios::fixed, ios::floatfield);

  setBranches(Run_Tree);
  double evtwt;
  auto lumi_weight = weightCalc(HistoTot, input);
  cout << "LumiWeight is " << lumi_weight << "\n";

  auto nentries_wtn = Run_Tree->GetEntries();
  cout << "nentries_wtn====" << nentries_wtn << "\n";

  // Declare Vectors for electron and tau from Z
  TLorentzVector El4Momentum, AntiIsoEle4Momentum, Tau4Momentum, Jet4Momentum;

  // Begin the event loop.
  for (Int_t i = 0; i < nentries_wtn; i++) {
    Run_Tree->GetEntry(i);
    evtwt = lumi_weight;

    // Variables needed later.
    bool OS(false), SS(false);
    vector<float> btag_pt;
    vector<int> good_ele_idx, good_tau_idx;

    // Give a progress report while running.
    if (i % 1000 == 0) {
      fprintf(stdout, "\r  Processed events: %8d of %8lld ", i, nentries_wtn);
      fflush(stdout);
    }

    // Apply muon trigger: HLT_Ele25_eta2p1_WPTight_Gsf_v
    // do this now so we don't waste time processing an event that will fail
    // definition: https://github.com/cmkuo/ggAnalysis/blob/master/ggNtuplizer/plugins/ggNtuplizer_globalEvent.cc#L132
    bool PassTrigger = (HLTEleMuX >> 0 & 1) == 1;
    if (!PassTrigger) {
      continue;
    }

    int numTau(0);
    for (int igen = 0; igen < nMC; igen++) {
      // tau: ID == 15 && Z: ID == 23
      if (fabs(mcPID->at(igen)) == 15 && mcMomPID->at(igen) == 23) {
        numTau++;
      }
    }

    // Z->ll must have 0 gen taus
    // Z->tautau must have 2 gen taus
    if (out.find("ToTauTau") != string::npos && numTau < 1) {
      continue;
    } else if (out.find("ToLL") != string::npos && numTau > 1) {
      continue;
    }

    // Loop over electrons in the event.
    for (auto iele = 0; iele < nEle; iele++) {
      // Electron kinematic selection.
      if (elePt->at(iele) < 15 || fabs(eleEta->at(iele)) > 2.4) {
        continue;
      }

      if (fabs(eleD0->at(iele)) > 0.045 || fabs(eleDz->at(iele)) > 0.2) {
        continue;
      }

      bool eleMVAId(false);
      if (fabs(eleSCEta->at(iele)) <= 0.8 && eleIDMVA->at(iele) > 0.941) {
        eleMVAId = true;
      } else if (fabs(eleSCEta->at(iele)) > 0.8 && fabs(eleSCEta->at(iele)) <= 1.5 && eleIDMVA->at(iele) > 0.899) {
        eleMVAId = true;
      } else if (fabs(eleSCEta->at(iele)) >= 1.5 && eleIDMVA->at(iele) > 0.758) {
        eleMVAId = true;
      }

      // Calculate electron isolation.
      float IsoEle = elePFChIso->at(iele) / elePt->at(iele);
      IsoEle += std::max(0., (elePFNeuIso->at(iele) + elePFPhoIso->at(iele) - 0.5 * elePFPUIso->at(iele)) / elePt->at(iele));

      // Apply electron isolation.
      if (IsoEle > 0.3) {
        AntiIsoEle4Momentum.SetPtEtaPhiE(elePt->at(iele), eleEta->at(iele), elePhi->at(iele), eleMass);
      } else {
        // Good electron, now add charge to vector and fill electron P4.
        good_ele_idx.push_back(iele);
        El4Momentum.SetPtEtaPhiM(elePt->at(iele), eleEta->at(iele), elePhi->at(iele), eleMass);
      }
    }  // End of electron loop

    // Apply dielectron veto and also make sure we found at least 1 good electron.
    if (good_ele_idx.size() > 1 && eleCharge->at(good_ele_idx.at(0)) * eleCharge->at(good_ele_idx.at(1)) < 0) {
      continue;
    } else if (good_ele_idx.size() == 0) {
      continue;
    }

    // Now apply tighter selection after we know there is only 1 good electron.
    auto idx = good_ele_idx.at(0);

    // Kinematic selection.
    if (elePt->at(idx) < 27 || fabs(eleEta->at(idx)) > 2.1) {
      continue;
    }

    // Transverse mass cut to remove W events.
    float ElMetTranverseMass = TMass_F(elePt->at(idx), elePt->at(idx) * cos(elePhi->at(idx)), elePt->at(idx) * sin(elePhi->at(idx)), pfMET, pfMETPhi);
    if (ElMetTranverseMass > 40) {
      continue;
    }

    // Loop over taus in the event.
    for (auto itau = 0; itau < nTau; itau++) {
      // Tau kinematic selection.
      if (tauPt->at(itau) < 30 || fabs(tauEta->at(itau)) > 2.3) {
        continue;
      }

      // Tau quality selection.
      if (taupfTausDiscriminationByDecayModeFinding->at(itau) < 0.5 ||
          !tauByTightMuonRejection3->at(itau) || !tauByMVA6LooseElectronRejection->at(itau) ||
          !tauByTightIsolationMVArun2v1DBoldDMwLT->at(itau)) {
        continue;
      }

      // Good tau, now add charge to vector and fill tau P4.
      good_tau_idx.push_back(itau);
      Tau4Momentum.SetPtEtaPhiM(tauPt->at(itau), tauEta->at(itau), tauPhi->at(itau), tauMass->at(itau));
    }  // End of tau loop

    // Apply ditau veto and also make sure we found at least 1 good tau.
    if (good_tau_idx.size() > 1 && tauCharge->at(good_tau_idx.at(0)) * tauCharge->at(good_tau_idx.at(1)) < 0) {
      continue;
    } else if (good_tau_idx.size() == 0) {
      continue;
    }

    // Set OS and SS.
    if (good_ele_idx.at(0) * good_tau_idx.at(0) < 0) {
      OS = true;
    } else {
      SS = true;
    }

    // Count the number of btagged jets.
    for (int ijet = 0; ijet < nJet; ijet++) {
      Jet4Momentum.SetPtEtaPhiE(jetPt->at(ijet), jetEta->at(ijet), jetPhi->at(ijet), jetEn->at(ijet));
      if (jetPt->at(ijet) > 20 && fabs(jetEta->at(ijet)) < 2.5 && jetCSV2BJetTags->at(ijet) > 0.8484 &&
          Jet4Momentum.DeltaR(Tau4Momentum) > 0.5 && Jet4Momentum.DeltaR(El4Momentum) > 0.5) {
        btag_pt.push_back(jetPt->at(ijet));
      }
    }

    // Apply b-tag veto to supress ttbar.
    if (btag_pt.size() > 0) {
      continue;
    }

    // Last, apply dR(mu, tau) > 0.5 selection.
    if (El4Momentum.DeltaR(Tau4Momentum) < 0.5) {
      continue;
    }

    // Do Pileup reweighting.
    if (!isData) {
      int puNUmmc = static_cast<int>(puTrue->at(0) * 10);
      int puNUmdata = static_cast<int>(puTrue->at(0) * 10);
      float PUMC_ = HistoPUMC->GetBinContent(puNUmmc + 1);
      float PUData_ = HistoPUData->GetBinContent(puNUmdata + 1);
      if (PUMC_ > 0) {
        evtwt *= (PUData_ / PUMC_);
      } else {
        cerr << "Divide by 0 found while PU reweighting" << endl;
      }
    }
    
    // Fill histograms.
    if (OS) {
      visibleMassOS->SetDefaultSumw2();
      visibleMassOS->Fill((El4Momentum + Tau4Momentum).M(), evtwt);
      visibleMassOSAntiIso->SetDefaultSumw2();
      visibleMassOSAntiIso->Fill((AntiIsoEle4Momentum + Tau4Momentum).M(), evtwt);
    } else if (SS) {
      visibleMassSS->SetDefaultSumw2();
      visibleMassSS->Fill((El4Momentum + Tau4Momentum).M(), evtwt);
      visibleMassSSAntiIso->SetDefaultSumw2();
      visibleMassSSAntiIso->Fill((AntiIsoEle4Momentum + Tau4Momentum).M(), evtwt);
    }
  }  // End Processing all entries

  // end of analysis code, close and write histograms/file
  fout->cd();
  fout->Write();
  fout->Close();
}

// Function to set all branch addresses.
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
