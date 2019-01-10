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
  auto evtwt = weightCalc(HistoTot, input);
  cout << "LumiWeight is " << evtwt << "\n";

  auto nentries_wtn = Run_Tree->GetEntries();
  cout << "nentries_wtn====" << nentries_wtn << "\n";

  // Declare Vectors for muon and tau from Z
  TLorentzVector Mu4Momentum, AntiIsoMu4Momentum, Tau4Momentum, Jet4Momentum;

  // Begin the event loop.
  for (Int_t i = 0; i < nentries_wtn; i++) {
    Run_Tree->GetEntry(i);

    // Variables needed later.
    bool OS(false), SS(false);
    vector<float> btag_pt;
    vector<int> good_muon_charge, good_tau_charge;

    // Give a progress report while running.
    if (i % 1000 == 0) {
      fprintf(stdout, "\r  Processed events: %8d of %8lld ", i, nentries_wtn);
      fflush(stdout);
    }

    // Apply muon trigger: HLT_IsoMu24_v
    // do this now so we don't waste time processing an event that will fail
    // definition: https://github.com/cmkuo/ggAnalysis/blob/master/ggNtuplizer/plugins/ggNtuplizer_globalEvent.cc#L151
    bool PassTrigger = (HLTEleMuX >> 19 & 1) == 1;
    if (!PassTrigger) {
      continue;
    }

    int numTau(0);
    for (int igen = 0; igen < nMC; igen++) {
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

    // Loop over muons in the event.
    for (auto imu = 0; imu < nMu; imu++) {
      // Muon kinematic selection.
      if (muPt->at(imu) < 30 || fabs(muEta->at(imu)) > 2.1) {
        continue;
      }

      if (fabs(muD0->at(imu)) > 0.045 || fabs(muDz->at(imu)) > 0.2) {
          continue;
        }

        // Apply medium muon ID.
        // definition: https://github.com/cmkuo/ggAnalysis/blob/master/ggNtuplizer/plugins/ggNtuplizer_muons.cc#L166-L171
        bool PassID = (muIDbit->at(imu) >> 1 & 1) == 1;
        if (!PassID) {
          continue;
      }

      // Transverse mass cut to remove W events.
      float MuMetTranverseMass = TMass_F(muPt->at(imu), muPt->at(imu) * cos(muPhi->at(imu)), muPt->at(imu) * sin(muPhi->at(imu)), pfMET, pfMETPhi);
      if (MuMetTranverseMass > 40) {
        continue;
      }

      // Calculate muon isolation.
      float IsoMu = muPFChIso->at(imu) / muPt->at(imu);
      if ((muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5 * muPFPUIso->at(imu)) > 0.0) {
        IsoMu += (muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5 * muPFPUIso->at(imu)) / muPt->at(imu);
      }

      // Apply muon isolation.
      if (IsoMu > 0.3) {
        AntiIsoMu4Momentum.SetPtEtaPhiE(muPt->at(imu), muEta->at(imu), muPhi->at(imu), MuMass);
      } else {
        // Good muon, now add charge to vector and fill muon P4.
        good_muon_charge.push_back(muCharge->at(imu));
        Mu4Momentum.SetPtEtaPhiE(muPt->at(imu), muEta->at(imu), muPhi->at(imu), MuMass);
      }
    }  // End of muon loop

    // Apply dimuon veto and also make sure we found at least 1 good muon.
    if (good_muon_charge.size() > 1 && good_muon_charge.at(0) * good_muon_charge.at(1) < 0) {
      continue;
    } else if (good_muon_charge.size() == 0) {
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
      good_tau_charge.push_back(tauCharge->at(itau));
      Tau4Momentum.SetPtEtaPhiM(tauPt->at(itau), tauEta->at(itau), tauPhi->at(itau), tauMass->at(itau));
    }  // End of tau loop

    // Apply ditau veto and also make sure we found at least 1 good tau.
    if (good_tau_charge.size() > 1 && good_tau_charge.at(0) * good_tau_charge.at(1) < 0) {
      continue;
    } else if (good_tau_charge.size() == 0) {
      continue;
    }

    // Set OS and SS.
    if (good_muon_charge.at(0) * good_tau_charge.at(0) < 0) {
      OS = true;
    } else {
      SS = true;
    }

    // Count the number of btagged jets.
    for (int ijet = 0; ijet < nJet; ijet++) {
      Jet4Momentum.SetPtEtaPhiE(jetPt->at(ijet), jetEta->at(ijet), jetPhi->at(ijet), jetEn->at(ijet));
      if (jetPt->at(ijet) > 20 && fabs(jetEta->at(ijet)) < 2.5 && jetCSV2BJetTags->at(ijet) > 0.8484 &&
          Jet4Momentum.DeltaR(Tau4Momentum) > 0.5 && Jet4Momentum.DeltaR(Mu4Momentum) > 0.5) {
        btag_pt.push_back(jetPt->at(ijet));
      }
    }

    // Apply b-tag veto to supress ttbar.
    if (btag_pt.size() > 0) {
      continue;
    }

    // Last, apply dR(mu, tau) > 0.5 selection.
    if (Mu4Momentum.DeltaR(Tau4Momentum) < 0.5) {
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
      visibleMassOS->Fill((Mu4Momentum + Tau4Momentum).M(), evtwt);
      visibleMassOSAntiIso->SetDefaultSumw2();
      visibleMassOSAntiIso->Fill((AntiIsoMu4Momentum + Tau4Momentum).M(), evtwt);
    } else if (SS) {
      visibleMassSS->SetDefaultSumw2();
      visibleMassSS->Fill((Mu4Momentum + Tau4Momentum).M(), evtwt);
      visibleMassSSAntiIso->SetDefaultSumw2();
      visibleMassSSAntiIso->Fill((AntiIsoMu4Momentum + Tau4Momentum).M(), evtwt);
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
