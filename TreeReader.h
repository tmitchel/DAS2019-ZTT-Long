#ifndef CODEX_ANALYZER_H
#define	CODEX_ANALYZER_H

#include "TROOT.h"
#include <sstream>
#include <stdio.h>
#include <vector>
#include <utility>
#include <iostream>
#include <map>
using namespace std;


// Declaration of leaf types
Int_t           run;
Long64_t        event;
Int_t           lumis;
Bool_t          isData;
Int_t           nVtx;
Int_t           nGoodVtx;
Int_t           nTrksPV;
Bool_t          isPVGood;
Float_t         vtx;
Float_t         vty;
Float_t         vtz;
Float_t         rho;
Float_t         rhoCentral;
ULong64_t       HLTEleMuX;
ULong64_t       HLTPho;
ULong64_t       HLTJet;
ULong64_t       HLTEleMuXIsPrescaled;
ULong64_t       HLTPhoIsPrescaled;
ULong64_t       HLTJetIsPrescaled;
vector<int>     *phoPrescale;
vector<float>   *pdf;
vector<float>   *pdfSystWeight;
vector<string>   *pdfSystWeightId;
Float_t         pdfWeight;
Float_t         pthat;
Float_t         processID;
Float_t         genWeight;
Float_t         genHT;
TString         *EventTag;
Int_t           nPUInfo;
vector<int>     *nPU;
vector<int>     *puBX;
vector<float>   *puTrue;
Int_t           nMC;
vector<int>     *mcPID;
vector<float>   *mcVtx;
vector<float>   *mcVty;
vector<float>   *mcVtz;
vector<float>   *mcPt;
vector<float>   *mcMass;
vector<float>   *mcEta;
vector<float>   *mcPhi;
vector<float>   *mcE;
vector<float>   *mcEt;
vector<int>     *mcGMomPID;
vector<int>     *mcMomPID;
vector<float>   *mcMomPt;
vector<float>   *mcMomMass;
vector<float>   *mcMomEta;
vector<float>   *mcMomPhi;
vector<unsigned short> *mcStatusFlag;
vector<int>     *mcParentage;
vector<int>     *mcStatus;
vector<float>   *mcCalIsoDR03;
vector<float>   *mcTrkIsoDR03;
vector<float>   *mcCalIsoDR04;
vector<float>   *mcTrkIsoDR04;
Float_t         genMET;
Float_t         genMETPhi;
Int_t           metFilters;
Float_t         pfMET;
Float_t         pfMETPhi;
Float_t         pfMETsumEt;
Float_t         pfMETmEtSig;
Float_t         pfMETSig;
Float_t         pfMET_T1JERUp;
Float_t         pfMET_T1JERDo;
Float_t         pfMET_T1JESUp;
Float_t         pfMET_T1JESDo;
Float_t         pfMET_T1UESUp;
Float_t         pfMET_T1UESDo;
Float_t         pfMETPhi_T1JESUp;
Float_t         pfMETPhi_T1JESDo;
Float_t         pfMETPhi_T1UESUp;
Float_t         pfMETPhi_T1UESDo;
Int_t           nPho;
vector<float>   *phoE;
vector<float>   *phoEt;
vector<float>   *phoEta;
vector<float>   *phoPhi;
vector<float>   *phoCalibE;
vector<float>   *phoCalibEt;
vector<float>   *phoSCE;
vector<float>   *phoSCRawE;
vector<float>   *phoESEn;
vector<float>   *phoESEnP1;
vector<float>   *phoESEnP2;
vector<float>   *phoSCEta;
vector<float>   *phoSCPhi;
vector<float>   *phoSCEtaWidth;
vector<float>   *phoSCPhiWidth;
vector<float>   *phoSCBrem;
vector<int>     *phohasPixelSeed;
vector<int>     *phoEleVeto;
vector<float>   *phoR9;
vector<float>   *phoHoverE;
vector<float>   *phoE1x3;
vector<float>   *phoE1x5;
vector<float>   *phoE2x2;
vector<float>   *phoE2x5Max;
vector<float>   *phoE5x5;
vector<float>   *phoESEffSigmaRR;
vector<float>   *phoSigmaIEtaIEtaFull5x5;
vector<float>   *phoSigmaIEtaIPhiFull5x5;
vector<float>   *phoSigmaIPhiIPhiFull5x5;
vector<float>   *phoE1x3Full5x5;
vector<float>   *phoE1x5Full5x5;
vector<float>   *phoE2x2Full5x5;
vector<float>   *phoE2x5MaxFull5x5;
vector<float>   *phoE5x5Full5x5;
vector<float>   *phoR9Full5x5;
vector<float>   *phoPFChIso;
vector<float>   *phoPFPhoIso;
vector<float>   *phoPFNeuIso;
vector<float>   *phoPFChWorstIso;
vector<float>   *phoCITKChIso;
vector<float>   *phoCITKPhoIso;
vector<float>   *phoCITKNeuIso;
vector<float>   *phoIDMVA;
vector<unsigned int> *phoFiredSingleTrgs;
vector<unsigned int> *phoFiredDoubleTrgs;
vector<unsigned int> *phoFiredL1Trgs;
vector<float>   *phoSeedTime;
vector<float>   *phoSeedEnergy;
vector<unsigned short> *phoxtalBits;
vector<unsigned short> *phoIDbit;
Int_t           nEle;
vector<int>     *eleCharge;
vector<int>     *eleChargeConsistent;
vector<float>   *eleEn;
vector<float>   *eleSCEn;
vector<float>   *eleESEn;
vector<float>   *eleESEnP1;
vector<float>   *eleESEnP2;
vector<float>   *eleD0;
vector<float>   *eleDz;
vector<float>   *eleSIP;
vector<float>   *elePt;
vector<float>   *eleEta;
vector<float>   *elePhi;
vector<float>   *eleR9;
vector<float>   *eleCalibPt;
vector<float>   *eleCalibEn;
vector<float>   *eleSCEta;
vector<float>   *eleSCPhi;
vector<float>   *eleSCRawEn;
vector<float>   *eleSCEtaWidth;
vector<float>   *eleSCPhiWidth;
vector<float>   *eleHoverE;
vector<float>   *eleEoverP;
vector<float>   *eleEoverPout;
vector<float>   *eleEoverPInv;
vector<float>   *eleBrem;
vector<float>   *eledEtaAtVtx;
vector<float>   *eledPhiAtVtx;
vector<float>   *eledEtaAtCalo;
vector<float>   *eleSigmaIEtaIEtaFull5x5;
vector<float>   *eleSigmaIPhiIPhiFull5x5;
vector<int>     *eleConvVeto;
vector<int>     *eleMissHits;
vector<float>   *eleESEffSigmaRR;
vector<float>   *elePFChIso;
vector<float>   *elePFPhoIso;
vector<float>   *elePFNeuIso;
vector<float>   *elePFPUIso;
vector<float>   *elePFClusEcalIso;
vector<float>   *elePFClusHcalIso;
vector<float>   *elePFMiniIso;
vector<float>   *eleIDMVA;
vector<float>   *eleIDMVAHZZ;
vector<float>   *eledEtaseedAtVtx;
vector<float>   *eleE1x5;
vector<float>   *eleE2x5;
vector<float>   *eleE5x5;
vector<float>   *eleE1x5Full5x5;
vector<float>   *eleE2x5Full5x5;
vector<float>   *eleE5x5Full5x5;
vector<float>   *eleR9Full5x5;
vector<int>     *eleEcalDrivenSeed;
vector<float>   *eleDr03EcalRecHitSumEt;
vector<float>   *eleDr03HcalDepth1TowerSumEt;
vector<float>   *eleDr03HcalDepth2TowerSumEt;
vector<float>   *eleDr03HcalTowerSumEt;
vector<float>   *eleDr03TkSumPt;
vector<float>   *elecaloEnergy;
vector<float>   *eleTrkdxy;
vector<float>   *eleKFHits;
vector<float>   *eleKFChi2;
vector<float>   *eleGSFChi2;
vector<vector<float> > *eleGSFPt;
vector<vector<float> > *eleGSFEta;
vector<vector<float> > *eleGSFPhi;
vector<vector<float> > *eleGSFCharge;
vector<vector<int> > *eleGSFHits;
vector<vector<int> > *eleGSFMissHits;
vector<vector<int> > *eleGSFNHitsMax;
vector<vector<float> > *eleGSFVtxProb;
vector<vector<float> > *eleGSFlxyPV;
vector<vector<float> > *eleGSFlxyBS;
vector<vector<float> > *eleBCEn;
vector<vector<float> > *eleBCEta;
vector<vector<float> > *eleBCPhi;
vector<vector<float> > *eleBCS25;
vector<vector<float> > *eleBCS15;
vector<vector<float> > *eleBCSieie;
vector<vector<float> > *eleBCSieip;
vector<vector<float> > *eleBCSipip;
vector<unsigned int> *eleFiredSingleTrgs;
vector<unsigned int> *eleFiredDoubleTrgs;
vector<unsigned int> *eleFiredL1Trgs;
vector<unsigned short> *eleIDbit;
Int_t           nMu;
vector<float>   *muPt;
vector<float>   *muEn;
vector<float>   *muEta;
vector<float>   *muPhi;
vector<int>     *muCharge;
vector<int>     *muType;
vector<unsigned short> *muIDbit;
vector<float>   *muD0;
vector<float>   *muDz;
vector<float>   *muSIP;
vector<float>   *muChi2NDF;
vector<float>   *muInnerD0;
vector<float>   *muInnerDz;
vector<int>     *muTrkLayers;
vector<int>     *muPixelLayers;
vector<int>     *muPixelHits;
vector<int>     *muMuonHits;
vector<int>     *muStations;
vector<int>     *muMatches;
vector<int>     *muTrkQuality;
vector<float>   *muIsoTrk;
vector<float>   *muPFChIso;
vector<float>   *muPFPhoIso;
vector<float>   *muPFNeuIso;
vector<float>   *muPFPUIso;
vector<float>   *muPFMiniIso;
vector<unsigned int> *muFiredTrgs;
vector<unsigned int> *muFiredL1Trgs;
vector<float>   *muInnervalidFraction;
vector<float>   *musegmentCompatibility;
vector<float>   *muchi2LocalPosition;
vector<float>   *mutrkKink;
vector<float>   *muBestTrkPtError;
vector<float>   *muBestTrkPt;
Int_t           nTau;
vector<bool>    *taupfTausDiscriminationByDecayModeFinding;
vector<bool>    *taupfTausDiscriminationByDecayModeFindingNewDMs;
vector<bool>    *tauByMVA6VLooseElectronRejection;
vector<bool>    *tauByMVA6LooseElectronRejection;
vector<bool>    *tauByMVA6MediumElectronRejection;
vector<bool>    *tauByMVA6TightElectronRejection;
vector<bool>    *tauByMVA6VTightElectronRejection;
vector<bool>    *tauByLooseMuonRejection3;
vector<bool>    *tauByTightMuonRejection3;
vector<bool>    *tauByLooseCombinedIsolationDeltaBetaCorr3Hits;
vector<bool>    *tauByMediumCombinedIsolationDeltaBetaCorr3Hits;
vector<bool>    *tauByTightCombinedIsolationDeltaBetaCorr3Hits;
vector<float>   *tauCombinedIsolationDeltaBetaCorrRaw3Hits;
vector<float>   *tauByIsolationMVArun2v1DBnewDMwLTraw;
vector<float>   *tauByIsolationMVArun2v1DBoldDMwLTraw;
vector<float>   *tauByIsolationMVArun2v1PWnewDMwLTraw;
vector<float>   *tauByIsolationMVArun2v1PWoldDMwLTraw;
vector<bool>    *tauByVTightIsolationMVArun2v1DBnewDMwLT;
vector<bool>    *tauByVTightIsolationMVArun2v1DBoldDMwLT;
vector<bool>    *tauByVTightIsolationMVArun2v1PWnewDMwLT;
vector<bool>    *tauByVTightIsolationMVArun2v1PWoldDMwLT;
vector<bool>    *tauByTightIsolationMVArun2v1DBnewDMwLT;
vector<bool>    *tauByTightIsolationMVArun2v1DBoldDMwLT;
vector<bool>    *tauByTightIsolationMVArun2v1PWnewDMwLT;
vector<bool>    *tauByTightIsolationMVArun2v1PWoldDMwLT;
vector<bool>    *tauByMediumIsolationMVArun2v1DBnewDMwLT;
vector<bool>    *tauByMediumIsolationMVArun2v1DBoldDMwLT;
vector<bool>    *tauByMediumIsolationMVArun2v1PWnewDMwLT;
vector<bool>    *tauByMediumIsolationMVArun2v1PWoldDMwLT;
vector<bool>    *tauByLooseIsolationMVArun2v1DBnewDMwLT;
vector<bool>    *tauByLooseIsolationMVArun2v1DBoldDMwLT;
vector<bool>    *tauByLooseIsolationMVArun2v1PWnewDMwLT;
vector<bool>    *tauByLooseIsolationMVArun2v1PWoldDMwLT;
vector<bool>    *tauByVLooseIsolationMVArun2v1DBnewDMwLT;
vector<bool>    *tauByVLooseIsolationMVArun2v1DBoldDMwLT;
vector<bool>    *tauByVLooseIsolationMVArun2v1PWnewDMwLT;
vector<bool>    *tauByVLooseIsolationMVArun2v1PWoldDMwLT;
vector<float>   *tauEta;
vector<float>   *tauPhi;
vector<float>   *tauPt;
vector<float>   *tauEt;
vector<float>   *tauCharge;
vector<float>   *tauP;
vector<float>   *tauPx;
vector<float>   *tauPy;
vector<float>   *tauPz;
vector<float>   *tauVz;
vector<float>   *tauEnergy;
vector<float>   *tauMass;
vector<float>   *tauDxy;
vector<float>   *tauZImpact;
vector<int>     *tauDecayMode;
vector<bool>    *tauLeadChargedHadronExists;
vector<float>   *tauLeadChargedHadronEta;
vector<float>   *tauLeadChargedHadronPhi;
vector<float>   *tauLeadChargedHadronPt;
vector<float>   *tauChargedIsoPtSum;
vector<float>   *tauNeutralIsoPtSum;
vector<float>   *tauPuCorrPtSum;
vector<int>     *tauNumSignalPFChargedHadrCands;
vector<int>     *tauNumSignalPFNeutrHadrCands;
vector<int>     *tauNumSignalPFGammaCands;
vector<int>     *tauNumSignalPFCands;
vector<int>     *tauNumIsolationPFChargedHadrCands;
vector<int>     *tauNumIsolationPFNeutrHadrCands;
vector<int>     *tauNumIsolationPFGammaCands;
vector<int>     *tauNumIsolationPFCands;
vector<float>   *taufootprintCorrection;
vector<float>   *tauphotonPtSumOutsideSignalCone;
vector<float>   *taudz;
vector<float>   *taudxy;
Int_t           nJet;
vector<float>   *jetPt;
vector<float>   *jetEn;
vector<float>   *jetEta;
vector<float>   *jetPhi;
vector<float>   *jetRawPt;
vector<float>   *jetRawEn;
vector<float>   *jetMt;
vector<float>   *jetArea;
vector<float>   *jetLeadTrackPt;
vector<float>   *jetLeadTrackEta;
vector<float>   *jetLeadTrackPhi;
vector<int>     *jetLepTrackPID;
vector<float>   *jetLepTrackPt;
vector<float>   *jetLepTrackEta;
vector<float>   *jetLepTrackPhi;
vector<float>   *jetCSV2BJetTags;
vector<float>   *jetJetProbabilityBJetTags;
vector<float>   *jetpfCombinedMVAV2BJetTags;
vector<int>     *jetPartonID;
vector<int>     *jetHadFlvr;
vector<float>   *jetGenJetEn;
vector<float>   *jetGenJetPt;
vector<float>   *jetGenJetEta;
vector<float>   *jetGenJetPhi;
vector<int>     *jetGenPartonID;
vector<float>   *jetGenEn;
vector<float>   *jetGenPt;
vector<float>   *jetGenEta;
vector<float>   *jetGenPhi;
vector<int>     *jetGenPartonMomID;
vector<float>   *jetP4Smear;
vector<float>   *jetP4SmearUp;
vector<float>   *jetP4SmearDo;
vector<bool>    *jetPFLooseId;
vector<int>     *jetID;
vector<float>   *jetPUID;
vector<float>   *jetJECUnc;
vector<unsigned int> *jetFiredTrgs;
vector<float>   *jetCHF;
vector<float>   *jetNHF;
vector<float>   *jetCEF;
vector<float>   *jetNEF;
vector<int>     *jetNCH;
vector<int>     *jetNNP;
vector<float>   *jetMUF;
vector<float>   *jetVtxPt;
vector<float>   *jetVtxMass;
vector<float>   *jetVtxNtrks;
vector<float>   *jetVtx3DVal;
vector<float>   *jetVtx3DSig;

#endif // #ifdef CODEX_ANALYZER_H
