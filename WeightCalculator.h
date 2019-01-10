#include <TLorentzVector.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"  // M_PI is in TMath
#include "TRandom3.h"
#include "TSystem.h"
#include "TTree.h"
using std::string;

float XSection(string name) {
  // https://docs.google.com/spreadsheets/d/1rWM3AlFKO8IJVaeoQkWZYWwSvicQ1QCXYSzH74QyZqE/edit?alt=json#gid=398123591

  // I could use a switch, but meh.
  if (name.find("DYJetsToLL_M-50_Inc") != string::npos) {
    return 5765.4;
  } else if (name.find("TTbar") != string::npos) {
    return 831.76;
  } else if (name.find("WJetsToLNu_Inc") != string::npos) {
    return 61526.;
  } else if (name.find("WW") != string::npos) {
    return 115.;
  } else if (name.find("WZ") != string::npos) {
    return 47.13;
  } else if (name.find("ZZ") != string::npos) {
    return 16.523;
  }
}

float weightCalc(TH1F *Histo, string name) {
  if (name.find("SingleElectron") != string::npos || name.find("SingleMuon") != string::npos) {
    return 1.;
  }

  float luminosity(35900.);
  if (Histo->GetBinContent(2) > 0) {
    return luminosity * XSection(name) * 1.0 / Histo->GetBinContent(2);
  } else {
    cerr << "No content in hEvents" << endl;
    return 0;
  }
}

float TMass_F(float pt3lep, float px3lep, float py3lep, float met, float metPhi) {
  return sqrt(pow(pt3lep + met, 2) - pow(px3lep + met * cos(metPhi), 2) - pow(py3lep + met * sin(metPhi), 2));
}
