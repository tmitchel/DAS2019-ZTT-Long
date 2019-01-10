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
using namespace std;

float XSection(std::string OutName) {
  //    https://docs.google.com/spreadsheets/d/1rWM3AlFKO8IJVaeoQkWZYWwSvicQ1QCXYSzH74QyZqE/edit?alt=json#gid=398123591

  map<string, float> cross_sections = {
      {"WJets", 61526},
      {"DYJetsToLL", 5765.4},
      {"TTJets", 831.76},
      {"WW", 115.0},
      {"WZ", 47.13},
      {"ZZ", 16.523}};

  return cross_sections[OutName];
}

float weightCalc(TH1F *Histo, std::string outputName) {

  if (outputName.find("SingleElectron") != string::npos || outputName.find("SingleMuon") != string::npos) {
    return 1.;
  }

  float luminosity(35900.);
  if (Histo->GetBinContent(2) > 0) {
    return luminosity * XSection(outputName) * 1.0 / Histo->GetBinContent(2);
  } else {
    cerr << "No content in hEvents" << endl;
    return 0;
  }
}

float TMass_F(float pt3lep, float px3lep, float py3lep, float met, float metPhi) {
  return sqrt(pow(pt3lep + met, 2) - pow(px3lep + met * cos(metPhi), 2) - pow(py3lep + met * sin(metPhi), 2));
}
