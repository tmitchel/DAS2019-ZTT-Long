# ZTT Long Exercise - DAS 2019 @ CMS LPC

This repository contains a reference implementation of the analyzers necessary for the ZTT long exercise at CMSDAS 2019. Students should be writing their own analyzer, but here is a completed version as a reference for the facilitators. All code except for the CombineHarvester package should be included here.

##### Table of Contents
[File Location](#loc) <br/>
[Compiling Analyzers](#compile) <br/>
[Running a single file](#single) <br/>
[Running on a directory](#run) <br/>

<a name='loc' />

## File Location
```
$ eosls /store/user/cmsdas/2019/long_exercises/ZTauTau
DYJetsToLL_M-50_Inc.root
SingleElectron.root
SingleMuon.root
TTbar.root
WJetsToLNu_Inc.root
WW.root
WZ.root
ZZ.root
```

<a name='compile' />

## Compiling Analyzers
Analyzers are compiled using the `Make.sh` script as shown below. The output will be a binary with the name `ZTT_XSection_el.exe`.
```
./Make.sh ZTT_XSection_el.cc
```

## Running a single file
The binary can be used to run over a single file to get a single output file as shown below.
```
./ZTT_XSection.exe TTJets.root root://cmseos.fnal.gov//store/user/cmsdas/2019/long_exercises/ZTauTau/TTbar.root
```

This will run on the file `TTbar.root` and produce an output file named `TTJets.root`

## Running on a directory
A python script is included to run an analyzer on all root files in a given directory shown below
```
python automate.py -e ZTT_XSection_mu.exe -i /store/user/tmitchel/DAS2019-ZTT-Long/ -o test_mu_v1
```
This will run the binary `ZTT_XSection_mu.exe` on all files in the directory `/store/user/tmitchel/DAS2019-ZTT-Long/`. The output file names will be the same as the input file name with the suffix `test_mu_v1` added immediately before `.root`.