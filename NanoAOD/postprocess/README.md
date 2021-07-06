# Bmm5 post processing scipts

## Introduction

Bmm5 analysis uses many different samples for a varity of studies and
measurements. It is more efficient to have samples optimized for a
specific task than using NanoAOD as universal input. A full set of
NanoAOD samples is of the order of 1TB in size and processing such
volume of data requires a batch system or long processing, which
introduce delays.

There are two main types of post processors:
* Skimmer - used to select a small subset of events keeping NanoAOD format
* FlatNtuple - simple tabular data format targetting a specific study containing small number of variables

## Installation

In order to postprocess (skim) NanoAOD we need https://github.com/cms-nanoAOD/nanoAOD-tools
```shell
cd $CMSSW_BASE/src
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
cd PhysicsTools/NanoAODTools
cmsenv
scram b
```
