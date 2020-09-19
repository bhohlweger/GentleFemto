# GentleFemto

This repository contains software to compute correlation functions from experimental data and theory using CATS

Befor the installation make sure
- to have ROOT6 installed (you can use the version that comes with AliROOT just do not enter the environment, see below for workaround) 
- to have CATS installed (https://github.com/dimihayl/DLM)
- to have RooUnfold installed (the AliROOT version is too old and not sufficient - https://gitlab.cern.ch/RooUnfold/RooUnfold) 

Global Variables best defined bashrc/bash_alias/profile:
- It is important to follow the naming scheme since these will be resolved by the CMake of GentleFemto
- CATS Related: have a global variable pointing to your CATS install directory (where the subfolder bin, include and lib is)
  - export CATS=/path/to/CATS_INSTALL/
  - append the CATS libraries to your LD_LIBRARY_PATH, i.e.

    export LD_LIBRARY_PATH="/path/to/CATS_INSTALL/lib:$LD_LIBRARY_PATH" or

    export LD_LIBRARY_PATH="${CATS}/lib:$LD_LIBRARY_PATH"

- RooUnfold related: have a global variable pointing to your RooUnfold install directory (libraries and includes should end up in one directory)
  - export ROOUNFOLD_ROOT="/path/to/RooUnfold_buildorinstall/"
  - append the RooUnfold libraries to your LD_LIBRARY_PATH, i.e.

    export LD_LIBRARY_PATH="/path/to/RooUnfold_buildorinstall/:$LD_LIBRARY_PATH" or

    export LD_LIBRARY_PATH="${ROOUNFOLD_ROOT}:$LD_LIBRARY_PATH"

# Read before installing anything with ROOT

Unfortunately the RooUnfold in AliROOT is missing several methods and is diverged from the actual RooUnfold master for quite a while. If you enter the AliROOT environment you will automatically load this version and paths will be fixed and set. Therefore we cannot do that and we need a 'clean' shell. In principal one can install a standalone version of ROOT, the installation that comes with AliROOT, however, is also fine and most probably already accessible to most users.
Independent of which way you choose, before using ROOT the 'thisroot.sh' of the corresponding installation needs to be loaded. In case of an AliROOT installation, the location of this file varies a bit depending on which system you have, the times you updated your alidist, the weather and the general mood of aliBuild. You should try to look for something like this: 

alice/sw/ubuntu1804_x86-64/ROOT/v6-18-04-1/bin/thisroot.sh

Else you can just briefly enter your AliROOT environment (make sure to not use this terminal to compile anything) and just 'echo ${ROOTSYS}'. It is recommended to define a seperate alias to load ROOT e.g.

alias justROOT='. /alice/sw/ubuntu1804_x86-64/ROOT/v6-18-04-1/bin/thisroot.sh'

In order to avoid conflicts with AliROOT it can then be used to load ROOT alone in case it is needed.

