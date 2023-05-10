# GasHydrateRecycling
This repository contains the source code used in the numerical simulations presented in the manuscript "Hidden periodic states in gas hydrate systems lead to spontaneous gas release without external triggers", authored by S. Gupta, E.B.-Galerne, C. Schmidt, and L. Rüpke.

## Instructions for installation
* System requirements: 
  * Ububtu LTS 18.04 or higher
  * packages: git, auto-conf, cmake (>=3.13), clang, g++, gcc, gfortran, superlu (dev), lapack, libblas, libboost, metis (dev), parmetis (dev), gmsh, paraview, openmpi 
* Execute following commands in terminal (in given order):
  * mkdir dune_2_8
  * cd mkdir dune_2_8
  * mkdir source
  * cd mkdir source
  * cp /downloads/installer.sh .
  * chmod 755 installer.sh
  * ./installer.sh dune
  * cd dune
  * ./buildmobules.sh
  * cd ../..
  * chmod 755 newproject.sh
  * ./newproject.sh GasHydrateRecycling
    * On prompt for "2) Which modules should this module depend on?", enter: dune-common dune-geometry dune-uggrid dune-grid dune-localfunctions dune-istl dune-typetree dune-functions dune-alugrid dune-pdelab
    * On prompt for "3) Project/Module version?", enter 1
    * On promt for "4) Maintainer's email address?", enter your email address.
  * chmod 755 buildproject.sh
  * ./buildproject.sh GasHydrateRecycling
  * cd GasHydrateRecycling/src/
  * rm -rf GasHydrateRecycling.cc
  * rm -rf CMakeLists.txt
  * cp \_all_source_files_in_repo\_ .
  * cd ../..
  * chmod 755 compile.sh
  * ./compile.sh GasHydrateRecycling

## To run the simulations:
* Execute following commands in terminal (in given order):
  * cd /dune_2_8/GasHydrateRecycling/release-build/src
  * ./problem \_your-user-name\_ \_input-file\_  
    * input files are located in the folder: /dune_2_8/GasHydrateRecycling/src/problem/inputs/cyclicity_study/
    * input files are the files with the extension ".ini". In the execution call, drop the .ini extension.
    * hint on \_your_user_name\_: The main executable looks for the following path: home/\_your_user_name\_/dune_2_8/GasHydrateRecycling/

## Files included in this repo:
* installation files:
  * installation/installer.sh
  * installation/newproject.sh
  * installation/buildproject.sh
  * installation/compile.sh
* source files 
  * duneincludes.hh
  * problem/main.cc
  * problem/initial.hh
  * problem/boundary.hh
  * problem/parameters.hh
  * problem/postprocess.hh
  * problem/inputs/cyclicity_study/scenarioXXX.ini 
    * Input files are set up for scenarios with each combination of the parameters sampled in the cyclicity study presented in this manuscript.
    * These parameters are: 1) permeability K, 2) kinetic rate of hydrate phase change kr, and 3) sediment burial rate vs.
  * operators/operations.hh
  * operators/localoperator.hh
  * operators/timeoperator.hh
  * operators/initialvaluefunction.hh
  * operators/boundaryvaluefunction.hh
  * operators/driver.hh
* outputs/problem/cyclicity_study
  *  Simulation outputs for the reference scenario is archived in this repository.
  *  "readme.txt" gives a quick overview on how to read and plot these.
