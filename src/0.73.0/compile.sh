#!/bin/bash
echo "Removing the previous version of the program along with the results...."
rm -R bin/
echo "Creating a new directory structure..."
mkdir bin
mkdir bin/Results
mkdir bin/Results/Surfaces
mkdir bin/Results/HeatMaps
cp ConfigSample.txt bin/
cp ConfigExperiment.txt bin/
echo "Compiling and linking z-lambda software..."
g++ main.cpp ExperimentData.cpp thermo.cpp -o zlambda -pthread -std=c++11 -O3
echo "Moving the program to the <bin> directory..."
mv zlambda bin/
