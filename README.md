<h1 align="center">Z-lambda software</h1>
<h3 align="center">for investigating heat transfer during z-scan experiments</h3>

Z-lambda is a software tool written in C++ for investigating heat transfer during Z-scan experiments with isotropic transparent solids. The program can calculate the function T(t,x,y,z) of the sample during Z-scan experiments, which makes it a valuable resource for researchers in the field of materials science and optics.

<h3 align="left">Installation</h3>

To use Z-lambda, you must have a C++ compiler installed on your system. The program can be compiled using at least C++11 compiler with pthreads. Simply download the source code and compile it on your system. For lunux and macOS simply execute bash file:

<pre><code class="language-bash">
./compile.sh
</code></pre>

The script will compile the programme and save it in the bin directory, along with the files needed to run it. Manually, compilation can be performed using the command:

<pre><code class="language-bash">
g++ main.cpp ExperimentData.cpp thermo.cpp -o zlambda -pthread -std=c++11 -O3
</code></pre>

<h3 align="left">Usage</h3>

To run Z-lambda, you must provide two input configuration files: ConfigExperiment.txt and ConfigSample.txt. These files specify the experimental parameters and the properties of the sample being tested, respectively.

ConfigExperiment.txt contains the following parameters:
Laser pulse energy
Laser pulse duration
Spot size
Number of pulses
Scanning step
Z-scan range
Scanning speed

ConfigSample.txt contains the following properties of the sample:

Absorption coefficient
Refractive index
Density
Specific heat capacity
Thermal conductivity
Once you have specified these parameters and properties in the input files, you can run the program by executing the following command:

<pre><code class="language-bash">
./zlambda.exe
</code></pre>
./zlambda.exe
The program will then simulate heat transfer using the specified parameters and properties, and output the temperature distribution of the sample during the Z-scan experiment.

<h3 align="left">Output</h3>

Z-lambda outputs the temperature distribution of the sample as a 3D matrix of values. These values can be visualized using any software capable of rendering 3D matrices. Additionally, the program can output a CSV file containing the temperature values for each point in the matrix.

<h3 align="left">Impact and Outlook</h3>

Z-lambda is a valuable resource for researchers investigating heat transfer during Z-scan (and beyond) experiments with isotropically transparent solids. Its versatility, ease of use, and open-ended nature make it an important tool for ongoing research in this field. Ongoing development efforts promise to further enhance its capabilities in the future.


<h3 align="left">Z-lambda source:</h3>
<p align="left"> <a href="https://www.w3schools.com/cpp/" target="_blank" rel="noreferrer"> <img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/cplusplus/cplusplus-original.svg" alt="cplusplus" width="40" height="40"/> </a> </p>