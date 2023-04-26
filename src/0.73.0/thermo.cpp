//  Z-lambda Software
//  Copyright(c) 2020 - 2023 by 
//  Rafał Miedziński (http://orcid.org/0000-0002-7526-089X) 
//  and
//  Izabela Fuks - Janczarek (http://orcid.org/0000-0002-5129-6783)
// 
//  MIT License
//
//  Created by Rafal Miedzinski & Izabela Fuks-Janczarek on 24/01/2020.
//

#include <cmath>
#include <fstream>
#include<iostream>
#include<string>
#include <algorithm>
#include <stdexcept>
#include"thermo.hpp"


using namespace std;

//using std::remove_if;
//q0 = 2^(5/2)*Beta*P / pi^(3/2) * r0^3

void SampleData::MyVersion()
{
    string major = " 0";
    string minor = "73";
    string path  = "0 ";

    cout << "<--                         Z-lambda(C) ver."
    <<major<<"."<<minor<<"."<<path
    <<"                         -->\n";
}

double SampleData::GetAlpha()
{
    if ((Density > 0) && (SpecificHeatCapacity > 0))
    {
        isThermalDiffusivity = true;
        return ThermalConductivity / (Density * SpecificHeatCapacity);
    }
    isThermalDiffusivity = false;
    return 0.0;
}

double SampleData::GetXmeters()
{
    return (dimensionX/1000);
}

double SampleData::GetYmeters()
{
    return (dimensionY/1000);
}

double SampleData::GetZmeters()
{
    return (dimensionZ/1000);
}

double SampleData::GetdX()
{
    //Function returs value in (m)
    dX = (dimensionX / 1000) / (N - 1);
    return dX;
}

double SampleData::GetdY()
{
    //Function returs value in (m)
    dY = (dimensionY / 1000) / (O - 1);
    return dY;
}

double SampleData::GetdZ()
{
    //Function returs value in (m)
    dZ = (dimensionZ / 1000) / (P - 1);
    return dZ;
}

double SampleData::GetSampleTramsmission()
{
    return SampleTrans;
}

double SampleData::GetSampleReflectance()
{
    return SampleReflectance;
}

double SampleData::GetSampleAbsorption()
{
    //Returns the sample absorption coefficient in [m^-1]
    
    if ((GetZmeters() > 0 ) && isSampleReflectance && (GetSampleTramsmission() > 0))
    {
        double Rsquare = GetSampleReflectance() * GetSampleReflectance();
        return (1 / GetZmeters()) * log((1 - Rsquare) / GetSampleTramsmission());
    }
    
    return 0.0;
}

double SampleData::GetSampleSpecificHeatCapacity()
{
    //       J/(kg * K)
    return SpecificHeatCapacity;
}

double SampleData::GetSampleDensity()
{
    // kg/m^3
    return Density;
}

double SampleData::GetThermalConductivity()
{
    return ThermalConductivity;
}

double SampleData::Getbeta_thermal()
{
    return beta_thermal;
}


double SampleData::Getalfa_xR()
{
    return alfa_xR;
}

double SampleData::Getalfa_xL()
{
    return alfa_xL;
}

double SampleData::Getalfa_yU()
{
    return alfa_yU;
}

double SampleData::Getalfa_yL()
{
    return alfa_yL;
}

double SampleData::Getalfa_zF()
{
    return alfa_zF;
}

double SampleData::Getalfa_zB()
{
    return alfa_zB;
}

double SampleData::GetTamb_zF()
{
    return Tamb_zF;
}

double SampleData::GetTamb_zB()
{
    return Tamb_zB;
}

double SampleData::GetTamb_xL()
{
    return Tamb_xL;
}

double SampleData::GetTamb_xR()
{
    return Tamb_xR;
}

double SampleData::GetTamb_yU()
{
    return Tamb_yU;
}

double SampleData::GetTamb_yD()
{
    return Tamb_yD;
}

unsigned int SampleData::GetN()
{
    return N;
}

unsigned int SampleData::GetO()
{
    return O;
}

unsigned int SampleData::GetP()
{
    return P;
}

double SampleData::GetTemp0()
{
    return Temp0;
}

double SampleData::Getepsilon_xL()
{
    return epsilon_xL;
}

double SampleData::Getepsilon_xR()
{
    return epsilon_xR;
}

double SampleData::Getepsilon_yU()
{
    return epsilon_yU;
}

double SampleData::Getepsilon_yL()
{
    return epsilon_yL;
}

double SampleData::Getepsilon_zF()
{
    return epsilon_zF;
}

double SampleData::Getepsilon_zB()
{
    return epsilon_zB;
}

void SampleData::SayHello()
{
    cout << "<------------------------------------------------------------------------------>" << "\n";
    MyVersion();
    //cout << "<--                         Z-lambda(C) ver. 0.71.1                          -->" << "\n";
    cout << "<--               https://github.com/rafalmiedzinski/Z-lambda                -->" << "\n";
    cout << "<--                  http://www.fizyka.ujd.edu.pl/zlambda                    -->" << "\n";
    cout << "<--       (c) 2020 - 2023 Rafal Miedzinski and Izabela Fuks-Janczarek        -->" << "\n";
    cout << "<--                                MIT License                               -->" << "\n";
    cout << "<--                                                                          -->" << "\n";
    cout << "<--                                                                          -->" << "\n";
    cout << "<-- Software Authors:                                                        -->" << "\n";
    cout << "<--      Rafal Miedzinski        (http://orcid.org/0000-0002-7526-089X)      -->" << "\n";
    cout << "<--      Izabela Fuks-Janczarek  (http://orcid.org/0000-0002-5129-6783)      -->" << "\n";
    cout << "<--                                                                          -->" << "\n";
    cout << "<------------------------------------------------------------------------------>" << "\n \n";

}


void SampleData::ShowMeASample()
{
    SayHello();
    GetAlpha();
    cout << "<------------------------------------------------------------------------------>" << "\n";
    cout << "<----------------     Physical parameters of the sample:     ------------------>" << "\n";
    cout << "<------------------------------------------------------------------------------>" << "\n \n";
    cout << " Density:                         rho = " << Density << " [kg / m^3]" << "\n";
    cout << " Specific heat capacity:           cw = " << SpecificHeatCapacity << " [J / kg K]" << "\n";
    cout << " Thermal conductivity:         lambda = " << ThermalConductivity << " [W / m K]" << "\n";
    cout << " Thermal l.beam absorptivity: beta_th = " << beta_thermal << " [1]" << "\n";
    cout << " Thermal diffusivity:           alpha = " << GetAlpha() << " [m^2 / s]" << "\n";
    cout << " Sample initial temperature:       T0 = " << GetTemp0() << " [K]" << "\n";
    cout << " Sample light transmission         Tl = " << GetSampleTramsmission() * 100 << " [%]" << "\n";
    cout << " Sample light reflectance          Rl = " << GetSampleReflectance() * 100 << " [%]" << "\n";
    cout << " Sample light absorption           Al = " << GetSampleAbsorption() / 100 << " [cm-1]" << "\n";
    cout << "<------------------------------------------------------------------------------>" << "\n \n";
    cout << " Sample dimension:            (x,y,z) = (" << dimensionX << ", " << dimensionY << ", " << dimensionZ << ") [mm]" << "\n";
    cout << " Number of the nodes:         (N,O,P) = (" << N << "," << O << "," << P << ")" << "\n";
    cout << " Space steps:              (dx,dy,dz) = (" << GetdX()*1000000 << "," << GetdY()*1000000 << "," << GetdZ()*1000000 << ") [um]" << "\n \n";
    cout << "<------------------------------------------------------------------------------>" << "\n \n";
    cout << " Thermal emissivity of the sample's wall:\n\n";
    cout << " Front wall                epsilon_zF = " << Getepsilon_zF() << endl;
    cout << " Back wall                 epsilon_zB = " << Getepsilon_zB() << endl<<endl;
    cout << " Left wall                 epsilon_xL = " << Getepsilon_xL() << endl;
    cout << " Right wall                epsilon_xR = " << Getepsilon_xR() << endl<<endl;
    cout << " Upper wall                epsilon_yU = " << Getepsilon_yU() << endl;
    cout << " Lower wall                epsilon_yL = " << Getepsilon_yL() << endl<<endl;

}

SampleData::SampleData()
{
    std::ifstream cFile("ConfigSample.txt");
    if (cFile.is_open())
    {
        std::string line;
        while (getline(cFile, line)) {
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace),
                line.end());
            if (line[0] == '#' || line.empty())
                continue;
            auto delimiterPos = line.find("=");
            auto name = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);
            try {
                if (name == "Density") {
                    Density = stod(value);
                    isDensity = true;
                }

                if (name == "SpecificHeatCapacity") {
                    SpecificHeatCapacity = stod(value);
                    isSpecificHeatCapacity = true;
                }
                if (name == "ThermalConductivity") { 
                    ThermalConductivity = stod(value);
                    isThermalConductivity = true;
                }

                if (name == "x") {
                    dimensionX = stod(value);
                    isdimensionX = true;
                }
                
                if (name == "y") {
                    dimensionY = stod(value);
                    isdimensionY = true;
                }
                
                if (name == "z") {
                    dimensionZ = stod(value);
                    isdimensionZ = true;
                    isdimensionZ = true;
                } 

                if (name == "alfa_zF") {
                    alfa_zF = stod(value);
                    isalfa_zF = true;
                }

                if (name == "alfa_zB") {
                    alfa_zB = stod(value);
                    isalfa_zB = true;
                }

                if (name == "alfa_xL") {
                    alfa_xL = stod(value);
                    isalfa_xL = true;
                }

                if (name == "alfa_xR") {
                    alfa_xR = stod(value);
                    isalfa_xR = true;
                }

                if (name == "alfa_yU") {
                    alfa_yU = stod(value);
                    isalfa_yU = true;
                }

                if (name == "alfa_yL") {
                    alfa_yL = stod(value);
                    isalfa_yL = true;
                }
                
                if (name == "epsilon_xL") {
                    epsilon_xL = stod(value);
                    isepsilon_xL = true;
                }
                
                if (name == "epsilon_xR") {
                    epsilon_xR = stod(value);
                    isepsilon_xR = true;
                }
                
                if (name == "epsilon_xL") {
                    epsilon_xL = stod(value);
                    isepsilon_xL = true;
                }
                
                if (name == "epsilon_yU") {
                    epsilon_yU = stod(value);
                    isepsilon_yU = true;
                }
                
                if (name == "epsilon_yL") {
                    epsilon_yL = stod(value);
                    isepsilon_yL = true;
                }
                
                if (name == "epsilon_zB") {
                    epsilon_zB = stod(value);
                    isepsilon_zB = true;
                }
                
                if (name == "epsilon_zF") {
                    epsilon_zF = stod(value);
                    isepsilon_zF = true;
                }

                if (name == "Tamb_zF") {
                    Tamb_zF = stod(value);
                    isTamb_zF = true;
                }

                if (name == "Tamb_zB") {
                    Tamb_zB = stod(value);
                    isTamb_zB = true;
                }

                if (name == "Tamb_xL") {
                    Tamb_xL = stod(value);
                    isTamb_xL = true;
                }

                if (name == "Tamb_xR") {
                    Tamb_xR = stod(value);
                    isTamb_xR = true;
                }

                if (name == "Tamb_yU") {
                    Tamb_yU = stod(value);
                    isTamb_yU = true;
                }

                if (name == "Tamb_yL") {
                    Tamb_yD = stod(value);
                    isTamb_yD = true;
                }

                if (name == "N") {
                    N = stoi(value);
                    isN = true;
                }

                if (name == "O") {
                    O = stoi(value);
                    isO = true;
                }
           
                if (name == "P") {
                    P = stoi(value);
                    if (P % 2 == 0) P++;
                    isP = true;
                }

                if (name == "Temp0") {
                    Temp0 = stod(value);
                    isTemp0 = true;
                }

                if (name == "Transmission") {
                    SampleTrans = stod(value);
                    isSampleTrans = true;
                }

                if (name == "Reflectance") {
                    SampleReflectance = stod(value);
                    isSampleReflectance = true;
                }

                if (name == "beta_thermal") {
                    beta_thermal = stod(value);
                    isbeta_thermal = true;
                }

 
            }
            catch (const std::invalid_argument & ia) {
                std::cerr << " Invalid value in: " << name << endl << "Error: " << ia.what() << endl;
                InitiateError = true;
            }
        }
    }
    else {
        std::cerr << " Error: Missed ConfigSample.txt file!\n";
        InitiateError = true;
    }
}

SampleData::~SampleData()
{
}


double OneD_Gaussian_Profile(double q0, double r0, double x)
{
    // q0 - wsp. przed eksponentą
    // r0 - laser beam radius
    //
    return q0*exp((-2*(pow(x, 2.0)))/(pow(r0, 2.0)));
}

double TwoD_Gaussian_Profile(double q0, double r0, double x, double y)
{
    // q0 - wsp. przed eksponentą
    // r0 - laser beam radius
    //
    return q0*exp((-2*(pow(x, 2.0)+pow(y, 2.0)))/(pow(r0, 2.0)));
}

double EfieldToIrradiance(double realE, double n)
{
    // Function converts electric field value (V/m) to power density (W/m^2)
    // realE - real part of electric field
    // imagE - imaginary part of electric field
    // c - speed of light (299792458 m/s)
    // eo - electric constant - (8.854187817e-12 F/m)
    
    return (299792458 * n * 8.854187817e-12 * 0.5 * realE * realE);
    
    // Results in: W/m2
    
}

double LambertBeer(double Io, double z, double alfa)
{
    // This function calculates the absorption of radiation using the Lambert-Beer law along the z 
    // path with absorption coefficient alpha.
    // 
    // Io - input power density (W/m^2)
    // z - length of the path taken by the light (m)
    // alpha - absorption coefficient (m^-1)
   
   return (Io*exp(-alfa*z));
}


double ThreeD_Gaussian_Profile(double omega_zero, double z_zero, double x, double y, double z)
{
    // q0 - coefficent ---->  q0 = [2^(5/2)*Beta*P] / [pi^(3/2) * r0^3]
    // r0 - laser beam radius
    //
    
    double omega_z_kwadrat = omega_zero*omega_zero*(1+pow((z/z_zero),2.0));
    return (omega_zero/sqrt(omega_z_kwadrat))*exp(- (x*x+y*y)/(omega_z_kwadrat) );
    
    
    //return q0*exp((-2*(pow(x, 2.0)+pow(y, 2.0)+pow(z, 2.0)))/(pow(r0, 2.0)));
}
