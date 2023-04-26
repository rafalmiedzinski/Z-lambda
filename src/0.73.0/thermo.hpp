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


#pragma once
#ifndef THERMO_H_INCLUDED
#define THERMO_H_INCLUDED

double OneD_Gaussian_Profile(double,double,double);
double TwoD_Gaussian_Profile(double, double, double, double);
double ThreeD_Gaussian_Profile(double, double, double, double, double);
double EfieldToIrradiance(double, double);
double LambertBeer(double, double, double);
//double CornerThreeD(double, double, double, double, double, double, double[3], double, double, double, double);





class SampleData
{
public:
    
    double GetAlpha(); //ThermalDiffusivity - a or alpha
    double GetXmeters();
    double GetYmeters();
    double GetZmeters();
    double GetdX();
    double GetdY();
    double GetdZ();

    double GetSampleTramsmission();

    double GetSampleReflectance();

    double GetSampleAbsorption();

    double GetSampleSpecificHeatCapacity();

    double GetSampleDensity();
  
    double GetThermalConductivity(); //lambda

    double Getbeta_thermal();

    

    //External convective heat transfer from the axis 
    // X - dXiR (Right) and dXiL (Left)
    // Y - dYiU (Upper) oraz dYiL (Lower)
    // Z - dZiF (Front) oraz dZiB (Back)
    double Getalfa_xR();
    double Getalfa_xL();
    double Getalfa_yU();
    double Getalfa_yL();
    double Getalfa_zF();
    double Getalfa_zB();
    
    //Emissivity epsilon of the sample's walls

    double Getepsilon_xL();
    double Getepsilon_xR();
    double Getepsilon_yU();
    double Getepsilon_yL();
    double Getepsilon_zF();
    double Getepsilon_zB();

    void SayHello();
    void MyVersion();

    double GetTamb_zF();
    double GetTamb_zB();

    double GetTamb_xL();
    double GetTamb_xR();

    double GetTamb_yU();
    double GetTamb_yD();


    unsigned int GetN();
    unsigned int GetO();
    unsigned int GetP();

    double GetTemp0();

    bool InitiateError = false;

    void ShowMeASample();
    SampleData();
    ~SampleData();

private:
    double Density = 1;                     //in [kg / m^3]
    bool isDensity = false;

    double SpecificHeatCapacity = 1;        //in [J/Kg K]
    bool isSpecificHeatCapacity = false;

    double ThermalConductivity = 1;         //in [W / K m] aka lambda
    bool isThermalConductivity = false;

    double ThermalDiffusivity = 1;           //in [m^2 / s] aka alfa
    bool isThermalDiffusivity = false;

    double SampleTrans = 1.0;               //Sample light transmission
    bool isSampleTrans = false;

    double SampleReflectance = 0.01;        //Sample light reflectance
    bool isSampleReflectance = false;

    double beta_thermal = 1.0;
    bool isbeta_thermal = false;

    double dimensionX = 1;                  //in [mm]
    bool isdimensionX = false;

    double dimensionY = 1;                  //in [mm]
    bool isdimensionY = false;

    double dimensionZ = 1;                  //in [mm]
    bool isdimensionZ = false;

    double dX = 1;                  //in [m] Space steps dx, dy, dz of the single node
    double dY = 1;                  //in [m]
    double dZ = 1;                  //in [m]

    double Temp0 = 293.15;          //Sample initial temperature [K]
    bool isTemp0 = false;

    unsigned int N = 10;
    bool isN = false;

    unsigned int O = 10;
    bool isO = false;

    unsigned int P = 10;
    bool isP = false;

    double Tamb_zF = 293.15;                  //in [K]
    bool isTamb_zF = false;

    double Tamb_zB = 293.15;                  //in [K]
    bool isTamb_zB = false;

    double Tamb_xR = 293.15;                  //in [K]
    bool isTamb_xR = false;

    double Tamb_xL = 293.15;                  //in [K]
    bool isTamb_xL = false;

    double Tamb_yU = 293.15;                  //in [K]
    bool isTamb_yU = false;

    double Tamb_yD = 293.15;                  //in [K]
    bool isTamb_yD = false;

    double alfa_zF = 5.0;                  //in  [W / (m^2*K)]
    bool isalfa_zF = false;

    double alfa_zB = 5.0;                  //in  [W / (m^2*K)]
    bool isalfa_zB = false;

    double alfa_xR = 5.0;                  //in  [W / (m^2*K)]
    bool isalfa_xR = false;

    double alfa_xL = 5.0;                  //in  [W / (m^2*K)]
    bool isalfa_xL = false;

    double alfa_yU = 5.0;                  //in  [W / (m^2*K)]
    bool isalfa_yU = false;

    double alfa_yL = 5.0;                  //in  [W / (m^2*K)]
    bool isalfa_yL = false;
    
    double epsilon_zF = 1.00;              //Thermal emissivity in [1]
    bool isepsilon_zF = false;
    
    double epsilon_zB = 1.00;              //Thermal emissivity in [1]
    bool isepsilon_zB = false;
    
    double epsilon_xL = 1.00;              //Thermal emissivity in [1]
    bool isepsilon_xL = false;
    
    double epsilon_xR = 1.00;              //Thermal emissivity in [1]
    bool isepsilon_xR = false;
    
    double epsilon_yU = 1.00;              //Thermal emissivity in [1]
    bool isepsilon_yU = false;
    
    double epsilon_yL = 1.00;              //Thermal emissivity in [1]
    bool isepsilon_yL = false;
    

    
    
};


template <class T> T ***Create3D(unsigned int N1, unsigned int N2, unsigned int N3)
{
    T *** array = new T ** [N1];

    array[0] = new T * [N1*N2];

    array[0][0] = new T [N1*N2*N3];

    unsigned int i,j;

    for( i = 0; i < N1; i++) {

        if (i < N1 -1 ) {

            array[0][(i+1)*N2] = &(array[0][0][(i+1)*N3*N2]);

            array[i+1] = &(array[0][(i+1)*N2]);

        }

        for( j = 0; j < N2; j++) {     
            if (j > 0) array[i][j] = array[i][j-1] + N3;
        }

    }

    //cout << endl;
    return array;
};

template <class T> void Delete3D(T ***array) {
    delete[] array[0][0]; 
    delete[] array[0];
    delete[] array;
};

#endif // thermo_H_INCLUDED
