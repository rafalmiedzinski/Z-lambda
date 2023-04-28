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
#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <thread>
#include <ctime>
#include <cstdlib>
#include "thermo.hpp"
#include "ExperimentData.hpp"


using namespace std;


ExperimentData::ExperimentData()
{
    std::ifstream cFile("ConfigExperiment.txt");
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
                if (name == "SimulationTime") {
                    SimulationTime = stod(value);
                    isSimulationTime = true;
                }
                
                if (name == "NumericalTimeStep") {
                    NumericalTimeStep = stod(value);
                    isNumericalTimeStep = true;
                }

                if (name == "LambdaOfLaser") {
                    LambdaOfLaser = stod(value);
                    isLambdaOfLaser = true;
                }

                if (name == "LaserPulseDuration") {
                    LaserPulseDuration = stod(value);
                    isLaserPulseDuration = true;
                    if (LaserPulseDuration == 0)
                    {
                        isCWLaser = true;
                        isPulsedLaser = false;
                    }

                    if (LaserPulseDuration > 0)
                    {
                        isCWLaser = false;
                        isPulsedLaser = true;
                    }
                }

                if (name == "LaserPulseRepetition") {
                    LaserPulseRepetition = stoi(value);
                    isLaserPulseRepetition = true;
                }

                if (name == "LaserCWPower") {
                    LaserCWPower = stod(value);
                    isLaserCWPower = true;
                }

                if (name == "LaserPulseEnergy") {
                    LaserPulseEnergy = stod(value);
                    isLaserPulseEnergy = true;
                }

                if (name == "ZScanRange") {
                    ZScanRange = stod(value);
                    ZScanRange /= 1000;
                    isZScanRange = true;
                }

                if (name == "SampleVelocity") {
                    SampleVelocity = stod(value);
                    isSampleVelocity = true;
                }

                if (name == "SingleShot") {
                    SingleShot = stoi(value);
                }

                if (name == "omegaZero") {
                    omegaZero = stod(value);
                    isomegaZero = true;
                }

                if (name == "SaveZPRInterval") {
                    SaveZPResultsInterval = stod(value);
                    isSaveZPResultsInterval = true;
                }

                if (name == "DoesSaveZPRInterval") {
                    DoSaveZeroPlane = stoi(value);
                    isDoSaveZeroPlane = true;
                }

                if (name == "SaveASInterval") {
                    SaveAllResultsInterval = stod(value);
                    isSaveAllResultsInterval = true;
                }

                if (name == "DoesSaveASInterval") {
                    DoSaveAllResults = stoi(value);
                    isDoSaveAllResults = true;
                }

                if (name == "Threads") {
                    Threads = stoi(value);
                    isThreads = true;
                }

            }
            catch (const std::invalid_argument & ia) {
                std::cerr << "Invalid value in: " << name << endl << "Error: " << ia.what() << endl;
                InitiateError = true;
            }
        }

    }
    else {
        std::cerr << "Error: Missed ConfigExperiment.txt file!\n";
        InitiateError = true;
    }
    setdeltat();
}

ExperimentData::~ExperimentData()
{
}

void ExperimentData::ShowMeExperiment()
{
    //Calculation of omega(z)
    double omegaZeroSquare = GetomegaZero() * GetomegaZero();
    double Zr = ((3.141592653589 * GetomegaZero() * GetomegaZero()) / (GetLaserWavelength() * 1e-9));
    double position = (GetZScanRange() / 2);
    double omegasquare = omegaZeroSquare * (1 + (position/ Zr) * (position / Zr));

    cout << "<------------------------------------------------------------------------------>" << "\n";
    cout << "<----------------        Parameters of the experiment         ----------------->" << "\n";
    cout << "<------------------------------------------------------------------------------>" << "\n \n";
    cout << " Simulation Time:                 t = " << SimulationTime << " [s]";
    cout << "\n Number of time iteration:       iT = " << GetTimeSteps();
    cout << "\n Time step duration:             dT = " << Getdeltat() << " [s]\n\n";
    cout << "<------------------------------------------------------------------------------>" << "\n";
    cout << "<----------------               LASER parameters              ----------------->" << "\n";
    cout << "<------------------------------------------------------------------------------>" << "\n \n";
    cout << " Wavelength:                 lambda = " << LambdaOfLaser << " [nm]";
    cout << "\n Beam diameter:         2 x omega_0 = " << (GetomegaZero() * 1e6) * 2 << " [um]";
    cout << "\n Beam diameter at:    2 x omega(-z) = " << (sqrt(omegasquare) * 1e3) << " [mm]";

    if (isCWLaser)
    {
        cout << "\n Laser type:                        = CW";
        cout << "\n Power of the Laser:             Pl = " << LaserCWPower << " [W]";
        cout << "\n Laser irradiance at z0:         I0 = " << GetLaserIrradianceI0inGWcm2() << " [kW / cm^2]";
    }

    if (isPulsedLaser)
    {
        cout << "\n Laser type:                        = Pulsed";
        cout << "\n Energy per pulse:               Ep = " << LaserPulseEnergy << " [J]";
        cout << "\n Laser irradiance at z0:         I0 = " << GetLaserIrradianceI0inGWcm2() << " [GW / cm^2]";
        cout << "\n Pulse duration:                tau = " << LaserPulseDuration << " [s]";
        cout << "\n Laser repetition frequency:     fr = " << LaserPulseRepetition << " [Hz]";
    }

}

int ExperimentData::GetSingleShot()
{
    return SingleShot;
}

double ExperimentData::GetZScanRange()
{
    //Result in (m)
    return ZScanRange;
}

double ExperimentData::GetSampleVeolocity()
{
    // value in (m/s)
    return (SampleVelocity/1000);
}

double ExperimentData::GetSimulationTime()
{
    return SimulationTime;
}

double ExperimentData::Getdeltat()
{
    return deltat;
}

double ExperimentData::GetdeltatSlowMotion()
{
    return deltatSlowMontion;
}

unsigned int ExperimentData::GetTimeSteps()
{
    if (isSimulationTime && isNumericalTimeStep)
    {
        return ((SimulationTime / NumericalTimeStep));
    }
    else
    {
        return 0;
    }  
}

double ExperimentData::GetLaserWavelength()
{
    return LambdaOfLaser;
}

double ExperimentData::GetLaserPulseDuration()
{
    if (LaserPulseDuration == 0)
    {
        isCWLaser = true;
        isPulsedLaser = false;
        return LaserPulseDuration;
    }
    else
    {
        isCWLaser = false;
        isPulsedLaser = true;
        return LaserPulseDuration;
    }
    
    
    return LaserPulseDuration;
}

unsigned int ExperimentData::GetLaserPulseRepetition()
{
    return LaserPulseRepetition;
}

double ExperimentData::GetLaserPulseEnergy()
{
    return LaserPulseEnergy;
}

double ExperimentData::GetLaserCWPower()
{
    return LaserCWPower;
}

double ExperimentData::GetomegaZero()
{
    //Results in (m)
    return omegaZero;
}

unsigned short int ExperimentData::GetThreads()
{
    if (Threads >= 64) return 64;
    return Threads;
}

double ExperimentData::GetLaserIrradianceI0inGWcm2()
{
    if (isPulsedLaser)
    {
        double Power = LaserPulseEnergy / LaserPulseDuration; //To change J to W
        Power /= 1e+9; // To change W to GW 
        double Area = (3.14159265358979323846 * omegaZero * omegaZero) * 10000; // to change m^2 to cm^2 
        return Power / Area; // results in [GW / cm^2]
    }
    
    if (isCWLaser)
    {
        double Power = LaserCWPower / 1e+3; // To change W to kW (in 1 sec. gives W)
        double Area = (3.14159265358979323846 * omegaZero * omegaZero) * 10000; // to change m^2 to cm^2 
        return Power / Area; // results in [kW / cm^2]
    }

    return 0.0;
}

double ExperimentData::GetSaveZPResultsInterval()
{
    if (isSaveZPResultsInterval) 
    {
        return SaveZPResultsInterval;
    }
    else
    {
        return std::numeric_limits<double>::max();
    }
}

double ExperimentData::GetSaveAllResultsInterval()
{
    if (isSaveAllResultsInterval)
    {
        return SaveAllResultsInterval;
    }
    else
    {
        return std::numeric_limits<double>::max();
    }
}




void ExperimentData::setdeltat()
{
    if (SimulationTime != 0 && NumericalTimeStep != 0) {
        deltat = NumericalTimeStep;
        deltatSlowMontion = LaserPulseDuration / 4; //This time is 1/4 of pulse duration - for better modelling of the laser pulse 
    }
    else
    {
        cout << "Some error occured: missing SimulationTime or/and TimeSteps in ConfigExperiment.txt file!\n";
        InitiateError = true;
    }
        
}

double SimulationData::GetOmega(double z)
{
    double omegasquare = omegaZeroSquare * (1 + ((z * z)/(GetZrSquare())));
    return sqrt(omegasquare);
}

bool SimulationData::CheckCalculationStability()
{
    return false;
}


void SimulationData::SetData(ExperimentData incomingexData, SampleData incomingsmData)
{
    
    exData = &incomingexData;
    spData = &incomingsmData;

    spData->ShowMeASample();
    exData->ShowMeExperiment();

    InitaiateCenterProfiles();           //Erase X,Y ZProfile.dat files
    InitiateCenterHeaatMapSurfaces(0);   //Erase CenterSurfaceXY.dat, CenterSurfaceXZ.dat, CenterSurfaceYZ.dat files
    InitiateSampleEnergyFile();          //Erase data in Results/SampleEnergy.dat
    
    if (incomingexData.GetSingleShot() == 1)
    {
        isSingleShot = true;
    }
    else
    {
        isSingleShot = false;
    }
    
    // Setting the saving intervals

    if (exData->GetSaveZPResultsInterval() < exData->GetSimulationTime())  
        SaveZPRIteration = exData->GetSaveZPResultsInterval() / exData->Getdeltat();

    if (exData->GetSaveAllResultsInterval() < exData->GetSimulationTime())
        SaveAllSampleIteration = exData->GetSaveAllResultsInterval() / exData->Getdeltat();

    double laserPeriod = (1.0 / exData->GetLaserPulseRepetition());
    LaserRepetitionIterations = laserPeriod / exData->Getdeltat();
   
    //double LPIExpTable[51];	//Laser Pulse Intensity exp(-0.5 * square) table 
                        //where square = ((currnetStep - 25.0) / 5.0) * ((currnetStep - 25.0) / 5.0);

    double someTempData = 0.0;
    for (int i = 0; i < 51; i++)
    {
        someTempData = ((i - 25.0) / 5.0) * ((i - 25.0) / 5.0);
        LPIExpTable[i] = exp(-0.5 * someTempData);
    }


    // -------------------------  Faster GetHeat Function -------------------------------
    //Exp Array of GetHeat function: exp(-alpha * sampleZnode * spData->GetdZ());
    //where  double alpha = spData->GetSampleAbsorption();
    GetHeat_expArray = new double[spData->GetP()];
    GetHeat_omegaArgument = new double[spData->GetP()];

    betaDivZmeters = (spData->Getbeta_thermal() / spData->GetZmeters());

    for (unsigned int i = 0; i < spData->GetP(); i++)
    {
        someTempData = -spData->GetSampleAbsorption() * i * spData->GetdZ();
        GetHeat_expArray[i] = exp(someTempData)* betaDivZmeters;
        GetHeat_omegaArgument[i] = (0.5 * spData->GetZmeters()) + (i * spData->GetdZ());
    }
    




    //SimulationTime = incomingexData.GetSimulationTime();
    SimulationTime = exData->GetSimulationTime();

    const double sigma = 5.670367e-8; //Stefan-Boltzman const.
    
    double dta = incomingexData.Getdeltat() * incomingsmData.GetAlpha();
    double dtasm = (exData->GetdeltatSlowMotion() * spData->GetAlpha()); //This time is 1/4 of pulse duration - for better modelling of the laser pulse 

    dFox = (dta) / pow(incomingsmData.GetdX(), 2);
    dFoy = (dta) / pow(incomingsmData.GetdY(), 2);
    dFoz = (dta) / pow(incomingsmData.GetdZ(), 2);

    dFoxSlowMotion = (dtasm) / pow(incomingsmData.GetdX(), 2);
    dFoySlowMotion = (dtasm) / pow(incomingsmData.GetdY(), 2);
    dFozSlowMotion = (dtasm) / pow(incomingsmData.GetdZ(), 2);

   
    //Set Number of nodes and T(n) and T(n+1) Temperature arrays;
    N = incomingsmData.GetN();
    O = incomingsmData.GetO();
    P = incomingsmData.GetP();

    halfN = N / 2;
    halfO = O / 2;
    halfP = P / 2;

    //Multithread settings
    if (exData->GetThreads() > 1)
    {
        isMultiThread = true;
    }
    else
    {
        isMultiThread = false;
    }
    ThreadNumber = exData->GetThreads();

    /*
    Maximum threads 

    The maximum number of threads supported by the program is 64. 
    The program can perform calculations for assumed number 
    of threads under conditon: 
    
    N/ThreadsNumber >=2. 
    
    If the condition is not met, the calculations will be performed 
    for a correspondingly smaller number of processor threads.

    */
    cout << endl << endl;
    cout << "<------------------------------------------------------------------------------>" << "\n";
    cout << "<----------------                  CPU threads                ----------------->" << "\n";
    cout << "<------------------------------------------------------------------------------>" << "\n \n";
    cout << endl << " User-requested number of CPU threads: " << ThreadNumber << endl;
    if (ThreadNumber > 64) ThreadNumber = 64;
    //Estimation of maximum ThreadNumber due to N-elements
    //Condition: N / ThreadNumber >=2, if not then ThreadNumber - 1
    while (((!(N / ThreadNumber)) >= 2) && (ThreadNumber >=2))
    {
        ThreadNumber--;
    }
    cout << " Number of processor threads determined by the Z-lambda: " << ThreadNumber << endl;

        ipart = ((N-2) / ThreadNumber);
        ipart1 = N / 4;
        ipart2 = 2 * ipart1;
        ipart3 = 3 * ipart1;


    Temp_n = Create3D<double>(N, O, P);
    Temp_np1 = Create3D<double>(N, O, P);
    
    

    // Setup the temperature arrays by
    for (unsigned int i = 0; i < N; i++) for (unsigned int j = 0; j < O; j++) for (unsigned int k = 0; k < P; k++)
    {
        Temp_n[i][j][k] = spData->GetTemp0();
        Temp_np1[i][j][k] = spData->GetTemp0();
    }

    //Information about RAM allocated by arrays
    double memory = (2 * N * O * P * sizeof(Temp_n) / 1024 / 1024);
    cout << "\n RAM allocated by arrays: "
        << memory << " MB \n \n";


    //External convective heat transfer from the axis 
    // X - dXiR (Right) and dXiL (Left)
    // Y - dYiU (Upper) oraz dYiL (Lower)
    // Z - dZiF (Front) oraz dZiB (Back)

    double lambda = incomingsmData.GetThermalConductivity();

    dXiL = (incomingsmData.Getalfa_xL() * incomingsmData.GetdX()) / lambda;     //double dXiP = (alfa_xP * dx) / lambda; 
    dXiR = (incomingsmData.Getalfa_xR() * incomingsmData.GetdX()) / lambda;     //double dXiL = (alfa_xL * dx) / lambda;

    dYiU = (incomingsmData.Getalfa_yU() * incomingsmData.GetdY()) / lambda;     //double dYiG = (alfa_yG * dy) / lambda; 
    dYiL = (incomingsmData.Getalfa_yL() * incomingsmData.GetdY()) / lambda;     //double dYiD = (alfa_yD * dy) / lambda;

    dZiF = (incomingsmData.Getalfa_zF() * incomingsmData.GetdZ()) / lambda;     //double dZiF = (alfa_zF * dz) / lambda; 
    dZiB = (incomingsmData.Getalfa_zB() * incomingsmData.GetdZ()) / lambda;     // double dZiB = (alfa_zB * dz) / lambda;

    Temp0 = (spData->GetTemp0());
    Tamb_zF = (spData->GetTamb_zF());
    Tamb_zB = (spData->GetTamb_zB());
    Tamb_xR = (spData->GetTamb_xR());
    Tamb_xL = (spData->GetTamb_xL());
    Tamb_yU = (spData->GetTamb_yU());
    Tamb_yD = (spData->GetTamb_yD());
    
    dXrR = (incomingsmData.Getepsilon_xR() * sigma * incomingsmData.GetdX()) / lambda;
    dXrL = (incomingsmData.Getepsilon_xL() * sigma * incomingsmData.GetdX()) / lambda;
    
    dYrU = (incomingsmData.Getepsilon_yU() * sigma * incomingsmData.GetdY()) / lambda;
    dYrL = (incomingsmData.Getepsilon_yL() * sigma * incomingsmData.GetdY()) / lambda;
    
    dZrF = (incomingsmData.Getepsilon_zF() * sigma * incomingsmData.GetdZ()) / lambda;
    dZrB = (incomingsmData.Getepsilon_zB() * sigma * incomingsmData.GetdZ()) / lambda;

    UnitCellVolume = incomingsmData.GetdX() * incomingsmData.GetdY() * incomingsmData.GetdZ();
 
    //Main calculation loop variable
    currentTime = 0;
    deltaSamplePositionSlowIteration = exData->GetSampleVeolocity() * exData->Getdeltat();
    deltaSamplePositionLaserIteration = exData->GetSampleVeolocity() * (exData->GetLaserPulseDuration() / 4); //?????? Sprawd to
    currentSamplePosition = -(exData->GetZScanRange() / 2); //Set sample position at -z

    //Faster calculations variables and pre-calculations
    omegaZeroSquare = exData->GetomegaZero()* exData->GetomegaZero();
    Zrsquare = GetZrSquare();

}

double SimulationData::GetDiameterAtZposition(double Z)
{
        double ZdivZr = (Z / GetZr());
        return sqrt(exData->GetomegaZero() * exData->GetomegaZero() * (1 + (ZdivZr * ZdivZr)));
}

double SimulationData::GetZr()
{
    //Returns Zr = pi* omegazerosqare / lambda (m)
    return ((pi * exData->GetomegaZero() * exData->GetomegaZero()) / (exData->GetLaserWavelength()*1e-9));
}

double SimulationData::GetZrSquare()
{
    return GetZr()*GetZr();
}



bool SimulationData::ShowMeSimulationData(ExperimentData incomingexData, SampleData incomingsmData)
{
    bool AllConditionOK = true;
    double T = incomingexData.GetTimeSteps();
    double t = incomingexData.GetSimulationTime();
    double dt = incomingexData.Getdeltat();
    double a = incomingsmData.GetAlpha();
    double result = 0;
    double T0pow3 = spData->GetTemp0() * spData->GetTemp0() * spData->GetTemp0();

    cout << "\n<------------------------------------------------------------------------------>" << "\n";
    cout << "<----------------   Checking the stability conditions        ------------------>" << "\n";
    cout << "<----------------       of numerical calculations.           ------------------>" << "\n";
    cout << "<------------------------------------------------------------------------------>" << "\n \n";

    
    cout << " Stability for core heat transfer:\n\n";
    cout << "  xi_x = " << dFox;
    if (dFox <= 0.5) cout << " <= 1/2  OK! \n";
    else
    {
        cout << " >= 1/2 Insufficient condition! [xi_x = (a*dt)/dx^2]\n";
        AllConditionOK = false;
    }

    cout << "  xi_x(sm) = " << dFoxSlowMotion;
    if (dFoxSlowMotion <= 0.5) cout << " <= 1/2  OK! \n\n";
    else
    {
        cout << " >= 1/2 Insufficient condition! [xi_x(sm) = (a*dt)/dx^2]\n";
        AllConditionOK = false;
    }

    cout << "  xi_y = " << dFoy;
    if (dFoy <= 0.5) cout << " <= 1/2  OK! \n";
    else
    {
        cout << " >= 1/2 Insufficient condition! [xi_y = (a*dt)/dy^2]\n";
        AllConditionOK = false;
    }

    cout << "  xi_y(sm) = " << dFoySlowMotion;
    if (dFoySlowMotion <= 0.5) cout << " <= 1/2  OK! \n\n";
    else
    {
        cout << " >= 1/2 Insufficient condition! [xi_y(sm) = (a*dt)/dy^2]\n";
        AllConditionOK = false;
    }

    cout << "  xi_z = " << dFoz;
    if (dFoz <= 0.5) cout << " <= 1/2  OK! \n";
    else
    {
        cout << " >= 1/2 Insufficient condition! [xi_z = (a*dt)/dz^2]\n";
        AllConditionOK = false;
    }

    cout << "  xi_z(sm) = " << dFozSlowMotion;
    if (dFozSlowMotion <= 0.5) cout << " <= 1/2  OK! \n \n";
    else
    {
        cout << " >= 1/2 Insufficient condition! [xi_z(sm) = (a*dt)/dz^2]\n";
        AllConditionOK = false;
    }

    
    cout << " Stability for faces heat transfer:\n";
    result = (dFox + dFoy + dFoz + dFox * dXiL + dFox * dXrL * T0pow3);
    cout << "  Left face:  (xi_x - xi_y - xi_z - xi_x * zeta_iL - xi_x * zeta_rL * T^3) = " << result;
    if (result <= 0.5) cout << " <= 1/2  OK! \n";
    else
    {
        cout << " >= 1/2 Insufficient condition! \n";
        AllConditionOK = false;
    }

    result = (dFox + dFoy + dFoz + dFox * dXiL + dFox * dXrL * T0pow3);
    cout << "  Right face: (xi_x - xi_y - xi_z - xi_x * zeta_iR - xi_x * zeta_rR * T^3) = " << result;
    if (result <= 0.5) cout << " <= 1/2  OK! \n\n";
    else
    {
        cout << " >= 1/2 Insufficient condition! \n";
        AllConditionOK = false;
    }
   
    result = (dFox + dFoy + dFoz + dFoy * dYiL + dFoy * dYrL * T0pow3);
    cout << "  Lower face: (xi_x - xi_y - xi_z - xi_y * zeta_iL - xi_y * zeta_rL * T^3) = " << result;
    if (result <= 0.5) cout << " <= 1/2  OK! \n";
    else
    {
        cout << " >= 1/2 Insufficient condition! \n";
        AllConditionOK = false;
    }

    result = (dFox + dFoy + dFoz + dFoy * dYiL + dFoy * dYrL * T0pow3);
    cout << "  Upper face: (xi_x - xi_y - xi_z - xi_y * zeta_iU - xi_y * zeta_rU * T^3) = " << result;
    if (result <= 0.5) cout << " <= 1/2  OK! \n\n";
    else
    {
        cout << " >= 1/2 Insufficient condition! \n";
        AllConditionOK = false;
    }
    
    result = (dFox + dFoy + dFoz + dFoz * dZiF + dFoz * dZrF * T0pow3);
    cout << "  Front face: (xi_x - xi_y - xi_z - xi_z * zeta_iF - xi_z * zeta_rF * T^3) = " << result;
    if (result <= 0.5) cout << " <= 1/2  OK! \n";
    else
    {
        cout << " >= 1/2 Insufficient condition! \n";
        AllConditionOK = false;
    }

    result = (dFox + dFoy + dFoz + dFoz * dZiB + dFoz * dZrB* T0pow3);
    cout << "  Back face:  (xi_x - xi_y - xi_z - xi_z * zeta_iB - xi_z * zeta_rB * T^3) = " << result;
    if (result <= 0.5) cout << " <= 1/2  OK! \n\n";
    else
    {
        cout << " >= 1/2 Insufficient condition! \n";
        AllConditionOK = false;
    }

    cout << " Stability for edges heat transfer:\n";
    result = (dFox + dFoy + dFoz + dFox * dXiL + dFoy * dYiU + dFox * dXrL * T0pow3 + dFoy * dYrU * T0pow3);
    cout << "  Top edges:  satbility result = " << result;
    if (result <= 0.5) cout << " <= 1/2  OK! \n";
    else
    {
        cout << " >= 1/2 Insufficient condition! \n";
        AllConditionOK = false;
    }

    
    cout << " \nTesting the simulation stability. The energy of the sample should be conserved and close to the value of 0.0000 mJ\n\n";


    double SomeTemperature = 300;
    Temp0   = SomeTemperature;
    Tamb_zF = SomeTemperature;
    Tamb_zB = SomeTemperature;
    Tamb_xR = SomeTemperature;
    Tamb_xL = SomeTemperature;
    Tamb_yU = SomeTemperature;
    Tamb_yD = SomeTemperature;
    
    // Setup the temperature arrays by
    for (unsigned int i = 0; i < N; i++) for (unsigned int j = 0; j < O; j++) for (unsigned int k = 0; k < P; k++)
    {
        Temp_n[i][j][k] = SomeTemperature;
        Temp_np1[i][j][k] = SomeTemperature;
    }

    uint64_t ms_start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    for (unsigned short int i = 0; i < 1001; i++)
    {
        CalculateAll();
        if ((i % 50) == 0)
        {
            cout << "  Iteration = " << setw(5) << i << "; Sample Energy = " <<
                (SimulationData::GetSampleCollectedEnergy() * 1000) << " (mJ)\n";
        }
    }
    uint64_t ms_stop = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    
    Time1000iteration = ms_stop - ms_start;

    

    if (GetSampleCollectedEnergy() == -nan("ind"))
    {
        cout << " \nSimulation calculation stability error! Check input data.\n\n";

        cout << " See the 0HeatMap.dat file in Result folder for more details. \n";

        cout << " Please wait, saving 0HeatMap.dat... ";
        
        SaveAllSAmpleHeatMap(currentTimeiteration, 0, HeatMapFileFormat);

        cout << "Done.\n";
            
        return false;
    }
    else
    {
        cout << "\n Calculations stability seems to be correct. Please veryfy the output files and check the sample energies for sure!\n";
        AllConditionOK = true;
    }

    // Setup the temperature arrays by
    for (unsigned int i = 0; i < N; i++) for (unsigned int j = 0; j < O; j++) for (unsigned int k = 0; k < P; k++)
    {
        Temp_n[i][j][k] = spData->GetTemp0();
        Temp_np1[i][j][k] = spData->GetTemp0();
    }


    Temp0 = (spData->GetTemp0());
    Tamb_zF = (spData->GetTamb_zF());
    Tamb_zB = (spData->GetTamb_zB());
    Tamb_xR = (spData->GetTamb_xR());
    Tamb_xL = (spData->GetTamb_xL());
    Tamb_yU = (spData->GetTamb_yU());
    Tamb_yD = (spData->GetTamb_yD());
    
    cout << "\n<------------------------------------------------------------------------------>" << "\n";
    cout << "<----------------          Other Simulation Summary          ------------------>" << "\n";
    cout << "<------------------------------------------------------------------------------>" << "\n \n";
    
    
    cout << " Sample will be irradiated by laser every: " << LaserRepetitionIterations << " iteration(s).\n";
    cout << " All sample heat map will be saved every: " << SaveAllSampleIteration << " iteration(s).\n";
    cout << " Center planes will be saved every: " << SaveZPRIteration << " iteration(s).\n";
    cout << " \n";

    return AllConditionOK;
}



void SimulationData::CalculateCornersTemperature()
{
    
    //Temperature calculation of the corner nodes [0,0,0][N,0,0][N,0,P][0,0,P][0,O,0][N,O,0][N,O,P][0,O,P]
    //double CornerThreeD( double dXi, double dYi, double dZi, double Tcz[3], double Tijk, double Ti_1, double Tj_1, double Tk_1)

    //External convective heat transfer from the axis 
    // X - dXiR (Right) and dXiL (Left)
    // Y - dYiU (Upper) oraz dYiL (Lower)
    // Z - dZiF (Front) oraz dZiB (Back)
    //{ Tcz_xL, Tcz_xD, Tcz_xF };

    double ambientT[3] = {Tamb_xL, Tamb_yD, Tamb_zF };
    Temp_np1[0][0][0] = CornerThreeD(dXiL, dYiL, dZiF, ambientT, Temp_n[0][0][0], Temp_n[1][0][0], Temp_n[0][1][0], Temp_n[0][0][1], dXrL, dYrL, dZrF);

    ambientT[0] = Tamb_xR;
    ambientT[1] = Tamb_yD;
    ambientT[2] = Tamb_zF;
    Temp_np1[N - 1][0][0] = CornerThreeD(dXiR, dYiL, dZiF, ambientT, Temp_n[N - 1][0][0], Temp_n[N - 2][0][0], Temp_n[N - 1][1][0], Temp_n[N - 1][0][1], dXrR, dYrL, dZrF);

    ambientT[0] = Tamb_xR;
    ambientT[1] = Tamb_yD;
    ambientT[2] = Tamb_zB;
    Temp_np1[N - 1][0][P - 1] = CornerThreeD(dXiR, dYiL, dZiB, ambientT, Temp_n[N - 1][0][P - 1], Temp_n[N - 2][0][P - 1], Temp_n[N - 1][1][P - 1], Temp_n[N - 1][0][P - 2], dXrR, dYrL, dZrB);

    ambientT[0] = Tamb_xL;
    ambientT[1] = Tamb_yD;
    ambientT[2] = Tamb_zB;
    Temp_np1[0][0][P - 1] = CornerThreeD(dXiL, dYiL, dZiB, ambientT, Temp_n[0][0][P - 1], Temp_n[1][0][P - 1], Temp_n[0][1][P - 1], Temp_n[0][0][P - 2], dXrL, dYrL, dZrB);

    ambientT[0] = Tamb_xL;
    ambientT[1] = Tamb_yU;
    ambientT[2] = Tamb_zF;
    Temp_np1[0][O - 1][0] = CornerThreeD(dXiL, dYiU, dZiF, ambientT, Temp_n[0][O - 1][0], Temp_n[1][O - 1][0], Temp_n[0][O - 2][0], Temp_n[0][O - 1][1], dXrL, dYrU, dZrF);

    ambientT[0] = Tamb_xR;
    ambientT[1] = Tamb_yU;
    ambientT[2] = Tamb_zF;
    Temp_np1[N - 1][O - 1][0] = CornerThreeD(dXiR, dYiU, dZiF, ambientT, Temp_n[N - 1][O - 1][0], Temp_n[N - 2][O - 1][0], Temp_n[N - 1][O - 2][0], Temp_n[N - 1][O - 1][1], dXrR, dYrU, dZrF);

    ambientT[0] = Tamb_xR;
    ambientT[1] = Tamb_yU;
    ambientT[2] = Tamb_zB;
    Temp_np1[N - 1][O - 1][P - 1] = CornerThreeD(dXiR, dYiU, dZiB, ambientT, Temp_n[N - 1][O - 1][P - 1], Temp_n[N - 2][O - 1][P - 1], Temp_n[N - 1][O - 2][P - 1], Temp_n[N - 1][O - 1][P - 2], dXrR, dYrU, dZrB);

    ambientT[0] = Tamb_xL;
    ambientT[1] = Tamb_yU;
    ambientT[2] = Tamb_zB;
    Temp_np1[0][O - 1][P - 1] = CornerThreeD(dXiL, dYiU, dZiB, ambientT, Temp_n[0][O - 1][P - 1], Temp_n[1][O - 1][P - 1], Temp_n[0][O - 2][P - 1], Temp_n[0][O - 1][P - 2], dXrL, dYrU, dZrB);


}

void SimulationData::CalculateWallsTemperature()
{
    //Six Walls (Front, Back, Left, Right, Upper, Lower)
    if (isMultiThread) 
    {
        //Front and Back Walls
        thread thread1{ [&]() { 
            for (unsigned int i = 1; i <= (N - 2); i++)
            {
                for (unsigned int j = 1; j <= (O - 2); j++)
                {
                    
                    // double SurfaceThreeD(double dFo_a, double dFo_b, double dFo_c, 
                    //double Tijk, 
                    //double Ta_1, double Tap1, 
                    //double Tb_1, double Tbp1, double Tcp1, 
                    //double dCi, double dCr, double Tcz);

                    //Front Wall
                    Temp_np1[i][j][0] =
                    SurfaceThreeD(dFox, dFoy, dFoz,
                                  Temp_n[i][j][0],
                                  Temp_n[i - 1][j][0], Temp_n[i + 1][j][0],
                                  Temp_n[i][j - 1][0], Temp_n[i][j + 1][0], Temp_n[i][j][1],
                                  dZiF, dZrF, Tamb_zF);

                    //Back Wall
                    Temp_np1[i][j][P - 1] =
                    SurfaceThreeD(dFox, dFoy, dFoz,
                                  Temp_n[i][j][P - 1], Temp_n[i - 1][j][P - 1], Temp_n[i + 1][j][P - 1],
                                  Temp_n[i][j - 1][P - 1], Temp_n[i][j + 1][P - 1], Temp_n[i][j][P - 2],
                                  dZiB, dZrB, Tamb_zB);

                }
            }
        } };

        thread thread2{ [&]()
        {

            for (unsigned int j = 1; j <= (O - 2); j++)
            {
                for (unsigned int k = 1; k <= (P - 2); k++)
                {

                    //Left Wall
                    Temp_np1[0][j][k] =
                        SurfaceThreeD(dFoy, dFoz, dFox,
                            Temp_n[0][j][k],
                            Temp_n[0][j - 1][k], Temp_n[0][j + 1][k],
                            Temp_n[0][j][k - 1], Temp_n[0][j][k + 1], Temp_n[1][j][k],
                            dXiL, dXrL, Tamb_xL);

                    //Right Wall
                    Temp_np1[N - 1][j][k] =
                        SurfaceThreeD(dFoy, dFoz, dFox,
                            Temp_n[N - 1][j][k],
                            Temp_n[N - 1][j - 1][k], Temp_n[N - 1][j + 1][k],
                            Temp_n[N - 1][j][k - 1], Temp_n[N - 1][j][k + 1], Temp_n[N - 2][j][k],
                            dXiR, dXrR, Tamb_xR);
                   
                }
            }
        } };

        thread thread3{ [&]()
        {
            for (unsigned int i = 1; i <= (N - 2); i++)
            {
                for (unsigned int k = 1; k <= (P - 2); k++)
                {

                    //Lower Wall
                    Temp_np1[i][0][k] =
                        SurfaceThreeD(dFox, dFoz, dFoy,
                            Temp_n[i][0][k],
                            Temp_n[i - 1][0][k], Temp_n[i + 1][0][k],
                            Temp_n[i][0][k - 1], Temp_n[i][0][k + 1], Temp_n[i][1][k],
                            dYiL, dYrL, Tamb_yD);

                    //Upper Wall
                    Temp_np1[i][O - 1][k] =
                        SurfaceThreeD(dFox, dFoz, dFoy,
                            Temp_n[i][O - 1][k],
                            Temp_n[i - 1][O - 1][k], Temp_n[i + 1][O - 1][k],
                            Temp_n[i][O - 1][k - 1], Temp_n[i][O - 1][k + 1], Temp_n[i][O - 2][k],
                            dYiU, dYrU, Tamb_yU);
                }
            }
        } };

        thread1.join();
        thread2.join();
        thread3.join();
    }
    else
    {
        //Front and Back Walls
        for (unsigned int i = 1; i <= (N - 2); i++)
        {
            for (unsigned int j = 1; j <= (O - 2); j++)
            {

                //Front Wall
                Temp_np1[i][j][0] =
                    SurfaceThreeD(dFox, dFoy, dFoz,
                        Temp_n[i][j][0],
                        Temp_n[i - 1][j][0], Temp_n[i + 1][j][0],
                        Temp_n[i][j - 1][0], Temp_n[i][j + 1][0], Temp_n[i][j][1],
                        dZiF, dZrF, Tamb_zF);

                //Back Wall
                Temp_np1[i][j][P - 1] =
                    SurfaceThreeD(dFox, dFoy, dFoz,
                        Temp_n[i][j][P - 1], Temp_n[i - 1][j][P - 1], Temp_n[i + 1][j][P - 1],
                        Temp_n[i][j - 1][P - 1], Temp_n[i][j + 1][P - 1], Temp_n[i][j][P - 2],
                        dZiB, dZrB, Tamb_zB);

            }
        }


        //Left and Rights Walls
        for (unsigned int j = 1; j <= (O - 2); j++)
        {
            for (unsigned int k = 1; k <= (P - 2); k++)
            {

                //Left Wall
                Temp_np1[0][j][k] =
                    SurfaceThreeD(dFoy, dFoz, dFox,
                        Temp_n[0][j][k],
                        Temp_n[0][j - 1][k], Temp_n[0][j + 1][k],
                        Temp_n[0][j][k - 1], Temp_n[0][j][k + 1], Temp_n[1][j][k],
                        dXiL, dXrL, Tamb_xL);

                //Right Wall
                Temp_np1[N - 1][j][k] =
                    SurfaceThreeD(dFoy, dFoz, dFox,
                        Temp_n[N - 1][j][k],
                        Temp_n[N - 1][j - 1][k], Temp_n[N - 1][j + 1][k],
                        Temp_n[N - 1][j][k - 1], Temp_n[N - 1][j][k + 1], Temp_n[N - 2][j][k],
                        dXiR, dXrR, Tamb_xR);
            }
        }

        //Upper and Lower Walls
        for (unsigned int i = 1; i <= (N - 2); i++)
        {
            for (unsigned int k = 1; k <= (P - 2); k++)
            {

                //Lower Wall
                Temp_np1[i][0][k] =
                    SurfaceThreeD(dFox, dFoz, dFoy,
                        Temp_n[i][0][k],
                        Temp_n[i - 1][0][k], Temp_n[i + 1][0][k],
                        Temp_n[i][0][k - 1], Temp_n[i][0][k + 1], Temp_n[i][1][k],
                        dYiL, dYrL, Tamb_yD);

                //Upper Wall
                Temp_np1[i][O - 1][k] =
                    SurfaceThreeD(dFox, dFoz, dFoy,
                        Temp_n[i][O - 1][k],
                        Temp_n[i - 1][O - 1][k], Temp_n[i + 1][O - 1][k],
                        Temp_n[i][O - 1][k - 1], Temp_n[i][O - 1][k + 1], Temp_n[i][O - 2][k],
                        dYiU, dYrU, Tamb_yU);
            }
        }
    }
}

void SimulationData::SetNormalTimeSteps()
{
    double dta = exData->Getdeltat() * spData->GetAlpha();

    dFox = (dta) / pow(spData->GetdX(), 2);
       
    dFoy = (dta) / pow(spData->GetdY(), 2);

    dFoz = (dta) / pow(spData->GetdZ(), 2);

}

void SimulationData::SetShortTimeSteps()
{
    double dtasm = (exData->GetdeltatSlowMotion()  * spData->GetAlpha()); //This time is 1/4 of pulse duration - for better modelling of the laser pulse 

    dFox = (dtasm) / pow(spData->GetdX(), 2);

    dFoy = (dtasm) / pow(spData->GetdY(), 2);

    dFoz = (dtasm) / pow(spData->GetdZ(), 2);
}

double SimulationData::ShootToThrill(double Z)
{
    //This function returns the time elapsed during the laser pulse
    // and simulate the irradiation of the sample at Z positions.
    double ElapsedTime = 0.0;
    double CurrentPower = 0.0;
    double PulseEnergy = exData->GetLaserPulseEnergy();

    double denominator = (spData->GetSampleSpecificHeatCapacity() * spData->GetSampleDensity());
    
    SetShortTimeSteps();

    double xdim = -spData->GetXmeters() / 2;
    double ydim = -spData->GetYmeters() / 2;
    double zdim = -spData->GetZmeters() / 2;

    if (isMultiThread)
    {

        auto function = [&](unsigned short int starti, unsigned short int endi)
        {
            double dtemperature = 0;
            for (unsigned int i = starti; i <= endi; i++)
                for (unsigned int j = 0; j < O; j++)
                    for (unsigned int k = 0; k < P; k++)
                    {
                        //Heat per unit vol. * Unit Vol. * time = Energy (W * m^-3 * s = J * m^-3)
                        //Below is Ved -> Volumetric energy density (J/m^3)
                        dtemperature = (GetHeat(xdim + i * spData->GetdX(),
                            ydim + j * spData->GetdY(),
                            Z,
                            CurrentPower, k) * exData->GetdeltatSlowMotion());//(J/m^3)

                        //deltaTemperature = Ved / (specific heat capacity * density), result in (K)
                        //Below is delta Temperature (K)
                        dtemperature /= denominator;
                        Temp_n[i][j][k] += dtemperature;
                    }
        };

        for (int dt = 0; dt < 51; dt++)
        {
        
            thread myThreads[64];

            istart = 0;
            iend = ipart;
            for (int i = 0; i < ThreadNumber; i++) {

                CurrentPower = LaserPulseIntensity(dt, PulseEnergy);
                
                myThreads[i] = thread(function, istart, iend);
                
                istart  += ipart;
                iend    += ipart;
                
                if ((i + 2) == ThreadNumber) iend = N - 2;
                // cout<<"istart = " << istart << "iend = " << iend << " ipart = " << ipart<< endl;
            }

            for (int i = 0; i < ThreadNumber; i++) myThreads[i].join();
        }       
    }
    else
    {
        for (int dt = 0; dt < 51; dt++)
        {

            CurrentPower = LaserPulseIntensity(dt, PulseEnergy);
            double dtemperature = 0;

            for (unsigned int i = 0; i < N; i++)
                for (unsigned int j = 0; j < O; j++)
                    for (unsigned int k = 0; k < P; k++)
                    {
                        //Heat per unit vol. * Unit Vol. * time = Energy (W * m^-3 * s = J * m^-3)
                        //Below is Ved -> Volumetric energy density (J/m^3)
                        dtemperature = (GetHeat(xdim + i * spData->GetdX(),
                            ydim + j * spData->GetdY(),
                            Z,
                            CurrentPower, k) * exData->GetdeltatSlowMotion());//(J/m^3)

                        //deltaTemperature = Ved / (specific heat capacity * density), result in (K)
                        //Below is delta Temperature (K)
                        dtemperature /= denominator;
                        Temp_n[i][j][k] += dtemperature;
                    }

            //CalculateAll();
            ElapsedTime += exData->GetdeltatSlowMotion();
        }
    }
    
    SetNormalTimeSteps();
    return ElapsedTime;
}

double SimulationData::GetHeat(double x, double y, double z, double LaserPower, unsigned int sampleZnode)
{
    //Function returns the heat per unit at x,y,z location of the sample i the depth of sampleZnode * dz
    //The result is in W/m^3

    double omega = GetOmega(z - GetHeat_omegaArgument[sampleZnode]);

    double a = (2 * LaserPower) / (pi * omega * omega);
    double square;

    square = (x * x + y * y) / (omega * omega);

    return a * exp( - 2 * square) * GetHeat_expArray[sampleZnode];
}

double SimulationData::GetSampleCollectedEnergy()
{
    double a = spData->GetSampleDensity() * spData->GetSampleSpecificHeatCapacity() 
                * spData->GetdX() * spData->GetdY() * spData->GetdY();
    double CollectedEnergy = 0.0;
    
    for (unsigned int i = 0; i<N; i++)
        for (unsigned int j = 0; j < O; j++)
            for (unsigned int k = 0; k < P; k++)
            {
                CollectedEnergy += ((Temp_n[i][j][k] - Temp0) * a);
            }
    
    
    return CollectedEnergy;
}

double SimulationData::LaserPulseIntensity(int currnetStep, double Amplitude)
{
    //Amplitude is the Pulse Energy in (J)
    //Full Gaussian pulse distribution is returned betwen interations 0 and 50.
    //The function returns the Power at the currentStep (J/s = W);
    //Array LPIExpTAble contains the exp(-0.5 * someconst) values
   
    return ((Amplitude / exData->GetLaserPulseDuration()) / pi) * LPIExpTable[currnetStep];
}

void SimulationData::CalculateAll()
{
    CalculateCornersTemperature();

    CalculateEdgeTemperatures();

    CalculateCoreTemperature();

    CalculateWallsTemperature();

    CopyArray();
}

void SimulationData::CalculateEdgeTemperatures()
{

    for (unsigned int i = 1; i <= (N - 2); i++)
        // 4 krawedzie w osi X (przednia gorna i dolna, tylna gorna i dolna)
        //Four edges of X axis (front upper and lower, back upper and lower)
    {
        //Front upper
        double ambientT[2] = {Tamb_yU ,Tamb_zF};
       
        //EdgeThreeD(double dFo_a, double dFo_b, double dFo_c, double dBi, double dCi, double Tcz[2], double Tijk, double Ta_1, double Tap1, double Tb_1, double Tc_1);
        Temp_np1[i][O - 1][0] = EdgeThreeD(dFox, dFoy, dFoz, dYiU, dZiF, ambientT,
            Temp_n[i][O - 1][0],
            Temp_n[i - 1][O - 1][0],
            Temp_n[i + 1][O - 1][0],
            Temp_n[i][O - 2][0],
            Temp_n[i][O - 1][1], dYrU, dZrF
        );

        //dFoy, dFoz, dYiG, dZiF, ambientT, Temp_n[i][0][0], Temp_n[i][O-2][0], Temp_n[i][1][1]);

        //Front lower
        ambientT[0] = Tamb_yD;
        ambientT[1] = Tamb_zF;
        Temp_np1[i][0][0] = EdgeThreeD(dFox, dFoy, dFoz, dYiL, dZiF, ambientT,
            Temp_n[i][0][0],
            Temp_n[i - 1][0][0],
            Temp_n[i + 1][0][0],
            Temp_n[i][1][0],
            Temp_n[i][0][1], dYrL, dZrF
        );


        //back upper
        ambientT[0] = Tamb_yU;
        ambientT[1] = Tamb_zB;
        //EdgeThreeD(double dFo_a, double dFo_b, double dFo_c, double dBi, double dCi, double Tcz[2], double Tijk, double Ta_1, double Tap1, double Tb_1, double Tc_1);
        Temp_np1[i][O - 1][P - 1] = EdgeThreeD(dFox, dFoy, dFoz, dYiU, dZiB, ambientT,
            Temp_n[i][O - 1][P - 1],
            Temp_n[i - 1][O - 1][P - 1],
            Temp_n[i + 1][O - 1][P - 1],
            Temp_n[i][O - 2][P - 1],
            Temp_n[i][O - 1][P - 2], dYrL, dZrB
        );


        //dFoy, dFoz, dYiG, dZiF, ambientT, Temp_n[i][0][0], Temp_n[i][O-2][0], Temp_n[i][1][1]);

        //Back lower
        ambientT[0] = Tamb_yD;
        ambientT[1] = Tamb_zB;
        Temp_np1[i][0][P - 1] = EdgeThreeD(dFox, dFoy, dFoz, dYiL, dZiB, ambientT,
            Temp_n[i][0][P - 1],
            Temp_n[i - 1][0][P - 1],
            Temp_n[i + 1][0][P - 1],
            Temp_n[i][1][P - 1],
            Temp_n[i][0][P - 2], dYrL, dZrB
        );

    }

    for (unsigned int j = 1; j <= (O - 2); j++)
        // 4 krawedzie w osi Y (lewe przednia i tylna, prawa przednia i tylna)
        //Four edges of Y axis (left front and back, right front and back)
    {
        double ambientT[2];
        //Left front
        ambientT[0] = Tamb_xL;
        ambientT[1] = Tamb_zF;
        //EdgeThreeD(double dFo_a, double dFo_b, double dFo_c, double dBi, double dCi, double Tcz[2], double Tijk, double Ta_1, double Tap1, double Tb_1, double Tc_1);
        Temp_np1[0][j][0] = EdgeThreeD(dFox, dFoy, dFoz, dXiL, dZiF, ambientT,
            Temp_n[0][j][0],
            Temp_n[0][j - 1][0],
            Temp_n[0][j + 1][0],
            Temp_n[1][j][0],
            Temp_n[0][j][1], dXrL, dZrF);


        //Left back
        ambientT[0] = Tamb_xL;
        ambientT[1] = Tamb_zB;
        //EdgeThreeD(double dFo_a, double dFo_b, double dFo_c, double dBi, double dCi, double Tcz[2], double Tijk, double Ta_1, double Tap1, double Tb_1, double Tc_1);
        Temp_np1[0][j][P - 1] = EdgeThreeD(dFox, dFoy, dFoz, dXiL, dZiB, ambientT,
            Temp_n[0][j][P - 1],
            Temp_n[0][j - 1][P - 1],
            Temp_n[0][j + 1][P - 1],
            Temp_n[1][j][P - 1],
            Temp_n[0][j][P - 2], dXrL, dZrB);

        //Right front
        ambientT[0] = Tamb_xR;
        ambientT[1] = Tamb_zF;
        //EdgeThreeD(double dFo_a, double dFo_b, double dFo_c, double dBi, double dCi, double Tcz[2], double Tijk, double Ta_1, double Tap1, double Tb_1, double Tc_1);
        Temp_np1[N - 1][j][0] = EdgeThreeD(dFox, dFoy, dFoz, dXiR, dZiF, ambientT,
            Temp_n[N - 1][j][0],
            Temp_n[N - 1][j - 1][0],
            Temp_n[N - 1][j + 1][0],
            Temp_n[N - 2][j][0],
            Temp_n[N - 1][j][1], dXrR, dZrF);

        //Right back
        ambientT[0] = Tamb_xR;
        ambientT[1] = Tamb_zB;
        //EdgeThreeD(double dFo_a, double dFo_b, double dFo_c, double dBi, double dCi, double Tcz[2], double Tijk, double Ta_1, double Tap1, double Tb_1, double Tc_1);
        Temp_np1[N - 1][j][P - 1] = EdgeThreeD(dFox, dFoy, dFoz, dXiR, dZiB, ambientT,
            Temp_n[N - 1][j][P - 1],
            Temp_n[N - 1][j - 1][P - 1],
            Temp_n[N - 1][j + 1][P - 1],
            Temp_n[N - 2][j][P - 1],
            Temp_n[N - 1][j][P - 2], dXrR, dZrB);
    }

    for (unsigned int k = 1; k <= (P - 2); k++)
        // 4 krawedzie w osi Y (lewe przednia i tylna, prawa przednia i tylna)
        //Four edges of Y axis (left front and back, right front and back)
    {
        double ambientT[2];
        //Left lower
        ambientT[0] = Tamb_yD;
        ambientT[1] = Tamb_xL;
        //EdgeThreeD(double dFo_a, double dFo_b, double dFo_c, double dBi, double dCi, double Tcz[2], double Tijk, double Ta_1, double Tap1, double Tb_1, double Tc_1);
        Temp_np1[0][0][k] = EdgeThreeD(dFox, dFoy, dFoz, dXiL, dZiF, ambientT,
            Temp_n[0][0][k],
            Temp_n[0][0][k - 1],
            Temp_n[0][0][k + 1],
            Temp_n[1][0][k],
            Temp_n[0][1][k], dYrL, dXrL);

        //Left upper
        ambientT[0] = Tamb_yU;
        ambientT[1] = Tamb_xL;
        //EdgeThreeD(double dFo_a, double dFo_b, double dFo_c, double dBi, double dCi, double Tcz[2], double Tijk, double Ta_1, double Tap1, double Tb_1, double Tc_1);
        Temp_np1[0][O - 1][k] = EdgeThreeD(dFox, dFoy, dFoz, dXiL, dZiB, ambientT,
            Temp_n[0][O - 1][k],
            Temp_n[0][O - 1][k - 1],
            Temp_n[0][O - 1][k + 1],
            Temp_n[1][O - 1][k],
            Temp_n[0][O - 2][k], dYrU, dXrL);

        //Right lower
        ambientT[0] = Tamb_yD;
        ambientT[1] = Tamb_xR;
        //EdgeThreeD(double dFo_a, double dFo_b, double dFo_c, double dBi, double dCi, double Tcz[2], double Tijk, double Ta_1, double Tap1, double Tb_1, double Tc_1);
        Temp_np1[N - 1][0][k] = EdgeThreeD(dFox, dFoy, dFoz, dXiR, dZiF, ambientT,
            Temp_n[N - 1][0][k],
            Temp_n[N - 1][0][k - 1],
            Temp_n[N - 1][0][k + 1],
            Temp_n[N - 2][0][k],
            Temp_n[N - 1][1][k], dYrL, dXrR);


        //Right upper
        ambientT[0] = Tamb_yU;
        ambientT[1] = Tamb_xR;
        //EdgeThreeD(double dFo_a, double dFo_b, double dFo_c, double dBi, double dCi, double Tcz[2], double Tijk, double Ta_1, double Tap1, double Tb_1, double Tc_1);
        Temp_np1[N - 1][O - 1][k] = EdgeThreeD(dFox, dFoy, dFoz, dXiR, dZiB, ambientT,
            Temp_n[N - 1][O - 1][k],
            Temp_n[N - 1][O - 1][k - 1],
            Temp_n[N - 1][O - 1][k + 1],
            Temp_n[N - 2][O - 1][k],
            Temp_n[N - 1][O - 2][k], dYrU, dXrR);
    }


    
}



void SimulationData::CalculateCoreTemperatureThreads(unsigned short int starti, unsigned short int endi)
{
    for (unsigned int k = 1; k <= (P - 2); k++) for (unsigned int i = starti; i <= endi; i++) for (unsigned int j = 1; j <= (O - 2); j++)
    {
        double Temptwotimes = 2 * Temp_n[i][j][k];
        Temp_np1[i][j][k] = (Temp_n[i - 1][j][k] + Temp_n[i + 1][j][k] - Temptwotimes) * dFox;
        Temp_np1[i][j][k] += (Temp_n[i][j - 1][k] + Temp_n[i][j + 1][k] - Temptwotimes) * dFoy;
        Temp_np1[i][j][k] += (Temp_n[i][j][k - 1] + Temp_n[i][j][k + 1] - Temptwotimes) * dFoz;
        Temp_np1[i][j][k] += Temp_n[i][j][k];
    }
}

void SimulationData::CalculateCoreTemperature()
{

//Temperature calculation inside the sample Nodes: [1 .. N-2][1 .. O-2][1 .. N-2]
    auto function = [&](unsigned short int starti, unsigned short int endi) {
        for (unsigned int k = 1; k <= (P - 2); k++) for (unsigned int i = starti; i <= endi; i++) for (unsigned int j = 1; j <= (O - 2); j++)
        {
            double Temptwotimes = 2 * Temp_n[i][j][k];
            Temp_np1[i][j][k] = (Temp_n[i - 1][j][k] + Temp_n[i + 1][j][k] - Temptwotimes) * dFox;
            Temp_np1[i][j][k] += (Temp_n[i][j - 1][k] + Temp_n[i][j + 1][k] - Temptwotimes) * dFoy;
            Temp_np1[i][j][k] += (Temp_n[i][j][k - 1] + Temp_n[i][j][k + 1] - Temptwotimes) * dFoz;
            Temp_np1[i][j][k] += Temp_n[i][j][k];
        }
    };


    if (isMultiThread)
    {
        
        thread myThreads[64];

        istart = 1;
        iend = ipart;
        for (int i = 0; i < ThreadNumber; i++) {
            myThreads[i] = thread(function, istart, iend);
            istart  += ipart;
            iend    += ipart;
            if ((i + 2) == ThreadNumber) iend = N - 2;
            
        }

        for (int i = 0; i < ThreadNumber; i++) myThreads[i].join();


    }
    else 
    {
        for (unsigned int k = 1; k <= (P - 2); k++) for (unsigned int i = 1; i <= (N - 2); i++) for (unsigned int j = 1; j <= (O - 2); j++)
        {
            double Temptwotimes = 2 * Temp_n[i][j][k];
            Temp_np1[i][j][k] = (Temp_n[i - 1][j][k] + Temp_n[i + 1][j][k] - Temptwotimes) * dFox;
            Temp_np1[i][j][k] += (Temp_n[i][j - 1][k] + Temp_n[i][j + 1][k] - Temptwotimes) * dFoy;
            Temp_np1[i][j][k] += (Temp_n[i][j][k - 1] + Temp_n[i][j][k + 1] - Temptwotimes) * dFoz;
            Temp_np1[i][j][k] += Temp_n[i][j][k];
        }

    }
    

}

void SimulationData::CopyArray()
{
    if ((isMultiThread) && (ThreadNumber >= 4))
    {
        //Copy the temperature nodes - preparation for the next iteration
        thread th1{ [&]()
        {
            for (unsigned int i = 0; i < ipart1; i++) for (unsigned int j = 0; j < O; j++) for (unsigned int k = 0; k < P; k++)
            {
                Temp_n[i][j][k] = Temp_np1[i][j][k];
            }
        } };

        thread th2{ [&]()
        {
            for (unsigned int i = ipart1; i < ipart2; i++) for (unsigned int j = 0; j < O; j++) for (unsigned int k = 0; k < P; k++)
            {
                Temp_n[i][j][k] = Temp_np1[i][j][k];
            }
        } };

        thread th3{ [&]()
        {
            for (unsigned int i = ipart2; i < ipart3; i++) for (unsigned int j = 0; j < O; j++) for (unsigned int k = 0; k < P; k++)
            {
                Temp_n[i][j][k] = Temp_np1[i][j][k];
            }
        } };

        thread th4{ [&]()
        {
            for (unsigned int i = ipart3; i < N; i++) for (unsigned int j = 0; j < O; j++) for (unsigned int k = 0; k < P; k++)
            {
                Temp_n[i][j][k] = Temp_np1[i][j][k];
            }
        } };

        th1.join();
        th2.join();
        th3.join();
        th4.join();
    }
    else
    {
        for (unsigned int i = 0; i < N; i++) for (unsigned int j = 0; j < O; j++) for (unsigned int k = 0; k < P; k++)
        {
            Temp_n[i][j][k] = Temp_np1[i][j][k];
        }
    }
    
    
}

void SimulationData::ShowCalculationProgress(unsigned int currentIteration, unsigned int TotalIterations)
{
    if (currentIteration == 0) cout << endl << char(218) << " 0%                     25%                     50%                     75%                    100% " << char(191) << endl<<" ";
    if (currentIteration % ((TotalIterations) / 100) == 0) //Progress bar
    {
        cout << "X"<<flush;
    }

}

SimulationData::SimulationData()
{

}

SimulationData::~SimulationData()
{
}

void SimulationData::StartCalculation(unsigned int NumberOfTimeIterations)
{
    unsigned int SaveAllCenterProfilesIndex = 0;
    unsigned int SaveAllCenterHeatMapSurfacesIndex = 0;
    unsigned int SaveAllSAmpleHeatMapIndex = 0;
 
    unsigned int SavePlanes = SaveZPRIteration;                 // Iteration "0" save the Planes
    unsigned int SaveHeatMap = SaveAllSampleIteration;          // Iteration "0" save Heat Map
    unsigned int ShotLaser = (LaserRepetitionIterations - 1);   // the first iteration is without pulse, second is with the laser pulse
    unsigned int SampleCollectedEnergyIterations = LaserRepetitionIterations / 4;
    unsigned int SaveSampleCollectedEnergy = SampleCollectedEnergyIterations;
    unsigned int NumberOfLaserShots = 1;
    bool stopshoot = false;
    float CalculationProgress = 0;

    ofstream file;


    if (isSingleShot)
    {
        
        double approxsimulationtime = 0;
        unsigned short int TestIteration = 1000;
        double ZScanRangeTime = (exData->GetZScanRange() / exData->GetSampleVeolocity());

        unsigned short int TotalNumberOfLaserPulses = ZScanRangeTime * exData->GetLaserPulseRepetition();
        unsigned short int TotalNumberOfSavingProfiles = NumberOfTimeIterations / SaveZPRIteration;
        unsigned short int TotalNumberOfHeatMapSavings = NumberOfTimeIterations / SaveAllSampleIteration;

        uint64_t LaserPulseGenerationTime = 0;
        uint64_t SavingDataProfilesAndMapTime = 0;
        uint64_t SavingHeatMapTime = 0;
        uint64_t Iteration1000Time = 0;
        uint64_t ms_start = 0;
        uint64_t ms_end = 0;

        cout << " in Single Shot Mode.\n\n";
        
        cout << " Irradiating the sample placed at z = 0 by laser: ";
        

        // ----------  Laser one shot
        ms_start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        currentTime += ShootToThrill(0);
        ms_end = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        LaserPulseGenerationTime = ms_end - ms_start;
        cout << " Done. \n Retrieving  the time of 1000 iterations: ";


        // ----------  1000 Iteration
        ms_start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        do
        {
            CalculateAll();
            ShowCalculationProgress(currentTimeiteration, TestIteration);
            currentTimeiteration++;   
        } while (currentTimeiteration < TestIteration);
        ms_end = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        
        Iteration1000Time =  ms_end - ms_start;

        cout << " Done. \n Saving center profiles and center heat maps: ";

        
        // ----------  Porfiles and Surfaces
        ms_start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        //SavePlanes function launch
        SaveAllCenterProfiles(SaveAllCenterProfilesIndex);
        SaveAllCenterHeatMapSurfaces(SaveAllCenterHeatMapSurfacesIndex);
        ms_end = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        SavingDataProfilesAndMapTime = ms_end - ms_start;
        cout << " Done. \n Saving full heat map of the sample. ";

        
        // ----------   Full Heat Maps
        ms_start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        SaveAllSAmpleHeatMap(currentTimeiteration, SaveAllSAmpleHeatMapIndex, HeatMapFileFormat);
        SaveAllSAmpleHeatMapIndex++;
        ms_end = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        SavingHeatMapTime = ms_end - ms_start;
        cout << " Done. \n";
     
        cout << endl<< endl;
        cout << " ---------------------------------------------------------------------------------------\n";
        cout << " -                                      -   Single  -              -                   -\n";
        cout << " -  Procedure                           -    Time   -   Quantity   -    Total Time     -\n";
        cout << " -                                      -    [ms]   -              -       [ms]        -\n";
        cout << " ---------------------------------------------------------------------------------------\n";
        cout << " - Single Laser Pulse Generation        - "
            << fixed << setw(9) << setprecision(1) << LaserPulseGenerationTime << " - "
            << fixed << setw(12) << setprecision(1) << TotalNumberOfLaserPulses << " - "
            << fixed << setw(17) << setprecision(1) << LaserPulseGenerationTime * TotalNumberOfLaserPulses << " -\n";

        cout << " - 1000 Iteration of Heat Transfer      - " 
            << fixed << setw(9) << setprecision(1) << Iteration1000Time << " - "
            << fixed << setw(12) << setprecision(1) << NumberOfTimeIterations << " - "
            << fixed << setw(17) << setprecision(1) << ((Iteration1000Time/1000) * NumberOfTimeIterations) << " -\n";
        
        cout << " - Full Heat Map Saving                 - " 
            << fixed << setw(9) << setprecision(1) << SavingHeatMapTime << " - "
            << fixed << setw(12) << setprecision(1) << TotalNumberOfHeatMapSavings << " - "
            << fixed << setw(17) << setprecision(1) << SavingHeatMapTime * TotalNumberOfHeatMapSavings << " -\n";

        cout << " - Profiles and Center Surfaces Saving  - "
            << fixed << setw(9) << setprecision(1) << SavingDataProfilesAndMapTime << " - "
            << fixed << setw(12) << setprecision(1) << TotalNumberOfSavingProfiles << " - "
            << fixed << setw(17) << setprecision(1) << SavingDataProfilesAndMapTime * TotalNumberOfSavingProfiles << " -\n";
        cout << " ---------------------------------------------------------------------------------------\n";

        approxsimulationtime = (Iteration1000Time/1000) * NumberOfTimeIterations;
        approxsimulationtime += LaserPulseGenerationTime * TotalNumberOfLaserPulses;
        approxsimulationtime += Iteration1000Time * NumberOfTimeIterations;
        approxsimulationtime += SavingHeatMapTime * TotalNumberOfHeatMapSavings;
        approxsimulationtime += SavingDataProfilesAndMapTime * TotalNumberOfSavingProfiles;
        approxsimulationtime /= 1000;

        cout << " -                                         Approximate total time: -"
            << fixed << setw(18) << setprecision(0) << approxsimulationtime << " -\n";

        approxsimulationtime /= 1000;//results in seconds
        approxsimulationtime /= 60; //results in minutes
        approxsimulationtime /= 60; //results in hours
        cout << " ---------------------------------------------------------------------------------------\n";
        cout << " - Estimated total calculation time     : "
            << fixed << setw(10) << setprecision(1) << approxsimulationtime << 
            " hours                             -\n";
        cout << " -                              which is: "
            << fixed << setw(10) << setprecision(2) << (approxsimulationtime /= 24) << 
            " days                              -\n";
        cout << " -                                                                                     -\n";

        if (isMultiThread) {
            approxsimulationtime = ThreadNumber;
        }
        else
        {
            approxsimulationtime = 1;
        } 
        cout << " - Number of Processor Threads used in calculations: "
            << fixed << setw(3) << setprecision(0) << approxsimulationtime << 
            "                               -\n";
        cout << " -                                                                                     -\n";
        cout << " - Energy collected by sample: " 
            << scientific << setw(15) << setprecision(5) << (SimulationData::GetSampleCollectedEnergy()*1000)
            << " (mJ)                                    -\n";
        cout << " ---------------------------------------------------------------------------------------\n";

        ms_start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        approxsimulationtime = GetSampleCollectedEnergy();
        ms_end = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();

        

    }
    else
    {
        do
        {
            CalculateAll();
            
            if ((ShotLaser == LaserRepetitionIterations) && !stopshoot)
            {
                currentTime += ShootToThrill(currentSamplePosition);
                CalculationProgress = (currentTimeiteration * 100.0) / NumberOfTimeIterations;
  
                cout << " Laser shots:" << setw(6) << NumberOfLaserShots 
                     << " t =" << setw(11) << setprecision(4) << currentTime
                     << " iter = " << setw(12) << currentTimeiteration
                     << " z = " << setw(12) << setprecision(6) << currentSamplePosition
                     << " Progress = " << setprecision(4) << CalculationProgress << " %\n";
                NumberOfLaserShots++;
                ShotLaser = 0;
            }
            
            if (SavePlanes == SaveZPRIteration)
            {
                //SavePlanes function launch
                SaveAllCenterProfiles(SaveAllCenterProfilesIndex);
                SaveAllCenterProfilesIndex++;

                SaveAllCenterHeatMapSurfaces(SaveAllCenterHeatMapSurfacesIndex);
                SaveAllCenterHeatMapSurfacesIndex++;

                CalculationProgress = (currentTimeiteration * 100.0) / NumberOfTimeIterations;
                        
                cout << " Saving all planes! t =" << setw(11) << setprecision(4) << currentTime 
                    << " iter = "  << setw(12) << currentTimeiteration
                    << " z = " << setw(12) << setprecision(6) << currentSamplePosition 
                    << " Progress = " << setprecision(4) << CalculationProgress << " %\n";
                SavePlanes = 0;
            }

            if (SaveHeatMap == SaveAllSampleIteration)
            {
                SaveAllSAmpleHeatMap(currentTimeiteration, SaveAllSAmpleHeatMapIndex, HeatMapFileFormat);
                SaveAllSAmpleHeatMapIndex++;

                CalculationProgress = (currentTimeiteration * 100.0) / NumberOfTimeIterations;

                cout << " Saving heat map!   t =" << setw(11) << setprecision(4) << currentTime 
                    << " iter = " << setw(12) <<currentTimeiteration
                    << " z = " << setw(12) << setprecision(6) <<currentSamplePosition
                    << " Progress = " << setprecision(4) << CalculationProgress << " %\n";
                SaveHeatMap = 0;
            }

            
            if (SaveSampleCollectedEnergy == SampleCollectedEnergyIterations)
            {
                SaveSampleEnergyFile();
                SaveSampleCollectedEnergy = 0;
            }
            
            SaveSampleCollectedEnergy++;
            ShotLaser++;
            SavePlanes++;
            SaveHeatMap++;

            currentTimeiteration++;
            currentSamplePosition += deltaSamplePositionSlowIteration;
            if (currentSamplePosition > (exData->GetZScanRange() / 2))
            {
                stopshoot = true;
                currentSamplePosition = (exData->GetZScanRange() / 2);
            }
                
            currentTime += exData->Getdeltat();

        } while (currentTimeiteration < NumberOfTimeIterations);
        


        //  -----------------------------
        //           Final saves
        // ------------------------------
        cout << endl << endl;
        cout << " Final saves:\n\n";
        CalculationProgress = (currentTimeiteration * 100.0) / NumberOfTimeIterations;

        cout << " Save all planes ! "
            << currentTime << " iter = "
            << currentTimeiteration << " z = "
            << currentSamplePosition
            << " Calculation Progress = " << setprecision(3) << CalculationProgress << " %\n";
        SavePlanes = 0;

        SaveAllCenterProfiles(SaveAllCenterProfilesIndex);
        SaveAllCenterProfilesIndex++;

        SaveAllCenterHeatMapSurfaces(SaveAllCenterHeatMapSurfacesIndex);
        SaveAllCenterHeatMapSurfacesIndex++;

        cout << " Save heat map ! t = "
            << currentTime << " iter = " << currentTimeiteration
            << " z = " << currentSamplePosition
            << " Calculation Progress = " << setprecision(3) << CalculationProgress << " %\n";
        SaveAllSAmpleHeatMap(currentTimeiteration, SaveAllSAmpleHeatMapIndex, HeatMapFileFormat);
        SaveAllSAmpleHeatMapIndex++;

        cout << "\n Iteration: " << currentTimeiteration 
             << "\n z = " << currentSamplePosition 
             << "\n Time = " << currentTime 
             << "\n\n Z-scan thermoanalysis complete. Calculations End.\n";

        // Write some funny text at the end of calculations:\

        srand(time(NULL));
        int num = rand() % 10 + 1;
        
        if (num == 1) {
            cout << " All your base are belong to us... and your calculation is complete!" << endl;
        }
        else if (num == 2) {
            cout << " May the source be with you, for your calculation is complete!" << endl;
        }
        else if (num == 3) {
            cout << " Congratulations! You have unlocked the achievement: 'Master of the Lambda'!" << endl;
        }
        else if (num == 4) {
            cout << " The cake may be a lie, but your calculation is not. It's complete!" << endl;
        }
        else if (num == 5) {
            cout << " To infinity and beyond! Your calculation has successfully reached its destination." << endl;
        }
        else if (num == 6) {
            cout << " It's dangerous to go alone! Take this... successfully completed calculation." << endl;
        }
        else if (num == 7) {
            cout << " Great Scott! Your calculation has successfully traveled through time and space." << endl;
        }
        else if (num == 8) {
            cout << " Congratulations, you have leveled up in Lambda Science! Your calculation is complete." << endl;
        }
        else if (num == 9) {
            cout << " It's official, you are a Lambda Ninja. Your calculation skills are unmatched!" << endl;
        }
        else {
            cout << " Congratulations! You have completed your calculation in record time. You're a wizard, Harry!" << endl;
        }



    }

     //DeallocateRAM(); - be carefull !!!
}

void SimulationData::DeallocateRAM()
{
    //Deallocate RAM
    cout << "\n\n Deallocate RAM: " << flush;
    Delete3D(Temp_n);
    Delete3D(Temp_np1);

    //delete[] GetHeat_expArray;
    //delete[] GetHeat_omegaArgument;
    //delete[] LPIExpTable;

    cout << " Done.\n \n" << flush;
}

double SimulationData::CornerThreeD(double dXi, double dYi, double dZi, double Tcz[3], double Tijk, double Ti_1, double Tj_1, double Tk_1, double dXr, double dYr, double dZr)
{
    double Two_dFox = 2 * dFox;
    double Two_dFoy = 2 * dFoy;
    double Two_dFoz = 2 * dFoz;
    
    double result = 0;
    result =  Two_dFox * (Ti_1 + dXi * Tcz[0]);
    result += Two_dFoy * (Tj_1 + dYi * Tcz[1]);
    result += Two_dFoz * (Tk_1 + dZi * Tcz[2]);
    result += Tijk * (1 - (Two_dFox) - (Two_dFoy) - (Two_dFoz) - (Two_dFoz * dZi) - (Two_dFoy * dYi) - (Two_dFox * dXi));
    

    return result;
}

double SimulationData::SurfaceThreeD(double dFo_a, double dFo_b, double dFo_c,
                     double Tijk, double Ta_1, double Tap1,
                     double Tb_1, double Tbp1, double Tcp1,
                     double dCi, double dCr, double Tcz)
{
    double result = 0;
    
    result  = Tijk * ( 1 - (2*dFo_a) - (2*dFo_b) - (2*dFo_c) - (2*dFo_c*dCi) - (2*dFo_c*dCr*Tijk*Tijk*Tijk));
    result += dFo_a*(Ta_1 + Tap1);
    result += dFo_b*(Tb_1 + Tbp1);
    result += 2*dFo_c*(dCi*Tcz + dCr*Tijk*Tijk*Tijk*Tijk + Tcp1);
   
    return result;
    
}

void SimulationData::InitaiateCenterProfiles()
{
    ofstream file;
    
    //First step - erase all data

    file.open("Results/XProfile.dat", std::ofstream::out | std::ofstream::trunc);
    file.close();

    file.open("Results/YProfile.dat", std::ofstream::out | std::ofstream::trunc);
    file.close();

    file.open("Results/ZProfile.dat", std::ofstream::out | std::ofstream::trunc);
    file.close();

    //Second step - open for append

    // file.open("XProfile.dat", std::ios_base::app); // append instead of overwrite
    // file.open("YProfile.dat", std::ios_base::app); // append instead of overwrite
    // file.open("ZProfile.dat", std::ios_base::app); // append instead of overwrite


}

void SimulationData::SaveAllCenterProfiles(unsigned int Index)
{
    ofstream fileX;
    ofstream fileY;
    ofstream fileZ;
    double dist;

    fileX.open("Results/XProfile.dat", std::ios_base::app); // append instead of overwrite
    fileY.open("Results/YProfile.dat", std::ios_base::app); // append instead of overwrite
    fileZ.open("Results/ZProfile.dat", std::ios_base::app); // append instead of overwrite

    if (!fileX) {
        cerr << "Can't save X profile in: XProfile.dat\n";
    }
    else
    {
        fileX << "# Index = " << Index << endl;
        fileX << "# Iteration = " << currentTimeiteration << endl;
        fileX << "# Time = " << currentTime << endl;
        fileX << "# Sample_position = " << currentSamplePosition << endl;
        fileX << "# Node     XDim       XCDim       Temperature\n";

        for (unsigned int i = 0; i < N; i++)
        {
            dist = (i * spData->GetdX());
            fileX << fixed << setw(8) << setprecision(1) << i;
            fileX << setw(12) << setprecision(8) << dist;
            fileX << setw(12) << setprecision(8) << dist - spData->GetXmeters() / 2;
            fileX << setw(18) << setprecision(12) << Temp_n[i][halfO][halfP] << endl;
        }
        fileX << endl << endl;
        fileX.close();
    }

    if (!fileY) {
        cerr << " Can't save Y profile in: YProfile.dat\n";
    }
    else
    {
        fileY << "# Index = " << Index << endl;
        fileY << "# Iteration = " << currentTimeiteration << endl;
        fileY << "# Time = " << currentTime << endl;
        fileY << "# Sample_position = " << currentSamplePosition << endl;
        fileY << "# Node     YDim       YCDim       Temperature\n";

        for (unsigned int i = 0; i < O; i++)
        {
            dist = (i * spData->GetdY());
            fileY << fixed << setw(8) << setprecision(1) << i;
            fileY << setw(12) << setprecision(8) << dist;
            fileY << setw(12) << setprecision(8) << dist - spData->GetYmeters() / 2;
            fileY << setw(18) << setprecision(12) << Temp_n[halfN][i][halfP] << endl;
        }

        fileY << endl << endl;
        fileY.close();
    }

    if (!fileZ) {
        cerr << " Can't save Z profile in: ZProfile.dat\n";
    }
    else
    {
        fileZ << "# Index = " << Index << endl;
        fileZ << "# Iteration = " << currentTimeiteration << endl;
        fileZ << "# Time = " << currentTime << endl;
        fileZ << "# Sample_position = " << currentSamplePosition << endl;
        fileZ << "# Node     ZDim       ZCDim       Temperature\n";

        for (unsigned int i = 0; i < P; i++)
        {
            dist = (i * spData->GetdZ());
            fileZ << fixed << setw(8) << setprecision(1) << i;
            fileZ << setw(12) << setprecision(8) << dist;
            fileZ << setw(12) << setprecision(8) << dist - spData->GetZmeters() / 2;
            fileZ << setw(18) << setprecision(12) << Temp_n[halfN][halfO][i] << endl;
        }

        fileZ << endl << endl;
        fileZ.close();
    }
}

void SimulationData::InitiateCenterHeaatMapSurfaces(unsigned int fileindex)
{
    
    //First step - erase all data
    ofstream file;
    string filename = "Results/Surfaces/CenterSurfaceXY_";
    filename.append(to_string(fileindex));
    filename.append(".dat");

    file.open(filename, std::ofstream::out | std::ofstream::trunc);
    file.close();


    filename = "Results/Surfaces/CenterSurfaceXZ_";
    filename.append(to_string(fileindex));
    filename.append(".dat");
    file.open(filename, std::ofstream::out | std::ofstream::trunc);
    file.close();

    
    filename = "Results/Surfaces/CenterSurfaceYZ_";
    filename.append(to_string(fileindex));
    filename.append(".dat");
    file.open(filename, std::ofstream::out | std::ofstream::trunc);
    file.close();
    
}

void SimulationData::SaveAllCenterHeatMapSurfaces(unsigned int Index)
{
    ofstream fileXY;
    ofstream fileXZ;
    ofstream fileYZ;
    string pathXY = "Results/Surfaces/CenterSurfaceXY_";
    string pathXZ = "Results/Surfaces/CenterSurfaceXZ_";
    string pathYZ = "Results/Surfaces/CenterSurfaceYZ_";

    if (((Index+1) % 100) == 0)
    {
        SavedCenterSurfacesCounter++;
        pathXY.append(to_string(SavedCenterSurfacesCounter));
        pathXZ.append(to_string(SavedCenterSurfacesCounter));
        pathYZ.append(to_string(SavedCenterSurfacesCounter));
        pathXY.append(".dat");
        pathXZ.append(".dat");
        pathYZ.append(".dat");
        InitiateCenterHeaatMapSurfaces(SavedCenterSurfacesCounter);


    }
    else
    {
        pathXY.append(to_string(SavedCenterSurfacesCounter));
        pathXZ.append(to_string(SavedCenterSurfacesCounter));
        pathYZ.append(to_string(SavedCenterSurfacesCounter));
        pathXY.append(".dat");
        pathXZ.append(".dat");
        pathYZ.append(".dat");
    }

    /*
    fileXY.open("Results/CenterSurfaceXY.dat", std::ios_base::app); // append instead of overwrite
    fileXZ.open("Results/CenterSurfaceXZ.dat", std::ios_base::app); // append instead of overwrite
    fileYZ.open("Results/CenterSurfaceYZ.dat", std::ios_base::app); // append instead of overwrite
    */

    fileXY.open(pathXY, std::ios_base::app); // append instead of overwrite
    fileXZ.open(pathXZ, std::ios_base::app); // append instead of overwrite
    fileYZ.open(pathYZ, std::ios_base::app); // append instead of overwrite

    if (!fileXY) {
        cerr << "Can't save XY surface in: CenterSurfaceXY.dat\n";
    }
    else
    {
        fileXY << "# Index = " << Index << endl;
        fileXY << "# Iteration = " << currentTimeiteration << endl;
        fileXY << "# Time = " << currentTime << endl;
        fileXY << "# Sample_position = " << currentSamplePosition << endl;
        fileXY << "# Number of nodes:\n";
        fileXY << "# X = " << N << " Y = " << O << endl;
        fileXY << "# Node sizes (m):\n";
        fileXY << "# dX = " << spData->GetdX() << " dY = " << spData->GetdY() << endl;
        fileXY << "# Rows -> X; Columns -> Y\n";

        for (unsigned int j = 0; j < O; j++)//reverse up and down
        {
            for (unsigned int i = 0; i < N; i++)
            {
                fileXY << fixed << setprecision(10) << Temp_n[i][j][halfP] << " ";
            }
            fileXY << endl;
        }
    }

    if (!fileXZ) {
        cerr << "Can't save XZ surface in: CenterSurfaceXZ.dat\n";
    }
    else
    {
        fileXZ << "# Index = " << Index << endl;
        fileXZ << "# Iteration = " << currentTimeiteration << endl;
        fileXZ << "# Time = " << currentTime << endl;
        fileXZ << "# Sample_position = " << currentSamplePosition << endl;
        fileXZ << "# Number of nodes:\n";
        fileXZ << "# X = " << N << " Z = " << P << endl;
        fileXZ << "# Node sizes (m):\n";
        fileXZ << "# dX = " << spData->GetdX() << " dZ = " << spData->GetdZ() << endl;
        fileXZ << "# Rows -> X; Columns -> Z\n";

        for (unsigned int k = 0; k < P; k++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                fileXZ <<  fixed << setprecision(10) << Temp_n[i][halfO][k] << " ";
            }
            fileXZ << endl;
        }

    }

    if (!fileYZ) {
        cerr << "Can't save YZ surface in: CenterSurfaceYZ.dat\n";
    }
    else
    {
        fileYZ << "# Index = " << Index << endl;
        fileYZ << "# Iteration = " << currentTimeiteration << endl;
        fileYZ << "# Time = " << currentTime << endl;
        fileYZ << "# Sample_position = " << currentSamplePosition << endl;
        fileYZ << "# Number of nodes:\n";
        fileYZ << "# Y = " << O << " Z = " << P << endl;
        fileYZ << "# Node sizes (m):\n";
        fileYZ << "# dY = " << spData->GetdY() << " dZ = " << spData->GetdZ() << endl;
        fileYZ << "# Rows -> Y; Columns -> Z\n";

        for (unsigned int k = 0; k < P; k++)
        {
            for (unsigned int j = 0; j < O; j++) //reverse up and down
            {
                fileYZ <<  fixed << setprecision(10) <<Temp_n[halfN][j][k] << " ";
            }
            fileYZ << endl;
        }

    }

    fileXY << endl << endl;
    fileXZ << endl << endl;
    fileYZ << endl << endl;

    fileXY.close();
    fileXZ.close();
    fileYZ.close();
}

void SimulationData::SaveAllSAmpleHeatMap(unsigned int IterationForFileName, unsigned int Index, unsigned int DataType)
{
    // DataType:
    // 1 - seven columns:  i,j,k x, y, z, T specebar separated (.dat)
    // 2 - four  columns: x, y, z T comma separated (.csv)
    string path = "Results/HeatMaps/";

    if (DataType == 1)
    {
        path.append(to_string(IterationForFileName));
        path.append("HeatMap.dat");
        ofstream HeatMapfile(path, ios::out);


        if (!HeatMapfile) {
            cerr << "Can't save Heat map in: HeatMap.dat\n";
        }
        else
        {
            double Xdist, Ydist, Zdist;
            double XdistHalf, YdistHalf, ZdistHalf;
            double XdistHalf_a, YdistHalf_a, ZdistHalf_a;

            XdistHalf_a = (spData->GetXmeters() / 2);
            YdistHalf_a = (spData->GetYmeters() / 2);
            ZdistHalf_a = (spData->GetZmeters() / 2);


            HeatMapfile << "# Index = " << Index << endl;
            HeatMapfile << "# Iteration = " << currentTimeiteration << endl;
            HeatMapfile << "# Time = " << currentTime << endl;
            HeatMapfile << "# Sample_position = " << currentSamplePosition << endl;
            HeatMapfile << "# Xnode Ynode Znode XDim YDim ZDim Temperature\n";

            for (unsigned int k = 0; k < P; k++) {

                Zdist = (k * spData->GetdZ());
                ZdistHalf = Zdist - ZdistHalf_a;

                for (unsigned int j = 0; j < O; j++) {

                    Ydist = (j * spData->GetdY());
                    YdistHalf = Ydist - YdistHalf_a;

                    for (unsigned int i = 0; i < N; i++)
                    {
                        Xdist = (i * spData->GetdX());
                        XdistHalf = Xdist - XdistHalf_a;

                        HeatMapfile << fixed << setw(5) << setprecision(1) << i;
                        HeatMapfile << fixed << setw(5) << setprecision(1) << j;
                        HeatMapfile << fixed << setw(5) << setprecision(1) << k;

                        //HeatMapfile << setw(12) << setprecision(8) << Xdist;
                        //HeatMapfile << setw(12) << setprecision(8) << Ydist;
                        //HeatMapfile << setw(12) << setprecision(8) << Ydist;

                        HeatMapfile << setw(9) << setprecision(4) << XdistHalf;
                        HeatMapfile << setw(9) << setprecision(4) << YdistHalf;
                        HeatMapfile << setw(9) << setprecision(4) << ZdistHalf;

                        HeatMapfile << setw(20) << setprecision(15) << Temp_n[i][j][k] << endl;
                    }

                }
            }
        }
        
        HeatMapfile.close();
    }

    if (DataType == 2)
    {
        path.append("HeatMap.csv.");
        path.append(to_string(SavedHeatMapsCounter));
        ofstream HeatMapfile(path, ios::out);

        if (!HeatMapfile) {
            cerr << "Can't save Heat map in: HeatMap.csv\n";
        }
        else
        {
            double Xdist, Ydist, Zdist;
            double XdistHalf, YdistHalf, ZdistHalf;
            double XdistHalf_a, YdistHalf_a, ZdistHalf_a;

            XdistHalf_a = (spData->GetXmeters() / 2);
            YdistHalf_a = (spData->GetYmeters() / 2);
            ZdistHalf_a = (spData->GetZmeters() / 2);


         
            HeatMapfile << "X" << N << ", Y" << O << ", Z" << P << ", Temperature, Sample Position \n";
   

            for (unsigned int k = 0; k < P; k++) {

                Zdist = (k * spData->GetdZ());
                ZdistHalf = Zdist - ZdistHalf_a;

                for (unsigned int j = 0; j < O; j++) {

                    Ydist = (j * spData->GetdY());
                    YdistHalf = Ydist - YdistHalf_a;

                    for (unsigned int i = 0; i < N; i++)
                    {

                        Xdist = (i * spData->GetdX());
                        XdistHalf = Xdist - XdistHalf_a;

                        //HeatMapfile << setw(12) << setprecision(8) << Xdist;
                        //HeatMapfile << setw(12) << setprecision(8) << Ydist;
                        //HeatMapfile << setw(12) << setprecision(8) << Ydist;

                        HeatMapfile << setw(11) << setprecision(4) << XdistHalf << ",";
                        HeatMapfile << setw(11) << setprecision(4) << YdistHalf << ",";
                        HeatMapfile << setw(11) << setprecision(4) << ZdistHalf << ",";

                        HeatMapfile << setw(20) << setprecision(15) << Temp_n[i][j][k] << ", ";
                        HeatMapfile << currentSamplePosition << endl;
                    }

                }
            }
        }
        SavedHeatMapsCounter++;
        HeatMapfile.close();
    }

}

void SimulationData::InitiateSampleEnergyFile()
{
    ofstream file;

    file.open("Results/SampleEnergy.dat", std::ofstream::out | std::ofstream::trunc);
    file << "#  Sample    Simulat.   Sample\n";
    file << "# Position     Time     Energy\n";
    file << "#   (mm)       (ms)      (mJ)\n";
    file.close();
}

void SimulationData::SaveSampleEnergyFile()
{
    ofstream file;
    file.open("Results/SampleEnergy.dat", std::ios_base::app); // append instead of overwrite
    double Energy;

    if (!file) {
        cerr << "Can't save Energy report in: Results/SampleEnergy.dat\n";
    }
    else
    {
        Energy = GetSampleCollectedEnergy();
        
        file << fixed << setw(10) << setprecision(4) << (currentSamplePosition*1000) << " ";
        file << fixed << setw(11) << setprecision(4) << (currentTime*1000) << " ";
        file << fixed << setw(12) << setprecision(8) << (Energy*1000) << "\n";
    }
}

bool SimulationData::SaveSurfaceX(unsigned int layernumberN, char* SurfaceFilename)
{
    //string path = ".\\results\\";
    string path = "";
    path.append(SurfaceFilename);
    path.append(".txt");
    ofstream surfacefile(path, ios::out);

    if (!surfacefile) {
        cerr << "Can't save X surface file: " << path << endl;
        return true;
    }
    else
    {
        for (unsigned int k = 0; k < P; k++)
        {
            for (unsigned int j = 0; j < O; j++)
            {
                surfacefile << Temp_n[layernumberN][j][k] << " ";
            }
            surfacefile << "\n";
        }
        surfacefile.close();
        return false;
    }
}

bool SimulationData::SaveSurfaceZ(unsigned int layernumberP, char* SurfaceFilename)
{
    //string path = ".\\results\\";
    string path = "";
    path.append(SurfaceFilename);
    path.append(".txt");
    ofstream surfacefile(path, ios::out);


    if (!surfacefile) {
        cerr << "Can't save Z surface file: " << path << endl;
        return true;
    }
    else
    {
        for (unsigned int j = 0; j < O; j++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                surfacefile << Temp_n[i][j][layernumberP] << " ";
            }
            surfacefile << "\n";
        }
        surfacefile.close();

        return false;

    }
}

bool SimulationData::SaveSurfaceY(unsigned int layernumberO, char* SurfaceFilename)
{
    //string path = ".\\results\\";
    string path = "";
    path.append(SurfaceFilename);

    path.append(".txt");
    ofstream surfacefile(path, ios::out);


    if (!surfacefile) {
        cerr << "Can't save Y surface file: " << path << endl;
        return true;
    }
    else
    {
        for (unsigned int k = 0; k < P; k++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                surfacefile << Temp_n[i][layernumberO][k] << " ";
            }
            surfacefile << "\n";
        }
        surfacefile.close();

        return false;

    }
}



double SimulationData::EdgeThreeD(double dFo_a, double dFo_b, double dFo_c, double dBi, double dCi, double Tcz[2], double Tijk, double Ta_1, double Tap1, double Tb_1, double Tc_1, double dBr, double dCr)
{
    double result = 0;
    double Tijk_3;
    double T_4_b, T_4_c;
    T_4_b = Tcz[0] * Tcz[0] * Tcz[0] * Tcz[0];
    T_4_c = Tcz[1] * Tcz[1] * Tcz[1] * Tcz[1];
    
   
    Tijk_3 = Tijk* Tijk* Tijk;


    result =  dFo_a * (Ta_1 + Tap1) + (2 * dFo_b * (Tb_1 + dBi * Tcz[0] + dBr * T_4_b));
    result += 2 * dFo_c * (Tb_1 + dCi * Tcz[1] +  dCr * T_4_c );
    result += Tijk * (1 - 2*dFo_a - 2*dFo_b - 2*dFo_c - (2*dFo_b*dBi) - (2*dFo_c*dCi) - (2 * dFo_b * dBr * Tijk_3) -   (2 * dFo_c * dCr * Tijk_3));
    
    if (result == -nan("ind"))
    {
        cout << currentTimeiteration << endl;
        double foo;
        cin >> foo;
    }

    return result;
   
}
