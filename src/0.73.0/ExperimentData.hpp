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
#ifndef EXPERIMENTDATA_H_INCLUDED
#define EXPERIMENTDATA_H_INCLUDED

class ExperimentData
{
public:
	ExperimentData();
	~ExperimentData();

	void ShowMeExperiment();
	
	int GetSingleShot();

	double GetZScanRange();
	double GetSampleVeolocity();
	
	double GetSimulationTime();
	double Getdeltat();
	double GetdeltatSlowMotion();


	unsigned int GetTimeSteps();

	double GetLaserWavelength();
	double GetLaserPulseDuration();
	unsigned int GetLaserPulseRepetition();
	double GetLaserPulseEnergy();
	double GetLaserCWPower();

	double GetomegaZero();

	unsigned short int GetThreads();

	double GetLaserIrradianceI0inGWcm2();

	double GetSaveZPResultsInterval();
	double GetSaveAllResultsInterval();

    bool InitiateError = false;

private:
	void setdeltat();

	int SingleShot = -1;
	

	double ZScanRange = 0;
	bool isZScanRange = false;

	double SampleVelocity = 0;
	bool isSampleVelocity = false;
	
	double SimulationTime = 0;
	bool isSimulationTime = false;



	double NumericalTimeStep = 0;
	bool isNumericalTimeStep = false;

	double deltat = 0;
	double deltatSlowMontion = 0;

	double LambdaOfLaser = 0;
	bool isLambdaOfLaser = false;

	double LaserPulseDuration = -10.0;
	bool isLaserPulseDuration = false;
	bool isPulsedLaser = false;			//Check that the laser is pulsed
	bool isCWLaser = false;				//or CW

	double LaserPulseEnergy = 0.0;		//(J)
	bool isLaserPulseEnergy = false;

	double LaserCWPower = 0.0;			//(W)
	bool isLaserCWPower = false;

	unsigned int LaserPulseRepetition = 0;
	bool isLaserPulseRepetition = false;

	double omegaZero = 1.0e-9;
	bool isomegaZero = false;

	double LaserIrradianceI0 = 0;

	unsigned short int Threads = 1;
	bool isThreads = false;
	
	double SaveZPResultsInterval = 0.01;
	bool isSaveZPResultsInterval = false;
	int DoSaveZeroPlane = false;
	bool isDoSaveZeroPlane = false;

	double SaveAllResultsInterval = 0.05;
	int isSaveAllResultsInterval = false;
	int DoSaveAllResults = false;
	bool isDoSaveAllResults = false;
};

class SimulationData
{
public:
	void SetData(ExperimentData incomingexData, SampleData incomingsmData);
	bool ShowMeSimulationData(ExperimentData incomingexData, SampleData incomingsmData);
	void StartCalculation(unsigned int NumberOfTimeIterations);

	void DeallocateRAM();
	
	bool SaveSurfaceX(unsigned int layernumberN, char* SurfaceFilename);
    bool SaveSurfaceY(unsigned int layernumberO, char* SurfaceFilename);
    bool SaveSurfaceZ(unsigned int layernumberP, char* SurfaceFilename);

	SimulationData();
	~SimulationData();

private:
	void CalculateCornersTemperature();
	void CalculateEdgeTemperatures();
	void CalculateCoreTemperature();
	void CalculateCoreTemperatureThreads(unsigned short int starti, unsigned short int endi);
    void CalculateWallsTemperature();

	const double pi = 3.14159265358979323846;

	double GetDiameterAtZposition(double Z);
	double GetZr();
	double GetZrSquare();
	double GetOmega(double z);

	bool CheckCalculationStability();

	void SetNormalTimeSteps();
	void SetShortTimeSteps();

	double ShootToThrill(double Z);
	double GetHeat(double x, double y, double z, double LaserPower, unsigned int sampleZnode);
	double GetSampleCollectedEnergy();

	double LaserPulseIntensity(int currnetStep, double Amplitude);

	void CalculateAll();

	void CopyArray();

	void ShowCalculationProgress(unsigned int currentIteration, unsigned int TotalIterations);

    double CornerThreeD(double dXi, double dYi, double dZi, double Tcz[3], double Tijk, double Ti_1, double Tj_1, double Tk_1, double dXr, double dYr, double dZr);
    double EdgeThreeD(double dFo_a, double dFo_b, double dFo_c, double dBi, double dCi, double Tcz[2], double Tijk, double Ta_1, double Tap1, double Tb_1, double Tc_1, double dBr, double dCr);
    double SurfaceThreeD(double dFo_a, double dFo_b, double dFo_c, double Tijk, double Ta_1, double Tap1, double Tb_1, double Tbp1, double Tcp1, double dCi, double dCr, double Tcz);
    //HeatMap file format
	unsigned int HeatMapFileFormat = 2; //1 - 7 columns *.dat files, 2 four columns *.csv files


	//Multithread settings
	bool isMultiThread = false;
	unsigned short int ThreadNumber = 1;
	unsigned short int ipart1, ipart2, ipart3;		//For 4threads fragmentation
	unsigned short int ipart = 0;					//For N threads
	unsigned short int istart;						//For N threads
	unsigned short int iend;						//For N threads


	//-----------------------------------------------------------------------------------
	//                             Function for data storage
	//-----------------------------------------------------------------------------------

	//Functions for data storages in X,Y and ZProfile.dat files
	//
	void InitaiateCenterProfiles();			 //Erase all data in X,Y and ZProfile.dat files
	void SaveAllCenterProfiles(unsigned int Index);			 //Append data to X,Y and ZProfile.dat files

	void InitiateCenterHeaatMapSurfaces(unsigned int fileindex);	//Erase all data in  X, Y A ZHMSurfaces.dat
	void SaveAllCenterHeatMapSurfaces(unsigned int Index);	//Append data to X, Y A ZHMSurfaces.dat

	//All Samples Node Saving Function
	void SaveAllSAmpleHeatMap(unsigned int IterationForFileName, unsigned int Index, unsigned int DataType);	//Save data in XXXXXHeatMap.dat file, where XXXXX - is IterationForFileName for DataType == 1
																												// or       in XXXXXHeatMap.csv file, where XXXXX - is IterationForFileName for DataType == 2
	unsigned int SavedHeatMapsCounter = 0;
	unsigned int SavedCenterSurfacesCounter = 0;


	void InitiateSampleEnergyFile();
	void SaveSampleEnergyFile();	//Save Gained or losses Energy by the sample  dada are stored in ene

	

	ExperimentData* exData;
	SampleData* spData;

	int SaveZPRIteration = 1;
	int SaveAllSampleIteration = 1;
	int LaserRepetitionIterations = 1;


	double dFox = 0;	//current time steps
	double dFoy = 0;
	double dFoz = 0;

	double dFoxNormalMotion = 0; //for normal calculations
	double dFoyNormalMotion = 0;
	double dFozNormalMotion = 0;


	double dFoxSlowMotion = 0; //for calculation in time steps around pulse width
	double dFoySlowMotion = 0;
	double dFozSlowMotion = 0;

	
	double Temp0 = 298;
	
	double Tamb_zF = 298;                  //in [K]
	double Tamb_zB = 298;                  //in [K]

	double Tamb_xR = 298;                  //in [K]
	double Tamb_xL = 298;                  //in [K]

	double Tamb_yU = 298;                  //in [K]
	double Tamb_yD = 298;                  //in [K]

	//External convective heat transfer from the axis 
	// X - dXiR (Right) and dXiL (Left)
	// Y - dYiU (Upper) oraz dYiL (Lower)
	// Z - dZiF (Front) oraz dZiB (Back)
	double dXiR = 0;
	double dXiL = 0;

	double dYiU = 0;
	double dYiL = 0;

	double dZiF = 0;
	double dZiB = 0;
    
    //Parameters related to radiadive heat transfer:
    double dXrR = 0;
    double dXrL = 0;

    double dYrU = 0;
    double dYrL = 0;

    double dZrF = 0;
    double dZrB = 0;

	double UnitCellVolume = 0;

	unsigned int currentTimeiteration = 0;

	double SimulationTime = 0;

	bool isSingleShot = false;
	
	char ProjectFolder;

	const char kPathSeparator =
	#if defined _WIN32 || defined __CYGWIN__
		'\\';
	#else
		'/';
	#endif

	//Arrays Dimensions (Number of Sample Nodes) in X, Y and Z axis
	unsigned int N, O, P;

	unsigned int halfN, halfO, halfP;

	double*** Temp_n;
	double*** Temp_np1;

	//Variablies for faster calculations

	double LPIExpTable[51];		//Laser Pulse Intensity exp(-0.5 * square) table 
								//where square = ((currnetStep - 25.0) / 5.0) * ((currnetStep - 25.0) / 5.0);
	int Time1000iteration = 0;
	
	// ----------------------- GetHeat Function some variables    ------------------------------
	double* GetHeat_expArray;			//Array of exp(-alpha * sampleZnode * spData->GetdZ()); in GetHeat function;
	double* GetHeat_omegaArgument;		//Array of (0.5 * spData->GetZmeters()) + sampleZnode * spData->GetdZ()) argument
	double betaDivZmeters = 0;  

	double omegaZeroSquare = 0.0;
	double Zrsquare = 0.0;

	//Variables for results data files



	//Variables used in void StartCalculation(unsigned int NumberOfTimeIterations) - function
	double currentTime = 0;
	double deltaSamplePositionSlowIteration = 0; 
	double deltaSamplePositionLaserIteration = 0;
	double currentSamplePosition = 0;
	
		
};


#endif

