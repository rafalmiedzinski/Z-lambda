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

#include "thermo.hpp"
#include "ExperimentData.hpp"
#include <iostream>
#include <thread> 
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

int main()
{
    //SampleData Sample;
    SampleData* pSample = new SampleData;

    //ExperimentData Experiment;
    ExperimentData* pExperiment = new ExperimentData;
    
    SimulationData Simulation;

    if (pSample->InitiateError) return 1;
    //pSample->ShowMeASample();

    if (pExperiment->InitiateError) return 1;

    Simulation.SetData (*pExperiment, *pSample);
    if(Simulation.ShowMeSimulationData(*pExperiment, *pSample)) {
        cout << " Starting the calculations";
    }
    else
    {
        cout << "\n Honestly - Something went wrong.\n";
        Simulation.DeallocateRAM();
        return 1;
    }
    
    Simulation.StartCalculation(pExperiment->GetTimeSteps());
    Simulation.DeallocateRAM();

    return 0;
}
