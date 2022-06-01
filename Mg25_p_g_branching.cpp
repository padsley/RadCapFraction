// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <cmath>
// #include <algorithm>

#include "RadCapFunctions.h"

// vector<double> vEx;
// 
// int nStates = 0;
// 
// bool VerboseFlag = false;
// void PrintIntensityVector(double *IntensityVector);
// void PrintIntensityMatrix(double **IntensityMatrix);
// 
// bool CheckFinished(double *IntensityVector);
// 
// void LoadLevelInformation();
// void LoadDecayInformation();
// void NormaliseDecayInformation(double **IntensityMatrix);
// 
// double **IntensityMatrix;
// double **sigIntensityMatrix;
// double *IntensityVector;
// double *nextIntensityVector;
// double *sigIntensityVector;
// double *nextsigIntensityVector;

// double CorrectForRecoil(double Egamma);

void Mg25_p_g_branching(int InitialStateIndex = 2)
{
    if(vEx.size()==0)
    {
        LoadLevelInformation();
    }
    
    //= new double*[nStates];
    IntensityMatrix    = new double*[nStates];
    sigIntensityMatrix = new double*[nStates];
    
    IntensityVector = new double[nStates];
    nextIntensityVector = new double[nStates];
    sigIntensityVector = new double[nStates];
    nextsigIntensityVector = new double[nStates];
    
    for(int i=0;i<nStates;i++){
        IntensityVector[i] = 0;
        nextIntensityVector[i] = 0;
        sigIntensityVector[i] = 0;
        nextsigIntensityVector[i] = 0;
        
        IntensityMatrix[i] = new double[nStates];
        sigIntensityMatrix[i] = new double[nStates];
        for(int j=0;j<nStates;j++)
        {
            IntensityMatrix[i][j] = 0;
            sigIntensityMatrix[i][j] = 0;
        }
    }
    
    IntensityVector[InitialStateIndex] = 1.0;
    
    if(IntensityMatrix[0][2]==0)
    {
        LoadDecayInformation();
    }
    
    NormaliseDecayInformation(IntensityMatrix);
    //     PrintIntensityMatrix(IntensityMatrix);
    
    //IntensityMatrix[initial state][final state]
    int count_steps = 0;
    while(!CheckFinished(IntensityVector,2))
    {
//         std::cout << "Step " << count_steps << endl;
        
                        if(VerboseFlag)
        {
            std::cout << "IntensityVector" << std::endl;
            PrintIntensityVector(IntensityVector,InitialStateIndex);
        }
                        if(VerboseFlag)
        {
            std::cout << "nextIntensityVector" << std::endl;
            PrintIntensityVector(nextIntensityVector,InitialStateIndex);
        }
        
        
        for(int j=2;j<nStates;j++)//this is the initial state
        {
            //             int j=0;
            //             if(i<2)j=2;
            //             else j=i+1;
            
//             if(IntensityVector[j]>0)std::cout << "Intensity vector: " << j << "\t" << IntensityVector[j] << "\t" << sigIntensityVector[j] << "\t" << 100*sigIntensityVector[j]/IntensityVector[j] << std::endl;
            
            for(int i=0;i<nStates;i++)//this is the final state
            {
                //                 if(IntensityVector[j]>0)std::cout << "i=" << i << "\t" << "j=" << j << "\t" << IntensityMatrix[j][i] << "\t" << IntensityVector[j] << std::endl;
                if(IntensityMatrix[j][i]>0 && IntensityVector[j]>0)
                {
//                     std::cout << "Final state: " << i << "\t" << "Initial state: " << j << "\t" << IntensityMatrix[j][i] << "\t" << sigIntensityMatrix[j][i] << std::endl;
                    
                    nextIntensityVector[i] += IntensityMatrix[j][i] * IntensityVector[j];
//                     nextsigIntensityVector[i] = nextIntensityVector[i] * sqrt(pow(sigIntensityVector[j]/IntensityVector[j],2.) + pow(sigIntensityMatrix[j][i]/IntensityMatrix[j][i],2.));
                    nextsigIntensityVector[i] = sqrt(pow(sigIntensityVector[i],2.) + pow(IntensityMatrix[j][i] * IntensityVector[j],2.)*(pow(sigIntensityMatrix[j][i]/IntensityMatrix[j][i],2.) + pow(sigIntensityVector[j]/IntensityVector[j],2.)));
                }
                
                
            }
            //             std::cout << std::endl;
            IntensityVector[j] = 0;
        }
        IntensityVector = nextIntensityVector;
        sigIntensityVector = nextsigIntensityVector;
        count_steps++;
    }
    
    std::cout << "State chosen has energy " << vEx.at(InitialStateIndex) << " keV" << std::endl;
    std::cout << "Ground-state fraction: " << IntensityVector[0] << "+-" << sigIntensityVector[0] << std::endl;
    std::cout << "Meta-state fraction: " << IntensityVector[1] << "+-" << sigIntensityVector[1] << std::endl;
    
    std::cout << "Sum check = " << IntensityVector[0] + IntensityVector[1] << std::endl;
}

