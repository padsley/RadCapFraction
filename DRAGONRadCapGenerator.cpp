#include "RadCapFunctions.h"



void DRAGONRadCapGenerator()
{
    
    //Standard function for loading variables and declaring stuff
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
    
    PrintLevels();
}
