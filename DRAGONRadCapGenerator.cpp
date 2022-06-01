#include "RadCapFunctions.h"

void DRAGONRadCapGenerator(int InitialStateIndex=2)
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
    
//     PrintLevels();
    
    int count_steps = 0;
    
    vector<double> ExcitationsOfLevelsIncluded;
    bool *LevelsHit = new bool[nStates];
    for(int i=0;i<nStates;i++)LevelsHit[i] = false;
    
    vector<int> ListOfLevelsHit;
    
    while(!CheckFinished(IntensityVector,2) && count_steps<3)
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

        
        //Find which levels have been touched by the decay cascade
        for(int j=0;j<nStates;j++)//loop over the initial states
        {
            for(int i=0;i<nStates;i++)//loop over the possible final states
            {
                if(IntensityMatrix[j][i]>0 && IntensityVector[j]>0) //check to see that there's something in the initial state and that there's a transition into the final state
                {
                    LevelsHit[i] = true;
                    
                    nextIntensityVector[i] += IntensityMatrix[j][i] * IntensityVector[j];
                    nextsigIntensityVector[i] = sqrt(pow(sigIntensityVector[i],2.) + pow(IntensityMatrix[j][i] * IntensityVector[j],2.)*(pow(sigIntensityMatrix[j][i]/IntensityMatrix[j][i],2.) + pow(sigIntensityVector[j]/IntensityVector[j],2.)));
                }
            }
            IntensityVector[j] = 0;
            
        }
        IntensityVector = nextIntensityVector;
        sigIntensityVector = nextsigIntensityVector;
        
//         count_steps++;
    }//end of the while condition
    
    for(int i=0;i<nStates;i++)
        if(LevelsHit[i])ListOfLevelsHit.push_back(i);//This makes a vector listing the levels hit in my calculation
        
        
    for(unsigned int i=0;i<ListOfLevelsHit.size();i++)
    {
        std::cout << "level(" << i << ") = " << vEx.at(ListOfLevelsHit.at(i)) << std::endl;
    }
        
        /*
        for(int j=2;j<nStates;j++)//this is the initial state
        {
            for(int i=0;i<nStates;i++)//this is the final state
            {
                if(IntensityMatrix[j][i]>0 && IntensityVector[j]>0)
                {                    
                    nextIntensityVector[i] += IntensityMatrix[j][i] * IntensityVector[j];
                    nextsigIntensityVector[i] = sqrt(pow(sigIntensityVector[i],2.) + pow(IntensityMatrix[j][i] * IntensityVector[j],2.)*(pow(sigIntensityMatrix[j][i]/IntensityMatrix[j][i],2.) + pow(sigIntensityVector[j]/IntensityVector[j],2.)));
                }
                
                
            }
            //             std::cout << std::endl;
            IntensityVector[j] = 0;
        }
        IntensityVector = nextIntensityVector;
        sigIntensityVector = nextsigIntensityVector;
        count_steps++;
        */
    
}
