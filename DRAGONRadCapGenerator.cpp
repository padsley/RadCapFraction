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
    LevelsHit[InitialStateIndex] = true;
    
    //     vector<int> ListOfLevelsHit;
    
//     double FinalLevels = 2;
    
    
    
    vector <int> ListOfLevelsHit = ProcessDecay(InitialStateIndex);
    while(ListOfLevelsHit.size()>LevelLimit)
    {
        std::cout << "Too many levels! Trying with the new minimum intensity of: " << MinIntensityAllowed << std::endl;
        ListOfLevelsHit = ProcessDecay(InitialStateIndex);        
    }
    
    if(ListOfLevelsHit.size()==0)std::cout << "Somehow have managed to get no levels ever populated..." << std::endl;
    
    for(unsigned int i=0;i<ListOfLevelsHit.size();i++)
        LevelsHit[i] = true;
            
        double SumIntensity = 0;
        for(int i=0;i<nStates;i++)
        {
            if(IntensityVector[i]>0 && VerboseFlag)std::cout << "Final population of level = " << IntensityVector[i] << std::endl;
            SumIntensity += IntensityVector[i];
        }
        
        std::cout << "Sum of intensities at the end = " << SumIntensity << std::endl;
        
        for(unsigned int i=0;i<ListOfLevelsHit.size();i++)
        {
            std::cout << "level(" << i << ") = " << vEx.at(ListOfLevelsHit.at(i))*1.e-3 << std::endl;
        }
        
        //now to loop over each of the levels hit and to work out the decay from each of the levels
        for(unsigned int i=0;i<ListOfLevelsHit.size();i++)
        {
            std::cout << std::endl;
            
            int DecayIndex = 1;
            
            int CountDecaysFromLevel = 0;
            
            bool TooManyDecays = false;
            
            vector<double> ListOfDecayIntensities;//This is here to store a list of the decay intensities
            ListOfDecayIntensities.clear();
            
            double WeakestDecayLimit = 0.0; //To be filled should there be a minimum decay intensity from the above step
            
            for(unsigned int j=0;j<ListOfLevelsHit.size();j++)
            {
                if(IntensityMatrix[ListOfLevelsHit.at(i)][ListOfLevelsHit.at(j)]>0)
                {
                    CountDecaysFromLevel++;
                    
                    ListOfDecayIntensities.push_back(IntensityMatrix[ListOfLevelsHit.at(i)][ListOfLevelsHit.at(j)]);
                }
            }
            
            if(CountDecaysFromLevel>DecayLimit)
            {
                std::cout << "Number of decays from this state exceeds the limit - trying to find the " << DecayLimit << " strongest decays to print instead" << std::endl;
                
                TooManyDecays = true;
                std::sort(ListOfDecayIntensities.begin(), ListOfDecayIntensities.end(), std::greater<double>());
                
                printVec(ListOfDecayIntensities);
                
                WeakestDecayLimit = ListOfDecayIntensities.at(DecayLimit);
                
                std::cout << "The weakest decay limit is: " << ListOfDecayIntensities.at(DecayLimit) << std::endl;
            }
            
            std::cout << "Have " << CountDecaysFromLevel << " decays from level " << i << std::endl;
            
            for(unsigned int j=0;j<ListOfLevelsHit.size();j++)
            {
                if(IntensityMatrix[ListOfLevelsHit[i]][ListOfLevelsHit[j]]>0 && !TooManyDecays)
                {
                    std::cout << "br(" << i << "," << DecayIndex << ") = " << IntensityMatrix[ListOfLevelsHit.at(i)][ListOfLevelsHit.at(j)]*100. << std::endl;
                    std::cout << "md(" << i << "," << DecayIndex << ") = " << j << std::endl;
                    
                    DecayIndex++;
                }
                else if(IntensityMatrix[ListOfLevelsHit[i]][ListOfLevelsHit[j]]>WeakestDecayLimit && TooManyDecays)
                {
                    
                    
                    std::cout << "br(" << i << "," << DecayIndex << ") = " << IntensityMatrix[ListOfLevelsHit.at(i)][ListOfLevelsHit.at(j)]*100. << std::endl;
                    std::cout << "md(" << i << "," << DecayIndex << ") = " << j << std::endl;
                    
                    DecayIndex++;
                }
            }
        }
        
        
        
        /*
         *       for(int j=2;j<nStates;j++)//this is the initial state
         *       {
         *           for(int i=0;i<nStates;i++)//this is the final state
         *           {
         *               if(IntensityMatrix[j][i]>0 && IntensityVector[j]>0)
         *               {                    
         *                   nextIntensityVector[i] += IntensityMatrix[j][i] * IntensityVector[j];
         *                   nextsigIntensityVector[i] = sqrt(pow(sigIntensityVector[i],2.) + pow(IntensityMatrix[j][i] * IntensityVector[j],2.)*(pow(sigIntensityMatrix[j][i]/IntensityMatrix[j][i],2.) + pow(sigIntensityVector[j]/IntensityVector[j],2.)));
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


std::vector<int> ProcessDecay(int InitialStateIndex)
{
    //         std::cout << "Step " << count_steps << endl;
    double TempMinIntensity = 1.; 
    vector<int> ListOfLevelsHit;
    ListOfLevelsHit.clear();
    
    for(int i=0;i<nStates;i++){
        IntensityVector[i] = 0;
        nextIntensityVector[i] = 0;
        sigIntensityVector[i] = 0;
        nextsigIntensityVector[i] = 0;
    }
    
    IntensityVector[InitialStateIndex] = 1.0;
    
    
    
    while(!CheckFinished(IntensityVector,FinalLevels))
    {
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
        for(int j=FinalLevels;j<nStates;j++)//loop over the initial states
        {
            for(int i=0;i<nStates;i++)//loop over the possible final states
            {
                if(IntensityMatrix[j][i]>MinIntensityAllowed && IntensityVector[j]>0) //check to see that there's something in the initial state and that there's a transition into the final state
                {
                    if(IntensityMatrix[j][i]<TempMinIntensity)TempMinIntensity=IntensityMatrix[j][i];
                    if(VerboseFlag)std::cout << "New minimum intensity = " << TempMinIntensity << std::endl;
                    
                    if(VerboseFlag)std::cout << "i = " << i << "\t" << "j = " << j << "\t IntensityMatrix = " << IntensityMatrix[j][i] << std::endl; 
                    
                    
                    
                    ListOfLevelsHit.push_back(i);//This makes a vector listing the levels hit in my calculation
                    
                    nextIntensityVector[i] += IntensityMatrix[j][i] * IntensityVector[j];
                    nextsigIntensityVector[i] = sqrt(pow(sigIntensityVector[i],2.) + pow(IntensityMatrix[j][i] * IntensityVector[j],2.)*(pow(sigIntensityMatrix[j][i]/IntensityMatrix[j][i],2.) + pow(sigIntensityVector[j]/IntensityVector[j],2.)));
                }
            }
            IntensityVector[j] = 0;
            
        }
        IntensityVector = nextIntensityVector;
        sigIntensityVector = nextsigIntensityVector;
        
        //         count_steps++;
    }//end of the "finished" while condition
    
    MinIntensityAllowed = TempMinIntensity;
    
    return ListOfLevelsHit;
}
