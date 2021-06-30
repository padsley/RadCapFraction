#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

vector<double> vEx;

int nStates = 0;

bool VerboseFlag = false;
void PrintIntensityVector(double *IntensityVector);
void PrintIntensityMatrix(double **IntensityMatrix);

bool CheckFinished(double *IntensityVector);

void LoadLevelInformation();
void LoadDecayInformation();
void NormaliseDecayInformation(double **IntensityMatrix);

double **IntensityMatrix;
double **sigIntensityMatrix;
double *IntensityVector;
double *nextIntensityVector;
double *sigIntensityVector;
double *nextsigIntensityVector;

double CorrectForRecoil(double Egamma);

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
    while(!CheckFinished(IntensityVector))
    {
//         std::cout << "Step " << count_steps << endl;
        
                        if(VerboseFlag)
        {
            std::cout << "IntensityVector" << std::endl;
            PrintIntensityVector(IntensityVector);
        }
                        if(VerboseFlag)
        {
            std::cout << "nextIntensityVector" << std::endl;
            PrintIntensityVector(nextIntensityVector);
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

void PrintIntensityVector(double *IntensityVector)
{
    for(int i=0;i<nStates;i++)
        std::cout << "IntensityVector[" << i << "] = " << IntensityVector[i] << std::endl;
    
    std::cout << std::endl;
}

void PrintIntensityMatrix(double **IntensityMatrix)
{
    for(int i=0;i<nStates;i++)
    {
        std::cout << "i = " << i << "\t";
        for(int j=0;j<nStates;j++)
        {
            std::cout << IntensityMatrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

bool CheckFinished(double *IntensityVector)
{
    bool result = true;
    
    for(int i=2;i<nStates;i++){
        if(IntensityVector[i]>0)result = false;
    }
    
    return result;
}

void NormaliseDecayInformation(double **IntensityMatrix)
{
    for(int i=0;i<nStates;i++)
    {
        double sum = 0., sigma_sum = 0;
        for(int j=0;j<nStates;j++)
        {
            if(IntensityMatrix[i][j]>0)
            {
                sum += IntensityMatrix[i][j];
                sigma_sum = sqrt(pow(sigma_sum,2.) + pow(sigIntensityMatrix[i][j],2.));
            }
        }
        
        for(int j=0;j<nStates;j++)
        {
            if(sum>0)
            {
                double dummy = IntensityMatrix[i][j];
                IntensityMatrix[i][j] /= sum;
                
//                 sigIntensityMatrix[i][j] /= sum;
//                 sigIntensityMatrix[i][j] = sqrt(pow(sigIntensityMatrix[i][j],2.) +
                sigIntensityMatrix[i][j] = IntensityMatrix[i][j] * sqrt(pow(sigIntensityMatrix[i][j]/dummy,2.) + pow(sigma_sum/sum,2.));
            }
                
            else if (sum==0)
            {
                IntensityMatrix[i][j] = 0;
                sigIntensityMatrix[i][j] = 0;
            }
            
            
//             if(IntensityMatrix[i][j]>0)std::cout << IntensityMatrix[i][j] << "\t" << sigIntensityMatrix[i][j] << "\t" << 100.*sigIntensityMatrix[i][j]/IntensityMatrix[i][j] << std::endl;
        }
    }
}

void LoadLevelInformation()
{
    ifstream input;
    input.open("Al26_ensdf.dat");
    
    if(input.is_open())
    {
        for(std::string line; std::getline(input,line);)
        {
            //          std::cout << line.substr(1,7) << std::endl;
            if(strcmp(line.substr(1,7).c_str(),"26AL  L")==0)
            {
                //              if(VerboseFlag)std::cout << "atof(line.substr(9,13).c_str() = " << atof(line.substr(9,13).c_str()) << std::endl;
                vEx.push_back(atof(line.substr(9,13).c_str()));
            }
        }
    }
    
    nStates = vEx.size();
    std::cout << "nStates = " << nStates << std::endl;
        
        input.close();
}

void PrintLevels()
{
    if(vEx.size()==0)LoadLevelInformation();
    
    for(unsigned int i=0;i<vEx.size();i++)
        std::cout << i << "\t" << vEx.at(i) << std::endl;
}

void LoadDecayInformation()
{
    fstream input;
    input.open("Al26_ensdf.dat");
    
    if(input.is_open())
    {
        std::vector<double>::iterator it;
        int StateIndex = 0; //state index is the index of the decaying state
        for(std::string line; std::getline(input,line);)
        {
            
            if(strcmp(line.substr(1,7).c_str(),"26AL  L")==0)
            {
                //              std::cout << "Current Energy is " << atof(line.substr(9,13).c_str()) << std::endl;
                
                it = find(vEx.begin(), vEx.end(), (double)atof(line.substr(9,13).c_str()));
                
                if(it!=vEx.end())
                {
                    StateIndex = std::distance(vEx.begin(),it);
                }
                else
                {
                    //                  std::cout << "Energy not found in vector" << std::endl;
                }
            }
            else if(strcmp(line.substr(1,7).c_str(),"26AL  G")==0)
            {             
                //              std::cout << "Gamma-ray transition found with energy " << atof(line.substr(9,13).c_str()) << std::endl;
                
                double Egamma = atof(line.substr(9,10).c_str());
                
                Egamma = CorrectForRecoil(Egamma);
                
                //              std::cout << "Start state energy: " << *it << std::endl;
                
                double FinalEnergy = *it - Egamma;
                //              std::cout << "FinalEnergy = " << FinalEnergy << std::endl;
                
                double Intensity = atof(line.substr(22,7).c_str());
                
                //work out the number of digits for the Intensity value
                std::string counter(line.substr(22,7));
                //tmp.erase(std::remove(tmp.begin(), tmp.end(), ' '), tmp.end());
                counter.erase(std::remove(counter.begin(), counter.end(), ' '), counter.end());
                int length = counter.length();
//                 std::cout << "Length = " << length << std::endl;
                
                
                int pos_decimal = 0, before_decimal = -1, after_decimal = -1;
                
                if(counter.find('.',0)!=string::npos)
                {
                
                    counter.find('.',0); //find the location of the decimal point
                
//                     std::cout << "position of decimal = " << pos_decimal << std::endl;
                
                    before_decimal = pos_decimal - 1;
        
                    if(before_decimal>3)std::cout << "I don't know how but there are more than three values before the decimal here and it should be normalised to 100.0! The string was: " << line.substr(22,7).c_str() << std::endl;
                
                    after_decimal = length - pos_decimal - 1;
                }
                else
                {
                    pos_decimal = 0;
                    before_decimal = length;
                    after_decimal = 0;
                }
                
                double sigIntensity = atof(line.substr(29,3).c_str()); //need to convert this to be in the last digits of the previous result
                
                std::string sig_counter(line.substr(29,3));
                sig_counter.erase(std::remove(sig_counter.begin(), sig_counter.end(), ' '), sig_counter.end());
                
                double multiplier = 0;
                
                if(after_decimal>0)
                {
                    multiplier = pow(10.,1-1.*after_decimal);
                }
                else
                {
//                     std::cout << "No points after the decimal" << std::endl;
                }
                
                sigIntensity *= multiplier; //this should take the sigma and make it the right value for the 
                
//                 std::cout << "Intensity string = " << counter.c_str() << "\t+-\t" << sig_counter.c_str() << std::endl;
//                 std::cout << "Intensity = " << Intensity << "\t+-\t" << sigIntensity << std::endl << std::endl;
                
                for(unsigned int i=0;i<vEx.size();i++)//looking for the final level so looping over the list here
                {
                    if(abs(FinalEnergy - vEx.at(i))<0.5)//only care about if the gamma energy matches a known level - if this isn't true ever then it's a problem I guess?
                    {
                        //                      std::cout << "Found matching state with index " << StateIndex << std::endl;
                        IntensityMatrix[StateIndex][i] = Intensity;
                        sigIntensityMatrix[StateIndex][i] = sigIntensity;
                    }
                    else
                    {
                        //                      std::cout << "Did NOT find matching state" << std::endl;
                    }
                    
                    
                }
                
            }
        }
    }
}

double CorrectForRecoil(double Egamma)
{
    double result = 0;
    
    double m_26Al = 24206.831449*1000; //keV/c/c
    
    double RecoilEnergy = 0.5 * pow(Egamma,2.) / m_26Al;
    
    return Egamma + RecoilEnergy;
    
}
