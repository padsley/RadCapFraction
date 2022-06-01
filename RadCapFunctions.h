#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <TString.h>

//**************************** THESE ARE THE THINGS TO CHANGE ***********************************
string ENSDF_input = "Ca40_ensdf.dat"; //change this file to get a different input - for some reason that I don't quite understand you need to have one blank line at the end of this file to get it not to crash horribly
int MassNumber = 40; //change this for the right mass number for the decay
string Symbol = "CA"; //change this for the right symbol for the decay NOTE THAT IT'S ALL CAPS!!!!

//40Ca
double MassRecoil = 37225.217616*1000; //This is the recoil mass in units of keV/c/c

//26Al:
// string ENSDF_input = "Al26_ensdf.dat";
// int MassNumber = 26;
// string Symbol = "AL";
// double MassRecoil = 24206.831449*1000;
//************************************************************************************************


vector<double> vEx;

int nStates = 0;

bool VerboseFlag = true;
void PrintIntensityVector(double *IntensityVector, int InitialStateIndex);
void PrintIntensityMatrix(double **IntensityMatrix);

bool CheckFinished(double *IntensityVector, int NumberFinalStates);

void LoadLevelInformation();
void LoadDecayInformation();
void NormaliseDecayInformation(double **IntensityMatrix);

double **IntensityMatrix;
double **sigIntensityMatrix;
double *IntensityVector;
double *nextIntensityVector;
double *sigIntensityVector;
double *nextsigIntensityVector;



double CorrectForRecoil(double Egamma, double mass_recoil)
{
    double result = 0;
    
//     double m_26Al = 24206.831449*1000; //keV/c/c
    
    double RecoilEnergy = 0.5 * pow(Egamma,2.) / mass_recoil;
    
    return Egamma + RecoilEnergy;
    
}

void PrintIntensityVector(double *IntensityVector,int InitialStateIndex)
{
    for(int i=0;i<nStates && i<=InitialStateIndex;i++)
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

bool CheckFinished(double *IntensityVector, int NumberFinalStates)//NumberFinalStates is the list of how many final states should be regarded to "finished" - for 26Al should be 2 since can be in the g.s. or the isomeric state, for other cases should be 1 for the ground state which should collect all of the strength
{
    bool result = true;
    
    for(int i=NumberFinalStates;i<nStates;i++){
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
    input.open(ENSDF_input);
    
    if(input.is_open())
    {
        for(std::string line; std::getline(input,line);)
        {
            std::string testStatement;
            testStatement = std::to_string(MassNumber);
            testStatement.append(Symbol);
            testStatement.append("  L");
            
            if(strcmp(line.substr(1,7).c_str(),testStatement.c_str())==0)
            {
                vEx.push_back(atof(line.substr(9,13).c_str()));
            }
        }
    }
    
    nStates = vEx.size();
    std::cout << "nStates = " << nStates << std::endl;
        
    input.close();
    
    std::cout << "Finished loading level information" << std::endl;
}

void PrintLevels()
{
    std::cout << "Printing Level Information" << std::endl;
    if(vEx.size()==0)LoadLevelInformation();
    
    for(unsigned int i=0;i<vEx.size();i++)
        std::cout << i << "\t" << vEx.at(i) << std::endl;
}

void LoadDecayInformation()
{
    std::cout << "Load Decay Information" << std::endl;
    
    fstream input;
    input.open(ENSDF_input);
    
    if(input.is_open())
    {
        std::vector<double>::iterator it;
        int StateIndex = 0; //state index is the index of the decaying state
        for(std::string line; std::getline(input,line);)
        {
            std::string testStatement;//FIND THE LEVELS
            testStatement = std::to_string(MassNumber);
            testStatement.append(Symbol);
            testStatement.append("  L");
            
            std::string testStatement2;//FIND THE GAMMAS
            testStatement2 = std::to_string(MassNumber);
            testStatement2.append(Symbol);
            testStatement2.append("  G");
            
            if(strcmp(line.substr(1,7).c_str(),testStatement.c_str())==0)
            {                
                it = find(vEx.begin(), vEx.end(), (double)atof(line.substr(9,13).c_str()));
                
                if(it!=vEx.end())
                {
                    StateIndex = std::distance(vEx.begin(),it);
                }
                else
                {
                    std::cout << "Energy not found in vector" << std::endl;
                }
            }
            else if(strcmp(line.substr(1,7).c_str(),testStatement2.c_str())==0)
            {             
//                 if(VerboseFlag)std::cout << "Finding gamma-ray information for state " << std::endl;
                
                double Egamma = atof(line.substr(9,10).c_str());
                
                Egamma = CorrectForRecoil(Egamma,MassRecoil);
                
                double FinalEnergy = *it - Egamma;
                
                double Intensity = atof(line.substr(22,7).c_str());
                
                //work out the number of digits for the Intensity value
                std::string counter(line.substr(22,7));
                counter.erase(std::remove(counter.begin(), counter.end(), ' '), counter.end());
                int length = counter.length();                
                
                int pos_decimal = 0, before_decimal = -1, after_decimal = -1;
                
                if(counter.find('.',0)!=string::npos)
                {
                    counter.find('.',0); //find the location of the decimal point
                
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
                
                for(unsigned int i=0;i<vEx.size();i++)//looking for the final level so looping over the list here
                {
//                     if(VerboseFlag)std::cout << abs(FinalEnergy - vEx.at(i)) << std::endl;
                    if(abs(FinalEnergy - vEx.at(i))<0.5)//only care about if the gamma energy matches a known level - if this isn't true ever then it's a problem I guess?
                    {
                        std::cout << "Found matching state with index " << StateIndex << std::endl;
                        IntensityMatrix[StateIndex][i] = Intensity;
                        sigIntensityMatrix[StateIndex][i] = sigIntensity;
                    }
                    else
                    {
                        //std::cout << "Did NOT find matching state" << std::endl;
                        //I don't think that this is actually useful because it's just saying "this isn't the right decay"
                    }
                }
            }
        }
    }
}