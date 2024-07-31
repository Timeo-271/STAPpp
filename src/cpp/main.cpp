/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Bar.h"
#include "T3.h"
#include "Outputter.h"
#include "Clock.h"

using namespace std;

int main(int argc, char *argv[])
{  

	if (argc != 2)
	{  
	    cout << "Usage: stap++ InputFileName\n";   
		exit(1); 
	}  
  
	string filename(argv[1]); 
    size_t found = filename.find_last_of('.'); 
  
    // 如果文件名包含扩展名  
    if (found != std::string::npos) {
        if (filename.substr(found) == ".dat")
            filename = filename.substr(0, found); 
        else {  
            cout << "*** Error *** Invalid file extension: " 
                 << filename.substr(found+1) << endl;   
            exit(1); 
        }  
    }  
  
    string InFile = filename + ".dat"; 
	string OutFile = filename + ".out";  
  
	CDomain* FEMData = CDomain::GetInstance();  
  
    Clock timer; 
    timer.Start();  
    
	if (!FEMData->ReadData(InFile, OutFile))
	{  
		cerr << "*** Error *** Data input failed!" << endl;
		exit(1);
	}  
      
    double time_input = timer.ElapsedTime(); 
    COutputter* Output = COutputter::GetInstance();   
  
    if (!FEMData->GetMODEX())  
    {  
        *Output << "Data check completed !" << endl << endl;  
        return 0;
    }  
  
	FEMData->AllocateMatrices(); // allocate storage
    
	FEMData->AssembleStiffnessMatrix();
      
    double time_assemble = timer.ElapsedTime(); 
    
	CLDLTSolver* Solver = new CLDLTSolver(FEMData->GetStiffnessMatrix()); 
      
    Solver->LDLT(); // 

    #ifdef _DEBUG_ // if defined DEBUG 
        Output->PrintStiffnessMatrix(); 
    #endif
         
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)  
    {  
        FEMData->AssembleForce(lcase + 1);  
        
        Solver->BackSubstitution(FEMData->GetForce());  
    
        // output case
        *Output << " LOAD CASE" << setw(5) << lcase + 1 << endl << endl << endl;  
    
    #ifdef _DEBUG_   
        Output->PrintDisplacement();  
    #endif  
    
        // output dis
        Output->OutputNodalDisplacement();  
    
        // calculate stress info
        Output->OutputElementStress();  
    }  
    
    double time_solution = timer.ElapsedTime();  
    
    timer.Stop();  
    
    // output time
    *Output << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"  
            << "     TIME FOR INPUT PHASE = " << time_input << endl  
            << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl  
            << "     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS = " << time_solution - time_assemble << endl << endl  
            << "     T O T A L   S O L U T I O N   T I M E = " << time_solution << endl << endl;  
    
    return 0;  
}