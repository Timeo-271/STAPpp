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
#include "Material.h"
#include "ElementGroup.h"
#include "gauss.h"
#include "T3.h"
#include "Q4.h"
#include <math.h>

using namespace std;

//	Clear an array
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()  
{  
    Title[0] = '0';  
    MODEX = 0;   
  
    NUMNP = 0;      
    NodeList = nullptr;
      
    NUMEG = 0;    
    EleGrpList = nullptr; 
      
    NLCASE = 0;     
    NLOAD = nullptr; 
    LoadCases = nullptr;  
      
    NEQ = 0;      
  
    Force = nullptr; 
    StiffnessMatrix = nullptr;  
}

//	Desconstructor
CDomain::~CDomain()  
{  
    delete [] NodeList;  
  
    delete [] EleGrpList;  
  
    delete [] NLOAD;
    delete [] LoadCases; 
  
    delete [] Force;
    delete StiffnessMatrix;

}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::GetInstance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::GetInstance(OutFile);

//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();

//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data
	if (ReadNodalPoints()){
        Output->OutputNodeInfo();}
    else{
        cerr << "*** Error *** Fail to read nodal point !" << endl;
        return false;
        }
    
//	Read load data
	if (ReadLoadCases())
        Output->OutputLoadInfo();
    else{
        cerr << "*** Error *** Fail to read load cases!" << endl;
        return false;
        }

//	Read element data
	if (ReadElements())
        Output->OutputElementInfo();
    else{
        cerr << "*** Error *** Fail to read element information !" << endl;
        return false;
        }

//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

	return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()  
{   
//	Read nodal point data lines  
    NodeList = new CNode[NUMNP];  
    //std::cout<<NUMNP;
//	Loop over for all nodal points  
    for (unsigned int np = 0; np < NUMNP; np++)  
    { 
        if (!NodeList[np].Read(Input))  
        {
            cerr << "*** Error *** Fail to read node list !" << endl;
            return false;  
        }  
  
        if (NodeList[np].NodeNumber != np + 1)  
        {
            cerr << "*** Error *** all node must be input in order !" << endl  
            << "   expect node : " << np + 1 << endl  
            << "   provided node : " << NodeList[np].NodeNumber << endl;  
  
            return false;  
        }  
    }  
    return true;  
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{

    NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < NodeList[np].NDF_; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode_[dof]) 
				NodeList[np].bcode_[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode_[dof] = NEQ;
			}
		}
	}
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	
    LoadCases = new CLoadCaseData[NLCASE];	// List all load cases
    
//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
    {
        unsigned int LL;
        
        Input >> LL;
        
        if (LL != lcase + 1)
        {
            cerr << "*** Error *** Load case must be inputted in order !" << endl
            << "   Expected load case : " << lcase + 1 << endl
            << "   Provided load case : " << LL << endl;
            
            return false;
        }

        LoadCases[lcase].Read(Input);
    }

	return true;
}

// Read element data 
bool CDomain::ReadElements()  
{  
    EleGrpList = new CElementGroup[NUMEG];  
  
    for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)  
    {  
        if (!EleGrpList[EleGrp].Read(Input))  
            return false;  
    }
    return true;  
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()  
{  
#ifdef _DEBUG_  
    COutputter* Output = COutputter::GetInstance();  
    *Output << setw(9) << "Ele = " << setw(22) << "Location Matrix" << endl;  
#endif  
  
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups  
    {  
        CElementGroup& ElementGrp = EleGrpList[EleGrp];  
        unsigned int NUME = ElementGrp.GetNUME();  
          
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp  
        {  
            CElement& Element = ElementGrp[Ele];  
  
            // Generate location matrix  
            Element.GenerateLocationMatrix();  
              
#ifdef _DEBUG_  
            unsigned int* LocationMatrix = Element.GetLocationMatrix();  
              
            *Output << setw(9) << Ele+1;  
            for (int i=0; i<Element.GetND(); i++)  
                *Output << setw(5) << LocationMatrix[i];  
            *Output << endl;  
#endif  
  
            StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());  
        }  
    }  
      
    StiffnessMatrix->CalculateMaximumHalfBandwidth();  
       
#ifdef _DEBUG_  
    *Output << endl;  
	Output->PrintColumnHeights();  
#endif  
  
}

//    Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//    and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
    //    Allocate for global force/displacement vector
    Force = new double[NEQ];
    
    //  Create the banded stiffness matrix
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);
    
    //    Calculate column heights
    CalculateColumnHeights();
    
    //    Calculate address of diagonal elements in banded matrix
    StiffnessMatrix->CalculateDiagnoalAddress();
    
    //    Allocate for banded global stiffness matrix
    StiffnessMatrix->Allocate();
    
    COutputter* Output = COutputter::GetInstance();
    Output->OutputTotalSystemData();
}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()  
{  
    // for all element 
    
    for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)  
    {  
        
        CElementGroup& ElementGrp = EleGrpList[EleGrp];  
           
        unsigned int NUME = ElementGrp.GetNUME();  
           
        unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
        //std::cout<<size<<endl;
        double* Matrix = new double[size]; 
          
        for (unsigned int Ele = 0; Ele < NUME; Ele++)  
        {  
            CElement& Element = ElementGrp[Ele];  
            
            Element.ElementStiffness(Matrix); 
            StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND()); 
        }  

        //std::cout<<Matrix[0]<<" "<<Matrix[1]<<" "<<Matrix[2]<<" "<<Matrix[3]<<" "<<Matrix[4]
        //<<" "<<Matrix[5]<<" "<<Matrix[6]<<" "<<Matrix[7]<<" "<<Matrix[8]
        //<<endl;
           
        delete[] Matrix;  
        Matrix = nullptr;  
    }  
       
#ifdef _DEBUG_  
    COutputter* Output = COutputter::GetInstance();  
    Output->PrintStiffnessMatrix();  
#endif  
}

//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)  
{  
	if (LoadCase > NLCASE)   
		return false;  
	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];  
   
    clear(Force, NEQ);
    for(unsigned int grpnum = 0; grpnum<NUMEG; grpnum++){
        
        CElementGroup& ElementGrp = EleGrpList[grpnum];
        ElementTypes ElementType = ElementGrp.GetElementType(); //get Element Type
        unsigned int NUME = ElementGrp.GetNUME();//get NUME
        switch (ElementType)  
        {  
            case ElementTypes::Bar://bar:only force
                for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)  
                {  
                    // get dof 
                    unsigned int dof = NodeList[LoadData->node_load[lnum] - 1].bcode_[LoadData->dof[lnum] - 1];  
                    
                    if(dof) // The DOF is activated  
                        Force[dof - 1] += LoadData->load[lnum];  
                } 
            break;

            case ElementTypes::T3: 
            //force
                for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)  
                {   
                    unsigned int dof = NodeList[LoadData->node_load[lnum] - 1].bcode_[LoadData->dof[lnum] - 1];  
                    
                    if(dof) // The DOF is activated  
                        Force[dof - 1] += LoadData->load[lnum];  
                }
            //naturla bc
                if(!(LoadData->EleGrp_num - grpnum - 1)){
                    for (unsigned int lnum = 0; lnum < LoadData->nnbc; lnum++) //2 points at a time 
                    {  
                        unsigned int dof1 = NodeList[LoadData->node_nbc[2*lnum] - 1].bcode_[LoadData->dof_nbc[lnum] - 1];  //left
                        unsigned int dof2 = NodeList[LoadData->node_nbc[2*lnum+1] - 1].bcode_[LoadData->dof_nbc[lnum] - 1]; //right
                        //length
                        double DX[2];
                        DX[0] = NodeList[LoadData->node_nbc[2*lnum] - 1].XYZ[0] - NodeList[LoadData->node_nbc[2*lnum+1] - 1].XYZ[0];
                        DX[1] = NodeList[LoadData->node_nbc[2*lnum] - 1].XYZ[1] - NodeList[LoadData->node_nbc[2*lnum+1] - 1].XYZ[1];
                        double DX2[2];	//  Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
                        DX2[0] = DX[0] * DX[0];
                        DX2[1] = DX[1] * DX[1];
                        double L2 = DX2[0] + DX2[1];
                        double L = sqrt(L2); //length
                    
                        CMaterial& Mat = ElementGrp.GetMaterial(grpnum);
                        C2DMaterial& CT3mat = static_cast<C2DMaterial&>(Mat);
                        double h = CT3mat.thickness;
                        //ElementGrp.GetMaterial
                        
                        if(dof1){ // The DOF is activated  
                            double t1 = LoadData->nbc[2*lnum];
                            double t2 = LoadData->nbc[2*lnum+1];
                            Force[dof1 - 1] += h*L*(2*t1+t2)/6; 
                        }
                        if(dof2){ // The DOF is activated  
                            double t1 = LoadData->nbc[2*lnum];
                            double t2 = LoadData->nbc[2*lnum+1]; 
                            Force[dof2 - 1] += h*L*(2*t2+t1)/6;
                        }
                        
                    } 
                }
            //body force
                for (unsigned int Ele = 0; Ele < NUME; Ele++)  
                {   
                    CElement& Element = ElementGrp[Ele]; //CElement
                    CT3& T3Element = dynamic_cast<CT3&>(Element); //turn to CT3
                    //std::cout<<size<<endl;
                    
                    unsigned int size = 9;
                    double* F = new double[size];
                    T3Element.ElementForce(F); //for every element
                    //std::cout<<F[0]<<F[1]<<endl;
                    CNode** ElementNode = Element.GetNodes(); //3 points in element
                    CNode* Node0 = ElementNode[0];//3 nodes
                    CNode* Node1 = ElementNode[1];
                    CNode* Node2 = ElementNode[2];
                    unsigned int dof[9] ={0};
                    dof[0] = NodeList[Node0->NodeNumber - 1].bcode_[0]; //1x
                    dof[1] = NodeList[Node0->NodeNumber - 1].bcode_[1]; 
                    dof[3] = NodeList[Node1->NodeNumber - 1].bcode_[0]; 
                    dof[4] = NodeList[Node1->NodeNumber - 1].bcode_[1]; 
                    dof[6] = NodeList[Node2->NodeNumber - 1].bcode_[0];
                    dof[7] = NodeList[Node2->NodeNumber - 1].bcode_[1];

                    for( int i = 0; i<9 ; i++){
                        if(dof[i]) 
                            Force[dof[i] - 1] += F[i]; 
                    }
                    delete[] F;
                    F = nullptr;
                }
            break;

            case ElementTypes::Q4: 
            //force
                for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)  
                {   
                    unsigned int dof = NodeList[LoadData->node_load[lnum] - 1].bcode_[LoadData->dof[lnum] - 1];  
                    
                    if(dof) // The DOF is activated  
                        Force[dof - 1] += LoadData->load[lnum];  
                }
            //naturla bc
                if(!(LoadData->EleGrp_num - grpnum - 1)){
                    for (unsigned int lnum = 0; lnum < LoadData->nnbc; lnum++) //2 points at a time 
                    {  
                        unsigned int dof1 = NodeList[LoadData->node_nbc[2*lnum] - 1].bcode_[LoadData->dof_nbc[lnum] - 1];  //left
                        unsigned int dof2 = NodeList[LoadData->node_nbc[2*lnum+1] - 1].bcode_[LoadData->dof_nbc[lnum] - 1]; //right
                        //length
                        double DX[2];
                        DX[0] = NodeList[LoadData->node_nbc[2*lnum] - 1].XYZ[0] - NodeList[LoadData->node_nbc[2*lnum+1] - 1].XYZ[0];
                        DX[1] = NodeList[LoadData->node_nbc[2*lnum] - 1].XYZ[1] - NodeList[LoadData->node_nbc[2*lnum+1] - 1].XYZ[1];
                        double DX2[2];	//  Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
                        DX2[0] = DX[0] * DX[0];
                        DX2[1] = DX[1] * DX[1];
                        double L2 = DX2[0] + DX2[1];
                        double L = sqrt(L2); //length
                        double J = L/2;
                        CMaterial& Mat = ElementGrp.GetMaterial(grpnum);
                        C2DMaterial& CQ4mat = static_cast<C2DMaterial&>(Mat);
                        unsigned int Elenum = LoadData->Ele_num;

                        CElement& Element = ElementGrp[Elenum-1]; //CElement
                        CQ4& Q4Element = dynamic_cast<CQ4&>(Element); //turn to CQ4

                        double* Fbc = new double[2];
                        
                        double t1 = LoadData->nbc[2*lnum];
                        double t2 = LoadData->nbc[2*lnum+1];

                        Q4Element.Calculate_NBC(Fbc, t1, t2);
                        //ElementGrp.GetMaterial
                        

                        if(dof1){ // The left DOF is activated  
                            Force[dof1 - 1] += J*Fbc[0]; 
                        }
                        if(dof2){ // The right DOF is activated  
                            Force[dof2 - 1] += J*Fbc[1];
                        }
                        delete[] Fbc;
                        Fbc = nullptr;
                        //std::cout<<Force[0]<<Force[1]<<Force[2]<<Force[3]<<endl;
                    } 
                }
            //body force
                for (unsigned int Ele = 0; Ele < NUME; Ele++)  
                {   
                    CElement& Element = ElementGrp[Ele]; //CElement
                    CQ4& Q4Element = dynamic_cast<CQ4&>(Element); //turn to CQ4
                    //std::cout<<size<<endl;
                    CMaterial& Mat = ElementGrp.GetMaterial(grpnum);
                    C2DMaterial& CQ4mat = static_cast<C2DMaterial&>(Mat);

                    unsigned int size = 12;
                    double* F = new double[size];

                    Q4Element.Calculate_BODY(F);

                    //std::cout<<F[0]<<F[1]<<endl;
                    CNode** ElementNode = Element.GetNodes(); //4 points in element
                    CNode* Node0 = ElementNode[0];//3 nodes
                    CNode* Node1 = ElementNode[1];
                    CNode* Node2 = ElementNode[2];
                    CNode* Node3 = ElementNode[3];
                    unsigned int dof[12] ={0};

                    dof[0] = NodeList[Node0->NodeNumber - 1].bcode_[0]; //1x
                    dof[1] = NodeList[Node0->NodeNumber - 1].bcode_[1]; 
                    dof[3] = NodeList[Node1->NodeNumber - 1].bcode_[0]; 
                    dof[4] = NodeList[Node1->NodeNumber - 1].bcode_[1]; 
                    dof[6] = NodeList[Node2->NodeNumber - 1].bcode_[0];
                    dof[7] = NodeList[Node2->NodeNumber - 1].bcode_[1];
                    dof[9] = NodeList[Node3->NodeNumber - 1].bcode_[0];
                    dof[10] = NodeList[Node3->NodeNumber - 1].bcode_[1];

                    for( int i = 0; i<12 ; i++){
                        if(dof[i]) 
                            Force[dof[i] - 1] += F[i]; 
                    }
                    delete[] F;
                    F = nullptr;
                }
                
                break;

            default: //error
                    std::cerr << "NDF Do not fit CNode::Read." << std::endl;
                    exit(5);
                    break;
        }
    }
	return true;  
}

