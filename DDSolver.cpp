/**
* @file DDSolver.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-8-29   19:14
* @note
* @todo
*/

#include "DDSolver.h"
#include "FDDomain.h"
#include "Material.h"
#include "SctmPhys.h"
using namespace SctmPhys;
using namespace SctmUtils;

DriftDiffusionSolver::DriftDiffusionSolver(FDDomain *domain): totalVertices(domain->GetVertices())
{
	getDDVertices(domain);
	//for testing the DD solver, the preparation is done in the DDTest class.
	//prepareSolver();
}

void DriftDiffusionSolver::SolveDD()
{
	
}

void DriftDiffusionSolver::prepareSolver()
{
	SctmUtils::UtilsMsg.PrintHeader("Solving Drift-Diffusion");
	
	int vertSize = this->vertices.size();
	this->rhsVector.resize(vertSize);
	this->elecDensity.resize(vertSize);

	this->temperature = 300; //temporarily used here TODO: modify to accord with the whole simulation
	buildVertexMap();
	setTimeStep();
	buildCoefficientMatrix();
	buildRhsVector();
	//the structure does not change, so the final coefficient matrix after refreshing does not change either
	refreshCoefficientMatrix();
}

void DriftDiffusionSolver::buildVertexMap()
{
	std::pair<VertexMapInt::iterator, bool> insertPairVertex;
	std::pair<VertexMapDouble::iterator, bool> insertPairPrpty;
	int vertID = 0;
	int equationID = 0;
	double phyValue = 0;
	FDVertex *currVert = NULL;

	//this map is filled in order to obtain the vertex index from the vertex internal id. This is useful
	//in setting up the equation, i.e. filling the matrix.
	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		vertID = currVert->GetID();

		equationID = iVert; //equation index is also the vertex index in the vertices container
		insertPairVertex = this->equationMap.insert(VertexMapInt::value_type(vertID, equationID));
		SCTM_ASSERT(insertPairVertex.second==true, 10011);

		insertPairPrpty = this->mobilityMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys.GetPhysPrpty(PhysProperty::eMobility)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->potentialMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys.GetPhysPrpty(PhysProperty::ElectrostaticPotential)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->lastElecDensMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys.GetPhysPrpty(PhysProperty::eDensity)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);
	}
}

void DriftDiffusionSolver::setBndCondCurrent(vector<double> &in_current)
{
}

void DriftDiffusionSolver::buildCoefficientMatrix()
{
	int matrixSize = this->vertices.size();
	matrixSolver.matrix.resize(matrixSize, matrixSize);
	matrixSolver.matrix.reserve(Eigen::VectorXd::Constant(matrixSize, 5));

	FDVertex *currVert = NULL;
	int vertID = 0;
	int indexEquation = 0;
	int indexCoefficient = 0;
	
	double mobility = 0;
	double coeff_adjacent = 0; // coefficient for adjacent vertex
	double coeff_center = 0; // coefficient for center vertex

	double deltaX = 0; // dJ/dx
	double deltaY = 0; // dJ/dy

	double KT_div_q = SctmPhys::k0 * this->temperature / SctmPhys::q;
	//we should consider the normalized current equation
	//for electron, J = -u[ n * p_phi/p_x - p_n/p_x ]

	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		indexEquation = equationMap[currVert->GetID()];
		SCTM_ASSERT(indexEquation==iVert, 100012);

		coeff_center = 0;
		// in the boundary cases, the length of one or two direction may be zero 
		deltaX = (currVert->EastLength + currVert->WestLength) / 2;
		deltaY = (currVert->NorthLength + currVert->SouthLength) / 2;
		SCTM_ASSERT(deltaX!=0 && deltaY!=0, 10015);

		//the initial filling method can be used for boundary vertices and non-boundary vertices
		if (   (currVert->EastVertex != NULL)
			 &&( (currVert->NortheastElem == NULL ? false : currVert->NortheastElem->Region->Type == FDRegion::Trapping)
			   ||(currVert->SoutheastElem == NULL ? false : currVert->SoutheastElem->Region->Type == FDRegion::Trapping) )
		   )
		{
			mobility = (mobilityMap[currVert->EastVertex->GetID()] + mobilityMap[currVert->GetID()]) / 2;
			indexCoefficient = equationMap[currVert->EastVertex->GetID()];
			coeff_adjacent = - mobility / currVert->EastLength / deltaX *
				( (potentialMap[currVert->EastVertex->GetID()] - potentialMap[currVert->GetID()]) / 2 - 1 );
			coeff_center += - mobility / currVert->EastLength / deltaX *
				( (potentialMap[currVert->EastVertex->GetID()] - potentialMap[currVert->GetID()]) / 2 + 1 );

			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
		}

		if (   (currVert->WestVertex != NULL)
		     &&( (currVert->NorthwestElem == NULL ? false : currVert->NorthwestElem->Region->Type == FDRegion::Trapping)
			   ||(currVert->SouthwestElem == NULL ? false : currVert->SouthwestElem->Region->Type == FDRegion::Trapping) )
		   )
		{
			mobility = (mobilityMap[currVert->WestVertex->GetID()] + mobilityMap[currVert->GetID()]) / 2;
			indexCoefficient = equationMap[currVert->WestVertex->GetID()];
			coeff_adjacent = - mobility / currVert->WestLength / deltaX *
				( (potentialMap[currVert->WestVertex->GetID()] - potentialMap[currVert->GetID()]) /2 - 1 );
			coeff_center += - mobility / currVert->WestLength / deltaX *
				( (potentialMap[currVert->WestVertex->GetID()] - potentialMap[currVert->GetID()]) /2 + 1 );

			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
		}

		if (   (currVert->NorthVertex != NULL)
			 &&( (currVert->NortheastElem == NULL ? false : currVert->NortheastElem->Region->Type == FDRegion::Trapping)
			   ||(currVert->NorthwestElem == NULL ? false : currVert->NorthwestElem->Region->Type == FDRegion::Trapping) )
		   )
		{
			mobility = (mobilityMap[currVert->NorthVertex->GetID()] + mobilityMap[currVert->GetID()]) / 2;
			indexCoefficient = equationMap[currVert->NorthVertex->GetID()];
			coeff_adjacent = - mobility / currVert->NorthLength / deltaY *
				( (potentialMap[currVert->NorthVertex->GetID()] - potentialMap[currVert->GetID()]) / 2 - 1 );
			coeff_center += - mobility / currVert->NorthLength / deltaY *
				( (potentialMap[currVert->NorthVertex->GetID()] - potentialMap[currVert->GetID()]) / 2 + 1 );

			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
		}

		if (   (currVert->SouthVertex != NULL)
			 &&( (currVert->SoutheastElem == NULL ? false : currVert->SoutheastElem->Region->Type == FDRegion::Trapping)
			   ||(currVert->SouthwestElem == NULL ? false : currVert->SouthwestElem->Region->Type == FDRegion::Trapping))
		   )
		{
			mobility = (mobilityMap[currVert->SouthVertex->GetID()] + mobilityMap[currVert->GetID()]) / 2;
			indexCoefficient = equationMap[currVert->SouthVertex->GetID()];
			coeff_adjacent = - mobility / currVert->SouthLength / deltaY *
				( (potentialMap[currVert->SouthVertex->GetID()] - potentialMap[currVert->GetID()]) / 2 - 1 );
			coeff_center += - mobility / currVert->SouthLength / deltaY *
				( (potentialMap[currVert->SouthVertex->GetID()] - potentialMap[currVert->GetID()]) / 2 + 1 );

			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
		}

		coeff_center += -1 / timeStep; // from p_n/p_t, p=partial differential

		indexCoefficient = equationMap[currVert->GetID()];
 		SCTM_ASSERT(indexCoefficient==indexEquation, 10012);
		//indexCoefficent = indexEquation = vertMap[currVert->GetID]
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_center;
	}
}

void DriftDiffusionSolver::buildRhsVector()
{
	FDVertex *currVert = NULL;
	double rhsVal = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		rhsVal = -lastElecDensMap[currVert->GetID()] / timeStep;
		this->rhsVector.at(iVert) = rhsVal;
	}
}

void DriftDiffusionSolver::refreshRhs()
{
	FDVertex *currVert = NULL;
	int equationID = 0;
	double q = SctmPhys::q;
	double bcVal = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		if (!currVert->BndCond.Valid(FDBoundary::eCurrentDensity))
			continue;

		//the vertex is at current density boundary condition
		equationID = equationMap[currVert->GetID()]; //equationIndex = iVert
		SCTM_ASSERT(currVert->BndCond.GetBCType(FDBoundary::eCurrentDensity) == FDBoundary::BC_Dirichlet, 10013);
		//if (currVert->BndCond.GetBCType(FDBoundary::eCurrentDensity) != FDBoundary::BC_Dirichlet)
		//	continue;

		//the vertex is at BC_Dirichlet boundary condition
		//only solve DD equation in trapping region
		if ( currVert->WestVertex == NULL
			||(  (currVert->WestVertex->NorthwestElem == NULL ? true : currVert->WestVertex->NorthwestElem->Region->Type != FDRegion::Trapping)
			   &&(currVert->WestVertex->SouthwestElem == NULL ? true : currVert->WestVertex->SouthwestElem->Region->Type != FDRegion::Trapping) )
		   )
		{
			bcVal = currVert->BndCond.GetBCValueWestEast(FDBoundary::eCurrentDensity);
			rhsVector.at(equationID) += bcVal / (currVert->EastLength / 2);
		}

		if ( currVert->EastVertex == NULL
			||(  (currVert->EastVertex->NortheastElem == NULL ? true : currVert->EastVertex->NortheastElem->Region->Type != FDRegion::Trapping)
			   &&(currVert->EastVertex->SoutheastElem == NULL ? true : currVert->EastVertex->SoutheastElem->Region->Type != FDRegion::Trapping) )
		   )
		{
			bcVal = currVert->BndCond.GetBCValueWestEast(FDBoundary::eCurrentDensity);
			rhsVector.at(equationID) += - bcVal / (currVert->WestLength / 2);
		}

		if ( currVert->SouthVertex == NULL
			||(  (currVert->SouthVertex->SoutheastElem == NULL ? true : currVert->SouthVertex->SoutheastElem->Region->Type != FDRegion::Trapping)
			   &&(currVert->SouthVertex->SouthwestElem == NULL ? true : currVert->SouthVertex->SouthwestElem->Region->Type != FDRegion::Trapping) )
		   )
		{
			bcVal = currVert->BndCond.GetBCValueSouthNorth(FDBoundary::eCurrentDensity);
			rhsVector.at(equationID) += bcVal / (currVert->NorthLength / 2);
		}

		if ( currVert->NorthVertex == NULL
			||(  (currVert->NorthVertex->NorthwestElem == NULL ? true : currVert->NorthVertex->NorthwestElem->Region->Type != FDRegion::Trapping)
			   &&(currVert->NorthVertex->NortheastElem == NULL ? true : currVert->NorthVertex->NortheastElem->Region->Type != FDRegion::Trapping) )
		   )
		{
			bcVal = currVert->BndCond.GetBCValueSouthNorth(FDBoundary::eCurrentDensity);
			rhsVector.at(equationID) += - bcVal / (currVert->SouthLength / 2);
		}
	}
}

void DriftDiffusionSolver::refreshCoefficientMatrix()
{
	//The boundary conditions are always BC_Dirichlet, so there is no need refreshing the matrix.
}

void DriftDiffusionSolver::getDDVertices(FDDomain *domain)
{
	FDVertex *currVert = NULL;

	for (size_t iVert = 0; iVert != this->totalVertices.size(); ++iVert)
	{
		currVert = this->totalVertices.at(iVert);
		if (
			((currVert->NorthwestElem != NULL) 
			&& (currVert->NorthwestElem->Region->Type == FDRegion::Trapping))
			
			||((currVert->NortheastElem != NULL) 
			&& (currVert->NortheastElem->Region->Type == FDRegion::Trapping))
			
			||((currVert->SoutheastElem != NULL) 
			&& (currVert->SoutheastElem->Region->Type == FDRegion::Trapping))
			
			||((currVert->SouthwestElem != NULL) 
			&& (currVert->SouthwestElem->Region->Type == FDRegion::Trapping))
		   )
		{
			this->vertices.push_back(currVert);
		}
	}
}

void DriftDiffusionSolver::setTimeStep()
{
	timeStep = UtilsTimeStep.NextTimeStep();
}

void DriftDiffusionSolver::fillBackElecDens()
{
	FDVertex *currVert = NULL;
	int equationID = 0;
	double edens = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		equationID = equationMap[currVert->GetID()];
		edens = this->elecDensity.at(equationID);
		currVert->Phys.SetPhysPrpty(PhysProperty::eDensity, edens);
	}
}

DDTest::DDTest(FDDomain *_domain) : DriftDiffusionSolver(_domain)
{
	prepareSolver();
}

void DDTest::prepareSolver()
{
	SctmUtils::UtilsMsg.PrintHeader("Solving Drift-Diffusion");

	int vertSize = this->vertices.size();
	this->rhsVector.resize(vertSize);
	this->elecDensity.resize(vertSize);

	this->temperature = 300; //temporarily used here TODO: modify to accord with the whole simulation
	buildVertexMap();
	setTimeStep();
	buildCoefficientMatrix();
	//UtilsDebug.PrintSparseMatrix(matrixSolver.matrix);
	buildRhsVector();
	UtilsDebug.PrintVector(rhsVector);
	//the structure does not change, so the final coefficient matrix after refreshing does not change either
	refreshCoefficientMatrix();
}

void DDTest::buildVertexMap()
{
	std::pair<VertexMapInt::iterator, bool> insertPairVertex;
	std::pair<VertexMapDouble::iterator, bool> insertPairPrpty;
	int vertID = 0;
	int equationID = 0;
	double phyValue = 0;
	FDVertex *currVert = NULL;

	//this map is filled in order to obtain the vertex index from the vertex internal id. This is useful
	//in setting up the equation, i.e. filling the matrix.
	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		vertID = currVert->GetID();

		equationID = iVert; //equation index is also the vertex index in the vertices container
		insertPairVertex = this->equationMap.insert(VertexMapInt::value_type(vertID, equationID));
		SCTM_ASSERT(insertPairVertex.second==true, 10011);

		insertPairPrpty = this->mobilityMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys.GetPhysPrpty(PhysProperty::eMobility)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->potentialMap.insert(VertexMapDouble::value_type(vertID, 0)); // suppose no electronic field exists
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->lastElecDensMap.insert(VertexMapDouble::value_type(vertID, 0)); // suppose the initial carrier density is 0
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);
	}
}

void DDTest::setBndCondCurrent()
{
	FDVertex *currVert = NULL;
	int equaID = 0;
	double bcVal = -1e-3; 
	// current density at boundary for test, in [A/cm^2]
	// note that the current direction is the reversed direction of electron flow

	Normalization norm = Normalization();

	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		//for the current tunneling from tunneling oxide
		if (  currVert->IsAtBoundary(FDBoundary::eCurrentDensity) 
			&&( (currVert->SoutheastElem == NULL ? false : currVert->SoutheastElem->Region->Type == FDRegion::Tunneling)
			||(currVert->SouthwestElem == NULL ? false : currVert->SouthwestElem->Region->Type == FDRegion::Tunneling) )
			)
		{
			currVert->BndCond.RefreshBndCondValue(FDBoundary::eCurrentDensity, 0, norm.PushCurrDens(bcVal));
		}
	}
}

void DDTest::SolveDD()
{
	setBndCondCurrent();
	refreshRhs();
	UtilsDebug.PrintVector(this->rhsVector, "right hand side vector");
	this->matrixSolver.SolveMatrix(rhsVector, this->elecDensity);
	UtilsDebug.PrintVector(this->elecDensity, "electron density");
	fillBackElecDens();

	SctmFileOperator write = SctmFileOperator("C:\\Users\\Lunzhy\\Desktop\\SctmTest\\DDTest\\eDensity.txt", SctmFileOperator::Write);
	write.WriteDDResultForOrigin(this->vertices, "electron density");
}
