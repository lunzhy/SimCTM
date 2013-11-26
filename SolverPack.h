/**
* @file SolverPack.h
* @brief
*
*
*
* @author
* @version 
* @date 2013-11-19   16:07
* @note
* @todo
*/
#include <map>

class FDDomain;
class TwoDimPoissonSolver;
class TunnelSolver;
class DriftDiffusionSolver;

typedef std::map<int, int> VertexMapInt; // <vertID, int>
typedef std::map<int, double> VertexMapDouble; // <vertID, double>

class SolverPack
{
public:
	SolverPack(FDDomain *_domain);
	void Run();
protected:
	FDDomain *domain;
	TwoDimPoissonSolver *poissonSolver;
	TunnelSolver *TunnelOxideSolver;
	TunnelSolver *BlockOxideSolver;
	DriftDiffusionSolver *ddSolver;

protected:
	void initialize();
	void callIteration();
	
	void fetchPoissonResult();
	void fetchTunnelOxideResult();
	void fetchBlockOxideResult();
	void fetchDDResult();

	void fakeFermiEnergy();
private:
	VertexMapDouble mapSiFermiAboveCBedge; // for input in the tunneling solver silicon fermi energy - silicon conduction band edge
	VertexMapDouble mapPotential;
	VertexMapDouble mapCurrDensFromTunnelLayer;
	VertexMapDouble mapCurrDensCoeff; // the coefficient to calculate current density for dd solver, in [A*cm]
};