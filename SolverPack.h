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
#ifndef _SOLVERPACK_H_
#define _SOLVERPACK_H_

#include <map>

class FDDomain;
class TwoDimPoissonSolver;
class TunnelSolver;
class DriftDiffusionSolver;
class TrapSolver;
class OneDimSubsSolver;

typedef std::map<int, int> VertexMapInt; // <vertID, int>
typedef std::map<int, double> VertexMapDouble; // <vertID, double>

class SolverPack
{
public:
	SolverPack(FDDomain *_domain);
	void Run();
protected:
	double temperature;
	string simStructure;
	FDDomain *domain;

	OneDimSubsSolver *subsSolver;
	TwoDimPoissonSolver *poissonSolver;
	TunnelSolver *tunnelOxideSolver;
	TunnelSolver *blockOxideSolver;
	TrapSolver *trappingSolver;
	DriftDiffusionSolver *ddSolver;

protected:
	void initialize();
	void callIteration();
	
	void fetchSubstrateResult();
	void fetchPoissonResult();
	void fetchTunnelOxideResult();
	void fetchBlockOxideResult();
	void fetchDDResult();
	void fetchTrappingResult();
	void readSubstrateFromFile();

private:
	VertexMapDouble mapChannelPotential; ///< the potential of channel vertices
	VertexMapDouble mapSiFermiAboveCBedge; ///< for input in the tunneling solver silicon fermi energy - silicon conduction band edge
	VertexMapDouble mapPotential; ///< map for potential
	VertexMapDouble mapCurrDens_Tunnel; ///< the tunneling current density for FN/DT tunneling, (Program, in [A/cm^2]) or tunneling coefficient (Retention, in [A*cm]) across the tunneling oxide
	VertexMapDouble mapCurrDensCoeff_Block; ///< the coefficient to calculate current density for dd solver, in [A*cm], in trap-to-gate tunnel solver
	VertexMapDouble mapCurrDensMFN; ///< the current density tunneling into trapping layer (Modified Fowler-Nordheim), in [A/cm^2], in subs-to-trapLayer tunnel solver
	VertexMapDouble mapCurrDensB2T; ///< tunneling current density in band-to-trap tunneling, in subs-to-trapLayer tunnel solver
	VertexMapDouble mapTransCoeffT2B_Block; ///< transmission coefficient in Trap-to-Band tunneling out, in trapLayer-to-gate tunnel solver
	VertexMapDouble mapTransCoeffT2B_Tunnel; ///< transmission coefficient in Trap-to-Band tunneling out from tunnling oxide. 
};

#endif