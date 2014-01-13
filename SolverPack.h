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

private:
	VertexMapDouble mapChannelPotential; ///< the potential of channel vertices
	VertexMapDouble mapSiFermiAboveCBedge; ///< for input in the tunneling solver silicon fermi energy - silicon conduction band edge
	VertexMapDouble mapPotential; ///< map for potential
	VertexMapDouble mapCurrDensFromTunnelLayer; ///< the tunneling current density across the tunneling layer, for FN/DT tunneling, in [A/cm^2]
	VertexMapDouble mapCurrDensCoeff; ///< the coefficient to calculate current density for dd solver, in [A*cm]
	VertexMapDouble mapCurrDensMFN; ///< the current density tunneling into trapping layer (Modified Fowler-Nordheim), in [A/cm^2]
	VertexMapDouble mapCurrDensB2T; ///< tunneling current density in band-to-trap tunneling
};

#endif