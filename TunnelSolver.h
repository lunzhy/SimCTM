/**
* @file TunnelSolver.h
* @brief This file contains the classes that account for solving tunneling problems
*
*
*
* @author
* @version 
* @date 2013-7-16   10:39
* @note
* @todo
*/
#ifndef _TUNNELSOLVER_H_
#define _TUNNELSOLVER_H_

#include <vector>
#include <map>

class FDDomain;
class FDVertex;
using std::vector;
typedef std::map<int, double> VertexMapDouble; // <vertID, value>, map of vertex physical value


/// @brief TunnelSolver is the base class dealing with tunneling mechanism
///
/// It should be note that all the variables in TunnelSolver and its derived class are in real value.
/// They are transformed to real values after being passed into this solver.
/// And the result after calculation is in U.I.
class TunnelSolver
{
	//friend class SctmUtils::SctmDebug;
public:
	TunnelSolver(FDDomain *_domain);
	void ReadInput(VertexMapDouble &fermi);
	virtual void SolveTunnel() = 0;
	virtual void ReturnResult(VertexMapDouble &ret) = 0;
	virtual void ReturnResult_MFN(VertexMapDouble &ret) {};
	virtual void ReturnResult_B2T(VertexMapDouble &ret) {};

protected:
	/// @brief initialize is used to initialize the solver
	/// 
	/// The initialization of the solver is called once in the initialization of the
	/// tunnel solver (base class) object, so only the properties which stays unchanged 
	/// during the following simulation is processed in this method.
	/// Initialization is used to find the starting vertices for tunneling, i.e., the 
	/// vertices at substrate. The other corresponding vertices is set when loading
	/// the band structures.
	///
	/// @pre
	/// @return void
	/// @note initialization is done in the base solver, so for each specific tunnel solver,
	/// the method to find tunneling start/end vertices are the same. That means vertices at
	/// channel and gate must appear in pairs.
	void initialize();
	virtual double getSupplyFunction(double energy);
	double getTransCoeff(double energy, vector<double> &deltax, vector<double> &emass, vector<double> &cbedge, int size = 0); //Transmission coefficient
	
	double calcDTFNtunneling(vector<double> &deltaX, vector<double> &emass, vector<double> &cbedge);
	double calcThermalEmission(vector<double> &deltaX, vector<double> &emass, vector<double> &cbedge);

	void loadBandStructure(FDVertex *startVert);

protected:
	double temperature;
	FDDomain *domain;
	vector<FDVertex *> vertsStart;
	vector<FDVertex *> vertsEnd;

	vector<FDVertex *> vertsTunnelOxideStart;
	vector<FDVertex *> vertsTunnelOxideEnd;
	vector<FDVertex *> vertsBlockOxideStart;
	vector<FDVertex *> vertsBlockOxideEnd;

	VertexMapDouble fermiAboveMap; // fermi energy - conduction band 

	//for tunneling layer
	vector<double> cbEdge_Tunnel;
	vector<double> eMass_Tunnel;
	vector<double> deltaX_Tunnel;

	//for trapping layer
	//the additional vertex for solving modified Fowler-Nordheim tunneling
	vector<double> cbEdge_Trap;
	vector<double> eMass_Trap;
	vector<double> deltaX_Trap;
	vector<double> eEnergyLevel_Trap; ///< trap energy level
	vector<FDVertex *> verts_Trap;

	//for blocking layer
	vector<double> cbEdge_Block;
	vector<double> eMass_Block;
	vector<double> deltaX_Block;

	vector<double> eCurrDens_DTFN; // the sequence is the same with the vertex in verticsTunnelStart

	double cbedgeTunnelFrom; ///< left electrode conduction band edge
	double cbedgeTunnelTo;
	double fermiEnergyTunnelFrom;
	double fermiEnergyTunnelTo;
	double effTunnelMass; ///< effective mass
	double eCurrDens; ///< the tunneling current density, in [A/m^2]
};

class SubsToTrapElecTunnel : public TunnelSolver
{
public:
	SubsToTrapElecTunnel(FDDomain *_domain);
	void SolveTunnel();
	void ReturnResult(VertexMapDouble &ret);
	void ReturnResult_MFN(VertexMapDouble &ret);
	void ReturnResult_B2T(VertexMapDouble &ret);
protected:
	double getSupplyFunction(double energy);
	void setSolver_DTFN(FDVertex *startVertex);
	/// @brief setSolver_Trap is used to set the vectors for solving mechanisms including trapping layer.
	/// 
	/// Setting solver of trapping layer is the preparation for solve MFN and B2T tunneling problems.
	/// 
	/// @param FDVertex * startVertex
	/// @pre
	/// @return void
	/// @note
	void setSolver_Trap(FDVertex *startVertex);
	void setTunnelTag();
	void calcCurrDens_MFN();
	void calcCurrDens_B2T();

	FDVertex *findTrapVertex_MFN(double energy, int &size);
	FDVertex *findTrapVertex_B2T(double energy, int &size);

	vector<double> cbEdge_TunnelTrap;
	vector<double> eMass_TunnelTrap;
	vector<double> deltaX_TunnelTrap;

	VertexMapDouble eCurrDensMap_MFN; ///< map for MFN tunneling current, in [A/cm^2]
	VertexMapDouble eCurrDensMap_B2T; ///< map for the electron current density from substrate in calculation of band-to-trap tunneling, in [A/cm^2]
};

class TrapToGateElecTunnel : public TunnelSolver
{
public:
	TrapToGateElecTunnel(FDDomain *_domain);
	void SolveTunnel();
	void ReturnResult(VertexMapDouble &ret);
protected:
	double getSupplyFunction(double energy);
	void setSolver_DTFN(FDVertex *endVertex);
	void setTunnelTag();

	vector<double> eCurrDensMap_T2B;
};

#endif