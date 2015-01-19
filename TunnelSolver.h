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
#include <string>
#include <algorithm>

class FDDomain;
class FDVertex;
using std::vector;
typedef std::map<int, double> VertexMapDouble; // <vertID, value>, map of vertex physical value


/// @brief TunnelSolver is the base class dealing with tunneling mechanism
///
/// It should be note that all the variables in TunnelSolver and its derived class are in real value.
/// They are transformed to real values after being passed into this solver.
/// And the result after calculation is in U.I.
class SubsToTrapElecTAT;
class TunnelSolver
{
	//friend class SctmUtils::SctmDebug;
public:
	enum TunnelMode
	{
		ElecTunnel,
		HoleTunnel,
	};
	enum TunnelDirection
	{
		North,
		South,
		East,
		West,
		NoTunnel,
	};
	TunnelSolver(FDDomain *_domain);
	void ReadInput(VertexMapDouble &fermi);
	virtual void SolveTunnel() = 0;
	virtual void ReturnResult(VertexMapDouble &ret) = 0;
	virtual void ReturnResult_MFN(VertexMapDouble &ret) {};
	virtual void ReturnResult_B2T(VertexMapDouble &ret) {};
	virtual void ReturnResult_T2B(VertexMapDouble &ret) {};

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
	double getTransCoeff(double energy, vector<double> &deltax, vector<double> &emass, vector<double> &cbedge, int size = 0, int startindex = 0); //Transmission coefficient
	
	double calcDTFNtunneling(vector<double> &deltaX, vector<double> &emass, vector<double> &cbedge, double cbedgeMax);
	double calcThermalEmission(vector<double> &deltaX, vector<double> &emass, vector<double> &cbedge, double cbedgeMin);

	void loadBandStructure(FDVertex *startVert);

	double supplyFunction_forCurrDens(double energy);
	double supplyFunction_forTunCoeff(double energy);

protected:
	double temperature;
	FDDomain *domain;
	bool tunnelTrapToGateEnable;
	TunnelMode tunMode;

	vector<FDVertex *> vertsTunnelOxideStart;
	vector<FDVertex *> vertsTunnelOxideEnd;
	vector<FDVertex *> vertsBlockOxideStart;
	vector<FDVertex *> vertsBlockOxideEnd;

	//the fermi above values are in normalized value
	VertexMapDouble efermiAboveMap; // fermi energy - conduction band 
	VertexMapDouble hfermiAboveMap; // hole fermi energy - valence band, offset value

	//cbEdge, vbEdge, eMass, hMass
	//for tunneling layer
	vector<double> bandEdge_Tunnel;
	vector<double> ehmass_Tunnel; ///< electron or hole mass
	vector<double> deltaX_Tunnel;

	//for trapping layer
	//the additional vertex for solving modified Fowler-Nordheim tunneling
	vector<double> bandEdge_Trap;
	vector<double> ehMass_Trap;
	vector<double> deltaX_Trap;
	vector<double> trapEnergyLevel; ///< trap energy level
	vector<FDVertex *> verts_Trap;

	//for blocking layer
	vector<double> bandEdge_Block;
	vector<double> ehMass_Block;
	vector<double> deltaX_Block;

	//these properties are used in the calculation of supply function
	double bandEdgeTunnelFrom; ///< left electrode conduction band edge
	double bandEdgeTunnelTo;

	double fermiEnergyTunnelFrom;
	double fermiEnergyTunnelTo;
	double effTunnelMass; ///< effective mass

	TunnelDirection tunDirection;
};

class SubsToTrapElecTunnel : public TunnelSolver
{
	friend class SubsToTrapElecTAT;
public:
	SubsToTrapElecTunnel(FDDomain *_domain);
	virtual void SolveTunnel();
	void ReturnResult(VertexMapDouble &ret);
	void ReturnResult_MFN(VertexMapDouble &ret);
	void ReturnResult_B2T(VertexMapDouble &ret);
	void ReturnResult_T2B(VertexMapDouble &ret);
protected:
	double getSupplyFunction(double energy);
	void setSolver_DTFN(FDVertex *startVertex); //This method is not currently used.
	/// @brief setSolver_Trap is used to set the vectors for solving mechanisms including trapping layer.
	/// 
	/// Setting solver of trapping layer is the preparation for solve MFN and B2T tunneling problems.
	/// 
	/// @pre the band edge / emass / deltax containers are filled.
	/// @return void
	/// @note
	void setSolver_Trap();
	virtual void setTunnelDirection(FDVertex *vertSubs, FDVertex *vertTrap);
	virtual void setTunnelTag();
	void calcCurrDens_MFN();
	void calcCurrDens_B2T();
	void calcTransCoeff_T2B();
	double updateTunCurrForCylindrical(double currdens, FDVertex* inVert, FDVertex* outVert);

	FDVertex *findTrapVertex_MFN(double energy, int &size);
	FDVertex *findTrapVertex_B2T(double energy, int &size);

	double eSubsBarrier;

	vector<double> cbEdge_TunnelTrap;
	vector<double> eMass_TunnelTrap;
	vector<double> deltaX_TunnelTrap;

	vector<double> eCurrDens_DTFN; // the sequence is the same with the vertex in verticsTunnelStart
	VertexMapDouble eCurrDensMap_MFN; ///< map for MFN tunneling current, in [A/cm^2]
	VertexMapDouble eCurrDensMap_B2T; ///< map for the electron current density from substrate in calculation of band-to-trap tunneling, in [A/cm^2]
	VertexMapDouble eTransCoeffMap_T2B; ///< map for Trap-to-Band tunneling out from trap site to substrate, especially in Retention.
	VertexMapDouble eCoeffMap_TAT2B; ///< map for Trap-Assisted Trap-to-Band tunneling out from trap site to substrate

	SubsToTrapElecTAT* solverTAT; 
};

class SubsToTrapHoleTunnel : public SubsToTrapElecTunnel
{
public:
	SubsToTrapHoleTunnel(FDDomain* _domain);
	void SolveTunnel();
protected:
	void setTunnelDirection(FDVertex* vertSubs, FDVertex* vertTrap);
	void setTunnelTag();
	void pretendToBeElecTun();

	double hSubsBarrier;
	vector<double> hCurrDens_DTFN;
};

class TrapToGateElecTunnel : public TunnelSolver
{
public:
	TrapToGateElecTunnel(FDDomain *_domain);
	virtual void SolveTunnel();
	void ReturnResult(VertexMapDouble &ret);
	void ReturnResult_T2B(VertexMapDouble &ret);
protected:
	double getSupplyFunction(double energy);
	void setSolver_DTFN(FDVertex *endVertex); //This method is not currently used.
	virtual void setTunnelDirection(FDVertex *vertTrap, FDVertex *vertGate);
	virtual void setTunnelTag();
	void setSolver_Trap();
	void calcTransCoeff_T2B();
	double updateTunCurrForCylindrical(double currdens, FDVertex* inVert, FDVertex* outVert);

	vector<double> cbEdge_TrapBlock;
	vector<double> eMass_TrapBlock;
	vector<double> deltaX_TrapBlock;

	vector<double> eCurrDens_DTFN; // the sequence is the same with the vertex in verticsTunnelStart
	VertexMapDouble eTransCoeffMap_T2B;
};

class TrapToGateHoleTunnel : public TrapToGateElecTunnel
{
public:
	TrapToGateHoleTunnel(FDDomain* _domain);
	void SolveTunnel();
protected:
	void setTunnelDirection(FDVertex *vertTrap, FDVertex *vertGate);
	void setTunnelTag();

	void pretendToBeElecTun();
	vector<double> hCurrDens_DTFN;
};

class SlopingTunnelTrapToGate
{
public:
	SlopingTunnelTrapToGate(FDDomain *_domain, FDVertex *_vert, TunnelSolver::TunnelMode _tunmode);
	void LoadBandStructureAlongPath(vector<double> &dx, vector<double> &cbedge, vector<double> &emass);
	FDVertex* GetGateVertex();

	static bool IsSlopingTunnel(FDVertex *vert);
protected:
	FDDomain *domain;
	double temperature;
	TunnelSolver::TunnelMode tunMode;
	FDVertex *vertStart;
	std::string vertRegName;

	double tanAngle;
	FDVertex *upmostVert;
	FDVertex *easternmostVert;
	FDVertex *westernmostVert;
	
	vector<double> dSlope; //the size of deltaX should be one less that the other two
	vector<double> bandEdge;
	vector<double> ehMass;

	void setBoundaryVerts();
	void setInterpolatedValues();
	void findLeftRigthVertex(FDVertex *vert, double xDist, FDVertex* &left, FDVertex* &right);
};

class SubsToTrapElecTAT
{
public:
	SubsToTrapElecTAT(SubsToTrapElecTunnel* _tunSolver);
	double SolveTAT(FDVertex* tunStart, FDVertex* tunEnd);
	double CalcCoeff_TrapAssistedT2B(FDVertex* trapvert);

protected:
	void setOxideTrapParam();
	void loadOxideVertices();
	bool calcTimeConstant(FDVertex* oxideVert, double& ctime, double& etime);
	double getTunCoeffTAT(FDVertex* startVert, FDVertex* endVert, std::string inout);
	double getDeltaX(FDVertex* vert);

	//for trap-assisted trap-to-band tunneling
	void calcTimeConstant_TrapAssistedT2B(FDVertex* oxideVert, FDVertex* trapVert, double& ctime, double& etime);
	double getTunCoeff_TrapAssistedT2B(FDVertex* oxideVert, FDVertex* trapVert);

protected:
	SubsToTrapElecTunnel* tunnelSolver;
	TunnelSolver::TunnelDirection tunnelDirection;
	FDVertex* tunStartVert;
	FDVertex* tunEndVert;
	vector<FDVertex*> oxideVertices;
	VertexMapDouble oxideTrapOccupation;

	double oxideTrapCrossSection;
	double oxideTrapDensity; // in [1/cm^2]
	double oxideTrapEnergyFromCond;
	double trapAssistedT2BTunnelFrequency;
};

#endif