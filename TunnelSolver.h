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

class TunnelSolver
{
	//friend class SctmUtils::SctmDebug;
public:
	TunnelSolver(FDDomain *_domain);
	void ReadInput(VertexMapDouble &fermi);
	virtual void SolveTunnel() = 0;
	virtual void ReturnResult(VertexMapDouble &ret) = 0;
	virtual void ReturnResult_MFN(VertexMapDouble &ret) {};

protected:
	virtual double getSupplyFunction(double energy);
	double getTransCoeff(double energy, vector<double> &deltax, vector<double> &emass, vector<double> &cbedge); //Transmission coefficient
	
	double calcDTFNtunneling();
	double calcThermalEmission();

	//virtual void setSolver_DTFN(FDVertex *startOrEndVertex) = 0;
	double calcCurrDens_DTFN(); //solver the tunneling current when solver is ready

protected:
	double temperature;
	FDDomain *domain;
	vector<FDVertex *> vertsStart;
	vector<FDVertex *> vertsEnd;

	VertexMapDouble fermiAboveMap; // fermi energy - conduction band 

	vector<double> cbEdge_Oxide;
	vector<double> eMass_Oxide;
	vector<double> deltaX_Oxide;
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
protected:
	/// @brief initialize is used to initialize the solver
	/// 
	/// The initialization of the solver is called once in the initialization of the
	/// tunnel solver object, so only the properties which stays unchanged during the
	/// following simulation is processed in this method.
	/// 
	/// @pre
	/// @return void
	/// @note initialization only exists in derived class.
	void initialize();
	double getSupplyFunction(double energy);
	void setSolver_DTFN(FDVertex *startVertex);

	void setSolver_MFN(FDVertex *startVertex);
	void calcCurrDens_MFN();

	FDVertex *findTrapVertex_MFN(double energy);

	//the additional vertex for solving modified Fowler-Nordheim tunneling
	vector<double> cbEdge_Trap;
	vector<double> eMass_Trap;
	vector<double> deltaX_Trap;
	vector<FDVertex *> verts_Trap;
	VertexMapDouble eCurrDensMap_MFN;
};

class TrapToGateElecTunnel : public TunnelSolver
{
public:
	TrapToGateElecTunnel(FDDomain *_domain);
	void SolveTunnel();
	void ReturnResult(VertexMapDouble &ret);
protected:
	void initialize();
	double getSupplyFunction(double energy);
	void setSolver_DTFN(FDVertex *endVertex);
};

#endif