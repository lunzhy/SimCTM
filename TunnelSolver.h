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

protected:
	virtual double getSupplyFunction(double energy);
	double getTransCoeff(double energy); //Transmission coefficient
	
	double calcDTFNtunneling();
	double calcThermalEmission();

	virtual void setSolver_Tunnel(FDVertex *startOrEndVertex) = 0;
	double solveCurrDens_Tunnel(); //solver the tunneling current when solver is ready

protected:
	FDDomain *domain;
	vector<FDVertex *> vertsStart_Tunnel;
	vector<FDVertex *> vertsEnd_Tunnel;

	vector<FDVertex *> vertsStart_Trap;
	vector<FDVertex *> vertsEnd_Trap;

	VertexMapDouble fermiAboveMap; // fermi energy - conduction band 

	vector<double> cbEdge;
	vector<double> elecMass;
	vector<double> deltaX;
	vector<double> eCurrDens_Tunnel; // the sequence is the same with the vertex in verticsTunnelStart

	double cbedgeTunnelFrom; ///< left electrode conduction band edge
	double cbedgeTunnelTo;
	double fermiEnergyTunnelFrom;
	double fermiEnergyTunnelTo;
	double effTunnelMass; ///< effective mass
	double temperature;
	double eCurrDens; ///< the tunneling current density, in [A/cm^2]
};

class SubsToTrapElecTunnel : public TunnelSolver
{
public:
	SubsToTrapElecTunnel(FDDomain *_domain);
	void SolveTunnel();
	void ReturnResult(VertexMapDouble &ret);
protected:
	void initialize(); ///< initialize only exists in derived class
	double getSupplyFunction(double energy);
	void setSolver_Tunnel(FDVertex *startVertex);
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
	void setSolver_Tunnel(FDVertex *endVertex);
};

#endif