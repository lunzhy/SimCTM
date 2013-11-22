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
typedef std::map<int, double> VertexMapDouble;
class TunnelSolver
{
	//friend class SctmUtils::SctmDebug;
public:
	TunnelSolver(FDDomain *_domain);
	void ReadInput(VertexMapDouble fermi);
	void SolveTunnel_Interface();

protected:
	virtual double getSupplyFunction(double energy);
	double getTransCoeff(double energy); //Transmission coefficient
	
	double calcDTFNtunneling();
	double calcThermalEmission();

	virtual void setSolver_Interface(FDVertex *startVertex) = 0;
	double solve_Interface(); // solver the interface tunneling current for specified vertex


protected:
	FDDomain *domain;
	vector<FDVertex *> vertsTunnelStart;
	
	vector<FDVertex *> vertsTunnelEnd_Interface;
	VertexMapDouble fermiAboveMap; // fermi energy - conduction band 

	vector<double> cbEdge;
	vector<double> elecMass;
	vector<double> deltaX;
	vector<double> eCurrDens_Interface; // the sequence is the same with the vertex in verticsTunnelStart

	double cbedgeTunnelFrom; ///< left electrode conduction band edge
	double cbedgeTunnelTo;
	double fermiEnergyTunnelFrom;
	double fermiEnergyTunnelTo;
	double effTunnelMass; ///< effective mass
	double temperature;
	double eCurrDens; ///< the tunneling current density, in [A/cm^2]
};

class SubsToGateEletronTunnel : public TunnelSolver
{
public:
	SubsToGateEletronTunnel(FDDomain *_domain);
protected:
	void initialize(); ///< initialize only exists in derived class
	double getSupplyFunction(double energy);
	void setSolver_Interface(FDVertex *startVertex);
};

#endif