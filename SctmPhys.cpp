/**
* @file SctmPhys.cpp
* @brief This file contains the detailed implementation of the physics problems in the simulation
*
*
*
* @author
* @version 
* @date 2013-8-5   15:22
* @note
* @todo
*/

#pragma once
#include "SctmPhys.h"
#include "Material.h"
#include "SctmUtils.h"
#include "DomainDetails.h"
#include "FDDomain.h"
#include "Normalization.h"
#include "SctmMath.h"

namespace SctmPhys
{
	const double &k0 = BoltzmanConstant;
	const double &h = PlanckConstant;
	const double &eps = VacuumDielectricConstant;
	const double &T0 = RoomTemperature;
	const double &q = ElementaryCharge;
	const double &ni = IntrinsicConcentration;
	
	double ReferencePotential;

	PhysProperty::PhysProperty(FDVertex *_vert)
	{
		vertSelf = _vert;
		bandgap = 0;
		electrostaticPotential = 0;
		//conductionBandEnergy = 0;
		//valenceBandEnergy = 0;
		electronAffinity = 0;
		e_mass = 0;
		h_mass = 0;
		//netCharge = 0;
		e_density = 0;
		e_mobility = 0;
		controlArea = 0;
	}

	void PhysProperty::SetPhysPrpty(Name prptyName, double prptyValue)
	{
		switch (prptyName)
		{
		case ElectrostaticPotential:
			electrostaticPotential = prptyValue;
			break;
		case eMass:
			e_mass = prptyValue;
			break;
		case hMass:
			h_mass = prptyValue;
			break;
		case ElectronAffinity:
			electronAffinity = prptyValue;
			break;
		case Bandgap:
			bandgap = prptyValue;
			break;
		case eMobility:
			e_mobility = prptyValue;
			break;
		case eDensity:
			e_density = prptyValue;
			break;
		case TunnelCoeff:
			tunnelCoeff = prptyValue;
			break;
		default:
			SCTM_ASSERT(SCTM_ERROR, 10019);
		}
	}

	double PhysProperty::GetPhysPrpty(Name prptyName) const
	{
		static double RefPotential = SctmPhys::ReferencePotential;

		double energy = 0;
		double affinity = 0;
		double pot = 0;
		double ret = 0;

		switch (prptyName)
		{
			case ElectrostaticPotential:
			{
				ret = electrostaticPotential;
				break;
			}
			case ConductionBandEnergy:
			{
				//Ec = -X-q(phi-phiRef)
				pot = this->electrostaticPotential;
				affinity = this->electronAffinity;
				energy = -affinity - (pot-RefPotential);
				ret = energy;
				break;
			}
			case ValenceBandEnergy:
			{
				//Ev = Ev-Eg = -X-q(phi-phiRef)-Eg
				ret =  GetPhysPrpty(ConductionBandEnergy) - bandgap;
				break;
			}
			case eMass:
			{
				ret = e_mass;
				break;
			}
			case hMass:
			{
				ret = h_mass;
				break;
			}	
			case ElectronAffinity:
			{
				ret = electronAffinity;
				break;
			}
			case Bandgap:
			{
				ret = bandgap;
				break;
			}
			case NetCharge:
			{
				//there should be four parts
				ret = - e_density;
				break;
			}	
			case eMobility:
			{
				ret = e_mobility;
				break;
			}
			case eDensity:
			{
				ret = e_density;
				break;
			}
			case DensityControlArea:
			{
				ret = controlArea;
				break;
			}
			case ElectricField_X:
			{
				double bndNormX = 0; // inner vertex
				if (vertSelf->IsAtBoundary(FDBoundary::Potential))
				{
					bndNormX = vertSelf->BndCond.GetBndDirection(FDBoundary::Potential).X();
				}
				//bnd_normX < 0, west boundary
				//bnd_normX > 0, east boundary
				//bnd_normX = 0, inner vertex or vertex at boundary but inside in x direction
				if ( bndNormX == 0 )
				{
					// for boundary vertex in terms of X direction
					double hw = vertSelf->WestLength;
					double he = vertSelf->EastLength;
					double fw = vertSelf->WestVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fe = vertSelf->EastVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fc = GetPhysPrpty(ElectrostaticPotential);
					SCTM_ASSERT(hw!=0 && he!=0, 10024);
					//Ex = - p_phi / p_x
					ret = - ( -he*he * fw + (he*he - hw*hw) * fc + hw*hw * fe ) / ( he*hw*(he + hw) );
					//ret = - ( fe - fw ) / ( he + hw );
				}
				else if (bndNormX < 0)
				{
					SCTM_ASSERT(vertSelf->EastVertex->EastVertex!=NULL, 10025);
					double he = vertSelf->EastLength;
					double hee = vertSelf->EastVertex->EastLength;
					double fe = vertSelf->EastVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fee = vertSelf->EastVertex->EastVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fc = GetPhysPrpty(ElectrostaticPotential);
					SCTM_ASSERT(hee!=0 && he!=0, 10024);
					ret = - ( -hee*(2*he + hee) * fc + (he + hee)*(he + hee) * fe - he*he * fee ) / ( he*hee*(he + hee) );
				}
				else //  (bndNormX > 0)
				{
					SCTM_ASSERT(vertSelf->WestVertex->WestVertex!=NULL, 10025);
					double hw = vertSelf->WestLength;
					double hww = vertSelf->WestVertex->WestLength;
					double fw = vertSelf->WestVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fww = vertSelf->WestVertex->WestVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fc = GetPhysPrpty(ElectrostaticPotential);
					SCTM_ASSERT(hw!=0 && hww!=0, 10024);
					ret = - ( hw*hw * fww - (hw + hww)*(hw + hww) * fw + hww*(2*hw + hww) * fc ) / ( hw*hww*(hw + hww) );
				}
				break;
			}
			case ElectricField_Y:
			{
				double bndNormY = 0;
				if (vertSelf->IsAtBoundary(FDBoundary::Potential))
				{
					bndNormY = vertSelf->BndCond.GetBndDirection(FDBoundary::Potential).Y();
				}
				//bnd_normY < 0, south boundary
				//bnd_normY > 0, north boundary
				//bnd_normY = 0, inner vertex
				if ( bndNormY == 0 )
				{
					double hs = vertSelf->SouthLength;
					double hn = vertSelf->NorthLength;
					double fs = vertSelf->SouthVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fn = vertSelf->NorthVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fc = GetPhysPrpty(ElectrostaticPotential);
					SCTM_ASSERT(hs!=0 && hn!=0, 10024);
					//Ey = - p_phi / p_y
					ret = - ( -hn*hn * fs + (hn*hn - hs*hs) * fc + hs*hs * fn ) / ( hn*hs*(hn + hs) );
				}
				else if ( bndNormY < 0 )
				{
					SCTM_ASSERT(vertSelf->NorthVertex->NorthVertex!=NULL, 10025);
					double hn = vertSelf->NorthLength;
					double hnn = vertSelf->NorthVertex->NorthLength;
					double fn = vertSelf->NorthVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fnn = vertSelf->NorthVertex->NorthVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fc = GetPhysPrpty(ElectrostaticPotential);
					SCTM_ASSERT(hnn!=0 && hn!=0, 10024);
					ret = - ( -hnn*(2*hn + hnn) * fc + (hn + hnn)*(hn + hnn) * fn - hn*hn * fnn ) / ( hn*hnn*(hn + hnn) );
				}
				else // ( bndNormX > 0 )
				{
					SCTM_ASSERT(vertSelf->SouthVertex->SouthVertex!=NULL, 10025);
					double hs = vertSelf->SouthLength;
					double hss = vertSelf->SouthVertex->SouthLength;
					double fs = vertSelf->SouthVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fss = vertSelf->SouthVertex->SouthVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fc = GetPhysPrpty(ElectrostaticPotential);
					SCTM_ASSERT(hs!=0 && hss!=0, 10024);
					ret = - ( hs*hs * fss - (hs + hss)*(hs + hss) * fs + hss*(2*hs + hss) * fc ) / ( hs*hss*(hs + hss) );
				}
				break;
			}
			case ElectricField:
			{
				double elecFieldX = GetPhysPrpty(ElectricField_X);
				double elecFieldY = GetPhysPrpty(ElectricField_Y);
				ret = SctmMath::sqrt( elecFieldX * elecFieldX + elecFieldY * elecFieldY );
				break;
			}
			case eCurrentDensity_X:
			{
				double bndNormX = 0;
				if ( vertSelf->IsAtBoundary(FDBoundary::eDensity) )
				{
					bndNormX = vertSelf->BndCond.GetBndDirection(FDBoundary::eDensity).X();
				}
				//bnd_normX < 0, west boundary
				//bnd_normX > 0, east boundary
				//bnd_normX = 0, inner vertex or vertex at boundary but inside in x direction
				if ( bndNormX == 0 )
				{
					double elecFieldX = GetPhysPrpty(ElectricField_X);

					double hw = vertSelf->WestLength;
					double he = vertSelf->EastLength;
					double nw = vertSelf->WestVertex->Phys->GetPhysPrpty(eDensity);
					double ne = vertSelf->EastVertex->Phys->GetPhysPrpty(eDensity);
					double nc = GetPhysPrpty(eDensity);
					SCTM_ASSERT(hw!=0 && he!=0, 10024);
					double pn_div_px = ( -he*he * nw + (he*he - hw*hw) * nc + hw*hw * ne ) / ( he*hw*(he + hw) );
					double mobility = GetPhysPrpty(eMobility);

					// J = -u( n * p_phi/p_x + p_n / p_x )
					ret = mobility * ( nc * elecFieldX - pn_div_px);
				}
				else
				{
					double bcNormX = vertSelf->BndCond.GetBCNormVector(FDBoundary::eDensity).X();

					if ( vertSelf->BndCond.GetBCTunnelTag() == FDBoundary::eTunnelIn )
					{
						double bcValue = vertSelf->BndCond.GetBCValue(FDBoundary::eDensity);
						double currDens = bcValue * bcNormX;
						ret = currDens;
					}
					else if (  vertSelf->BndCond.GetBCTunnelTag() == FDBoundary::eTunnelOut )
					{
						double tunCoeff = GetPhysPrpty(TunnelCoeff);
						double dens = GetPhysPrpty(eDensity);
						ret = tunCoeff * dens * bcNormX;
					}
					else
					{
						ret = 0;
					}
				}
				break;
			}
			case eCurrentDensity_Y:
			{
				double bndNormY = 0;
				if ( vertSelf->IsAtBoundary(FDBoundary::eDensity) )
				{
					bndNormY = vertSelf->BndCond.GetBndDirection(FDBoundary::eDensity).Y();
				}
				//bnd_normY < 0, south boundary
				//bnd_normY > 0, north boundary
				//bnd_normY = 0, inner vertex
				if ( bndNormY == 0 )
				{
					double elecFieldY = GetPhysPrpty(ElectricField_Y);

					double hs = vertSelf->SouthLength;
					double hn = vertSelf->NorthLength;
					double ns = vertSelf->SouthVertex->Phys->GetPhysPrpty(eDensity);
					double nn = vertSelf->NorthVertex->Phys->GetPhysPrpty(eDensity);
					double nc = GetPhysPrpty(eDensity);
					SCTM_ASSERT(hs!=0 && hn!=0, 10024);
					double pn_div_py = ( -hn*hn * ns + (hn*hn - hs*hs) * nc + hs*hs * nn ) / ( hn*hs*(hn + hs) );
					double mobility = GetPhysPrpty(eMobility);

					ret = mobility * ( nc * elecFieldY - pn_div_py );
				}
				else
				{
					double bcNormY = vertSelf->BndCond.GetBCNormVector(FDBoundary::eDensity).Y();

					if ( vertSelf->BndCond.GetBCTunnelTag() == FDBoundary::eTunnelIn )
					{
						double bcValue = vertSelf->BndCond.GetBCValue(FDBoundary::eDensity);
						double currDens = bcValue * bcNormY;
						ret = currDens;
					}
					else if ( vertSelf->BndCond.GetBCTunnelTag() == FDBoundary::eTunnelOut )
					{
						double tunCoeff = GetPhysPrpty(TunnelCoeff);
						double dens = GetPhysPrpty(eDensity);
						ret = tunCoeff * dens * bcNormY;
					}
					else
					{
						ret = 0;
					}
				}
				break;
			}
			case eCurrentDensity:
			{
				double eCurrDensX = GetPhysPrpty(eCurrentDensity_X);
				double eCurrDensY = GetPhysPrpty(eCurrentDensity_Y);
				ret = SctmMath::sqrt( eCurrDensX * eCurrDensX + eCurrDensY * eCurrDensY );
				break;
			}
			case TunnelCoeff:
			{
				ret = tunnelCoeff;
				break;
			}
			default:
			{
				// use SCTM_ASSERT for non-existed property
				SCTM_ASSERT(SCTM_ERROR, 10001);
			}
		}

		return ret;
	}

	void PhysProperty::FillVertexPhysUsingMatPropty(PhysProperty::Name vertexPhys,
		MaterialDB::MatProperty::Name matPrpty)
	{
		//TODO : need to solve the problem of mutual including.
		//the problem is solved but with some unknowns.
		double tot = 0; //total area
		double sum = 0; //sum corresponds to the integral
		double physValue = 0; //the final physical value related to vertex

		using MaterialDB::GetMatPrpty;
		FDElement *currElem = NULL;
		currElem = vertSelf->NortheastElem;
		tot += ( currElem != NULL ) ? currElem->Area : 0;
		sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		currElem = vertSelf->NorthwestElem;
		tot += ( currElem != NULL ) ? currElem->Area : 0;
		sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		currElem = vertSelf->SoutheastElem;
		tot += ( currElem != NULL ) ? currElem->Area : 0;
		sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		currElem = vertSelf->SouthwestElem;
		tot += ( currElem != NULL ) ? currElem->Area : 0;
		sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		SCTM_ASSERT(tot>=0, 10004);
		physValue = sum / tot;

		vertSelf->Phys->SetPhysPrpty(vertexPhys, physValue);
	}		

	void PhysProperty::FillVertexPhysUsingMatPropty(PhysProperty::Name vertexPhys, 
		MaterialDB::MatProperty::Name matPrpty, FDRegion::RegionType rType)
	{
		//TODO : need to solve the problem of mutual including.
		//the problem is solved but with some unknowns.
		double tot = 0; //total area
		double sum = 0; //sum corresponds to the integral
		double physValue = 0; //the final physical value related to vertex

		using MaterialDB::GetMatPrpty;
		FDElement *currElem = NULL;
		currElem = vertSelf->NortheastElem;
		tot += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? currElem->Area : 0;
		sum += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		currElem = vertSelf->NorthwestElem;
		tot += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? currElem->Area : 0;
		sum += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		currElem = vertSelf->SoutheastElem;
		tot += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? currElem->Area : 0;
		sum += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		currElem = vertSelf->SouthwestElem;
		tot += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? currElem->Area : 0;
		sum += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		if (tot == 0)
		{
			physValue = 0;
		}
		else
		{
			physValue = sum / tot;
		}
		//SCTM_ASSERT(tot>=0, 10004);
		vertSelf->Phys->SetPhysPrpty(vertexPhys, physValue);
	}

	void PhysProperty::CalculateDensityControlArea()
	{
		double area = 0;
		FDElement *currElem = NULL;

		currElem = vertSelf->NorthwestElem;
		if ( (currElem != NULL) && (currElem->Region->Type == FDRegion::Trapping) )
		{
			area += 0.25 * currElem->Area;
		}
		currElem = vertSelf->NortheastElem;
		if ( (currElem != NULL) && (currElem->Region->Type == FDRegion::Trapping) )
		{
			area += 0.25 * currElem->Area;
		}
		currElem = vertSelf->SouthwestElem;
		if ( (currElem != NULL) && (currElem->Region->Type == FDRegion::Trapping) )
		{
			area += 0.25 * currElem->Area;
		}
		currElem = vertSelf->SoutheastElem;
		if ( (currElem != NULL) && (currElem->Region->Type == FDRegion::Trapping) )
		{
			area += 0.25 * currElem->Area;
		}
		//for vertex not related to trapping layer, the density control area is 0
		this->controlArea = area;
	}

	void PhysProperty::UpdateValue(Name prptyName, double val)
	{
		//only several properties are updated
		SetPhysPrpty(prptyName, val);
	}

	void SetPhysConstant()
	{
		//the value returned by GetMatPrpty is normalized value
		double mp = MaterialDB::GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_HoleMass);
		double mn = MaterialDB::GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_ElectronMass);
		double bandgap = MaterialDB::GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_Bandgap);
		double affinity = MaterialDB::GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_ElectronAffinity);
		double temp = T0;
		using SctmUtils::Normalization;
		Normalization norm = Normalization();
		//CAUTION: Is reference potential related to temperature
		//the value of reference potential is
		//phi(electron affinity) + Eg/2 + 3/4*kT/q*ln(mp/mn)
		//affinity and bandgap is in eV, so ReferencePotential is in [eV]
		//the reference potential should be normalized.
		SctmPhys::ReferencePotential = norm.PullEnergy(affinity) + norm.PullEnergy(bandgap / 2) - 3/4 * k0 * temp / q * SctmMath::ln( mp / mn );
		SctmPhys::ReferencePotential = norm.PushPotential(SctmPhys::ReferencePotential);
	}

}