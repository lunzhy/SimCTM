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

#include "SctmPhys.h"
#include "Material.h"
#include "SctmUtils.h"
#include "DomainDetails.h"
#include "FDDomain.h"
#include "Normalization.h"
#include "SctmMath.h"
#include <vector>
#include "TripleCells.h"

using SctmUtils::SctmData;

namespace SctmPhys
{
	const double &k0 = BoltzmanConstant;
	const double &h = PlanckConstant;
	const double &eps = VacuumDielectricConstant;
	const double &T0 = RoomTemperature;
	const double &q = ElementaryCharge;
	const double &ni = IntrinsicConcentration;
	const double &m0 = ElectronMass;
	
	double ReferencePotential;

	PhysProperty::PhysProperty(FDVertex *_vert)
	{
		temperature = SctmUtils::SctmGlobalControl::Get().Temperature;

		vertSelf = _vert;
		bandgap = 0;
		electrostaticPotential = 0;
		densControlArea = 0;
		epsilon = 0;
		
		electronAffinity = 0;
		e_DOSmass = 0;
		e_mass = 0;
		e_density = 0;
		e_mobility = 0;

		e_tunnelCoeff = 0;
		e_currdensMFN_X = 0;
		e_currdensMFN_Y = 0;
		e_subsCurrDensB2T = 0;

		e_emissionTime = 0;
		e_captureTime = 0;

		h_mass = 0;
		h_DOSmass = 0;
		h_mobility = 0;
		h_density = 0;
		h_tunnelCoeff = 0;

		//all the vectors and maps are initialized with 0 size
	}

	void PhysProperty::SetPhysPrpty(Name prptyName, double prptyValue)
	{
		switch (prptyName)
		{
		case ElectrostaticPotential:
			electrostaticPotential = prptyValue;
			break;
		case  eDOSMass:
			e_DOSmass = prptyValue;
			break;
		case eMass:
			e_mass = prptyValue;
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
		case eTunnelCoeff:
			e_tunnelCoeff = prptyValue;
			break;
		case DielectricConstant:
			epsilon = prptyValue;
			break;
		case eCurrDensMFN_X:
			e_currdensMFN_X = prptyValue;
			break;
		case eCurrDensMFN_Y:
			e_currdensMFN_Y = prptyValue;
			break;
		case eSubsCurrDensB2T:
			e_subsCurrDensB2T = prptyValue;
			break;
		case eEmissionTime:
			e_emissionTime = prptyValue;
			break;
		case eCaptureTime:
			e_captureTime = prptyValue;
			break;
		//for holes
		case hMass:
			h_mass = prptyValue;
			break;
		case hDOSMass:
			h_DOSmass = prptyValue;
			break;
		case hMobility:
			h_mobility = prptyValue;
			break;
		case hDensity:
			h_density = prptyValue;
			break;
		case hTunnelCoeff:
			h_tunnelCoeff = prptyValue;
			break;
		default:
			SCTM_ASSERT(SCTM_ERROR, 10019);
		}
	}

	double PhysProperty::GetPhysPrpty(Name prptyName, MaterialDB::Mat::Name matName /* = Mat::ErrorMaterial */) const
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
				//affinity = this->electronAffinity;
				affinity = GetPhysPrpty(ElectronAffinity, matName);
				energy = -affinity - (pot-RefPotential);
				ret = energy;
				break;
			}
			case ValenceBandEnergy:
			{
				//Ev = Ev-Eg = -X-q(phi-phiRef)-Eg
				
				ret =  GetPhysPrpty(ConductionBandEnergy, matName) - GetPhysPrpty(Bandgap, matName);
				break;
			}
			case ElectronAffinity:
			{
				if (matName == MaterialDB::Mat::ErrorMaterial)
				{
					ret = this->electronAffinity;
				}
				else
				{
					ret = getMultiPrptyValue(ElectronAffinity, matName);
				}
				//ret = electronAffinity;
				break;
			}
			case Bandgap:
			{
				if (matName == MaterialDB::Mat::ErrorMaterial)
				{
					ret = this->bandgap;
				}
				else
				{
					ret = getMultiPrptyValue(Bandgap, matName);
				}
				//ret = bandgap;
				break;
			}
			case NetFreeCharge:
			{
				//there should be four parts
				ret = -GetPhysPrpty(eDensity) + GetPhysPrpty(hDensity);
				break;
			}
			case eDOSMass:
			{
				ret = e_DOSmass;
				break;
			}
			case eMass:
			{
				ret = e_mass;
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
				ret = densControlArea;
				break;
			}
			case DielectricConstant:
			{
				ret = epsilon;
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
			case ElectricFieldTrap_X:
			{
				double bcNormX = 0; // inner vertex
				if (vertSelf->IsAtBoundary(FDBoundary::eDensity))
				{
					bcNormX = vertSelf->BndCond.GetBndDirection(FDBoundary::eDensity).X();
				}
				//bnd_normX < 0, west boundary
				//bnd_normX > 0, east boundary
				//bnd_normX = 0, inner vertex or vertex at boundary but inside in x direction
				if (bcNormX == 0)
				{
					// for boundary vertex in terms of X direction
					ret = GetPhysPrpty(PhysProperty::ElectricField_X);
				}
				else if (bcNormX < 0)
				{
					SCTM_ASSERT(vertSelf->EastVertex->EastVertex != NULL, 10025);
					double he = vertSelf->EastLength;
					double hee = vertSelf->EastVertex->EastLength;
					double fe = vertSelf->EastVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fee = vertSelf->EastVertex->EastVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fc = GetPhysPrpty(ElectrostaticPotential);
					SCTM_ASSERT(hee != 0 && he != 0, 10024);
					ret = -(-hee*(2 * he + hee) * fc + (he + hee)*(he + hee) * fe - he*he * fee) / (he*hee*(he + hee));
				}
				else //  (bndNormX > 0)
				{
					SCTM_ASSERT(vertSelf->WestVertex->WestVertex != NULL, 10025);
					double hw = vertSelf->WestLength;
					double hww = vertSelf->WestVertex->WestLength;
					double fw = vertSelf->WestVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fww = vertSelf->WestVertex->WestVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fc = GetPhysPrpty(ElectrostaticPotential);
					SCTM_ASSERT(hw != 0 && hww != 0, 10024);
					ret = -(hw*hw * fww - (hw + hww)*(hw + hww) * fw + hww*(2 * hw + hww) * fc) / (hw*hww*(hw + hww));
				}
				break;
			}
			case ElectricFieldTrap_Y:
			{
				double bcNormY = 0;
				if (vertSelf->IsAtBoundary(FDBoundary::eDensity))
				{
					bcNormY = vertSelf->BndCond.GetBCNormVector(FDBoundary::eDensity).Y();
				}
				//bnd_normY < 0, south boundary
				//bnd_normY > 0, north boundary
				//bnd_normY = 0, inner vertex
				if (bcNormY == 0)
				{
					ret = GetPhysPrpty(PhysProperty::ElectricField_Y);
				}
				else if (bcNormY < 0)
				{
					SCTM_ASSERT(vertSelf->NorthVertex->NorthVertex != NULL, 10025);
					double hn = vertSelf->NorthLength;
					double hnn = vertSelf->NorthVertex->NorthLength;
					double fn = vertSelf->NorthVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fnn = vertSelf->NorthVertex->NorthVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fc = GetPhysPrpty(ElectrostaticPotential);
					SCTM_ASSERT(hnn != 0 && hn != 0, 10024);
					ret = -(-hnn*(2 * hn + hnn) * fc + (hn + hnn)*(hn + hnn) * fn - hn*hn * fnn) / (hn*hnn*(hn + hnn));
				}
				else // ( bndNormX > 0 )
				{
					SCTM_ASSERT(vertSelf->SouthVertex->SouthVertex != NULL, 10025);
					double hs = vertSelf->SouthLength;
					double hss = vertSelf->SouthVertex->SouthLength;
					double fs = vertSelf->SouthVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fss = vertSelf->SouthVertex->SouthVertex->Phys->GetPhysPrpty(ElectrostaticPotential);
					double fc = GetPhysPrpty(ElectrostaticPotential);
					SCTM_ASSERT(hs != 0 && hss != 0, 10024);
					ret = -(hs*hs * fss - (hs + hss)*(hs + hss) * fs + hss*(2 * hs + hss) * fc) / (hs*hss*(hs + hss));
				}
				break;
			}
			case ElectricFieldTrap:
			{
				double elecFieldX = GetPhysPrpty(ElectricFieldTrap_X);
				double elecFieldY = GetPhysPrpty(ElectricFieldTrap_Y);
				ret = SctmMath::sqrt(elecFieldX * elecFieldX + elecFieldY * elecFieldY);
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

					// J = u (-n * p_phi/p_x + p_n / p_x )
					ret = mobility * ( nc * elecFieldX + pn_div_px);
				}
				else
				{
					double bcNormX = vertSelf->BndCond.GetBCNormVector(FDBoundary::eDensity).X();

					if ( vertSelf->BndCond.GetElecTunnelTag() == FDBoundary::eTunnelIn )
					{
						double bcValue = vertSelf->BndCond.GetBCValue(FDBoundary::eDensity);
						double currDens = bcValue * bcNormX;
						ret = currDens;
					}
					else if (  vertSelf->BndCond.GetElecTunnelTag() == FDBoundary::eTunnelOut )
					{
						double tunCoeff = GetPhysPrpty(eTunnelCoeff);
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

					ret = mobility * ( nc * elecFieldY + pn_div_py );
				}
				else
				{
					double bcNormY = vertSelf->BndCond.GetBCNormVector(FDBoundary::eDensity).Y();

					if ( vertSelf->BndCond.GetElecTunnelTag() == FDBoundary::eTunnelIn )
					{
						double bcValue = vertSelf->BndCond.GetBCValue(FDBoundary::eDensity);
						double currDens = bcValue * bcNormY;
						ret = currDens;
					}
					else if ( vertSelf->BndCond.GetElecTunnelTag() == FDBoundary::eTunnelOut )
					{
						double tunCoeff = GetPhysPrpty(eTunnelCoeff);
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
			case eTunnelCoeff:
			{
				ret = e_tunnelCoeff;
				break;
			}
			case eEffDOS:
			{
				double mass = GetPhysPrpty(PhysProperty::eDOSMass) * SctmPhys::m0;
				double kT = SctmPhys::k0 * GetPhysPrpty(PhysProperty::Temperature); // the simulation temperature
				double h = SctmPhys::h;
				double per_m3_in_per_cm3 = 1 / SctmMath::pow((1 / SctmPhys::cm_in_m), 3);

				// 3/2 = 1, 3.0/2 = 1.5
				ret = 2 * SctmMath::pow((2*SctmMath::PI*mass*kT), 3.0/2) / h / h / h; // in [m^-3]
				ret = ret * per_m3_in_per_cm3; // convert to [cm^-3]

				using SctmUtils::Normalization;
				Normalization norm = Normalization(GetPhysPrpty(PhysProperty::Temperature));
				ret = norm.PushDensity(ret);
				break;
			}
			case eThermalVelocity:
			{
				double mass = GetPhysPrpty(PhysProperty::eDOSMass) * SctmPhys::m0;
				double kT = SctmPhys::k0 * GetPhysPrpty(PhysProperty::Temperature);
				double m_in_cm = 1/ SctmPhys::cm_in_m;

				ret = SctmMath::sqrt(3 * kT / mass); // in [m/s]
				ret = ret * m_in_cm;
			
				using SctmUtils::Normalization;
				Normalization norm = Normalization(GetPhysPrpty(PhysProperty::Temperature));
				ret = norm.PushVelocity(ret);
				break;
			}
			case Temperature:
			{
				ret = temperature;
				break;
			}
			case eCurrDensMFN_X:
			{
				ret = e_currdensMFN_X;
				break;
			}
			case eCurrDensMFN_Y:
			{
				ret = e_currdensMFN_Y;
				break;
			}
			case eSubsCurrDensB2T:
			{
				ret = e_subsCurrDensB2T;
				break;
			}
			case eEmissionTime:
			{
				ret = e_emissionTime;
				break;
			}
			case eCaptureTime:
			{
				ret = e_captureTime;
				break;
			}
			case eOccupationTAT:
			{
				ret = GetPhysPrpty(eEmissionTime) / (GetPhysPrpty(eEmissionTime) + GetPhysPrpty(eCaptureTime));
				break;
			}
			//for holes
			case hMass:
			{
				ret = h_mass;
				break;
			}
			case hDOSMass:
			{
				ret = h_DOSmass;
				break;
			}
			case hMobility:
			{
				ret = h_mobility;
				break;
			}
			case hDensity:
			{
				ret = h_density;
				break;
			}
			case hTunnelCoeff:
			{
				ret = h_tunnelCoeff;
				break;
			}
			case hCurrentDensity_X:
			{
				double bndNormX = 0;
				//hDensity boundary condition is the same as eDensity boundary condition, except for the value
				if (vertSelf->IsAtBoundary(FDBoundary::eDensity))
				{
					bndNormX = vertSelf->BndCond.GetBndDirection(FDBoundary::eDensity).X();
				}
				//bnd_normX < 0, west boundary
				//bnd_normX > 0, east boundary
				//bnd_normX = 0, inner vertex or vertex at boundary but inside in x direction
				if (bndNormX == 0)
				{
					double elecFieldX = GetPhysPrpty(ElectricField_X);

					double hw = vertSelf->WestLength;
					double he = vertSelf->EastLength;
					double pw = vertSelf->WestVertex->Phys->GetPhysPrpty(hDensity);
					double pe = vertSelf->EastVertex->Phys->GetPhysPrpty(hDensity);
					double pc = GetPhysPrpty(hDensity);
					SCTM_ASSERT(hw != 0 && he != 0, 10024);
					double pp_div_px = (-he*he * pw + (he*he - hw*hw) * pc + hw*hw * pe) / (he*hw*(he + hw));
					double mobility = GetPhysPrpty(hMobility);

					// J = -u (p * p_phi/p_x + p_p / p_x )
					ret = -mobility * (-pc * elecFieldX + pp_div_px);
				}
				else
				{
					double bcNormX = vertSelf->BndCond.GetBCNormVector(FDBoundary::eDensity).X();

					if (vertSelf->BndCond.GetHoleTunnelTag() == FDBoundary::hTunnelIn)
					{
						double bcValue = vertSelf->BndCond.GetBCValue(FDBoundary::hDensity);
						double currDens = bcValue * bcNormX;
						ret = currDens;
					}
					else if (vertSelf->BndCond.GetHoleTunnelTag() == FDBoundary::hTunnelOut)
					{
						double tunCoeff = GetPhysPrpty(hTunnelCoeff);
						double dens = GetPhysPrpty(hDensity);
						ret = tunCoeff * dens * bcNormX;
					}
					else
					{
						ret = 0;
					}
				}
				break;
			}
			case hCurrentDensity_Y:
			{
				double bndNormY = 0;
				//hDensity boundary condition is the same as eDensity boundary condition, except for the value
				if (vertSelf->IsAtBoundary(FDBoundary::eDensity))
				{
					bndNormY = vertSelf->BndCond.GetBndDirection(FDBoundary::eDensity).Y();
				}
				//bnd_normY < 0, south boundary
				//bnd_normY > 0, north boundary
				//bnd_normY = 0, inner vertex
				if (bndNormY == 0)
				{
					double elecFieldY = GetPhysPrpty(ElectricField_Y);

					double hs = vertSelf->SouthLength;
					double hn = vertSelf->NorthLength;
					double ps = vertSelf->SouthVertex->Phys->GetPhysPrpty(hDensity);
					double pn = vertSelf->NorthVertex->Phys->GetPhysPrpty(hDensity);
					double pc = GetPhysPrpty(hDensity);
					SCTM_ASSERT(hs != 0 && hn != 0, 10024);
					double pp_div_py = (-hn*hn * ps + (hn*hn - hs*hs) * pc + hs*hs * pn) / (hn*hs*(hn + hs));
					double mobility = GetPhysPrpty(hMobility);

					// J = -u (p * p_phi/p_x + p_p / p_x )
					ret = -mobility * (-pc * elecFieldY + pp_div_py);
				}
				else
				{
					double bcNormY = vertSelf->BndCond.GetBCNormVector(FDBoundary::eDensity).Y();

					if (vertSelf->BndCond.GetHoleTunnelTag() == FDBoundary::hTunnelIn)
					{
						double bcValue = vertSelf->BndCond.GetBCValue(FDBoundary::hDensity);
						double currDens = bcValue * bcNormY;
						ret = currDens;
					}
					else if (vertSelf->BndCond.GetHoleTunnelTag() == FDBoundary::hTunnelOut)
					{
						double tunCoeff = GetPhysPrpty(hTunnelCoeff);
						double dens = GetPhysPrpty(hDensity);
						ret = tunCoeff * dens * bcNormY;
					}
					else
					{
						ret = 0;
					}
				}
				break;
			}
			case hCurrentDensity:
			{
				double hCurrDensX = GetPhysPrpty(hCurrentDensity_X);
				double hCurrDensY = GetPhysPrpty(hCurrentDensity_Y);
				ret = SctmMath::sqrt(hCurrDensX * hCurrDensX + hCurrDensY * hCurrDensY);
				break;
			}
			case hEffDOS:
			{
				double mass = GetPhysPrpty(PhysProperty::hDOSMass) * SctmPhys::m0;
				double kT = SctmPhys::k0 * GetPhysPrpty(PhysProperty::Temperature); // the simulation temperature
				double h = SctmPhys::h;
				double per_m3_in_per_cm3 = 1 / SctmMath::pow((1 / SctmPhys::cm_in_m), 3);

				// 3/2 = 1, 3.0/2 = 1.5
				ret = 2 * SctmMath::pow((2 * SctmMath::PI*mass*kT), 3.0 / 2) / h / h / h; // in [m^-3]
				ret = ret * per_m3_in_per_cm3; // convert to [cm^-3]

				using SctmUtils::Normalization;
				Normalization norm = Normalization(GetPhysPrpty(PhysProperty::Temperature));
				ret = norm.PushDensity(ret);
				break;
			}
			case hThermalVelocity:
			{
				double mass = GetPhysPrpty(PhysProperty::hDOSMass) * SctmPhys::m0;
				double kT = SctmPhys::k0 * GetPhysPrpty(PhysProperty::Temperature);
				double m_in_cm = 1 / SctmPhys::cm_in_m;

				ret = SctmMath::sqrt(3 * kT / mass); // in [m/s]
				ret = ret * m_in_cm;

				using SctmUtils::Normalization;
				Normalization norm = Normalization(GetPhysPrpty(PhysProperty::Temperature));
				ret = norm.PushVelocity(ret);
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

	void PhysProperty::FillVertexPhysUsingMatPrpty(PhysProperty::Name vertexPhys,
		MaterialDB::MatProperty::Name matPrpty, bool onlyTrapRegion/* = false */)
	{
		double tot = 0; //total area
		double sum = 0; //sum corresponds to the integral
		double physValue = 0; //the final physical value related to vertex
		double temp = 0;

		using MaterialDB::GetMatPrpty;
		FDElement *currElem = NULL;

		if (!onlyTrapRegion)
		{
			tot = 0;
			sum = 0;
			currElem = vertSelf->NortheastElem;
			tot += (currElem != NULL) ? currElem->Area : 0;
			sum += (currElem != NULL) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

			currElem = vertSelf->NorthwestElem;
			tot += (currElem != NULL) ? currElem->Area : 0;
			sum += (currElem != NULL) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

			currElem = vertSelf->SoutheastElem;
			tot += (currElem != NULL) ? currElem->Area : 0;
			sum += (currElem != NULL) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

			currElem = vertSelf->SouthwestElem;
			tot += (currElem != NULL) ? currElem->Area : 0;
			sum += (currElem != NULL) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;
		}
		else //only for trapping region
		{
			tot = 0;
			sum = 0;
			currElem = vertSelf->NortheastElem;
			if (!FDDomain::isNotTrappingElem(currElem))
			{
				temp = GetMatPrpty(currElem->Region->Mat, matPrpty);
				tot += currElem->Area;
				sum += currElem->Area * GetMatPrpty(currElem->Region->Mat, matPrpty);
			}

			currElem = vertSelf->NorthwestElem;
			if (!FDDomain::isNotTrappingElem(currElem))
			{
				temp = GetMatPrpty(currElem->Region->Mat, matPrpty);
				tot += currElem->Area;
				sum += currElem->Area * GetMatPrpty(currElem->Region->Mat, matPrpty);
			}

			currElem = vertSelf->SoutheastElem;
			if (!FDDomain::isNotTrappingElem(currElem))
			{
				temp = GetMatPrpty(currElem->Region->Mat, matPrpty);
				tot += currElem->Area;
				sum += currElem->Area * GetMatPrpty(currElem->Region->Mat, matPrpty);
			}

			currElem = vertSelf->SouthwestElem;
			if (!FDDomain::isNotTrappingElem(currElem))
			{
				temp = GetMatPrpty(currElem->Region->Mat, matPrpty);
				tot += currElem->Area;
				sum += currElem->Area * GetMatPrpty(currElem->Region->Mat, matPrpty);
			}
		}

		SCTM_ASSERT(tot >= 0, 10004);
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
		if (!FDDomain::isNotTrappingElem(currElem))
		{
			area += 0.25 * currElem->Area;
		}
		currElem = vertSelf->NortheastElem;
		if (!FDDomain::isNotTrappingElem(currElem))
		{
			area += 0.25 * currElem->Area;
		}
		currElem = vertSelf->SouthwestElem;
		if (!FDDomain::isNotTrappingElem(currElem))
		{
			area += 0.25 * currElem->Area;
		}
		currElem = vertSelf->SoutheastElem;
		if (!FDDomain::isNotTrappingElem(currElem))
		{
			area += 0.25 * currElem->Area;
		}
		//for vertex not related to trapping layer, the density control area is 0
		this->densControlArea = area;
	}

	void PhysProperty::UpdateValue(Name prptyName, double val)
	{
		//only several properties are updated
		SetPhysPrpty(prptyName, val);
	}

	void PhysProperty::SetMultiPrpty(PhysProperty::Name vertPhy, MaterialDB::MatProperty::Name matPrpty)
	{
		using namespace MaterialDB;
		FDElement *currElem = NULL;
		Mat::Name currMatName = Mat::ErrorMaterial;

		//load all the elements
		//Set south elements first. This order will determines the order of band data at vertex belonging to different material.
		//And it will affect data plotting in PySimFig
		vector<FDElement *> elems;
		currElem = vertSelf->SouthwestElem;
		if (currElem != NULL)
		{
			elems.push_back(currElem);
		}
		currElem = vertSelf->SoutheastElem;
		if (currElem != NULL)
		{
			elems.push_back(currElem);
		}
		currElem = vertSelf->NorthwestElem;
		if (currElem != NULL)
		{
			elems.push_back(currElem);
		}
		currElem = vertSelf->NortheastElem;
		if (currElem != NULL)
		{
			elems.push_back(currElem);
		}

		for (size_t iv = 0; iv != elems.size(); ++iv)
		{
			currMatName = elems.at(iv)->Region->Mat->MatName();

			//add the material name to related material name
			if (std::find(relatedMatName.begin(), relatedMatName.end(), currMatName) == relatedMatName.end())
			{
				this->relatedMatName.push_back(currMatName);
			}

			switch (vertPhy)
			{
				case PhysProperty::ElectronAffinity:
				{
					if (multiElectronAffinity.find(currMatName) == multiElectronAffinity.end())
					{
						multiElectronAffinity[currMatName] = GetMatPrpty(GetMaterial(currMatName), matPrpty);
					}
					break;
				}
				case PhysProperty::Bandgap:
				{
					if (multiBandgap.find(currMatName) == multiBandgap.end())
					{
						multiBandgap[currMatName] = GetMatPrpty(GetMaterial(currMatName), matPrpty);
					}
					break;
				}
				case PhysProperty::DielectricConstant:
				{
					if (multiDielectricConstant.find(currMatName) == multiDielectricConstant.end())
					{
						multiDielectricConstant[currMatName] = GetMatPrpty(GetMaterial(currMatName), matPrpty);
					}
					break;
				}
				default:
					SCTM_ASSERT(SCTM_ERROR, 10038);
					break;
			}
		}
	}

	double PhysProperty::getMultiPrptyValue(PhysProperty::Name verPhy, MaterialDB::Mat::Name matName) const
	{
		switch (verPhy)
		{
			case ElectronAffinity:
			{
				if (multiElectronAffinity.find(matName) != multiElectronAffinity.end())
				{
					return multiElectronAffinity.at(matName);
				}
				else
				{
					SCTM_ASSERT(SCTM_ERROR, 10039);
				}
				break;
			}
			case Bandgap:
			{
				if (multiBandgap.find(matName) != multiBandgap.end())
				{
					return multiBandgap.at(matName);
				}
				else
				{
					SCTM_ASSERT(SCTM_ERROR, 10039);
				}
				break;
			}
			default:
				SCTM_ASSERT(SCTM_ERROR, 10038);
				break;
		}
		return 0;
	}

	bool PhysProperty::HasMultiPrpty(PhysProperty::Name prptyName) const
	{
		switch (prptyName)
		{
			case ElectronAffinity:
			case ConductionBandEnergy:
			case ValenceBandEnergy:
			case Bandgap:
			{
				if (!(multiElectronAffinity.size() == multiBandgap.size() &&
					multiElectronAffinity.size() == relatedMatName.size()))
				{
					SCTM_ASSERT(SCTM_ERROR, 10040);
				}
				return (relatedMatName.size() != 1);
				break;
			}
			default:
				return false;
		}
	}

	std::vector<MaterialDB::Mat::Name> const& PhysProperty::GetRelatedMaterialNames()
	{
		return relatedMatName;
	}

	void SetPhysConstant()
	{
		using namespace MaterialDB;
		using SctmUtils::SctmGlobalControl;
		//the value returned by GetMatPrpty is normalized value
		double mp = MaterialDB::GetMatPrpty(GetMaterial(Mat::Silicon), MaterialDB::MatProperty::Mat_HoleDOSMass);
		double mn = MaterialDB::GetMatPrpty(GetMaterial(Mat::Silicon), MaterialDB::MatProperty::Mat_ElecDOSMass);
		double bandgap = MaterialDB::GetMatPrpty(GetMaterial(Mat::Silicon), MaterialDB::MatProperty::Mat_Bandgap);
		double affinity = MaterialDB::GetMatPrpty(GetMaterial(Mat::Silicon), MaterialDB::MatProperty::Mat_ElectronAffinity);
		double temperature = SctmGlobalControl::Get().Temperature;
		
		using SctmUtils::Normalization;
		Normalization norm = Normalization(SctmGlobalControl::Get().Temperature);
		//CAUTION: Is reference potential related to temperature ?
		//the value of reference potential is
		//phi(electron affinity) + Eg/2 + 3/4*kT/q*ln(mp/mn)
		//affinity and bandgap is in eV, so ReferencePotential is in [eV]
		//the reference potential should be normalized.
		SctmPhys::ReferencePotential = norm.PullEnergy(affinity) + norm.PullEnergy(bandgap / 2.0) - 0.75 * k0 * temperature / q * SctmMath::ln( mp / mn );
		SctmPhys::ReferencePotential = norm.PushPotential(SctmPhys::ReferencePotential);
	}


	TrapProperty::TrapProperty(FDVertex *_vert)
	{
		vertSelf = _vert;
		highFrqEpsilon = 0;
		trapDensity = 0;

		e_trapDensity = 0;
		e_trapped = 0;
		e_crossSection = 0;
		e_trapCrossSection = 0;
		e_energyFromCondBand = 0;
		e_frequencyT2B = 0;
		e_frequencyPF = 0;
		e_transCoeffT2B = 0;

		//for holes
		h_trapDensity = 0;
		h_trapped = 0;
		h_crossSection = 0;
		h_trapCrossSection = 0;
		h_frequencyT2B = 0;
		h_frequencyPF = 0;
		h_energyFromValeBand = 0;
	}

	void TrapProperty::FillTrapPrptyUsingMatPrpty(TrapProperty::Name trapPrpty, MaterialDB::MatProperty::Name matPrpty)
	{
		FDElement *currElem = NULL;
		double tot = 0; //total area
		double sum = 0; //sum corresponds to the integral
		double trapValue = 0; //the final trap property value related to vertex

		using MaterialDB::GetMatPrpty;
		double temp = 0;
		currElem = vertSelf->NorthwestElem;
		if ( !FDDomain::isNotTrappingElem(currElem) )
		{
			temp = GetMatPrpty(currElem->Region->Mat, matPrpty);
			tot += currElem->Area;
			sum += currElem->Area * GetMatPrpty(currElem->Region->Mat, matPrpty);
		}
		
		currElem = vertSelf->NortheastElem;
		if ( !FDDomain::isNotTrappingElem(currElem) )
		{
			temp = GetMatPrpty(currElem->Region->Mat, matPrpty);
			tot += currElem->Area;
			sum += currElem->Area * GetMatPrpty(currElem->Region->Mat, matPrpty);
		}
		
		currElem = vertSelf->SouthwestElem;
		if ( !FDDomain::isNotTrappingElem(currElem) )
		{
			temp = GetMatPrpty(currElem->Region->Mat, matPrpty);
			tot += currElem->Area;
			sum += currElem->Area * GetMatPrpty(currElem->Region->Mat, matPrpty);
		}
		
		currElem = vertSelf->SoutheastElem;
		if ( !FDDomain::isNotTrappingElem(currElem) )
		{
			temp = GetMatPrpty(currElem->Region->Mat, matPrpty);
			tot += currElem->Area;
			sum += currElem->Area * GetMatPrpty(currElem->Region->Mat, matPrpty);
		}

		if (tot == 0)
		{
			trapValue = 0;
		}
		else
		{
			trapValue = sum / tot;
		}
		vertSelf->Trap->SetTrapPrpty(trapPrpty, trapValue);
	}

	void TrapProperty::SetTrapPrpty(TrapProperty::Name trapPrpty, double val)
	{
		switch (trapPrpty)
		{
		case HighFrqEpsilon:
			highFrqEpsilon = val;
			break;
		case TrapDensity:
			trapDensity = val;
			break;
		case eCrossSection:
			e_crossSection = val;
			break;
		case eTrapCrossSection:
			e_trapCrossSection = val;
			break;
		case eEnergyFromCondBand:
			e_energyFromCondBand = val;
			break;
		case eTrapDensity:
			e_trapDensity = val;
			break;
		case eTrapped:
			e_trapped = val;
			break;
		case eTransCoeffT2B:
			e_transCoeffT2B = val;
			break;
		case eFrequency_T2B:
			e_frequencyT2B = val;
			break;
		case eFrequency_PF:
			e_frequencyPF = val;
			break;
		//for holes
		case hTrapped:
			h_trapped = val;
			break;
		case hTrapDensity:
			h_trapDensity = val;
			break;
		case hCrossSection:
			h_crossSection = val;
			break;
		case hTrapCrossSection:
			h_trapCrossSection = val;
			break;
		case hEnergyFromValeBand:
			h_energyFromValeBand = val;
			break;
		case hFrequency_T2B:
			h_frequencyT2B = val;
			break;
		case hFrequency_PF:
			h_frequencyPF = val;
			break; 
		default:
			SCTM_ASSERT(SCTM_ERROR, 10028);
		}
	}

	double TrapProperty::GetTrapPrpty(TrapProperty::Name trapPrpty) const
	{
		double ret = 0;
		switch (trapPrpty)
		{
			case HighFrqEpsilon:
			{
				ret = highFrqEpsilon;
				break;
			}
			case NetTrappedCharge:
			{
				// there should be 2 parts here
				ret = -GetTrapPrpty(eTrapped) + GetTrapPrpty(hTrapped);
				break;
			}
			case TrapDensity:
			{
				ret = trapDensity;
				break;
			}
			case eTrapped:
			{
				ret = e_trapped;
				break;
			}
			case eTrapDensity:
			{
				ret = e_trapDensity;
				break;
			}
			case eOccupation:
			{
				//ret = e_trapped / e_trapDensity;
				ret = e_trapped / GetTrapPrpty(TrapProperty::TrapDensity);
				break;
			}
			case eEnergyFromCondBand:
			{
				ret = e_energyFromCondBand;
				break;
			}
			case eCrossSection:
			{
				ret = e_crossSection;
				break;
			}
			case eTrapCrossSection:
			{
				ret = e_trapCrossSection;
				break;
			}
			case eEmptyTrapDens:
			{
				ret = e_trapDensity - e_trapped;
				break;
			}
			case eCaptureCoeff_J_Model:
			{
				double mobility = 0;
				double elecField = 0;
				mobility = vertSelf->Phys->GetPhysPrpty(PhysProperty::eMobility);
				elecField = vertSelf->Phys->GetPhysPrpty(PhysProperty::ElectricField);
				ret = GetTrapPrpty(eCrossSection) * mobility * elecField;
				break;
			}
			case eCaptureCoeff_V_Model:
			{
				double eVelocity = vertSelf->Phys->GetPhysPrpty(PhysProperty::eThermalVelocity);
				ret = eVelocity * GetTrapPrpty(eCrossSection);
				break;
			}
			case eTrappedCapCoeff_V_Model:
			{
				// this coefficient describes trapped electrons capturing holes
				double eVelocity = vertSelf->Phys->GetPhysPrpty(PhysProperty::hThermalVelocity);
				ret = eVelocity * GetTrapPrpty(eTrapCrossSection);
				break;
			}
			case eTrappedCapCoeff_J_Model:
			{
				double mobility = 0;
				double elecField = 0;
				mobility = vertSelf->Phys->GetPhysPrpty(PhysProperty::hMobility);
				elecField = vertSelf->Phys->GetPhysPrpty(PhysProperty::ElectricField);
				ret = GetTrapPrpty(eTrapCrossSection) * mobility * elecField;
				break;
			}
			case eEmissionCoeff_BasicSRH:
			{
				double eVelocity = 0;
				double eEffectiveDOS = 0;

				eVelocity = vertSelf->Phys->GetPhysPrpty(PhysProperty::eThermalVelocity);
				eEffectiveDOS = vertSelf->Phys->GetPhysPrpty(PhysProperty::eEffDOS);

				double trapDepth = GetTrapPrpty(TrapProperty::eEnergyFromCondBand);

				using SctmUtils::SctmGlobalControl;
				static string pfModel = SctmGlobalControl::Get().PhysicsPFModel;
				
				if (pfModel == "None" || pfModel == "Frequency") //ordinary temperature-based emission
				{
					ret = GetTrapPrpty(eCrossSection) * eVelocity * eEffectiveDOS *
						SctmMath::exp(-trapDepth); // kT/q will disappear with normalized energy
				}
				else if (pfModel == "EtDecrease")
				{
					double pfDecrease = GetTrapPrpty(TrapProperty::eTrapEnergyDecreasePF);

					using SctmUtils::Normalization;
					Normalization norm = Normalization(SctmGlobalControl::Get().Temperature);
					double freq = norm.PullFrequency(GetTrapPrpty(eCrossSection) * eVelocity * eEffectiveDOS);
					
					ret = GetTrapPrpty(eCrossSection) * eVelocity * eEffectiveDOS *
						SctmMath::exp(-(trapDepth - pfDecrease)); // kT/q will disappear with normalized energy
					//TODO: warning when pfDecrease > trapEnergy
					if (pfDecrease > trapDepth)
					{
						SctmData::Get().WritePooleFrenkelInfo();
					}
				}
				else
				{
					SCTM_ASSERT(SCTM_ERROR, 10045);
				}
				
				break;
			}
			case eCoeff_B2T:
			{
				double q = SctmPhys::q;
				double eCurrDens = vertSelf->Phys->GetPhysPrpty(PhysProperty::eSubsCurrDensB2T);
				double eXsecion = GetTrapPrpty(eCrossSection);
 
				ret =  eCurrDens * eXsecion;
				break;
			}
			case eTransCoeffT2B:
			{
				ret = e_transCoeffT2B;
				break;
			}
			case eFrequency_T2B:
			{
				ret = e_frequencyT2B;
				break;
			}
			case eEmissionCoeff_T2B:
			{
				double frequency = GetTrapPrpty(eFrequency_T2B);
				double transCoeff = GetTrapPrpty(eTransCoeffT2B);
				ret = frequency * transCoeff;
				break;
			}
			case eFrequency_PF:
			{
				ret = e_frequencyPF;
				break;
			}
			case eTrapEnergyDecreasePF:
			{
				double temperature = vertSelf->Phys->GetPhysPrpty(PhysProperty::Temperature);
				using SctmUtils::Normalization;
				Normalization norm = Normalization(temperature);
				double q = SctmPhys::q;

				double elecField = norm.PullElecField(vertSelf->Phys->GetPhysPrpty(PhysProperty::ElectricFieldTrap)); // in [V/cm], real value
				double eps = GetTrapPrpty(TrapProperty::HighFrqEpsilon) *
					SctmPhys::VacuumDielectricConstant / (1 / SctmPhys::cm_in_m); // in [F/cm]


				double deltaEt = SctmMath::sqrt(q * q * q / SctmMath::PI / eps) * SctmMath::sqrt(elecField) / q; // in [eV], real value.
				return norm.PushEnergy(deltaEt);
			}
			case eEmissionCoeff_PF:
			{
				double deltaEt = GetTrapPrpty(TrapProperty::eTrapEnergyDecreasePF); // in [eV], normalized value
				double frequency = GetTrapPrpty(TrapProperty::eFrequency_PF);
				double trapDepth = GetTrapPrpty(TrapProperty::eEnergyFromCondBand);

				double coeff = frequency * SctmMath::exp(-(trapDepth - deltaEt));
				if (deltaEt > trapDepth)
				{
					SctmData::Get().WritePooleFrenkelInfo();
				}
				return coeff;
			}
			//for holes
			case hTrapped:
			{
				ret = h_trapped;
				break;
			}
			case hTrapDensity:
			{
				ret = h_trapDensity;
				break;
			}
			case hOccupation:
			{
				//ret = h_trapped / h_trapDensity;
				ret = h_trapped / GetTrapPrpty(TrapProperty::TrapDensity);
				break;
			}
			case hEmptyTrapDens:
			{
				ret = h_trapDensity - h_trapped;
				break;
			}
			case hEnergyFromValeBand:
			{
				ret = h_energyFromValeBand;
				break;
			}
			case hCrossSection:
			{
				ret = h_crossSection;
				break;
			}
			case hTrapCrossSection:
			{
				ret = h_trapCrossSection;
				break;
			}
			case hFrequency_T2B:
			{
				ret = h_frequencyT2B;
				break;
			}
			case hFrequency_PF:
			{
				ret = h_frequencyPF;
				break;
			}
			case hCaptureCoeff_J_Model:
			{
				double mobility = 0;
				double elecField = 0;
				mobility = vertSelf->Phys->GetPhysPrpty(PhysProperty::hMobility);
				elecField = vertSelf->Phys->GetPhysPrpty(PhysProperty::ElectricField);
				ret = GetTrapPrpty(hCrossSection) * mobility * elecField;
				break;
			}
			case hCaptureCoeff_V_Model:
			{
				double eVelocity = vertSelf->Phys->GetPhysPrpty(PhysProperty::hThermalVelocity);
				ret = eVelocity * GetTrapPrpty(hCrossSection);
				break;
			}
			case hTrappedCapCoeff_V_Model:
			{
				// this coefficient describes trapped holes capturing free electrons
				double eVelocity = vertSelf->Phys->GetPhysPrpty(PhysProperty::eThermalVelocity);
				ret = eVelocity * GetTrapPrpty(hTrapCrossSection);
				break;
			}
			case hTrappedCapCoeff_J_Model:
			{
				double mobility = 0;
				double elecField = 0;
				mobility = vertSelf->Phys->GetPhysPrpty(PhysProperty::eMobility); // Jn * sigma_h
				elecField = vertSelf->Phys->GetPhysPrpty(PhysProperty::ElectricField);
				ret = GetTrapPrpty(hTrapCrossSection) * mobility * elecField;
				break;
			}
			case hEmissionCoeff_BasicSRH:
			{
				double hVelocity = 0;
				double hEffectiveDOS = 0;

				hVelocity = vertSelf->Phys->GetPhysPrpty(PhysProperty::hThermalVelocity);
				hEffectiveDOS = vertSelf->Phys->GetPhysPrpty(PhysProperty::hEffDOS);

				double trapDepth = GetTrapPrpty(TrapProperty::hEnergyFromValeBand);

				using SctmUtils::SctmGlobalControl;
				static string pfModel = SctmGlobalControl::Get().PhysicsPFModel;

				if (pfModel == "None" || pfModel == "Frequency") //ordinary temperature-based emission
				{
					ret = GetTrapPrpty(hCrossSection) * hVelocity * hEffectiveDOS *
						SctmMath::exp(-trapDepth); // kT/q will disappear with normalized energy
				}
				else if (pfModel == "EtDecrease")
				{
					double pfDecrease = GetTrapPrpty(TrapProperty::hTrapEnergyDecreasePF);
					ret = GetTrapPrpty(hCrossSection) * hVelocity * hEffectiveDOS *
						SctmMath::exp(-(trapDepth - pfDecrease)); // kT/q will disappear with normalized energy
					//TODO: warning when pfDecrease > trapEnergy
					if (pfDecrease > trapDepth)
					{
						SctmData::Get().WritePooleFrenkelInfo();
					}
				}
				else
				{
					SCTM_ASSERT(SCTM_ERROR, 10045);
				}

				break;
			}
			case hTrapEnergyDecreasePF:
			{
				double temperature = vertSelf->Phys->GetPhysPrpty(PhysProperty::Temperature);
				using SctmUtils::Normalization;
				Normalization norm = Normalization(temperature);
				double q = SctmPhys::q;

				double elecField = norm.PullElecField(vertSelf->Phys->GetPhysPrpty(PhysProperty::ElectricFieldTrap)); // in [V/cm], real value
				double eps = GetTrapPrpty(TrapProperty::HighFrqEpsilon) *
					SctmPhys::VacuumDielectricConstant / (1 / SctmPhys::cm_in_m); // in [F/cm]


				double deltaEt = SctmMath::sqrt(q * q * q / SctmMath::PI / eps) * SctmMath::sqrt(elecField) / q; // in [eV], real value.
				return norm.PushEnergy(deltaEt);
			}
			case hEmissionCoeff_PF:
			{
				double deltaEt = GetTrapPrpty(TrapProperty::hTrapEnergyDecreasePF); // in [eV], normalized value
				double frequency = GetTrapPrpty(TrapProperty::hFrequency_PF);
				double trapDepth = GetTrapPrpty(TrapProperty::hEnergyFromValeBand);

				double coeff = frequency * SctmMath::exp(-(trapDepth - deltaEt));
				if (deltaEt > trapDepth)
				{
					SctmData::Get().WritePooleFrenkelInfo();
				}
				return coeff;
			}
			default:
			{
				SCTM_ASSERT(SCTM_ERROR, 10029);
			}
		}

		return ret;
	}

	double CalculateFlatbandShift_domain(FDDomain *domain)
	{
		//TODO: this was a temporary method used in one dimensional problems to get the start vertex in the calculation
		//now the domain flatband voltage shift calculation method is 
		//to compute the weighted-average flat voltage shift of each slice with respect to its control length.
		//However, this method is still a structure-dependent method. The pre-knowledge of the structure has to been known.
		static FDContact *subsContact = domain->GetContact("Channel");
		static vector<FDVertex *> &channelVerts = subsContact->GetContactVerts();

		double VfbShift_slice = 0;
		double controlLength = 0;
		double sum = 0;
		double length = 0;
		FDVertex *vert = NULL;

		FDVertex *eastVert = NULL;
		FDVertex *westVert = NULL;
		for (size_t iVert = 0; iVert != channelVerts.size(); ++iVert)
		{
			vert = channelVerts.at(iVert);
			VfbShift_slice = CalculateFlatbandShift_slice_for1D(vert);
			
			eastVert = vert->EastVertex;
			westVert = vert->WestVertex;
			
			controlLength = 0;
			if (eastVert != NULL)
			{
				if (std::find(channelVerts.begin(), channelVerts.end(), eastVert) != channelVerts.end())
				{
					controlLength += vert->EastLength / 2;
				}
			}
			if (westVert != NULL)
			{
				if (std::find(channelVerts.begin(), channelVerts.end(), westVert) != channelVerts.end())
				{
					controlLength += vert->WestLength / 2;
				}
			}
			length += controlLength;

			sum += VfbShift_slice * controlLength;
		}
		
		return sum / length; // in normalized value
	}

	double CalculateFlatbandShift_slice_for1D(FDVertex *channelVert)
	{
		FDVertex *startVert = channelVert;

		using namespace SctmUtils;
		Normalization norm = Normalization(SctmGlobalControl::Get().Temperature);

		FDVertex *currVert = NULL;
		FDVertex *vertForCap = NULL; // the vertex pointer for calculation of capacitance
		double densCtrlArea = 0;

		double eFreeDens = 0;
		double eTrappedDens = 0;
		double eLineDens = 0;
		double hFreeDens = 0;
		double hTrappedDens = 0;
		double hLineDens = 0;

		double freeLineCharge = 0;
		double trapLineCharge = 0;

		double cap_reciprocal = 0;
		double epsilon = 0;
		double wide = 0;
		double delta_d = 0;
		double VfbShift = 0;
		
		//begin to calculate flatband voltage shift in the specific slice.
		currVert = startVert;
		while (currVert != NULL)
		{
			//for vertex in the trapping layer
			if (currVert->Trap != NULL)
			{
				densCtrlArea = currVert->Phys->GetPhysPrpty(PhysProperty::DensityControlArea);

				freeLineCharge = densCtrlArea * currVert->Phys->GetPhysPrpty(PhysProperty::NetFreeCharge);
				trapLineCharge = densCtrlArea * currVert->Trap->GetTrapPrpty(TrapProperty::NetTrappedCharge);

// 				eFreeDens = currVert->Phys->GetPhysPrpty(PhysProperty::eDensity);
// 				eTrappedDens = currVert->Trap->GetTrapPrpty(TrapProperty::eTrapped);
// 				eLineDens = densCtrlArea * (eFreeDens + eTrappedDens);
// 				hFreeDens = currVert->Phys->GetPhysPrpty(PhysProperty::hDensity);
// 				hTrappedDens = currVert->Trap->GetTrapPrpty(TrapProperty::hTrapped);
// 				hLineDens = densCtrlArea * (hFreeDens + hTrappedDens);

				vertForCap = currVert;
				cap_reciprocal = 0; // the reciprocal of capacitance

				while (vertForCap != NULL)
				{
					wide = (vertForCap->EastLength + vertForCap->WestLength) / 2;

					if (SctmGlobalControl::Get().Coordinate == "Cylindrical")
					{
						delta_d = vertForCap->R * SctmMath::ln((vertForCap->R + vertForCap->NorthLength / 2) / 
							(vertForCap->R - vertForCap->SouthLength / 2));
					}
					else // SctmGlobalControl::Get().Coordinate == "Cartesian"
					{
						delta_d = (vertForCap->SouthLength + vertForCap->NorthLength) / 2;
					}

					epsilon = vertForCap->Phys->GetPhysPrpty(PhysProperty::DielectricConstant);
					cap_reciprocal += delta_d / epsilon / wide;

					vertForCap = vertForCap->NorthVertex;
				}
				VfbShift += -(freeLineCharge + trapLineCharge) * cap_reciprocal;
			}
			currVert = currVert->NorthVertex;
		}

		return VfbShift;
	}

	double CalculateFlatbandShift_slice_for2D(FDVertex *channelVert)
	{
		FDVertex *startVert = channelVert;

		using namespace SctmUtils;
		Normalization norm = Normalization(SctmGlobalControl::Get().Temperature);

		FDVertex *currVert = NULL;
		FDVertex *vertForCap = NULL; // the vertex pointer for calculation of capacitance
		
		double densCtrlArea = 0;
		double eFreeDens = 0;
		double eTrappedDens = 0;
		double eLineDens = 0;
		double cap_reciprocal = 0;
		double epsilon = 0;
		double wide = 0;
		double delta_d = 0;
		double VfbShift = 0;

		//begin to calculate flatband voltage shift in the specific slice.
		currVert = startVert;
		while (currVert != NULL)
		{
			//for vertex in the trapping layer
			if (currVert->Trap != NULL)
			{
				densCtrlArea = currVert->Phys->GetPhysPrpty(PhysProperty::DensityControlArea);
				eFreeDens = currVert->Phys->GetPhysPrpty(PhysProperty::eDensity);
				eTrappedDens = currVert->Trap->GetTrapPrpty(TrapProperty::eTrapped);
				eLineDens = densCtrlArea * (eFreeDens + eTrappedDens);

				vertForCap = currVert;
				cap_reciprocal = 0; // the reciprocal of capacitance

				while (vertForCap != NULL)
				{
					if (SctmGlobalControl::Get().Structure == "Triple")
					{
						if (TripleCells::IsEndOfEffectiveCapacitor(vertForCap))
						{
							break;
						}
					}
					else if (SctmGlobalControl::Get().Structure == "TripleFull")
					{
						if (TripleCellsFull::IsEndOfEffectiveCapacitor(vertForCap))
						{
							break;
						}
					}
					wide = (vertForCap->EastLength + vertForCap->WestLength) / 2;
					delta_d = (vertForCap->SouthLength + vertForCap->NorthLength) / 2;
					epsilon = vertForCap->Phys->GetPhysPrpty(PhysProperty::DielectricConstant);
					cap_reciprocal += delta_d / epsilon / wide;

					vertForCap = vertForCap->NorthVertex;
				}
				VfbShift += eLineDens * cap_reciprocal;
			}
			currVert = currVert->NorthVertex;
		}

		return VfbShift;
	}

}
