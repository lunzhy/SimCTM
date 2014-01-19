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
		//conductionBandEnergy = 0;
		//valenceBandEnergy = 0;
		electronAffinity = 0;
		e_mass = 0;
		h_mass = 0;
		//netCharge = 0;
		e_density = 0;
		e_mobility = 0;
		controlArea = 0;
		epsilon = 0;

		tunnelCoeff = 0;
		e_currdensMFN_X = 0;
		e_currdensMFN_Y = 0;
		e_subsCurrDensB2T = 0;

		//all the vectors and maps are initialized with 0 size
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
			case eEffDOS:
			{
				double mass = GetPhysPrpty(PhysProperty::eMass) * SctmPhys::m0;
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
				double mass = GetPhysPrpty(PhysProperty::eMass) * SctmPhys::m0;
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
			default:
			{
				// use SCTM_ASSERT for non-existed property
				SCTM_ASSERT(SCTM_ERROR, 10001);
			}
		}

		return ret;
	}

	void PhysProperty::FillVertexPhysUsingMatPrpty(PhysProperty::Name vertexPhys,
		MaterialDB::MatProperty::Name matPrpty)
	{
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

	void PhysProperty::FillVertexPhysUsingMatPrpty(PhysProperty::Name vertexPhys, 
		MaterialDB::MatProperty::Name matPrpty, FDRegion::TypeName rType)
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
						multiElectronAffinity[currMatName] = GetMatPrpty(MaterialMap(currMatName), matPrpty);
					}
					break;
				}
				case PhysProperty::Bandgap:
				{
					if (multiBandgap.find(currMatName) == multiBandgap.end())
					{
						multiBandgap[currMatName] = GetMatPrpty(MaterialMap(currMatName), matPrpty);
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
		double mp = MaterialDB::GetMatPrpty(MaterialMap(Mat::Silicon), MaterialDB::MatProperty::Mat_HoleMass);
		double mn = MaterialDB::GetMatPrpty(MaterialMap(Mat::Silicon), MaterialDB::MatProperty::Mat_ElectronMass);
		double bandgap = MaterialDB::GetMatPrpty(MaterialMap(Mat::Silicon), MaterialDB::MatProperty::Mat_Bandgap);
		double affinity = MaterialDB::GetMatPrpty(MaterialMap(Mat::Silicon), MaterialDB::MatProperty::Mat_ElectronAffinity);
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
		e_trapDensity = 0;
		e_trapped = 0;
		e_crossSection = 0;
		energyFromCondBand = 0;
		e_frequencyT2B = 0;
		e_transCoeffT2B = 0;
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
		case eCrossSection:
			e_crossSection = val;
			break;
		case EnergyFromCondBand:
			energyFromCondBand = val;
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
		case eFrequencyT2B:
			e_frequencyT2B = val;
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
				ret = e_trapped / e_trapDensity;
				break;
			}
			case EnergyFromCondBand:
			{
				ret = energyFromCondBand;
				break;
			}
			case eCrossSection:
			{
				ret = e_crossSection;
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
			case eEmissionCoeff_BasicSRH:
			{
				double eVelocity = 0;
				double eEffectiveDOS = 0;
				double temperature = vertSelf->Phys->GetPhysPrpty(PhysProperty::Temperature);

				double kT_div_q = SctmPhys::k0 * temperature / SctmPhys::q;
				eVelocity = vertSelf->Phys->GetPhysPrpty(PhysProperty::eThermalVelocity);
				eEffectiveDOS = vertSelf->Phys->GetPhysPrpty(PhysProperty::eEffDOS);

				using SctmUtils::Normalization;
				Normalization norm = Normalization(temperature);
				double trapEnergy = norm.PullEnergy(GetTrapPrpty(TrapProperty::EnergyFromCondBand));

				ret = GetTrapPrpty(eCrossSection) * eVelocity * eEffectiveDOS * 
					SctmMath::exp( - trapEnergy / kT_div_q);
				
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
			case NetCharge:
			{
				// there should be 2 parts here
				ret = - GetTrapPrpty(eTrapped);
				break;
			}
			case eTransCoeffT2B:
			{
				ret = e_transCoeffT2B;
				break;
			}
			case eFrequencyT2B:
			{
				ret = e_frequencyT2B;
				break;
			}
			case eEmissionCoeff_T2B:
			{
				double frequency = GetTrapPrpty(eFrequencyT2B);
				double transCoeff = GetTrapPrpty(eTransCoeffT2B);
				ret = frequency * transCoeff;
				break;
			}
			default:
			{
				SCTM_ASSERT(SCTM_ERROR, 10029);
			}
		}

		return ret;
	}

	double CalculateFlatbandShift(FDDomain *domain)
	{
		static FDContact *subsContact = domain->GetContact("Channel");
		//TODO: this is a temporary method used in one dimensional problems to get the start vertex in the calculation
		static FDVertex *startVert = subsContact->GetContactVerts().at(0);

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

		return VfbShift; // in normalized value
	}

}