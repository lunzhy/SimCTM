/**
* @file Normalization.h
* @brief This file deals with the normalization problems in the simulation.
*
*
*
* @author
* @version 
* @date 2013-7-9   9:38
* @note
* @todo
*/
#ifndef _NORMALIZATION_H_
#define _NORMALIZATION_H_
#include <vector>
#include <string>

/// @brief This namespace contains the common utilities used in the simulation
/// 
/// 
namespace SctmUtils
{	
	//the initialization of the const for normalization is not used.
	/*
	const double EPSILON = SctmPhys::eps / ( 1 / SctmPhys::cm_in_m); // in [F/cm]
	const double CHARGE = SctmPhys::ElementaryCharge; // in [C]
	const double BOLTZMAN = SctmPhys::BoltzmanConstant; // in [J/K]
	const double INTRINSIC_CONC_SI = 1.0e10;// in [cm-3]
	const double ROOM_TEMP = SctmPhys::RoomTemperature;// in [K]
	*/
	
	//the const double cannot be initialized here with other method.
	//const double RELATIVE_EPSILON_SI = MaterialDB::GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_DielectricConstant);
	
	/// @brief This class defines the normalization methods for the values used in the simulation
	///
	///
	class Normalization
	{
		friend class SctmTimeStep;
	public:
		/// @brief physical parameter normalization conversion direction
		///
		/// used in conversing the parameter vectors
		enum ConverseDirection
		{
			Push,///< normalize the read-in parameter
			Pull///< transfer the normalized parameter to real value
		};
		/// @brief Normalization a construction method for this class
		/// 
		/// This is used to initialize the normalization class with given temperature.
		/// 
		/// @param double temperature
		/// @pre
		/// @return 
		/// @note
		Normalization(double temperature);
		/// @brief Normalization is a construction method for this class
		/// 
		/// This method is used to initialize the class without temperature will considerate the room temperature.
		/// 
		/// @param void
		/// @pre
		/// @return 
		/// @note
		Normalization(void);
		~Normalization(void);
	private:
		/// the factor of normalization and conversion
		double lengthFactor;///< factor of length normalization
		double potentialFactor;///< factor of potential normalization
		double elecFieldFactor;///< factor of electric field normalization
		double densityFactor;///< factor of concentration/density normalization, including charge and carrier density
		double temperature;///< factor of the system temperature
		double diffusionFactor;///< factor of diffusion coefficient
		double mobilityFactor;///< factor of carrier mobility
		double timeFactor;///factor of time step
		double currDensFactor;///factor of current density
	public:

		/// @brief The method for normalization of the parameters.
		/// 
		/// Push/Pull methods are used for parameters normalization. Push is used to transform the real value into the 
		/// value that is normalized for simplifying the following calculation. Push is used to convert the values used
		/// in the simulation to real values.
		/// The following methods are designed for different parameters.
		/// 
		/// @param double length/...
		/// @pre temperature is predefined when the class is constructed.
		/// @return double
		/// @note
		//in [K]
		inline double PushLength(double length)
		{
			return length / lengthFactor;
		}
		inline double PullLength(double length)
		{
			return length * lengthFactor;
		}
		/// @brief This method is used to convert the parameter vector. 
		/// 
		/// By using this method, the normalized values and real values are converted according to the given
		/// conversion direction.
		/// The same is with the following similar method.
		/// 
		/// @param std::vector<double> & real real values
		/// @param std::vector<double> & norm normalized values used for the calculation
		/// @param ConverseDirection direction
		/// @pre
		/// @return void
		/// @note
		void ConverseLengthVector(std::vector<double> &real, std::vector<double> &norm, ConverseDirection direction);

		//potential, in [V]
		inline double PushPotential(double potential)
		{
			return potential / potentialFactor;
		}
		inline double PullPotential(double potential)
		{
			return potential * potentialFactor;
		}
		void ConversePotentialVector(std::vector<double> &real, std::vector<double> &norm, ConverseDirection direction);
		
		//in electric field, in [V*cm^-1]
		inline double PushElecField(double elecField)
		{
			return elecField / elecFieldFactor;
		}
		inline double PullElecField(double elecField)
		{
			return elecField * elecFieldFactor;
		}
		void ConverseElecFieldVector(std::vector<double> &real, std::vector<double> &norm, ConverseDirection direction);
		
		//concentration/density, in [cm^-3]
		inline double PushDensity(double density)
		{
			return density / densityFactor;
		}
		inline double PullDensity(double density)
		{
			return density * densityFactor;
		}
		void ConveseDensityVector(std::vector<double> &real, std::vector<double> &norm, ConverseDirection direction);
		
		//diffusion coefficient, in [cm^2/s]
		inline double PushDiffusion(double diffusion)
		{
			return diffusion / diffusionFactor;
		}
		inline double PullDiffusion(double diffusion)
		{
			return diffusion * diffusionFactor;
		}

		//mobility, in [cm^2/V/s]
		inline double PushMobility(double mobility)
		{
			return mobility / mobilityFactor;
		}
		inline double PullMobility(double mobility)
		{
			return mobility * mobilityFactor;
		}

		//time step, in [s]
		inline double PushTime(double time)
		{
			return time / timeFactor;
		}
		inline double PullTime(double time)
		{
			return time * timeFactor;
		}

		//current density, in [A/cm^2]
		inline double PushCurrDens(double currdens)
		{
			return currdens / currDensFactor;
		}
		inline double PullCurrDens(double currdens)
		{
			return currdens * currDensFactor;
		}
		
		//area, in [cm^2]
		inline double PullLineDensity(double lineDensity)
		{
			return lineDensity * densityFactor * lengthFactor * lengthFactor;
		}

		//energy is always used with q, in [eV]
		inline double PushEnergy(double energy)
		{
			return energy / potentialFactor;
		}
		inline double PullEnergy(double energy)
		{
			return energy * potentialFactor;
		}

		//tunneling probability, the coefficient for density to calculate current density.
		//J = P * n
		inline double PushTunCoeff(double coeff)
		{
			return coeff / (currDensFactor / densityFactor);
		}
		inline double PullTunCoeff(double coeff)
		{
			return coeff * (currDensFactor * densityFactor);
		}
	private:
		/// @brief initFactors is used to initialize the normalization factors.
		/// 
		///
		/// 
		/// @pre
		/// @return void
		/// @note
		void initFactors();
	};
}

#endif