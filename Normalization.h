/**
* @file Normalization.h
* @brief
*
*
*
* @author
* @version 
* @date 2013-7-9   9:38
* @note
* @todo
*/
#pragma once
#include <vector>
#include <string>


/// @brief This namespace contains the common utilities used in the simulation
/// 
/// 
namespace SctmUtils
{
	const double EPSILON = 8.854188e-14;// in [F*cm-1]
	const double CHARGE = 1.602177e-19;// in [C]
	const double BOLTZMAN = 1.380662e-23;// in [J*K-1]
	const double INTRINSIC_CONC_SI = 1.0e10;// in [cm-3]
	const double ROOM_TEMP = 300;// in [K]
	const double RELATIVE_EPSILON_SI = 11.6;
	
	
	
	/// @brief This class defines the normalization methods for the values used in the simulation
	///
	///
	class Normalization
	{
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
		double concFactor;///< factor of concentration normalization, including charge and carrier concentration
		double temperature;///< factor of the system temperature
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

		//in [V]
		inline double PushPotential(double potential)
		{
			return potential / potentialFactor;
		}
		inline double PullPotential(double potential)
		{
			return potential * potentialFactor;
		}
		void ConversePotentialVector(std::vector<double> &real, std::vector<double> &norm, ConverseDirection direction);
		//in [V*cm-1]
		inline double PushElecField(double elecField)
		{
			return elecField / elecFieldFactor;
		}
		inline double PullElecField(double elecField)
		{
			return elecField * elecFieldFactor;
		}
		void ConverseElecFieldVector(std::vector<double> &real, std::vector<double> &norm, ConverseDirection direction);
		//in [cm-3]
		inline double PushConcentration(double conc)
		{
			return conc / concFactor;
		}
		inline double PullConcentration(double conc)
		{
			return conc * concFactor;
		}
		void ConveseConcVector(std::vector<double> &real, std::vector<double> &norm, ConverseDirection direction);
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