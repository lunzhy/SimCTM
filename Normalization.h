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

namespace Utility
{
	const double EPSILON = 8.854188e-14;// in [F*cm-1]
	const double CHARGE = 1.602177e-19;// in [C]
	const double BOLTZMAN = 1.380662e-23;// in [J*K-1]
	const double INTRINSIC_CONC_SI = 1.0e10;// in [cm-3]
	const double ROOM_TEMP = 300;// in [K]
	const double RELATIVE_EPSILON_SI = 11.6;
	
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
		Normalization(double temperature);
		Normalization(void);//initialize the class without temperature will considerate the room temperature
		~Normalization(void);
	private:
		/// the factor of normalization and conversion
		double lengthFactor;///< factor of length normalization
		double potentialFactor;///< factor of potential normalization
		double elecFieldFactor;///< factor of electric field normalization
		double concFactor;///< factor of concentration normalization, including charge and carrier concentration
		double temperature;///< factor of the system temperature
	public:

		/// @brief ��һ���������
		/// 
		/// Push/Pull ���ڹ�һ�����������PushΪ���������ֵ��һ��Ϊ����ʹ��ֵ�� Pull������ʹ��ֵת��Ϊ��ֵ
		/// ���¶���������������Ĵ�����ͬ���ڴ�initial�ļ������ʱ�򣬶��ڶ����ֵ����һ��һ��д���������
		/// 
		/// @param double length
		/// @pre Ҫ�����¶�
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
		/// @brief ��ֵ�����Լ���һֵ�������໥ת�� 
		/// 
		/// ����������ʱ��������Ҫ����һֵ����ת����ֵ���������ʱ��ֱ�ӰѼ����һ����ת��Ϊ��ֵ����
		/// 
		/// @param std::vector<double> & real ʵ��ֵ����ֵ
		/// @param std::vector<double> & norm ��һֵ��������
		/// @param ConverseDirection direction ת������
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
		void initFactors();
	};
}