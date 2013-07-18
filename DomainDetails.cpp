/**
* @file DomainDetails.cpp
* @brief This file contains the detailed implementation of the methods declared in the header file.
*
*
*
* @author
* @version 
* @date 2013-7-4   19:41
* @note
* @todo
*/

#include "DomainDetails.h"
#include "GeneralMath.h"

double FDVertex::Distance( FDVertex *vertex1, FDVertex *vertex2 )
{
	return GeneralMath::sqrt(GeneralMath::square(vertex1->X - vertex2->X) + GeneralMath::square(vertex1->Y - vertex2->Y));
}
