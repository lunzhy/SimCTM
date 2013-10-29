/**
* @file SimpleONO.h
* @brief This file deals with a simple finite differential domain problem, which is consisted of three layers.
*
*
*
* @author
* @version 
* @date 2013-7-5   11:15
* @note
* @todo
*/
#ifndef _SIMPLEONO_H_
#define _SIMPLEONO_H_
#include "FDDomain.h"
class FDRegion;

/// @brief SimpleONO is used for set the domain for finite differential method
///
/// This class is derived from SimpleONO. It can build a rectangle domain by using FDDomainHelper.
class SimpleONO : public FDDomain
{
	friend class DDTest;
	friend class DriftDiffusionSolver;
public:
	/// @brief BuildDomain is used to build the mesh of rectangle domain defined in this class
	/// 
	/// This method is called in outer class to build domain without pass-in parameters, this is
	/// because the parameters and details of the domain are defined and set in the class.
	/// In the mesh structure of this rectangle domain, the grid length in x direction is uniform
	/// And the grid length of y direction is uniform in each region.
	/// 
	/// @pre
	/// @return void
	/// @note
	void BuildDomain();
protected:
	double xLength; ///< the length in x direction of the rectangle structure.
	double yLengthTunnel; ///< the length in y direction of tunneling oxide
	double yLengthTrap; ///< the length in y direction of trapping layer
	double yLengthBlock; ///< the length in y direction of blocking layer

	int xCntVertex; ///< vertex count in x direction
	int yCntVertexTunnel; ///< vertex count in y direction of tunneling oxide
	int yCntVertexTrap; ///< vertex count in y direction of trapping layer
	int yCntVertexBlock; ///< vertex count in y direction of blocking layer
	int yCntTotalVertex; ///< total vertex count in y direction

	double xGrid; ///< the grid length in x direction
	double yGridTunnel; ///< the grid length in y direction of tunneling oxide
	double yGridTrap; ///< the grid length in y direction of trapping layer
	double yGridBlock; ///< the grid length in y direction of blocking layer

	double gatePotential; ///< the gate potential, related to FDContact Gate
	double channelPotential; ///< the channel potential 

private:
	/// @brief setParameters is used to set the structure parameters of the rectangle domain
	/// 
	/// Currently, the parameters used to build the mesh of this domain are defined and
	/// and can be revised here.
	/// 
	/// @pre
	/// @return void
	/// @note
	void setParameters();
	/// @brief setDomainDetails is used to define the domain details in the structure.
	/// 
	/// The regions, vertices and elements are defined in this method.
	/// 
	/// @pre
	/// @return void
	/// @note
	void setDomainDetails();
	/// @brief setAdjacency is used to set adjacent domain details for each vertex and element
	/// 
	/// Adjacent domain details are the vertices and elements near a specified vertex and element.
	/// It is used to link the vertices and elements together and build the data structure of the domain. 
	///
	/// @pre
	/// @return void
	/// @note
	void setAdjacency();
	/// @brief prepareStructure is used to prepare the domain structure for the following process.
	/// 
	/// Currently, setParameter is the only method called in this method. It is designed to include all
	/// the process used to prepare the simulation structure.
	/// 
	/// @pre
	/// @return void
	/// @note If current vertex/element doesn't have an adjacent neighbor in specific direction, the pointer is set NULL. 
	/// If current vertex doesn't have a west/east/south/north edge, the length is set 0.
	void prepareStructure();
	/// @brief printStructure is used to print the processing result of the designed structure.
	/// 
	/// This is a temporary method used to check the results.
	/// 
	/// @pre
	/// @return void
	/// @note
	void printStructure();
	
	/// @brief yNextGridLength is used to determine the next grid length in y direction
	/// 
	/// This is used in setting the coordinates when initializing the vertices at given vertex
	/// For the last vertex, this method returns 0.  The unit is [cm].
	/// 
	/// @param int vertexY, the current vertex index in y direction
	/// @pre
	/// @return double
	/// @note
	double yNextGridLength( int vertexY );
	/// @brief xNextGridLength is used to determine the next grid length in y direction
	/// 
	/// This is used in setting the coordinates when initializing the vertices at given vertex
	/// For the last vertex, this method returns 0. The unit is [cm].
	/// 
	/// @param int vertexX, the current vertex index in x direction
	/// @pre
	/// @return double
	/// @note
	double xNextGridLength( int vertexX );
	/// @brief thisRegion is used to determine the region with given element index in y direction.
	/// 
	/// This method is used to set the region sequence of the designed structure.
	/// Change the definition of this class to define new structure.
	/// 
	/// @param int elemY
	/// @pre
	/// @return FDRegion *
	/// @note the initial region sequence is tunneling oxide (region[0]), 
	/// trapping layer (region[1]) and blocking oxide (region[2]).
	FDRegion * thisRegion( int elemY);
	/// @brief stuffPotential is used to stuff potential obtained from calculation results
	/// 
	/// This is a temporary method for initializing the potential. The potential is calculated with
	/// channel potential and the electric field along the gate stack. These values are obtained from
	/// Sentaurus results. Here, we use same reference potential as in Sentaurus.
	/// 
	/// @pre
	/// @return void
	/// @note
	void stuffPotential();
	/// @brief setVertexPhysics is used to set the physical value related to the specific vertex.
	/// 
	/// These physical values are calculated from the material-related properties using the area of the adjacent elements 
	/// Electron affinity, bandgap, electron mass are calculated and filled into the physical value of the vertex.
	/// 
	/// @pre
	/// @return void
	/// @note
	void setVertexPhysics();
	/// @brief refreshBandEnergy is used to refresh band energy profile using the potential
	/// 
	/// The conduction/valence band energy is computed through the use of the obtained potential profile.
	/// 
	/// @pre
	/// @return void
	/// @note
	void refreshBandEnergy();
	/// @brief setBoundaryCondition is used to set the boundary condition for each vertex.
	/// 
	/// For each vertex, a first judgment on whether it is a vertex related to contact is made.
	/// And thus the boundary condition is set within this vertex. Secondly, if the contact is at
	/// the boundary, another type of boundary condition will be applied to this vertex.
	/// 
	/// @pre
	/// @return void
	/// @note
	void setBoundaryCondition();
	void setBoundaryCondition(bool fake);
	static bool isValidElem(FDElement *elem);
	static bool isNotElemOf(FDRegion::RegionType rType, FDElement *elem);
	void setVertBC_eDensity(FDVertex *vert);
	void setVertBC_Potential(FDVertex *vert);
};

#endif