/*
  Copyright (C) 2012,2013,2017(H)
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "log4espp.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "NodeGrid.hpp"

#include <boost/bind.hpp>
#include "mpi.hpp"

namespace espressopp {
  namespace storage {
    LOG4ESPP_LOGGER(NodeGrid::logger, "DomainDecomposition.NodeGrid");
  
    NodeGridIllegal::NodeGridIllegal()
      : std::invalid_argument("node grid dimensions have to be positive") {}

    NodeGrid::
    NodeGrid(const Int3D& grid,
         const longint nodeId,
         const Real3D& domainSize,
         const boost::python::list& _neiListx,
         const boost::python::list& _neiListy,
         const boost::python::list& _neiListz)
      : Grid(grid)
    {
      if (grid[0] <= 0 || grid[1] <= 0 || grid[2] <= 0) {
        throw NodeGridIllegal();
      }
      for (int i = 0; i < boost::python::len(_neiListx); ++i) {
        neiListx.push_back(boost::python::extract<int>(_neiListx[i]));
      }
      for (int i = 0; i < boost::python::len(_neiListy); ++i) {
        neiListy.push_back(boost::python::extract<int>(_neiListy[i]));
      }
      for (int i = 0; i < boost::python::len(_neiListz); ++i) {
        neiListz.push_back(boost::python::extract<int>(_neiListz[i]));
      }

      init(nodeId, domainSize);
    }

    NodeGrid::
    NodeGrid(const Int3D& grid,
         const longint nodeId,
         const Real3D& domainSize,
         const std::vector<int>& _neiListx,
         const std::vector<int>& _neiListy,
         const std::vector<int>& _neiListz)
      : Grid(grid)
    {
      if (grid[0] <= 0 || grid[1] <= 0 || grid[2] <= 0) {
        throw NodeGridIllegal();
      }
      // ################################   H   ############################################
      // Domain size static changes individual Neighbor Lists X, Y and Z
      // ###################################################################################
      neiListx = _neiListx;
      neiListy = _neiListy;
      neiListz = _neiListz;

      init(nodeId, domainSize);
    }


    void NodeGrid::init(const longint nodeId,
                        const Real3D& domainSize)
    {
    maxDomainSizeInCells[0]=neiListx[neiListx.size()-1];
    maxDomainSizeInCells[1]=neiListy[neiListy.size()-1];
    maxDomainSizeInCells[2]=neiListz[neiListz.size()-1];


    mapIndexToPosition(localNodePos, mpiWorld->rank());  // check this  |  mapIndexToPosition(nodePos, node);
    localBoxSize[0] = (static_cast<real>((neiListx[localNodePos[0]+1])-(neiListx[localNodePos[0]]))/static_cast<real>(maxDomainSizeInCells[0]))*domainSize[0];
    invLocalBoxSize[0] = 1.0/localBoxSize[0];

    localBoxSize[1] = (static_cast<real>((neiListy[localNodePos[1]+1])-(neiListy[localNodePos[1]]))/static_cast<real>(maxDomainSizeInCells[1]))*domainSize[1];
    invLocalBoxSize[1] = 1.0/localBoxSize[1];

    localBoxSize[2] = (static_cast<real>((neiListz[localNodePos[2]+1])-(neiListz[localNodePos[2]]))/static_cast<real>(maxDomainSizeInCells[2]))*domainSize[2];
    invLocalBoxSize[2] = 1.0/localBoxSize[2];

  smallestLocalBoxDiameter = std::min(std::min(localBoxSize[0], localBoxSize[1]), localBoxSize[2]);

    std::cout << "NODE GRID INIT, id " << nodeId << ", localBoxSize " << localBoxSize << ", LEFT: " << getMyLeft() << ", RIGHT " << getMyRight() <<  ", maxdomainsizeincells " << maxDomainSizeInCells[0] << " " << maxDomainSizeInCells[1] << " " << maxDomainSizeInCells[2] << std::endl;

  calcNodeNeighbors(nodeId);
    }


    // ################################   H   ############################################
    // Check for Changes!!!   mapPositionToNodeClipped
    // ###################################################################################
    longint NodeGrid::
    mapPositionToNodeClipped(const Real3D& pos) const
    {
      Int3D cpos, npos;
      real kCellGridSize[3];
      //kCellGridSize=localBoxSize[0]/(static_cast<real>((neiListx[localNodePos[0]+1])-(neiListx[localNodePos[0]]))/static_cast<real>(maxDomainSizeInCells[0]))/static_cast<real>(maxDomainSizeInCells[0]);
      
      std::vector<std::vector<int> > neiList = {neiListx, neiListy, neiListz};

      for (int dir = 0; dir < 3; ++dir) {
          kCellGridSize[dir]=localBoxSize[dir]/(static_cast<real>(neiList[dir][localNodePos[dir]+1])-(neiList[dir][localNodePos[dir]]));
          cpos[dir] = static_cast< int >(pos[dir]/kCellGridSize[dir]);
        if (cpos[dir] <= 0) {
          npos[dir] = 0;
        }
        //else if (cpos[dir] >= getGridSize(dir)) {
        else if (cpos[dir] >= maxDomainSizeInCells[dir]) {
          npos[dir] = getGridSize(dir) - 1;
        }
        else {
          for (int j=0;j<neiList[dir].size(); ++j) {
            if (cpos[dir]>=(neiList[dir][j]) && cpos[dir]<(neiList[dir][j+1])) {
                npos[dir]=j;
            }
          }
        }
      }

//      std::cerr << "loc B size: " << localBoxSize << ", loc node pos: " << localNodePos << ", max domain size in cells: " 
//          << maxDomainSizeInCells[0] << ", " <<maxDomainSizeInCells[1] << ", " <<  maxDomainSizeInCells[2] << ", kCellGridSize: " 
//          << kCellGridSize[0] << ", " <<kCellGridSize[1] << ", " <<  kCellGridSize[2] << ", pos: " << pos << ", cpos: " << cpos << ", npos: " << npos << std::endl;
      return mapPositionToIndex(npos);
    }

    void NodeGrid::calcNodeNeighbors(longint node)
    {
      Int3D nPos;
  
      mapIndexToPosition(nodePos, node);

      LOG4ESPP_DEBUG(logger, "my position: "
		     << node << " -> "
		     << nodePos[0] << " "
		     << nodePos[1] << " "
		     << nodePos[2]);

      for(int dir = 0; dir < 3; ++dir) {
	for(int j = 0; j < 3; ++j) {
	  nPos[j] = nodePos[j];
	}

	// left neighbor in direction dir
	nPos[dir] = nodePos[dir] - 1;
	if(nPos[dir] < 0) {
	  nPos[dir] += getGridSize(dir);
	}
	nodeNeighbors[2*dir] = mapPositionToIndex(nPos);
	LOG4ESPP_DEBUG(logger, "left neighbor in dir " << dir << ": "
		       << getNodeNeighborIndex(2*dir) << " <-> "
		       << nPos[0] << " "
		       << nPos[1] << " "
		       << nPos[2]);

	// right neighbor in direction dir
	nPos[dir] = nodePos[dir] + 1;
	if(nPos[dir] >= getGridSize(dir)) {
	  nPos[dir] -= getGridSize(dir);
	}
	nodeNeighbors[2*dir + 1] = mapPositionToIndex(nPos);

	LOG4ESPP_DEBUG(logger, "right neighbor in dir " << dir << ": "
		       << getNodeNeighborIndex(2*dir+1) << " <-> "
		       << nPos[0] << " "
		       << nPos[1] << " "
		       << nPos[2]);


    // todo: probably has to be altered... if(nodePos[dir] < frame....
	// left or right boundary ?
	boundaries[2*dir] = (nodePos[dir] == 0) ? ToRight : 0;
	boundaries[2*dir + 1] = (nodePos[dir] == getGridSize(dir) - 1) ? ToLeft : 0;

	LOG4ESPP_DEBUG(logger, "boundaries in dir " << dir << ": "
		       << getBoundary(2*dir) << " "
		       << getBoundary(2*dir + 1));
      }
    }

  }
}
