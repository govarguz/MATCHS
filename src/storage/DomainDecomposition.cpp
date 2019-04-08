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


#include "python.hpp"

#include <memory>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include "log4espp.hpp"
#include "System.hpp"

#include "Real3D.hpp"
#include "DomainDecomposition.hpp"
#include "bc/BC.hpp"
#include "Int3D.hpp"
#include "Buffer.hpp"

#include "iterator/CellListIterator.hpp"
#include "esutil/Error.hpp"
#include "esutil/Assert.hpp"

#include "boost/serialization/vector.hpp"
#include <vector>
#include <boost/python.hpp>


using namespace boost;
using namespace std;

namespace espressopp { 
  namespace storage {


  const int DD_COMM_TAG = 0xab;

  LOG4ESPP_LOGGER(DomainDecomposition::logger, "DomainDecomposition");

  std::string formatMismatchMessage(const Int3D& gridRequested,
                    int nodesAvailable) {
    std::ostringstream out;
    out << "requested node grid (" << gridRequested
            << ") does not match number of nodes in the communicator (" << nodesAvailable << ")";
    return out.str();
  }

  NodeGridMismatch::
  NodeGridMismatch(const Int3D& gridRequested, int nodesAvailable)
    : std::invalid_argument
  (formatMismatchMessage(gridRequested, nodesAvailable))
  {}

  DomainDecomposition::
  DomainDecomposition(shared_ptr< System > _system,
          int _halfCellInt,
          const Int3D& _nodeGrid,
          const Int3D& _cellGrid,
          const boost::python::list& neiListx,
          const boost::python::list& neiListy,
          const boost::python::list& neiListz)
    : Storage(_system), exchangeBufferSize(0), halfCellInt(_halfCellInt) {
    LOG4ESPP_INFO(logger, "node grid = "
          << _nodeGrid[0] << "x" << _nodeGrid[1] << "x" << _nodeGrid[2]
          << " cell grid = "
          << _cellGrid[0] << "x" << _cellGrid[1] << "x" << _cellGrid[2]
          << " Neighbor List size = "
          << boost::python::len(neiListx) << "x" << boost::python::len(neiListy) << "x" << boost::python::len(neiListz));

    ///TODO:
    ///TODO: unify whether number of cells in cell grid is multiplied by halfCellInt!!!
    ///TODO:
    nodeGrid = NodeGrid(_nodeGrid, getSystem()->comm->rank(), getSystem()->bc->getBoxL(), neiListx, neiListy, neiListz);
    if (nodeGrid.getNumberOfCells() != getSystem()->comm->size()) {
      throw NodeGridMismatch(_nodeGrid, getSystem()->comm->size());
    }
    construct(_cellGrid); // Hcomment: this inlcludes CreateCG, initCellInts and prepGhostComm
    LOG4ESPP_DEBUG(logger, "done");
  }

  DomainDecomposition::
  DomainDecomposition(shared_ptr< System > _system,
          int _halfCellInt,
          const Int3D& _nodeGrid,
          const Int3D& _cellGrid,
          const std::vector<int>& neiListx,
          const std::vector<int>& neiListy,
          const std::vector<int>& neiListz)
    : Storage(_system), exchangeBufferSize(0), halfCellInt(_halfCellInt) {
    LOG4ESPP_INFO(logger, "node grid = "
          << _nodeGrid[0] << "x" << _nodeGrid[1] << "x" << _nodeGrid[2]
          << " cell grid = "
          << _cellGrid[0] << "x" << _cellGrid[1] << "x" << _cellGrid[2]
          << " Neighbor List size = "
          << neiListx.size() << "x" << neiListy.size() << "x" << neiListz.size());

    ///TODO:
    ///TODO: unify whether number of cells in cell grid is multiplied by halfCellInt!!!
    ///TODO:
    nodeGrid = NodeGrid(_nodeGrid, getSystem()->comm->rank(), getSystem()->bc->getBoxL(), neiListx, neiListy, neiListz);
    if (nodeGrid.getNumberOfCells() != getSystem()->comm->size()) {
      throw NodeGridMismatch(_nodeGrid, getSystem()->comm->size());
    }
    construct(_cellGrid);
    LOG4ESPP_DEBUG(logger, "done");
  }

  DomainDecomposition::
  DomainDecomposition(shared_ptr< System > _system,
          int _halfCellInt,
          const Int3D& _nodeGrid,
          const std::vector<int>& neiListx,
          const std::vector<int>& neiListy,
          const std::vector<int>& neiListz)
    : Storage(_system), exchangeBufferSize(0), halfCellInt(_halfCellInt){
    nodeGrid = NodeGrid(_nodeGrid, getSystem()->comm->rank(), getSystem()->bc->getBoxL(), neiListx, neiListy, neiListz);
    if (nodeGrid.getNumberOfCells() != getSystem()->comm->size()) {
      throw NodeGridMismatch(_nodeGrid, getSystem()->comm->size());
    }

    ///TODO:
    ///TODO: unify whether number of cells in cell grid is multiplied by halfCellInt!!!
    ///TODO:
    Int3D cellGrid(neiListx[nodeGrid.getNodePosition(0) + 1] - neiListx[nodeGrid.getNodePosition(0)],
                   neiListy[nodeGrid.getNodePosition(1) + 1] - neiListy[nodeGrid.getNodePosition(1)],
                   neiListz[nodeGrid.getNodePosition(2) + 1] - neiListz[nodeGrid.getNodePosition(2)]);

    LOG4ESPP_INFO(logger, "node grid = "
          << _nodeGrid[0] << "x" << _nodeGrid[1] << "x" << _nodeGrid[2]
          << " cell grid = "
          << cellGrid[0] << "x" << cellGrid[1] << "x" << cellGrid[2]
          << " Neighbor List size = "
          << neiListx.size() << "x" << neiListy.size() << "x" << neiListz.size());

    construct(cellGrid);
    LOG4ESPP_DEBUG(logger, "done");
  }


  DomainDecomposition::
  DomainDecomposition(shared_ptr< System > _system,
          int _halfCellInt,
          const Int3D& _nodeGrid,
          const boost::python::list& neiListx,
          const boost::python::list& neiListy,
          const boost::python::list& neiListz)
    : DomainDecomposition(_system, _halfCellInt, _nodeGrid,
                          std::vector<int>(boost::python::stl_input_iterator<int>(neiListx),boost::python::stl_input_iterator<int>()),
                          std::vector<int>(boost::python::stl_input_iterator<int>(neiListy),boost::python::stl_input_iterator<int>()),
                          std::vector<int>(boost::python::stl_input_iterator<int>(neiListz),boost::python::stl_input_iterator<int>()))
  {
  }

  DomainDecomposition::
  DomainDecomposition(shared_ptr< System > _system,
          const Int3D& _nodeGrid,
          const std::vector<int>& neiListx,
          const std::vector<int>& neiListy,
          const std::vector<int>& neiListz,
          const std::vector<int>& halfCellMaskx,
          const std::vector<int>& halfCellMasky,
          const std::vector<int>& halfCellMaskz)
    : Storage(_system), exchangeBufferSize(0) { // todo: move halfCellInt from Storage to DD?
    nodeGrid = NodeGrid(_nodeGrid, getSystem()->comm->rank(), getSystem()->bc->getBoxL(), neiListx, neiListy, neiListz);
    if (nodeGrid.getNumberOfCells() != getSystem()->comm->size()) {
      throw NodeGridMismatch(_nodeGrid, getSystem()->comm->size());
    }

    ESPR_ASSERT(neiListx.size() == halfCellMaskx.size() + 1);
    ESPR_ASSERT(neiListy.size() == halfCellMasky.size() + 1);
    ESPR_ASSERT(neiListz.size() == halfCellMaskz.size() + 1);

    ///TODO:
    ///TODO: unify whether number of cells in cell grid is multiplied by halfCellInt!!!
    ///TODO:
    /// In this particular constructor, we do multiply it
    /// It means, that neiLists are refering to number of cells in the coarse sense
    /// so if neiList[i+1] - neiList[i] = a, processor i has a*halfCellInt cells,
    /// where halfCellInt is determined from the halfCellLists.
    ///
    Int3D nodePos = nodeGrid.getNodePosition();

    halfCellInt = min(halfCellMaskx[nodePos[0]], min(halfCellMasky[nodePos[1]], halfCellMaskz[nodePos[2]]));
    Int3D cellGrid((neiListx[nodePos[0]+ 1] - neiListx[nodePos[0]]) * halfCellInt,
                   (neiListy[nodePos[1] + 1] - neiListy[nodePos[1]]) * halfCellInt,
                   (neiListz[nodePos[2] + 1] - neiListz[nodePos[2]]) * halfCellInt);

    LOG4ESPP_INFO(logger, "node grid = "
          << _nodeGrid[0] << "x" << _nodeGrid[1] << "x" << _nodeGrid[2]
          << " cell grid = "
          << cellGrid[0] << "x" << cellGrid[1] << "x" << cellGrid[2]
          << " Neighbor List size = "
          << neiListx.size() << "x" << neiListy.size() << "x" << neiListz.size());

    construct(cellGrid);
    LOG4ESPP_DEBUG(logger, "done");
  }

  DomainDecomposition::
  DomainDecomposition(shared_ptr< System > _system,
          const Int3D& _nodeGrid,
          const boost::python::list& neiListx,
          const boost::python::list& neiListy,
          const boost::python::list& neiListz,
          const boost::python::list& halfCellMaskx,
          const boost::python::list& halfCellMasky,
          const boost::python::list& halfCellMaskz)
    : DomainDecomposition(_system, _nodeGrid,
                          std::vector<int>(boost::python::stl_input_iterator<int>(neiListx),boost::python::stl_input_iterator<int>()),
                          std::vector<int>(boost::python::stl_input_iterator<int>(neiListy),boost::python::stl_input_iterator<int>()),
                          std::vector<int>(boost::python::stl_input_iterator<int>(neiListz),boost::python::stl_input_iterator<int>()),
                          std::vector<int>(boost::python::stl_input_iterator<int>(halfCellMaskx),boost::python::stl_input_iterator<int>()),
                          std::vector<int>(boost::python::stl_input_iterator<int>(halfCellMasky),boost::python::stl_input_iterator<int>()),
                          std::vector<int>(boost::python::stl_input_iterator<int>(halfCellMaskz),boost::python::stl_input_iterator<int>()))
  {
  }


  void DomainDecomposition::construct(const Int3D &_cellGrid)
  {
    createCellGrid(_cellGrid);
    initCellInteractions();
    prepareGhostCommunication();
  }

  void DomainDecomposition::createCellGrid(const Int3D& _cellGrid)
  {
    real myLeft[3];
    real myRight[3];

    LOG4ESPP_INFO(logger, "my node grid position: "
          << nodeGrid.getNodePosition(0) << " "
          << nodeGrid.getNodePosition(1) << " "
          << nodeGrid.getNodePosition(2) << " -> "
          << getSystem()->comm->rank());

    LOG4ESPP_DEBUG(logger, "my neighbors: "
           << nodeGrid.getNodeNeighborIndex(0) << "<->"
           << nodeGrid.getNodeNeighborIndex(1) << ", "
           << nodeGrid.getNodeNeighborIndex(2) << "<->"
           << nodeGrid.getNodeNeighborIndex(3) << ", "
           << nodeGrid.getNodeNeighborIndex(4) << "<->"
           << nodeGrid.getNodeNeighborIndex(5));

    for (int i = 0; i < 3; ++i) {
      myLeft[i] = nodeGrid.getMyLeft(i);
      myRight[i] = nodeGrid.getMyRight(i);
    }

    cellGrid = CellGrid(_cellGrid, myLeft, myRight, halfCellInt);

    LOG4ESPP_INFO(logger, "local box "
          << myLeft[0] << "-" << myRight[0] << ", "
          << myLeft[1] << "-" << myRight[1] << ", "
          << myLeft[2] << "-" << myRight[2]);

    longint nLocalCells = 1;
    longint nRealCells = 1;
    for (int i = 0; i < 3; ++i) {
      nRealCells *= cellGrid.getGridSize(i);
      nLocalCells *= cellGrid.getFrameGridSize(i);
    }

    resizeCells(nLocalCells);

    realCells.reserve(nRealCells);
    ghostCells.reserve(nLocalCells - nRealCells);

    markCells();

    LOG4ESPP_DEBUG(logger, "total # cells=" << nLocalCells
           << ", # real cells=" << nRealCells
           << ", frame cell grid = (" << cellGrid.getFrameGridSize(0)
           << ", " << cellGrid.getFrameGridSize(1)
           << ", " << cellGrid.getFrameGridSize(2)
           << ")");
  }

  void DomainDecomposition::markCells() {
    realCells.resize(0);
    ghostCells.resize(0);
    realBoundaryCells.resize(0);

    for (int o = 0; o < cellGrid.getFrameGridSize(2); ++o) {
      for (int n = 0; n < cellGrid.getFrameGridSize(1); ++n) {
        for (int m = 0; m < cellGrid.getFrameGridSize(0); ++m) {
          Cell *cur = &cells[cellGrid.mapPositionToIndex(m, n, o)];
#ifdef CELL_EXTRA_DATA
          cur->id = cellGrid.mapPositionToIndex(m, n, o);
          Int3D intIdx(m, n, o);
          cur->grid_pos = intIdx;
          for(int i = 0; i < 3; i++) {
            cur->myLeft[i] = nodeGrid.getMyLeft(i) + (intIdx[i] - halfCellInt) * cellGrid.getCellSize(i);
            cur->myRight[i] = nodeGrid.getMyLeft(i) + (intIdx[i] - halfCellInt + 1) * cellGrid.getCellSize(i);
          }
#endif
          if (cellGrid.isInnerCell(m, n, o)) {
            LOG4ESPP_TRACE(logger, "cell " << (cur - &cells[0]) << " is inner cell (" << m << ", " << n << ", " << o << ")");
            realCells.push_back(cur);
            if (!cellGrid.isInnerNoCommCell(m, n, o)){
                realBoundaryCells.push_back(cur);
            }
          } else {
            LOG4ESPP_TRACE(logger, "cell " << (cur - &cells[0]) << " is ghost cell (" << m << ", " << n << ", " << o << ")");
            ghostCells.push_back(cur);
          }
        }
      }
    }
  }

  // TODO one should take care of rc and system size
  /* NOTE to E++ Developers: There are two possible ways of 'taking care' of rc
 the first one is by redistributing the cell grid. Which can be done within the HeSpaDDA
framework by restarting the DomainDecomposition. The second would be by an upcoming development
where different cell sizes can be dynamically arranged. For details (horacio.v.g@gmail.com)*/
  /** scale position coordinates of all real particles by factor s */
  void DomainDecomposition::scaleVolume(real s, bool particleCoordinates){
	if(particleCoordinates) Storage::scaleVolume( s );
    
    real maxCut = getSystem() -> maxCutoff;
    real skinL = getSystem() -> getSkin();
    real cs = maxCut + skinL;
    if( cs > s*cellGrid.getSmallestCellDiameter() ){
      Real3D Li = getSystem() -> bc -> getBoxL(); // getting the system size
      real minL = min(Li[0], min(Li[1],Li[2]));
      if(cs > minL){
        esutil::Error err(getSystemRef().comm);
        stringstream msg;
        msg<<"Error. The current system size "<< minL <<" smaller then cutoff+skin "<< cs;
        err.setException( msg.str() );
      }
      else{
        cellAdjust();
      }
    }
    else{
      cellGrid.scaleVolume( s );
      nodeGrid.scaleVolume( s );
    }
  }
  // anisotropic version
  void DomainDecomposition::scaleVolume(Real3D s, bool particleCoordinates){
	if(particleCoordinates) Storage::scaleVolume( s );
    
    real maxCut = getSystem() -> maxCutoff;
    real skinL = getSystem() -> getSkin();
    real cs = maxCut + skinL;
    real cellD = cellGrid.getSmallestCellDiameter();
    
    real r0 = s[0]*cellD;
    real r1 = s[1]*cellD;
    real r2 = s[2]*cellD;
    
    if( cs > min( min( r0, r1), r2 ) ){
      Real3D Li = getSystem() -> bc -> getBoxL(); // getting the system size
      real minL = min(Li[0], min(Li[1],Li[2]));
      if(cs > minL){
        esutil::Error err(getSystemRef().comm);
        stringstream msg;
        msg<<"Error. The current system size "<< minL <<" smaller then cutoff+skin "<< cs;
        err.setException( msg.str() );
      }
      else
        cellAdjust();
    }
    else{
      cellGrid.scaleVolume(s);
      nodeGrid.scaleVolume(s);
    }
  }
  
  Int3D DomainDecomposition::getInt3DCellGrid(){
    return Int3D( cellGrid.getGridSize(0),
                  cellGrid.getGridSize(1),
                  cellGrid.getGridSize(2)
                );
  }
  Int3D DomainDecomposition::getInt3DNodeGrid(){
    return Int3D( nodeGrid.getGridSize(0),
                  nodeGrid.getGridSize(1),
                  nodeGrid.getGridSize(2)
                );
  }

  void DomainDecomposition::cellAdjust(){   
    /* // create an appropriate cell grid
    Real3D box_sizeL = getSystem() -> bc -> getBoxL();
    real skinL = getSystem() -> getSkin();
    real maxCutoffL = getSystem() -> maxCutoff;
            
    // nodeGrid is already defined
    Int3D _nodeGrid(nodeGrid.getGridSize());
    // new cellGrid
    real rc_skin = maxCutoffL + skinL;
    int ix = (int)(box_sizeL[0] / (rc_skin * _nodeGrid[0]));
    int iy = (int)(box_sizeL[1] / (rc_skin * _nodeGrid[1]));
    int iz = (int)(box_sizeL[2] / (rc_skin * _nodeGrid[2]));
    Int3D _newCellGrid(ix, iy, iz);

    // save all particles to temporary vector
    std::vector<ParticleList> tmp_pl;
    size_t _N = realCells.size();
    tmp_pl.reserve( _N );
    for(CellList::Iterator it(realCells); it.isValid(); ++it) {
      tmp_pl.push_back((*it)->particles);
    }
    */
    // reset all cells info
    invalidateGhosts();
    cells.clear();
    localCells.clear();
    realCells.clear();
    ghostCells.clear();
    realBoundaryCells.clear();
    for(int i=0; i<6; i++){
      commCells[i].reals.clear();
      commCells[i].ghosts.clear();
    }
    /*
    // ################################   H   ############################################
    // Check for variable pressure simulations!!!
    // createCellGrid(_nodeGrid, _newCellGrid);
    // Has been deactivated
    // ###################################################################################

    // creating new grids
    createCellGrid(_nodeGrid, _newCellGrid);
    initCellInteractions();
    prepareGhostCommunication();
    
    // pushing the particles back to the empty cell ("do we have to check particles?")
    for(int i=0; i<tmp_pl.size(); i++){
      for (size_t p = 0; p < tmp_pl[i].size(); ++p) {
        Particle& part = tmp_pl[i][p];
        const Real3D& pos = part.position();
        Cell *sortCell = mapPositionToCellClipped(pos);
        appendUnindexedParticle(sortCell->particles, part);
      }
    }

    for(CellList::Iterator it(realCells); it.isValid(); ++it) {
      updateLocalParticles((*it)->particles);
    }
    
    exchangeGhosts();
    onParticlesChanged();
    */
  }

  void DomainDecomposition::initCellInteractions() {
    LOG4ESPP_DEBUG(logger, "setting up neighbors for " << cells.size() << " cells");
    //TODO, sadly, I cannot use getSystem() -> maxCutoff here, since DomainDecomposition is constructed BEFORE the interactions are registered
    //TODO thus I introduced this new paraemter for the constructor. It is not nice, but needed here
    //real maxCut = getSystem() -> maxCutoff;
    real maxCut = getSystem() -> expectedMaxCutoff;
    std::cout << "expectedMaxCutoff: " << getSystem()->expectedMaxCutoff << std::endl;
    if(maxCut < 0.0) {
        throw std::runtime_error("system.expectedMaxCutoff has not been set. It needs to be set in order to allow skiping some cells from neighbors");
    }
    real skinL = getSystem() -> getSkin();
    real cs = (maxCut + skinL);
    for(int i = 0; i < 3; i++)
        if(cellGrid.getInnerCellsBegin(i) != halfCellInt)
            throw std::runtime_error("DomainDecomposition:InitCellInteractions: cell grid not complying with halfCell approach");
    for (int o = cellGrid.getInnerCellsBegin(2); o < cellGrid.getInnerCellsEnd(2); ++o) {
      // if it was a full cell grid, this cell woudl be contained in cell with o : 
      int o_coarse = o / halfCellInt;
      for (int n = cellGrid.getInnerCellsBegin(1); n < cellGrid.getInnerCellsEnd(1); ++n) {
        int n_coarse = n / halfCellInt;
        for (int m = cellGrid.getInnerCellsBegin(0); m < cellGrid.getInnerCellsEnd(0); ++m) {
          int m_coarse = m / halfCellInt;
          longint cellIdx = cellGrid.mapPositionToIndex(m, n, o);
          Cell *cell = &cells[cellIdx];

          // index of a "suppercell" in fictional grid wit full cells
          // it is needed to create ordering of cells consistent over nodes with different halfCellInt
          longint cellIdx_coarse = cellGrid.mapPositionToIndex(m_coarse, n_coarse, o_coarse);

          LOG4ESPP_TRACE(logger, "setting up neighbors for cell " << cell - getFirstCell()
                << " @ " << m << " " << n << " " << o);

          cell->neighborCells.reserve((2*halfCellInt + 1) * (2*halfCellInt + 1) * (2*halfCellInt + 1) - 1);

          // loop all neighbor cells
          for (int p = o - halfCellInt; p <= o + halfCellInt; ++p) {
            int p_coarse = p / halfCellInt;
            for (int q = n - halfCellInt; q <= n + halfCellInt; ++q) {
              int q_coarse = q / halfCellInt;
              for (int r = m - halfCellInt; r <= m + halfCellInt; ++r) {
                int r_coarse = r / halfCellInt;

                // do not consider cell itself as neigbor
                if (p == o && q == n && r == m)
                  continue;

                // in the halfcell (or smaller) case, it is possible, that some of the possibly neigboring
                // cells are so far away, that they do not need to be considered
                Int3D cell_dist_int = Int3D(m, n, o) - Int3D(r, q, p);
                Real3D cell_dist;
                for (int i = 0; i < 3; i++){
                  if(cell_dist_int[i] != 0){
                    //this will give us the minimal distance between the two cells in i direction 
                    //(measured in multiples of cell size)
                    cell_dist_int[i] = abs(cell_dist_int[i]) - 1;
                  }
                  cell_dist[i] = cell_dist_int[i] * cellGrid.getCellSize(i);
                }
                if(cell_dist.abs() > cs){
                  std::cout << "dist " << cell_dist.abs() << ", cutoff+skin " << cs << ", cdi " <<  cell_dist_int << "  , cellsize " <<cellGrid.getCellSize(0)
                                      << ", " << cellGrid.getCellSize(1)<< ", " << cellGrid.getCellSize(2) << std::endl;
                  continue;
                }

                longint cell2Idx = cellGrid.mapPositionToIndex(r, q, p);
                Cell *cell2 = &cells[cell2Idx];

                // this is a sort of global ordering of cells. It has to be consistent over nodes,
                // in order that iterator over particle pairs in creation of verlet list (and others)
                // consideres each GLOBAL pair of particles exactly once (it is tricky especially by 
                // the border, where 1 particle is real on 1 processor and 2. particle on another. Than
                // it has to be clear which processor takes care of this pair.
                // It is even more tricky, when neigboring processors have different halfCellInt. Thus
                // this implementation...
                longint cell2Idx_coarse = cellGrid.mapPositionToIndex(r_coarse, q_coarse, p_coarse);
                bool neighbor_prefered;
                if(cell2Idx_coarse != cellIdx_coarse){
                    // if cells belong to different suppercells, do as it would be done for supercells
                    neighbor_prefered = (cell2Idx_coarse < cellIdx_coarse);
                }
                else{
                    // if they belong to the same supercell, it is not possible, that cell2 is ghost
                    // (and cell is never ghost). Both cells are thus real at one processor and thus 
                    // order is not really important (does not need to be consistent with anything other
                    neighbor_prefered = (cell2Idx < cellIdx);
                }

                cell->neighborCells.push_back(NeighborCellInfo(cell2, (neighbor_prefered)));
                //cell->neighborCells.push_back(NeighborCellInfo(cell2, (cell2Idx<cellIdx)));

                LOG4ESPP_TRACE(logger, "neighbor cell " << cell2 - getFirstCell()
                      << " @ " << r << " " << q << " " << p << ((cell2Idx<cellIdx) ? " is" : " is not") << " taken" );
              }
            }
          }

           int expected = (2*halfCellInt + 1) * (2*halfCellInt + 1) * (2*halfCellInt + 1) - 1;
           if(cell->neighborCells.size() < expected)
               std::cout << cell->neighborCells.size() << "neighbor cells instead of " << expected << std::endl;
        }
      }
    }

    LOG4ESPP_DEBUG(logger, "done");
  }

  Cell *DomainDecomposition::mapPositionToCell(const Real3D& pos) {
    return &cells[cellGrid.mapPositionToCell(pos)];
  }

  Cell *DomainDecomposition::mapPositionToCellWithGhosts(const Real3D& pos) {
    return &cells[cellGrid.mapPositionToCellWithGhosts(pos)];
  }

  Cell *DomainDecomposition::mapPositionToCellClipped(const Real3D& pos) {
    return &cells[cellGrid.mapPositionToCellClipped(pos)];
  }

  Cell *DomainDecomposition::mapPositionToCellChecked(const Real3D& pos) {
    longint c = cellGrid.mapPositionToCellChecked(pos);
    if (c == CellGrid::noCell) {
      return 0;
    }
    else{
      return &cells[c];
    }
  }

  longint DomainDecomposition::mapPositionToNodeClipped(const Real3D& pos) {
    return nodeGrid.mapPositionToNodeClipped(pos);
  }

  bool DomainDecomposition::checkIsRealParticle(longint id, const Real3D& pos) {
    bool isReal = getSystem()->comm->rank() == mapPositionToNodeClipped(pos);
    bool isInDomain = true;
    for (int i = 0; i < 3; i++) {
      if(pos[i] < nodeGrid.getMyLeft()[i]) isInDomain = false;
      if(pos[i] > nodeGrid.getMyRight()[i]) isInDomain = false;
    }
    if(isReal != isInDomain) {
        std::cout << "WARNING: is real: " << pos << ", node: " + std::string(nodeGrid.getMyLeft()) + "<->" + std::string(nodeGrid.getMyRight())
          << ", answer: " << getSystem()->comm->rank() << "==" << mapPositionToNodeClipped(pos) << std::endl;
    }

    return isReal;
  }

  bool DomainDecomposition::appendParticles(ParticleList &l, int dir) {
    bool outlier = false;

    LOG4ESPP_DEBUG(logger, "got " << l.size() << " particles");

    for (ParticleList::iterator it = l.begin(), end = l.end(); it != end; ++it) {
      Real3D& pos = it->position();

      if (nodeGrid.getBoundary(dir) != 0) {
          getSystem()->bc->foldCoordinate(pos, it->image(), nodeGrid.convertDirToCoord(dir));
          LOG4ESPP_TRACE(logger, "folded coordinate " << nodeGrid.convertDirToCoord(dir)
                  << " of particle " << it->id());
      }

      longint cell;
      if (cellGrid.mapPositionToCellCheckedAndClipped(cell, pos)) {
          LOG4ESPP_TRACE(logger, "particle " << it->id()
                  << " @ " << pos << " is not inside node domain");
          outlier = true;
      }

      LOG4ESPP_TRACE(logger, "append part " << it->id() << " to cell "
              << cell);

      appendIndexedParticle(cells[cell].particles, *it);
    }
    return outlier;
  }

  void DomainDecomposition::decomposeRealParticles() {

     //std::cout << getSystem()->comm->rank() << ": " << " decomposeRealParticles\n";
     //std::cout << getSystem()->comm->rank() << "localCells size: " <<  localCells.size() << std::endl;

    LOG4ESPP_DEBUG(logger, "starting, expected comm buffer size " << exchangeBufferSize);

    // allocate send/recv buffers. We use the size as we need maximally so far, to avoid reallocation
    // TODO: This might be a problem when all particles are created on a single node initially!
    ParticleList sendBufL;
    sendBufL.reserve(exchangeBufferSize);
    ParticleList sendBufR;
    sendBufR.reserve(exchangeBufferSize);
    ParticleList recvBufL;
    recvBufL.reserve(exchangeBufferSize);
    ParticleList recvBufR;
    recvBufR.reserve(exchangeBufferSize);

    bool allFinished;
    do {
      bool finished = true;

      for (int coord = 0; coord < 3; ++coord) {
        LOG4ESPP_DEBUG(logger, "starting with direction " << coord);

        if (nodeGrid.getGridSize(coord) > 1) {
          for (std::vector<Cell*>::iterator it = realCells.begin(),
            end = realCells.end(); it != end; ++it) {

            Cell &cell = **it;

            // do not use an iterator here, since we need to take out particles during the loop
            for (size_t p = 0; p < cell.particles.size(); ++p) {
              Particle &part = cell.particles[p];
              const Real3D& pos = part.position();

              // check whether the particle is now "left" of the local domain
              if (pos[coord] - cellGrid.getMyLeft(coord) < -ROUND_ERROR_PREC) {
                LOG4ESPP_TRACE(logger, "send particle left " << part.id());
                moveIndexedParticle(sendBufL, cell.particles, p);
                // redo same particle since we took one out here, so it's a new one
                --p;
              }
              // check whether the particle is now "right" of the local domain
              else if (pos[coord] - cellGrid.getMyRight(coord) >= ROUND_ERROR_PREC) {
                LOG4ESPP_TRACE(logger, "send particle right " << part.id());
                moveIndexedParticle(sendBufR, cell.particles, p);
                --p;
              }
              // Sort particles in cells of this node during last direction
              else if (coord == 2) {
                const Real3D& pos = part.position();
                Cell *sortCell = mapPositionToCellChecked(pos);
                if (sortCell != &cell) {
                  if (sortCell == 0) {
                    // particle is not in the local domain
                    LOG4ESPP_DEBUG(logger, "take another loop: particle " << part.id()
                            << " @ " << pos <<
                            " is not inside node domain after neighbor exchange");
                    // isnan function is C99 only, x != x is only true if x == nan
                    if (pos[0] != pos[0] || pos[1] != pos[1] || pos[2] != pos[2]) {
                      // TODO: error handling
                      LOG4ESPP_ERROR(logger, "particle " << part.id() <<
                              " has moved to outer space (one or more coordinates are nan)");
                    } else {
                      // particle stays where it is, and will be sorted in the next round
                      finished = false;
                    }
                  } else {
                    // particle is in the local domain
                    moveIndexedParticle(sortCell->particles, cell.particles, p);
                    --p;
                  }
                }
              }

            }
          }

          // Exchange particles, odd-even rule
          if (nodeGrid.getNodePosition(coord) % 2 == 0) {
            sendParticles(sendBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
            recvParticles(recvBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
            sendParticles(sendBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
            recvParticles(recvBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
          } else {
            recvParticles(recvBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
            sendParticles(sendBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
            recvParticles(recvBufL, nodeGrid.getNodeNeighborIndex(2 * coord));
            sendParticles(sendBufR, nodeGrid.getNodeNeighborIndex(2 * coord + 1));
          }

          // sort received particles to cells
          if (appendParticles(recvBufL, 2 * coord) && coord == 2)
              finished = false;
          if (appendParticles(recvBufR, 2 * coord + 1) && coord == 2)
              finished = false;

          // reset send/recv buffers
          sendBufL.resize(0);
          sendBufR.resize(0);
          recvBufL.resize(0);
          recvBufR.resize(0);


        } else {
          /* Single node direction case (no communication)
            Fold particles that have left the box */
          for (std::vector< Cell* >::iterator it = realCells.begin(),
            end = realCells.end(); it != end; ++it) {
            Cell &cell = **it;
            // do not use an iterator here, since we have need to take out particles during the loop
            for (size_t p = 0; p < cell.particles.size(); ++p) {
                Particle &part = cell.particles[p];
                getSystem()->bc->foldCoordinate(part.position(), part.image(), coord);
                LOG4ESPP_TRACE(logger, "folded coordinate " << coord << " of particle " << part.id());

                if (coord == 2) {
                    Cell *sortCell = mapPositionToCellChecked(part.position());

                    if (sortCell != &cell) {
                        if (sortCell == 0) {
                            LOG4ESPP_DEBUG(logger, "take another loop: particle " << part.id()
                                    << " @ " << part.position()
                                    << " is not inside node domain after neighbor exchange");
                            const Real3D& pos = part.position();
                            // isnan function is C99 only, x != x is only true if x == nan
                            if (pos[0] != pos[0] || pos[1] != pos[1] || pos[2] != pos[2]) {
                                LOG4ESPP_ERROR(logger, "particle " << part.id() <<
                                        " has moved to outer space (one or more coordinates are nan)");
                            } else {
                                // particle stays where it is, and will be sorted in the next round
                                finished = false;
                            }
                        } else {
                            moveIndexedParticle(sortCell->particles, cell.particles, p);
                            --p;
                        }
                    }
                }
            }
          }
        }

        LOG4ESPP_DEBUG(logger, "done with direction " << coord);
      }

      // Communicate wether particle exchange is finished
      mpi::all_reduce(*getSystem()->comm, finished, allFinished, std::logical_and<bool>());
    } while (!allFinished);

    exchangeBufferSize = std::max(exchangeBufferSize,
                  std::max(sendBufL.capacity(),
                       std::max(sendBufR.capacity(),
                            std::max(recvBufL.capacity(),
                                 recvBufR.capacity()))));

    LOG4ESPP_DEBUG(logger, "finished exchanging particles, new send/recv buffer size " << exchangeBufferSize);

    LOG4ESPP_DEBUG(logger, "done");
  }

  void DomainDecomposition::exchangeGhosts() {
    LOG4ESPP_DEBUG(logger, "exchangeGhosts -> ghost communication sizes first, real->ghost");
    doGhostCommunication(true, true, dataOfExchangeGhosts);
  }

  void DomainDecomposition::updateGhosts() {
    LOG4ESPP_DEBUG(logger, "updateGhosts -> ghost communication no sizes, real->ghost");
    doGhostCommunication(false, true, dataOfUpdateGhosts);
  }

  void DomainDecomposition::updateGhostsV() {
    LOG4ESPP_DEBUG(logger, "updateGhostsV -> ghost communication no sizes, real->ghost velocities");
    doGhostCommunication(false, true, 2); // 2 is the bitflag for particle momentum
  }

  void DomainDecomposition::collectGhostForces() {
    LOG4ESPP_DEBUG(logger, "collectGhosts -> ghost communication no sizes, ghost->real");
    doGhostCommunication(false, false);
  }

  void DomainDecomposition::fillCells(std::vector<Cell *> &cv,
                  const int leftBoundary[3],
                  const int rightBoundary[3]) {
    LOG4ESPP_DEBUG(logger, "filling: "
           << leftBoundary[0] << "-" << (rightBoundary[0] - 1) << " "
           << leftBoundary[1] << "-" << (rightBoundary[1] - 1) << " "
           << leftBoundary[2] << "-" << (rightBoundary[2] - 1));

    longint total = 1;
    for (int i = 0; i < 3; ++i) {
        if (leftBoundary[i] < 0 || leftBoundary[i] > cellGrid.getFrameGridSize(i) ||
                rightBoundary[i] < 0 || rightBoundary[i] > cellGrid.getFrameGridSize(i) ||
                leftBoundary[i] >= rightBoundary[i]) {
            throw std::runtime_error("DomainDecomposition::fillCells: wrong cell grid specified internally");
        }
        total *= (rightBoundary[i] - leftBoundary[i]);
    }
    cv.reserve(total);

    for (int o = leftBoundary[0]; o < rightBoundary[0]; ++o) {
        for (int n = leftBoundary[1]; n < rightBoundary[1]; ++n) {
            for (int m = leftBoundary[2]; m < rightBoundary[2]; ++m) {
                int i = cellGrid.mapPositionToIndex(o, n, m);
                LOG4ESPP_TRACE(logger, "add cell " << i);
                cv.push_back(&cells[i]);
            }
        }
    }

    LOG4ESPP_DEBUG(logger, "expected " << total << " cells, filled with " << cv.size());
  }

  void DomainDecomposition::prepareGhostCommunication() {
    // direction loop: x, y, z
    for (int coord = 0; coord < 3; ++coord) {
        // boundaries of area to send
        int leftBoundary[3], rightBoundary[3];
        /* boundaries perpendicular directions are the same for left/right send.
        We also send the ghost frame that we have already, so the data amount
        increase with each cycle.

        For a direction that was done already, i.e. is smaller than dir,
        we take the full ghost frame, otherwise only the inner frame.  */
        for (int offset = 1; offset <= 2; ++offset) {
            int otherCoord = (coord + offset) % 3;
            if (otherCoord < coord) {
                leftBoundary[otherCoord] = 0;
                rightBoundary[otherCoord] = cellGrid.getFrameGridSize(otherCoord);
            } else {
                leftBoundary[otherCoord] = cellGrid.getInnerCellsBegin(otherCoord);
                rightBoundary[otherCoord] = cellGrid.getInnerCellsEnd(otherCoord);
            }
        }

        //  lr loop: left right - loop
        for (int lr = 0; lr < 2; ++lr) {
            int dir = 2 * coord + lr;

            /* participating real particles from this node */
            LOG4ESPP_DEBUG(logger, "direction " << dir << " reals");

            if (lr == 0) {
                leftBoundary[coord] = cellGrid.getInnerCellsBegin(coord);
                rightBoundary[coord] = cellGrid.getInnerCellsBegin(coord) + cellGrid.getFrameWidth();
            } else {
                leftBoundary[coord] = cellGrid.getInnerCellsEnd(coord) - cellGrid.getFrameWidth();
                rightBoundary[coord] = cellGrid.getInnerCellsEnd(coord);
            }
            fillCells(commCells[dir].reals, leftBoundary, rightBoundary);

            /* participating ghosts from this node */
            LOG4ESPP_DEBUG(logger, "direction " << dir << " ghosts");

            if (lr == 0) {
                leftBoundary[coord] = cellGrid.getInnerCellsEnd(coord);
                rightBoundary[coord] = cellGrid.getInnerCellsEnd(coord) + cellGrid.getFrameWidth();
            } else {
                leftBoundary[coord] = cellGrid.getInnerCellsBegin(coord) - cellGrid.getFrameWidth();
                rightBoundary[coord] = cellGrid.getInnerCellsBegin(coord);
            }
            fillCells(commCells[dir].ghosts, leftBoundary, rightBoundary);
        }
    }
  }
  
  void DomainDecomposition::
  doGhostCommunication(bool sizesFirst, bool realToGhosts, int extradata) {
    LOG4ESPP_DEBUG(logger, "do ghost communication " << (sizesFirst ? "with sizes " : "")
           << (realToGhosts ? "reals to ghosts " : "ghosts to reals ") << extradata);

    /* direction loop: x, y, z.
   Here we could in principle build in a one sided ghost
   communication, simply by taking the lr loop only over one
   value. */
    for (int _coord = 0; _coord < 3; ++_coord) {
      /* inverted processing order for ghost force communication,
        since the corner ghosts have to be collected via several
        nodes. We now add back the corner ghost forces first again
        to ghost forces, which only eventually go back to the real
        particle.
      */
      int coord = realToGhosts ? _coord : (2 - _coord);
      real curCoordBoxL = getSystem()->bc->getBoxL()[coord];

      // lr loop: left right
      for (int lr = 0; lr < 2; ++lr) {
        int dir         = 2 * coord + lr;
        int oppositeDir = 2 * coord + (1 - lr);

        Real3D shift(0, 0, 0);

        shift[coord] = nodeGrid.getBoundary(dir) * curCoordBoxL;

        LOG4ESPP_DEBUG(logger, "direction " << dir);

        if (nodeGrid.getGridSize(coord) == 1) {
          LOG4ESPP_DEBUG(logger, "local communication");

          // copy operation, we have to receive as many cells as we send
          if (commCells[dir].ghosts.size() != commCells[dir].reals.size()) {
            throw std::runtime_error("DomainDecomposition::doGhostCommunication: send/recv cell structure mismatch during local copy");
          }

          for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i) {
            if (realToGhosts) {
              copyRealsToGhosts(*commCells[dir].reals[i], *commCells[dir].ghosts[i], extradata, shift);
            } else {
              addGhostForcesToReals(*commCells[dir].ghosts[i], *commCells[dir].reals[i]);
            }
          }
        }
        else {
          // exchange size information, if necessary
          if (sizesFirst) {
            LOG4ESPP_DEBUG(logger, "exchanging ghost cell sizes");

            // prepare buffers
            std::vector<longint> sendSizes, recvSizes;
            sendSizes.reserve(commCells[dir].reals.size());
            for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i) {
              sendSizes.push_back(commCells[dir].reals[i]->particles.size());
            }
            recvSizes.resize(commCells[dir].ghosts.size());

            // exchange sizes, odd-even rule
            if (nodeGrid.getNodePosition(coord) % 2 == 0) {
              LOG4ESPP_DEBUG(logger, "sending to node " << nodeGrid.getNodeNeighborIndex(dir)
                        << ", then receiving from node " << nodeGrid.getNodeNeighborIndex(oppositeDir));
              getSystem()->comm->send(nodeGrid.getNodeNeighborIndex(dir), DD_COMM_TAG, &(sendSizes[0]), sendSizes.size());
              getSystem()->comm->recv(nodeGrid.getNodeNeighborIndex(oppositeDir), DD_COMM_TAG, &(recvSizes[0]), recvSizes.size());
            }
            else {
              LOG4ESPP_DEBUG(logger, "receiving from node " << nodeGrid.getNodeNeighborIndex(oppositeDir)
                        << ", then sending to node " << nodeGrid.getNodeNeighborIndex(dir));
              getSystem()->comm->recv(nodeGrid.getNodeNeighborIndex(oppositeDir), DD_COMM_TAG, &(recvSizes[0]), recvSizes.size());
              getSystem()->comm->send(nodeGrid.getNodeNeighborIndex(dir), DD_COMM_TAG, &(sendSizes[0]), sendSizes.size());
            }

            // resize according to received information
            for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i) {
              commCells[dir].ghosts[i]->particles.resize(recvSizes[i]);
            }
            LOG4ESPP_DEBUG(logger, "exchanging ghost cell sizes done");
          }

          // prepare send and receive buffers
          longint receiver, sender;
          outBuffer.reset();
          if (realToGhosts) {
            receiver = nodeGrid.getNodeNeighborIndex(dir);
            sender = nodeGrid.getNodeNeighborIndex(oppositeDir);
            for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i) {
              packPositionsEtc(outBuffer, *commCells[dir].reals[i], extradata, shift);
            }
          }
          else {
            receiver = nodeGrid.getNodeNeighborIndex(oppositeDir);
            sender = nodeGrid.getNodeNeighborIndex(dir);
            for (int i = 0, end = commCells[dir].ghosts.size(); i < end; ++i) {
              packForces(outBuffer, *commCells[dir].ghosts[i]);
            }
          }

          // exchange particles, odd-even rule
          if (nodeGrid.getNodePosition(coord) % 2 == 0) {
            outBuffer.send(receiver, DD_COMM_TAG);
            inBuffer.recv(sender, DD_COMM_TAG);
          } else {
            inBuffer.recv(sender, DD_COMM_TAG);
            outBuffer.send(receiver, DD_COMM_TAG);
          }

          // unpack received data
          if (realToGhosts) {
            for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i) {
              unpackPositionsEtc(*commCells[dir].ghosts[i], inBuffer, extradata);
            }
          }
          else {
            for (int i = 0, end = commCells[dir].reals.size(); i < end; ++i) {
              unpackAndAddForces(*commCells[dir].reals[i], inBuffer);
            }
          }
        }
      }
    }
    LOG4ESPP_DEBUG(logger, "ghost communication finished");
  }

  void DomainDecomposition::printCellParticles(Cell* cell, std::string where, bool isReal, std::vector<int> only_particles)
  {
      std::cout << std::setprecision(1) << where << " CELL INFO "
#ifdef CELL_EXTRA_DATA
          << cell->id << " [" << cell->grid_pos << "]"
#endif
          << std::endl;
      for (auto nci : cell->neighborCells){
          std::cout << where << "    NEIGHBOR " 
#ifdef CELL_EXTRA_DATA
          << nci.cell->id << " [" << nci.cell->grid_pos << "]"
#endif
          << ", useForAllPairs " << nci.useForAllPairs << std::endl;
      }

      for (Particle p : cell->particles) {
        if (only_particles.empty() || std::find(only_particles.begin(), only_particles.end(), p.id()) != only_particles.end()){
#ifdef CELL_EXTRA_DATA
           bool in_cell = true;
           for (int i = 0; i < 3; i++)
              if((p.position()[i] < cell->myLeft[i]) || (p.position()[i] > cell->myRight[i]))
                 in_cell = false;
#endif
            std::cout << std::setprecision(3) << std::setw(5) << std::fixed << where
              << ": P " << (isReal ? "REAL  " : "GHOST ") << p.id() << " (" << p.position()
              << std::setprecision(2) << "), V(" << p.velocity() << "), F(" << p.force()
              << std::setprecision(1) << "), MPI " << getSystem()->comm->rank() << "(" << this->nodeGrid.getMyLeft() << "-" << this->nodeGrid.getMyRight() << ")"
#ifdef CELL_EXTRA_DATA
              << ", CELL " << cell->id << "[" << cell->grid_pos << "], (" << cell->myLeft << "<->" << cell->myRight << ")"
              << ", INSIDE " << in_cell
#endif
              << std::endl;
        }
      }
  }

  void DomainDecomposition::printAllParticles(std::string where, std::vector<int> only_particles)
  {
    //std::cout << "PRINTING PARTICLES" << std::endl;
    for (auto cell : realCells){
        printCellParticles(cell, where, true, only_particles);
    }
    for (auto cell : ghostCells){
        printCellParticles(cell, where, false, only_particles);
    }
  }

  void DomainDecomposition::printNodeParticles(int rank)
  {
      if(getSystem()->comm->rank() == rank)
          printAllParticles(std::string("RANK ") + std::to_string(rank), std::vector<int>());
  }



  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void DomainDecomposition::registerPython() {
    using namespace espressopp::python;
    class_< DomainDecomposition, bases< Storage >, boost::noncopyable >
    ("storage_DomainDecomposition", init< shared_ptr< System >,
          //int,
          const Int3D&,
          //const Int3D&,
          const boost::python::list& ,const boost::python::list&, const boost::python::list&,
          const boost::python::list& ,const boost::python::list&, const boost::python::list& >())
    .def("mapPositionToNodeClipped", &DomainDecomposition::mapPositionToNodeClipped)
    .def("getCellGrid", &DomainDecomposition::getInt3DCellGrid)
    .def("getNodeGrid", &DomainDecomposition::getInt3DNodeGrid)
    .def("cellAdjust", &DomainDecomposition::cellAdjust)
    .def("printNodeParticles", &DomainDecomposition::printNodeParticles)
    ;
  }


  }
}
