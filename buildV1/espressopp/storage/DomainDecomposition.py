#  Copyright (C) 2012,2013,2017(H)
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


r"""
**************************************
espressopp.storage.DomainDecomposition
**************************************


.. function:: espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid,,neiListx,neiListy,neiListz)

		:param system: 
		:param nodeGrid: 
                :param cellGrid:
                :param neiListx:
                :param neiListy:
                :param neiListz:
		:type system: 
                :type Real3D:
                :type Real3D:
                :type boost:vector:
                :type boost:vector:
                :type boost:vector:

.. function:: espressopp.storage.DomainDecomposition.getCellGrid()

		:rtype: 

.. function:: espressopp.storage.DomainDecomposition.getNodeGrid()

		:rtype: 
"""
from espressopp import pmi
from espressopp.esutil import cxxinit
from _espressopp import storage_DomainDecomposition
from espressopp import Int3D, toInt3DFromVector
from espressopp.tools import decomp
from espressopp import check
import mpi4py.MPI as MPI

from espressopp.storage.Storage import *

class DomainDecompositionLocal(StorageLocal, storage_DomainDecomposition):
    'The (local) DomainDecomposition.'
    def __init__(self, system, nodeGrid, neiListx, neiListy, neiListz, halfCellMaskx, halfCellMasky, halfCellMaskz):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            p1 = pmi._MPIcomm.rank % nodeGrid[0]
            aux1 =pmi._MPIcomm.rank/nodeGrid[0]  # HeSpaDDA comment: Refers to the order in which processors are given within the Linked-Cell-List
            p2 = aux1 % nodeGrid[1]
            aux2 = aux1/nodeGrid[1]
            p3 = aux2  # HeSpaDDA comment: The processors triplet (p1,p2,p3) have been extracted and are ready for the construction of the cells neighbor list
            cxxinit(self, storage_DomainDecomposition, system, nodeGrid,  neiListx, neiListy, neiListz,
                    halfCellMaskx, halfCellMasky, halfCellMaskz)

    def getCellGrid(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getCellGrid(self)

    def getNodeGrid(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getNodeGrid(self)

    def printNodeParticles(self, rank):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.printNodeParticles(self, rank)

if pmi.isController:
    class DomainDecomposition(Storage):
        pmiproxydefs = dict(
          cls = 'espressopp.storage.DomainDecompositionLocal',  
          pmicall = ['getCellGrid', 'getNodeGrid', 'cellAdjust', 'printNodeParticles']
        )
        def __init__(self, system,
                     nodeGrid='auto', 
                     neiListx='auto',
                     neiListy='auto',
                     neiListz='auto',
                     halfCellMaskx='auto',
                     halfCellMasky='auto',
                     halfCellMaskz='auto',
                     nocheck=False):
            # do sanity checks for the system first
            if nocheck:
              self.next_id = 0
              self.pmiinit(system,  nodeGrid, neiListx, neiListy, neiListz,
                           halfCellMaskx, halfCellMasky, halfCellMaskz)
            else:
              if check.System(system, 'bc'):
                if nodeGrid == 'auto':
                  nodeGrid = decomp.nodeGrid(system.comm.rank)
                else:
                  nodeGrid = toInt3DFromVector(nodeGrid)
#                if cellGrid == 'auto':
#                  cellGrid = Int3D(2,2,2)
#                else:
#                  cellGrid = cellGrid
                if neiListx == 'auto':
                  neiListx = neiListx
                else:
                  neiListx = neiListx
                if neiListy == 'auto':
                  neiListy = neiListy
                else:
                  neiListy = neiListy
                if neiListz == 'auto':
                  neiListz = neiListz
                else:
                  neiListz = neiListz
                if halfCellMaskx == 'auto':
                  halfCellMaskx = [1] * (len(neiListx) - 1)
                if halfCellMasky == 'auto':
                  halfCellMasky = [1] * (len(neiListy) - 1)
                if halfCellMaskz == 'auto':
                  halfCellMaskz = [1] * (len(neiListz) - 1)
                # minimum image convention check:
#                for k in range(3):
#                  if nodeGrid[k]*cellGrid[k]== 1 :
#                    print(("Warning! cellGrid[{}] has been "
#                             "adjusted to 2 (was={})".format(k, cellGrid[k])))
#                    cellGrid[k] = 2
                self.next_id = 0
                self.pmiinit(system, nodeGrid, neiListx, neiListy, neiListz,
                             halfCellMaskx, halfCellMasky, halfCellMaskz)
              else:
                print 'Error: could not create DomainDecomposition object'
