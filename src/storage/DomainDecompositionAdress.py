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
********************************************
espressopp.storage.DomainDecompositionAdress
********************************************

The DomainDecompositionAdress is the Domain Decomposition for AdResS and H-
AdResS simulations. It makes sure that tuples (i.e. a coarse-grained particle
and its corresponding atomistic particles) are always stored together on one CPU.
When setting DomainDecompositionAdress you have to provide the system as well as
the nodegrid and the cellgrid.

Example - setting DomainDecompositionAdress:

>>> system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid, neiListx, neiListy, neiListz)


.. function:: espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid, neiListx, neiListy, neiListz)

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

"""

from espressopp import pmi
from espressopp.esutil import cxxinit
from _espressopp import storage_DomainDecompositionAdress
from espressopp import toInt3DFromVector
from espressopp.tools import decomp
from espressopp import check

from espressopp.storage.Storage import *

class DomainDecompositionAdressLocal(StorageLocal,storage_DomainDecompositionAdress):

    def __init__(self, system, halfCellInt, nodeGrid, cellGrid, neiListx, neiListy, neiListz):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            p1 = pmi._MPIcomm.rank % nodeGrid[0]
            aux1 =pmi._MPIcomm.rank/nodeGrid[0]  # HDD comment: Getting the order of processors
            p2 = aux1 % nodeGrid[1]
            aux2 = aux1/nodeGrid[1]
            p3 = aux2  # HDD comment: Obtaining the processors per axes (x,y,z)
            cellGrid = Int3D(neiListx[p1+1]-neiListx[p1],neiListy[p2+1]-neiListy[p2],neiListz[p3+1]-neiListz[p3])
            #print "real Cells are",cellGrid
            cxxinit(self, storage_DomainDecompositionAdress, system, halfCellInt, nodeGrid, cellGrid, neiListx, neiListy, neiListz)
if pmi.isController:
    class DomainDecompositionAdress(Storage):
        pmiproxydefs = dict(
            cls = 'espressopp.storage.DomainDecompositionAdressLocal',
            pmicall = ['getCellGrid', 'cellAdjust']
            )
        def __init__(self, system,
                     halfCellInt,
                     nodeGrid='auto',
                     cellGrid='auto',
                     neiListx='auto',
                     neiListy='auto',
                     neiListz='auto'):
            if nodeGrid == 'auto':
                nodeGrid = Int3D(system.comm.rank, 1, 1)
            else:
                nodeGrid = toInt3DFromVector(nodeGrid)

            if cellGrid == 'auto':
                # TODO: Implement
                raise 'Automatic cell size calculation not yet implemented'
            else:
                cellGrid = toInt3DFromVector(cellGrid)
            if neiListx == 'auto':
              neiListx = neiListx
            else:
              neiListx = neiListx
            if neiListy == 'auto':
              neiListy = neiListx
            else:
              neiListy = neiListy
            if neiListz == 'auto':
              neiListz = neiListx
            else:
              neiListz = neiListz
            self.next_id = 0
            self.pmiinit(system, halfCellInt, nodeGrid, cellGrid, neiListx, neiListy, neiListz)
