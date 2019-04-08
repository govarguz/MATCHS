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
*************************************************
espressopp.storage.DomainDecompositionNonBlocking
*************************************************


.. function:: espressopp.storage.DomainDecompositionNonBlocking(system, nodeGrid, cellGrid, neiListx, neiListy, neiListz)

		:param system: 
		:param nodeGrid: 
		:param cellGrid: 
		:type system: 
                :type Real3D:
                :type real3D:
                :type boost:vector:
                :type boost:vector:
                :type boost:vector:
"""
from espressopp import pmi
from espressopp.esutil import cxxinit
from _espressopp import storage_DomainDecomposition
from _espressopp import storage_DomainDecompositionNonBlocking
from espressopp import Int3D, toInt3DFromVector
import mpi4py.MPI as MPI

from espressopp.storage.DomainDecomposition import *

class DomainDecompositionNonBlockingLocal(DomainDecompositionLocal, storage_DomainDecompositionNonBlocking):
    'The (local) DomainDecompositionNonBlocking.'
    def __init__(self, system, nodeGrid, cellGrid, neiListx, neiListy, neiListz):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            p1 = pmi._MPIcomm.rank % nodeGrid[0]
            aux1 =pmi._MPIcomm.rank/nodeGrid[0]  # HDD comment: Getting the order of processors
            p2 = aux1 % nodeGrid[1]
            aux2 = aux1/nodeGrid[1]
            p3 = aux2  # HDD comment: Obtaining the processors per axes (x,y,z)
            cellGrid = Int3D(neiListx[p1+1]-neiListx[p1],neiListy[p2+1]-neiListy[p2],neiListz[p3+1]-neiListz[p3])
            cxxinit(self, storage_DomainDecompositionNonBlocking, system, nodeGrid, cellGrid, neiListx, neiListy, neiListz)
if pmi.isController:
    class DomainDecompositionNonBlocking(DomainDecomposition):
        pmiproxydefs = dict(
          cls = 'espressopp.storage.DomainDecompositionNonBlockingLocal'
        )
        def __init__(self, system,
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
              neiListx =0
            else:
              neiListx = neiListx
            if neiListy == 'auto':
              neiListy =0
            else:
              neiListy = neiListy
            if neiListz == 'auto':
              neiListz =0
            else:
              neiListz = neiListz
            self.next_id = 0
            self.pmiinit(system, nodeGrid, cellGrid, neiListx, neiListy, neiListz)
