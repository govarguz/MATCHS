#!/usr/bin/env python2
#  Copyright (C) 2012-2017(H)
#      Max Planck Institute for Polymer Research
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

###########################################################################
#                                                                         #
#  ESPResSo++ Python script for a Polymer Melt System including           #
#  runtime details
#                                                                         #
###########################################################################

import time
import sys
import espressopp
sys.path.append('../../testsuite/polymer_melt/')
from check import check_results

nsteps          = 20
isteps          = 200
maxcutoff = rc  = 1.12246  #     = pow(2.0, 1.0/6.0)
skin            = 0.3
dt = timestep   = 0.005


# set temperature to None for NVE-simulations
temperature = 1.0


duplicate_int    = int(sys.argv[1])
domain_dec    = sys.argv[2]

halfCellInt_coarse = int(sys.argv[3])
halfCellInt_dense = int(sys.argv[4])

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
def get_positions(storage, num_particles):
    positions = []
    for pid in range(1, num_particles+1):
        part = storage.getParticle(pid)
        positions.append(part.pos)
    return positions

def get_velocities(storage, num_particles):
    velocities = []
    for pid in range(1, num_particles+1):
        part = storage.getParticle(pid)
        velocities.append(part.v)
    return velocities

def get_forces(storage, num_particles):
    forces = []
    for pid in range(1, num_particles+1):
        part = storage.getParticle(pid)
        forces.append(part.f)
    return forces



######################################################################

print espressopp.Version().info()
print 'Setting up simulation ...'
bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.lammps.read('polymer_melt.lammps')
#bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.lammps.read('minisys.lammps')
bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.replicate(bonds, angles, x, y, z, Lx, Ly, Lz, xdim=duplicate_int, ydim=duplicate_int, zdim=duplicate_int)
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
box = (Lx, Ly, Lz)

# <previous> system, integrator = espressopp.standard_system.Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)
# <instead of that, the following taken from free movement>
rng            = espressopp.esutil.RNG()
bc             = espressopp.bc.OrthorhombicBC(rng, box)

system = espressopp.System()
system.bc      = bc
system.rng     = rng
system.skin    = skin
system.expectedMaxCutoff = maxcutoff

nodeGrid       = espressopp.tools.decomp.nodeGrid(box, maxcutoff, skin, espressopp.MPI.COMM_WORLD.size)
#cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, maxcutoff, skin, halfCellInt)

print ("node grid: ", nodeGrid)#, ", cell grid: ", cellGrid)

neiListx,neiListy,neiListz=espressopp.tools.decomp.neiListHom(nodeGrid, box, maxcutoff, skin)  # If no multiscale or Inhomogeneous then use neiListHom

num_coarse = (len(neiListx) - 1) / 2
num_dense = len(neiListx) - 1 - num_coarse
halfCellMaskx = [halfCellInt_coarse] * num_coarse + [halfCellInt_dense] * num_dense
halfCellMasky = [halfCellInt_dense] * (len(neiListy) - 1)
halfCellMaskz = [halfCellInt_dense] * (len(neiListz) - 1)


print("neilist X: ", neiListx)
print("neilist Y: ", neiListy)
print("neilist Z: ", neiListz)
print("halfcellmask X: ", halfCellMaskx)
print("halfcellmask Y: ", halfCellMasky)
print("halfcellmask Z: ", halfCellMaskz)
#ddstorage = espressopp.storage.DomainDecompositionReference(system, nodeGrid, neiListx, neiListy, neiListz,
if(domain_dec == "dd"):
    ddstorage = espressopp.storage.DomainDecomposition(system, nodeGrid, neiListx, neiListy, neiListz,
                                                   halfCellMaskx, halfCellMasky, halfCellMaskz)
elif(domain_dec == "dd_duplicate"):
    ddstorage = espressopp.storage.DomainDecompositionDuplicate(system, nodeGrid, neiListx, neiListy, neiListz,
                                                   halfCellMaskx, halfCellMasky, halfCellMaskz)
elif(domain_dec == "dd_reference"):
    ddstorage = espressopp.storage.DomainDecompositionReference(system, nodeGrid, neiListx, neiListy, neiListz,
                                                   halfCellMaskx, halfCellMasky, halfCellMaskz)
else:
    sys.exit("invalid domain decomposition method")


system.storage = ddstorage

integrator     = espressopp.integrator.VelocityVerlet(system)
integrator.dt  = dt



#   <original polymer melt starts here>
print 'box...',box
# add particles to the system and then decompose
# do this in chunks of 1000 particles to speed it up
props = ['id', 'type', 'mass', 'pos']
new_particles = []
for i in range(num_particles):
  part = [i + 1, 0, 1.0, espressopp.Real3D(x[i], y[i], z[i])]
  new_particles.append(part)
  if i % 1000 == 0:
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()
    new_particles = []
system.storage.addParticles(new_particles, *props)
system.storage.decompose()

# Lennard-Jones with Verlet list
vl      = espressopp.VerletList(system, cutoff = rc)
potLJ   = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)
interLJ = espressopp.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interLJ)

# FENE bonds
fpl = espressopp.FixedPairList(system.storage)
fpl.addBonds(bonds)
#print(bonds)
potFENE = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espressopp.interaction.FixedPairListFENE(system, fpl, potFENE)
system.addInteraction(interFENE)

# Cosine with FixedTriple list
ftl = espressopp.FixedTripleList(system.storage)
ftl.addTriples(angles)
potCosine = espressopp.interaction.Cosine(K=1.5, theta0=3.1415926)
interCosine = espressopp.interaction.FixedTripleListCosine(system, ftl, potCosine)
system.addInteraction(interCosine)

# print simulation parameters
print ''
print 'number of particles = ', num_particles
print 'density             = ', density
print 'rc                  = ', rc
print 'dt                  = ', integrator.dt
print 'skin                = ', system.skin
print 'temperature         = ', temperature
print 'nsteps              = ', nsteps
print 'isteps              = ', isteps
print 'NodeGrid            = ', system.storage.getNodeGrid()
print 'CellGrid            = ', system.storage.getCellGrid()
print ''

# espressopp.tools.decomp.tuneSkin(system, integrator)
PRINT_RANK = 0
#system.storage.printNodeParticles(PRINT_RANK)
initial_verlet_list = vl.getAllPairs()

initial_positions = get_positions(system.storage, num_particles)
espressopp.tools.analyse.info(system, integrator)

start_time = time.clock()
for k in range(nsteps):
  integrator.run(isteps)
  if(k == 4):
      pos_5x200_steps = get_positions(system.storage, num_particles)
      vel_5x200_steps = get_velocities(system.storage, num_particles)
      for_5x200_steps = get_forces(system.storage, num_particles)

  espressopp.tools.analyse.info(system, integrator)

end_time = time.clock()
espressopp.tools.analyse.info(system, integrator)
espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)

pressure = espressopp.analysis.Pressure(system).compute()
ekin = espressopp.analysis.EnergyKin(system).compute()
epot0 = system.getInteraction(0).computeEnergy()
epot1 = system.getInteraction(1).computeEnergy()
epot2 = system.getInteraction(2).computeEnergy()
final_observables = (pressure, ekin, epot0, epot1, epot2)

correct = check_results(initial_verlet_list, initial_positions, pos_5x200_steps, vel_5x200_steps, for_5x200_steps, final_observables)
if(correct):
    print("everything within tolerance")
else:
    print("results not within tolerance")
    raise(Exception("Test failed: tolerance not reached"))
