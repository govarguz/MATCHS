#!/usr/bin/env python                                                               
# -*- coding: iso-8859-1 -*-                                                        

import time
import espressopp
import logging
from mpi4py import MPI
import sys

# arguments -- halfCell in coarse and dense region
halfCellInt_coarse = int(sys.argv[1])
halfCellInt_dense = int(sys.argv[2])

# Please modify follwoing 4 parameters to create your system
density1           = 0.85
density2           = 0.45
L                  = 40
monomers_per_chain = 50
##############

num_chains1        = int(density1*L**3/monomers_per_chain)
num_chains2        = int(density2*L**3/monomers_per_chain)
num_chains         = num_chains1 + num_chains2

# set temperature to None for NVE-simulations
temperature = 1.0
seed        = 6543215 # seed for random

a1 = 2
a2 = 1
a3 = 1

# set the potential for fine-grained polymer properties
# LJ
epsilon = 1.0
sigma   = 1.0
rc_lj   = pow(2.0, 1.0/6.0)
# FENE
K_fene    = 30.0
r0_fene   =  0.0
rmax_fene =  1.5
rc_fene   = 1.6

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################

spline  = 2                                # spline interpolation type (1, 2, 3)

tabfileLJ = "pot-lj.txt"
tabfileFENE = "pot-fene.txt"

# writes the tabulated file
def writeTabFile(pot, name, N, low=0.0, high=2.5, body=2):
    outfile = open(name, "w")
    delta = (high - low) / (N - 1)
     
    for i in range(N):
        r = low + i * delta
        energy = pot.computeEnergy(r)
        if body == 2:# this is for 2-body potentials
            force = pot.computeForce(espressopp.Real3D(r, 0.0, 0.0))[0]
        else: # this is for 3- and 4-body potentials
            force = pot.computeForce(r)
        outfile.write("%15.8g %15.8g %15.8g\n"%(r, energy, force))
     
    outfile.close()

nsteps      = 200
isteps      = 200
rc          = max(rc_lj, rc_fene)
maxcutoff   = rc
skin        = 0.3
timestep    = 0.005
box         = (a1*L, a2*L, a3*L)



#random.seed(seed)

print espressopp.Version().info()
print 'Setting up simulation ...'

#logging.getLogger("SteepestDescent").setLevel(logging.INFO)

system         = espressopp.System()
system.rng     = espressopp.esutil.RNG()
system.rng.seed(seed)
system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin
system.expectedMaxCutoff = maxcutoff

nodeGrid       = espressopp.tools.decomp.nodeGrid(box, maxcutoff, skin, espressopp.MPI.COMM_WORLD.size)
#cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
#nodeGrid = (4, 1, 1)

print ("node grid: ", nodeGrid) #, ", cell grid: ", cellGrid)

neiListx,neiListy,neiListz=espressopp.tools.decomp.neiListHom(nodeGrid, box, maxcutoff, skin)  # If no multiscale or Inhomogeneous then use neiListHom
num_coarse = (len(neiListx) - 1) / 2
num_dense = len(neiListx) - 1 - num_coarse
halfCellMaskx = [halfCellInt_coarse] * num_coarse + [halfCellInt_dense] * num_dense
halfCellMasky = [halfCellInt_dense] * (len(neiListy) - 1)
halfCellMaskz = [halfCellInt_dense] * (len(neiListz) - 1)

#
print("neilist X: ", neiListx)
print("neilist Y: ", neiListy)
print("neilist Z: ", neiListz)
print("halfcellmask X: ", halfCellMaskx)
print("halfcellmask Y: ", halfCellMasky)
print("halfcellmask Z: ", halfCellMaskz)

#   ddstorage = espressopp.storage.DomainDecomposition(system, nodeGrid, neiListx, neiListy, neiListz,
#                                                    halfCellMaskx, halfCellMasky, halfCellMaskz)

ddstorage = espressopp.storage.DomainDecompositionDuplicate(system, nodeGrid, neiListx, neiListy, neiListz,
                                                   halfCellMaskx, halfCellMasky, halfCellMaskz)

print("ddstorage created")
system.storage = ddstorage
print("ddstorage assigned")
# --- end

#system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

integrator     = espressopp.integrator.VelocityVerlet(system)
integrator.dt  = timestep
thermostat     = espressopp.integrator.LangevinThermostat(system)
thermostat.gamma  = 0.5
thermostat.temperature = temperature
integrator.addExtension(thermostat)

# set the polymer properties
bondlen            = 0.97

props    = ['id', 'type', 'mass', 'pos', 'v']
vel_zero = espressopp.Real3D(0.0, 0.0, 0.0)

bondlist  = espressopp.FixedPairList(system.storage)
#anglelist = espressopp.FixedTripleList(system.storage)
pid      = 1
type     = 0
mass     = 5.0

res_file = open('microscopic_nb1_no40.res')
# add particles to the system and then decompose
# do this in chunks of 1000 particles to speed it up
chain = []
for i in range(num_chains):
  #print "Check_Num", i
  #if i < num_chains/2 or i > num_chains/2 + 100 - 1:
  if i >= num_chains1:
    mass = 1.0
  startpos = system.bc.getRandomPos()
  positions, bonds, angles = espressopp.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen, True)
  j = 0
  while j < monomers_per_chain:
    line = res_file.readline()
    parameters = line.split()
    i_diff = 0
    if (len(parameters) < 15): # originaly 12
      i_diff = 1
    if parameters[0] == "ATOM":
      res_positions = espressopp.Real3D((float(parameters[6 - i_diff]) + a1*L)%(a1*L),
                    	                (float(parameters[7 - i_diff]) + a2*L)%(a2*L),
                                        (float(parameters[8 - i_diff]) + a3*L)%(a3*L))
      print pid + j, mass, res_positions
      part = [pid + j, type, mass, res_positions, vel_zero]
      chain.append(part)
      j += 1
  pid += monomers_per_chain
  system.storage.addParticles(chain, *props)
  system.storage.decompose()
  chain = []
  bondlist.addBonds(bonds)
  #anglelist.addTriples(angles)
system.storage.addParticles(chain, *props)
system.storage.decompose()

num_particles = num_chains * monomers_per_chain
density = num_particles * 1.0 / (L * L * L * a1 * a2 * a3)

# Lennard-Jones with Verlet list
print "#set LJ"
vl      = espressopp.VerletList(system, cutoff = rc_lj)
potLJ   = espressopp.interaction.LennardJones(epsilon, sigma, cutoff=rc_lj, shift=0)
interLJ = espressopp.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
#print 'Generating potential files ... (%2s)\n' % (tabfileLJ)
#writeTabFile(potLJ, tabfileLJ, N=257, low=0.01, high=rc_lj)
#potTabLJ = espressopp.interaction.Tabulated(itype=spline, filename=tabfileLJ, cutoff=rc_lj)
#interLJ = espressopp.interaction.VerletListTabulated(vl)
#interLJ.setPotential(type1=0, type2=0, potential=potTabLJ)
system.addInteraction(interLJ)

# FENE bonds
print "#set FENE"
potFENE = espressopp.interaction.FENECapped(K=K_fene, r0=r0_fene, rMax=rmax_fene, cutoff=rc_fene, caprad=1.4999)
interFENE = espressopp.interaction.FixedPairListFENECapped(system, bondlist, potFENE)
#print 'Generating potential files ... (%2s)\n' % (tabfileFENE)
#writeTabFile(potFENE, tabfileFENE, N=513, low=0.0001, high=potFENE.cutoff)
#potTabFENE = espressopp.interaction.Tabulated(itype=spline, filename=tabfileFENE)
#interFENE = espressopp.interaction.FixedPairListTabulated(system, bondlist, potTabFENE)
system.addInteraction(interFENE)

# Cosine with FixedTriple list
#potCosine = espressopp.interaction.Cosine(K=1.5, theta0=0.)
#interCosine = espressopp.interaction.FixedTripleListCosine(system, anglelist, potCosine)
#system.addInteraction(interCosine)

filename = "polymer_melt.pdb"
#espressopp.tools.pdb.pdbwrite(filename, system, monomers_per_chain, False)
#espressopp.tools.pdb.fastwritepdb(filename, system, monomers_per_chain, False)

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
#print 'CellGrid            = ', system.storage.getCellGrid()
print ''

# espressopp.tools.decomp.tuneSkin(system, integrator)

start_time = time.clock()
for k in range(nsteps):
  integrator.run(isteps)
  #integrator.run(1)
  espressopp.tools.analyse.info(system, integrator)
  #espressopp.tools.pdb.fastwritepdb(filename, system, monomers_per_chain, True)
end_time = time.clock()
#espressopp.tools.analyse.info(system, integrator)
espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)
espressopp.tools.pdb.pdbwrite(filename, system, monomers_per_chain, False)
