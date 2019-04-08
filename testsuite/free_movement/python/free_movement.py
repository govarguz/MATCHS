import sys
import espressopp
from espressopp import Real3D

#script paramters: halfCellInt_outer, halfCellInt_center
#script is designed to run with 27 mpi tasks and prescribe the 
#middle processor of the 3x3x3 grid different (higher) halfCellInt
#than to the processors around

halfCellInt_outer    = int(sys.argv[1])
halfCellInt_center    = int(sys.argv[2])
if(halfCellInt_outer > halfCellInt_center):
    raise ValueError("incorrect setup")

print("got parameters half cell outer " + str(halfCellInt_outer))

L              = 10.0
maxcutoff      = 0.1
#maxcutoff      = 0.2 * pow(2.0, 1.0/6.0)
skin           = 0.4

num_particles = 1000

dt             = 0.005
nsteps         = 400

TOL            = 1e-5

# creates particles with initial random position and velocity 
# expected end position is calculated from that
# assuming periodic boundary conditions, thus if particle should fly out of the box,
# it is assumed that it returns from the other side
def random_particles_periodic(system, num_part, time, L):
  particle_list = []
  for pid in range(1, num_part + 1):
    pos  = system.bc.getRandomPos()
    v    = Real3D(system.rng(), system.rng(), system.rng())

    end_pos = pos + time * v
    for dim in range(3):
      while(end_pos[dim] > L):
        end_pos[dim] -= L
      while(end_pos[dim] < 0.0):
        end_pos[dim] += L

    part = [pid, pos, v, end_pos]
    particle_list.append(part)

  return particle_list

# creates particles with random initial and end position within the box
# the velocity is then calculated from that
# In this case, thus, particle is not supposed to leave the box
def random_particles_no_leave(system, num_part, time, L):
  particle_list = []
  for pid in range(1, num_part + 1):
    pos     = system.bc.getRandomPos()
    end_pos = system.bc.getRandomPos()

    v = (end_pos - pos) / time

    part = [pid, pos, v, end_pos]
    particle_list.append(part)

  return particle_list

def single_particle(num_part, time):
  pos = Real3D(6.2171742242955474, 2.097406695157744, 5.632429607545326)
  end_pos = Real3D(2.960239482582736, 7.745615077405823, 0.6471872473406748)

  #v = Real3D(-6.513869483425623, 11.296416764496158, -9.970484720409303)
  v = (end_pos - pos) / time

  particle_list = []
  for pid in range(1, num_part + 1):
    particle_list.append([pid, pos, v, end_pos])
  return particle_list

def check_particles(system, my_particle_list, num_part):
  for pid in range(1, num_part + 1):
    p = system.storage.getParticle(pid)
    if(p == None):
        continue
    #expected = my_particle_list[p.id-1][-1]
    expected = my_particle_list[pid-1][-1]
    #print("checking particle pid ", pid, ", internal id ", p.id, " ghost ", p.isGhost)

    if((expected - p.pos).abs() > TOL):
      print("Wrong position of particle ", p.id, ", position: ", p.pos, ", expected: ", expected)
      print("PARTICLE PROPERTIES: ", my_particle_list[pid-1])
      raise(Exception("does not match"))
  print("Tested " + str(num_part) +" particles, all OK")

box            = (L,L,L)
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
halfCellMaskx = [halfCellInt_outer] * (len(neiListx) - 1)
halfCellMasky = [halfCellInt_outer] * (len(neiListy) - 1)
halfCellMaskz = [halfCellInt_outer] * (len(neiListz) - 1)

halfCellMaskx[1] = halfCellInt_center
halfCellMasky[1] = halfCellInt_center
halfCellMaskz[1] = halfCellInt_center

print("neilist X: ", neiListx)
print("neilist Y: ", neiListy)
print("neilist Z: ", neiListz)
print("halfcellmask X: ", halfCellMaskx)
print("halfcellmask Y: ", halfCellMasky)
print("halfcellmask Z: ", halfCellMaskz)
#ddstorage = espressopp.storage.DomainDecompositionReference(system, nodeGrid, neiListx, neiListy, neiListz,
ddstorage = espressopp.storage.DomainDecompositionDuplicate(system, nodeGrid, neiListx, neiListy, neiListz,
                                                   halfCellMaskx, halfCellMasky, halfCellMaskz)

system.storage = ddstorage

integrator     = espressopp.integrator.VelocityVerlet(system)
integrator.dt  = dt

my_particle_list = random_particles_periodic(system, num_particles, dt*nsteps, L)
#my_particle_list = random_particles_no_leave(system, num_particles, dt*nsteps, L)
#my_particle_list = single_particle(num_particles, dt*nsteps)

masses = [system.rng() for i in range(num_particles)]
types = [system.rng(2) for i in range(num_particles)]
particle_list = [[pid, pos, type, v, mass] for ([pid, pos, v, end_pos], mass, type) in zip(my_particle_list, masses, types)]
system.storage.addParticles(particle_list, 'id', 'pos', 'type', 'v', 'mass')

system.storage.decompose()

integrator.run(nsteps)

#for i in range(nsteps):
#    integrator.run(1)
#    for pid in range(1, num_particles + 1):
#      p = system.storage.getParticle(pid)
#      print("POSITION ", pid, ", internal id ", p.id, "position: ", p.pos, "v: ", p.v)


check_particles(system, my_particle_list, num_particles)
