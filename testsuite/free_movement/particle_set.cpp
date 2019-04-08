#include "storage/Storage.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "particle_set.h"
#include "esutil/Assert.hpp"

//namespace espressopp{
    ParticleSet::ParticleSet(double total_time, double particle_mass)
        : total_time(total_time), particle_mass(particle_mass)
    {

    }

    ParticleSet::~ParticleSet()
    {
    }

    int ParticleSet::num_particles() const
    {
        ESPR_ASSERT(initial_pos.size() == final_pos.size());
        ESPR_ASSERT(initial_pos.size() == velocity.size());

        return initial_pos.size();
    }

    void ParticleSet::add_random_particles(shared_ptr<System> system, int num_add_particles)
    {
        // we are adding on all processors at the moment
        // -> we cannot randomly initialize the generator
//        srand(time(NULL));

//        if(system->comm->rank() == 0) {
            for(int part = 0; part < num_add_particles; part++) {
                Real3D init, final;

                for(int i = 0; i < 3; i++) {                    
                    init[i] = (rand() * system->bc->getBoxL()[i]) / RAND_MAX;
                    final[i] = (rand() * system->bc->getBoxL()[i]) / RAND_MAX;
                }
                add_particle(init, final);
            }
  //      }
    }

    void ParticleSet::add_random_periodic_particles(shared_ptr<System> system, int num_add_particles)
    {
//        if(system->comm->rank() == 0) {
            for(int p = 0; p < num_add_particles; p++) {
                Real3D pos;
                for(int i = 0; i < 3; i++) {
                    pos[i] = (rand() * system->bc->getBoxL()[i]) / RAND_MAX;
                }
                int direction = rand() % 3;
                add_periodic_particle(system, pos, direction);
            }
//        }
    }

    void ParticleSet::add_particle(Real3D init_p, Real3D final_p, Real3D vel)
    {
        initial_pos.push_back(init_p);
        final_pos.push_back(final_p);
        velocity.push_back(vel);
    }

    void ParticleSet::add_particle(Real3D init_p, Real3D final_p)
    {
        add_particle(init_p, final_p, (final_p - init_p) / total_time);
    }

    void ParticleSet::add_periodic_particle(shared_ptr<System> system, Real3D initial_and_final_pos, int direction)
    {
        Real3D vel(0,0,0);
        vel[direction] = system->bc->getBoxL()[direction] / total_time;
        add_particle(initial_and_final_pos, initial_and_final_pos, vel);
    }


    void ParticleSet::register_particles(shared_ptr<System> system)
    {
        // this is of course not a scalable way to do it, just simple....
//        boost::mpi::broadcast(*system->comm, initial_pos, 0);
//        boost::mpi::broadcast(*system->comm, final_pos, 0);
//        boost::mpi::broadcast(*system->comm, velocity, 0);

        for(int partID = 0; partID < num_particles(); partID++)
        {
            Particle* particle = system->storage->addParticle(partID, initial_pos[partID]);
            if(particle)
            {
                std::cout << "really adding, proc " << system->comm->rank() << ", ID: " << partID << ", "
                          << initial_pos[partID] <<  ", vel: " << velocity[partID] << std::endl;
                particle->setV(velocity[partID]);
                if(particle_mass > 0.0)
                        particle->setMass(particle_mass);
            }
        }
    }

    int ParticleSet::check_particles(shared_ptr<System> system) const
    {
        int errors = 0;
        const real tol = 1e-10;
        int particles_found[num_particles()];
        int particles_found_glob[num_particles()];

        for(int partID = 0; partID < num_particles(); partID++)
        {
            particles_found[partID] = particles_found_glob[partID] = 0;

            CellList &cellList = system->storage->getRealCells();
            for(int ci = 0; ci < cellList.size(); ci++) {
                Cell* cell = cellList[ci];
                for(int pi = 0; pi < cell->particles.size(); pi++)
                {
                    Particle part = cell->particles[pi];
                    auto partID = part.getId();
            //                if(part.getId() > num_particles)
            //                    std::cout << "RANK " << system->comm->rank() << ": CHECK: PARTICLE " << part.getId() << std::endl;

                    if(particles_found[partID] > 0)
                        std::cout << "RANK " << system->comm->rank() << ": WRONG: PARTICLE " << partID
                              << " multiplicity even in one processor" << std::endl;
                    particles_found[partID]++;
                    Real3D pos = part.getPos();
                    Real3D diff = pos - final_pos[partID];
                    if (diff.abs() > tol){
                        errors += 1;
                        std::cout <<  "RANK " << system->comm->rank() << ": PARTICLE " << partID
                               << ": WRONG position " << pos << ", expected " << final_pos[partID] << std::endl;
                    }
                }
            }
        }

        boost::mpi::all_reduce(*system->comm, particles_found, num_particles(), particles_found_glob, std::plus<int>());

/// todo: why this does not work?
/// todo: there are sometimes some strangely present but invalid particles!!


//        for(int partID = 0; partID < num_particles; partID++)
//        {
//            particles_found[partID] = particles_found_glob[partID] = 0;
//            if(Particle* particle = system->storage->lookupRealParticle(partID)){
//                if(particle->ghost())
//                    continue;
////                if(particle->getId() != partID)
////                    continue;
//                particles_found[partID] = 1;
//                Real3D pos = particle->getPos();
//                Real3D diff = pos - final_pos[partID];
//                if (diff.abs() > tol){
//                    errors += 1;
//                    std::cout <<  "RANK " << system->comm->rank() << ": PARTICLE " << partID << ", ID: " << particle->getId()
//                              << ": WRONG position " << pos << ", expected " << final_pos[partID] << std::endl;
//                }
//            }

//        }

        boost::mpi::all_reduce(*system->comm, particles_found, num_particles(), particles_found_glob, std::plus<int>());

        for(int partID = 0; partID < num_particles(); partID++)
        {
            if(particles_found_glob[partID] != 1)
                errors += 1;
        }

        return errors;
    }

//}
