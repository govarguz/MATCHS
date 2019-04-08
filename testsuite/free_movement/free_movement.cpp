#include "storage/Storage.hpp"
#include "storage/DomainDecomposition.hpp"
#include "storage/DomainDecompositionAdress.hpp"
#include "esutil/RNG.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "main/espressopp_common.hpp"
#include "integrator/VelocityVerlet.hpp"

#include "particle_set.h"

using namespace espressopp;

int main(int argc, char **argv)
{
    initMPIEnv(argc, argv);

    const int num_steps = 20;
    const double total_time = 2.0;
    const double time_step = total_time / num_steps;
    const int num_particles = 1000;

    std::vector<int> neiListx, neiListy, neiListz;
    neiListx.push_back(0); neiListx.push_back(4); neiListx.push_back(9);
    neiListy.push_back(0); neiListy.push_back(3);  neiListy.push_back(10); neiListy.push_back(15);
    neiListz.push_back(0); neiListz.push_back(5); neiListz.push_back(8);

    Int3D nodeGrid(2,3,2);

    int halfCellInt = 1;
    boost::shared_ptr<class System > system(new System());
    boost::shared_ptr<class esutil::RNG > rng(new esutil::RNG());
    system->rng = rng;
    Real3D boxL(5., 5., 5.);
    boost::shared_ptr<class bc::OrthorhombicBC> bc(new bc::OrthorhombicBC(rng, boxL));
    system->bc = bc;

    boost::shared_ptr< mpi::communicator > newcomm = boost::make_shared< mpi::communicator >();
    system->comm = newcomm;

    boost::shared_ptr<storage::DomainDecomposition> storage(
                  new storage::DomainDecomposition(system, halfCellInt,
                                                   nodeGrid, //cellGrid,
                                                   neiListx, neiListy, neiListz));
    system->storage = storage;

    ParticleSet particle_set(total_time);
    particle_set.add_random_particles(system, num_particles);
    particle_set.register_particles(system);

    integrator::VelocityVerlet integrator(system);
    integrator.setTimeStep(time_step);

    integrator.run(num_steps);

    //std::cout << "rank " << system->comm->rank() << ", num local particles " << storage->getNLocalParticles() << ", ghost part " << storage->getNGhostParticles() << std:: endl;
    int num_errors = particle_set.check_particles(system);
    finalizeMPIEnv();
    std::cout << "num errors " << num_errors << std::endl;
    return num_errors;
}
