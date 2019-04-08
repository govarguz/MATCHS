#include "storage/Storage.hpp"
#include "storage/DomainDecomposition.hpp"
#include "storage/DomainDecompositionAdress.hpp"
#include "esutil/RNG.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "main/espressopp_common.hpp"
#include "integrator/MDIntegrator.hpp"
#include "integrator/VelocityVerlet.hpp"

#include "particle_set.h"

using namespace espressopp;


class Probe : public SystemAccess
{
public:

    Probe(boost::shared_ptr<System> system, boost::shared_ptr<integrator::MDIntegrator> _integrator, int _particle_id) :
        SystemAccess(system), integrator(_integrator), particle_id(_particle_id)
    {
        conn_init = integrator->runInit.connect(
                boost::bind(&Probe::probe_particle, this));
//        conn_recalc1 = integrator->recalc1.connect(
//                boost::bind(&Probe::probe_particle, this));
        conn_aftIntV = integrator->aftIntV.connect(
                boost::bind(&Probe::probe_particle, this));
    }

    void probe_particle()
    {
//        std::cout << "   PROBE   " << std::endl;
        CellList &cellList = getSystem()->storage->getRealCells();
        for(int ci = 0; ci < cellList.size(); ci++) {
            Cell* cell = cellList[ci];
//            for(auto cid : cell->neighborCells)
//                std::cout << cid << ", ";
            for(int pi = 0; pi < cell->particles.size(); pi++)
            {
                Particle part = cell->particles[pi];
                if(part.getId() != particle_id)
                    continue;

                std::cout << "PROBE (" << getSystem()->comm->rank() << ") time step " << integrator->getStep() << ", REAL CELL ID " << ci
                          << ", position " << part.getPos() << std::endl;

            }
        }

    }

    boost::signals2::connection conn_init, conn_recalc1, conn_aftIntV;
    boost::shared_ptr<integrator::MDIntegrator> integrator;
    int particle_id;
};

int main(int argc, char **argv)
{
    initMPIEnv(argc, argv);


    const int num_steps = 20;
    const double total_time = 2.0;
    const double time_step = total_time / num_steps;

    const int num_particles = 1000;
    const double skin = 0.1;

    std::vector<int> neiListx = {0,3,6,9};
    std::vector<int> neiListy = {0,3,6,9,12};
    std::vector<int> neiListz = {0,3,6,9};

    std::vector<int> halfCellx = {1,2,1};
    std::vector<int> halfCelly = {1,2,2,1};
    std::vector<int> halfCellz = {1,2,1};

    Int3D nodeGrid(1,2,1);


    boost::shared_ptr<class System > system(new System());
    boost::shared_ptr<class esutil::RNG > rng(new esutil::RNG());
    system->rng = rng;
    Real3D boxL(5., 5., 5.);
    boost::shared_ptr<class bc::OrthorhombicBC> bc(new bc::OrthorhombicBC(rng, boxL));
    system->bc = bc;
    system->setSkin(skin);

    boost::shared_ptr< mpi::communicator > newcomm = boost::make_shared< mpi::communicator >();
    system->comm = newcomm;

    boost::shared_ptr<storage::DomainDecomposition> storage(
                  new storage::DomainDecomposition(system, //halfCellInt,
                                                   nodeGrid, //cellGrid,
                                                   neiListx, neiListy, neiListz,
                                                   halfCellx, halfCelly, halfCellz));
    system->storage = storage;

    ParticleSet particle_set(total_time);
    particle_set.add_random_particles(system, num_particles / 2);
//      particle_set.add_particle(Real3D(3.27684, 0.905806, 1.52619), Real3D(1.30785, 4.28778, 1.70677));

//    particle_set.add_particle(Real3D(2.5, 0.5, 2.5), Real3D(2.5, 4.5, 2.5));
//    particle_set.add_periodic_particle(system, Real3D(2.5, 0.5, 2.5), 1);
//    particle_set.add_random_periodic_particles(system, max_num_particles / 2);
    particle_set.register_particles(system);

    shared_ptr<class integrator::VelocityVerlet> integrator(new integrator::VelocityVerlet(system));
    integrator->setTimeStep(time_step);

//    Probe probe(system, integrator, 0);

    integrator->run(num_steps);
//    std::cout << "rank " << system->comm->rank() << ", num local particles " << storage->getNLocalParticles()
//              << ", real part " << storage->getNLocalParticles() - storage->getNGhostParticles()
//              << ", ghost part " << storage->getNGhostParticles() << std:: endl;
    int num_errors = particle_set.check_particles(system);
    finalizeMPIEnv();
    std::cout << "num errors " << num_errors << std::endl;
    return num_errors;
}
