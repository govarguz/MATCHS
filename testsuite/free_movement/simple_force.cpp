#include "storage/Storage.hpp"
#include "storage/DomainDecomposition.hpp"
#include "storage/DomainDecompositionAdress.hpp"
#include "esutil/RNG.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "main/espressopp_common.hpp"
#include "integrator/MDIntegrator.hpp"
#include "integrator/VelocityVerlet.hpp"
#include "interaction/SingleParticleInteractionTemplate.hpp"
#include "interaction/FixedPairListInteractionTemplate.hpp"
#include "interaction/FixedPairDistListInteractionTemplate.hpp"
#include "FixedPairList.hpp"
#include "FixedPairDistList.hpp"
#include "interaction/GravityTruncated.hpp"

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

class TrivialPotential : public interaction::SingleParticlePotentialTemplate<TrivialPotential>
{
public:
    TrivialPotential(double force) : const_force(force)
    {

    }

    espressopp::real _computeEnergyRaw(const Particle& p, const bc::BC& bc) const
    {
        return 0;
    }

    bool _computeForceRaw(Real3D& force, const Particle& p, const bc::BC& bc) const
    {
        force = Real3D(const_force,0,0);
        return true;
    }

    virtual espressopp::real getMaxCutoff()
    {
        return 1;
    }

private:
    double const_force;

};

class ConstantPairForcePotential : public interaction::PotentialTemplate<ConstantPairForcePotential>
{
public:
    ConstantPairForcePotential(double force) : const_force(force)
    {

    }

    espressopp::real _computeEnergySqrRaw(espressopp::real x ) const
    {
        return 0;
    }

    bool _computeForce(Real3D& force, const Real3D& dist) const
    {
        return _computeForceRaw(force, dist, dist.sqr());
    }

    bool _computeForceRaw(Real3D& force, const Real3D& dist, espressopp::real distSqr) const
    {
        force = - const_force * dist / sqrt(distSqr);
    }

//    virtual espressopp::real getMaxCutoff()
//    {
//        return 1;
//    }

private:
    double const_force;
};

int main(int argc, char **argv)
{
    initMPIEnv(argc, argv);

    const int num_steps = 100;
    const double total_time = 2.0;
    const double skin = 0.1;
    const double const_force = 2.5;
    const double particle_mass = 10;
    const double time_step = total_time / num_steps;

    std::vector<int> neiListx = {0,3};
    std::vector<int> neiListy = {0,3,6,9,12,15};
    std::vector<int> neiListz = {0,3};

    std::vector<int> halfCellx = {1};
    std::vector<int> halfCelly = {1,1,1,1,1};
    std::vector<int> halfCellz = {1};
//    std::vector<int> neiListx = {0,3,6,9};
//    std::vector<int> neiListy = {0,3,6,9,12};
//    std::vector<int> neiListz = {0,3,6,9};

//    std::vector<int> halfCellx = {1,1,1};
//    std::vector<int> halfCelly = {1,1,1,1};
//    std::vector<int> halfCellz = {1,1,1};

    Int3D nodeGrid(1,5,1);


    boost::shared_ptr<class System > system(new System());
    boost::shared_ptr<class esutil::RNG > rng(new esutil::RNG());
    system->rng = rng;
    Real3D boxL(5., 5., 5.);
    boost::shared_ptr<class bc::OrthorhombicBC> bc(new bc::OrthorhombicBC(rng, boxL));
    system->bc = bc;
    system->setSkin(skin);

    auto newcomm = boost::make_shared< mpi::communicator >();
    system->comm = newcomm;

    auto storage = boost::make_shared<storage::DomainDecomposition>(system, //halfCellInt,
                                                   nodeGrid, //cellGrid,
                                                   neiListx, neiListy, neiListz,
                                                   halfCellx, halfCelly, halfCellz);

    system->storage = storage;

    double a = const_force / particle_mass; //(F/m)
    double pos_change_by_force = 0.5 * a * total_time * total_time;
    std::cout << "Each particle should be moved by force in x and z directions by " << pos_change_by_force << std::endl;

    ParticleSet particle_set(total_time, particle_mass);
    // in x, postition of particles is changed by TrivialPotential, force allways acts in direction (1,0,0)
    // in y, position changed by free movement, velocity is specified as third parameter of add_particle
    // in z, position changed by force by which 2 particles apply to each other. It is in z direction, since
    // booth particles have the same x and y positions and its strenght is independent of the distance of the particles,
    // as defined in ConstantPairForcePotential
    particle_set.add_particle(Real3D(2.0, 0.5, 1.5), Real3D(2.0 + pos_change_by_force, 4.5, 1.5 + pos_change_by_force), Real3D(0, 4./total_time, 0));
    particle_set.add_particle(Real3D(2.0, 0.5, 3.5), Real3D(2.0 + pos_change_by_force, 4.5, 3.5 - pos_change_by_force), Real3D(0, 4./total_time, 0));
    particle_set.register_particles(system);

    auto trivial_potential = boost::make_shared<TrivialPotential>(const_force);
    auto direction_force = boost::make_shared<interaction::SingleParticleInteractionTemplate<TrivialPotential> >(system, trivial_potential);

    boost::shared_ptr<FixedPairList> pair_list = boost::make_shared<FixedPairList>(storage);
    pair_list->add(0,1);
    auto constant_pair_force_potential = boost::make_shared<ConstantPairForcePotential>(const_force);
    auto constant_pair_force =
            boost::make_shared<interaction::FixedPairListInteractionTemplate<ConstantPairForcePotential> >
            (system, pair_list, constant_pair_force_potential);

    system->addInteraction(constant_pair_force);
    system->addInteraction(direction_force);

    boost::shared_ptr<class integrator::VelocityVerlet> integrator(new integrator::VelocityVerlet(system));
    integrator->setTimeStep(time_step);

    //Probe probe(system, integrator, 0);

    integrator->run(num_steps);
//    std::cout << "rank " << system->comm->rank() << ", num local particles " << storage->getNLocalParticles()
//              << ", real part " << storage->getNLocalParticles() - storage->getNGhostParticles()
//              << ", ghost part " << storage->getNGhostParticles() << std:: endl;
    int num_errors = particle_set.check_particles(system);
    finalizeMPIEnv();
    std::cout << "num errors " << num_errors << std::endl;
    return num_errors;
}
