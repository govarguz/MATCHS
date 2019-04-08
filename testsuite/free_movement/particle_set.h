#ifndef PARTICLE_SET_H
#define PARTICLE_SET_H

#include "include/types.hpp"
#include "main/espressopp_common.hpp"
#include "storage/Storage.hpp"

using namespace espressopp;

//namespace espressopp{
    struct ParticleSet
    {
        ParticleSet(double total_time, double particle_mass = 0.0);
        ~ParticleSet();

        // velocity calculated from initial and final position, assuming, that there is no force
        // applied and thus the motion is linear
        void add_particle(Real3D initial_pos, Real3D final_pos);

        // specify all initial position, final position and velocity
        void add_particle(Real3D initial_pos, Real3D final_pos, Real3D velocity);

        void add_random_particles(boost::shared_ptr<class System > system, int num_add_particles);

        // add particle, that starts in a position, goes in direction of one of the axis, goes
        // through the (periodic) boundary and returns to its initial position
        void add_periodic_particle(boost::shared_ptr<class System > system, Real3D initial_and_final_pos, int direction);
        void add_random_periodic_particles(boost::shared_ptr<class System > system, int num_add_particles);

        void register_particles(boost::shared_ptr<System> system);
        int check_particles(boost::shared_ptr<System> system) const;

        int num_particles() const;

        std::vector<Real3D> initial_pos, final_pos, velocity;
        double total_time;
        double particle_mass;
    };
//}

#endif // PARTICLE_SET_H
