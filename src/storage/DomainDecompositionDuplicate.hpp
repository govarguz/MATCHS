#ifndef _STORAGE_DOMAINDECOMPOSITIONDUPLICATE_HPP
#define _STORAGE_DOMAINDECOMPOSITIONDUPLICATE_HPP

#include <storage/DomainDecomposition.hpp>

namespace espressopp{
  namespace storage{
    template <typename DD>
    class DomainDecompositionDuplicateTemplate : public DD
    {
      public:
      DomainDecompositionDuplicateTemplate(shared_ptr< System > system,
              const Int3D& _nodeGrid,
              const boost::python::list& neiListx,
              const boost::python::list& neiListy,
              const boost::python::list& neiListz,
              const boost::python::list& halfCellMaskx,
              const boost::python::list& halfCellMasky,
              const boost::python::list& halfCellMaskz);


      // follow the routines, which have to be redefined
      // they perform communication, and thus they will be performed using dd_comm
      //--------------------------------------------------------
      public:
      virtual void updateGhosts();
      virtual void updateGhostsV();
      virtual void collectGhostForces();
      protected:
      virtual void decomposeRealParticles();
      virtual void exchangeGhosts();
      //--------------------------------------------------------

      protected:

      void transRealParticlesBeforeDecompose();
      void transRealParticlesAfterDecompose();

      void transGhostParticlesBeforeExchangeGhosts();
      void transGhostParticlesAfterExchangeGhosts();

      // TODO: reals could be optimized: only cells close to the border have to be dealt with
      void copyParticlePropertiesToComm(CellList computeSource, CellList commTarget, std::vector<Cell*> cellOrder, bool positions);
      void copyParticlePropertiesToCompute(CellList commSource, bool positions);

      // transfer REAL particles from COMPUTE to COMM layer before calling update ghosts
      void transParticlesBeforeUpdate();
      // transfer ALL particles from COMPUTE to COMM layer before calling collect ghost forces
      void transParticlesBeforeCollect();
      // transfer GHOST particles from COMM to COMPUTE layer after calling update ghosts or collect ghost forces
      void transParticlesAfterUpdate();
      // transfer REAL particles from COMM to COMPUTE layer after calling update ghosts or collect ghost forces
      void transParticlesAfterCollect();

      void check_layer_consistency(std::string msg);

      shared_ptr<DD> dd_comm;

      // all this is extra data and it would be nicer, if we could get rid of that
      std::map<Cell*, Cell*> fine_to_coarse;
      std::map<size_t, Cell*> cell_in_dd_comm;
      std::map<size_t, Cell*> cell_in_dd_compute;
      std::vector<Cell*> order_of_real_cells;
      std::vector<Cell*> order_of_ghost_cells;

    };

    class DomainDecompositionDuplicate : public DomainDecompositionDuplicateTemplate<DomainDecomposition>
    {
      public:

      DomainDecompositionDuplicate(shared_ptr< System > system,
              const Int3D& _nodeGrid,
              const boost::python::list& neiListx,
              const boost::python::list& neiListy,
              const boost::python::list& neiListz,
              const boost::python::list& halfCellMaskx,
              const boost::python::list& halfCellMasky,
              const boost::python::list& halfCellMaskz)
        : DomainDecompositionDuplicateTemplate<DomainDecomposition>(system, _nodeGrid, neiListx, neiListy, neiListz,
                            halfCellMaskx, halfCellMasky, halfCellMaskz)
      {
      }

      static void registerPython();
    };

    class DomainDecompositionReference : public DomainDecomposition
    {
      public:

      DomainDecompositionReference(shared_ptr< System > system,
              const Int3D& _nodeGrid,
              const boost::python::list& neiListx,
              const boost::python::list& neiListy,
              const boost::python::list& neiListz,
              const boost::python::list& halfCellMaskx,
              const boost::python::list& halfCellMasky,
              const boost::python::list& halfCellMaskz)
        : DomainDecomposition(system, _nodeGrid, neiListx, neiListy, neiListz,
                            halfCellMaskx, halfCellMasky, halfCellMaskz)
        {
        }
      // follow the routines, which have to be redefined
      // they perform communication, and thus they will be performed using dd_comm
      //--------------------------------------------------------
      public:
      virtual void updateGhosts();
      virtual void updateGhostsV();
      virtual void collectGhostForces();
      protected:
      virtual void decomposeRealParticles();
      virtual void exchangeGhosts();
      //--------------------------------------------------------
      public:
      static void registerPython();
    };

  } // namespace
} // namespace

#endif
