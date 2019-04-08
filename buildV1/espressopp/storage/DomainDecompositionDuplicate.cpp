#include <storage/Storage.hpp>
#include <storage/DomainDecomposition.hpp>
#include <storage/DomainDecompositionDuplicate.hpp>

//#define DEBUG_OUTPUT
std::vector<int> only_cells = {23599, 23598, 23597};

namespace espressopp{
  namespace storage{

    template<typename DD>
    DomainDecompositionDuplicateTemplate<DD>::DomainDecompositionDuplicateTemplate(
              shared_ptr< System > system,
              const Int3D& _nodeGrid,
              const boost::python::list& neiListx,
              const boost::python::list& neiListy,
              const boost::python::list& neiListz,
              const boost::python::list& halfCellMaskx,
              const boost::python::list& halfCellMasky,
              const boost::python::list& halfCellMaskz)
      :DD(system, _nodeGrid, neiListx, neiListy, neiListz, halfCellMaskx, halfCellMasky, halfCellMaskz)
    {
      dd_comm = make_shared<DD>(system, 1, _nodeGrid, neiListx, neiListy, neiListz);

      for(int x_fine = 0; x_fine < this->cellGrid.getFrameGridSize(0); x_fine++)
      {
        int x_coarse = x_fine / this->halfCellInt;
        for(int y_fine = 0; y_fine < this->cellGrid.getFrameGridSize(1); y_fine++)
        {
          int y_coarse = y_fine / this->halfCellInt;
          for(int z_fine = 0; z_fine < this->cellGrid.getFrameGridSize(2); z_fine++)
          {
            int z_coarse = z_fine / this->halfCellInt;
            longint index_fine = this->cellGrid.mapPositionToIndex(x_fine, y_fine, z_fine);
            longint index_coarse = dd_comm->getCellGrid().mapPositionToIndex(x_coarse, y_coarse, z_coarse);
//            std::cout << "ABC " << this->cellGrid.getFrameGridSize(0) << ", " << this->cellGrid.getFrameGridSize(1) << ", "
//                      << this->cellGrid.getFrameGridSize(2) << ", " << index_fine << " (" << x_fine << ", " << y_fine << ", " << z_fine << "), "
//                      << index_coarse << " (" << x_coarse << ", " << y_coarse << ", " << z_coarse << ")" << std::endl;
            fine_to_coarse.insert({&this->cells[index_fine], &(this->dd_comm->cells[index_coarse])});
          }
        }
      }


      // connect signals
      //this->afterRecvParticles.connect(dd_comm->afterRecvParticles);
      dd_comm->afterRecvParticles.connect(this->afterRecvParticles);
      dd_comm->beforeSendParticles.connect(this->beforeSendParticles);

      // zero the timer
      this->timerOverheads.reset();
      this->timeOverheadUpdateGhosts = 0.0;
      this->timeOverheadCollectForces = 0.0;
      this->timeOverheadDecomposeReal = 0.0;
      this->timeOverheadExchangeGhosts = 0.0;

    }

    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::check_layer_consistency(std::string msg)
    {
        if(this->getNRealParticles() != dd_comm->getNRealParticles() || this->getNGhostParticles() != dd_comm->getNGhostParticles()){
            std::cout << "LAYERS INCONSISTENT " << msg << ": COMPUTE: " << this->getNLocalParticles() << ", " << this->getNRealParticles() << ", "
                                                                        << this->getNRealBoundaryParticles() << ", " << this->getNGhostParticles()
                         << ", COMM: " << dd_comm->getNLocalParticles() << ", " << dd_comm->getNRealParticles() << ", " << dd_comm->getNRealBoundaryParticles()
                                                                        << ", " << dd_comm->getNGhostParticles() << std::endl;
        }

    }


    //************************************************
    // transfer REAL particles from COMPUTE to COMM layer before calling decompose
    // ***********************************************
    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::transRealParticlesBeforeDecompose()
    {
        dd_comm->removeAllParticles();

        for (Cell* cell : this->realCells)
        {
            for (Particle p : cell->particles)
            {
                Cell* comm_cell = dd_comm->mapPositionToCellClipped(p.position());

//                if(comm_cell != dd_comm->mapPositionToCell(p.position())){
//                    Cell* comm_cell_2 = dd_comm->mapPositionToCell(p.position());
//                    throw(std::runtime_error("transRealParticlesBeforeDecompose: cell discrepancy " +
//                                            std::to_string(comm_cell->id) + " : " + std::to_string(comm_cell_2->id) +
//                                            ", particle position " + std::to_string(p.position().at(0)) + " " +
//                                            std::to_string(p.position().at(1)) + " " + std::to_string(p.position().at(2))));
//                }
                if(std::find(dd_comm->realCells.begin(), dd_comm->realCells.end(), comm_cell) == dd_comm->realCells.end())
                    throw(std::runtime_error("transRealParticlesBeforeDecompose: not a reall cell"));

                comm_cell->particles.push_back(p);
            }
        }

        for (Cell* cell : dd_comm->realCells)
        {
            dd_comm->updateLocalParticles(cell->particles);
        }

        this->removeAllParticlesNoSignal();
    }


    //************************************************
    // transfer REAL particles from COMM to COMPUTE layer after calling decompose
    // ***********************************************
    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::transRealParticlesAfterDecompose()
    {
        // real particles have been decomposed, so we should fix to which cell do they belong in the comm DD.
        this->cell_in_dd_comm.clear();
        for (Cell* cell : dd_comm->realCells)
        {
            for (Particle p : cell->particles)
            {
               this->cell_in_dd_comm[p.id()] = cell;
            }
        }

        this->removeAllParticlesNoSignal();

        this->order_of_real_cells.clear();

        for (Cell* cell : dd_comm->realCells)
        {
            bool isBoundaryCell = (std::find(dd_comm->realBoundaryCells.begin(), dd_comm->realBoundaryCells.end(), cell)
                    != dd_comm->realBoundaryCells.end());
            for (Particle p : cell->particles)
            {
                Cell* comm_cell = dd_comm->mapPositionToCellChecked(p.position());

                if ((comm_cell == 0) || (comm_cell != cell))
                    throw(std::runtime_error("transRealParticlesAfterDecompose, comm, Projected cell not the same as expected " + std::to_string((long int)comm_cell)));
                if(std::find(dd_comm->realCells.begin(), dd_comm->realCells.end(), comm_cell) == dd_comm->realCells.end())
                    throw(std::runtime_error("transRealParticlesAfterDecompose: not a reall cell"));

                Cell* compute_cell = this->mapPositionToCellChecked(p.position());

                if ((compute_cell == 0) || (fine_to_coarse[compute_cell] != comm_cell))
                    throw(std::runtime_error("transRealParticlesAfterDecompose, compute, Projected cell not the same as expected: "
#ifdef CELL_EXTRA_DATA
                                             + std::to_string((unsigned long int) comm_cell->id) + ", "
                                             + std::to_string((unsigned long int) fine_to_coarse[compute_cell]->id)));
#else
                                             ));
#endif

                if(std::find(this->realCells.begin(), this->realCells.end(), compute_cell) == this->realCells.end())
                    throw(std::runtime_error("transRealParticlesAfterDecompose: not a reall cell"));

                compute_cell->particles.push_back(p);
                if(isBoundaryCell)
                    order_of_real_cells.push_back(compute_cell);
            }
        }

        // after the particles were moved to the compute layer (compute DD), we should updateLocalParticles
        // and also fix to which cell do particles belong in the compute DD
        this->cell_in_dd_compute.clear();
        for (Cell* cell : this->realCells)
        {
            this->updateLocalParticles(cell->particles);
            for (Particle p : cell->particles)
            {
               this->cell_in_dd_compute[p.id()] = cell;
            }
        }

        //cannot remove here, since I would also remove the ghost particles!!!
        // !!!!!
        //dd_comm->removeAllParticles();
    }

    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::transGhostParticlesBeforeExchangeGhosts()
    {
        // actually no need to do anything, assuming, that ghost particles have been already
        // cleared in transRealParticlesBeforeDecompose()
    }

    //************************************************
    // transfer GHOST particles from COMM to COMPUTE layer after calling exchange ghosts
    // ***********************************************
    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::transGhostParticlesAfterExchangeGhosts()
    {
        // this->removeAllParticles removed this-> ghost cells when calling transRealParticlesAfterDecompose()
        //
        this->order_of_ghost_cells.clear();

        // ghost particles have been exchanged, so we should fix to which cell do they belong in the comm DD.
        //this->cell_in_dd_comm.clear(); .. this has already been cleared for reals
        for (Cell* cell : dd_comm->ghostCells)
        {
            for (Particle p : cell->particles)
            {
               if(this->cell_in_dd_comm.find(p.id()) != this->cell_in_dd_comm.end())
                   throw(std::runtime_error("transGhostParticlesAfterExchangeGhosts, cell_id_dd_comm already contains index " + std::to_string(p.id())));
               this->cell_in_dd_comm[p.id()] = cell;
            }
        }

        for (Cell* cell : dd_comm->ghostCells)
        {
            //std::cout << "GOING through ghost cells, " << cell->id << std::endl;
            for (Particle p : cell->particles)
            {
                Cell* comm_cell = dd_comm->mapPositionToCell(p.position());
                //TODO shouldn't it be this? it is a bit messy...
                //Cell* comm_cell = dd_comm->mapPositionToCellWithGhosts(p.position());
                //std::cout << "projecting PARTICLE back to calc, " << p.id() << ", " << cell->id << ", " << comm_cell->id << std::endl;
                // particles in dd_comm have been resorted, so they should be in the correct cell (given the position)
                if ((comm_cell == 0) || (comm_cell != cell))
                    throw(std::runtime_error(
                    //std::cout << ((
                                "transGhostParticlesAfterExchangeGhosts, comm, Projected cell not the same as expected, "
                           + std::to_string(p.id()) + ": " + std::string(p.position())
#ifdef CELL_EXTRA_DATA
          //<< ", CELL " << cell->id << "[" << cell->grid_pos << "], (" << cell->myLeft << "<->" << cell->myRight << ")"
                                             + ", comm_cell: " + std::to_string(comm_cell->id) + ": " + std::string(comm_cell->myLeft) + "<->" + std::string(comm_cell->myRight)
                                             + ", cell: " + std::to_string(cell->id) + ": " + std::string(cell->myLeft) + "<->" + std::string(cell->myRight)));
#else
                                             ));
#endif
                if(std::find(dd_comm->ghostCells.begin(), dd_comm->ghostCells.end(), comm_cell) == dd_comm->ghostCells.end()){
                    //std::cout <<  ((
                    throw(std::runtime_error(
                            "transGhostParticlesAfterExchangeGhosts: not a ghost cell in COMM layer"));
                }

                // carefull - ghost cells in compute layer have to be treated differently
                Cell* compute_cell = this->mapPositionToCellWithGhosts(p.position());
                if ((compute_cell == 0) || (fine_to_coarse[compute_cell] != comm_cell))
                {
                    throw(std::runtime_error("transGhostParticlesAfterExchangeGhosts, compute, Projected cell not the same as expected"));
                }
                if(std::find(this->ghostCells.begin(), this->ghostCells.end(), compute_cell) == this->ghostCells.end()) {
                    //std::cout << ((
                    throw(std::runtime_error(
                    "transGhostParticlesAfterExchangeGhosts: not a ghost cell in COMPUTE layer"));
                }

                compute_cell->particles.push_back(p);
                order_of_ghost_cells.push_back(compute_cell);

            }
        }

        for (Cell* cell : this->ghostCells)
        {
            this->updateLocalParticles(cell->particles);
            for (Particle p : cell->particles)
            {
               if(this->cell_in_dd_compute.find(p.id()) != this->cell_in_dd_compute.end())
                   throw(std::runtime_error("transGhostParticlesAfterExchangeGhosts, cell_id_dd_compute already contains index " + std::to_string(p.id())));
               this->cell_in_dd_compute[p.id()] = cell;
            }
        }
    }

    //****************************************

    struct IteratorMap : public std::map<Cell*, ParticleList::iterator> {
        public:
            IteratorMap(CellList list) {
                num_particles = 0;
                for (Cell* cell : list){
                    (*this)[cell] = cell->particles.begin();
                    num_particles += cell->particles.size();
                }
            }

            void checkExhausted() {
                for (const auto &x : (*this)){
                    if(x.second != x.first->particles.end())
                        throw(std::runtime_error("transParticlesBeforeUpdateOrCollect, some of real cell particle lists not exhausted"));
                }
            }

            int num_particles;
    };


    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>:: copyParticlePropertiesToComm(CellList computeSource, CellList commTarget, std::vector<Cell*> cellOrder, bool positions)
    {
       // init iterators
       IteratorMap compute_iterators(computeSource);
       IteratorMap comm_iterators(commTarget);

       if(compute_iterators.num_particles != comm_iterators.num_particles)
           throw(std::runtime_error("transParticlesBeforeUpdateOrCollect, different number of COMM and COMPUTE particles"));

       // insert particles to comm layer in the correct order
       for(Cell* compute_cell_ptr: cellOrder)
       {
           Particle* p_compute = &*(compute_iterators[compute_cell_ptr]++);
           Cell* comm_cell_ptr = cell_in_dd_comm[p_compute->id()];
           Particle* p_comm = &*(comm_iterators[comm_cell_ptr]++);

           if(p_compute->id() != p_comm->id())
               throw(std::runtime_error("transParticlesBeforeUpdateOrCollect, particles comm and compute do not match"));

           //copy 
           if(positions)
               p_comm->position() = p_compute->position();
           else
               p_comm->force() = p_compute->force();
       }

       // make sure that all particles from compute layer have been used
       compute_iterators.checkExhausted();
       comm_iterators.checkExhausted();
    }

    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::copyParticlePropertiesToCompute(CellList commSource, bool positions)
    {
        for(Cell* cell: commSource)
        {
            for(Particle p_comm : cell->particles)
            {
                Particle* p_compute = this->lookupLocalParticle(p_comm.id());
                if(p_compute == 0)
                    throw(std::runtime_error("transParticlesAfterUpdateOrCollect, such particle does not exist in compute layer"));
                if(positions)
                    p_compute->position() = p_comm.position();
                else
                    p_compute->force() = p_comm.force();
            }
        }
    }

    // transfer REAL particles from COMPUTE to COMM layer before calling update ghosts
    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::transParticlesBeforeUpdate()
    {
       // befor the update ghosts, it is enough to copy positions of real particles, since the old
       // values of positions in the ghost layer are rewritten
       copyParticlePropertiesToComm(this->realBoundaryCells, dd_comm->realBoundaryCells, this->order_of_real_cells, true);
    }

    // transfer ALL particles from COMPUTE to COMM layer before calling collect ghost forces
    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::transParticlesBeforeCollect()
    {
       // before the collect ghost forces, both real and ghost layer have to be updated. That is because
       // forces are additive and as received from the ghosts, they are ADDED to the values in real particles
       copyParticlePropertiesToComm(this->realBoundaryCells, dd_comm->realBoundaryCells, this->order_of_real_cells, false);
       copyParticlePropertiesToComm(this->ghostCells, dd_comm->ghostCells, this->order_of_ghost_cells, false);
    }

    // transfer GHOST particles from COMM to COMPUTE layer after calling update ghosts or collect ghost forces
    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::transParticlesAfterUpdate()
    {
        copyParticlePropertiesToCompute(dd_comm->ghostCells, true);
    }

    // transfer REAL particles from COMM to COMPUTE layer after calling update ghosts or collect ghost forces
    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::transParticlesAfterCollect()
    {
        copyParticlePropertiesToCompute(dd_comm->realBoundaryCells, false);
    }


    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::updateGhosts()
    {
      real time;
#ifdef DEBUG_OUTPUT
      std::cout << "DDD updateGhosts" << std::endl;
      this->printAllParticles(   "DD_COMPUTATION  1 update gho", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 1 update gho", only_cells);
      this->getSystem()->comm->barrier();
#endif
      time = this->timerOverheads.getElapsedTime();
      transParticlesBeforeUpdate();
      this->timeOverheadUpdateGhosts += this->timerOverheads.getElapsedTime() - time;
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  2 update gho", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 2 update gho", only_cells);
      this->getSystem()->comm->barrier();
#endif

      // do the actuall operation
      dd_comm->updateGhosts();

#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  3 update gho", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 3 update gho", only_cells);
      this->getSystem()->comm->barrier();
#endif
      time = this->timerOverheads.getElapsedTime();
      transParticlesAfterUpdate();
      this->timeOverheadUpdateGhosts += this->timerOverheads.getElapsedTime() - time;
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  4 update gho", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 4 update gho", only_cells);
      std::cout << "DDD END updateGhosts" << std::endl;
#endif
    }

    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::updateGhostsV()
    {
     real time;
#ifdef DEBUG_OUTPUT
      std::cout << "DDD updateGhostsV" << std::endl;
#endif
      time = this->timerOverheads.getElapsedTime();
      transParticlesBeforeUpdate();
      this->timeOverheadUpdateGhosts += this->timerOverheads.getElapsedTime() - time;

      // do the actuall operation
      dd_comm->updateGhostsV();

      time = this->timerOverheads.getElapsedTime();
      transParticlesAfterUpdate();
      this->timeOverheadUpdateGhosts += this->timerOverheads.getElapsedTime() - time;
#ifdef DEBUG_OUTPUT
      std::cout << "DDD END updateGhostsV" << std::endl;
#endif
    }

    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::collectGhostForces()
    {
        real time;

#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      std::cout << "DDD collectGhostForces" << std::endl;
      this->printAllParticles(   "DD_COMPUTATION  1 col gh for", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 1 col gh for", only_cells);
      this->getSystem()->comm->barrier();
#endif
      time = this->timerOverheads.getElapsedTime();
      transParticlesBeforeCollect();
      this->timeOverheadCollectForces += this->timerOverheads.getElapsedTime() - time;
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  2 col gh for", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 2 col gh for", only_cells);
      this->getSystem()->comm->barrier();
#endif

      // do the actuall operation
      dd_comm->collectGhostForces();

#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  3 col gh for", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 3 col gh for", only_cells);
      this->getSystem()->comm->barrier();
#endif
      time = this->timerOverheads.getElapsedTime();
      transParticlesAfterCollect();
      this->timeOverheadCollectForces += this->timerOverheads.getElapsedTime() - time;
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  4 col gh for", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 4 col gh for", only_cells);
      std::cout << "DDD END collectGhostForces" << std::endl;
#endif
    }

    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::decomposeRealParticles()
    {
      real time;
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      std::cout << "DDD decomposeRealParticles" << std::endl;
      this->printAllParticles(   "DD_COMPUTATION  1 dec real p", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 1 dec real p", only_cells);
      this->getSystem()->comm->barrier();
#endif
      time = this->timerOverheads.getElapsedTime();
      transRealParticlesBeforeDecompose();
      this->timeOverheadDecomposeReal += this->timerOverheads.getElapsedTime() - time;
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  2 dec real p", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 2 dec real p", only_cells);
      this->getSystem()->comm->barrier();
#endif

      // do the actuall operation
      dd_comm->decomposeRealParticles();

#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  3 dec real p", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 3 dec real p", only_cells);
      this->getSystem()->comm->barrier();
#endif
      time = this->timerOverheads.getElapsedTime();
      transRealParticlesAfterDecompose();
      this->timeOverheadDecomposeReal += this->timerOverheads.getElapsedTime() - time;
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  4 dec real p", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 4 dec real p", only_cells);
      std::cout << "DDD END decomposeRealParticles" << std::endl;
#endif
    }

    template<typename DD>
    void DomainDecompositionDuplicateTemplate<DD>::exchangeGhosts()
    {
      real time;
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      std::cout << "DDD exchangeGhosts" << std::endl;
      this->printAllParticles(   "DD_COMPUTATION  1 exch ghost", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 1 exch ghost", only_cells);
      this->getSystem()->comm->barrier();
#endif
      //transRealParticlesBeforeDecompose();
      //transGhostParticlesBeforeExchangeGhosts();
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  2 exch ghost", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 2 exch ghost", only_cells);
      this->getSystem()->comm->barrier();
#endif

      // do the actuall operation
      dd_comm->exchangeGhosts();

#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  3 exch ghost", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 3 exch ghost", only_cells);
      this->getSystem()->comm->barrier();
#endif
      time = this->timerOverheads.getElapsedTime();
      //transRealParticlesAfterDecompose();
      transGhostParticlesAfterExchangeGhosts();
      this->timeOverheadExchangeGhosts += this->timerOverheads.getElapsedTime() - time;
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  4 exch ghost", only_cells);
      dd_comm->printAllParticles("DD_COMUNICATION 4 exch ghost", only_cells);
      std::cout << "DDD END exchangeGhosts" << std::endl;
#endif
      this->check_layer_consistency("after exchange ghosts");
    }


    ///***********************************************************************


    void DomainDecompositionReference::updateGhosts()
    {
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      std::cout << "DDD updateGhosts" << std::endl;
      this->printAllParticles(   "DD_COMPUTATION  1 update gho", only_cells);
      this->getSystem()->comm->barrier();
#endif

      DomainDecomposition::updateGhosts();

#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  4 update gho", only_cells);
      std::cout << "DDD END updateGhosts" << std::endl;
      this->getSystem()->comm->barrier();
#endif
    }

    void DomainDecompositionReference::updateGhostsV()
    {
      std::cout << "DDD updateGhostsV" << std::endl;
      DomainDecomposition::updateGhostsV();
      std::cout << "DDD END updateGhostsV" << std::endl;
    }

    void DomainDecompositionReference::collectGhostForces()
    {
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      std::cout << "DDD collectGhostForces" << std::endl;
      this->printAllParticles(   "DD_COMPUTATION  1 col gh for", only_cells);
      this->getSystem()->comm->barrier();
#endif

      DomainDecomposition::collectGhostForces();

#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  4 col gh for", only_cells);
      std::cout << "DDD END collectGhostForces" << std::endl;
      this->getSystem()->comm->barrier();
#endif
    }

    void DomainDecompositionReference::decomposeRealParticles()
    {
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      std::cout << "DDD decomposeRealParticles" << std::endl;
      this->printAllParticles(   "DD_COMPUTATION  1 dec real p", only_cells);
      this->getSystem()->comm->barrier();
#endif

      DomainDecomposition::decomposeRealParticles();

#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  4 dec real p", only_cells);
      std::cout << "DDD END decomposeRealParticles" << std::endl;
      this->getSystem()->comm->barrier();
#endif
    }

   void DomainDecompositionReference::exchangeGhosts()
   {
#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      std::cout << "DDD exchangeGhosts" << std::endl;
      this->printAllParticles(   "DD_COMPUTATION  1 exch ghost", only_cells);
      this->getSystem()->comm->barrier();
#endif

      this->DomainDecomposition::exchangeGhosts();

#ifdef DEBUG_OUTPUT
      this->getSystem()->comm->barrier();
      this->printAllParticles(   "DD_COMPUTATION  4 exch ghost", only_cells);
      std::cout << "DDD END exchangeGhosts" << std::endl;
      this->getSystem()->comm->barrier();
#endif

   }

  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void DomainDecompositionDuplicate::registerPython() {
    using namespace espressopp::python;
    class_< DomainDecompositionDuplicate, bases< Storage >, boost::noncopyable >
    ("storage_DomainDecompositionDuplicate", init< shared_ptr< System >,
          //int,
          const Int3D&,
          //const Int3D&,
          const boost::python::list& ,const boost::python::list&, const boost::python::list&,
          const boost::python::list& ,const boost::python::list&, const boost::python::list& >())
    .def("mapPositionToNodeClipped", &DomainDecompositionDuplicate::mapPositionToNodeClipped)
    .def("getCellGrid", &DomainDecompositionDuplicate::getInt3DCellGrid)
    .def("getNodeGrid", &DomainDecompositionDuplicate::getInt3DNodeGrid)
    .def("cellAdjust", &DomainDecompositionDuplicate::cellAdjust)
    .def("printNodeParticles", &DomainDecomposition::printNodeParticles)
    ;
  }

  void DomainDecompositionReference::registerPython() {
    using namespace espressopp::python;
    class_< DomainDecompositionReference, bases< Storage >, boost::noncopyable >
    ("storage_DomainDecompositionReference", init< shared_ptr< System >,
          //int,
          const Int3D&,
          //const Int3D&,
          const boost::python::list& ,const boost::python::list&, const boost::python::list&,
          const boost::python::list& ,const boost::python::list&, const boost::python::list& >())
    .def("mapPositionToNodeClipped", &DomainDecompositionReference::mapPositionToNodeClipped)
    .def("getCellGrid", &DomainDecompositionReference::getInt3DCellGrid)
    .def("getNodeGrid", &DomainDecompositionReference::getInt3DNodeGrid)
    .def("cellAdjust", &DomainDecompositionReference::cellAdjust)
    .def("printNodeParticles", &DomainDecomposition::printNodeParticles)
    ;
  }

  } //namespace
} //namespace
