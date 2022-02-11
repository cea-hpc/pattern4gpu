#include "Pattern4GPUModule.h"

#define P4GPU_PROFILING // Pour activer le profiling
#include "P4GPUTimer.h"

#include "Pattern4GPU4Kokkos.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/* Implémentation d'origine de computeCqsAndVector()                         */
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
_computeCqsAndVector_Vori() {

  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans _computeCqsAndVector_Vori";

  constexpr Real k025 = 0.25;

  CellGroup active_cells = defaultMesh()->cellFamily()->findGroup("active_cells");

  P4GPU_DECLARE_TIMER(subDomain(), ComputeCqs); P4GPU_START_TIMER(ComputeCqs);
  UniqueArray<Real3> pos(8);
  ENUMERATE_CELL (cell_i,  allCells()) {
    for (Integer ii = 0; ii < 8; ++ii) {
      pos[ii] = m_node_coord_bis[cell_i->node(ii)];
    }

    m_cell_cqs[cell_i][0] = -k025*math::vecMul(pos[4]-pos[3], pos[1]-pos[3]);
    m_cell_cqs[cell_i][1] = -k025*math::vecMul(pos[0]-pos[2], pos[5]-pos[2]);
    m_cell_cqs[cell_i][2] = -k025*math::vecMul(pos[1]-pos[3], pos[6]-pos[3]);
    m_cell_cqs[cell_i][3] = -k025*math::vecMul(pos[7]-pos[2], pos[0]-pos[2]);
    m_cell_cqs[cell_i][4] = -k025*math::vecMul(pos[5]-pos[7], pos[0]-pos[7]);
    m_cell_cqs[cell_i][5] = -k025*math::vecMul(pos[1]-pos[6], pos[4]-pos[6]);
    m_cell_cqs[cell_i][6] = -k025*math::vecMul(pos[5]-pos[2], pos[7]-pos[2]);
    m_cell_cqs[cell_i][7] = -k025*math::vecMul(pos[6]-pos[3], pos[4]-pos[3]);
  }
  P4GPU_STOP_TIMER(ComputeCqs);

  P4GPU_DECLARE_TIMER(subDomain(), NodeVectorUpdate); P4GPU_START_TIMER(NodeVectorUpdate);
  ENUMERATE_NODE (node_i, allNodes()) {
    m_node_vector[node_i].assign(0., 0., 0.);
  }

  // Calcul du gradient de pression
  ENUMERATE_CELL (cell_i, active_cells) {
    ENUMERATE_NODE (node_i, cell_i->nodes()) {
      m_node_vector[node_i] += (m_cell_arr1[cell_i] + m_cell_arr2[cell_i]) *
        m_cell_cqs[cell_i][node_i.index()];
    }
  }

  m_node_vector.synchronize();
  P4GPU_STOP_TIMER(NodeVectorUpdate);

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Implémentation CPU Arcane multi-thread                                    */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_computeCqsAndVector_Vmt() {

  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans _computeCqsAndVector_Vmt";

  ParallelLoopOptions options;

  // On calcule les CQs sur les mailles
  options.setPartitioner(ParallelLoopOptions::Partitioner::Auto);
  arcaneParallelForeach(allCells(), options, [&](CellVectorView cells) {
    constexpr Real k025 = 0.25;
    ENUMERATE_CELL (cell_i, cells) {
      // Recopie les coordonnées locales (pour le cache)
      Real3 pos[8];
      for (Integer ii = 0; ii < 8; ++ii) {
        pos[ii] = m_node_coord_bis[cell_i->node(ii)];
      }

      m_cell_cqs[cell_i][0] = -k025*math::vecMul(pos[4]-pos[3], pos[1]-pos[3]);
      m_cell_cqs[cell_i][1] = -k025*math::vecMul(pos[0]-pos[2], pos[5]-pos[2]);
      m_cell_cqs[cell_i][2] = -k025*math::vecMul(pos[1]-pos[3], pos[6]-pos[3]);
      m_cell_cqs[cell_i][3] = -k025*math::vecMul(pos[7]-pos[2], pos[0]-pos[2]);
      m_cell_cqs[cell_i][4] = -k025*math::vecMul(pos[5]-pos[7], pos[0]-pos[7]);
      m_cell_cqs[cell_i][5] = -k025*math::vecMul(pos[1]-pos[6], pos[4]-pos[6]);
      m_cell_cqs[cell_i][6] = -k025*math::vecMul(pos[5]-pos[2], pos[7]-pos[2]);
      m_cell_cqs[cell_i][7] = -k025*math::vecMul(pos[6]-pos[3], pos[4]-pos[3]);
    }
  });

  auto node_index_in_cells = m_acc_env->nodeIndexInCells();
  const Integer max_node_cell = m_acc_env->maxNodeCell();

  // Puis, on applique les CQs sur les noeuds
  arcaneParallelForeach(allNodes(), options, [&](NodeVectorView nodes) {
    ENUMERATE_NODE (node_i, nodes) {
      Int32 first_pos = node_i.localId() * max_node_cell;
      Real3 node_vec = Real3::zero();
      ENUMERATE_CELL(cell_i, node_i->cells()) {
        if (m_is_active_cell[cell_i]) { // la maille ne contribue que si elle est active
          Int16 node_index = node_index_in_cells[first_pos + cell_i.index()];
          node_vec += (m_cell_arr1[cell_i]+m_cell_arr2[cell_i])
            * m_cell_cqs[cell_i][node_index];
        }
      }
      m_node_vector[node_i] = node_vec;
    }
  });

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Implémentation CPU version 2                                              */
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
_computeCqsAndVector_Vmt_v2() {

  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans _computeCqsAndVector_Vmt_v2";

  ParallelLoopOptions options;
  options.setPartitioner(ParallelLoopOptions::Partitioner::Auto);

  arcaneParallelForeach(allNodes(), options, [&](NodeVectorView nodes) {
  ENUMERATE_NODE (node_i, nodes) {
    m_node_vector[node_i].assign(0., 0., 0.);
  }
  });
//  m_node_vector.fill(Real3::zero());

  constexpr Real k025 = 0.25;

  CellGroup active_cells = defaultMesh()->cellFamily()->findGroup("active_cells");

  const Integer ipos[8][4] = {
    {4, 3, 1, 3},
    {0, 2, 5, 2},
    {1, 3, 6, 3},
    {7, 2, 0, 2},
    {5, 7, 0, 7},
    {1, 6, 4, 6},
    {5, 2, 7, 2},
    {6, 3, 4, 3}
  };
  for(Integer idx_node=0 ; idx_node<8 ; ++idx_node) {

    arcaneParallelForeach(active_cells, options, [&](CellVectorView cells) {
    const Integer ip0 = ipos[idx_node][0];
    const Integer ip1 = ipos[idx_node][1];
    const Integer ip2 = ipos[idx_node][2];
    const Integer ip3 = ipos[idx_node][3];

    ENUMERATE_CELL (cell_i, cells) {

      const Real3& nd0 = m_node_coord_bis[cell_i->node(ip0)];
      const Real3& nd1 = m_node_coord_bis[cell_i->node(ip1)];
      const Real3& nd2 = m_node_coord_bis[cell_i->node(ip2)];
      const Real3& nd3 = m_node_coord_bis[cell_i->node(ip3)];

      Real3 cell_cqs = -k025*math::vecMul(nd0-nd1, nd2-nd3);

      m_node_vector[cell_i->node(idx_node)] += 
        (m_cell_arr1[cell_i] + m_cell_arr2[cell_i]) * cell_cqs;
    }
    });
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Implémentation API GPU Arcane version 1                                   */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline void computeCQs(Real3 pos[8], Span<Real3> out_cqs) {
  constexpr Real k025 = 0.25;
  Real3 p0 = pos[0];
  Real3 p1 = pos[1];
  Real3 p2 = pos[2];
  Real3 p3 = pos[3];
  Real3 p4 = pos[4];
  Real3 p5 = pos[5];
  Real3 p6 = pos[6];
  Real3 p7 = pos[7];

  out_cqs[0] = -k025*math::vecMul(p4-p3, p1-p3);
  out_cqs[1] = -k025*math::vecMul(p0-p2, p5-p2);
  out_cqs[2] = -k025*math::vecMul(p1-p3, p6-p3);
  out_cqs[3] = -k025*math::vecMul(p7-p2, p0-p2);
  out_cqs[4] = -k025*math::vecMul(p5-p7, p0-p7);
  out_cqs[5] = -k025*math::vecMul(p1-p6, p4-p6);
  out_cqs[6] = -k025*math::vecMul(p5-p2, p7-p2);
  out_cqs[7] = -k025*math::vecMul(p6-p3, p4-p3);
}

void Pattern4GPUModule::
_computeCqsAndVector_Varcgpu_v1() {

  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans _computeCqsAndVector_Varcgpu_v1";

  P4GPU_DECLARE_TIMER(subDomain(), ComputeCqs); P4GPU_START_TIMER(ComputeCqs);

  bool is_sync_cell_cqs=!(options()->getCcavCqsSyncVersion() == CCAV_CS_nosync);
  CellGroup cell_group=(is_sync_cell_cqs ? ownCells() : allCells());
  {
    auto queue = m_acc_env->newQueue();
    auto command = makeCommand(queue);

    auto in_node_coord_bis = ax::viewIn(command,m_node_coord_bis);
    auto out_cell_cqs = ax::viewInOut(command,m_cell_cqs);

    auto cnc = m_acc_env->connectivityView().cellNode();

    command << RUNCOMMAND_ENUMERATE(Cell, cid, cell_group){
      // Recopie les coordonnées locales (pour le cache)
      Real3 pos[8];
      Int32 index=0;
      for( NodeLocalId nid : cnc.nodes(cid) ){
        pos[index]=in_node_coord_bis[nid];
        ++index;
      }
      // Calcule les résultantes aux sommets
      computeCQs(pos, out_cell_cqs[cid]);
    };
  }
  if (is_sync_cell_cqs) {
    PROF_ACC_BEGIN("syncCellCQs");
    if (options()->getCcavCqsSyncVersion() == CCAV_CS_bulksync_std) {
      m_cell_cqs.synchronize();
    } else {
      ARCANE_ASSERT(options()->getCcavCqsSyncVersion() == CCAV_CS_bulksync_sync,
          ("CCAV_CS_bulksync_sync obligatoire"));
      m_acc_env->vsyncMng()->globalSynchronize(m_cell_cqs);
    }
    PROF_ACC_END;
  }
  P4GPU_STOP_TIMER(ComputeCqs);

  P4GPU_DECLARE_TIMER(subDomain(), NodeVectorUpdate); P4GPU_START_TIMER(NodeVectorUpdate);

  auto async_node_vector_update = [&](NodeGroup node_group, bool high_priority=false) -> Ref<RunQueue> {
    // On inverse boucle Cell <-> Node car la boucle originelle sur les mailles n'est parallélisable
    // Du coup, on boucle sur les Node
    // On pourrait construire un groupe de noeuds des mailles active_cells et boucler sur ce groupe
    // Mais ici, on a pré-construit un tableau global au mailles qui indique si une maille fait partie
    // du groupe active_cells ou pas (m_is_active_cell)
    // Rem : en décomp. de dom., pour la plupart des sous-dom. on aura : active_cells = allCells
    // Ainsi, en bouclant sur tous les noeuds allNodes(), pour un noeud donné :
    //   1- on initialise le vecteur à 0
    //   2- on calcule les constributions des mailles connectées au noeud uniquement si elles sont actives
    ax::RunQueueBuildInfo bi;
    if (high_priority) {
      // 0 = priorité par défaut
      // Plus la valeur de priorité est faible, plus la queue sera prioritaire
      bi.setPriority(-10); 
    }
    auto ref_queue = makeQueueRef(m_acc_env->runner(), bi);
    ref_queue->setAsync(true);
    auto command = makeCommand(ref_queue.get());

    auto in_cell_arr1 = ax::viewIn(command, m_cell_arr1);
    auto in_cell_arr2 = ax::viewIn(command, m_cell_arr2);
    auto in_cell_cqs  = ax::viewIn(command, m_cell_cqs);
    auto in_is_active_cell = ax::viewIn(command, m_is_active_cell);

    auto out_node_vector = ax::viewOut(command, m_node_vector);

    auto node_index_in_cells = m_acc_env->nodeIndexInCells();
    const Integer max_node_cell = m_acc_env->maxNodeCell();

    auto nc_cty = m_acc_env->connectivityView().nodeCell();

    command << RUNCOMMAND_ENUMERATE(Node,nid,node_group) {
      Int32 first_pos = nid.localId() * max_node_cell;
      Integer index = 0;
      Real3 node_vec = Real3::zero();
      for( CellLocalId cid : nc_cty.cells(nid) ){
        if (in_is_active_cell[cid]) { // la maille ne contribue que si elle est active
          Int16 node_index = node_index_in_cells[first_pos + index];
          node_vec += (in_cell_arr1[cid]+in_cell_arr2[cid]) 
            * in_cell_cqs[cid][node_index];
        }
        ++index;
      }
      out_node_vector[nid] = node_vec;
    }; // non bloquant

    return ref_queue;
  };

  if (options()->getCcavVectorSyncVersion() == CCAV_VS_bulksync_std ||
      options()->getCcavVectorSyncVersion() == CCAV_VS_bulksync_queue) 
  {
    // Calcul sur tous les noeuds "own"
    Ref<RunQueue> ref_queue = async_node_vector_update(ownNodes());
    ref_queue->barrier();

    // Puis comms
    PROF_ACC_BEGIN("syncNodeVector");
    if (options()->getCcavVectorSyncVersion() == CCAV_VS_bulksync_std) {
      m_node_vector.synchronize();
    } else {
      ARCANE_ASSERT(options()->getCcavVectorSyncVersion()==CCAV_VS_bulksync_queue,
          ("Ici, option differente de bulksync_queue"));
      auto vsync = m_acc_env->vsyncMng();
      //vsync->globalSynchronize(m_node_vector);
      vsync->globalSynchronizeQueue(ref_queue, m_node_vector);
      //vsync->globalSynchronizeDevThr(m_node_vector);
      //vsync->globalSynchronizeDevQueues(m_node_vector);
    }
    PROF_ACC_END;

  } 
  else if (options()->getCcavVectorSyncVersion() == CCAV_VS_overlap_evqueue ||
      options()->getCcavVectorSyncVersion() == CCAV_VS_overlap_evqueue_d) 
  {
    auto vsync = m_acc_env->vsyncMng();
    SyncItems<Node>* sync_nodes = vsync->getSyncItems<Node>();

    // On amorce sur le DEVICE le calcul sur les noeuds "own" sur le bord du
    // sous-domaine sur la queue ref_queue_bnd (_bnd = boundary)
    Ref<RunQueue> ref_queue_bnd = async_node_vector_update(sync_nodes->sharedItems(), /*high_priority=*/true);
    // ici, le calcul n'est pas terminé sur le DEVICE

    // On amorce sur le DEVICE le calcul des noeuds intérieurs dont ne dépendent
    // pas les comms sur la queue ref_queue_inr (_inr = inner)
    Ref<RunQueue> ref_queue_inr = async_node_vector_update(sync_nodes->privateItems());

    // Sur la même queue de bord ref_queue_bnd, on amorce le packing des données
    // puis les comms MPI sur CPU, puis unpacking des données et on synchronise 
    // la queue ref_queue_bnd
//    vsync->globalSynchronizeQueue(ref_queue_bnd, m_node_vector);
    if (options()->getCcavVectorSyncVersion() == CCAV_VS_overlap_evqueue)
      vsync->globalSynchronizeQueueEvent(ref_queue_bnd, m_node_vector);
    else if (options()->getCcavVectorSyncVersion() == CCAV_VS_overlap_evqueue_d)
      vsync->globalSynchronizeQueueEventD(ref_queue_bnd, m_node_vector);
    // ici, après cet appel, ref_queue_bnd est synchronisée

    // On attend la terminaison des calculs intérieurs
    ref_queue_inr->barrier();
  } 
  else if (options()->getCcavVectorSyncVersion() == CCAV_VS_overlap_iqueue) 
  {
    auto vsync = m_acc_env->vsyncMng();
    SyncItems<Node>* sync_nodes = vsync->getSyncItems<Node>();

    // On amorce sur le DEVICE le calcul sur les noeuds "own" sur le bord du
    // sous-domaine sur la queue ref_queue_bnd (_bnd = boundary)
    Ref<RunQueue> ref_queue_bnd = async_node_vector_update(sync_nodes->sharedItems(), /*high_priority=*/true);
    // ici, le calcul n'est pas terminé sur le DEVICE
    
    // Sur la même queue de bord ref_queue_bnd, on amorce le packing des données
    auto ref_sync_req = vsync->iGlobalSynchronizeQueue(ref_queue_bnd, m_node_vector);

    // On amorce sur le DEVICE le calcul des noeuds intérieurs dont ne dépendent
    // pas les comms sur la queue ref_queue_inr (_inr = inner)
    Ref<RunQueue> ref_queue_inr = async_node_vector_update(sync_nodes->privateItems());

    // Une fois le packing des données terminé sur ref_queue_bnd
    // on effectue les comms MPI sur CPU, 
    // puis unpacking des données et on synchronise la queue ref_queue_bnd
    ref_sync_req->wait();
    // ici, après cet appel, ref_queue_bnd est synchronisée

    // On attend la terminaison des calculs intérieurs
    ref_queue_inr->barrier();
  } 
  else
  {
    ARCANE_ASSERT(options()->getCcavVectorSyncVersion()==CCAV_VS_nosync,
        ("Ici, pas de synchro"));
    // Pas de synchro
    Ref<RunQueue> ref_queue = async_node_vector_update(allNodes());
    ref_queue->barrier();
  }

  P4GPU_STOP_TIMER(NodeVectorUpdate);
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Implémentation API GPU Arcane version 5 par GG                            */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_computeCqsAndVector_Varcgpu_v5()
{
  constexpr Arcane::Real k025 = 0.25;

  auto queue = m_acc_env->newQueue();
  const Integer nb_cell = allCells().size();
  const Integer nb_node = allNodes().size();
  {
    auto command = makeCommand(queue);

    auto in_node_coord_bis = viewIn(command,m_node_coord_bis);
    //auto out_cell_cqs = viewInOut(command,m_numarray_cqs);
    auto out_cell_cqs = Real3_View8(*m_numarray_cqs);

    auto cnc = m_acc_env->connectivityView().cellNode();

    command << RUNCOMMAND_LOOP1(iter,nb_cell){
      auto [cell_i] = iter();
      CellLocalId cid{(Int32)cell_i};//Int32 cell_i = cid.localId();
      std::array<Real3,8> pos;
      auto nodes = cnc.nodes(cid);
      for( Integer index = 0; index<8; ++index ){
        pos[index] = in_node_coord_bis[nodes[index]];
      }

      out_cell_cqs(0, cell_i) = -k025 * Arcane::math::cross(pos[4] - pos[3], pos[1] - pos[3]);
      out_cell_cqs(1, cell_i) = -k025 * Arcane::math::cross(pos[0] - pos[2], pos[5] - pos[2]);
      out_cell_cqs(2, cell_i) = -k025 * Arcane::math::cross(pos[1] - pos[3], pos[6] - pos[3]);
      out_cell_cqs(3, cell_i) = -k025 * Arcane::math::cross(pos[7] - pos[2], pos[0] - pos[2]);
      out_cell_cqs(4, cell_i) = -k025 * Arcane::math::cross(pos[5] - pos[7], pos[0] - pos[7]);
      out_cell_cqs(5, cell_i) = -k025 * Arcane::math::cross(pos[1] - pos[6], pos[4] - pos[6]);
      out_cell_cqs(6, cell_i) = -k025 * Arcane::math::cross(pos[5] - pos[2], pos[7] - pos[2]);
      out_cell_cqs(7, cell_i) = -k025 * Arcane::math::cross(pos[6] - pos[3], pos[4] - pos[3]);
    };
  }
  {
    // On inverse boucle Cell <-> Node car la boucle originelle sur les mailles n'est parallélisable
    // Du coup, on boucle sur les Node
    // On pourrait construire un groupe de noeuds des mailles active_cells et boucler sur ce groupe
    // Mais ici, on a pré-construit un tableau global au mailles qui indique si une maille fait partie
    // du groupe active_cells ou pas (m_is_active_cell)
    // Rem : en décomp. de dom., pour la plupart des sous-dom. on aura : active_cells = allCells
    // Ainsi, en bouclant sur tous les noeuds allNodes(), pour un noeud donné :
    //   1- on initialise le vecteur à 0
    //   2- on calcule les constributions des mailles connectées au noeud uniquement si elles sont actives
    auto command = makeCommand(queue);

    auto in_cell_arr1 = viewIn(command, m_cell_arr1);
    auto in_cell_arr2 = viewIn(command, m_cell_arr2);
    //auto in_cell_cqs  = viewIn(command, m_numarray_cqs);
    auto in_cell_cqs = Real3_View8(*m_numarray_cqs);
    auto in_is_active_cell = viewIn(command, m_is_active_cell);

    auto out_node_vector = ax::viewOut(command, m_node_vector);

    auto node_index_in_cells = m_acc_env->nodeIndexInCells();
    const Integer max_node_cell = m_acc_env->maxNodeCell();

    auto nc_cty = m_acc_env->connectivityView().nodeCell();

    command << RUNCOMMAND_LOOP1(iter,nb_node){
      auto [node_i] = iter();
      NodeLocalId nid{(Int32)node_i};
      Int32 first_pos = node_i * 8;
      Real3 node_vec;
      auto cells = nc_cty.cells(nid);
      for( Integer index = 0; index<8; ++index ){
        Int32 cid_as_int = cells[index];
        CellLocalId cid { cid_as_int };
        if (in_is_active_cell[cid]) { // la maille ne contribue que si elle est active
          Int16 node_index = node_index_in_cells[first_pos + index];
          node_vec += (in_cell_arr1[cid]+in_cell_arr2[cid]) * in_cell_cqs(node_index,cid_as_int);
        }
      }
      out_node_vector[nid] = node_vec;
    };
  }
}

/*---------------------------------------------------------------------------*/
/* Implémentation Kokkos                                                     */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_computeCqsAndVector_Vkokkos() {

  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans _computeCqsAndVector_Vkokkos";

  m_kokkos_wrapper->computeCqsAndVectorV2();

  PROF_ACC_END;

  // A decommenter pour effectuer les comparaisons numeriques STDENV_VERIF=WRITE / READ
  // m_kokkos_wrapper->syncHostData(allCells(), allNodes(), m_node_vector, m_node_coord_bis, m_cell_cqs, m_cell_arr1, m_cell_arr2);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
computeCqsAndVector() {

  PROF_ACC_BEGIN(__FUNCTION__);

  switch (options()->getComputeCqsVectorVersion()) {
    case CCVV_ori: _computeCqsAndVector_Vori(); break;
    case CCVV_mt: _computeCqsAndVector_Vmt(); break;
    case CCVV_mt_v2: _computeCqsAndVector_Vmt_v2(); break;
    case CCVV_arcgpu_v1: _computeCqsAndVector_Varcgpu_v1(); break;
    case CCVV_arcgpu_v5: _computeCqsAndVector_Varcgpu_v5(); break;
    case CCVV_kokkos: _computeCqsAndVector_Vkokkos(); break;
    default: break;
  };

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

