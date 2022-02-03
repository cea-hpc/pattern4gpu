#include "Pattern4GPUModule.h"
#include <arcane/materials/MeshBlockBuildInfo.h>
#include <arcane/materials/MeshEnvironmentBuildInfo.h>
#include <arcane/materials/IMeshBlock.h>
#include <arcane/materials/MeshMaterialModifier.h>
#include <arcane/materials/MatItemEnumerator.h>
#include <arcane/materials/ComponentPartItemVectorView.h>
#include <arcane/materials/MaterialVariableBuildInfo.h>
#include <arcane/IMesh.h>
#include <arcane/IParallelMng.h>
#include <arcane/IItemFamily.h>
#include <arcane/utils/ArcaneGlobal.h>
#include <arcane/utils/StringBuilder.h>
#include <arcane/AcceleratorRuntimeInitialisationInfo.h>
#include <arcane/ServiceBuilder.h>

#include "msgpass/VarSyncMng.h"
#include <arcane/IVariableSynchronizer.h>

#define P4GPU_PROFILING // Pour activer le profiling
#include "P4GPUTimer.h"

#include "Pattern4GPU4Kokkos.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Pattern4GPUModule::
Pattern4GPUModule(const ModuleBuildInfo& mbi)
  : ArcanePattern4GPUObject(mbi), 
  m_mesh_material_mng(IMeshMaterialMng::getReference(defaultMesh())),
  m_compxx(MaterialVariableBuildInfo(
        m_mesh_material_mng, "Compxx", IVariable::PTemporary | IVariable::PExecutionDepend)),
  m_compxy(MaterialVariableBuildInfo(
        m_mesh_material_mng, "Compxy", IVariable::PTemporary | IVariable::PExecutionDepend)),
  m_compyy(MaterialVariableBuildInfo(
        m_mesh_material_mng, "Compyy", IVariable::PTemporary | IVariable::PExecutionDepend)),
  m_tmp1(VariableBuildInfo(mesh(), "Tmp1", IVariable::PTemporary | IVariable::PExecutionDepend)),
  m_kokkos_wrapper(nullptr)
{
}

Pattern4GPUModule::
~Pattern4GPUModule() {
  delete m_allenvcell_converter;
  if (m_kokkos_wrapper) {
    delete m_kokkos_wrapper;
    KokkosWrapper::end();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
accBuild()
{
  PROF_ACC_BEGIN(__FUNCTION__);

  m_acc_env = ServiceBuilder<IAccEnv>(subDomain()).getSingleton();
  m_acc_env->initAcc();

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initKokkosWrapper()
{
  // Tout est déporté dans le wrapper pour s'assurer de la séparation de l'environnement Kokkos et de l'environnement Arcane
  m_kokkos_wrapper = new KokkosWrapper();
  m_kokkos_wrapper->init(allCells(), allNodes(), m_is_active_cell, m_acc_env->nodeIndexInCells());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initP4GPU()
{
  PROF_ACC_START_CAPTURE; // la capture du profiling commence réellement ici

  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initP4GPU";

  // On peut créer maintenant l'objet car la composition des environnements
  // est connue car le le pt d'entree GeomEnv.InitGeomEnv a été appelé
  m_allenvcell_converter=new CellToAllEnvCellConverter(m_mesh_material_mng);

  // On impose un pas de temps (pour l'instant, non paramétrable)
  m_global_deltat = 1.e-3;

  // Pour accélérateur
  m_acc_env->initMesh(mesh());

  // init kokkos
  if (options()->getWithKokkos())
    initKokkosWrapper();

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initTensor()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initTensor";

  // Initialise m_tensor de telle sorte qu'il soit symétrique
  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();
  if (options()->getInitTensorVersion() == ITV_ori) 
  {
    ENUMERATE_ENV(ienv, m_mesh_material_mng) {
      IMeshEnvironment* env = *ienv;
      Integer env_id = env->id();
      const Real d=3*env_id;
      ENUMERATE_ENVCELL (envcell_i, env) {
        EnvCell envcell=(*envcell_i);
        Cell cell=envcell.globalCell();
        const Node& first_node=cell.node(0);
        const Real3& c=node_coord[first_node];
        const Real dd=d+0.5*sin(1+c.x+c.y+c.z);

        Real3x3& tens3x3 = m_tensor[envcell_i];
        tens3x3.x.x=dd+1.;       tens3x3.x.y=dd+1.5;      tens3x3.x.z=dd+1.75;
        tens3x3.y.x=tens3x3.x.y; tens3x3.y.y=dd+2.;       tens3x3.y.z=dd+2.5;
        tens3x3.z.x=tens3x3.x.z; tens3x3.z.y=tens3x3.y.z; tens3x3.z.z=-tens3x3.x.x-tens3x3.y.y;
      }
    }
    CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
    const Real3x3 zero_r3x3(Real3x3::zero());

    ENUMERATE_CELL (cell_i, allCells()) {
      AllEnvCell allenvcell = allenvcell_converter[*cell_i];

      Real3x3 tens3x3 = zero_r3x3;
      if (m_volume[cell_i]>0.) {
        ENUMERATE_CELL_ENVCELL (envcell_i, allenvcell) {
          tens3x3 += m_volume[envcell_i]*m_tensor[envcell_i]; // Real3x3 += Real * Real3x3
        }
        tens3x3 /= m_volume[cell_i];
      }
      m_tensor[cell_i] = tens3x3;
    }
  }
  else if (options()->getInitTensorVersion() == ITV_arcgpu_v1)
  {
    auto queue = m_acc_env->newQueue();
    {
      auto command = makeCommand(queue);

      auto in_node_coord  = ax::viewIn(command, node_coord);
      auto in_volume_g    = ax::viewIn(command, m_volume.globalVariable());
      auto out_tensor_g   = ax::viewOut(command, m_tensor.globalVariable());

      MultiEnvVar<Real> menv_volume(m_volume, m_mesh_material_mng);
      auto in_volume(menv_volume.span());

      MultiEnvVar<Real3x3> menv_tensor(m_tensor, m_mesh_material_mng);
      auto inout_tensor(menv_tensor.span());

      // Pour décrire l'accés multi-env sur GPU
      auto in_menv_cell(m_acc_env->multiEnvCellStorage()->viewIn(command));

      auto cnc = m_acc_env->connectivityView().cellNode();

      command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {

        NodeLocalId first_nid(cnc.nodes(cid)[0]);
        Real3 c=in_node_coord[first_nid];
        const Real dc=0.5*sin(1+c.x+c.y+c.z);

        Real vol_glob = in_volume_g[cid];
        Real3x3 tens_glob = Real3x3::zero();

        for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
          auto evi = in_menv_cell.envCell(cid,ienv);
          Integer env_id = in_menv_cell.envId(cid,ienv);
          const Real d=3*env_id;
          const Real dd=d+dc;

          Real3x3& tens3x3 = inout_tensor.ref(evi); // référence sur la valeur partielle

          tens3x3.x.x=dd+1.;       tens3x3.x.y=dd+1.5;      tens3x3.x.z=dd+1.75;
          tens3x3.y.x=tens3x3.x.y; tens3x3.y.y=dd+2.;       tens3x3.y.z=dd+2.5;
          tens3x3.z.x=tens3x3.x.z; tens3x3.z.y=tens3x3.y.z; tens3x3.z.z=-tens3x3.x.x-tens3x3.y.y;

          if (vol_glob>0) {
            Real vol_part = in_volume[evi];
            tens_glob += vol_part*tens3x3; // Real3x3 += Real * Real3x3
          }
        }
        if (vol_glob>0) {
          out_tensor_g[cid] = tens_glob/vol_glob;
        }
      };
    }
  }
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initNodeVector()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initNodeVector";

  // chaque noeud aura un vecteur de norme 1 mais dans des directions
  // différentes
  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();
  if (options()->getInitNodeVectorVersion() == INVV_ori) 
  {
    ENUMERATE_NODE(node_i, allNodes()) {
      const Real3& c=node_coord[node_i];
      Real cos_th=cos(c.x+c.y+c.z); // garantit une valeur dans [-1,+1]
      Real sin_th=math::sqrt(1-cos_th*cos_th); // garantit une valeur dans [0,+1]
      Real phi=(c.x+1)*(c.y+1)*(c.z+1);
      m_node_vector[node_i]=Real3(sin_th*cos(phi), sin_th*sin(phi),cos_th);
    }
  }
  else if (options()->getInitNodeVectorVersion() == INVV_mt) 
  {
    ParallelLoopOptions options;
    options.setPartitioner(ParallelLoopOptions::Partitioner::Static);
    arcaneParallelForeach(allNodes(), options, [&](NodeVectorView nodes){
      ENUMERATE_NODE(node_i, nodes) {
        const Real3& c=node_coord[node_i];
        Real cos_th=cos(c.x+c.y+c.z); // garantit une valeur dans [-1,+1]
        Real sin_th=math::sqrt(1-cos_th*cos_th); // garantit une valeur dans [0,+1]
        Real phi=(c.x+1)*(c.y+1)*(c.z+1);
        m_node_vector[node_i]=Real3(sin_th*cos(phi), sin_th*sin(phi),cos_th);
      }
    });
  }
  else if (options()->getInitNodeVectorVersion() == INVV_arcgpu_v1)
  {
    auto queue = m_acc_env->newQueue();
    auto command = makeCommand(queue);

    auto in_node_coord = ax::viewIn(command, node_coord);
    auto out_node_vector = ax::viewOut(command, m_node_vector);
    
    command << RUNCOMMAND_ENUMERATE(Node, nid, allNodes()) {
      const Real3 c=in_node_coord[nid];
      Real cos_th=cos(c.x+c.y+c.z); // garantit une valeur dans [-1,+1]
      Real sin_th=math::sqrt(1-cos_th*cos_th); // garantit une valeur dans [0,+1]
      Real phi=(c.x+1)*(c.y+1)*(c.z+1);
      out_node_vector[nid]=Real3(sin_th*cos(phi), sin_th*sin(phi),cos_th);
    };
  }
  else if (options()->getInitNodeVectorVersion() == INVV_kokkos)
  {
    m_kokkos_wrapper->initNodeVector(node_coord, allNodes());
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Synchronise les noeuds fantômes pour node_vector                          */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
syncNodeVector() {
  PROF_ACC_BEGIN(__FUNCTION__);

  auto queue = makeQueueRef(m_acc_env->runner());
  queue->setAsync(true);
  auto vsync = m_acc_env->vsyncMng();
  vsync->globalSynchronizeQueue(queue, m_node_vector);

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initNodeCoordBis()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initNodeCoordBis";

  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();

  if (options()->getInitNodeCoordBisVersion() == INCBV_ori)
  {
    ENUMERATE_NODE(node_i, allNodes()) {
      m_node_coord_bis[node_i]=node_coord[node_i];
    }
  }
  else if (options()->getInitNodeCoordBisVersion() == INCBV_arcgpu_v1)
  {
    auto queue = m_acc_env->newQueue();
    auto command = makeCommand(queue);

    auto in_node_coord = ax::viewIn(command, node_coord);
    auto out_node_coord_bis = ax::viewOut(command, m_node_coord_bis);

    command << RUNCOMMAND_ENUMERATE(Node, nid, allNodes()) {
      out_node_coord_bis[nid]=in_node_coord[nid];
    };
  }
  else if (options()->getInitNodeCoordBisVersion() == INCBV_kokkos)
  {
    m_kokkos_wrapper->initNodeCoordBis(node_coord, allNodes());
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initCellArr12()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initCellArr12";

  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();
  if (options()->getInitCellArr12Version() == IA12V_ori) 
  {
    ENUMERATE_CELL(cell_i, allCells()) {
      const Node& first_node=(*cell_i).node(0);
      const Real3& c=node_coord[first_node];
      m_cell_arr1[cell_i]=1.+math::abs(sin(c.x+1)*cos(c.y+1)*sin(c.z+2));
      m_cell_arr2[cell_i]=2.+math::abs(cos(c.x+2)*sin(c.y+1)*cos(c.z+1));
    }
  }
  else if (options()->getInitCellArr12Version() == IA12V_arcgpu_v1)
  {
    auto async_init_cell_arr12 = [&](CellGroup cell_group, bool high_priority=false) -> Ref<RunQueue> {
      ax::RunQueueBuildInfo bi;
      if (high_priority) {
        // 0 = priorité par défaut
        // Plus la valeur de priorité est faible, plus la queue sera prioritaire
        bi.setPriority(-10); 
      }
      auto ref_queue = makeQueueRef(m_acc_env->runner(), bi);
      ref_queue->setAsync(true);
      auto command = makeCommand(ref_queue.get());

      auto in_node_coord = ax::viewIn(command, node_coord);

      auto out_cell_arr1 = ax::viewOut(command, m_cell_arr1);
      auto out_cell_arr2 = ax::viewOut(command, m_cell_arr2);

      auto cnc = m_acc_env->connectivityView().cellNode();

      command << RUNCOMMAND_ENUMERATE(Cell, cid, ownCells()) {
        NodeLocalId first_nid(cnc.nodes(cid)[0]);
        Real3 c=in_node_coord[first_nid];
        out_cell_arr1[cid]=1.+math::abs(sin(c.x+1)*cos(c.y+1)*sin(c.z+2));
        out_cell_arr2[cid]=2.+math::abs(cos(c.x+2)*sin(c.y+1)*cos(c.z+1));
      };

      return ref_queue;
    };

    Integer sync_cell_sync=2; // 0 = pas de synchro
    if (sync_cell_sync==1)
    {
      // Calcul que sur les mailles intérieures
      auto ref_queue = async_init_cell_arr12(ownCells());
      ref_queue->barrier();

      PROF_ACC_BEGIN("syncCellArr12");
#if 0
      m_cell_arr1.synchronize();
      m_cell_arr2.synchronize();
#else
      auto vsync = m_acc_env->vsyncMng();
      vsync->globalSynchronize(m_cell_arr1);
      vsync->globalSynchronize(m_cell_arr2);
#endif
      PROF_ACC_END;
    }
    else if (sync_cell_sync==2)
    {
      auto vsync = m_acc_env->vsyncMng();
      SyncItems<Cell>* sync_cells = vsync->getSyncItems<Cell>();

      auto ref_queue_bnd = async_init_cell_arr12(sync_cells->sharedItems(), /*high_priority=*/true);
      auto ref_queue_inr = async_init_cell_arr12(sync_cells->privateItems());
      // TODO : aggréger les comms de m_cell_arr1 et m_cell_arr2
      vsync->globalSynchronizeQueue(ref_queue_bnd, m_cell_arr1);
      vsync->globalSynchronizeQueue(ref_queue_bnd, m_cell_arr2);
      ref_queue_inr->barrier();
    }
    else
    {
      // Pas de synchro avec les voisins
      auto ref_queue = async_init_cell_arr12(allCells());
      ref_queue->barrier();
    }
  }
  else if (options()->getInitCellArr12Version() == IA12V_kokkos)
  {
    m_kokkos_wrapper->initCellArr12(allCells(), defaultMesh()->nodesCoordinates());
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Initialisation à 0 de la CQS
 */
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initCqs()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initCqs";

  // Valable en 3D, 8 noeuds par maille
  m_cell_cqs.resize(8);

  if (options()->getInitCqsVersion() == ICQV_ori)
  { 
    ENUMERATE_CELL (cell_i, allCells()) {
      for(Integer inode(0) ; inode<8 ; ++inode) {
        m_cell_cqs[cell_i][inode] = Real3::zero();
      }
    }
  }
  else if (options()->getInitCqsVersion() == ICQV_arcgpu_v1)
  {
    auto queue = m_acc_env->newQueue();
    auto command = makeCommand(queue);

    auto out_cell_cqs = ax::viewOut(command, m_cell_cqs);

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
      for(Integer inode(0) ; inode<8 ; ++inode) {
        out_cell_cqs[cid][inode] = Real3::zero();
      }
    };
  }
  else if (options()->getInitCqsVersion() == ICQV_arcgpu_v5)
  {
      m_numarray_cqs = new NumArray<Real3,2>();
      m_numarray_cqs->resize(8,allCells().size());
  }
  else if (options()->getInitCqsVersion() == ICQV_kokkos)
  {
    m_kokkos_wrapper->initCqs();
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Initialisation de la CQS avec des valeurs non nulles
 */
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initCqs1()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initCqs1";

  // Valable en 3D, 8 noeuds par maille
  m_cell_cqs.resize(8);

  if (options()->getInitCqs1Version() == ICQ1V_ori)
  { 
    Real cos_inode[8], sin_inode[8];
    for(Integer inode(0) ; inode<8 ; ++inode) {
      cos_inode[inode] = cos(1+inode);
      sin_inode[inode] = sin(1+inode);
    }
    ENUMERATE_CELL (cell_i, allCells()) {
      Integer cid = cell_i.localId();
      Real sin_cid = sin(1+cid);
      Real cos_cid = cos(1+cid);
      for(Integer inode(0) ; inode<8 ; ++inode) {
        Real cx = 1e-3*(1+math::abs(sin_cid*cos_inode[inode]));
        Real cy = 1e-3*(1+math::abs(sin_cid*sin_inode[inode]));
        Real cz = 1e-3*(1+math::abs(cos_cid*sin_inode[inode]));
        m_cell_cqs[cell_i][inode] = Real3(cx,cy,cz);
      }
    }
  }
  else if (options()->getInitCqs1Version() == ICQ1V_arcgpu_v1)
  {
    auto queue = m_acc_env->newQueue();
    auto command = makeCommand(queue);

    NumArray<Real,1> cos_inode(8);
    NumArray<Real,1> sin_inode(8);
    Span<Real> out_cos_inode(cos_inode.to1DSpan());
    Span<Real> out_sin_inode(sin_inode.to1DSpan());
    for(Integer inode(0) ; inode<8 ; ++inode) {
      out_cos_inode[inode] = cos(1+inode);
      out_sin_inode[inode] = sin(1+inode);
    }

    auto in_cos_inode = ax::viewIn(command, cos_inode);
    auto in_sin_inode = ax::viewIn(command, sin_inode);
    auto out_cell_cqs = ax::viewOut(command, m_cell_cqs);

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
      Real sin_cid = sin(1+cid.localId());
      Real cos_cid = cos(1+cid.localId());
      for(Integer inode(0) ; inode<8 ; ++inode) {
        Real cx = 1e-3*(1+math::abs(sin_cid*in_cos_inode(inode)));
        Real cy = 1e-3*(1+math::abs(sin_cid*in_sin_inode(inode)));
        Real cz = 1e-3*(1+math::abs(cos_cid*in_sin_inode(inode)));
        out_cell_cqs[cid][inode] = Real3(cx,cy,cz);
      }
    };
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Effectue un travail bidon juste pour charger les connectivités sur GPU    */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
burnConnectivity() {
  PROF_ACC_BEGIN(__FUNCTION__);
  auto queue = m_acc_env->newQueue();
  {
    auto command = makeCommand(queue);

    VariableNodeInteger tmp1(VariableBuildInfo(mesh(), "TemporaryTmp1"));
    auto out_tmp1 = ax::viewOut(command, tmp1);

    auto node_index_in_cells = m_acc_env->nodeIndexInCells();
    const Integer max_node_cell = m_acc_env->maxNodeCell();
    auto nc_cty = m_acc_env->connectivityView().nodeCell();

    command << RUNCOMMAND_ENUMERATE(Node, nid, allNodes()) {
      Int32 first_pos = nid.localId() * max_node_cell;
      Integer index = 0;
      Integer sum = 0;
      for( CellLocalId cid : nc_cty.cells(nid) ){
        Int16 node_index = node_index_in_cells[first_pos + index];
        sum += node_index;
        ++index;
      }
      out_tmp1[nid] = sum;
    };
  }
  {
    auto command = makeCommand(queue);

    VariableCellInteger tmp2(VariableBuildInfo(mesh(), "TemporaryTmp2"));
    auto out_tmp2 = ax::viewOut(command, tmp2);

    auto cnc = m_acc_env->connectivityView().cellNode();

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
      Integer index = 0;
      for( NodeLocalId nid : cnc.nodes(cid) ){
        ++index;
      }
      out_tmp2[cid] = index;
    };
  }
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Effectue un travail bidon juste pour charger is_active_cell sur GPU    */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
burnIsActiveCell() {
  PROF_ACC_BEGIN(__FUNCTION__);
  auto queue = m_acc_env->newQueue();
  {
    auto command = makeCommand(queue);

    VariableCellInteger tmp2(VariableBuildInfo(mesh(), "TemporaryTmp2"));
    auto in_is_active_cell = ax::viewIn(command, m_is_active_cell);
    auto out_tmp2 = ax::viewOut(command, tmp2);

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
      out_tmp2[cid] = in_is_active_cell[cid];
    };
  }
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
updateVectorFromTensor() {

  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans updateVectorFromTensor";

  if (options()->getUpdateVectorFromTensorVersion() == UVTV_ori)
  {
    Integer nb_blocks = m_mesh_material_mng->blocks().size();

    for (Integer i = 0; i < nb_blocks; i++) {
      // on récupère le groupe de mailles associé au bloc
      IMeshBlock* b = (m_mesh_material_mng->blocks())[i];
      CellGroup cell_group = b->cells();

      Real3x3 cell_tensor;
      ENUMERATE_CELL (cell_i, cell_group) {
        // boucle sur les noeuds de la maille
        ENUMERATE_NODE (node_i, cell_i->nodes()) {
          cell_tensor = m_tensor[cell_i];
          m_node_vector[node_i] -= math::prodTensVec(cell_tensor,
              m_cell_cqs[cell_i][node_i.index()]);
        }
      }
    }  // end iblock loop
  }
  else if (options()->getUpdateVectorFromTensorVersion() == UVTV_mt)
  {
    Integer nb_blocks = m_mesh_material_mng->blocks().size();

    for (Integer i = 0; i < nb_blocks; i++) {
      // on récupère le groupe de mailles associé au bloc
      IMeshBlock* b = (m_mesh_material_mng->blocks())[i];
      CellGroup cell_group = b->cells();

      auto node_index_in_cells = m_acc_env->nodeIndexInCells();
      const Integer max_node_cell = m_acc_env->maxNodeCell();

      ParallelLoopOptions options;
      options.setPartitioner(ParallelLoopOptions::Partitioner::Auto);

      NodeGroup node_group = allNodes(); // TODO : passer aux noeuds du blocks
      arcaneParallelForeach(node_group, options, [&](NodeVectorView nodes) {
      ENUMERATE_NODE (node_i, nodes) {
        Int32 first_pos = node_i.localId() * max_node_cell;
        ENUMERATE_CELL(cell_i, node_i->cells()) {
          if (true) { // TODO : vrai ssi cell_i est dans cell_group
            Int16 node_index = node_index_in_cells[first_pos + cell_i.index()];
            const Real3x3& cell_tensor = m_tensor[cell_i];
            m_node_vector[node_i] -= math::prodTensVec(cell_tensor,
              m_cell_cqs[cell_i][node_index]);
          }
        }
      }
      });
    }  // end iblock loop
  }
  PROF_ACC_END;
}

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
//  bool is_sync_cell_cqs=true;
  bool is_sync_cell_cqs=false;
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
#if 0
    m_cell_cqs.synchronize();
#else
    m_acc_env->vsyncMng()->globalSynchronize(m_cell_cqs);
#endif
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

  Integer sync_node_vector=2; // 0 = pas de synchro
  if (sync_node_vector == 1) 
  {
    // Calcul sur tous les noeuds "own"
    Ref<RunQueue> ref_queue = async_node_vector_update(ownNodes());
    ref_queue->barrier();

    // Puis comms
    PROF_ACC_BEGIN("syncNodeVector");
#if 0
    m_node_vector.synchronize();
#else
    auto vsync = m_acc_env->vsyncMng();
    //vsync->globalSynchronize(m_node_vector);
    vsync->globalSynchronizeQueue(ref_queue, m_node_vector);
    //vsync->globalSynchronizeDevThr(m_node_vector);
    //vsync->globalSynchronizeDevQueues(m_node_vector);
#endif
    PROF_ACC_END;

  } 
  else if (sync_node_vector == 2) 
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
    vsync->globalSynchronizeQueueEvent(ref_queue_bnd, m_node_vector);
    // ici, après cet appel, ref_queue_bnd est synchronisée

    // On attend la terminaison des calculs intérieurs
    ref_queue_inr->barrier();
  } 
  else if (sync_node_vector == 3) 
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

ARCANE_REGISTER_MODULE_PATTERN4GPU(Pattern4GPUModule);

