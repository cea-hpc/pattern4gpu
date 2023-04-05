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

#include "msgpass/VarSyncMng.h"
#include <arcane/IVariableSynchronizer.h>

#define P4GPU_PROFILING // Pour activer le profiling
#include "P4GPUTimer.h"

#include "Pattern4GPU4Kokkos.h"

#include "accenv/SingletonIAccEnv.h"


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
  m_kokkos_wrapper(nullptr), m_all_cells_2_all_env_cells(nullptr)
{
}

Pattern4GPUModule::
~Pattern4GPUModule() {
  delete m_allenvcell_converter;
  if (m_kokkos_wrapper) {
    delete m_kokkos_wrapper;
    KokkosWrapper::end();
  }

  delete m_buf_addr_mng;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
accBuild()
{
  PROF_ACC_BEGIN(__FUNCTION__);

  m_acc_env = SingletonIAccEnv::accEnv(subDomain());

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initKokkosWrapper()
{
  // Tout est déporté dans le wrapper pour s'assurer de la séparation de l'environnement Kokkos et de l'environnement Arcane
  m_kokkos_wrapper = new KokkosWrapper();
  m_kokkos_wrapper->init(allCells(), allNodes(), 
      m_acc_env->multiEnvMng()->isActiveCell(), m_acc_env->nodeIndexInCells());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initP4GPU()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initP4GPU";

  // On peut créer maintenant l'objet car la composition des environnements
  // est connue car le le pt d'entree GeomEnv.InitGeomEnv a été appelé
  m_allenvcell_converter=new CellToAllEnvCellConverter(m_mesh_material_mng);

  // Les groupes actives_cells et active_nodes ont été créés par GeomEnv.InitGeomEnv
  CellGroup active_cells = defaultMesh()->cellFamily()->findGroup("active_cells");
  NodeGroup active_nodes = defaultMesh()->nodeFamily()->findGroup("active_nodes");
  m_acc_env->multiEnvMng()->setActiveItemsFromGroups(active_cells, active_nodes);

  // On impose un pas de temps (pour l'instant, non paramétrable)
  m_global_deltat = 1.e-3;

  // init kokkos
  if (options()->getWithKokkos())
    initKokkosWrapper();

  // TEST : pour amortir le cout des allocs pour GPU
  if (!m_buf_addr_mng) {
    m_buf_addr_mng = new BufAddrMng(m_acc_env->runner(), m_mesh_material_mng);
  }

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
      auto in_menv_cell(m_acc_env->multiEnvMng()->viewIn(command));

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
  vsync->globalSynchronize(queue, m_node_vector);

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

  NodeGroup node_group = (options()->getCcavCqsSyncVersion() == VS_nosync ?
      allNodes() : ownNodes());

  if (options()->getInitNodeCoordBisVersion() == INCBV_ori)
  {
    ENUMERATE_NODE(node_i, node_group) {
      m_node_coord_bis[node_i]=node_coord[node_i];
    }
  }
  else if (options()->getInitNodeCoordBisVersion() == INCBV_arcgpu_v1)
  {
    auto queue = m_acc_env->newQueue();
    auto command = makeCommand(queue);

    auto in_node_coord = ax::viewIn(command, node_coord);
    auto out_node_coord_bis = ax::viewOut(command, m_node_coord_bis);

    command << RUNCOMMAND_ENUMERATE(Node, nid, node_group) {
      out_node_coord_bis[nid]=in_node_coord[nid];
    };
  }
  else if (options()->getInitNodeCoordBisVersion() == INCBV_kokkos)
  {
    m_kokkos_wrapper->initNodeCoordBis(node_coord, node_group);
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
    auto async_init_cell_arr12 = [&](CellGroup cell_group, RunQueue* async_queue) {
      auto command = makeCommand(async_queue);

      auto in_node_coord = ax::viewIn(command, node_coord);

      auto out_cell_arr1 = ax::viewOut(command, m_cell_arr1);
      auto out_cell_arr2 = ax::viewOut(command, m_cell_arr2);

      auto cnc = m_acc_env->connectivityView().cellNode();

      command << RUNCOMMAND_ENUMERATE(Cell, cid, cell_group) {
        NodeLocalId first_nid(cnc.nodes(cid)[0]);
        Real3 c=in_node_coord[first_nid];
        out_cell_arr1[cid]=1.+math::abs(sin(c.x+1)*cos(c.y+1)*sin(c.z+2));
        out_cell_arr2[cid]=2.+math::abs(cos(c.x+2)*sin(c.y+1)*cos(c.z+1));
      };
    };

    MeshVariableSynchronizerList mvsl(m_acc_env->vsyncMng());

    mvsl.add(m_cell_arr1);
    mvsl.add(m_cell_arr2);

    m_acc_env->vsyncMng()->computeAndSync(
	ownCells(),
	async_init_cell_arr12,
	mvsl);
  }
  else if (options()->getInitCellArr12Version() == IA12V_kokkos)
  {
    m_kokkos_wrapper->initCellArr12(allCells(), defaultMesh()->nodesCoordinates());
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Initialisation à 0 de la CQS                                              */
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
    m_numarray_cqs = new NumArray<Real3,MD_Dim2>();
    m_numarray_cqs->resize(8,allCells().size());
    
    auto queue = m_acc_env->newQueue();
    auto command = makeCommand(queue);

    auto out_cell_cqs = Real3_View8(*m_numarray_cqs);

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
      Int32 cid_as_int = cid.localId();
      for(Integer inode(0) ; inode<8 ; ++inode) {
        out_cell_cqs(inode, cid_as_int) = Real3::zero();
      }
    };
  }
  else if (options()->getInitCqsVersion() == ICQV_kokkos)
  {
    m_kokkos_wrapper->initCqs();
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Initialisation de la CQS avec des valeurs non nulles                      */
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

    NumArray<Real,MD_Dim1> cos_inode(8);
    NumArray<Real,MD_Dim1> sin_inode(8);
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
	(void)cid; // pour éviter un warning du type "unused-but-set-variable"
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
	(void)nid; // pour éviter un warning du type "unused-but-set-variable"
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
    auto in_is_active_cell = ax::viewIn(command, m_acc_env->multiEnvMng()->isActiveCell());
    auto out_tmp2 = ax::viewOut(command, tmp2);

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
      out_tmp2[cid] = in_is_active_cell[cid];
    };
  }
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_PATTERN4GPU(Pattern4GPUModule);

