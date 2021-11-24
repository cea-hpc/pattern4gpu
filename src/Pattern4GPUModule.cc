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
  m_tmp1(VariableBuildInfo(mesh(), "Tmp1", IVariable::PTemporary | IVariable::PExecutionDepend))
{
}

Pattern4GPUModule::
~Pattern4GPUModule() {
  delete m_allenvcell_converter;
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
initP4GPU()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initP4GPU";

  // On peut créer maintenant l'objet car la composition des environnements
  // est connue car le le pt d'entree GeomEnv.InitGeomEnv a été appelé
  m_allenvcell_converter=new CellToAllEnvCellConverter(m_mesh_material_mng);

  // On impose un pas de temps (pour l'instant, non paramétrable)
  m_global_deltat = 1.e-3;

  // Pour accélérateur
  m_acc_env->initMesh(mesh());

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

    Real3x3& tens3x3 = m_tensor[cell_i];
    tens3x3 = zero_r3x3;
    if (m_volume[cell_i]>0.) {
      ENUMERATE_CELL_ENVCELL (envcell_i, allenvcell) {
        tens3x3 += m_volume[envcell_i]*m_tensor[envcell_i]; // Real3x3 += Real * Real3x3
      }
      tens3x3 /= m_volume[cell_i];
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
    auto queue = m_acc_env->newQueue();
    auto command = makeCommand(queue);

    auto in_node_coord = ax::viewIn(command, node_coord);

    auto out_cell_arr1 = ax::viewOut(command, m_cell_arr1);
    auto out_cell_arr2 = ax::viewOut(command, m_cell_arr2);

    auto cnc = m_acc_env->connectivityView().cellNode();

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
      NodeLocalId first_nid(cnc.nodes(cid)[0]);
      Real3 c=in_node_coord[first_nid];
      out_cell_arr1[cid]=1.+math::abs(sin(c.x+1)*cos(c.y+1)*sin(c.z+2));
      out_cell_arr2[cid]=2.+math::abs(cos(c.x+2)*sin(c.y+1)*cos(c.z+1));
    };
  }
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
_updateVariable(const MaterialVariableCellReal& volume, MaterialVariableCellReal& f)
{
  PROF_ACC_BEGIN(__FUNCTION__);
  CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
  ENUMERATE_CELL (cell_i, allCells()) {
    AllEnvCell allenvcell = allenvcell_converter[*cell_i];

    if (allenvcell.nbEnvironment() <= 1) continue;

    Real sum = 0.;
    ENUMERATE_CELL_ENVCELL (envcell_i, allenvcell) {
      sum += f[envcell_i];
    }
    f[cell_i] = sum / volume[cell_i];
  }

  ENUMERATE_ENV (env_i, m_mesh_material_mng) {
    ENUMERATE_ENVCELL (envcell_i, *env_i) {
      f[envcell_i] /= volume[envcell_i];
    }
  }
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Ne met à jour que les valeurs partielles                                  */
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
_updateVariableV2(const MaterialVariableCellReal& volume, MaterialVariableCellReal& f)
{
  PROF_ACC_BEGIN(__FUNCTION__);
  
  ENUMERATE_ENV (env_i, m_mesh_material_mng) {
    ENUMERATE_ENVCELL (envcell_i, *env_i) {
      f[envcell_i] /= volume[envcell_i];
    }
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
/* UpdateTensor implementation arcgpu_v2 pour les mailles pures              */
/*---------------------------------------------------------------------------*/
Ref<ax::RunQueue> Pattern4GPUModule::
_asyncUpdateVariableV2Pur(const char* kernel_name,
    const MaterialVariableCellReal& volume, MaterialVariableCellReal& f)
{
  auto queue_ref = makeQueueRef(m_acc_env->runner());
  queue_ref->setAsync(true);
  {
      auto command = makeCommand(queue_ref.get());

      auto in_env_id     = ax::viewIn(command, m_env_id);
      // suffixe _p = _pure
      auto in_volume_p   = ax::viewIn(command, volume.globalVariable());
      auto inout_f_p     = ax::viewInOut(command, f.globalVariable());

      command.addKernelName(kernel_name) << RUNCOMMAND_ENUMERATE(Cell, cid, allCells())
      {
        if (in_env_id[cid]>=0) { // vrai ssi cid maille pure
          Real upd_f = inout_f_p[cid]/in_volume_p[cid];
          inout_f_p[cid] = upd_f; // operateur /= non supporté
        }
      };
  }

  return queue_ref;
}

void Pattern4GPUModule::
_updateTensorPure_arcgpu_v2()
{
  /* Les recopies sont indépendantes par environnement et mailles pures/mixtes
   * On peut donc récupérer toutes les valeurs pures de tous les environnements
   * puis les traiter
   */
  auto queue_pur = m_acc_env->newQueue();
  queue_pur.setAsync(true);
  // Remplir les variables composantes
  {
    auto command = makeCommand(queue_pur);

    auto in_env_id     = ax::viewIn(command, m_env_id);
    // suffixe _p = _pure
    auto in_tensor_p   = ax::viewIn(command, m_tensor.globalVariable());
    auto out_compxx_p  = ax::viewOut(command, m_compxx.globalVariable());
    auto out_compxy_p  = ax::viewOut(command, m_compxy.globalVariable());
    auto out_compyy_p  = ax::viewOut(command, m_compyy.globalVariable());

    command.addKernelName("tens2scal_pur") << RUNCOMMAND_ENUMERATE(Cell, cid, allCells())
    {
      if (in_env_id[cid]>=0) { // vrai ssi cid maille pure
        const Real3x3& real3x3 = in_tensor_p[cid];
        out_compxx_p[cid] = real3x3.x.x;
        out_compxy_p[cid] = real3x3.x.y;
        out_compyy_p[cid] = real3x3.y.y;
      }
    };
  }

  queue_pur.barrier();

  // On doit attendre que les recopies à partir du tenseur soient terminées
  auto queue_xx = _asyncUpdateVariableV2Pur("updxx_pur", m_volume, m_compxx);
  auto queue_xy = _asyncUpdateVariableV2Pur("updxy_pur", m_volume, m_compxy);
  auto queue_yy = _asyncUpdateVariableV2Pur("updyy_pur", m_volume, m_compyy);
  queue_xx->barrier();
  queue_xy->barrier();
  queue_yy->barrier();

  // Maintenant, on peut remettre à jour m_tensor pour les mailles pures
  {
    auto command = makeCommand(queue_pur);

    auto in_env_id     = ax::viewIn(command, m_env_id);
    // suffixe _p = _pure
    auto in_compxx_p    = ax::viewIn(command, m_compxx.globalVariable());
    auto in_compxy_p    = ax::viewIn(command, m_compxy.globalVariable());
    auto in_compyy_p    = ax::viewIn(command, m_compyy.globalVariable());
    auto inout_tensor_p = ax::viewInOut(command, m_tensor.globalVariable());

    command.addKernelName("scal2tens_pur") << RUNCOMMAND_ENUMERATE(Cell, cid, allCells())
    {
      if (in_env_id[cid]>=0) { // vrai ssi cid maille pure
        Real3x3 real3x3 = inout_tensor_p[cid];

        real3x3.x.y = in_compxy_p[cid];
        real3x3.y.x = in_compxy_p[cid];

        real3x3.x.x = in_compxx_p[cid];
        real3x3.y.y = in_compyy_p[cid];
        real3x3.z.z = - in_compxx_p[cid] - in_compyy_p[cid];

        inout_tensor_p[cid] = real3x3;
      }
    };
  }
  queue_pur.barrier();

  if (defaultMesh()->dimension() == 3)
  {
    // Remplir les variables composantes (on utilise les memes variables tampon que pour les champs en 2D)
    // On "lance" les recopies des valeurs pures du tenseur
    {
      auto command = makeCommand(queue_pur);

      auto in_env_id     = ax::viewIn(command, m_env_id);
      // suffixe _p = _pure
      auto in_tensor_p   = ax::viewIn(command, m_tensor.globalVariable());
      auto out_compxx_p  = ax::viewOut(command, m_compxx.globalVariable());
      auto out_compyy_p  = ax::viewOut(command, m_compyy.globalVariable());

      command.addKernelName("tens2scal3d_pur") << RUNCOMMAND_ENUMERATE(Cell, cid, allCells())
      {
        if (in_env_id[cid]>=0) { // vrai ssi cid maille pure
          const Real3x3& real3x3 = in_tensor_p[cid];
          out_compxx_p[cid] = real3x3.x.z;
          out_compyy_p[cid] = real3x3.y.z;
        }
      };
    }
    queue_pur.barrier();

    // Projeter les variables
    auto queue_xz = _asyncUpdateVariableV2Pur("updxz_pur", m_volume, m_compxx);
    auto queue_yz = _asyncUpdateVariableV2Pur("updyz_pur", m_volume, m_compyy);
    queue_xz->barrier();
    queue_yz->barrier();

    // Ranger les variables composantes
    {
      auto command = makeCommand(queue_pur);

      auto in_env_id     = ax::viewIn(command, m_env_id);
      // suffixe _p = _pure
      auto in_compxx_p    = ax::viewIn(command, m_compxx.globalVariable());
      auto in_compyy_p    = ax::viewIn(command, m_compyy.globalVariable());
      auto inout_tensor_p = ax::viewInOut(command, m_tensor.globalVariable());

      command.addKernelName("scal2tens3d_pur") << RUNCOMMAND_ENUMERATE(Cell, cid, allCells())
      {
        if (in_env_id[cid]>=0) { // vrai ssi cid maille pure
          Real3x3 real3x3 = inout_tensor_p[cid];

          real3x3.x.z = in_compxx_p[cid];
          real3x3.y.z = in_compyy_p[cid];
          real3x3.z.x = in_compxx_p[cid];
          real3x3.z.y = in_compyy_p[cid];

          inout_tensor_p[cid] = real3x3;
        }
      };
    }
    queue_pur.barrier();
  }
}

/*---------------------------------------------------------------------------*/
/* UpdateTensor implementation arcgpu_v2 pour les mailles mixtes             */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_asyncUpdateVariableV2Mix(IMeshEnvironment* env,
    MaterialVariableCellReal& volume, MaterialVariableCellReal& f)
{
  auto menv_queue = m_acc_env->multiEnvQueue();
  auto command = makeCommand(menv_queue->queue(env->id()));

  Span<const Real> in_volume(envView(volume, env));
  Span<Real>       inout_f  (envView(f     , env));

  // Nombre de mailles impures (mixtes) de l'environnement
  Integer nb_imp = env->impureEnvItems().nbItem();

  command << RUNCOMMAND_LOOP1(iter, nb_imp) {
    auto [imix] = iter(); // imix \in [0,nb_imp[

    inout_f[imix] /= in_volume[imix];

  }; // asynchrone par rapport au CPU et aux autres environnements
}

void Pattern4GPUModule::
_updateTensorImpure_arcgpu_v2()
{
  // Les calculs des mailles mixtes par environnement sont indépendants
  // Remplir les variables composantes
  auto menv_queue = m_acc_env->multiEnvQueue();
  ENUMERATE_ENV(ienv,m_mesh_material_mng){
    IMeshEnvironment* env = *ienv;

    auto command = makeCommand(menv_queue->queue(env->id()));

    Span<const Real3x3> in_tensor(envView(m_tensor, env));
    Span<Real>         out_compxx(envView(m_compxx, env));
    Span<Real>         out_compxy(envView(m_compxy, env));
    Span<Real>         out_compyy(envView(m_compyy, env));

    // Nombre de mailles impures (mixtes) de l'environnement
    Integer nb_imp = env->impureEnvItems().nbItem();

    command << RUNCOMMAND_LOOP1(iter, nb_imp) {
      auto [imix] = iter(); // imix \in [0,nb_imp[

      const Real3x3& real3x3 = in_tensor[imix];
      out_compxx[imix] = real3x3.x.x;
      out_compxy[imix] = real3x3.x.y;
      out_compyy[imix] = real3x3.y.y;

    }; // asynchrone par rapport au CPU et aux autres environnements
  }
  menv_queue->waitAllQueues();

  ENUMERATE_ENV(ienv,m_mesh_material_mng){
    IMeshEnvironment* env = *ienv;

    _asyncUpdateVariableV2Mix(env, m_volume, m_compxx);
    _asyncUpdateVariableV2Mix(env, m_volume, m_compxy);
    _asyncUpdateVariableV2Mix(env, m_volume, m_compyy);
  }
  menv_queue->waitAllQueues();

  // Maintenant, on peut remettre à jour m_tensor pour les mailles mixtes
  ENUMERATE_ENV(ienv,m_mesh_material_mng){
    IMeshEnvironment* env = *ienv;

    auto command = makeCommand(menv_queue->queue(env->id()));

    Span<const Real>  in_compxx   (envView(m_compxx, env));
    Span<const Real>  in_compxy   (envView(m_compxy, env));
    Span<const Real>  in_compyy   (envView(m_compyy, env));
    Span<Real3x3>     inout_tensor(envView(m_tensor, env));

    // Nombre de mailles impures (mixtes) de l'environnement
    Integer nb_imp = env->impureEnvItems().nbItem();

    command << RUNCOMMAND_LOOP1(iter, nb_imp) {
      auto [imix] = iter(); // imix \in [0,nb_imp[

      Real3x3& real3x3 = inout_tensor[imix];

      real3x3.x.y = in_compxy[imix];
      real3x3.y.x = in_compxy[imix];

      real3x3.x.x = in_compxx[imix];
      real3x3.y.y = in_compyy[imix];
      real3x3.z.z = - in_compxx[imix] - in_compyy[imix];

    }; // asynchrone par rapport au CPU et aux autres environnements
  }
  menv_queue->waitAllQueues();

  if (defaultMesh()->dimension() == 3)
  {
    // Remplir les variables composantes (on utilise les memes variables tampon que pour les champs en 2D)
    // On "lance" les recopies des valeurs mixtes du tenseur
    ENUMERATE_ENV(ienv,m_mesh_material_mng){
      IMeshEnvironment* env = *ienv;

      auto command = makeCommand(menv_queue->queue(env->id()));

      Span<const Real3x3> in_tensor(envView(m_tensor, env));
      Span<Real>         out_compxx(envView(m_compxx, env));
      Span<Real>         out_compyy(envView(m_compyy, env));

      // Nombre de mailles impures (mixtes) de l'environnement
      Integer nb_imp = env->impureEnvItems().nbItem();

      command << RUNCOMMAND_LOOP1(iter, nb_imp) {
        auto [imix] = iter(); // imix \in [0,nb_imp[

        const Real3x3& real3x3 = in_tensor[imix];
        out_compxx[imix] = real3x3.x.z;
        out_compyy[imix] = real3x3.y.z;

      }; // asynchrone par rapport au CPU et aux autres environnements
    }
    menv_queue->waitAllQueues();

    // Projeter les variables
    ENUMERATE_ENV(ienv,m_mesh_material_mng){
      IMeshEnvironment* env = *ienv;

      _asyncUpdateVariableV2Mix(env, m_volume, m_compxx);
      _asyncUpdateVariableV2Mix(env, m_volume, m_compyy);
    }
    menv_queue->waitAllQueues();

    // Ranger les variables composantes
    ENUMERATE_ENV(ienv,m_mesh_material_mng){
      IMeshEnvironment* env = *ienv;

      auto command = makeCommand(menv_queue->queue(env->id()));

      Span<const Real>  in_compxx   (envView(m_compxx, env));
      Span<const Real>  in_compyy   (envView(m_compyy, env));
      Span<Real3x3>     inout_tensor(envView(m_tensor, env));

      // Nombre de mailles impures (mixtes) de l'environnement
      Integer nb_imp = env->impureEnvItems().nbItem();

      command << RUNCOMMAND_LOOP1(iter, nb_imp) {
        auto [imix] = iter(); // imix \in [0,nb_imp[

        Real3x3& real3x3 = inout_tensor[imix];

        real3x3.x.z = in_compxx[imix];
        real3x3.y.z = in_compyy[imix];
        real3x3.z.x = in_compxx[imix];
        real3x3.z.y = in_compyy[imix];

      }; // asynchrone par rapport au CPU et aux autres environnements
    }
    menv_queue->waitAllQueues();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
updateTensor()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans updateTensor";

  if (options()->getUpdateTensorVersion() == UVV_ori)
  {
    // Remplir les variables composantes
    ENUMERATE_ENV (env_i, m_mesh_material_mng) {
      ENUMERATE_ENVCELL (envcell_i, *env_i) {
        const Real3x3& real3x3 = m_tensor[envcell_i];
        m_compxx[envcell_i] = real3x3.x.x;
        m_compxy[envcell_i] = real3x3.x.y;
        m_compyy[envcell_i] = real3x3.y.y;
      }
    }

    _updateVariable(m_volume, m_compxx);
    _updateVariable(m_volume, m_compxy);
    _updateVariable(m_volume, m_compyy);

    // Ranger les variables composantes
    ENUMERATE_ENV (env_i, m_mesh_material_mng) {
      ENUMERATE_ENVCELL (envcell_i, *env_i) {
        Real3x3& real3x3 = m_tensor[envcell_i];

        real3x3.x.y = m_compxy[envcell_i];
        real3x3.y.x = m_compxy[envcell_i];

        real3x3.x.x = m_compxx[envcell_i];
        real3x3.y.y = m_compyy[envcell_i];
        real3x3.z.z = - m_compxx[envcell_i] - m_compyy[envcell_i];
      }
    }

    if (defaultMesh()->dimension() == 3) {
      // Remplir les variables composantes (on utilise les memes variables tampon que pour les champs en 2D)
      ENUMERATE_ENV (env_i, m_mesh_material_mng) {
        ENUMERATE_ENVCELL (envcell_i, *env_i) {
          m_compxx[envcell_i] = m_tensor[envcell_i].x.z;
          m_compyy[envcell_i] = m_tensor[envcell_i].y.z;
        }
      }

      // Projeter les variables
      _updateVariable(m_volume, m_compxx);
      _updateVariable(m_volume, m_compyy);

      // Ranger les variables composantes
      ENUMERATE_ENV (env_i, m_mesh_material_mng) {
        ENUMERATE_ENVCELL (envcell_i, *env_i) {
          Real3x3& real3x3 = m_tensor[envcell_i];

          real3x3.x.z = m_compxx[envcell_i];
          real3x3.y.z = m_compyy[envcell_i];
          real3x3.z.x = m_compxx[envcell_i];
          real3x3.z.y = m_compyy[envcell_i];
        }
      }
    }  // end if (dim == 3)
  }
  else if (options()->getUpdateTensorVersion() == UVV_ori_v2)
  {
    // Même résultats numériques que ori
    // Mais _updateVariableV2 ne calcule que les valeurs partielles

    // Remplir les variables composantes
    ENUMERATE_ENV (env_i, m_mesh_material_mng) {
      ENUMERATE_ENVCELL (envcell_i, *env_i) {
        const Real3x3& real3x3 = m_tensor[envcell_i];
        m_compxx[envcell_i] = real3x3.x.x;
        m_compxy[envcell_i] = real3x3.x.y;
        m_compyy[envcell_i] = real3x3.y.y;
      }
    }

    _updateVariableV2(m_volume, m_compxx);
    _updateVariableV2(m_volume, m_compxy);
    _updateVariableV2(m_volume, m_compyy);

    // Ranger les variables composantes
    ENUMERATE_ENV (env_i, m_mesh_material_mng) {
      ENUMERATE_ENVCELL (envcell_i, *env_i) {
        Real3x3& real3x3 = m_tensor[envcell_i];

        real3x3.x.y = m_compxy[envcell_i];
        real3x3.y.x = m_compxy[envcell_i];

        real3x3.x.x = m_compxx[envcell_i];
        real3x3.y.y = m_compyy[envcell_i];
        real3x3.z.z = - m_compxx[envcell_i] - m_compyy[envcell_i];
      }
    }

    if (defaultMesh()->dimension() == 3) {
      // Remplir les variables composantes (on utilise les memes variables tampon que pour les champs en 2D)
      ENUMERATE_ENV (env_i, m_mesh_material_mng) {
        ENUMERATE_ENVCELL (envcell_i, *env_i) {
          m_compxx[envcell_i] = m_tensor[envcell_i].x.z;
          m_compyy[envcell_i] = m_tensor[envcell_i].y.z;
        }
      }

      // Projeter les variables
      _updateVariableV2(m_volume, m_compxx);
      _updateVariableV2(m_volume, m_compyy);

      // Ranger les variables composantes
      ENUMERATE_ENV (env_i, m_mesh_material_mng) {
        ENUMERATE_ENVCELL (envcell_i, *env_i) {
          Real3x3& real3x3 = m_tensor[envcell_i];

          real3x3.x.z = m_compxx[envcell_i];
          real3x3.y.z = m_compyy[envcell_i];
          real3x3.z.x = m_compxx[envcell_i];
          real3x3.z.y = m_compyy[envcell_i];
        }
      }
    }  // end if (dim == 3)
  }
  else if (options()->getUpdateTensorVersion() == UVV_arcgpu_v2)
  {
    // Même résultats numériques que ori
    // Mais _updateVariableV2 ne calcule que les valeurs partielles

    m_acc_env->checkMultiEnvGlobalCellId(m_mesh_material_mng);

    _updateTensorPure_arcgpu_v2();
    _updateTensorImpure_arcgpu_v2();
  }
  else if (options()->getUpdateTensorVersion() == UVV_ori_v3)
  {
    // On met aussi à jour la valeur moyenne sur les mailles mixtes
    // Donne des résultats numériques différents que ori et ori_v2

    // Remplir les variables composantes
    ENUMERATE_ENV (env_i, m_mesh_material_mng) {
      ENUMERATE_ENVCELL (envcell_i, *env_i) {
        const Real3x3& real3x3 = m_tensor[envcell_i];
        m_compxx[envcell_i] = real3x3.x.x;
        m_compxy[envcell_i] = real3x3.x.y;
        m_compyy[envcell_i] = real3x3.y.y;
      }
    }

    _updateVariable(m_volume, m_compxx);
    _updateVariable(m_volume, m_compxy);
    _updateVariable(m_volume, m_compyy);

    // Ranger les variables composantes
    ENUMERATE_ENV (env_i, m_mesh_material_mng) {
      ENUMERATE_ENVCELL (envcell_i, *env_i) {
        Real3x3& real3x3 = m_tensor[envcell_i];

        real3x3.x.y = m_compxy[envcell_i];
        real3x3.y.x = m_compxy[envcell_i];

        real3x3.x.x = m_compxx[envcell_i];
        real3x3.y.y = m_compyy[envcell_i];
        real3x3.z.z = - m_compxx[envcell_i] - m_compyy[envcell_i];
      }
    }
    // Recopie des valeurs moyennes sur les mailles mixtes
    CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
    ENUMERATE_CELL (cell_i, allCells()) {
      AllEnvCell allenvcell = allenvcell_converter[*cell_i];
      if (allenvcell.nbEnvironment() > 1) {
        Real3x3& real3x3 = m_tensor[cell_i];

        real3x3.x.y = m_compxy[cell_i];
        real3x3.y.x = m_compxy[cell_i];

        real3x3.x.x = m_compxx[cell_i];
        real3x3.y.y = m_compyy[cell_i];
        real3x3.z.z = - m_compxx[cell_i] - m_compyy[cell_i];
      }
    }

    if (defaultMesh()->dimension() == 3) {
      // Remplir les variables composantes (on utilise les memes variables tampon que pour les champs en 2D)
      ENUMERATE_ENV (env_i, m_mesh_material_mng) {
        ENUMERATE_ENVCELL (envcell_i, *env_i) {
          m_compxx[envcell_i] = m_tensor[envcell_i].x.z;
          m_compyy[envcell_i] = m_tensor[envcell_i].y.z;
        }
      }

      // Projeter les variables
      _updateVariable(m_volume, m_compxx);
      _updateVariable(m_volume, m_compyy);

      // Ranger les variables composantes
      ENUMERATE_ENV (env_i, m_mesh_material_mng) {
        ENUMERATE_ENVCELL (envcell_i, *env_i) {
          Real3x3& real3x3 = m_tensor[envcell_i];

          real3x3.x.z = m_compxx[envcell_i];
          real3x3.y.z = m_compyy[envcell_i];
          real3x3.z.x = m_compxx[envcell_i];
          real3x3.z.y = m_compyy[envcell_i];
        }
      }
      // Recopie des valeurs moyennes sur les mailles mixtes
      ENUMERATE_CELL (cell_i, allCells()) {
        AllEnvCell allenvcell = allenvcell_converter[*cell_i];
        if (allenvcell.nbEnvironment() > 1) {
          Real3x3& real3x3 = m_tensor[cell_i];

          real3x3.x.z = m_compxx[cell_i];
          real3x3.y.z = m_compyy[cell_i];
          real3x3.z.x = m_compxx[cell_i];
          real3x3.z.y = m_compyy[cell_i];
        }
      }
    }  // end if (dim == 3)
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

  {
    auto queue = m_acc_env->newQueue();
    auto command = makeCommand(queue);

    auto in_node_coord_bis = ax::viewIn(command,m_node_coord_bis);
    auto out_cell_cqs = ax::viewInOut(command,m_cell_cqs);

    auto cnc = m_acc_env->connectivityView().cellNode();

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()){
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
    auto queue = m_acc_env->newQueue();
    auto command = makeCommand(queue);

    auto in_cell_arr1 = ax::viewIn(command, m_cell_arr1);
    auto in_cell_arr2 = ax::viewIn(command, m_cell_arr2);
    auto in_cell_cqs  = ax::viewIn(command, m_cell_cqs);
    auto in_is_active_cell = ax::viewIn(command, m_is_active_cell);

    auto out_node_vector = ax::viewOut(command, m_node_vector);

    auto node_index_in_cells = m_acc_env->nodeIndexInCells();
    const Integer max_node_cell = m_acc_env->maxNodeCell();

    auto nc_cty = m_acc_env->connectivityView().nodeCell();

    command << RUNCOMMAND_ENUMERATE(Node,nid,allNodes()) {
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
    };
  }

  PROF_ACC_END;
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
  };

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_PATTERN4GPU(Pattern4GPUModule);

