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

using namespace Arcane;
using namespace Arcane::Materials;

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
/* UpdateTensor implementation arcgpu_v2a pour les mailles pures             */
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
_updateTensorPure_arcgpu_v2a()
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
/* UpdateTensor implementation arcgpu_v2a pour les mailles mixtes            */
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
_updateTensorImpure_arcgpu_v2a()
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
/* Implem updateTensor pour arcgpu_v2b en 3D                                 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_updateTensor3D_arcgpu_v2b()
{
  auto queue_pur = m_acc_env->newQueue();
  queue_pur.setAsync(true);
  {
    auto command = makeCommand(queue_pur);

    auto in_env_id       = ax::viewIn   (command, m_env_id);
    // suffixe _p = _pure
    auto in_volume_p     = ax::viewIn   (command, m_volume.globalVariable());
    auto inout_tensor_p  = ax::viewInOut(command, m_tensor.globalVariable());

    command.addKernelName("pure") << RUNCOMMAND_ENUMERATE(Cell, cid, allCells())
    {
      if (in_env_id[cid]>=0) { // vrai ssi cid maille pure
        Real3x3 real3x3 = inout_tensor_p[cid];

        Real vol = in_volume_p[cid];
        real3x3.x.x /= vol;
        real3x3.x.y /= vol;
        real3x3.x.z /= vol;

        real3x3.y.x = real3x3.x.y;
        real3x3.y.y /= vol;
        real3x3.y.z /= vol;

        real3x3.z.x = real3x3.x.z;
        real3x3.z.y = real3x3.y.z;
        real3x3.z.z = -real3x3.x.x - real3x3.y.y;

        inout_tensor_p[cid] = real3x3;
      }
    };
  }

  auto menv_queue = m_acc_env->multiEnvQueue();
  ENUMERATE_ENV(ienv,m_mesh_material_mng){
    IMeshEnvironment* env = *ienv;

    auto command = makeCommand(menv_queue->queue(env->id()));

    Span<const Real>  in_volume   (envView(m_volume, env));
    Span<Real3x3>     inout_tensor(envView(m_tensor, env));

    // Nombre de mailles impures (mixtes) de l'environnement
    Integer nb_imp = env->impureEnvItems().nbItem();

    command << RUNCOMMAND_LOOP1(iter, nb_imp) {
      auto [imix] = iter(); // imix \in [0,nb_imp[

      Real3x3& real3x3 = inout_tensor[imix];

      Real vol = in_volume[imix];
      real3x3.x.x /= vol;
      real3x3.x.y /= vol;
      real3x3.x.z /= vol;

      real3x3.y.x = real3x3.x.y;
      real3x3.y.y /= vol;
      real3x3.y.z /= vol;

      real3x3.z.x = real3x3.x.z;
      real3x3.z.y = real3x3.y.z;
      real3x3.z.z = -real3x3.x.x - real3x3.y.y;

    }; // asynchrone par rapport au CPU et aux autres environnements
  }

  queue_pur.barrier();
  menv_queue->waitAllQueues();
}

/*---------------------------------------------------------------------------*/
/* Implem updateTensor pour arcgpu_v3b en 3D                                 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_updateTensor3D_arcgpu_v3b()
{
  auto queue = m_acc_env->newQueue();
  {
    auto command = makeCommand(queue);

    auto in_volume_g  = ax::viewIn (command, m_volume.globalVariable());
    auto out_tensor_g = ax::viewOut(command, m_tensor.globalVariable());

    MultiEnvVar<Real> menv_volume(m_volume, m_mesh_material_mng);
    auto in_volume(menv_volume.span());

    MultiEnvVar<Real3x3> menv_tensor(m_tensor, m_mesh_material_mng);
    auto inout_tensor(menv_tensor.span());

    // Pour décrire l'accés multi-env sur GPU
    auto in_menv_cell(m_acc_env->multiEnvCellStorage()->viewIn(command));

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {

      Real sum_xx=0., sum_xy=0., sum_xz=0.;
      Real            sum_yy=0., sum_yz=0.;
      for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
        auto evi = in_menv_cell.envCell(cid,ienv);

        Real vol = in_volume[evi];
        Real3x3& real3x3 = inout_tensor.ref(evi); // référence sur la valeur partielle

        sum_xx += real3x3.x.x;
        sum_xy += real3x3.x.y;
        sum_xz += real3x3.x.z;

        sum_yy += real3x3.y.y;
        sum_yz += real3x3.y.z;

        real3x3.x.x /= vol;
        real3x3.x.y /= vol;
        real3x3.x.z /= vol;

        real3x3.y.x = real3x3.x.y;
        real3x3.y.y /= vol;
        real3x3.y.z /= vol;

        real3x3.z.x = real3x3.x.z;
        real3x3.z.y = real3x3.y.z;
        real3x3.z.z = -real3x3.x.x - real3x3.y.y;
      }

      // Valeurs moyennes uniquement sur les mailles mixtes
      if (in_menv_cell.nbEnv(cid)>1) {
        Real vol_glob = in_volume_g[cid];
        Real3x3 tens_glob; 

        tens_glob.x.x = sum_xx/vol_glob;
        tens_glob.x.y = sum_xy/vol_glob;
        tens_glob.x.z = sum_xz/vol_glob;

        tens_glob.y.x = tens_glob.x.y;
        tens_glob.y.y = sum_yy/vol_glob;
        tens_glob.y.z = sum_yz/vol_glob;

        tens_glob.z.x = tens_glob.x.z;
        tens_glob.z.y = tens_glob.y.z;
        tens_glob.z.z = -tens_glob.x.x - tens_glob.y.y;

        out_tensor_g[cid] = tens_glob;
      }
    };
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
  else if (options()->getUpdateTensorVersion() == UVV_arcgpu_v2a)
  {
    // Même résultats numériques que ori
    // Mais _updateVariableV2 ne calcule que les valeurs partielles

    m_acc_env->checkMultiEnvGlobalCellId(m_mesh_material_mng);

    _updateTensorPure_arcgpu_v2a();
    _updateTensorImpure_arcgpu_v2a();
  }
  else if (options()->getUpdateTensorVersion() == UVV_arcgpu_v2b)
  {
    // Même résultats numériques que ori
    // On n'utilise pas des tableaux temporaires,
    // on travaille directement sur m_tensor

    m_acc_env->checkMultiEnvGlobalCellId(m_mesh_material_mng);

    if (defaultMesh()->dimension() == 3) 
    {
      _updateTensor3D_arcgpu_v2b();
    } 
    else 
    {
      fatal() << "UVV_arcgpu_v2b non implemente pour dim != 2";
    }
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
  else if (options()->getUpdateTensorVersion() == UVV_arcgpu_v3b)
  {
    // Même résultats numériques que ori_v3
    // On n'utilise pas des tableaux temporaires,
    // on travaille directement sur m_tensor

    m_acc_env->checkMultiEnvGlobalCellId(m_mesh_material_mng);

    if (defaultMesh()->dimension() == 3) 
    {
      _updateTensor3D_arcgpu_v3b();
    } 
    else 
    {
      fatal() << "UVV_arcgpu_v3b non implemente pour dim != 2";
    }
  }
  
  PROF_ACC_END;
}

