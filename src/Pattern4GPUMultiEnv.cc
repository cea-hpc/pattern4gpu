#include "Pattern4GPUModule.h"

#include <arcane/materials/ComponentPartItemVectorView.h>
#include <arcane/materials/MeshMaterialVariableSynchronizerList.h>
#include <arcane/AcceleratorRuntimeInitialisationInfo.h>

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/* INITIALISATION DES VARIABLES MULTI-ENVIRONNMENT                           */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Ecriture m_menv_var[123] dans m_menv_var[123]_visu pour visualisation     */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_dumpVisuMEnvVar() {
  if (options()->visuMEnvVar()) {
    m_menv_var1_visu.fill(0.);
    m_menv_var2_visu.fill(0.);
    m_menv_var3_visu.fill(0.);

    ENUMERATE_ENV(ienv, m_mesh_material_mng) {
      IMeshEnvironment* env = *ienv;
      Integer env_id = env->id();

      ENUMERATE_ENVCELL (envcell_i, env) {
        Cell cell = (*envcell_i).globalCell();
        m_menv_var1_visu[cell][env_id] = m_menv_var1[envcell_i];
        m_menv_var2_visu[cell][env_id] = m_menv_var2[envcell_i];
        m_menv_var3_visu[cell][env_id] = m_menv_var3[envcell_i];
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/* Initialisation des variables multi-envrionnement                          */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
initMEnvVar() {
  PROF_ACC_BEGIN(__FUNCTION__);

  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();

  if (options()->getInitMenvVarVersion() == IMVV_ori) {
    bool to_sync = true;
    //bool to_sync = false;
    CellGroup cell_group = (to_sync ? ownCells() : allCells());
    CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
    ENUMERATE_CELL(icell, cell_group) {
      Cell cell = * icell;
      const Node& first_node=cell.node(0);
      const Real3& c=node_coord[first_node];
      // grandeur moyenne (si maille mixte) et/ou pure (si maille pure)
      m_menv_var1[icell]=1.+math::abs(sin(c.x+1)*cos(c.y+1)*sin(c.z+2)); 
      m_menv_var2[icell]=2.+math::abs(cos(c.x+2)*sin(c.y+1)*cos(c.z+1));
      m_menv_var3[icell]=3.+math::abs(sin(c.x+2)*sin(c.y+2)*cos(c.z+1));
      m_menv_iv1[icell] = icell.localId();

      AllEnvCell all_env_cell = allenvcell_converter[cell];
      if (all_env_cell.nbEnvironment() !=1) { // uniquement mailles mixtes
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;
          m_menv_var1[ev] = 0.; 
          m_menv_var2[ev] = m_frac_vol[ev] * m_menv_var2[icell];
          m_menv_var3[ev] = m_frac_vol[ev] * m_menv_var3[icell];
          m_menv_iv1[ev] = icell.localId()+ev.environmentId();
        }
      }
    }
    if (to_sync) {
#if 0
      MeshMaterialVariableSynchronizerList mmvsl(m_mesh_material_mng);
      m_menv_var1.synchronize(mmvsl);
      m_menv_var2.synchronize(mmvsl);
      m_menv_var3.synchronize(mmvsl);
      mmvsl.apply();
#elif 0
      auto ref_queue = m_acc_env->refQueueAsync();
      m_acc_env->vsyncMng()->multiMatSynchronize(m_menv_var1, ref_queue);
      m_acc_env->vsyncMng()->multiMatSynchronize(m_menv_var2, ref_queue);
      m_acc_env->vsyncMng()->multiMatSynchronize(m_menv_var3, ref_queue);
#else
      MeshVariableSynchronizerList mvsl(m_acc_env->vsyncMng());
      mvsl.add(m_menv_var1);
      mvsl.add(m_menv_iv1);
      mvsl.add(m_menv_var2);
      mvsl.add(m_tensor);
      mvsl.add(m_node_vector); // grandeur globale
      mvsl.add(m_menv_var3);
      auto ref_queue = m_acc_env->refQueueAsync();
      m_acc_env->vsyncMng()->synchronize(mvsl, ref_queue);
#endif
    }
  }
  else if (options()->getInitMenvVarVersion() == IMVV_arcgpu_v1) 
  {
    auto queue = m_acc_env->newQueue();
    {
      auto command = makeCommand(queue);

      auto in_node_coord = ax::viewIn(command, node_coord);
      auto out_menv_var1_g = ax::viewOut(command, m_menv_var1.globalVariable());
      auto out_menv_var2_g = ax::viewOut(command, m_menv_var2.globalVariable());
      auto out_menv_var3_g = ax::viewOut(command, m_menv_var3.globalVariable());
      auto out_menv_iv1_g  = ax::viewOut(command, m_menv_iv1 .globalVariable());

      auto cnc = m_acc_env->connectivityView().cellNode();

      command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
        NodeLocalId first_nid(cnc.nodes(cid)[0]);
        Real3 c=in_node_coord[first_nid];
        out_menv_var1_g[cid]=1.+math::abs(sin(c.x+1)*cos(c.y+1)*sin(c.z+2));
        out_menv_var2_g[cid]=2.+math::abs(cos(c.x+2)*sin(c.y+1)*cos(c.z+1));
        out_menv_var3_g[cid]=3.+math::abs(sin(c.x+2)*sin(c.y+2)*cos(c.z+1));
        out_menv_iv1_g[cid]=cid;
      };
    }

    // Les calculs des mailles mixtes par environnement sont indépendants
    auto menv_queue = m_acc_env->multiEnvQueue();
    ENUMERATE_ENV(ienv,m_mesh_material_mng){
      IMeshEnvironment* env = *ienv;
      Integer env_id = env->id();

      auto command = makeCommand(menv_queue->queue(env->id()));

      auto in_menv_var2_g = ax::viewIn(command, m_menv_var2.globalVariable());
      auto in_menv_var3_g = ax::viewIn(command, m_menv_var3.globalVariable());

      Span<const Real>    in_frac_vol   (envView(m_frac_vol , env));
      Span<const Integer> in_global_cell(envView(m_global_cell, env));

      Span<Real>          out_menv_var1 (envView(m_menv_var1, env));
      Span<Real>          out_menv_var2 (envView(m_menv_var2, env));
      Span<Real>          out_menv_var3 (envView(m_menv_var3, env));
      Span<Integer>       out_menv_iv1  (envView(m_menv_iv1 , env));

      // Nombre de mailles impures (mixtes) de l'environnement
      Integer nb_imp = env->impureEnvItems().nbItem();

      command << RUNCOMMAND_LOOP1(iter, nb_imp) {
        auto [imix] = iter(); // imix \in [0,nb_imp[
        CellLocalId cid(in_global_cell[imix]); // on récupère l'identifiant de la maille globale

        out_menv_var1[imix] = 0;
        out_menv_var2[imix] = in_frac_vol[imix] * in_menv_var2_g[cid];
        out_menv_var3[imix] = in_frac_vol[imix] * in_menv_var3_g[cid];
        out_menv_iv1[imix] = cid+env_id;

      }; // asynchrone par rapport au CPU et aux autres environnements
    }
    menv_queue->waitAllQueues();

    MeshVariableSynchronizerList mvsl(m_acc_env->vsyncMng());
    mvsl.add(m_menv_var1);
    mvsl.add(m_menv_iv1);
    mvsl.add(m_menv_var2);
    mvsl.add(m_tensor);
    mvsl.add(m_menv_var3);
    auto ref_queue = m_acc_env->refQueueAsync();
    m_acc_env->vsyncMng()->synchronize(mvsl, ref_queue);
  }

  // Sortie des variables multi-environnement pour la visu
  if (options()->visuMEnvVar()) {
    Integer nb_env = m_mesh_material_mng->environments().size();
    m_menv_var1_visu.resize(nb_env);
    m_menv_var2_visu.resize(nb_env);
    m_menv_var3_visu.resize(nb_env);
  }
  _dumpVisuMEnvVar();
  
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*                               PATTERN 1                                   */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* On calcule seulement les grandeurs partielles sur les mailles mixtes      */
/*   for(env : ) {                                                           */
/*     for(ev : env.impure()) { // Mailles mixtes                            */
/*       gpart[ev] = calcpart(ev, ev.globalCell());                          */
/*     }                                                                     */
/*   }                                                                       */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
partialImpureOnly() {
  PROF_ACC_BEGIN(__FUNCTION__);

  if (options()->getPartialImpureOnlyVersion() == PIOV_ori)
  {
    CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
    ENUMERATE_CELL(icell, allCells()){
      Cell cell = * icell;
      AllEnvCell all_env_cell = allenvcell_converter[cell];
      if (all_env_cell.nbEnvironment() !=1) { // uniquement mailles mixtes
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;
          m_menv_var1[ev] = m_frac_vol[ev] * m_menv_var1[cell];
        }
      }
    }
  }
  else if (options()->getPartialImpureOnlyVersion() == PIOV_arcgpu_v1)
  {
    // Les calculs des mailles mixtes par environnement sont indépendants
    auto menv_queue = m_acc_env->multiEnvQueue();
    ENUMERATE_ENV(ienv,m_mesh_material_mng){
      IMeshEnvironment* env = *ienv;

      auto command = makeCommand(menv_queue->queue(env->id()));

      auto in_menv_var1_g = ax::viewIn(command, m_menv_var1.globalVariable());
      Span<const Real>    in_frac_vol   (envView(m_frac_vol , env));
      Span<const Integer> in_global_cell(envView(m_global_cell, env));
      Span<Real>          out_menv_var1 (envView(m_menv_var1, env));

      // Nombre de mailles impures (mixtes) de l'environnement
      Integer nb_imp = env->impureEnvItems().nbItem();

      command << RUNCOMMAND_LOOP1(iter, nb_imp) {
        auto [imix] = iter(); // imix \in [0,nb_imp[
        CellLocalId cid(in_global_cell[imix]); // on récupère l'identifiant de la maille globale

        out_menv_var1[imix] = in_frac_vol[imix] * in_menv_var1_g[cid];

      }; // asynchrone par rapport au CPU et aux autres environnements
    }
    menv_queue->waitAllQueues();
  }

  _dumpVisuMEnvVar();

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*                               PATTERN 2                                   */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* On calcule les grandeurs partielles sur les mailles de chaque environnement
 * (mailles pures et mixtes)
 *
 *  for(env : ) {
 *    for(ev : env) { // Mailles pures et mixtes  
 *      gpart[ev] = calcpart(ev);  
 *    } 
 *  }
 *                                                                           */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
partialOnly() {
  PROF_ACC_BEGIN(__FUNCTION__);

  if (options()->getPartialOnlyVersion() == POV_ori)
  {
    ENUMERATE_ENV(ienv, m_mesh_material_mng) {
      IMeshEnvironment* env = *ienv;
      ENUMERATE_ENVCELL(iev,env)
      {
        m_menv_var1[iev] = math::sqrt(m_menv_var2[iev]/m_menv_var3[iev]);
      }
    }
  }
  else if (options()->getPartialOnlyVersion() == POV_arcgpu_v1)
  {
    m_acc_env->checkMultiEnvGlobalCellId(m_mesh_material_mng);

    // On traite séquentiellement les environnements, 
    // chaque environnement étant calculé en parallèle
    ENUMERATE_ENV(ienv, m_mesh_material_mng) {
      IMeshEnvironment* env = *ienv;

      // Mailles pures
      auto queue_pur = m_acc_env->newQueue();
      queue_pur.setAsync(true);
      {
        auto command = makeCommand(queue_pur);

        // Nombre de mailles pures de l'environnement
        Integer nb_pur = env->pureEnvItems().nbItem();

        // Pour les mailles pures, valueIndexes() est la liste des ids locaux des mailles
        Span<const Int32> in_cell_id(env->pureEnvItems().valueIndexes());

        // suffixe _p = _pure
        auto in_menv_var2_p         = ax::viewIn(command, m_menv_var2.globalVariable());
        auto in_menv_var3_p         = ax::viewIn(command, m_menv_var3.globalVariable());
        auto out_menv_var1_p        = ax::viewOut(command, m_menv_var1.globalVariable());

        command << RUNCOMMAND_LOOP1(iter, nb_pur) {
          auto [ipur] = iter(); // ipur \in [0,nb_pur[
          CellLocalId cid(in_cell_id[ipur]); // accés indirect à la valeur de la maille

          out_menv_var1_p[cid] = math::sqrt(in_menv_var2_p[cid]/in_menv_var3_p[cid]);

        }; // non-bloquant et asynchrone par rapport au CPU et autres queues
      }

      // Mailles mixtes
      auto queue_mix = m_acc_env->newQueue();
      queue_mix.setAsync(true);
      {
        auto command = makeCommand(queue_mix);

        // Nombre de mailles impures (mixtes) de l'environnement
        Integer nb_imp = env->impureEnvItems().nbItem();

        // suffixe _i = _impure
        Span<const Real> in_menv_var2_i (envView(m_menv_var2, env));
        Span<const Real> in_menv_var3_i (envView(m_menv_var3, env));

        Span<Real> out_menv_var1_i    (envView(m_menv_var1, env));

        command << RUNCOMMAND_LOOP1(iter, nb_imp) {
          auto [imix] = iter(); // imix \in [0,nb_imp[

          out_menv_var1_i[imix] = math::sqrt(in_menv_var2_i[imix]/in_menv_var3_i[imix]);

        }; // non-bloquant et asynchrone par rapport au CPU et autres queues
      }
      queue_pur.barrier();
      queue_mix.barrier();
    }
  }
  else if (options()->getPartialOnlyVersion() == POV_arcgpu_v2)
  {
    m_acc_env->checkMultiEnvGlobalCellId(m_mesh_material_mng);

    // Comme le traitement est uniforme sur tous les environnements
    // on peut traiter toutes les mailles pures de tous les env simultanément
    // On lance de manière asynchrone les calculs des valeurs pures sur GPU sur queue_glob
    auto queue_glob = m_acc_env->newQueue();
    queue_glob.setAsync(true);
    {
      auto command = makeCommand(queue_glob);

      auto in_env_id              = ax::viewIn(command, m_env_id);
      // suffixe _p = _pure
      auto in_menv_var2_p         = ax::viewIn(command, m_menv_var2.globalVariable());
      auto in_menv_var3_p         = ax::viewIn(command, m_menv_var3.globalVariable());
      auto out_menv_var1_p        = ax::viewOut(command, m_menv_var1.globalVariable());

      command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells())
      {
        if (in_env_id[cid]>=0) { // vrai ssi cid maille pure
          out_menv_var1_p[cid] = math::sqrt(in_menv_var2_p[cid]/in_menv_var3_p[cid]);
        }
      };
    }

    // On traite en concurrence les mailles mixtes des environnements, 
    // chaque environnement étant calculé en parallèle
    auto menv_queue = m_acc_env->multiEnvQueue();
    ENUMERATE_ENV(ienv, m_mesh_material_mng) {
      IMeshEnvironment* env = *ienv;

      // Mailles mixtes
      auto command = makeCommand(menv_queue->queue(env->id()));

      // Nombre de mailles impures (mixtes) de l'environnement
      Integer nb_imp = env->impureEnvItems().nbItem();

      // suffixe _i = _impure
      Span<const Real> in_menv_var2_i (envView(m_menv_var2, env));
      Span<const Real> in_menv_var3_i (envView(m_menv_var3, env));

      Span<Real> out_menv_var1_i    (envView(m_menv_var1, env));

      command << RUNCOMMAND_LOOP1(iter, nb_imp) {
        auto [imix] = iter(); // imix \in [0,nb_imp[

        out_menv_var1_i[imix] = math::sqrt(in_menv_var2_i[imix]/in_menv_var3_i[imix]);

      }; // non-bloquant et asynchrone par rapport au CPU et autres queues
    }

    queue_glob.barrier();
    menv_queue->waitAllQueues();
  }
  else if (options()->getPartialOnlyVersion() == POV_arcgpu_v3)
  {
    m_acc_env->checkMultiEnvGlobalCellId(m_mesh_material_mng);
    
    auto queue = m_acc_env->newQueue();
    queue.setAsync(true);
    
    MultiEnvVarHD<Real> menv_menv_var1(m_menv_var1, m_buf_addr_mng);
    MultiEnvVarHD<Real> menv_menv_var2(m_menv_var2, m_buf_addr_mng);
    MultiEnvVarHD<Real> menv_menv_var3(m_menv_var3, m_buf_addr_mng);

    auto end_cpy_event = m_buf_addr_mng->asyncCpyHToD(queue);
    UniqueArray<Ref<ax::RunQueueEvent>> events{end_cpy_event};

    auto out_menv_var1  (menv_menv_var1.spanD());
    auto in_menv_var2   (menv_menv_var2.spanD());
    auto in_menv_var3   (menv_menv_var3.spanD());

    // Description du traitement pour un environnement
    auto comp_var1 = [&](IMeshEnvironment* env, 
        ConstArrayView<EnvVarIndex> levis, RunQueue* async_queue) {
      auto command = makeCommand(async_queue);

      // Mailles de l'environnement
      Span<const EnvVarIndex> in_levis(levis);
      Integer nb_evis = in_levis.size();

      command << RUNCOMMAND_LOOP1(iter, nb_evis) {
        auto [i] = iter(); // i \in [0,nb_evis[
        const auto evi = in_levis[i];

        out_menv_var1.setValue(evi, 
            math::sqrt(in_menv_var2[evi]/in_menv_var3[evi]));
      }; 
    }; // fin lambda comp_var1

    // Effectue à la fois le calcul sur tous les environnements + synchronisation
    m_acc_env->vsyncMng()->enumerateEnvAndSyncOnEvents(events,
        comp_var1, m_menv_var1,
        options()->ponlyVar1SyncVersion()
        );
  }

  _dumpVisuMEnvVar();

  PROF_ACC_END;
}


/*---------------------------------------------------------------------------*/
/*                               PATTERN 3                                   */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* On calcule les grandeurs partielles sur les mailles de chaque environnement
 * (mailles pures et mixtes) et on moyennise sur les mailles mixtes
 *
 *  for(env : ) {
 *    for(ev : env) { // Mailles pures et mixtes  
 *      gpart[ev] = calcpart(ev, ev.globalCell());
 *    } 
 *  }
 *  for(cell : mixtes) {
 *    gmoy[cell] = 0
 *    for(ev : cell.allEnvCell()) {
 *      gmoy[cell] += pond(ev) * gpart[ev];
 *    }
 *  }
 *                                                                           */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
partialAndMean() {
  PROF_ACC_BEGIN(__FUNCTION__);

  if (options()->getPartialAndMeanVersion() == PMV_ori)
  {
    // On calcule les grandeurs partielles env par env
    ENUMERATE_ENV(ienv, m_mesh_material_mng) {
      IMeshEnvironment* env = *ienv;
      ENUMERATE_ENVCELL(iev,env)
      {
        m_menv_var1[iev] = math::sqrt(m_menv_var2[iev]/m_menv_var3[iev]);
      }
    }

    // Puis on effectue la moyenne sur chacune des mailles mixtes
    CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
    ENUMERATE_CELL(icell, allCells()){
      Cell cell = * icell;
      AllEnvCell all_env_cell = allenvcell_converter[cell];
      if (all_env_cell.nbEnvironment() !=1) { // uniquement mailles mixtes
        m_menv_var1[icell] = 0.;
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          m_menv_var1[icell] += m_frac_vol[ienvcell] * m_menv_var1[ienvcell];
        }
      }
    }
  }
  else if (options()->getPartialAndMeanVersion() == PMV_ori_v2)
  {
    CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
    ENUMERATE_CELL(icell, allCells()){
      Cell cell = * icell;
      AllEnvCell all_env_cell = allenvcell_converter[cell];
      // Calcul des valeurs partielles
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        m_menv_var1[ienvcell] = math::sqrt(m_menv_var2[ienvcell]/m_menv_var3[ienvcell]);
      }
      // Puis on moyennise uniquement sur les mailles mixtes
      if (all_env_cell.nbEnvironment() !=1) { // uniquement mailles mixtes
        m_menv_var1[icell] = 0.;
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          m_menv_var1[icell] += m_frac_vol[ienvcell] * m_menv_var1[ienvcell];
        }
      }
    }
  }
  else if (options()->getPartialAndMeanVersion() == PMV_arcgpu_v1)
  {
    // On boucle sur l'intégralité du maillage pour faire 2 choses :
    //  - init à 0 de la grandeur globale (moyenne) pour toutes les mailles
    //  - puis calcul de la valeur partielle pour les mailles pures
    auto queue = m_acc_env->newQueue();
    {
      auto command = makeCommand(queue);

      auto in_env_id         = ax::viewIn(command, m_env_id);
      auto in_menv_var2_g    = ax::viewIn(command, m_menv_var2.globalVariable());
      auto in_menv_var3_g    = ax::viewIn(command, m_menv_var3.globalVariable());
      auto out_menv_var1_g   = ax::viewOut(command, m_menv_var1.globalVariable());

      command << RUNCOMMAND_ENUMERATE(Cell,cid,allCells()) {

        out_menv_var1_g[cid] = 0.; // pour préparer la moyenne

        if (in_env_id[cid]>=0) { // Uniquement maille pure
          out_menv_var1_g[cid] = math::sqrt(in_menv_var2_g[cid]/in_menv_var3_g[cid]);
        }
      };
    }

    // Puis, environnement par environnement, 
    //  - calcul de la valeur partielle sur chaque maille mixte
    //  - contribution de la maille mixte dans la maille globale
    ENUMERATE_ENV(ienv,m_mesh_material_mng){
      IMeshEnvironment* env = *ienv;

      // Les kernels sont lancés environnement par environnement les uns après les autres
      auto command = makeCommand(queue);

      Span<const Integer> in_global_cell    (envView(m_global_cell, env));
      Span<const Real>    in_frac_vol       (envView(m_frac_vol,    env)); 
      Span<const Real>    in_menv_var2      (envView(m_menv_var2,   env)); 
      Span<const Real>    in_menv_var3      (envView(m_menv_var3,   env)); 
      Span<Real>          inout_menv_var1   (envView(m_menv_var1,   env)); 

      auto out_menv_var1_g = ax::viewOut(command, m_menv_var1.globalVariable());

      // Nombre de mailles impures (mixtes) de l'environnement
      Integer nb_imp = env->impureEnvItems().nbItem();

      command << RUNCOMMAND_LOOP1(iter, nb_imp) {
        auto [imix] = iter(); // imix \in [0,nb_imp[
        CellLocalId cid(in_global_cell[imix]); // on récupère l'identifiant de la maille globale

        // Calcul de la valeur partielle
        inout_menv_var1[imix] = math::sqrt(in_menv_var2[imix]/in_menv_var3[imix]);

        // Contribution à la moyenne (globale)
        out_menv_var1_g[cid] += in_frac_vol[imix] * inout_menv_var1[imix];
      };
    }
#if 1
    auto ref_queue = m_acc_env->refQueueAsync();
    m_acc_env->vsyncMng()->multiMatSynchronize(m_menv_var1, ref_queue);
#else
    m_menv_var1.synchronize();
#endif
  }
  else if (options()->getPartialAndMeanVersion() == PMV_arcgpu_v2)
  {
    auto queue = m_acc_env->newQueue();
    queue.setAsync(true);
#if 0
    MultiEnvVar<Real> menv_menv_var1(m_menv_var1, m_mesh_material_mng);
    auto inout_menv_var1(menv_menv_var1.span());

    MultiEnvVar<Real> menv_menv_var2(m_menv_var2, m_mesh_material_mng);
    auto in_menv_var2(menv_menv_var2.span());

    MultiEnvVar<Real> menv_menv_var3(m_menv_var3, m_mesh_material_mng);
    auto in_menv_var3(menv_menv_var3.span());

    MultiEnvVar<Real> menv_frac_vol(m_frac_vol, m_mesh_material_mng);
    auto in_frac_vol(menv_frac_vol.span());
#else
    MultiEnvVarHD<Real> menv_menv_var1(m_menv_var1, m_buf_addr_mng);
    MultiEnvVarHD<Real> menv_menv_var2(m_menv_var2, m_buf_addr_mng);
    MultiEnvVarHD<Real> menv_menv_var3(m_menv_var3, m_buf_addr_mng);
    MultiEnvVarHD<Real> menv_frac_vol (m_frac_vol,  m_buf_addr_mng);

    auto end_cpy_event = m_buf_addr_mng->asyncCpyHToD(queue);
    UniqueArray<Ref<ax::RunQueueEvent>> events{end_cpy_event};

    auto inout_menv_var1(menv_menv_var1.spanD());
    auto in_menv_var2   (menv_menv_var2.spanD());
    auto in_menv_var3   (menv_menv_var3.spanD());
    auto in_frac_vol    (menv_frac_vol .spanD());
#endif

    auto comp_var1 = [&](CellGroup cell_group, RunQueue* async_queue) {
      auto command = makeCommand(async_queue);

      auto in_env_id       = ax::viewIn(command, m_env_id);
      auto out_menv_var1_g = ax::viewOut(command, m_menv_var1.globalVariable());

      // Pour décrire l'accés multi-env sur GPU
      auto in_menv_cell(m_acc_env->multiEnvCellStorage()->viewIn(command));

      command << RUNCOMMAND_ENUMERATE(Cell, cid, cell_group) {

        // Calcul des valeurs partielles pour tous les environnements de la maille
        for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
          auto evi = in_menv_cell.envCell(cid,ienv);

          inout_menv_var1.setValue(evi, 
              math::sqrt(in_menv_var2[evi]/in_menv_var3[evi]));
        }

        // Puis on moyennise uniquement sur les mailles mixtes
        if (in_env_id[cid]<0) { // Maille mixte ou vide
          Real sum_var1=0.;
          for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
            auto evi = in_menv_cell.envCell(cid,ienv);

            sum_var1 += in_frac_vol[evi] * inout_menv_var1[evi];
          }
          out_menv_var1_g[cid] = sum_var1;
        }
      };
    }; // fin lambda comp_var1

#if 0
    m_acc_env->vsyncMng()->computeMatAndSyncOnEvents(events,
        comp_var1, m_menv_var1,
        options()->getPmeanVar1SyncVersion());
#else
    MeshVariableSynchronizerList mvsl(m_acc_env->vsyncMng());
    mvsl.add(m_menv_var1);
//    mvsl.add(m_menv_var2);
//    mvsl.add(m_menv_var3);

    m_acc_env->vsyncMng()->computeMatAndSyncOnEvents(events,
        comp_var1, mvsl,
        options()->getPmeanVar1SyncVersion());
#endif

    queue.barrier();
  }

  _dumpVisuMEnvVar();

  PROF_ACC_END;
}


/*---------------------------------------------------------------------------*/
/*                               PATTERN 4                                   */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* On combine des sommes de grandeurs partielles pour mettre à jour des
 * grandeurs globales
 *
 *  for(cell : allCells()) {
 *    som2 = 0;
 *    for(ev : cell.allEnvCell()) { // ttes les mailles partielles de la maille
 *      all_part2[ev.envId()] = gpart2[ev]/gmoy2[cell];
 *      som2 += all_part2[ev.envId()];
 *    }
 *    som3 = 0; 
 *    for(ev : cell.allEnvCell()) { // ttes les mailles partielles de la maille
 *      all_part2[ev.envId()] *= (som2+1);
 *      gpart3[ev] = all_part2[ev.envId()] * gmoy3[cell];
 *      som3 += all_part2[ev.envId()];
 *    }
 *    gmoy1[cell] = som3;
 *  }
 *                                                                           */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
partialAndMean4() {
  PROF_ACC_BEGIN(__FUNCTION__);

  if (options()->getPartialAndMean4Version() == PM4V_ori)
  {
    debug() << "PM4V_ori";
    Integer nb_env = m_mesh_material_mng->environments().size();
    RealUniqueArray all_part2(nb_env);

    CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
    ENUMERATE_CELL(icell, allCells()) {
      Cell cell = * icell;
      AllEnvCell all_env_cell = allenvcell_converter[cell];

      Real sum2=0.;
      ENUMERATE_CELL_ENVCELL(iev,all_env_cell) {
        Integer env_id = (*iev).environmentId();
        all_part2[env_id] = m_menv_var2[iev]/m_menv_var2[icell];
        sum2 += all_part2[env_id];
      }

      Real sum3=0.;
      ENUMERATE_CELL_ENVCELL(iev,all_env_cell) {
        Integer env_id = (*iev).environmentId();
        all_part2[env_id] *= (sum2+1.);
        m_menv_var3[iev] = all_part2[env_id] * m_menv_var3[icell];
        sum3 += all_part2[env_id];
      }

      m_menv_var1[icell] = sum3;
    }
  }
  else if (options()->getPartialAndMean4Version() == PM4V_ori_v2)
  {
    debug() << "PM4V_ori_v2";

    CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
    ENUMERATE_CELL(icell, allCells()) {
      Cell cell = * icell;
      AllEnvCell all_env_cell = allenvcell_converter[cell];

      Real sum2=0.;
      ENUMERATE_CELL_ENVCELL(iev,all_env_cell) {
        sum2 += m_menv_var2[iev]/m_menv_var2[icell];
      }

      Real sum3=0.;
      ENUMERATE_CELL_ENVCELL(iev,all_env_cell) {
        Real contrib2 = (m_menv_var2[iev]/m_menv_var2[icell])*(sum2+1.);
        m_menv_var3[iev] = contrib2 * m_menv_var3[icell];
        sum3 += contrib2;
      }

      m_menv_var1[icell] = sum3;
    }
  }
  else if (options()->getPartialAndMean4Version() == PM4V_arcgpu_v1)
  {
    debug() << "PM4V_arcgpu_v1";

    // La variable temporaire sum2 va se transformer en tableau global 
    // temporaire au maillage m_tmp1.
    // On va faire une première passe sur toutes les mailles pures/mixtes
    // pour calculer m_tmp1
    auto queue = m_acc_env->newQueue();
    {
      auto command = makeCommand(queue);

      auto in_env_id         = ax::viewIn(command, m_env_id);
      auto in_menv_var2_g    = ax::viewIn(command, m_menv_var2.globalVariable());
      auto out_sum2_g        = ax::viewOut(command, m_tmp1);

      command << RUNCOMMAND_ENUMERATE(Cell,cid,allCells()) {

        out_sum2_g[cid] = 0.; // pour préparer la moyenne

        if (in_env_id[cid]>=0) { // Uniquement maille pure
          // Ici, valeur partielle pure = valeur moyenne
          // je sais, c'est égal à 1 mais ça aurait pu être une autre formule
          out_sum2_g[cid] = in_menv_var2_g[cid]/in_menv_var2_g[cid];
        }
      };
    }

    // Puis, environnement par environnement, 
    //  - calcul de la valeur partielle sur chaque maille mixte
    //  - contribution de la maille mixte dans la maille globale
    ENUMERATE_ENV(ienv,m_mesh_material_mng){
      IMeshEnvironment* env = *ienv;

      // Les kernels sont lancés environnement par environnement les uns après les autres
      auto command = makeCommand(queue);

      Span<const Integer> in_global_cell    (envView(m_global_cell, env));
      Span<const Real>    in_menv_var2      (envView(m_menv_var2,   env)); 

      auto in_menv_var2_g    = ax::viewIn(command, m_menv_var2.globalVariable());
      auto out_sum2_g        = ax::viewOut(command, m_tmp1);

      // Nombre de mailles impures (mixtes) de l'environnement
      Integer nb_imp = env->impureEnvItems().nbItem();

      command << RUNCOMMAND_LOOP1(iter, nb_imp) {
        auto [imix] = iter(); // imix \in [0,nb_imp[
        CellLocalId cid(in_global_cell[imix]); // on récupère l'identifiant de la maille globale

        // Contribution à la moyenne (globale)
        out_sum2_g[cid] += in_menv_var2[imix]/in_menv_var2_g[cid];
      };
    }
    // Fin de calcul m_tmp1 (sum2)
    
    // Maintenant, on peut utiliser m_tmp1 (sum2) pour calculer les valeurs globales de m_env_var1
    // On reprend le même pattern
    {
      auto command = makeCommand(queue);

      auto in_env_id         = ax::viewIn(command, m_env_id);
      auto in_sum2_g         = ax::viewIn(command, m_tmp1);
      auto in_menv_var2_g    = ax::viewIn(command, m_menv_var2.globalVariable());
      auto out_menv_var1_g   = ax::viewOut(command, m_menv_var1.globalVariable());
      auto inout_menv_var3_g = ax::viewInOut(command, m_menv_var3.globalVariable());

      command << RUNCOMMAND_ENUMERATE(Cell,cid,allCells()) {

        out_menv_var1_g[cid] = 0.; // pour préparer la moyenne

        if (in_env_id[cid]>=0) { // Uniquement maille pure
          Real sum2 = in_sum2_g[cid];
          Real contrib2 = (in_menv_var2_g[cid]/in_menv_var2_g[cid])*(sum2+1.);
          inout_menv_var3_g[cid] = contrib2 * inout_menv_var3_g[cid];

          out_menv_var1_g[cid] = contrib2;
        }
      };
    }

    // Puis, environnement par environnement, 
    //  - calcul de la valeur partielle sur chaque maille mixte
    //  - contribution de la maille mixte dans la maille globale
    ENUMERATE_ENV(ienv,m_mesh_material_mng){
      IMeshEnvironment* env = *ienv;

      // Les kernels sont lancés environnement par environnement les uns après les autres
      auto command = makeCommand(queue);

      Span<const Integer> in_global_cell  (envView(m_global_cell, env));
      Span<const Real>    in_menv_var2    (envView(m_menv_var2,   env)); 
      Span<Real>          out_menv_var3   (envView(m_menv_var3,   env)); 

      auto in_menv_var2_g    = ax::viewIn(command, m_menv_var2.globalVariable());
      auto in_menv_var3_g    = ax::viewIn(command, m_menv_var3.globalVariable());
      auto in_sum2_g         = ax::viewIn(command, m_tmp1);
      auto inout_menv_var1_g = ax::viewInOut(command, m_menv_var1.globalVariable());


      // Nombre de mailles impures (mixtes) de l'environnement
      Integer nb_imp = env->impureEnvItems().nbItem();

      command << RUNCOMMAND_LOOP1(iter, nb_imp) {
        auto [imix] = iter(); // imix \in [0,nb_imp[
        CellLocalId cid(in_global_cell[imix]); // on récupère l'identifiant de la maille globale

        // Calcul de valeur partielle et de la contibution à la somme
        Real sum2 = in_sum2_g[cid];
        Real contrib2 = (in_menv_var2[imix]/in_menv_var2_g[cid])*(sum2+1.);
        out_menv_var3[imix] = contrib2 * in_menv_var3_g[cid];

        // Contribution à la moyenne (globale)
        inout_menv_var1_g[cid] += contrib2;
      };
    }
  }
  else if (options()->getPartialAndMean4Version() == PM4V_arcgpu_v2)
  {
    debug() << "PM4V_arcgpu_v2";

    auto queue = m_acc_env->newQueue();
    auto command = makeCommand(queue);

    MultiEnvVar<Real> menv_menv_var2(m_menv_var2, m_mesh_material_mng);
    auto in_menv_var2(menv_menv_var2.span());

    MultiEnvVar<Real> menv_menv_var3(m_menv_var3, m_mesh_material_mng);
    auto out_menv_var3(menv_menv_var3.span());

    auto in_menv_var2_g    = ax::viewIn(command, m_menv_var2.globalVariable());
    auto in_menv_var3_g    = ax::viewIn(command, m_menv_var3.globalVariable());
    auto out_menv_var1_g   = ax::viewOut(command, m_menv_var1.globalVariable());

    // Pour décrire l'accés multi-env sur GPU
    auto in_menv_cell(m_acc_env->multiEnvCellStorage()->viewIn(command));


    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {

      Real sum2=0.;
      for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
        auto evi = in_menv_cell.envCell(cid,ienv);

        sum2 += in_menv_var2[evi]/in_menv_var2_g[cid];
      }

      Real sum3=0.;
      for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
        auto evi = in_menv_cell.envCell(cid,ienv);

        Real contrib2 = (in_menv_var2[evi]/in_menv_var2_g[cid])*(sum2+1.);
        out_menv_var3.setValue(evi, contrib2 * in_menv_var3_g[cid]);
        sum3 += contrib2;
      }

      out_menv_var1_g[cid] = sum3;
    };
  }

  _dumpVisuMEnvVar();

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
