#include "Pattern4GPUModule.h"

#include <arcane/materials/ComponentPartItemVectorView.h>
#include <arcane/AcceleratorRuntimeInitialisationInfo.h>

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/* UTILITAIRES POUR PREPARER LES CALCULS MULTI-ENVIRONNEMENT SUR GPU         */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Calcul des cell_id globaux : permet d'associer à chaque maille impure (mixte) */
/* l'identifiant de la maille globale                                        */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_computeMultiEnvGlobalCellId() {
  PROF_ACC_BEGIN(__FUNCTION__);

  // Calcul des cell_id globaux 
  CellToAllEnvCellConverter& all_env_cell_converter=*m_allenvcell_converter;
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Integer cell_id = cell.localId();
    m_global_cell[cell] = cell_id;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    if (all_env_cell.nbEnvironment() !=1) {
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;
        m_global_cell[ev] = cell_id;
      }
      // Maille mixte, contient l'opposé du nombre d'environnements
      m_env_id[icell] = -all_env_cell.nbEnvironment();
    } else {
      // Maille pure, cette boucle est de taille 1
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;
        // Cette affectation n'aura lieu qu'une fois
        m_env_id[icell] = ev.environmentId();
      }
    }
  }

  _checkMultiEnvGlobalCellId();
  PROF_ACC_END;
}

void Pattern4GPUModule::
_checkMultiEnvGlobalCellId() {
#ifdef ARCANE_DEBUG
  debug() << "_checkMultiEnvGlobalCellId";

  // Vérification
  ENUMERATE_ENV(ienv, m_mesh_material_mng) {
    IMeshEnvironment* env = *ienv;
    Integer env_id = env->id();
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;
      Cell cell(ev.globalCell());
      ARCANE_ASSERT(cell.localId()==m_global_cell[ev], ("lid differents"));
      AllEnvCell all_env_cell(ev.allEnvCell());
      if (all_env_cell.nbEnvironment()==1) {
        ARCANE_ASSERT(m_env_id[cell]==env_id, ("cell pure : environnement id incorrect dans m_env_id[cell]"));
      } else {
        ARCANE_ASSERT(m_env_id[cell]==-all_env_cell.nbEnvironment(), ("cell mixte : m_env_id[cell] différent de -nbEnvironment()"));
      }
    }
  }
#endif
}

/*---------------------------------------------------------------------------*/
/* A appeler après l'initialisation de la carte des environnements           */ 
/* pour préparer  traitement des environnements sur accélérateur             */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_initEnvForAcc() {

  m_menv_queue = new MultiAsyncRunQueue(m_runner, m_mesh_material_mng->environments().size());

  // construit le tableau multi-env m_global_cell_id et le tableau global m_env_id
  _computeMultiEnvGlobalCellId();

  _updateEnvForAcc();
}

/*---------------------------------------------------------------------------*/
/* Préparer les données multi-envronnement pour l'accélérateur               */
/* A appeler quand la carte des environnements change                        */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_updateEnvForAcc() {
  // "Conseils" accés mémoire
  ENUMERATE_ENV(ienv,m_mesh_material_mng){
    IMeshEnvironment* env = *ienv;
    // La liste des mailles pures évoluant, il faut donner les nouveaux pointeur et taille
    m_acc_mem_adv->setReadMostly(env->pureEnvItems().valueIndexes());
  }
}


/*---------------------------------------------------------------------------*/
/* INITIALISATION DES VARIABLES MULTI-ENVIRONNMENT                           */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Ecriture m_menv_var1 dans m_menv_var1_visu pour visualisation             */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_dumpVisuMEnvVar() {
  if (options()->visuMEnvVar()) {
    m_menv_var1_visu.fill(0.);
#if 1
    ENUMERATE_ENV(ienv, m_mesh_material_mng) {
      IMeshEnvironment* env = *ienv;
      Integer env_id = env->id();

      ENUMERATE_ENVCELL (envcell_i, env) {
        Cell cell = (*envcell_i).globalCell();
        m_menv_var1_visu[cell][env_id] = m_menv_var1[envcell_i];
      }
    }
#else
    CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
    ENUMERATE_CELL(icell, allCells()) {
      Cell cell = * icell;
      AllEnvCell all_env_cell = allenvcell_converter[cell];
      if (all_env_cell.nbEnvironment() !=1) { // uniquement mailles mixtes
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;
          Integer env_id = ev.environmentId();
          m_menv_var1_visu[cell][env_id] = m_menv_var1[ev];
        }
      }
    }
#endif    
  }
}

/*---------------------------------------------------------------------------*/
/* Initialisation des variables multi-envrionnement                          */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
initMEnvVar1() {
  PROF_ACC_BEGIN(__FUNCTION__);

  _initEnvForAcc(); // TODO : en faire un point d'entrée ?

  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();

  if (options()->getInitMenvVar1Version() == IMV1V_ori) {
    CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
    ENUMERATE_CELL(icell, allCells()) {
      Cell cell = * icell;
      const Node& first_node=cell.node(0);
      const Real3& c=node_coord[first_node];
      m_menv_var1[icell]=1.+math::abs(sin(c.x+1)*cos(c.y+1)*sin(c.z+2)); // grandeur moyenne (si maille mixte) et/ou pure (si maille pure)

      AllEnvCell all_env_cell = allenvcell_converter[cell];
      if (all_env_cell.nbEnvironment() !=1) { // uniquement mailles mixtes
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;
          m_menv_var1[ev] = 0.; 
        }
      }
    }
  }
  else if (options()->getInitMenvVar1Version() == IMV1V_arcgpu_v1) 
  {
    auto queue = makeQueue(m_runner);
    queue.setAsync(true);
    {
      auto command = makeCommand(queue);

      auto in_node_coord = ax::viewIn(command, node_coord);
      auto out_menv_var1_g = ax::viewOut(command, m_menv_var1.globalVariable());

      auto cnc = m_connectivity_view.cellNode();

      command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
        NodeLocalId first_nid(cnc.nodes(cid)[0]);
        Real3 c=in_node_coord[first_nid];
        out_menv_var1_g[cid]=1.+math::abs(sin(c.x+1)*cos(c.y+1)*sin(c.z+2));
      };
    }

    // Les calculs des mailles mixtes par environnement sont indépendants
    ENUMERATE_ENV(ienv,m_mesh_material_mng){
      IMeshEnvironment* env = *ienv;

      auto command = makeCommand(m_menv_queue->queue(env->id()));

      Span<Real>          out_menv_var1 (envView(m_menv_var1, env));

      // Nombre de mailles impures (mixtes) de l'environnement
      Integer nb_imp = env->impureEnvItems().nbItem();

      command << RUNCOMMAND_LOOP1(iter, nb_imp) {
        auto [imix] = iter(); // imix \in [0,nb_imp[

        out_menv_var1[imix] = 0;

      }; // asynchrone par rapport au CPU et aux autres environnements
    }
    m_menv_queue->waitAllQueues();
    queue.barrier();
  }

  // Sortie des variables multi-environnement pour la visu
  if (options()->visuMEnvVar()) {
    Integer nb_env = m_mesh_material_mng->environments().size();
    m_menv_var1_visu.resize(nb_env);
  }
  _dumpVisuMEnvVar();
  
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*                               PATTERN 1                                   */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* On calcule seulement les grandeurs partielles sur les mailles mixte       */
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
    ENUMERATE_ENV(ienv,m_mesh_material_mng){
      IMeshEnvironment* env = *ienv;

      auto command = makeCommand(m_menv_queue->queue(env->id()));

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
    m_menv_queue->waitAllQueues();
  }

  _dumpVisuMEnvVar();

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
