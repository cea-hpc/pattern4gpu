#ifndef VAR_SYNC_MNG_H
#define VAR_SYNC_MNG_H

#include "accenv/AcceleratorUtils.h"
#include <arcane/IMesh.h>
#include <arcane/utils/MultiArray2.h>
#include <arcane/IParallelMng.h>

#include "msgpass/SyncItems.h"
#include "msgpass/SyncBuffers.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/* Encapsule une requête de synchronisation sur une variable globale         */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
class GlobalSyncRequest {
 public:
  GlobalSyncRequest(
    SyncItems<ItemType>* sync_items,
    SyncBuffers* sync_buffers,
    IParallelMng* pm,
    Int32ConstArrayView neigh_ranks,
    Ref<RunQueue> ref_queue,
    MeshVarRefT<ItemType, DataType> var);

  virtual ~GlobalSyncRequest();

  // Termine les comms de m_requests, synchronise la queue m_ref_queue, 
  // met à jour les items fantômes de var
  void wait();

 protected:
  bool m_is_over=false;  //! Indique si la requête est terminée
  IParallelMng* m_pm=nullptr;  //! pour effectuer les send/receive proprement dit
  Int32ConstArrayView m_neigh_ranks;  //! Liste des rangs des voisins
  UniqueArray<Parallel::Request> m_requests;  //! les requetes send/receive proprement dites

  ConstMultiArray2View<Integer> m_ghost_item_idx_pn;  //! par voisin, items fantômes à mettre à jour
  
  // Les buffers pour tous les voisins envois/réceptions sur host et device
  MultiBufView m_buf_snd_h;
  MultiBufView m_buf_snd_d;
  MultiBufView m_buf_rcv_h;
  MultiBufView m_buf_rcv_d;

  Ref<RunQueue> m_ref_queue;  //! la queue sur laquelle effectuer les opérations sur GPU
  MeshVarRefT<ItemType, DataType> m_var;  //! variable Arcane dont il faut mettre les items fantômes à jour
};

/*---------------------------------------------------------------------------*/
/* Gère les synchronisations des mailles fantômes par Message Passing        */
/*---------------------------------------------------------------------------*/
class VarSyncMng {
 public:
  VarSyncMng(IMesh* mesh, ax::Runner& runner, AccMemAdviser* acc_mem_adv);
  virtual ~VarSyncMng();

  //! Retourne vrai si un GPU est dispo pour exécuter les calculs
  bool isAcceleratorAvailable() const;

  // Retourne l'instance de SyncItems<T> en fonction de T
  template<typename ItemType>
  SyncItems<ItemType>* getSyncItems();

  // Equivalent à un var.synchronize() où var est une variable globale (i.e. non multi-mat)
  template<typename MeshVariableRefT>
  void globalSynchronize(MeshVariableRefT var);

  // Equivalent à un globalSynchronize pour lequel les données de var sont sur le DEVice
  // La queue asynchrone ref_queue est synchronisé en fin d'appel
  template<typename MeshVariableRefT>
  void globalSynchronizeQueue(Ref<RunQueue> ref_queue, MeshVariableRefT var);

  // Equivalent à un globalSynchronize pour lequel les données de var sont sur le DEVice
  template<typename MeshVariableRefT>
  void globalSynchronizeDevThr(MeshVariableRefT var);

  // Equivalent à un globalSynchronize pour lequel les données de var sont sur le DEVice
  template<typename MeshVariableRefT>
  void globalSynchronizeDevQueues(MeshVariableRefT var);

  // Amorce un globalSynchronizeQueue
  template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
  Ref<GlobalSyncRequest<ItemType, DataType, MeshVarRefT> > iGlobalSynchronizeQueue(Ref<RunQueue> ref_queue, MeshVarRefT<ItemType, DataType> var);

 protected:
  
  // Pré-allocation des buffers de communication pour miniser le nb de réallocations
  void _preAllocBuffers();

 protected:

  ax::Runner& m_runner;

  IParallelMng* m_pm;  //! pour effectuer les send/receive proprement dit

  Integer m_nb_nei;  //! Nb de voisins (m_neigh_ranks.size())
  Int32ConstArrayView m_neigh_ranks;  //! Liste des rangs des voisins

  SyncItems<Cell>* m_sync_cells=nullptr;
  SyncItems<Node>* m_sync_nodes=nullptr;

  SyncBuffers* m_sync_buffers=nullptr;

  MultiAsyncRunQueue* m_neigh_queues=nullptr;  //! Pour gérer plusieurs queues pour les voisins
};

/*---------------------------------------------------------------------------*/
/* Spécialisations pour retourner le nb de fois un type élémentaire DataType */
/* est répété pour un représenter ItemType                                   */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
Integer get_var_degree(MeshVarRefT<ItemType, DataType> var) {
  // Ne devrait jamais être appelé
  ARCANE_ASSERT(false, ("get_var_degree à spécifialiser"));
  return 0;
}
template<typename ItemType, typename DataType>
Integer get_var_degree(MeshVariableScalarRefT<ItemType, DataType> var) {
  return 1;
}

template<typename ItemType, typename DataType>
Integer get_var_degree(MeshVariableArrayRefT<ItemType, DataType> var) {
  return var.arraySize();
}

#endif

