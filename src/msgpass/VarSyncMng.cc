#include "msgpass/VarSyncMng.h"
#include "msgpass/PackTransfer.h"

#include <arcane/IVariableSynchronizer.h>
#include <arcane/IItemFamily.h>
#include <arcane/IParallelMng.h>

/*---------------------------------------------------------------------------*/
/* Gère les synchronisations des mailles fantômes par Message Passing        */
/*---------------------------------------------------------------------------*/
VarSyncMng::VarSyncMng(IMesh* mesh, ax::Runner& runner) :
  m_runner (runner)
{
  IItemFamily* cell_family = mesh->cellFamily();
  IVariableSynchronizer* var_sync = cell_family->allItemsSynchronizer();

  m_pm = mesh->parallelMng();

  // Hypothèse, la liste des voisins est la même quelle que soit le type d'item
  // Donc, je peux récupérer celle issue des mailles
  m_neigh_ranks = var_sync->communicatingRanks();
  m_nb_nei = m_neigh_ranks.size();

  m_sync_cells = new SyncItems<Cell>(mesh,m_neigh_ranks);
  m_sync_nodes = new SyncItems<Node>(mesh,m_neigh_ranks);
  m_sync_buffers = new SyncBuffers(isAcceleratorAvailable());
  m_neigh_queues = new MultiAsyncRunQueue(m_runner, m_nb_nei, /*unlimited=*/true);
}

VarSyncMng::~VarSyncMng() {
  delete m_sync_cells;
  delete m_sync_nodes;
  delete m_sync_buffers;
  delete m_neigh_queues;
}

/*---------------------------------------------------------------------------*/
/* Retourne vrai si un GPU est dispo pour exécuter les calculs               */
/*---------------------------------------------------------------------------*/
bool VarSyncMng::isAcceleratorAvailable() const {
  return ax::impl::isAcceleratorPolicy(m_runner.executionPolicy());
}

/*---------------------------------------------------------------------------*/
/* Spécialisations pour retourner l'instance de SyncItems<T> en fonction de T*/
/*---------------------------------------------------------------------------*/
template<>
SyncItems<Cell>* VarSyncMng::getSyncItems() {
  return m_sync_cells;
}

template<>
SyncItems<Node>* VarSyncMng::getSyncItems() {
  return m_sync_nodes;
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/

#define INST_VAR_SYNC_MNG_GET_SYNC_ITEMS(__ItemType__) \
  template SyncItems<__ItemType__>* VarSyncMng::getSyncItems()

INST_VAR_SYNC_MNG_GET_SYNC_ITEMS(Cell);
INST_VAR_SYNC_MNG_GET_SYNC_ITEMS(Node);

