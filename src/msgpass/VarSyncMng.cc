#include "msgpass/VarSyncMng.h"
#include "msgpass/PackTransfer.h"

#include <arcane/IVariableSynchronizer.h>
#include <arcane/IItemFamily.h>
#include <arcane/IParallelMng.h>

// Définie ailleurs
bool is_comm_device_aware();

/*---------------------------------------------------------------------------*/
/* Gère les synchronisations des mailles fantômes par Message Passing        */
/*---------------------------------------------------------------------------*/
VarSyncMng::VarSyncMng(IMesh* mesh, ax::Runner& runner, AccMemAdviser* acc_mem_adv) :
  m_mesh   (mesh),
  m_runner (runner)
{
  m_is_device_aware = is_comm_device_aware();

  IItemFamily* cell_family = m_mesh->cellFamily();
  IVariableSynchronizer* var_sync = cell_family->allItemsSynchronizer();

  m_pm = m_mesh->parallelMng();

  // Hypothèse, la liste des voisins est la même quelle que soit le type d'item
  // Donc, je peux récupérer celle issue des mailles
  m_neigh_ranks = var_sync->communicatingRanks();
  m_nb_nei = m_neigh_ranks.size();

  m_sync_cells = new SyncItems<Cell>(m_mesh,m_neigh_ranks, acc_mem_adv);
  m_sync_nodes = new SyncItems<Node>(m_mesh,m_neigh_ranks, acc_mem_adv);
  m_sync_buffers = new SyncBuffers(isAcceleratorAvailable());
  m_neigh_queues = new MultiAsyncRunQueue(m_runner, m_nb_nei, /*unlimited=*/true);

  _preAllocBuffers();

  m_pack_events.resize(m_nb_nei);
  m_transfer_events.resize(m_nb_nei);
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    if (isAcceleratorAvailable()) {
      m_pack_events[inei] = new AcceleratorEvent();
      m_transfer_events[inei] = new AcceleratorEvent();
    } else {
      m_pack_events[inei] = nullptr;
      m_transfer_events[inei] = nullptr;
    }
  }
  ax::RunQueueBuildInfo bi;
  bi.setPriority(-10); // la priorité doit être la même que celle de la queue qui servira au pack/unpack des buffers de comms
  m_ref_queue_data = makeQueueRef(m_runner, bi);
  m_ref_queue_data->setAsync(true);
}

VarSyncMng::~VarSyncMng() {
  delete m_sync_cells;
  delete m_sync_nodes;
  delete m_sync_buffers;
  delete m_neigh_queues;

  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    delete m_pack_events[inei];
    delete m_transfer_events[inei];
  }
}

/*---------------------------------------------------------------------------*/
/* Retourne vrai si un GPU est dispo pour exécuter les calculs               */
/*---------------------------------------------------------------------------*/
bool VarSyncMng::isAcceleratorAvailable() const {
  return AcceleratorUtils::isAvailable(m_runner);
}

/*---------------------------------------------------------------------------*/
/* Retourne vrai si on peut utiliser les adresses dans DEVICE pour les comms */
/*---------------------------------------------------------------------------*/
bool VarSyncMng::isDeviceAware() const {
  return m_is_device_aware;
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
/* Effectue une première allocation des buffers pour les communications      */
/* Ceci est une pré-allocation pour miniser le nb de réallocations           */
/*---------------------------------------------------------------------------*/
void VarSyncMng::_preAllocBuffers() {
  auto sync_items = getSyncItems<Cell>();
  
  auto owned_item_idx_pn = sync_items->ownedItemIdxPn();
  auto ghost_item_idx_pn = sync_items->ghostItemIdxPn();

  // Alloue pour un buffer de 24 Real3 par maille
  Integer degree = 24;

  m_sync_buffers->resetBuf();
  // On prévoit une taille max du buffer qui va contenir tous les messages
  m_sync_buffers->addEstimatedMaxSz<Real3>(owned_item_idx_pn, degree);
  m_sync_buffers->addEstimatedMaxSz<Real3>(ghost_item_idx_pn, degree);
  // Le buffer de tous les messages est réalloué si pas assez de place
  m_sync_buffers->allocIfNeeded();
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/

#define INST_VAR_SYNC_MNG_GET_SYNC_ITEMS(__ItemType__) \
  template SyncItems<__ItemType__>* VarSyncMng::getSyncItems()

INST_VAR_SYNC_MNG_GET_SYNC_ITEMS(Cell);
INST_VAR_SYNC_MNG_GET_SYNC_ITEMS(Node);

