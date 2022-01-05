#include "msgpass/VarSyncMng.h"

#include <arcane/IVariableSynchronizer.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ISubDomain.h>
#include <arcane/IParallelMng.h>
#include <arcane/MeshVariableScalarRef.h>

// Retourne le eItemKind en fonction de ItemType
template<typename ItemType>
eItemKind get_item_kind() {
  return IK_Unknown;
}

template<>
eItemKind get_item_kind<Cell>() {
  return IK_Cell;
}

template<>
eItemKind get_item_kind<Node>() {
  return IK_Node;
}

// Retourne le groupe de tous les items d'un ItemType donné
template<typename ItemType>
ItemGroupT<ItemType> get_all_items(IMesh* mesh) {
  return ItemGroupT<ItemType>();
}

template<>
ItemGroupT<Cell> get_all_items(IMesh* mesh) {
  return mesh->allCells();
}

template<>
ItemGroupT<Node> get_all_items(IMesh* mesh) {
  return mesh->allNodes();
}


/*---------------------------------------------------------------------------*/
/* Encapsule la liste des items à envoyer/recevoir pour un type d'item donné */
/*---------------------------------------------------------------------------*/
template<typename ItemType>
SyncItems<ItemType>::SyncItems(IMesh* mesh, Int32ConstArrayView neigh_ranks) 
{
  eItemKind item_kind = get_item_kind<ItemType>();
  IItemFamily* item_family = mesh->itemFamily(item_kind);
  IVariableSynchronizer* var_sync = item_family->allItemsSynchronizer();
  
  Integer nb_nei = neigh_ranks.size();

  // "shared" ou "owned" : les items intérieurs au sous-domaine et qui doivent être envoyés
  // "ghost" : les items fantômes pour lesquels on va recevoir des informations
  m_nb_owned_item_pn.resize(nb_nei);
  m_nb_ghost_item_pn.resize(nb_nei);

  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    m_nb_owned_item_pn[inei] = var_sync->sharedItems(inei).size();
    m_nb_ghost_item_pn[inei] = var_sync->ghostItems(inei).size();
  }

  // On va récupérer les identifiants proprement dits
  ItemGroupT<ItemType> item_group = get_all_items<ItemType>(mesh);
  GroupIndexTable& lid_to_index = *item_group.localIdToIndex().get();

  auto lids2itemidx = [&](Int32ConstArrayView lids, Int32ArrayView item_idx)
  {
    for(Integer ilid=0 ; ilid<lids.size() ; ++ilid) {
      Int32 lid = lids[ilid];
      Integer idx_in_group = lid_to_index[lid];
      ARCANE_ASSERT(idx_in_group != -1, ("idx_in_group == -1"));
      item_idx[ilid] = idx_in_group;
    }
  };

  m_owned_item_idx_pn.resize(m_nb_owned_item_pn);
  m_ghost_item_idx_pn.resize(m_nb_ghost_item_pn);

  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    lids2itemidx(var_sync->sharedItems(inei), m_owned_item_idx_pn[inei]);
    lids2itemidx(var_sync->ghostItems(inei) , m_ghost_item_idx_pn[inei]);
  }
}

/*---------------------------------------------------------------------------*/
/* Encapsule des vues sur plusieurs buffers de communication                 */
/*---------------------------------------------------------------------------*/
template<typename DataType>
MultiBufView<DataType>::MultiBufView()
{ }

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
template<typename DataType>
MultiBufView<DataType>::MultiBufView(ArrayView<Byte*> ptrs, Int64ConstArrayView sizes) :
  m_multi_buf (sizes.size())
{
  ARCANE_ASSERT(ptrs.size()==sizes.size(), ("ptrs.size()!=sizes.size()"));
  Integer nb_buf=sizes.size();
  for(Integer i=0 ; i<nb_buf ; ++i) {

    // ptrs[i] doit être aligné sur alignof(DataType)
    ARCANE_ASSERT(reinterpret_cast<size_t>(ptrs[i])%alignof(DataType)==0, 
        ("L'adresse ptrs[i] n'est pas aligne sur alignof(DataType)"));

    // sizes[i] doit être un multiple de sizeof(DataType)
    ARCANE_ASSERT(sizes[i]%sizeof(DataType)==0, 
        ("sizes[i] n'est pas un multiple de sizeof(DataType)"));

    m_multi_buf[i] = 
      ArrayView<DataType>(static_cast<Integer>(sizes[i]/sizeof(DataType)),
          reinterpret_cast<DataType*>(ptrs[i]));
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
template<typename DataType>
MultiBufView<DataType>::MultiBufView(const MultiBufView<DataType>& rhs) :
  m_multi_buf (rhs.m_multi_buf)
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Accès en lecture/écriture au i-ème buffer
template<typename DataType>
ArrayView<DataType> MultiBufView<DataType>::operator[](Integer i) {
  return m_multi_buf[i];
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Retourne [beg_ptr, end_ptr[ qui contient tous les buffers (peut-être espacés de trous)
template<typename DataType>
Span<Byte> MultiBufView<DataType>::rangeSpan() {
  if (m_multi_buf.size()==0) {
    return Span<Byte>();
  } else {
    Byte* beg_ptr=reinterpret_cast<Byte*>(m_multi_buf[0].data());
    Integer last = m_multi_buf.size()-1;
    Byte* end_ptr=reinterpret_cast<Byte*>(m_multi_buf[last].data()+m_multi_buf[last].size());
    Int64 sz = end_ptr-beg_ptr;
    return Span<Byte>(beg_ptr, sz);
  }
}

/*---------------------------------------------------------------------------*/
/* A partir des items à communiquer, estime une borne sup de la taille du    */ 
/* buffer en octets                                                          */
/*---------------------------------------------------------------------------*/
template<typename DataType>
Int64 SyncBuffers::estimatedMaxBufSz(IntegerConstArrayView item_sizes) {
  Integer nb_nei = item_sizes.size(); // nb de voisins
  Int64 estim_max_buf_sz = 0;
  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    // Dans le pire des cas, le décalage est d'au plus alignof(DataType)-1 octets
    estim_max_buf_sz += (item_sizes[inei]*sizeof(DataType) + alignof(DataType)-1);
  }
  return estim_max_buf_sz;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void SyncBuffers::resetBuf() {
  m_buf_estim_sz = 0;
  m_first_av_pos = 0;
  m_buf_bytes.resize(0); // ou bien .clear() ?
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
template<typename DataType>
void SyncBuffers::addEstimatedMaxSz(ConstMultiArray2View<Integer> item_idx_pn) {
  m_buf_estim_sz += estimatedMaxBufSz<DataType>(item_idx_pn.dim2Sizes());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void SyncBuffers::allocIfNeeded() {
  m_buf_bytes.resize(m_buf_estim_sz);
}

/*---------------------------------------------------------------------------*/
/* A partir du nb d'items par voisin item_sizes et d'un buffer de données déjà
 * alloué buf_bytes,
 * retourne une vue par voisin des buffers
 */
/*---------------------------------------------------------------------------*/
template<typename DataType>
MultiBufView<DataType> SyncBuffers::_multiBufView(
    IntegerConstArrayView item_sizes,
    Span<Byte> buf_bytes) {

  if (estimatedMaxBufSz<DataType>(item_sizes)>buf_bytes.size()) {
    // Il y a un risque que le buffer déjà alloué ne soit pas assez grand
    return MultiBufView<DataType>();
  }

  Integer nb_nei = item_sizes.size(); // nb de voisins
  UniqueArray<Byte*> ptrs(nb_nei); // le pointeur de base du buffer par voisin
  Int64UniqueArray sizes_in_bytes(nb_nei); // la taille en octets du buffer par voisin

  Byte* cur_ptr{buf_bytes.data()};
  size_t available_space = buf_bytes.size();
  Integer inei;

  for(inei=0 ; available_space>0 && inei<nb_nei ; ++inei) {
    // Par voisin, le tableau de valeurs doit être aligné sur alignof(DataType)
    void* cur_ptr_v = static_cast<void*>(cur_ptr);
    if (std::align(alignof(DataType), sizeof(DataType), cur_ptr_v, available_space)) {

      cur_ptr = static_cast<Byte*>(cur_ptr_v); // cur_ptr_v a été potentiellement modifié

      // Ici, cur_ptr a été modifié et est aligné sur alignof(DataType)
      // available_space a été diminué du nb d'octets = cur_ptr(après appel) - cur_ptr(avant appel)

      // Calcul en octets de l'occupation des valeurs pour le voisin inei
      size_t sz_nei_in_bytes = item_sizes[inei]*sizeof(DataType);

      ptrs[inei] = cur_ptr;
      sizes_in_bytes[inei] = sz_nei_in_bytes;

      cur_ptr += sz_nei_in_bytes; // ici, cur_ptr n'est plus forcement aligné avec alignof(T)
      if (sz_nei_in_bytes <= available_space) {
        available_space -= sz_nei_in_bytes;
      } else {
        ARCANE_ASSERT(false, ("Espace insuffisant pour aligner les données dans le buffer, available_space va devenir négatif"));
        break; // available_space ne pourra jamais être négatif car size_t est non signé
      }
    } else {
      ARCANE_ASSERT(false, ("Espace insuffisant pour aligner les données dans le buffer d'après std::align"));
      break;
    }
  }

  if (inei==nb_nei) {
    MultiBufView<DataType> mb(ptrs, sizes_in_bytes);
    return mb;
  } else {
    // On ne devait jamais arriver là
    ARCANE_ASSERT(false, ("On ne devrait pas etre la"));
    return MultiBufView<DataType>();
  }
}

/*---------------------------------------------------------------------------*/
/* */
/*---------------------------------------------------------------------------*/
template<typename DataType>
MultiBufView<DataType> SyncBuffers::multiBufView(
    ConstMultiArray2View<Integer> item_idx_pn) {

  Int64 av_space = m_buf_bytes.size()-m_first_av_pos;
  Span<Byte> buf_bytes{m_buf_bytes.subView(m_first_av_pos, av_space)};
  auto mb = _multiBufView<DataType>(item_idx_pn.dim2Sizes(), buf_bytes);
  auto rg{mb.rangeSpan()}; // Encapsule [beg_ptr, end_ptr[
  Byte* end_ptr = rg.data()+rg.size();
  m_first_av_pos = (end_ptr - m_buf_bytes.data());
  return mb;
}

/*---------------------------------------------------------------------------*/
/* Gère les synchronisations des mailles fantômes par Message Passing        */
/*---------------------------------------------------------------------------*/
VarSyncMng::VarSyncMng(IMesh* mesh) {
  IItemFamily* cell_family = mesh->cellFamily();
  IVariableSynchronizer* var_sync = cell_family->allItemsSynchronizer();

  m_pm = mesh->parallelMng();

  // Hypothèse, la liste des voisins est la même quelle que soit le type d'item
  // Donc, je peux récupérer celle issue des mailles
  m_neigh_ranks = var_sync->communicatingRanks();
  m_nb_nei = m_neigh_ranks.size();

  m_sync_cells = new SyncItems<Cell>(mesh,m_neigh_ranks);
  m_sync_nodes = new SyncItems<Node>(mesh,m_neigh_ranks);
  m_sync_buffers = new SyncBuffers();
}

VarSyncMng::~VarSyncMng() {
  delete m_sync_cells;
  delete m_sync_nodes;
  delete m_sync_buffers;
}


/*---------------------------------------------------------------------------*/
/* Spécialisations pour retourner l'instance de SyncItems<T> en fonction de T*/
/*---------------------------------------------------------------------------*/
template<>
SyncItems<Cell>* VarSyncMng::_getSyncItems() {
  return m_sync_cells;
}

template<>
SyncItems<Node>* VarSyncMng::_getSyncItems() {
  return m_sync_nodes;
}

/*---------------------------------------------------------------------------*/
/* Equivalent à un var.synchronize() où var est une variable globale         */ 
/* (i.e. non multi-mat)                                                      */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType>
void VarSyncMng::globalSynchronize(MeshVariableScalarRefT<ItemType, DataType> var)
{
  /*
  // Distinguer la lecture de var/construction des buffers des messages proprement dits

  // Remplir les buffers de communications à partir de var (pack ?)
  //  ==> il faut un type pour encapsuler les buffers de communications
  //  Idealement, les données pour tous les voisins devraient être dans un seul buffer contigu
  //  pour faciliter les transferts H<->D
  // Si un seul buffer, attention à avoir des données qui sont alignées (std::align ?)

  // Envoyer/recevoir les données avec MPI

  // Lire les buffers MPI reçus pour écrire les données dans var
  */

  SyncItems<ItemType>* sync_items = _getSyncItems<ItemType>();

  auto owned_item_idx_pn = sync_items->ownedItemIdxPn();
  auto ghost_item_idx_pn = sync_items->ghostItemIdxPn();

  m_sync_buffers->resetBuf();
  // On prévoit une taille max du buffer qui va contenir tous les messages
  m_sync_buffers->addEstimatedMaxSz<DataType>(owned_item_idx_pn);
  m_sync_buffers->addEstimatedMaxSz<DataType>(ghost_item_idx_pn);
  // Le buffer de tous les messages est réalloué si pas assez de place
  m_sync_buffers->allocIfNeeded();

  // On récupère les adresses et tailles des buffers d'envoi et de réception
  auto buf_snd = m_sync_buffers->multiBufView<DataType>(owned_item_idx_pn);
  auto buf_rcv = m_sync_buffers->multiBufView<DataType>(ghost_item_idx_pn);

  // L'échange proprement dit des valeurs de var
  UniqueArray<Parallel::Request> requests;
  auto var_arr = var.asArray();

  // On amorce les réceptions
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce la réception
    auto req_rcv = m_pm->recv(buf_rcv[inei], rank_nei, /*blocking=*/false);
    requests.add(req_rcv);
  }

  // On remplit les buffers sur CPU, TODO : sur GPU
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {

    // On lit les valeurs de var pour les recopier dans le buffer d'envoi
    auto owned_item_idx_nei = owned_item_idx_pn[inei];
    Integer nb_owned_item_idx_nei = owned_item_idx_nei.size();

    ArrayView<DataType> snd_vals = buf_snd[inei];

    for(Integer i=0 ; i<nb_owned_item_idx_nei ; ++i) {
      LocalIdType lid{owned_item_idx_nei[i]};
      snd_vals[i] = var_arr[lid];
    }
  }

  // On amorce les envois
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce l'envoi
    auto req_snd = m_pm->send(buf_snd[inei], rank_nei, /*blocking=*/false);
    requests.add(req_snd);
  }

  m_pm->waitAllRequests(requests);
  requests.clear();

  // On recopie les valeurs reçues dans les buffers dans var
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    auto ghost_item_idx_nei = ghost_item_idx_pn[inei];
    Integer nb_ghost_item_idx_nei = ghost_item_idx_nei.size();

    ArrayView<DataType> rcv_vals = buf_rcv[inei];

    for(Integer i=0 ; i<nb_ghost_item_idx_nei ; ++i) {
      LocalIdType lid{ghost_item_idx_nei[i]};
      var_arr[lid] = rcv_vals[i];
    }
  }
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/
#include <arcane/utils/Real3.h>

#define INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE(__ItemType__, __DataType__) \
  template void VarSyncMng::globalSynchronize(MeshVariableScalarRefT<__ItemType__,__DataType__> var)

INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE(Cell, Real);
INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE(Node, Real3);

