#include "msgpass/VarSyncMng.h"

#include <arcane/IVariableSynchronizer.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ISubDomain.h>
#include <arcane/IParallelMng.h>
#include <arcane/MeshVariableScalarRef.h>
#include <arcane/MeshVariableArrayRef.h>

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
MultiBufView::MultiBufView()
{ }

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
MultiBufView::MultiBufView(ArrayView<Byte*> ptrs, Int64ConstArrayView sizes) :
  m_ptrs  (ptrs),
  m_sizes (sizes)
{
  ARCANE_ASSERT(ptrs.size()==sizes.size(), ("ptrs.size()!=sizes.size()"));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
MultiBufView::MultiBufView(const MultiBufView& rhs) :
  m_ptrs  (rhs.m_ptrs),
  m_sizes (rhs.m_sizes)
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Convertit un buffer d'octets en buffer de DataType
template<typename DataType>
ArrayView<DataType> MultiBufView::valBuf(ArrayView<Byte> buf) {
  // buf.data() doit être aligné sur alignof(DataType)
  ARCANE_ASSERT(reinterpret_cast<size_t>(buf.data())%alignof(DataType)==0, 
      ("L'adresse buf.data() n'est pas aligne sur alignof(DataType)"));

  // buf.size() doit être un multiple de sizeof(DataType)
  ARCANE_ASSERT(buf.size()%sizeof(DataType)==0, 
      ("buf.size() n'est pas un multiple de sizeof(DataType)"));

  return ArrayView<DataType>(
      static_cast<Integer>(buf.size()/sizeof(DataType)), 
      reinterpret_cast<DataType*>(buf.data()));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Convertit un buffer d'octets en buffer 2D de DataType dont la taille dans la deuxième dimension est dim2_size
template<typename DataType>
Array2View<DataType> MultiBufView::valBuf2(ArrayView<Byte> buf, Integer dim2_size) {
  // buf.data() doit être aligné sur alignof(DataType)
  ARCANE_ASSERT(reinterpret_cast<size_t>(buf.data())%alignof(DataType)==0, 
      ("L'adresse buf.data() n'est pas aligne sur alignof(DataType)"));

  // buf.size() doit être un multiple de dim2_size*sizeof(DataType)
  ARCANE_ASSERT(buf.size()%(dim2_size*sizeof(DataType))==0, 
      ("buf.size() n'est pas un multiple de dim2_size*sizeof(DataType)"));

  return Array2View<DataType>(reinterpret_cast<DataType*>(buf.data()),
      static_cast<Integer>(buf.size()/(dim2_size*sizeof(DataType))), 
      dim2_size);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Accès en lecture/écriture au i-ème buffer d'octets
ArrayView<Byte> MultiBufView::byteBuf(Integer i) {
  return ArrayView<Byte>(m_sizes[i], m_ptrs[i]);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Retourne [beg_ptr, end_ptr[ qui contient tous les buffers (peut-être espacés de trous)
Span<Byte> MultiBufView::rangeSpan() {
  if (m_ptrs.size()==0) {
    return Span<Byte>();
  } else {
    Byte* beg_ptr=m_ptrs[0];
    Integer last = m_ptrs.size()-1;
    Byte* end_ptr=m_ptrs[last]+m_sizes[last];
    Int64 sz = end_ptr-beg_ptr;
    return Span<Byte>(beg_ptr, sz);
  }
}

/*---------------------------------------------------------------------------*/
/* A partir des items à communiquer, estime une borne sup de la taille du    */ 
/* buffer en octets                                                          */
/*---------------------------------------------------------------------------*/
template<typename DataType>
Int64 SyncBuffers::estimatedMaxBufSz(IntegerConstArrayView item_sizes, 
    Integer degree) {
  Integer nb_nei = item_sizes.size(); // nb de voisins
  Int64 estim_max_buf_sz = 0;
  Int64 sizeof_item = sizeof(DataType)*degree;
  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    // Dans le pire des cas, le décalage est d'au plus alignof(DataType)-1 octets
    estim_max_buf_sz += (item_sizes[inei]*sizeof_item + alignof(DataType)-1);
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
void SyncBuffers::addEstimatedMaxSz(ConstMultiArray2View<Integer> item_idx_pn,
    Integer degree) {
  m_buf_estim_sz += estimatedMaxBufSz<DataType>(item_idx_pn.dim2Sizes(), degree);
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
MultiBufView SyncBuffers::_multiBufView(
    IntegerConstArrayView item_sizes, Integer degree,
    Span<Byte> buf_bytes) {

  if (estimatedMaxBufSz<DataType>(item_sizes, degree)>buf_bytes.size()) {
    // Il y a un risque que le buffer déjà alloué ne soit pas assez grand
    return MultiBufView();
  }

  Integer nb_nei = item_sizes.size(); // nb de voisins
  UniqueArray<Byte*> ptrs(nb_nei); // le pointeur de base du buffer par voisin
  Int64UniqueArray sizes_in_bytes(nb_nei); // la taille en octets du buffer par voisin

  Byte* cur_ptr{buf_bytes.data()};
  size_t available_space = buf_bytes.size();
  size_t sizeof_item = sizeof(DataType)*degree;
  Integer inei;

  for(inei=0 ; available_space>0 && inei<nb_nei ; ++inei) {
    // Par voisin, le tableau de valeurs doit être aligné sur alignof(DataType)
    void* cur_ptr_v = static_cast<void*>(cur_ptr);
    if (std::align(alignof(DataType), sizeof(DataType), cur_ptr_v, available_space)) {

      cur_ptr = static_cast<Byte*>(cur_ptr_v); // cur_ptr_v a été potentiellement modifié

      // Ici, cur_ptr a été modifié et est aligné sur alignof(DataType)
      // available_space a été diminué du nb d'octets = cur_ptr(après appel) - cur_ptr(avant appel)

      // Calcul en octets de l'occupation des valeurs pour le voisin inei
      size_t sz_nei_in_bytes = item_sizes[inei]*sizeof_item;

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
    MultiBufView mb(ptrs, sizes_in_bytes);
    return mb;
  } else {
    // On ne devait jamais arriver là
    ARCANE_ASSERT(false, ("On ne devrait pas etre la"));
    return MultiBufView();
  }
}

/*---------------------------------------------------------------------------*/
/* */
/*---------------------------------------------------------------------------*/
template<typename DataType>
MultiBufView SyncBuffers::multiBufView(
    ConstMultiArray2View<Integer> item_idx_pn, Integer degree) {

  Int64 av_space = m_buf_bytes.size()-m_first_av_pos;
  Span<Byte> buf_bytes{m_buf_bytes.subView(m_first_av_pos, av_space)};
  auto mb = _multiBufView<DataType>(item_idx_pn.dim2Sizes(), degree, buf_bytes);
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

/*---------------------------------------------------------------------------*/
/* cpy_var2buf */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
void cpy_var2buf(IntegerConstArrayView item_idx, 
    MeshVarRefT<ItemType, DataType> var,
    ArrayView<Byte> buf) 
{
  // Ne devrait jamais être appelé
  ARCANE_ASSERT(false, ("cpy_var2buf à spécifialiser"));
}

// Spécialisation pour MeshVariable***Scalar***RefT
template<typename ItemType, typename DataType>
void cpy_var2buf(IntegerConstArrayView item_idx, 
    MeshVariableScalarRefT<ItemType, DataType> var,
    ArrayView<Byte> buf) 
{
  auto var_arr = var.asArray();

  ArrayView<DataType> buf_vals(MultiBufView::valBuf<DataType>(buf));

  Integer nb_item_idx = item_idx.size();

  for(Integer i=0 ; i<nb_item_idx ; ++i) {
    LocalIdType lid{item_idx[i]};
    buf_vals[i] = var_arr[lid];
  }
}

// Spécialisation pour MeshVariable***Array***RefT
template<typename ItemType, typename DataType>
void cpy_var2buf(IntegerConstArrayView item_idx, 
    MeshVariableArrayRefT<ItemType, DataType> var,
    ArrayView<Byte> buf) 
{
  auto var_arr = var.asArray();

  Integer degree = var.arraySize();
  // Vue sur tableau 2D [nb_item][degree]
  Array2View<DataType> buf_vals(MultiBufView::valBuf2<DataType>(buf, degree));

  Integer nb_item_idx = item_idx.size();

  for(Integer i=0 ; i<nb_item_idx ; ++i) {
    LocalIdType lid{item_idx[i]};
    Span<const DataType> in_var_arr  (var_arr[lid]);
    Span<DataType>       out_buf_vals(buf_vals[i]);
    for(Integer j=0 ; j<degree ; ++j) {
      out_buf_vals[j] = in_var_arr[j];
    }
  }
}

/*---------------------------------------------------------------------------*/
/* cpy_buf2var */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
void cpy_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVarRefT<ItemType, DataType> var) 
{
  // Ne devrait jamais être appelé
  ARCANE_ASSERT(false, ("cpy_var2buf à spécifialiser"));
}

// Spécialisation pour MeshVariable***Scalar***RefT
template<typename ItemType, typename DataType>
void cpy_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVariableScalarRefT<ItemType, DataType> var)
{
  auto var_arr = var.asArray();

  ConstArrayView<DataType> buf_vals(MultiBufView::valBuf<DataType>(buf));

  Integer nb_item_idx = item_idx.size();

  for(Integer i=0 ; i<nb_item_idx ; ++i) {
    LocalIdType lid{item_idx[i]};
    var_arr[lid] = buf_vals[i];
  }
}

// Spécialisation pour MeshVariable***Array***RefT
template<typename ItemType, typename DataType>
void cpy_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVariableArrayRefT<ItemType, DataType> var)
{
  auto var_arr = var.asArray();

  Integer degree = var.arraySize();
  // Vue sur tableau 2D [nb_item][degree]
  ConstArray2View<DataType> buf_vals(MultiBufView::valBuf2<DataType>(buf, degree));


  Integer nb_item_idx = item_idx.size();

  for(Integer i=0 ; i<nb_item_idx ; ++i) {
    LocalIdType lid{item_idx[i]};
    Span<const DataType> in_buf_vals(buf_vals[i]);
    Span<DataType>       out_var_arr(var_arr[lid]);
    for(Integer j=0 ; j<degree ; ++j) {
      out_var_arr[j] = in_buf_vals[j];
    }
  }
}

/*---------------------------------------------------------------------------*/
/* Equivalent à un var.synchronize() où var est une variable globale         */ 
/* (i.e. non multi-mat)                                                      */
/*---------------------------------------------------------------------------*/
template<typename MeshVariableRefT>
void VarSyncMng::globalSynchronize(MeshVariableRefT var)
{
  using ItemType = typename MeshVariableRefT::ItemType;
  using DataType = typename MeshVariableRefT::DataType;
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

  // Pour un ItemType donné, combien de DataType sont utilisés ? => degree
  Integer degree = get_var_degree(var);

  m_sync_buffers->resetBuf();
  // On prévoit une taille max du buffer qui va contenir tous les messages
  m_sync_buffers->addEstimatedMaxSz<DataType>(owned_item_idx_pn, degree);
  m_sync_buffers->addEstimatedMaxSz<DataType>(ghost_item_idx_pn, degree);
  // Le buffer de tous les messages est réalloué si pas assez de place
  m_sync_buffers->allocIfNeeded();

  // On récupère les adresses et tailles des buffers d'envoi et de réception
  auto buf_snd = m_sync_buffers->multiBufView<DataType>(owned_item_idx_pn, degree);
  auto buf_rcv = m_sync_buffers->multiBufView<DataType>(ghost_item_idx_pn, degree);

  // L'échange proprement dit des valeurs de var
  UniqueArray<Parallel::Request> requests;

  // On amorce les réceptions
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce la réception
    auto byte_buf_rcv = buf_rcv.byteBuf(inei); // le buffer de réception pour inei
    auto req_rcv = m_pm->recv(byte_buf_rcv, rank_nei, /*blocking=*/false);
    requests.add(req_rcv);
  }

  // On remplit les buffers sur CPU, TODO : sur GPU
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {

    // On lit les valeurs de var pour les recopier dans le buffer d'envoi
    auto buf_snd_inei = buf_snd.byteBuf(inei); // buffer dans lequel on va écrire
    // "buf_snd[inei] <= var"
    cpy_var2buf(owned_item_idx_pn[inei], var, buf_snd_inei);
  }

  // On amorce les envois
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce l'envoi
    auto byte_buf_snd = buf_snd.byteBuf(inei); // le buffer d'envoi pour inei
    auto req_snd = m_pm->send(byte_buf_snd, rank_nei, /*blocking=*/false);
    requests.add(req_snd);
  }

  m_pm->waitAllRequests(requests);
  requests.clear();

  // On recopie les valeurs reçues dans les buffers dans var
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    auto buf_rcv_inei = buf_rcv.byteBuf(inei); // buffer duquel on va lire les données
    // "var <= buf_rcv[inei]"
    cpy_buf2var(ghost_item_idx_pn[inei], buf_rcv_inei, var);
  }
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/
#include <arcane/utils/Real3.h>

#define INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE(__MeshVariableRefT__) \
  template void VarSyncMng::globalSynchronize(__MeshVariableRefT__ var)

INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE(VariableCellReal);
INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE(VariableNodeReal3);

INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE(VariableCellArrayReal3);
