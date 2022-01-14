#include "msgpass/VarSyncMng.h"

#include <arcane/IVariableSynchronizer.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/ISubDomain.h>
#include <arcane/IParallelMng.h>
#include <arcane/MeshVariableScalarRef.h>
#include <arcane/MeshVariableArrayRef.h>
#include <arcane/accelerator/IRunQueueStream.h>

#include <thread>
#include <mpi.h>

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
SyncItems<ItemType>::SyncItems(IMesh* mesh, Int32ConstArrayView neigh_ranks) :
  m_buf_owned_item_idx    (platform::getAcceleratorHostMemoryAllocator()),
  m_indexes_owned_item_pn (platform::getAcceleratorHostMemoryAllocator()),
  m_nb_owned_item_pn      (platform::getAcceleratorHostMemoryAllocator()),
  m_buf_ghost_item_idx    (platform::getAcceleratorHostMemoryAllocator()),
  m_indexes_ghost_item_pn (platform::getAcceleratorHostMemoryAllocator()),
  m_nb_ghost_item_pn      (platform::getAcceleratorHostMemoryAllocator())
{
  eItemKind item_kind = get_item_kind<ItemType>();
  IItemFamily* item_family = mesh->itemFamily(item_kind);
  IVariableSynchronizer* var_sync = item_family->allItemsSynchronizer();
  
  Integer nb_nei = neigh_ranks.size();

  // "shared" ou "owned" : les items intérieurs au sous-domaine et qui doivent être envoyés
  // "ghost" : les items fantômes pour lesquels on va recevoir des informations
  m_indexes_owned_item_pn.resize(nb_nei);
  m_indexes_ghost_item_pn.resize(nb_nei);
  m_nb_owned_item_pn.resize(nb_nei);
  m_nb_ghost_item_pn.resize(nb_nei);

  Integer accu_nb_owned=0;
  Integer accu_nb_ghost=0;
  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    m_nb_owned_item_pn[inei] = var_sync->sharedItems(inei).size();
    m_nb_ghost_item_pn[inei] = var_sync->ghostItems(inei).size();

    m_indexes_owned_item_pn[inei] = accu_nb_owned;
    m_indexes_ghost_item_pn[inei] = accu_nb_ghost;
    
    accu_nb_owned += m_nb_owned_item_pn[inei];
    accu_nb_ghost += m_nb_ghost_item_pn[inei];
  }
  m_buf_owned_item_idx.resize(accu_nb_owned);
  m_buf_ghost_item_idx.resize(accu_nb_ghost);

  // On construit des multi-vues sur des zones allouées en mémoire managées
  MultiArray2View<Integer> owned_item_idx_pn(m_buf_owned_item_idx.view(),
      m_indexes_owned_item_pn.constView(), m_nb_owned_item_pn.constView());

  MultiArray2View<Integer> ghost_item_idx_pn(m_buf_ghost_item_idx.view(),
      m_indexes_ghost_item_pn.constView(), m_nb_ghost_item_pn.constView());

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

  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    lids2itemidx(var_sync->sharedItems(inei), owned_item_idx_pn[inei]);
    lids2itemidx(var_sync->ghostItems(inei) , ghost_item_idx_pn[inei]);
  }
}

/*---------------------------------------------------------------------------*/
/* Encapsule des vues sur plusieurs buffers de communication                 */
/*---------------------------------------------------------------------------*/
MultiBufView::MultiBufView()
{ }

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
MultiBufView::MultiBufView(ArrayView<Byte*> ptrs, Int64ConstArrayView sizes,
    eLocMem loc_mem) :
  m_ptrs    (ptrs),
  m_sizes   (sizes),
  m_loc_mem (loc_mem)
{
  ARCANE_ASSERT(ptrs.size()==sizes.size(), ("ptrs.size()!=sizes.size()"));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
MultiBufView::MultiBufView(const MultiBufView& rhs) :
  m_ptrs    (rhs.m_ptrs),
  m_sizes   (rhs.m_sizes),
  m_loc_mem (rhs.m_loc_mem)
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
/* SyncBuffers                                                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
SyncBuffers::SyncBuffers(bool is_acc_avl) :
  m_is_accelerator_available (is_acc_avl)
{
  for(Integer imem(0) ; imem<2 ; ++imem) {
    m_buf_mem[imem].m_ptr=nullptr;
    m_buf_mem[imem].m_size=0;
    m_buf_mem[imem].m_first_av_pos=0;
    m_buf_mem[imem].m_loc_mem=LM_HostMem;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
SyncBuffers::~SyncBuffers() {
  for(Integer imem(0) ; imem<2 ; ++imem) {
    if (m_buf_mem[imem].m_loc_mem == LM_HostMem) {
#ifdef ARCANE_COMPILING_CUDA
      if (m_is_accelerator_available) {
        cudaFreeHost(reinterpret_cast<void*>(m_buf_mem[imem].m_ptr));
      } else {
        delete[] m_buf_mem[imem].m_ptr;
      }
#else
      delete[] m_buf_mem[imem].m_ptr;
#endif
    }
#ifdef ARCANE_COMPILING_CUDA
    else if (m_buf_mem[imem].m_ptr)
    {
      ARCANE_ASSERT(m_buf_mem[imem].m_loc_mem == LM_DevMem, 
          ("Impossible de libérer de la mémoire qui n'est pas sur le device"));
      cudaFree(reinterpret_cast<void*>(m_buf_mem[imem].m_ptr));
    }
#endif
  }
}

/* A partir des items à communiquer, estime une borne sup de la taille du    */ 
/* buffer en octets                                                          */
/*---------------------------------------------------------------------------*/
template<typename DataType>
Int64 SyncBuffers::estimatedMaxBufSz(IntegerConstArrayView item_sizes, 
    Integer degree) {
  // HYPOTHESE 1 : même valeur de sizeof(DataType) sur CPU et GPU
  // HYPOTHESE 2 : même valeur de alignof(DataType) sur CPU et GPU
  // TODO : comment le vérifier ?
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
  for(Integer imem(0) ; imem<2 ; ++imem) {
    m_buf_mem[imem].m_first_av_pos=0;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
template<typename DataType>
void SyncBuffers::addEstimatedMaxSz(ConstMultiArray2View<Integer> item_idx_pn,
    Integer degree) {
  m_buf_estim_sz += estimatedMaxBufSz<DataType>(item_idx_pn.dim2Sizes(), degree);
}

/*---------------------------------------------------------------------------*/
/* Reallocation dans la mémoire hôte */
/*---------------------------------------------------------------------------*/
void SyncBuffers::BufMem::reallocIfNeededOnHost(Int64 wanted_size, bool is_acc_avl) {
  m_loc_mem = LM_HostMem;
  // S'il n'y a pas assez d'espace, on réalloue (peu importe les données précédentes)
  if (m_size < wanted_size) {
#ifdef ARCANE_COMPILING_CUDA
    if (is_acc_avl) {
      cudaFreeHost(reinterpret_cast<void*>(m_ptr));
      void* h_ptr;
      cudaMallocHost(&h_ptr, wanted_size);
      m_ptr = reinterpret_cast<Byte*>(h_ptr);
    } else {
      delete[] m_ptr;
      m_ptr = new Byte[wanted_size];
    }
#else
    delete[] m_ptr;
    m_ptr = new Byte[wanted_size];
#endif
    m_size = wanted_size;
  }
}

/*---------------------------------------------------------------------------*/
/* Reallocation dans la mémoire device */
/*---------------------------------------------------------------------------*/
#ifdef ARCANE_COMPILING_CUDA
void SyncBuffers::BufMem::reallocIfNeededOnDevice(Int64 wanted_size) {
  m_loc_mem = LM_DevMem;
  // S'il n'y a pas assez d'espace, on réalloue (peu importe les données précédentes)
  if (m_size < wanted_size) {
    cudaFree(reinterpret_cast<void*>(m_ptr));
    void* d_ptr;
    cudaMalloc(&d_ptr, wanted_size);
    m_ptr = reinterpret_cast<Byte*>(d_ptr);
    m_size = wanted_size;
  }
}
#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void SyncBuffers::allocIfNeeded() {
  // D'abord l'hote
  m_buf_mem[0].reallocIfNeededOnHost(m_buf_estim_sz, m_is_accelerator_available);

  // Puis le device si celui-ci existe
  if (m_is_accelerator_available) {
#ifdef ARCANE_COMPILING_CUDA
    m_buf_mem[1].reallocIfNeededOnDevice(m_buf_estim_sz);
#endif
  }
  if (!m_is_accelerator_available) {
    // Pour débugger, le buffer "device" se trouve dans la mémoire hôte
    m_buf_mem[1].reallocIfNeededOnHost(m_buf_estim_sz, m_is_accelerator_available);
  }
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
    ConstMultiArray2View<Integer> item_idx_pn, Integer degree, Integer imem) {

  auto& buf_mem = m_buf_mem[imem];
  Byte* new_ptr = buf_mem.m_ptr+buf_mem.m_first_av_pos;
  Int64 av_space = buf_mem.m_size-buf_mem.m_first_av_pos;
  Span<Byte> buf_bytes(new_ptr, av_space);

  auto mb = _multiBufView<DataType>(item_idx_pn.dim2Sizes(), degree, buf_bytes);
  mb.locMem() = buf_mem.m_loc_mem;

  auto rg{mb.rangeSpan()}; // Encapsule [beg_ptr, end_ptr[
  Byte* end_ptr = rg.data()+rg.size();
  buf_mem.m_first_av_pos = (end_ptr - buf_mem.m_ptr);
  return mb;
}

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
    const MeshVarRefT<ItemType, DataType>& var,
    ArrayView<Byte> buf) 
{
  // Ne devrait jamais être appelé
  ARCANE_ASSERT(false, ("cpy_var2buf à spécifialiser"));
}

// Spécialisation pour MeshVariable***Scalar***RefT
template<typename ItemType, typename DataType>
void cpy_var2buf(IntegerConstArrayView item_idx, 
    const MeshVariableScalarRefT<ItemType, DataType>& var,
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
    const MeshVariableArrayRefT<ItemType, DataType>& var,
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
    MeshVarRefT<ItemType, DataType> &var) 
{
  // Ne devrait jamais être appelé
  ARCANE_ASSERT(false, ("cpy_var2buf à spécifialiser"));
}

// Spécialisation pour MeshVariable***Scalar***RefT
template<typename ItemType, typename DataType>
void cpy_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVariableScalarRefT<ItemType, DataType> &var)
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
    MeshVariableArrayRefT<ItemType, DataType> &var)
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
  auto buf_snd = m_sync_buffers->multiBufView<DataType>(owned_item_idx_pn, degree, 0);
  auto buf_rcv = m_sync_buffers->multiBufView<DataType>(ghost_item_idx_pn, degree, 0);

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
/* async_cpy */
/*---------------------------------------------------------------------------*/
void async_cpy(Span<Byte> dst_buf, eLocMem dst_loc_mem,
    Span<const Byte> src_buf, eLocMem src_loc_mem, RunQueue& queue) {
  ARCANE_ASSERT(!(src_loc_mem==LM_DevMem && dst_loc_mem==LM_DevMem), ("Copie de mémoire device à mémoire device non supportée"));

  ARCANE_ASSERT(src_buf.size()==dst_buf.size(), ("Les buffers src et dst n'ont pas la meme taille"));

  if (src_loc_mem==LM_HostMem && dst_loc_mem==LM_HostMem)
  {
    std::memcpy(dst_buf.data(), src_buf.data(), src_buf.size());
  }
  else
  {
#ifdef ARCANE_COMPILING_CUDA
    auto* rq = queue._internalStream();
    cudaStream_t* s = reinterpret_cast<cudaStream_t*>(rq->_internalImpl());

    cudaMemcpyKind mem_kind = (src_loc_mem==LM_DevMem && dst_loc_mem==LM_HostMem ?
        cudaMemcpyDeviceToHost : cudaMemcpyHostToDevice);
#if 0
    static const char* str_mem_kind[5] = {"cudaMemcpyHostToHost", 
      "cudaMemcpyHostToDevice", "cudaMemcpyDeviceToHost", 
      "cudaMemcpyDeviceToDevice", "cudaMemcpyDefault"};
    std::cout << str_mem_kind[mem_kind] << " : " << src_buf.size() << " bytes" << std::endl;
#endif
    cudaMemcpyAsync(dst_buf.data(), src_buf.data(), src_buf.size(), mem_kind, *s);
#else
    ARCANE_ASSERT(false, ("src_buf ou dst_buf est déclarée en mémoire device mais support CUDA non disponible"));
#endif
  }
}

/*---------------------------------------------------------------------------*/
/* async_cpy */
/*---------------------------------------------------------------------------*/
void async_cpy(MultiBufView out_buf, MultiBufView in_buf, RunQueue& queue) {
  async_cpy(out_buf.rangeSpan(), out_buf.locMem(), 
      in_buf.rangeSpan(), in_buf.locMem(), queue);
}

/*---------------------------------------------------------------------------*/
/* async_cpy_var2buf */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
void async_cpy_var2buf(IntegerConstArrayView item_idx, 
    const MeshVarRefT<ItemType, DataType>& var,
    ArrayView<Byte> buf, RunQueue& queue) 
{
  // Ne devrait jamais être appelé
  ARCANE_ASSERT(false, ("async_cpy_var2buf à spécifialiser"));
}

// Spécialisation pour MeshVariable***Scalar***RefT
template<typename ItemType, typename DataType>
void async_cpy_var2buf(IntegerConstArrayView item_idx, 
    const MeshVariableScalarRefT<ItemType, DataType>& var,
    ArrayView<Byte> buf, RunQueue& queue) 
{
  using ItemIdType = typename ItemType::LocalIdType;

  auto command = makeCommand(queue);

  auto in_var = ax::viewIn(command, var);

  Span<const Integer> in_item_idx(item_idx);
  Span<DataType> buf_vals(MultiBufView::valBuf<DataType>(buf));

  Integer nb_item_idx = item_idx.size();

  command << RUNCOMMAND_LOOP1(iter, nb_item_idx) {
    auto [i] = iter();
    ItemIdType lid(in_item_idx[i]); 
    buf_vals[i] = in_var[lid];
  }; // asynchrone
}

/*---------------------------------------------------------------------------*/
/* async_cpy_buf2var */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
void async_cpy_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVarRefT<ItemType, DataType> &var, RunQueue& queue) 
{
  // Ne devrait jamais être appelé
  ARCANE_ASSERT(false, ("async_cpy_buf2var à spécifialiser"));
}

// Spécialisation pour MeshVariable***Scalar***RefT
template<typename ItemType, typename DataType>
void async_cpy_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVariableScalarRefT<ItemType, DataType> &var, RunQueue& queue)
{
  using ItemIdType = typename ItemType::LocalIdType;

  auto command = makeCommand(queue);

  Span<const Integer>  in_item_idx(item_idx);
  Span<const DataType> buf_vals(MultiBufView::valBuf<DataType>(buf));

  auto out_var = ax::viewOut(command, var);

  Integer nb_item_idx = item_idx.size();

  command << RUNCOMMAND_LOOP1(iter, nb_item_idx) {
    auto [i] = iter();
    ItemIdType lid(in_item_idx[i]);
    out_var[lid] = buf_vals[i];
  }; // asynchrone
}

/*---------------------------------------------------------------------------*/
/* Equivalent à un var.synchronize() où var est une variable globale         */ 
/* (i.e. non multi-mat) dont les données sont présentes sur GPU              */
/*---------------------------------------------------------------------------*/
template<typename MeshVariableRefT>
void VarSyncMng::globalSynchronizeDev(MeshVariableRefT var)
{
  if (m_nb_nei==0) {
    return;
  }

  using ItemType = typename MeshVariableRefT::ItemType;
  using DataType = typename MeshVariableRefT::DataType;

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
  // sur l'HOTE (_h et LM_HostMem)
  auto buf_snd_h = m_sync_buffers->multiBufView<DataType>(owned_item_idx_pn, degree, 0);
  auto buf_rcv_h = m_sync_buffers->multiBufView<DataType>(ghost_item_idx_pn, degree, 0);

  // On récupère les adresses et tailles des buffers d'envoi et de réception 
  // sur le DEVICE (_d et LM_DevMem)
  auto buf_snd_d = m_sync_buffers->multiBufView<DataType>(owned_item_idx_pn, degree, 1);
  auto buf_rcv_d = m_sync_buffers->multiBufView<DataType>(ghost_item_idx_pn, degree, 1);

  // L'échange proprement dit des valeurs de var
  UniqueArray<Parallel::Request> requests;

  // On amorce les réceptions
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce la réception sur l'HOTE
    auto byte_buf_rcv = buf_rcv_h.byteBuf(inei); // le buffer de réception pour inei
    auto req_rcv = m_pm->recv(byte_buf_rcv, rank_nei, /*blocking=*/false);
    requests.add(req_rcv);
  }

  auto queue = makeQueue(m_runner);
  queue.setAsync(true);

  // On enchaine sur le device : 
  //    copie de var_dev dans buf_dev 
  //    puis transfert buf_dev => buf_hst

  // On remplit les buffers sur le DEVICE
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {

    // On lit les valeurs de var pour les recopier dans le buffer d'envoi
    auto buf_snd_inei = buf_snd_d.byteBuf(inei); // buffer dans lequel on va écrire
    // "buf_snd[inei] <= var"
    async_cpy_var2buf(owned_item_idx_pn[inei], var, buf_snd_inei, queue);
  }

  // transfert buf_snd_d => buf_snd_h
  async_cpy(buf_snd_h, buf_snd_d, queue);
  queue.barrier(); // attendre que la copie sur l'hôte soit terminée pour envoyer les messages

  // On amorce les envois sur l'HOTE
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce l'envoi
    auto byte_buf_snd = buf_snd_h.byteBuf(inei); // le buffer d'envoi pour inei
    auto req_snd = m_pm->send(byte_buf_snd, rank_nei, /*blocking=*/false);
    requests.add(req_snd);
  }

  m_pm->waitAllRequests(requests);
  requests.clear();

  // transfert buf_rcv_h => buf_rcv_d
  async_cpy(buf_rcv_d, buf_rcv_h, queue);

  // On recopie les valeurs reçues dans les buffers dans var sur le DEVICE
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    auto buf_rcv_inei = buf_rcv_d.byteBuf(inei); // buffer duquel on va lire les données
    // "var <= buf_rcv_d[inei]"
    async_cpy_buf2var(ghost_item_idx_pn[inei], buf_rcv_inei, var, queue);
  }
  queue.barrier(); // attendre que les copies soient terminées sur GPU
}

/*---------------------------------------------------------------------------*/
/* Equivalent à un var.synchronize() où var est une variable globale         */ 
/* (i.e. non multi-mat) en utilisant des comms bloquantes dans des taches    */
/*---------------------------------------------------------------------------*/
template<typename MeshVariableRefT>
void VarSyncMng::globalSynchronizeDevThr(MeshVariableRefT var)
{
  if (m_nb_nei==0) {
    return;
  }

  using ItemType = typename MeshVariableRefT::ItemType;
  using DataType = typename MeshVariableRefT::DataType;

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
  // sur l'HOTE (_h et LM_HostMem)
  auto buf_snd_h = m_sync_buffers->multiBufView<DataType>(owned_item_idx_pn, degree, 0);
  auto buf_rcv_h = m_sync_buffers->multiBufView<DataType>(ghost_item_idx_pn, degree, 0);

  // On récupère les adresses et tailles des buffers d'envoi et de réception 
  // sur le DEVICE (_d et LM_DevMem)
  auto buf_snd_d = m_sync_buffers->multiBufView<DataType>(owned_item_idx_pn, degree, 1);
  auto buf_rcv_d = m_sync_buffers->multiBufView<DataType>(ghost_item_idx_pn, degree, 1);

#define USE_MPI_REQUEST

#ifdef USE_MPI_REQUEST
//#warning "USE_MPI_REQUEST"
  using RequestType = MPI_Request;
#else
  using RequestType = Parallel::Request;
#endif

  // L'échange proprement dit des valeurs de var
  Integer tag=1000;
  UniqueArray<RequestType> requests(2*m_nb_nei);
  IntegerUniqueArray msg_types(2*m_nb_nei); // nature des messages 

  // On amorce les réceptions
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce la réception
    auto byte_buf_rcv_h = buf_rcv_h.byteBuf(inei); // le buffer de réception pour inei
#ifdef USE_MPI_REQUEST
    MPI_Irecv(byte_buf_rcv_h.data(), byte_buf_rcv_h.size(), MPI_BYTE, rank_nei, tag, 
        MPI_COMM_WORLD, &(requests[inei]));
#else
    requests[inei] = m_pm->recv(byte_buf_rcv_h, rank_nei, /*blocking=*/false);
#endif
    msg_types[inei] = inei+1; // >0 pour la réception
  }

  // La tâche à effectuer pour un voisin
  auto lbd_sender = [&](Integer inei, RunQueue& queue) {

    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    auto byte_buf_snd_d = buf_snd_d.byteBuf(inei); // le buffer d'envoi pour inei sur le DEVICE
    auto byte_buf_snd_h = buf_snd_h.byteBuf(inei); // le buffer d'envoi pour inei sur l'HOTE

    // "byte_buf_snd_d <= var"
    async_cpy_var2buf(owned_item_idx_pn[inei], var, byte_buf_snd_d, queue);

    // transfert buf_snd_d[inei] => buf_snd_h[inei]
    async_cpy(byte_buf_snd_h, buf_snd_h.locMem(),
        byte_buf_snd_d, buf_snd_d.locMem(), queue);

    // attendre que la copie sur l'hôte soit terminée pour envoyer les messages
    queue.barrier(); 

    // On amorce les envois
#ifdef USE_MPI_REQUEST
    MPI_Isend(byte_buf_snd_h.data(), byte_buf_snd_h.size(), MPI_BYTE, rank_nei, tag,
        MPI_COMM_WORLD, &(requests[m_nb_nei+inei]));
#else
    requests[m_nb_nei+inei] = m_pm->send(byte_buf_snd_h, rank_nei, /*blocking=*/false);
#endif
    msg_types[m_nb_nei+inei] = -inei-1; // <0 pour l'envoi
  };

//#define USE_THR_SENDER
#ifdef USE_THR_SENDER
#warning "USE_THR_SENDER"
  UniqueArray<std::thread*> thr_sender(m_nb_nei);

  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    thr_sender[inei] = 
      new std::thread(lbd_sender, inei, std::ref(m_neigh_queues->queue(inei)));
  }

  // On attend la fin de tous les threads
  for(auto thr : thr_sender) {
    thr->join();
    delete thr;
  }
#else
  // On lance en SEQUENTIEL
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    lbd_sender(inei, m_neigh_queues->queue(inei));
  }
#endif

  // Tache qui unpack les données du buffer reçues par un voisin inei
  //    copie de var_dev dans buf_dev 
  //    puis transfert buf_dev => buf_hst
  auto lbd_unpacker = [&](Integer inei, RunQueue& queue) {
    auto byte_buf_rcv_h = buf_rcv_h.byteBuf(inei); // buffer des données reçues sur l'HOTE
    auto byte_buf_rcv_d = buf_rcv_d.byteBuf(inei); // buffer des données reçues à transférer sur le DEVICE

    // transfert buf_rcv_h[inei] => buf_rcv_d[inei]
    async_cpy(byte_buf_rcv_d, buf_rcv_d.locMem(),
        byte_buf_rcv_h, buf_rcv_h.locMem(), queue);

    // "var <= buf_rcv_d[inei]"
    async_cpy_buf2var(ghost_item_idx_pn[inei], byte_buf_rcv_d, var, queue);

    // attendre que les copies soient terminées sur GPU
    queue.barrier(); 
  };
  UniqueArray<std::thread*> thr_unpacker;

  ARCANE_ASSERT(2*m_nb_nei==requests.size(), 
      ("Le nb de requetes n'est pas egal à 2 fois le nb de voisins"));

  UniqueArray<RequestType> requests2(2*m_nb_nei);
  IntegerUniqueArray msg_types2(2*m_nb_nei); 
  UniqueArray<bool> is_done_req(2*m_nb_nei);
#ifdef USE_MPI_REQUEST
  IntegerUniqueArray array_of_indices(2*m_nb_nei);
#endif

  // On utilise des vues pour éviter de réallouer en permanence des tableaux
  ArrayView<RequestType> pending_requests(requests.view());
  ArrayView<Integer> pending_types(msg_types.view());

  ArrayView<RequestType> upd_pending_requests(requests2.view());
  ArrayView<Integer> upd_pending_types(msg_types2.view());

  ArrayView<RequestType> tmp_pending_requests;
  ArrayView<Integer> tmp_pending_types;

  Integer nb_iter_wait_some = 0;
  Integer nb_pending_rcv = m_nb_nei;

  while(nb_pending_rcv>0) {

    Integer nb_pending_req = pending_requests.size();

    // On dimenensionne is_done_requests au nb de requêtes d'avant waitSomeRequests
    // et on initialise à false
    ArrayView<bool> is_done_requests(is_done_req.subView(0, nb_pending_req));
    for(Integer ireq=0 ; ireq<nb_pending_req ; ++ireq) {
      is_done_requests[ireq]=false;
    }

    // Attente de quelques requetes
#ifdef USE_MPI_REQUEST
    Integer nb_req_done=0;
    MPI_Waitsome(pending_requests.size(), pending_requests.data(),
        &nb_req_done, array_of_indices.data(), MPI_STATUSES_IGNORE);
    IntegerArrayView done_indexes(array_of_indices.subView(0, nb_req_done));
#else
    IntegerUniqueArray done_indexes = m_pm->waitSomeRequests(pending_requests);
#endif

    for(Integer idone_req : done_indexes) {
      if (pending_types[idone_req] > 0) { // >0 signifie que c'est une requête de reception

        nb_pending_rcv--; // on une requete de reception en moins

        // On récupère l'indice du voisin
        Integer inei = pending_types[idone_req]-1;
        ARCANE_ASSERT(inei>=0 && inei<m_nb_nei, ("Mauvais indice de voisin"));

        // Maintenant qu'on a reçu le buffer pour le inei-ième voisin, 
        // on unpack les donnees dans un thread
//#define USE_THR_UNPACKER
#ifdef USE_THR_UNPACKER
#warning "USE_THR_UNPACKER"
        thr_unpacker.add(
            new std::thread(lbd_unpacker, inei, std::ref(m_neigh_queues->queue(inei)))
            );
#else
        lbd_unpacker(inei, m_neigh_queues->queue(inei));
#endif
      }
      is_done_requests[idone_req] = true;
    }
    // Il faut créer le nouveau tableau de requêtes pending dans upd_*
    Integer upd_nb_pending_req=0;
    for(Integer ireq=0 ; ireq<nb_pending_req ; ++ireq) {
      if (!is_done_requests[ireq]) {
        upd_pending_requests[upd_nb_pending_req]=pending_requests[ireq];
        upd_pending_types   [upd_nb_pending_req]=pending_types[ireq];
        upd_nb_pending_req++;
      }
    }

    // On échange les vues pour qu'à l'itération suivante 
    // pending_requests pointe vers upd_pending_types
    tmp_pending_requests = upd_pending_requests.subView(0, upd_nb_pending_req);
    upd_pending_requests = pending_requests;
    pending_requests = tmp_pending_requests;

    tmp_pending_types = upd_pending_types.subView(0, upd_nb_pending_req);
    upd_pending_types = pending_types;
    pending_types = tmp_pending_types;

    nb_iter_wait_some++;
#if 0
    std::ostringstream ostr;
    ostr << "P=" << m_pm->commRank() 
      << ", iter_wait_some=" << nb_iter_wait_some
      << ", nb_done=" << done_indexes.size();
    std::cout << ostr.str() << std::endl;
#endif
  }

  // Ici, toutes les requetes de receptions sont forcement terminées 
  // (condition de la boucle while précédente)
  // Mais il peut rester encore des requetes d'envoi en cours
  if (pending_requests.size()) {
    // Normalement, il ne reste que des requêtes d'envois
    ARCANE_ASSERT(pending_requests.size()<=m_nb_nei, 
        ("Il ne peut pas rester un nb de requetes d'envoi supérieur au nb de voisins"));
    for(Integer msg_type : pending_types) {
      ARCANE_ASSERT(msg_type<0, 
          ("Un message d'envoi doit avoir un type négatif ce qui n'est pas le cas"));
    }
#if 0
    std::ostringstream ostr;
    ostr << "P=" << m_pm->commRank() 
      << ", WaitAll pending_requests.size()=" << pending_requests.size();
    std::cout << ostr.str() << std::endl;
#endif
#ifdef USE_MPI_REQUEST
    MPI_Waitall(pending_requests.size(), pending_requests.data(), MPI_STATUSES_IGNORE);
#else
    m_pm->waitAllRequests(pending_requests);
#endif
  }

#ifdef USE_THR_UNPACKER
#warning "USE_THR_UNPACKER : join"
  // On attend la fin de tous les threads unpackers
  for(auto thr : thr_unpacker) {
    thr->join();
    delete thr;
  }
#endif
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

#define INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_DEV(__MeshVariableRefT__) \
  template void VarSyncMng::globalSynchronizeDev(__MeshVariableRefT__ var)

INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_DEV(VariableCellReal);
INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_DEV(VariableNodeReal3);

#define INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_DEV_THR(__MeshVariableRefT__) \
  template void VarSyncMng::globalSynchronizeDevThr(__MeshVariableRefT__ var)

INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_DEV_THR(VariableCellReal);
INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_DEV_THR(VariableNodeReal3);
