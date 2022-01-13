#ifndef VAR_SYNC_MNG_H
#define VAR_SYNC_MNG_H

#include "accenv/AcceleratorUtils.h"
#include <arcane/IMesh.h>
#include <arcane/utils/MultiArray2.h>

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/* Encapsule la liste des items à envoyer/recevoir */
/*---------------------------------------------------------------------------*/
template<typename ItemType>
class SyncItems {
 public:
  SyncItems(IMesh* mesh, Int32ConstArrayView neigh_ranks);
  virtual ~SyncItems() {}

  auto nbOwnedItemIdxPn() const {
    return m_nb_owned_item_pn.constView();
  }

  auto nbGhostItemIdxPn() const {
    return m_nb_ghost_item_pn.constView();
  }

  auto ownedItemIdxPn() const {
    return ConstMultiArray2View<Integer>(m_buf_owned_item_idx.constView(),
        m_indexes_owned_item_pn.constView(), m_nb_owned_item_pn.constView());
  }

  auto ghostItemIdxPn() const {
    return ConstMultiArray2View<Integer>(m_buf_ghost_item_idx.constView(),
        m_indexes_ghost_item_pn.constView(), m_nb_ghost_item_pn.constView());
  }

 protected:
  // "shared" ou "owned" : les items intérieurs au sous-domaine et qui doivent être envoyés
  // "ghost" : les items fantômes pour lesquels on va recevoir des informations
  // _pn : _per_neigh, info par voisin
  
  // Nb d'items par voisin et listes des indexes de ces items
  // Owned ou Shared
  IntegerUniqueArray m_buf_owned_item_idx;
  IntegerUniqueArray m_indexes_owned_item_pn;
  IntegerUniqueArray m_nb_owned_item_pn;

  // Ghost
  IntegerUniqueArray m_buf_ghost_item_idx;
  IntegerUniqueArray m_indexes_ghost_item_pn;
  IntegerUniqueArray m_nb_ghost_item_pn;
};

/*---------------------------------------------------------------------------*/
/* Emplacement de la mémoire                                                 */
/*---------------------------------------------------------------------------*/
enum eLocMem {
  LM_HostMem=0,
  LM_DevMem,
  MAX_LocMem
};

/*---------------------------------------------------------------------------*/
/* Encapsule les buffers des valeurs à envoyer/recevoir                      */
/*---------------------------------------------------------------------------*/
class MultiBufView {
 public:
  MultiBufView();
  MultiBufView(ArrayView<Byte*> ptrs, Int64ConstArrayView sizes, eLocMem loc_mem=LM_HostMem);
  MultiBufView(const MultiBufView& rhs);

  MultiBufView& operator=(const MultiBufView& rhs);

  //! Convertit un buffer d'octets en buffer de DataType
  template<typename DataType>
  static ArrayView<DataType> valBuf(ArrayView<Byte> buf);

  //! Convertit un buffer d'octets en buffer 2D de DataType dont la taille dans la deuxième dimension est dim2_size
  template<typename DataType>
  static Array2View<DataType> valBuf2(ArrayView<Byte> buf, Integer dim2_size);

  //! Accès en lecture/écriture au i-ème buffer d'octets
  ArrayView<Byte> byteBuf(Integer i);

  //! Retourne [beg_ptr, end_ptr[ qui contient tous les buffers (peut-être espacés de trous)
  Span<Byte> rangeSpan();

  //! Emplacement de la mémoire
  eLocMem& locMem() { return m_loc_mem; }
  eLocMem locMem() const { return m_loc_mem; }

 protected:
  SharedArray<Byte*> m_ptrs;
  SharedArray<Int64> m_sizes;
  eLocMem m_loc_mem=LM_HostMem;  //! emplacement de la mémoire
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
class SyncBuffers {
 public:
  SyncBuffers(bool is_acc_avl);
  virtual ~SyncBuffers();

  //
  void resetBuf();

  // 
  template<typename DataType>
  void addEstimatedMaxSz(ConstMultiArray2View<Integer> item_idx_pn, Integer degree);

  // 
  void allocIfNeeded();

  /*!
   * \brief A partir des nb d'items à communiquer, estime une borne sup de la taille du buffer en octets
   */
  template<typename DataType>
  Int64 estimatedMaxBufSz(IntegerConstArrayView item_sizes, Integer degree);

  /*!
   * \brief A partir de la vue sur m_buf_bytes, construit une vue par voisin des buffers
   */
  template<typename DataType>
  MultiBufView multiBufView(
    ConstMultiArray2View<Integer> item_idx_pn, Integer degree, Integer imem);

 protected:
  /*!
   * \brief A partir de la vue sur un buffer déjà alloué, construit une vue par voisin des buffers
   */
  template<typename DataType>
  MultiBufView _multiBufView(
      IntegerConstArrayView item_sizes, Integer degree,
      Span<Byte> buf_bytes);

 protected:
  struct BufMem {
    Byte* m_ptr=nullptr;  //! Adresse de base du buffer
    Int64 m_size=0;  //! Taille totale (en octets) du buffer
    Int64 m_first_av_pos=0;  //! Première position disponible
    eLocMem m_loc_mem=LM_HostMem;  //! Emplacement de la mémoire
    void reallocIfNeededOnHost(Int64 wanted_size, bool is_acc_avl);
#ifdef ARCANE_COMPILING_CUDA
    void reallocIfNeededOnDevice(Int64 wanted_size);
#endif
  };
 protected:
  bool m_is_accelerator_available=false;  //! Vrai si un GPU est disponible pour les calculs
  Int64 m_buf_estim_sz=0;  //! Taille qui va servir à allouer
  // Pour gérer les buffers sur l'hote et le device
  BufMem m_buf_mem[2];
};

/*---------------------------------------------------------------------------*/
/* Gère les synchronisations des mailles fantômes par Message Passing        */
/*---------------------------------------------------------------------------*/
class VarSyncMng {
 public:
  VarSyncMng(IMesh* mesh, ax::Runner& runner);
  virtual ~VarSyncMng();

  //! Retourne vrai si un GPU est dispo pour exécuter les calculs
  bool isAcceleratorAvailable() const;

  // Equivalent à un var.synchronize() où var est une variable globale (i.e. non multi-mat)
  template<typename MeshVariableRefT>
  void globalSynchronize(MeshVariableRefT var);

  // Equivalent à un globalSynchronize pour lequel les données de var sont sur le DEVice
  template<typename MeshVariableRefT>
  void globalSynchronizeDev(MeshVariableRefT var);

  // Equivalent à un globalSynchronize pour lequel les données de var sont sur le DEVice
  template<typename MeshVariableRefT>
  void globalSynchronizeDevThr(MeshVariableRefT var);

 protected:

  // Retourne l'instance de SyncItems<T> en fonction de T
  template<typename ItemType>
  SyncItems<ItemType>* _getSyncItems();

 protected:

  ax::Runner& m_runner;

  IParallelMng* m_pm;  //! pour effectuer les send/receive proprement dit

  Integer m_nb_nei;  //! Nb de voisins (m_neigh_ranks.size())
  Int32ConstArrayView m_neigh_ranks;  //! Liste des rangs des voisins

  SyncItems<Cell>* m_sync_cells=nullptr;
  SyncItems<Node>* m_sync_nodes=nullptr;

  SyncBuffers* m_sync_buffers=nullptr;
};

#endif

