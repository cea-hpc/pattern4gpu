#include "msgpass/SyncBuffers.h"

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
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/

#define INST_SYNC_BUFFERS(__DataType__) \
  template ArrayView<__DataType__> MultiBufView::valBuf<__DataType__>(ArrayView<Byte> buf); \
  template Array2View<__DataType__> MultiBufView::valBuf2<__DataType__>(ArrayView<Byte> buf, Integer dim2_size); \
  template MultiBufView SyncBuffers::multiBufView<__DataType__>(ConstMultiArray2View<Integer> item_idx_pn, Integer degree, Integer imem); \
  template void SyncBuffers::addEstimatedMaxSz<__DataType__>(ConstMultiArray2View<Integer> item_idx_pn, Integer degree)

INST_SYNC_BUFFERS(Real);
INST_SYNC_BUFFERS(Real3);

