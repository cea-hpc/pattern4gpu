#include "msgpass/PackTransfer.h"

/*---------------------------------------------------------------------------*/
/* async_transfer */
/*---------------------------------------------------------------------------*/
void async_transfer(Span<Byte> dst_buf, eLocMem dst_loc_mem,
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
/* async_transfer */
/*---------------------------------------------------------------------------*/
void async_transfer(MultiBufView out_buf, MultiBufView in_buf, RunQueue& queue) {
  async_transfer(out_buf.rangeSpan(), out_buf.locMem(), 
      in_buf.rangeSpan(), in_buf.locMem(), queue);
}

