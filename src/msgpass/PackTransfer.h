#ifndef MSG_PASS_PACK_TRANSFER_H
#define MSG_PASS_PACK_TRANSFER_H

#include "accenv/AcceleratorUtils.h"
#include "msgpass/SyncBuffers.h"

#include <arcane/accelerator/IRunQueueStream.h>

using namespace Arcane;

void async_transfer(Span<Byte> dst_buf, eLocMem dst_loc_mem,
    Span<const Byte> src_buf, eLocMem src_loc_mem, RunQueue& queue);

void async_transfer(MultiBufView out_buf, MultiBufView in_buf, RunQueue& queue);

/*---------------------------------------------------------------------------*/
/* async_pack_var2buf */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
void async_pack_var2buf(IntegerConstArrayView item_idx, 
    const MeshVarRefT<ItemType, DataType>& var,
    ArrayView<Byte> buf, RunQueue& queue) 
{
  // Ne devrait jamais être appelé
  ARCANE_ASSERT(false, ("async_pack_var2buf à spécifialiser"));
}

// Spécialisation pour MeshVariable***Scalar***RefT
template<typename ItemType, typename DataType>
void async_pack_var2buf(IntegerConstArrayView item_idx, 
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
/* async_unpack_buf2var */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
void async_unpack_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVarRefT<ItemType, DataType> &var, RunQueue& queue) 
{
  // Ne devrait jamais être appelé
  ARCANE_ASSERT(false, ("async_unpack_buf2var à spécifialiser"));
}

// Spécialisation pour MeshVariable***Scalar***RefT
template<typename ItemType, typename DataType>
void async_unpack_buf2var(IntegerConstArrayView item_idx, 
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

#endif

