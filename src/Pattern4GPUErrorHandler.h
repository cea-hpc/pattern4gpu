#ifndef PATTERN_4_GPU_ERROR_HANDLER_H
#define PATTERN_4_GPU_ERROR_HANDLER_H

#include "Pattern4GPUUtils.h"

using namespace Arcane;


class Pattern4GPUErrorHandler
{
 public:
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Pattern4GPUErrorHandler():
  min_cid_on_err(MemoryAllocationOptions(platform::getAcceleratorHostMemoryAllocator(),eMemoryLocationHint::MainlyHost),1)
{
  min_cid_on_err[0] = std::numeric_limits<Int32>::max();
}

Pattern4GPUErrorHandler(const Pattern4GPUErrorHandler& rhs):
 min_cid_on_err(rhs.min_cid_on_err),
  utils(rhs.utils)
{
}

~Pattern4GPUErrorHandler(){
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
 int* Get_min_ptr(){
   return min_cid_on_err.data();
 }

Pattern4GPUUtils*Get_utils_ptr(){
  return &utils;
}

int P4GPUGetError() const{
  return P4GPUIsError() ?min_cid_on_err[0]: -1;
  }

bool P4GPUIsError() const{
  return min_cid_on_err[0] != std::numeric_limits<Int32>::max();
  }
 private:
  UniqueArray<Integer> min_cid_on_err;
  Pattern4GPUUtils utils;

 };

class Pattern4GPUErrorHandlerView
{
  public:
  Pattern4GPUErrorHandlerView(Pattern4GPUErrorHandler& hand)
  {
    out_min_cid_on_err = hand.Get_min_ptr();
    utils_ptr = hand.Get_utils_ptr();
  }
  
  ARCCORE_HOST_DEVICE Pattern4GPUErrorHandlerView(const Pattern4GPUErrorHandlerView& hand):
    out_min_cid_on_err(hand.out_min_cid_on_err),
    utils_ptr(hand.utils_ptr)
  {
  }
  ARCCORE_HOST_DEVICE ~Pattern4GPUErrorHandlerView()
  {
  }
  
  ARCCORE_HOST_DEVICE void P4GPUThrowError(int cid) const{
  utils_ptr->P4GPU_atomicMin(out_min_cid_on_err,cid);
  }


  private: 
  mutable int*out_min_cid_on_err;
  Pattern4GPUUtils*utils_ptr;

};

#endif
