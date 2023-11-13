#ifndef PATTERN_4_GPU_UTILS_H
#define PATTERN_4_GPU_UTILS_H

#include <mutex>

#ifdef ARCCORE_DEVICE_CODE
#include <cuda/atomic>
#endif
using namespace Arcane;


class Pattern4GPUUtils
{
 public:
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE Pattern4GPUUtils():
  min_mutex()
{
}

ARCCORE_HOST_DEVICE Pattern4GPUUtils(const Pattern4GPUUtils& rhs):
  min_mutex()
{
}

ARCCORE_HOST_DEVICE ~Pattern4GPUUtils(){};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE void P4GPU_atomicMin(int*old, int val) const{
#ifdef ARCCORE_DEVICE_CODE
   atomicMin(old,val);
#else
   min_mutex.lock();
   *old = (val < *old)?val:*old;
   min_mutex.unlock();
#endif
}

ARCCORE_HOST_DEVICE void P4GPU_atomicAdd(Real3*old, Real3 val) const{
#ifdef ARCCORE_DEVICE_CODE
   atomicAdd(&(old->x),val.x);
   atomicAdd(&(old->y),val.y);
   atomicAdd(&(old->z),val.z);
#else
   min_mutex.lock();
   old->x += val.x;
   old->y += val.y;
   old->z += val.z;
   min_mutex.unlock();
#endif
}

 private:
   mutable std::mutex min_mutex;
};

#endif
