#ifndef ACCELERATOR_UTILS_H
#define ACCELERATOR_UTILS_H

#include "arcane/IApplication.h"
#include "arcane/accelerator/Reduce.h"
#include "arcane/accelerator/Runner.h"
#include "arcane/accelerator/Views.h"
#include "arcane/accelerator/Accelerator.h"
#include "arcane/accelerator/RunCommandLoop.h"
#include "arcane/accelerator/RunCommandEnumerate.h"

/*---------------------------------------------------------------------------*/
/* Pour les accélérateurs                                                    */
/*---------------------------------------------------------------------------*/

namespace ax = Arcane::Accelerator;

/*---------------------------------------------------------------------------*/
/* Pour le profiling sur accélérateur                                        */
/*---------------------------------------------------------------------------*/

#if defined(ARCANE_COMPILING_CUDA) && defined(PROF_ACC)
#define USE_PROF_ACC
#endif

#if defined(USE_PROF_ACC)

#warning "PROF_ACC : instrumentation avec nvtx"
#include <nvtx3/nvToolsExt.h>

#ifndef PROF_ACC_BEGIN
#define PROF_ACC_BEGIN(__name__) nvtxRangePushA(__name__)
#endif

#ifndef PROF_ACC_END
#define PROF_ACC_END nvtxRangePop()
#endif

#else

//#warning "Pas d'instrumentation"
#ifndef PROF_ACC_BEGIN
#define PROF_ACC_BEGIN(__name__)
#endif

#ifndef PROF_ACC_END
#define PROF_ACC_END 
#endif

#endif

/*---------------------------------------------------------------------------*/
/* "Conseiller" mémoire, permet de caractériser des accès mémoire            */
/*---------------------------------------------------------------------------*/
class AccMemAdviser {
 public:
  AccMemAdviser(bool enable=true) : 
    m_enable(enable) {
#ifdef ARCANE_HAS_CUDA
    if (m_enable) {
      cudaGetDevice(&m_device);
    }
#endif
  }
  ~AccMemAdviser() {}

  bool enable() const { return m_enable; }

  template<typename ViewType>
  void setReadMostly(ViewType view) {
#ifdef ARCANE_HAS_CUDA
    if (m_enable && view.size()) {
      cudaMemAdvise (view.data(), view.size(), cudaMemAdviseSetReadMostly,m_device);
    }
#endif
  }
 private:
  bool m_enable=true;
  int m_device=-1;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
