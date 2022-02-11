#include <arcane/accelerator/IRunQueueStream.h>

#include "accenv/AcceleratorEvent.h"

AcceleratorEvent::AcceleratorEvent() {
#ifdef ARCANE_COMPILING_CUDA
  cudaEventCreateWithFlags(&(m_event), cudaEventDisableTiming);
#endif
}

AcceleratorEvent::~AcceleratorEvent() {
#ifdef ARCANE_COMPILING_CUDA
  cudaEventDestroy(m_event);
#endif
}

/*---------------------------------------------------------------------------*/
/* Enregistre un evenement sur la queue                                      */
/*---------------------------------------------------------------------------*/
void AcceleratorEvent::record(RunQueue& queue) {
#ifdef ARCANE_COMPILING_CUDA
  auto* rq = queue._internalStream();
  cudaStream_t* s = reinterpret_cast<cudaStream_t*>(rq->_internalImpl());

  cudaEventRecord(m_event, *s);
#endif
}

/*---------------------------------------------------------------------------*/
/* Bloque l'hôte tant que l'événement n'arrive pas                           */
/*---------------------------------------------------------------------------*/
void AcceleratorEvent::wait() {
#ifdef ARCANE_COMPILING_CUDA
  // L'hote bloque jusqu'à ce que l'evenement arrive
  cudaEventSynchronize(m_event);
#endif
}

/*---------------------------------------------------------------------------*/
/* Bloque les kernels suivant sur la queue tant que l'événement n'arrive pas */
/*---------------------------------------------------------------------------*/
void AcceleratorEvent::queueWait(RunQueue& queue) {
#ifdef ARCANE_COMPILING_CUDA
  auto* rq = queue._internalStream();
  cudaStream_t* s = reinterpret_cast<cudaStream_t*>(rq->_internalImpl());

  cudaStreamWaitEvent(*s, m_event);
#endif
}

