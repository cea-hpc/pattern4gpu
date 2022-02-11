#ifndef ACC_ENV_ACCELERATOR_EVENT_H
#define ACC_ENV_ACCELERATOR_EVENT_H

#include "accenv/AcceleratorUtils.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/* Encapsule un événement sur le DEVICE                                      */
/*---------------------------------------------------------------------------*/
class AcceleratorEvent {
 public:
  AcceleratorEvent();
  virtual ~AcceleratorEvent();

  // Enregistre l'événement au sein d'une RunQueue
  void record(ax::RunQueue& queue);

  // Attend l'occurence de l'événement
  void wait();

  // Bloque les kernels suivant sur la queue tant que l'événement n'arrive pas
  void queueWait(ax::RunQueue& queue);

 protected:
#ifdef ARCANE_COMPILING_CUDA
  cudaEvent_t m_event;
#endif  
};

#endif
