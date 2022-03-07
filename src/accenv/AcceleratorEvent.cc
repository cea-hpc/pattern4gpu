#include "accenv/AcceleratorEvent.h"

AcceleratorEvent::AcceleratorEvent(ax::Runner& runner) :
  m_event (makeEvent(runner))
{
}

AcceleratorEvent::~AcceleratorEvent() {
}

/*---------------------------------------------------------------------------*/
/* Enregistre un evenement sur la queue                                      */
/*---------------------------------------------------------------------------*/
void AcceleratorEvent::record(RunQueue& queue) {
  queue.recordEvent(m_event);
}

/*---------------------------------------------------------------------------*/
/* Bloque l'hôte tant que l'événement n'arrive pas                           */
/*---------------------------------------------------------------------------*/
void AcceleratorEvent::wait() {
  m_event.wait();
}

/*---------------------------------------------------------------------------*/
/* Bloque les kernels suivant sur la queue tant que l'événement n'arrive pas */
/*---------------------------------------------------------------------------*/
void AcceleratorEvent::queueWait(RunQueue& queue) {
  queue.waitEvent(m_event);
}

