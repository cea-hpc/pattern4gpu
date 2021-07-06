// -*- coding: utf-8 -*-
#ifndef _P4GPU_P4GPUTIMER_H_
#define _P4GPU_P4GPUTIMER_H_
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include "arcane/ISubDomain.h"
#include "arcane/ITimeStats.h"
#include "arcane/Timer.h"
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * @brief Timer, comme Arcane::Timer::Action avec une fonction start() et stop()
 *
 * @note Ne devrait pas etre utilisee directement, mais via les macros definies
 * plus bas, uniquement actives si la macro P4GPU_PROFILING est definie :
 *   @a P4GPU_FUNCTION_TIMER
 *   @a P4GPU_DECLARE_TIMER
 *   @a P4GPU_START_TIMER
 *   @a P4GPU_STOP_TIMER
 */
class P4GPUTimer
{
	public:
	P4GPUTimer(ISubDomain * sd, const String & name)
	: m_sd(sd), m_name(name), m_is_started(false)
	{}

	~P4GPUTimer()
	{ if ( m_is_started ) m_sd->timeStats()->endAction(m_name, false); }

	void start()
	{
		if ( m_is_started ) return;
		m_is_started = true;
		m_sd->timeStats()->beginAction(m_name);
	}

	void stop()
	{
		if ( m_is_started ) m_sd->timeStats()->endAction(m_name, false);
		m_is_started = false;
	}

	private:
	ISubDomain * m_sd;
	String m_name;
	bool m_is_started;
};
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#if defined(P4GPU_PROFILING)
	/**
	 * @brief Declare un timer qui mesure le temps passe dans une fonction.
	 *
	 * Mesure le temps passe entre la construction de ce timer et sa destruction.
	 * Si cette macro est placee au debut d'une fonction, il mesure le temps
	 * passe dans cette fonction.
	 */
	#define P4GPU_FUNCTION_TIMER(sub_domain) Arcane::Timer::Action p4gpu_function_timer(sub_domain, __func__);
	/**
	 * @brief Declare un timer qui mesure le temps entre 2 instants du calcul.
	 *
	 * Cette macro declare un timer. Il peut etre demarre avec @a P4GPU_START_TIMER
	 * et arrete avec @a P4GPU_STOP_TIMER. Il est automatiquement arrete lors
	 * de sa destruction.
	 */
	#define P4GPU_DECLARE_TIMER(sub_domain, timer_name) P4GPUTimer p4gpu_timer_##timer_name(sub_domain, #timer_name);
	/**
	 * @brief Demarre un timer declare avec @a P4GPU_DECLARE_TIMER.
	 */
	#define P4GPU_START_TIMER(timer_name) p4gpu_timer_##timer_name.start();
	/**
	 * @brief Arrete un timer declare avec @a P4GPU_DECLARE_TIMER.
	 */
	#define P4GPU_STOP_TIMER(timer_name) p4gpu_timer_##timer_name.stop();
#else
	#define P4GPU_FUNCTION_TIMER(sub_domain)
	#define P4GPU_DECLARE_TIMER(sub_domain, timer_name)
	#define P4GPU_START_TIMER(timer_name)
	#define P4GPU_STOP_TIMER(timer_name)
#endif
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#endif
