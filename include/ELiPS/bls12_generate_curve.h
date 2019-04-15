#ifndef BLS12_GENERATE_CURVE_H
#define BLS12_GENERATE_CURVE_H

#include <ELiPS/mpn.h>

/**
 * @brief generate X.
 *
 */
extern void BLS12_generate_X();

/**
 * @brief generate prime.
 *
 */
extern int  BLS12_generate_prime();

/**
 * @brief generate order.
 *
 */
extern int  BLS12_generate_order();

/**
 * @brief generate trace.
 *
 */
extern void BLS12_generate_trace();

/**
 * @brief weil's theorem.
 *
 */
extern void BLS12_weil();

#endif
