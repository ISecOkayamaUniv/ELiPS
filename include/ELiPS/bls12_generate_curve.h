#ifndef BLS12_GENERATE_CURVE_H
#define BLS12_GENERATE_CURVE_H

#include <ELiPS/mpn.h>

/**
 * @brief generate X.
 *
 */
void BLS12_generate_X();

/**
 * @brief generate prime.
 *
 */
int  BLS12_generate_prime();

/**
 * @brief generate order.
 *
 */
int  BLS12_generate_order();

/**
 * @brief generate trace.
 *
 */
void BLS12_generate_trace();

/**
 * @brief weil's theorem.
 *
 */
void BLS12_weil();

#endif
