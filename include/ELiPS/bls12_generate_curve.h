#ifndef BLS12_GENERATE_CURVE_H
#define BLS12_GENERATE_CURVE_H

#include <ELiPS/mpn.h>

/**
 * @brief generate X.
 *
 */
extern void bls12_generate_X();

/**
 * @brief generate prime.
 *
 */
extern int  bls12_generate_prime();

/**
 * @brief generate order.
 *
 */
extern int  bls12_generate_order();

/**
 * @brief generate trace.
 *
 */
extern void bls12_generate_trace();

/**
 * @brief weil's theorem.
 *
 */
extern void bls12_weil();

#endif
