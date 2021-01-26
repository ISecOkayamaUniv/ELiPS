#ifndef BLS12_GENERATE_PARAMETER_H
#define BLS12_GENERATE_PARAMETER_H

#include <ELiPS/fp2.h>


/**
 * @brief get 1 cube root.
 *
 */
extern void bls12_get_epsilon();

/**
 * @brief get inverse element of 2.
 *
 */
extern void bls12_get_Two_inv();

/**
 * @brief set basis.
 *
 */
extern void bls12_set_basis();

/**
 * @brief set frobenius constat.
 *
 */
extern void bls12_set_frobenius_constant();

/**
 * @brief set b (y^2=ax^3+[b]).
 *
 */
extern void bls12_set_curve_parameter();

/**
 * @brief set square root of 2.
 *
 */
extern void bls12_set_root2();

/**
 * @brief set montgomery reduction paramater.
 *
 */
extern void bls12_set_montgomery();

#endif
