#ifndef BLS12_GENERATE_PARAMETER_H
#define BLS12_GENERATE_PARAMETER_H

#include <ELiPS/Fp2.h>


/**
 * @brief get 1 cube root.
 *
 */
void BLS12_get_epsilon();

/**
 * @brief get inverse element of 2.
 *
 */
void BLS12_get_Two_inv();

/**
 * @brief set basis.
 *
 */
void BLS12_set_basis();

/**
 * @brief set frobenius constat.
 *
 */
void BLS12_set_frobenius_constant();

/**
 * @brief set b (y^2=ax^3+[b]).
 *
 */
void BLS12_set_curve_parameter();

/**
 * @brief set square root of 2.
 *
 */
void BLS12_set_root2();

#endif
