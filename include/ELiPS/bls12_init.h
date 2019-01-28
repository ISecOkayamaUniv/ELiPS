#ifndef BLS12_INIT_H
#define BLS12_INIT_H

#include <ELiPS/bls12_generate_curve.h>
#include <ELiPS/bls12_generate_parameter.h>
#include <ELiPS/Fp12.h>
#include <ELiPS/EFp12.h>


/**
 * @brief Initializes a BLS curve.
 *
 */
void BLS12_init();

/**
 * @brief Initializes a BLS curve parameters.
 *
 */
void BLS12_init_parameters();

/**
 * @brief Clear a BLS curve parameters.
 *
 */
void BLS12_clear();

/**
 * @brief Print a BLS curve parameters.
 *
 */
void BLS12_print_parameters();

#endif
