#ifndef BLS12_INIT_H
#define BLS12_INIT_H

#include <ELiPS/bls12_generate_curve.h>
#include <ELiPS/bls12_generate_parameter.h>
#include <ELiPS/efp12.h>
#include <ELiPS/fr.h>


/**
 * @brief Initializes a bls curve.
 *
 */
extern void bls12_init();

/**
 * @brief Initializes a bls curve parameters.
 *
 */
extern void bls12_init_parameters();

/**
 * @brief Clear a bls curve parameters.
 *
 */
extern void bls12_clear();

/**
 * @brief Print a bls curve parameters.
 *
 */
extern void bls12_print_parameters();

#endif
