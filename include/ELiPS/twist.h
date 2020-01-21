#ifndef TWIST_H
#define TWIST_H

#include <ELiPS/efp12.h>

/**
 * @brief Set efp12_t on G2 to efp2_t.
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set on G2.
 */
extern void efp12_to_efp2(efp2_t *ANS,efp12_t *A);

/**
 * @brief Set efp2_t to efp12_t.
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp2_to_efp12(efp12_t *ANS,efp2_t *A);

/**
 * @brief Set efp12_t on G1 to efp_t.
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set on G1.
 */
extern void efp12_to_efp(efp_t *ANS,efp12_t *A);

/**
 * @brief Set efp_t to efp12_t.
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp_to_efp12(efp12_t *ANS,efp_t *A);

#endif
