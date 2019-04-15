#ifndef TWIST_H
#define TWIST_H

#include <ELiPS/EFp12.h>

/**
 * @brief Set EFp12 on G2 to EFp2.
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set on G2.
 */
extern void EFp12_to_EFp2(EFp2 *ANS,EFp12 *A);

/**
 * @brief Set EFp2 to EFp12.
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp2_to_EFp12(EFp12 *ANS,EFp2 *A);

/**
 * @brief Set EFp12 on G1 to EFp.
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set on G1.
 */
extern void EFp12_to_EFp(EFp *ANS,EFp12 *A);

/**
 * @brief Set EFp to EFp12.
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp_to_EFp12(EFp12 *ANS,EFp *A);

#endif
