#include "ELiPS/bls12.h"
#include "ELiPS/fr.h"
/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void) {
  bls12_init();
  //test_fp_montgomery(100);
  int n=10000;
  g1_test(n);
  g2_test(n);
  g3_test(n);
  debug_pairing(n);
  g1_set_random_test(n);
  g2_set_random_test(n);

  return 0;
}
