#include <ELiPS/bls_g2.h>

void bls_generate_g2_for_efpm2(efpm2_t *p) {
  mpz_t exp;
  efpm2_t p1, p2;
  //init
  mpz_init(exp);
  efpm2_init(&p1);
  efpm2_init(&p2);
  mpz_tdiv_q(exp, efpm2_total, order_z);
  mpz_tdiv_q(exp, exp, order_z);
  efpm2_rational_point(&p1);
  efpm2_scm(&p1, &p1, exp);
  fpm2_frobenius(&p2.x, &p1.x);
  fpm2_frobenius(&p2.y, &p1.y);
  fpm2_set_neg(&p1.y, &p1.y);
  efpm2_eca(p, &p2, &p1);
  //clear
  mpz_clear(exp);
}