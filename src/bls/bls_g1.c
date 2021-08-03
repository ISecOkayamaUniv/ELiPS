#include <ELiPS/bls_g1.h>

void bls_generate_g1_for_efpm2(efpm2_t *p) {
  efp_t p1, p2;
  mpz_t tmp1;
  //init
  efp_init(&p1);
  efp_init(&p2);
  mpz_init(tmp1);
  mpz_tdiv_q(tmp1, efp_total, order_z);

  efp_set_random(&p1);
  efp_scm(&p2, &p1, tmp1);
  //set
  fpm2_set_fp(&p->x, &p2.x);
  fpm2_set_fp(&p->y, &p2.y);
  mpz_clear(tmp1);
}
