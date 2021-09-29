#include <ELiPS/bls12.h>
#include <gmp.h>
#include <openssl/sha.h>
#include <stdio.h>
//#if ARCBIT == 64
#define sha512LIMB 4
#define G2_STR_SIZE (FPLIMB * 64 + 10)       //(4 * FPLIMB * 64 / 4 + 10)
#define SHA_STR_SIZE (sha512LIMB * 16 + 10)  //(sha512LIMB * 64 / 4 + 10)
//#endif

// #if ARCBIT == 32
// #define sha512LIMB 16
// #define G2_STR_SIZE (FPLIMB * 32 + 10)   //(4 * FPLIMB * 32 / 4 + 10)
// #define SHA_STR_SIZE (sha512LIMB * 8 + 10)  //(sha512LIMB * 32 / 4 + 10)
// #endif

g1_t g1_base;
g2_t g2_base;
g3_t g3_base;
typedef struct {
  mp_limb_t x0[sha512LIMB];
} sha512_t;

void sha512_init(sha512_t *A) {
  mpn_zero(A->x0, sha512LIMB);
}

void sha512_print(char *str, sha512_t *A) {
  gmp_printf("%s%Nu", str, A->x0, sha512LIMB);
  //gmp_printf("%s%Nx", str, A->x0, sha512LIMB);
}

void sha512_println(char *str, sha512_t *A) {
  gmp_printf("%s%Nx\n", str, A->x0, sha512LIMB);
  //gmp_printf("%s%Nx\n", str, A->x0, sha512LIMB);
}

int sha512_cmp(sha512_t *A, sha512_t *B) {
  if (mpn_cmp(A->x0, B->x0, sha512LIMB) == 0) {
    return 0;
  } else {
    return 1;
  }
}

void sha512_mul_256(sha512_t *ANS, sha512_t *A) {
  mpn_lshift(ANS->x0, A->x0, sha512LIMB, 8);
}

void sha512_add_ui(sha512_t *ANS, sha512_t *A, unsigned long int UI) {
  mpn_add_ui(ANS->x0, A->x0, sha512LIMB, UI);
}

void sha512_from_str(sha512_t *ans, char *msg) {
  unsigned char digest[64];
  SHA512_CTX sha_ctx;
  SHA512_Init(&sha_ctx);                      //コンテキストを初期化
  SHA512_Update(&sha_ctx, msg, strlen(msg));  //msgを入力にする
  SHA512_Final(digest, &sha_ctx);             //digestに出力

#if DEBUG_PRINT
  printf("msg = %s\n", msg);
  for (int i = 0; i < 64; i++) printf("%02x ", digest[i]);
  printf("\n");
#endif

  sha512_init(ans);
  for (int i = 0; i < 64; i++) {
    sha512_mul_256(ans, ans);  //unsigned cahr = 1byte = 8bit = 2^8 = 256
    sha512_add_ui(ans, ans, (unsigned long int)digest[i]);
  }

#if DEBUG_PRINT
  sha512_println("hash_msg = ", ans);
  printf("\n\n");
#endif
}

void fr_set_sha512(fr_t *ANS, sha512_t *A) {
  if (sha512LIMB < FRLIMB) {
    mpn_copyd(ANS->x0, A->x0, sha512LIMB);
  } else {
    mp_limb_t dumy[sha512LIMB];
    mpn_tdiv_qr(dumy, ANS->x0, 0, A->x0, sha512LIMB, order, FRLIMB);
  }
}

void bls12_hash_init() {
  // for determined point
  gmp_randstate_t s;
  gmp_randinit_default(s);
  gmp_randseed_ui(s, 1);
  g1_init(&g1_base);
  g2_init(&g2_base);
  g3_init(&g3_base);

  g1_set_random(&g1_base, s);
  g2_set_random(&g2_base, s);
  g1g2_to_g3_pairing(&g3_base, &g1_base, &g2_base);

  g1_println("g1_base=", &g1_base);
  g2_println("g2_base=", &g2_base);
  g3_println("g3_base=", &g3_base);
}

void fr_hash_from_str(fr_t *ANS, char *m) {
  sha512_t s;
  sha512_init(&s);
  sha512_from_str(&s, m);
  fr_set_sha512(ANS, &s);
}

void g1_hash_from_str(g1_t *ANS, char *m) {
  fr_t s_fr;
  fr_init(&s_fr);
  fr_hash_from_str(&s_fr, m);
  g1_scm(ANS, &g1_base, &s_fr);
}

void g2_hash_from_str(g2_t *ANS, char *m) {
  fr_t s_fr;
  fr_init(&s_fr);
  fr_hash_from_str(&s_fr, m);
  g2_scm(ANS, &g2_base, &s_fr);
}

void g3_hash_from_str(g3_t *ANS, char *m) {
  fr_t s_fr;
  fr_init(&s_fr);
  fr_hash_from_str(&s_fr, m);
  g3_exp(ANS, &g3_base, &s_fr);
}

void hash_set_g1(char *ANS, g1_t *A) {
  if (A->infinity)
    printf("hash set g1: infinity error\n");
  else {
    mpz_t total;
    mpz_t tmp;
    mpz_init(total);
    mpz_init(tmp);
    mpz_set_mpn_size(tmp, A->x.x0, FPLIMB_BITS, FPLIMB);
    mpz_add(total, total, tmp);
    mpz_mul(total, total, prime_z);
    mpz_set_mpn_size(tmp, A->y.x0, FPLIMB_BITS, FPLIMB);
    mpz_add(total, total, tmp);
    mpz_mul(total, total, prime_z);
    char *str;
    str = (char *)malloc(mpz_sizeinbase(total, 10) + 2);
    mpz_get_str(str, 10, total);
    //printf("str = %s\n", str);
    sha512_from_str(&ANS, str);
    //sha512_println("a = ", ANS);

    mpz_clear(total);
    mpz_clear(tmp);
  }
}

void hash_set_g2(char *ANS, g2_t *A) {
  if (A->infinity)
    printf("hash set g2: infinity error\n");
  else {
    mpz_t total;
    mpz_t tmp;
    mpz_init(total);
    mpz_init(tmp);

    mpz_set_mpn_size(tmp, A->x.x0.x0, FPLIMB_BITS, FPLIMB);
    mpz_add(total, total, tmp);
    mpz_mul(total, total, prime_z);
    mpz_set_mpn_size(tmp, A->x.x1.x0, FPLIMB_BITS, FPLIMB);
    mpz_add(total, total, tmp);
    mpz_mul(total, total, prime_z);
    mpz_set_mpn_size(tmp, A->y.x0.x0, FPLIMB_BITS, FPLIMB);
    mpz_add(total, total, tmp);
    mpz_mul(total, total, prime_z);
    mpz_set_mpn_size(tmp, A->y.x1.x0, FPLIMB_BITS, FPLIMB);
    mpz_add(total, total, tmp);
    mpz_mul(total, total, prime_z);
    char *str;
    str = (char *)malloc(mpz_sizeinbase(total, 10) + 2);
    mpz_get_str(str, 10, total);
    // printf("str = %s\n", str);
    sha512_from_str(ANS, str);
    // sha_println(ans);

    mpz_clear(total);
    mpz_clear(tmp);
  }
}

void hash_set_g3(char *ANS, g3_t *A) {
  mpz_t total;
  mpz_t tmp;
  mpz_init(total);
  mpz_init(tmp);
  sha512_t ans_sha;

  mpz_set_mpn_size(tmp, A->x0.x0.x0.x0, FPLIMB_BITS, FPLIMB);
  mpz_add(total, total, tmp);
  mpz_mul(total, total, prime_z);
  mpz_set_mpn_size(tmp, A->x0.x0.x1.x0, FPLIMB_BITS, FPLIMB);
  mpz_add(total, total, tmp);
  mpz_mul(total, total, prime_z);
  mpz_set_mpn_size(tmp, A->x0.x1.x0.x0, FPLIMB_BITS, FPLIMB);
  mpz_add(total, total, tmp);
  mpz_mul(total, total, prime_z);
  mpz_set_mpn_size(tmp, A->x0.x1.x1.x0, FPLIMB_BITS, FPLIMB);
  mpz_add(total, total, tmp);
  mpz_mul(total, total, prime_z);
  mpz_set_mpn_size(tmp, A->x0.x2.x0.x0, FPLIMB_BITS, FPLIMB);
  mpz_add(total, total, tmp);
  mpz_mul(total, total, prime_z);
  mpz_set_mpn_size(tmp, A->x0.x2.x1.x0, FPLIMB_BITS, FPLIMB);
  mpz_add(total, total, tmp);
  mpz_mul(total, total, prime_z);
  mpz_set_mpn_size(tmp, A->x1.x0.x0.x0, FPLIMB_BITS, FPLIMB);
  mpz_add(total, total, tmp);
  mpz_mul(total, total, prime_z);
  mpz_set_mpn_size(tmp, A->x1.x0.x1.x0, FPLIMB_BITS, FPLIMB);
  mpz_add(total, total, tmp);
  mpz_mul(total, total, prime_z);
  mpz_set_mpn_size(tmp, A->x1.x1.x0.x0, FPLIMB_BITS, FPLIMB);
  mpz_add(total, total, tmp);
  mpz_mul(total, total, prime_z);
  mpz_set_mpn_size(tmp, A->x1.x1.x1.x0, FPLIMB_BITS, FPLIMB);
  mpz_add(total, total, tmp);
  mpz_mul(total, total, prime_z);
  mpz_set_mpn_size(tmp, A->x1.x2.x0.x0, FPLIMB_BITS, FPLIMB);
  mpz_add(total, total, tmp);
  mpz_mul(total, total, prime_z);
  mpz_set_mpn_size(tmp, A->x1.x2.x1.x0, FPLIMB_BITS, FPLIMB);
  mpz_add(total, total, tmp);

  char *str;
  str = (char *)malloc(mpz_sizeinbase(total, 10) + 2);
  mpz_get_str(str, 10, total);
  // printf("str = %s\n", str);
  sha512_from_str(ANS, str);
  // sha512_println("a = ", ANS);

  mpz_clear(total);
  mpz_clear(tmp);
}