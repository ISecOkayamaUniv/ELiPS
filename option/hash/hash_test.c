#include "elips_sha512.h"  // hash header.
int main() {
  char m[1000] = "Decode test is success!";
  char s1s2F_hash[64];   //sha512 -> 256bit -> 64byte
  char _s1s2F_hash[64];  //sha512 -> 256bit -> 64byte
  char enc_m[1000] = "Decode test is success!";
  char dec_m[1000] = "Decode test is success!";

  bls12_init();
  bls12_hash_init();  //hash init. It is necessary.
  g1_t P, s1P;
  g2_t Q, s2Q;
  g3_t F, s1s2F, _s1s2F;
  fr_t s1, s2, s1s2;

  g1_init(&P);
  g2_init(&Q);
  g3_init(&F);

  g1_hash_from_str(&P, "I'm g1 hash function. Please input string here.");  //g1's hash from char
  g2_hash_from_str(&Q, "I'm g2 hash function. Please input string here.");  //g2's hash from char

  fr_hash_from_str(&s1, "I'm fr hash function. Please input string here.");
  fr_hash_from_str(&s2, "I'm fr hash function. Please input string here.");
  fr_mul(&s1s2, &s1, &s2);

  g1g2_to_g3_pairing(&F, &P, &Q);
  g3_exp(&s1s2F, &F, &s1s2);
  hash_set_g3(s1s2F_hash, &s1s2F);
  for (int i = 0; i < 64; i++) {
    enc_m[i] = m[i] ^ s1s2F_hash[i];
  }
  printf("enc_message=%s\n", enc_m);

  g1_scm(&s1P, &P, &s1);
  g2_scm(&s2Q, &Q, &s2);
  g1g2_to_g3_pairing(&_s1s2F, &s1P, &s2Q);
  hash_set_g3(_s1s2F_hash, &s1s2F);
  for (int i = 0; i < 64; i++) {
    dec_m[i] = enc_m[i] ^ _s1s2F_hash[i];
  }
  printf("dec_message=%s\n", dec_m);

  return 0;
}