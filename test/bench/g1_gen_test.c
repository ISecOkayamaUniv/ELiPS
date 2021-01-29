#include "bench.h"
//mpz_t sqrt_power_z;
//mpz_t g1_power;

// void power_init_a(){
//     mpz_t tmp;
//     mpz_init(tmp);
  
//     mpz_init(sqrt_power_z);
//     mpz_init(g1_power);

//     //set sqrt_power_z
//     mpz_sub_ui(sqrt_power_z,prime_z,3);
//     mpz_div_ui(sqrt_power_z,sqrt_power_z,4);

//     //set g1_power
//     mpz_tdiv_q(g1_power,efp_total,order_z);
// }

void fp_pow_montgomery(fp_t *ANS,fp_t *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fp_t tmp;
    fp_init(&tmp);//not need?
    fp_set(&tmp,A);

    for(i=1;i<length; i++){
        fp_mulmod_montgomery(&tmp,&tmp,&tmp);
        if(binary[i]=='1'){
            fp_mulmod_montgomery(&tmp,A,&tmp);
        }
    }
    fp_set(ANS,&tmp);

}






// void fp_sqrt_montgomery(fp_t *ANS,fp_t *A){
//     fp_pow_montgomery(ANS,A,sqrt_power_z);
// }

// int fp_legendre_montgomery(fp_t *A){
//     int i;
//     fp_t tmp1_fp;
//     fp_init(&tmp1_fp);
//     fp_pow_montgomery(&tmp1_fp,A,legendre_power_z);
//     if(mpn_cmp(tmp1_fp.x0,RmodP,FPLIMB)==0)        i=1;
//     else if(mpn_cmp_ui(tmp1_fp.x0,FPLIMB,0)==0)    i=0;
//     else                    i=-1;
//     return i;
// }

// void efp_set_random_montgomery(efp_t *P,gmp_randstate_t state){
//     fp_t tmp1,tmp2,tmp_x;
//     fp_init(&tmp1);
//     fp_init(&tmp2);
//     fp_init(&tmp_x);

//     while(1){
//         fp_set_random_montgomery(&P->x,state);
//         fp_mulmod_montgomery(&tmp1,&P->x,&P->x);
//         fp_mulmod_montgomery(&tmp2,&tmp1,&P->x);
//         fp_add_mpn(&tmp_x,&tmp2,curve_b);
//         if(fp_legendre_montgomery(&tmp_x)==1){
//             fp_sqrt_montgomery(&P->y,&tmp_x);
//             break;
//         }
//     }
//     P->infinity=0;
// }

// void efp_scm_montgomery(efp_t *ANS,efp_t *P,mpz_t scalar){
//      if(mpz_cmp_ui(scalar,0)==0){
//         ANS->infinity=1;
//         return;
//     }else if(mpz_cmp_ui(scalar,1)==0){
//         efp_set(ANS,P);
//         return;
//     }

//     efp_t Tmp_P,Next_P;
//     efp_init(&Tmp_P);
//     efp_set(&Tmp_P,P);
//     efp_init(&Next_P);
//     int i,length;
//     length=(int)mpz_sizeinbase(scalar,2);
//     char binary[length+1];
//     mpz_get_str(binary,2,scalar);

//     efp_set(&Next_P,&Tmp_P);
//     for(i=1;i<length; i++){
//         efp_ecd(&Next_P,&Next_P);
//         if(binary[i]=='1'){
//             efp_eca(&Next_P,&Next_P,&Tmp_P);
//         }
//     }
//     efp_set(ANS,&Next_P);
// }

// void g1_set_random_montgomery(g1_t *P,gmp_randstate_t state){
//     g1_t tmp;
//     g1_init(&tmp);
//     efp_set_random_montgomery(&tmp,state);
//     efp_scm_montgomery(P,&tmp,g1_power);
// }

// void g1_set_random_montgomery(g1_t *ANS, gmp_randstate_t state){
//     efp12_t P;
//     g1_t P_twisted;
//     bls12_generate_g1_montgomery(&P);
//     efp12_to_efp(&P_twisted,&P);
//     //efp_to_montgomery(ANS,&P_twisted);
//     ANS->infinity=0;
// }

void sqrt_test(int cnt){
  fp_t A,A_sqr,tmp,A_neg;
  mpz_t two;
  mpz_init(two);
  mpz_set_ui(two,2);
  fp_set_random(&A,state);
  fp_pow(&A_sqr,&A,two);
  fp_set_neg(&A_neg,&A);
  fp_println("A=",&A);
  fp_println("A_sqr=",&A_sqr);
  int legendre=fp_legendre(&A_sqr);
  printf("legendre %d\n",legendre);

  //normal
  fp_sqrt(&tmp,&A_sqr);
  if(fp_cmp(&tmp,&A)==0 || fp_cmp(&tmp,&A_neg)==0) printf("ok!\n");
  else printf("wrong!\n");

  //fast
  fp_legendre_sqrt(&tmp,&A_sqr);
  if(fp_cmp(&tmp,&A)==0 || fp_cmp(&tmp,&A_neg)==0) printf("ok!\n");
  else printf("wrong!\n");

  //fast_montgomery
//   fp_to_montgomery(&A_sqr,&A_sqr);
//   fp_sqrt_montgomery(&tmp,&A_sqr);
//   fp_mod_montgomery(&tmp,&tmp);
//   if(fp_cmp(&tmp,&A)==0 || fp_cmp(&tmp,&A_neg)==0) printf("ok!\n");
//   else printf("wrong!\n");

//   int i;
//   for(i=0;i<cnt;i++){

//   }
}


int main(){
  bls12_init();
  //power_init();
  bls12_print_parameters();
  //fr_order_init();
  //power_init_a();
  efp_t A;
  efp_set_random_fast(&A,state);
  efp12_t PP;
  bls12_generate_g1_fast(&PP);

  g1_set_random_test(100);
  // fp_t A,A_sqr,tmp,A_neg;
  // mpz_t two;
  // mpz_init(two);
  // mpz_set_ui(two,2);
  // fp_set_random(&A,state);
  // fp_pow(&A_sqr,&A,two);
  // fp_set_neg(&A_neg,&A);
  // fp_println("A=",&A);
  // fp_println("A_sqr=",&A_sqr);
  // int legendre=fp_legendre(&A_sqr);
  // printf("legendre %d\n",legendre);

  // //normal
  // fp_sqrt(&tmp,&A_sqr);
  // if(fp_cmp(&tmp,&A)==0 || fp_cmp(&tmp,&A_neg)==0) printf("ok!\n");
  // else printf("wrong!\n");

  // //fast
  // fp_sqrt_fast(&tmp,&A_sqr);
  // if(fp_cmp(&tmp,&A)==0 || fp_cmp(&tmp,&A_neg)==0) printf("ok!\n");
  // else printf("wrong!\n");
//   g1_t A;
//   fr_t sca;
//   //efp_set_random_montgomery(&A);
//   g1_set_random_montgomery(&A,state);
//   fr_set_mpn(&sca,order);
//   efp_mod_montgomery(&A,&A);
//   efp_scm(&A,&A,order_z);
//   efp_printf("A=",&A);


  return 0;
}
