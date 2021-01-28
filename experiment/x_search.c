#include "ELiPS/bls12.h"
#include<stdlib.h>
fp_t alpha2;
fp2_t beta3;
fp6_t ganma2;
FILE *fp;
int curve_b_int;
int curve_b_13_flag;
int curve_b_23_flag;
int curve_b_3_pow;
fp_t curve_b_13,curve_b_23;
fp_t digit3;
int X_mod72(){
  mpz_t m;
  int ans;
  mpz_init(m);
  mpz_mod_ui(m,X_z,72);
  ans = (int)mpz_get_ui(m);
  mpz_clear(m);
  return ans;
}
void bls12_init_without_x_binary(){
    int i,j;

    mpz_init(X_z);
    mpz_init(prime_z);
    mpz_init(order_z);
    mpz_init(trace_z);
    mpz_init(root_X);
    mpz_init(root_2);

    mpz_init(efp_total);
    mpz_init(efp12_total);
    mpn_zero(curve_b,FPLIMB);

    // for(i=0; i<bls12_X_length+1; i++){
    //     bls12_X_binary[i]=0;
    // }

    // for(i=0; i<bls12_X2_length+1; i++){
    //     bls12_X2_binary[i]=0;
    // }

    mpn_zero(epsilon1,FPLIMB);
    mpn_zero(epsilon2,FPLIMB);
    mpz_init(Two_inv_z);
    fp2_init(&Alpha_1);
    fp2_init(&Alpha_1_inv);

    for(i=0; i<12; i++){
        for(j=0; j<6; j++){
            fp2_init(&frobenius_constant[i][j]);
        }
        for(j=0; j<2; j++){
            fp2_init(&skew_frobenius_constant[i][j]);
        }
    }
}
void bls12_clear_without_x_binary(){
  int i,j;

    mpz_clear(X_z);
    mpz_clear(prime_z);
    mpz_clear(order_z);
    mpz_clear(trace_z);
    mpz_clear(root_X);
    mpz_clear(root_2);

    mpz_clear(efp_total);
    mpz_clear(efp12_total);

    mpz_clear(Two_inv_z);
}
int p_test(){
  //gmp_randinit_default (state);
  //gmp_randseed_ui(state,(unsigned long)time(NULL));
  int i;
  static mpz_t buf;
  mpz_init(buf);
  mpz_set_ui(X_z,0);
  for(i=bls12_X_length; i>=0; i--){
    if(bls12_X_binary[i]==1){
      mpz_ui_pow_ui(buf,2,i);
      mpz_add(X_z,X_z,buf);
    }else if(bls12_X_binary[i]==-1){
      mpz_ui_pow_ui(buf,2,i);
      mpz_sub(X_z,X_z,buf);
      }
    }
  mpz_clear(buf);
  if(bls12_generate_prime()==1 && bls12_generate_order()==1){
    bls12_generate_trace();
    bls12_weil();
    bls12_get_epsilon();
    bls12_get_Two_inv();
    bls12_set_basis();
    bls12_set_frobenius_constant();
    //bls12_set_curve_parameter();
    bls12_set_root2();
	  //bls12_ry();
	  pre_montgomery();
  }else{
    //bls12_clear();
    return -1;
    //printf("error : prime\nexit\n");
  }
  return 1;
}
int fpx_test(){
  printf("fpx_test\n");
  //initialize irreducible polynomial
  static fp_t tmp;
  fp_init(&tmp);
  fp_init(&alpha2);
  fp2_init(&beta3);
  fp6_init(&ganma2);
  fp_set_ui(&tmp, 0);
  //alpha^2=-1
  fp_sub_ui(&alpha2,&tmp,1);
  //beta^3=alpha + 1
  fp_set_ui(&beta3.x0,1);
  fp_set_ui(&beta3.x1,1);
  //ganma^2=beta
  fp_set_ui(&ganma2.x1.x0,1);

  //check irreducible
  if(fp_legendre(&alpha2)!=-1) return -1;
  if(fp2_isCNR(&beta3)!=-1) return -1;
  if(fp6_legendre(&ganma2)!=-1) return -1;
  return 1;
}
int twist_type(){
  efp12_t P;
  efp12_init(&P);
  bls12_generate_g1(&P);
  efp12_scm(&P, &P, order_z);
  if (P.infinity==0)
    return 9999;
  /*
  1-------elips
  2-------relic
  other---error
  */
  efp12_t Q;
  efp12_init(&Q);
  printf("generate g2");
  bls12_generate_g2(&Q);
  if(fp_cmp_ui(&Q.x.x0.x1.x0,0)==0){
    return 1;
  }else if(fp_cmp_ui(&Q.x.x0.x2.x0,0)==0){
    return 2;
  }else{
    return 9999;
  }
}
int curve_b_test(){
  efp_t ANS,test_P;
  efp_init(&ANS);
  efp_init(&test_P);
  efp_set_random(&test_P);
  efp_scm(&ANS,&test_P,efp_total);
  if(ANS.infinity) return 1;
  else return -1;
}
int curve_b_make(){
  //generate bls12( b is minimum )
  curve_b_int=0;
  curve_b_13_flag=0;
  curve_b_23_flag=0;
  //1/3,2/3
  fp_set_ui(&digit3,3);
  fp_inv(&curve_b_13,&digit3);
  fp_inv(&curve_b_23,&digit3);
  fp_mul_ui(&curve_b_23,&curve_b_23,2);

  mpn_copyd(curve_b,curve_b_13.x0,FPLIMB);
  if(curve_b_test()==1) curve_b_13_flag=1;
  mpn_copyd(curve_b,curve_b_23.x0,FPLIMB);
  if(curve_b_test()==1) curve_b_23_flag=1;

  for(mpn_zero(curve_b,FPLIMB);mpn_cmp(curve_b,prime,FPLIMB);mpn_add_ui(curve_b,curve_b,FPLIMB,1)){

    //if(curve_b_test()==1) return 1;
    if(curve_b_test()==1){
      //printf("curveb\n");
      //gmp_printf("b=%Nu\n",curve_b,FPLIMB);
      //getchar();
      return 1;
    }
    curve_b_int++;
  }
  return -1;
}
// int curve_b_make_divided3(){
//   //generate bls12( b is minimum )
//   fp_t digit3,digit3_tmp;
//   int curve_int=1;
//   //1/3,2/3
//   fp_set_ui(&digit3,3);
//   fp_inv(&digit3,&digit3);
//   fp_set(&digit3_tmp,&digit3);

//   mpn_copyd(curve_b,digit3.x0,FPLIMB);
//   while(1){
//     //if(curve_b_test()==1) return 1;
//     if(curve_b_test()==1){
//       printf("curveb\n");
//       gmp_printf("b=%Nu\n",curve_b,FPLIMB);
//       printf("curve_b=%d/3",curve_int);
//       getchar();
//       //return 1;
//     }
//     fp_add(&digit3,&digit3,&digit3_tmp);
//     mpn_copyd(curve_b,digit3.x0,FPLIMB);
//     curve_int++;
//   }
//   return -1;
// }


int curve_b_make_divided3(){
  //find b such 2^n/3
  fp_t digit3,digit3_tmp;
  int curve_int=1;
  //1/3,2/3
  fp_set_ui(&digit3,3);
  fp_inv(&digit3,&digit3);
  fp_set(&digit3_tmp,&digit3);

  mpn_copyd(curve_b,digit3.x0,FPLIMB);
  while(1){
    if(curve_b_test()==1){
      printf("curveb\n");
      gmp_printf("b=%Nu\n",curve_b,FPLIMB);
      printf("curve_b=%d/3",curve_int);
      return 1;
    }
    fp_mul_ui(&digit3,&digit3,2);
    //mpn_copyd(curve_b,digit3.x0,FPLIMB);
    mpn_copyd(curve_b,digit3.x0,FPLIMB);
    curve_int=curve_int*2;
  }
  return -1;
}
int alltest_with_twist_type(){
  bls12_clear_without_x_binary();
  bls12_init_without_x_binary();
  int type=0;
  if(p_test()==-1) return -1;
  printf("p_test:success\n");
  if(fpx_test()==-1) return -1;
  printf("fpx_test:success\n");
  //
  if(curve_b_make_divided3()==-1) return -1;
  printf("curve_b_make:success\n");
  //
  if(curve_b_make()==-1) return -1;
  printf("curve_b_make:success\n");
  bls12_print_parameters();
  type = twist_type();
  if(type==9999) return -1;
  printf("twist_type:success\n");
  return type;
}
void X_init(){
  int v=0;
  for(v=bls12_X_length;v>=0;v--){
    bls12_X_binary[v]=0;
  }
}
void X_print(){
  int v=0;
  printf("X_print\n");
  for(v=bls12_X_length;v>=0;v--){
    if(bls12_X_binary[v]!=0) printf("bls12_X_binary[%d]=%d\n",v,bls12_X_binary[v]);
  }
  printf("curveb\n");
  gmp_printf("b=%Nu\n",curve_b,FPLIMB);

}
void X_test_and_print(){
  int test=0,v=0;
  test=alltest_with_twist_type();

  if(test!=-1){
    for(v=bls12_X_length;v>=0;v--){
      if(bls12_X_binary[v]!=0) fprintf(fp,"x[%d]=%d ",v,bls12_X_binary[v]);
    }
    //fprintf(fp,",%d,%d,%d,%d,%d\n",(int)mpz_sizeinbase(prime_z,2),curve_b_int,test,curve_b_13_flag,curve_b_23_flag);
    printf("twist type=%d\n",test);
    X_print();
  }
}
void X_test_and_print_tofile(){
  int test=0,v=0;
  test=alltest_with_twist_type();
  int xmod72;
  //mpz_init(xmod72);

  if(test!=-1){
    for(v=bls12_X_length;v>=0;v--){
      if(bls12_X_binary[v]!=0) fprintf(fp,"x[%d]=%d ",v,bls12_X_binary[v]);
    }
    xmod72 = X_mod72(xmod72);
    fprintf(fp,",%d,%d,%d,%d,%d,%d\n",(int)mpz_sizeinbase(prime_z,2),curve_b_int,test,curve_b_13_flag,curve_b_23_flag,xmod72);
    printf("twist type=%d\n",test);
    X_print();
  }
  //mpz_init(xmod72);
}
void X_search_return(int v_before, int hamming_weight)
{
  int v,u=0,p,tmp;
  static int binary[2]={1,-1};//element of X
  tmp = hamming_weight-1;
  for(u=0;u<2;u++){
    for(v=v_before-1;v>=0;v--){
      bls12_X_binary[v]=binary[u];
      if(tmp>0) X_search_return(v,tmp);
      //generated x
      if(tmp<=0) X_test_and_print_tofile();
      bls12_X_binary[v]=0;
    }
  }
}
void X_search_hamming_4(){
  X_init();
  bls12_X_binary[77]=1;
  X_search_return(77,3);
}
void X_search_hamming_4_minus(){
  X_init();
  bls12_X_binary[77]=-1;
  X_search_return(77,3);
}
void X_search_hamming_4_76(){
  X_init();
  bls12_X_binary[76]=1;
  X_search_return(76,3);
}
void X_search_hamming_3_76(){
  X_init();
  bls12_X_binary[76]=1;
  X_search_return(76,2);
}
void X_search_hamming_3_76_minus(){
  X_init();
  bls12_X_binary[76]=-1;
  X_search_return(76,2);
}

void montgomery_test(){
  fp_t A,B,C,C_2;

  fp_init(&A);
  fp_init(&B);
  fp_init(&C);
  fp_init(&C_2);

  fp_set_random(&A,state);
  fp_set_random(&B,state);

  fp_mul(&C,&A,&B);

  //montgomery
  fp_to_montgomery(&A,&A);
  fp_to_montgomery(&B,&B);
  fp_mulmod_montgomery(&C_2,&A,&B);
  fp_mod_montgomery(&C_2,&C_2);

  if(fp_cmp(&C_2,&C)==0) printf("montgomery ok!\n");
  else printf("unko\n");
}
int main(){
    int type=0;
  bls12_init();
  //montgomery_test();
  char *fname = "test.csv";
  char *s1 = "X,prime_bit,curve_b,twist_type";

  fp = fopen( fname, "w" );
  if( fp == NULL ){
    printf( "%sファイルが開けません¥n", fname );
    return -1;
  }

  fprintf( fp, "%s\n", s1);
  X_search_hamming_4();
  X_search_hamming_4_minus();

  //X_search_hamming_4_76_minus();

}
