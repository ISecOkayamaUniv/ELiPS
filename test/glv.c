#include "ELiPS/bls12.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
void fp12_frobenius_test(){
    printf("\n==============fp12 frobenius test===========\n");
    fp12_t A_fp12,A_fro,A_fro2;
    fp12_init(&A_fp12);
    fp12_init(&A_fro);
    fp12_set_random(&A_fp12,state);
    fp12_pow(&A_fro,&A_fp12,prime_z);
    fp12_println("A=",&A_fro);
    fp12_frobenius_map_p1(&A_fro2,&A_fp12);
    fp12_println("A=",&A_fro2);
    if(fp12_cmp(&A_fro,&A_fro2)==0) printf("ok!\n");
    else printf("ng!\n");
}
void g2_frobenius_test(){
    printf("\n=============G2 frobenius test===========\n");
    efp12_t A;
    efp2_t A_twisted,Q_fro,Q_fro2;
    mpz_t buf;
    mpz_init(buf);
#ifdef X_MINUS
    mpz_neg(buf,X_z);
#endif
#ifdef X_PLUS
    mpz_set(buf,X_z);
#endif
    efp12_init(&A);
    efp2_init(&Q_fro);
    efp2_init(&Q_fro2);
    bls12_generate_g2(&A);
    efp12_to_efp2(&A_twisted,&A);
    efp2_scm(&Q_fro,&A_twisted,prime_z);
    efp2_skew_frobenius_map_p1(&Q_fro2,&A_twisted);
// #ifdef X_MINUS
//     efp2_set_neg(&Q_fro2,&A_fro2);
// #endif
    efp2_println("xQ=",&Q_fro);
    efp2_println("xQ=",&Q_fro2);
    if(efp2_cmp(&Q_fro,&Q_fro2)==0) printf("ok!\n");
    else printf("ng!\n");
}
void efp12_frobenius_test(){
    printf("\n=============G2 frobenius test===========\n");
    efp12_t Q_fro,Q_fro2,Q;
    efp12_init(&Q_fro);
    efp12_init(&Q_fro2);
    efp12_init(&Q);

    bls12_generate_g2(&Q);
    efp12_println("A=",&Q);

    efp12_scm(&Q_fro,&Q,prime_z);
    fp12_frobenius_map_p1(&Q_fro2.x,&Q.x);
    fp12_frobenius_map_p1(&Q_fro2.y,&Q.y);
    efp12_println("A=",&Q_fro);
    efp12_println("A=",&Q_fro2);
    if(efp12_cmp(&Q_fro,&Q_fro2)==0) printf("ok!\n");
    else printf("ng!\n");

}
int main(void){
    bls12_init();
    bls12_print_parameters();
    /*fp12 frobenius test*/
    fp12_frobenius_test();
    getchar();
    g2_frobenius_test();
    getchar();
    efp12_frobenius_test();


    
    /*efp12 frobenius test*/
    

   

    //test_field(0,1000,0,0,0);
    //test_fp(1000);
    //test_fp2(10000);
    //test_fp6(100);
    //test_mod(1000000);
    //test_Frobenius_map();
    //test_skew_frobenius_map();
    //test_twist();
    //test_mod(10000000,100000);

    //test_efp(100,100);
    //test_efp2(1000,1000);
    //test_efp12(10,10,10);

    //bls12_test_rational_point();
    //bls12_test_opt_ate_pairing(100);
    //bls12_test_g1_scm(10000);
    //bls12_test_g2_scm(10000);
    //bls12_test_g3_exp(100);

    //test_All();
    return 0;
}
