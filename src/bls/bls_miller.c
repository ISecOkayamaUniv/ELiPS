#include <ELiPS/bls_miller.h>

void bls_miller_for_optate_basic(fpm2_t *ANS, efpm2_t *p, efpm2_t *q){
    fpm2_t tmp1,tmp2;
    fpm2_t f;
    efpm2_t p1,q1,t;
    
    int i,length;
    int binary[512];
    length=(int)mpz_sizeinbase(order_z,2);
    printf("length = %d\n",length);
    //init
    efpm2_init(&p1);
    efpm2_init(&q1);
    efpm2_init(&t);
    for(i=0;i<512;i++){
        binary[i] = 0;
    }

    //set
    efpm2_set(&p1,p);
    efpm2_set(&q1,q);

    mpz_t s;
    mpz_init(s);
    mpz_set(s,order_z);
    for(i=0;mpz_popcount(s);i++){
        if(mpz_scan1(s,0)==i){
            binary[i] = 1;
            mpz_clrbit(s,i);
        }
    }
    printf("binary");
    for(i=0;i<256;i++){
        printf("%d ",binary[i]);
    }


    //step1
    fpm2_set_ui(&f,1);
    efpm2_set(&t,&p1);

    for(i=308-1;i>=0;i--){
        fpm2_mul(&tmp1,&f,&f);
        bls_f_ltt(&tmp2,&q1,&t);
        fpm2_mul(&f,&tmp1,&tmp2);
        efpm2_ecd(&t,&t);
        if(binary[i]){
            bls_f_ltp(&tmp1,&p1,&q1,&t);
            fpm2_mul(&f,&f,&tmp1);
            efpm2_eca(&t,&t,&p1);
        }
    }
    //check
    efpm2_println("t =",&t);
    //set
    fpm2_set(ANS,&f);
    /*
    for(i=length-1; i>=0; i--){
        switch(bls_X_binary[i]){
            case 0:

                ff_ltt(&f,&T,&mapped_P,&L);
                break;
            case 1:
                ff_ltt(&f,&T,&mapped_P,&L);
                f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
                break;
            case -1:
                ff_ltt(&f,&T,&mapped_P,&L);
                f_ltq(&f,&T,&mapped_Q_neg,&mapped_P,&L);
                break;
            default:
                break;
        }
    }
    */


}