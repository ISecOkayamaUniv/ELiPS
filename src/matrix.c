#include <ELiPS/matrix.h>


void matrix_init(matrix_t *a){
  int i,j;
  for(i=0;i<MATRIX_SIZE;i++){
    for(j=0;j<MATRIX_SIZE;j++){
      fp_init(&a->x[i][j]);
    }
  }
}
void matrix_set(matrix_t *b, matrix_t *a){
  int i,j;
  for(i=0;i<MATRIX_SIZE;i++){
    for(j=0;j<MATRIX_SIZE;j++){
      fp_set(&b->x[i][j],&a->x[i][j]);
    }
  }
}
void matrix_set_unit(matrix_t *a){
  int i,j;
  for(i=0;i<MATRIX_SIZE;i++){
    for(j=0;j<MATRIX_SIZE;j++){
      if(i==j)
        fp_set_ui(&a->x[i][j],1);
      else
        fp_init(&a->x[i][j]);
    }
  }
}
void matrix_set_random(matrix_t *a, gmp_randstate_t state){
  int i,j;
  for(i=0;i<MATRIX_SIZE;i++){
    for(j=0;j<MATRIX_SIZE;j++){
      fp_set_random(&a->x[i][j],state);
    }
  }
}
int matrix_cmp(matrix_t *a, matrix_t *b){
  int i,j;
  int result = 0;

  for(i=0;i<MATRIX_SIZE;i++){
    for(j=0;j<MATRIX_SIZE;j++){
      if(fp_cmp(&a->x[i][j],&b->x[i][j])!=0){
        result = 1;
      }
    }
  }
  return result;
}
void matrix_printf(char *str, matrix_t *a){
  int i,j;
  for(i=0;i<MATRIX_SIZE;i++){
    for(j=0;j<MATRIX_SIZE;j++){
      printf("Row %d Column %d ",i,j);
      fp_printf("",&a->x[i][j]);
      printf("\n");
    }
  }
}
void matrix_mul(matrix_t *c, matrix_t *a, matrix_t *b){
  int i,j,k;
  fp_t tmp1;
  matrix_t ans;
  matrix_init(&ans);
  for(i=0;i<MATRIX_SIZE;i++){
    for(j=0;j<MATRIX_SIZE;j++){
      for(k=0;k<MATRIX_SIZE;k++){
        fp_mul(&tmp1,&a->x[i][k],&b->x[k][j]);
        fp_add(&ans.x[i][j],&ans.x[i][j],&tmp1);
      }
    }
  }
  matrix_set(c,&ans);
}
void matrix_inv(matrix_t *b, matrix_t *a){
  int i,j,k;
  fp_t tmp1,tmp2;
  matrix_t t;
  matrix_t t_inv;

  matrix_set(&t,a);
  matrix_set_unit(&t_inv);
  for(i=0;i<MATRIX_SIZE;i++){
    fp_inv(&tmp1,&t.x[i][i]);
    for(j=0;j<MATRIX_SIZE;j++){
      fp_mul(&t.x[i][j],&t.x[i][j],&tmp1);
      fp_mul(&t_inv.x[i][j],&t_inv.x[i][j],&tmp1);
    }
    for(j=0;j<MATRIX_SIZE;j++){
      if(i!=j){
        fp_set(&tmp1,&t.x[j][i]);
        for(k=0;k<MATRIX_SIZE;k++){
          fp_mul(&tmp2,&t.x[i][k],&tmp1);
          fp_sub(&t.x[j][k],&t.x[j][k],&tmp2);
          fp_mul(&tmp2,&t_inv.x[i][k],&tmp1);
          fp_sub(&t_inv.x[j][k],&t_inv.x[j][k],&tmp2);
        }
      }
    }
  }
  matrix_set(b,&t_inv);
}
/*
void matrix_build_fp12_and_fp12cv(){
  int i,j;
  matrix_t matrix_fp12,matrix_fp12_inv;
  matrix_t matrix_fp12cv,matrix_fp12cv_inv;
  fp12_t tmp[12],buf;
  fp12cv_t tmp_cv[12],buf_cv;
  mpz_t exp;
  mpz_init(exp);
  mpz_pow_ui(exp,prime_z,12);
  mpz_sub_ui(exp,exp,1);
  mpz_tdiv_q_ui(exp,exp,13);//exp = (p^12-1)/13

  matrix_t a1,a2,e;
  matrix_init(&a1);
  matrix_init(&a2);
  matrix_set_unit(&e);

  //matrix karatsuba
  while(1){
    fp12_set_random(&buf,state);
    fp12_pow(&buf,&buf,exp);
    if(fp12_cmp_one(&buf)!=0){
      break;
    }
  }
  fp12_set(&tmp[0],&buf);
  for(i=1;i<12;i++){
    fp12_frobenius_map_p1(&tmp[i],&tmp[i-1]);
  }
  for(i=0;i<12;i++){
    fp_set(&matrix_fp12.x[i][0],&tmp[i].x0.x0.x0);
    fp_set(&matrix_fp12.x[i][1],&tmp[i].x0.x0.x1);
    fp_set(&matrix_fp12.x[i][4],&tmp[i].x0.x1.x0);
    fp_set(&matrix_fp12.x[i][5],&tmp[i].x0.x1.x1);
    fp_set(&matrix_fp12.x[i][8],&tmp[i].x0.x2.x0);
    fp_set(&matrix_fp12.x[i][9],&tmp[i].x0.x2.x1);
    fp_set(&matrix_fp12.x[i][2],&tmp[i].x1.x0.x0);
    fp_set(&matrix_fp12.x[i][3],&tmp[i].x1.x0.x1);
    fp_set(&matrix_fp12.x[i][6],&tmp[i].x1.x1.x0);
    fp_set(&matrix_fp12.x[i][7],&tmp[i].x1.x1.x1);
    fp_set(&matrix_fp12.x[i][10],&tmp[i].x1.x2.x0);
    fp_set(&matrix_fp12.x[i][11],&tmp[i].x1.x2.x1);
  }
  matrix_inv(&matrix_fp12_inv,&matrix_fp12);

  //matrix cvma
  while(1){
    fp12cv_set_random(&buf_cv,state);
    fp12cv_pow(&buf_cv,&buf_cv,exp);
    if(fp12cv_cmp_one(&buf_cv)!=0){
      break;
    }
  }
  
  fp12cv_set(&tmp_cv[0],&buf_cv);
  for(i=1;i<12;i++){
    fp12cv_frobenius(&tmp_cv[i],&tmp_cv[i-1]);
    //fp12cv_mul(&tmp_cv[i],&tmp_cv[i-1],&buf_cv);
  }
  for(i=0;i<12;i++){
    fp_set(&matrix_fp12cv.x[i][0],&tmp_cv[i].x0.x0);
    fp_set(&matrix_fp12cv.x[i][1],&tmp_cv[i].x0.x1);
    fp_set(&matrix_fp12cv.x[i][2],&tmp_cv[i].x0.x2);
    fp_set(&matrix_fp12cv.x[i][3],&tmp_cv[i].x0.x3);
    fp_set(&matrix_fp12cv.x[i][4],&tmp_cv[i].x1.x0);
    fp_set(&matrix_fp12cv.x[i][5],&tmp_cv[i].x1.x1);
    fp_set(&matrix_fp12cv.x[i][6],&tmp_cv[i].x1.x2);
    fp_set(&matrix_fp12cv.x[i][7],&tmp_cv[i].x1.x3);
    fp_set(&matrix_fp12cv.x[i][8],&tmp_cv[i].x2.x0);
    fp_set(&matrix_fp12cv.x[i][9],&tmp_cv[i].x2.x1);
    fp_set(&matrix_fp12cv.x[i][10],&tmp_cv[i].x2.x2);
    fp_set(&matrix_fp12cv.x[i][11],&tmp_cv[i].x2.x3);
  }
  
  matrix_inv(&matrix_fp12cv_inv,&matrix_fp12cv);
  matrix_mul(&a1,&matrix_fp12cv,&matrix_fp12cv_inv);
  matrix_mul(&a2,&matrix_fp12cv_inv,&matrix_fp12cv);
  printf("matrix_fp12cv cmp = %d\n",matrix_cmp(&a1,&a2));

  

  matrix_mul(&matrix_of_fp12_to_fp12cv,&matrix_fp12_inv,&matrix_fp12cv);
  //matrix_mul(&matrix_of_fp12_to_fp12cv,&matrix_fp12cv,&matrix_fp12_inv);
  matrix_inv(&matrix_of_fp12cv_to_fp12,&matrix_of_fp12_to_fp12cv);
 
  
  matrix_mul(&a1,&matrix_of_fp12_to_fp12cv,&matrix_of_fp12cv_to_fp12);
  //matrix_printf("a1 = ",&a1);
  printf("m1 * m2 = e cmp = %d\n",matrix_cmp(&a1,&e));
  

  mpz_clear(exp);
}
*/
void matrix_build_fp12_and_fpm2(){
  printf("\n------build matrix fp12 and fpm2------\n");
  int i;
  fp12_t a[12];
  fpm2_t b[12];
  matrix_t ma,ma_inv;
  matrix_t mb,mb_inv;
  mpz_t exp;
  matrix_t x1,x2;

  mpz_init(exp);
  mpz_pow_ui(exp,prime_z,12);
  mpz_sub_ui(exp,exp,1);
  mpz_tdiv_q_ui(exp,exp,13);//exp = (p^12-1)/13

  //karatsuba
  while(1){
    fp12_set_random(&a[0],state);
    fp12_pow(&a[0],&a[0],exp);
    if(fp12_cmp_one(&a[0])!=0){
      break;
    }
  }
  
  for(i=1;i<12;i++){
    fp12_frobenius_map_p1(&a[i],&a[i-1]);
  }

/*
  for(i=0;i<12;i++){
    fp_set(&ma.x[i][0],&a[i].x0.x0.x0);
    fp_set(&ma.x[i][1],&a[i].x0.x0.x1);
    fp_set(&ma.x[i][4],&a[i].x0.x1.x0);
    fp_set(&ma.x[i][5],&a[i].x0.x1.x1);
    fp_set(&ma.x[i][8],&a[i].x0.x2.x0);
    fp_set(&ma.x[i][9],&a[i].x0.x2.x1);
    fp_set(&ma.x[i][2],&a[i].x1.x0.x0);
    fp_set(&ma.x[i][3],&a[i].x1.x0.x1);
    fp_set(&ma.x[i][6],&a[i].x1.x1.x0);
    fp_set(&ma.x[i][7],&a[i].x1.x1.x1);
    fp_set(&ma.x[i][10],&a[i].x1.x2.x0);
    fp_set(&ma.x[i][11],&a[i].x1.x2.x1);
  }
  */
  
  for(i=0;i<12;i++){
    fp_set(&ma.x[i][0],&a[i].x0.x0.x0);
    fp_set(&ma.x[i][1],&a[i].x0.x0.x1);
    fp_set(&ma.x[i][2],&a[i].x0.x1.x0);
    fp_set(&ma.x[i][3],&a[i].x0.x1.x1);
    fp_set(&ma.x[i][4],&a[i].x0.x2.x0);
    fp_set(&ma.x[i][5],&a[i].x0.x2.x1);
    fp_set(&ma.x[i][6],&a[i].x1.x0.x0);
    fp_set(&ma.x[i][7],&a[i].x1.x0.x1);
    fp_set(&ma.x[i][8],&a[i].x1.x1.x0);
    fp_set(&ma.x[i][9],&a[i].x1.x1.x1);
    fp_set(&ma.x[i][10],&a[i].x1.x2.x0);
    fp_set(&ma.x[i][11],&a[i].x1.x2.x1);
  }
  
  matrix_inv(&ma_inv,&ma);
  matrix_mul(&x1,&ma,&ma_inv);
  matrix_mul(&x2,&ma_inv,&ma);
  printf("a * a^-1 cmp a^-1 * a = %d\n",matrix_cmp(&x1,&x2));

  //fpm2
  while(1){
    fpm2_set_random(&b[0],state);
    fpm2_pow(&b[0],&b[0],exp);
    if(fpm2_cmp_ui(&b[0],1)!=0){
      break;
    }
  }
  for(i=1;i<12;i++){
    fpm2_frobenius(&b[i],&b[i-1]);
  }
  
  for(i=0;i<12;i++){
    fp_set(&mb.x[i][0],&b[i].x[0].x[0]);
    fp_set(&mb.x[i][1],&b[i].x[0].x[1]);
    fp_set(&mb.x[i][2],&b[i].x[0].x[2]);
    fp_set(&mb.x[i][3],&b[i].x[0].x[3]);
    fp_set(&mb.x[i][4],&b[i].x[1].x[0]);
    fp_set(&mb.x[i][5],&b[i].x[1].x[1]);
    fp_set(&mb.x[i][6],&b[i].x[1].x[2]);
    fp_set(&mb.x[i][7],&b[i].x[1].x[3]);
    fp_set(&mb.x[i][8],&b[i].x[2].x[0]);
    fp_set(&mb.x[i][9],&b[i].x[2].x[1]);
    fp_set(&mb.x[i][10],&b[i].x[2].x[2]);
    fp_set(&mb.x[i][11],&b[i].x[2].x[3]);
  }
  
  matrix_inv(&mb_inv,&mb);
  matrix_mul(&x1,&mb,&mb_inv);
  matrix_mul(&x2,&mb_inv,&mb);
  printf("b * b^-1 cmp b^-1 * b = %d\n",matrix_cmp(&x1,&x2));
  
  matrix_init(&matrix_of_fp12_to_fpm2);
  matrix_init(&matrix_of_fpm2_to_fp12);
  matrix_mul(&matrix_of_fp12_to_fpm2,&ma_inv,&mb);
  matrix_inv(&matrix_of_fpm2_to_fp12,&matrix_of_fp12_to_fpm2);
  //check
  matrix_mul(&x1,&matrix_of_fp12_to_fpm2,&matrix_of_fpm2_to_fp12);
  matrix_mul(&x1,&matrix_of_fpm2_to_fp12,&matrix_of_fp12_to_fpm2);
  printf("M *M^-1 cmp M^-1 * M = %d\n",matrix_cmp(&x1,&x2));
}
void matrix_build_fp12_and_fp12cv(){
  printf("\n------build matrix fp12 and fp12cv------\n");
  int i,j;
  fp12_t a[12];
  fp12cv_t b[12];
  matrix_t ma,ma_inv;
  matrix_t mb,mb_inv;
  mpz_t exp;
  matrix_t x1,x2;

  mpz_init(exp);
  mpz_pow_ui(exp,prime_z,12);
  mpz_sub_ui(exp,exp,1);
  mpz_tdiv_q_ui(exp,exp,13);//exp = (p^12-1)/13

  //karatsuba
  while(1){
    fp12_set_random(&a[0],state);
    fp12_pow(&a[0],&a[0],exp);
    if(fp12_cmp_one(&a[0])!=0){
      break;
    }
  }
  
  for(i=1;i<12;i++){
    fp12_frobenius_map_p1(&a[i],&a[i-1]);
  }

  for(i=0;i<12;i++){
    fp_set(&ma.x[i][0],&a[i].x0.x0.x0);
    fp_set(&ma.x[i][1],&a[i].x0.x0.x1);
    fp_set(&ma.x[i][2],&a[i].x0.x1.x0);
    fp_set(&ma.x[i][3],&a[i].x0.x1.x1);
    fp_set(&ma.x[i][4],&a[i].x0.x2.x0);
    fp_set(&ma.x[i][5],&a[i].x0.x2.x1);
    fp_set(&ma.x[i][6],&a[i].x1.x0.x0);
    fp_set(&ma.x[i][7],&a[i].x1.x0.x1);
    fp_set(&ma.x[i][8],&a[i].x1.x1.x0);
    fp_set(&ma.x[i][9],&a[i].x1.x1.x1);
    fp_set(&ma.x[i][10],&a[i].x1.x2.x0);
    fp_set(&ma.x[i][11],&a[i].x1.x2.x1);
  }
  
  matrix_inv(&ma_inv,&ma);
  matrix_mul(&x1,&ma,&ma_inv);
  matrix_mul(&x2,&ma_inv,&ma);
  printf("a * a^-1 cmp a^-1 * a = %d\n",matrix_cmp(&x1,&x2));

  //fp12cv
  while(1){
    fp12cv_set_random(&b[0],state);
    fp12cv_pow(&b[0],&b[0],exp);
    if(fp12cv_cmp_ui(&b[0],1)!=0){
      break;
    }
  }
  for(i=1;i<12;i++){
    fp12cv_frobenius(&b[i],&b[i-1]);
  }
  for(i=0;i<12;i++){
    fp_set(&mb.x[i][0],&b[i].x0.x0);
    fp_set(&mb.x[i][1],&b[i].x0.x1);
    fp_set(&mb.x[i][2],&b[i].x0.x2);
    fp_set(&mb.x[i][3],&b[i].x0.x3);
    fp_set(&mb.x[i][4],&b[i].x1.x0);
    fp_set(&mb.x[i][5],&b[i].x1.x1);
    fp_set(&mb.x[i][6],&b[i].x1.x2);
    fp_set(&mb.x[i][7],&b[i].x1.x3);
    fp_set(&mb.x[i][8],&b[i].x2.x0);
    fp_set(&mb.x[i][9],&b[i].x2.x1);
    fp_set(&mb.x[i][10],&b[i].x2.x2);
    fp_set(&mb.x[i][11],&b[i].x2.x3);
  }
  matrix_inv(&mb_inv,&mb);
  matrix_mul(&x1,&mb,&mb_inv);
  matrix_mul(&x2,&mb_inv,&mb);
  printf("b * b^-1 cmp b^-1 * b = %d\n",matrix_cmp(&x1,&x2));
 
  matrix_mul(&matrix_of_fp12_to_fp12cv,&ma_inv,&mb);
  matrix_mul(&matrix_of_fp12cv_to_fp12,&mb_inv,&ma);
  //check
  matrix_mul(&x1,&matrix_of_fp12_to_fp12cv,&matrix_of_fp12cv_to_fp12);
  matrix_mul(&x2,&matrix_of_fp12cv_to_fp12,&matrix_of_fp12_to_fp12cv);
  printf("M *M^-1 cmp M^-1 * M = %d\n",matrix_cmp(&x1,&x2));
}
void fp12_to_fp12cv(fp12cv_t *b, fp12_t *a){
  int i;
  fp_t tmp1;
  fp_t ans[12];

  for(i=0;i<12;i++){
    fp_mul(&ans[i],&a->x0.x0.x0,&matrix_of_fp12_to_fp12cv.x[0][i]);
    fp_mul(&tmp1,&a->x0.x0.x1,&matrix_of_fp12_to_fp12cv.x[1][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x1.x0,&matrix_of_fp12_to_fp12cv.x[2][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x1.x1,&matrix_of_fp12_to_fp12cv.x[3][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x2.x0,&matrix_of_fp12_to_fp12cv.x[4][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x2.x1,&matrix_of_fp12_to_fp12cv.x[5][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x0.x0,&matrix_of_fp12_to_fp12cv.x[6][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x0.x1,&matrix_of_fp12_to_fp12cv.x[7][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x1.x0,&matrix_of_fp12_to_fp12cv.x[8][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x1.x1,&matrix_of_fp12_to_fp12cv.x[9][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x2.x0,&matrix_of_fp12_to_fp12cv.x[10][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x2.x1,&matrix_of_fp12_to_fp12cv.x[11][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
  }
  fp_set(&b->x0.x0,&ans[0]);
  fp_set(&b->x0.x1,&ans[1]);
  fp_set(&b->x0.x2,&ans[2]);
  fp_set(&b->x0.x3,&ans[3]);
  fp_set(&b->x1.x0,&ans[4]);
  fp_set(&b->x1.x1,&ans[5]);
  fp_set(&b->x1.x2,&ans[6]);
  fp_set(&b->x1.x3,&ans[7]);
  fp_set(&b->x2.x0,&ans[8]);
  fp_set(&b->x2.x1,&ans[9]);
  fp_set(&b->x2.x2,&ans[10]);
  fp_set(&b->x2.x3,&ans[11]);
}
void fp12cv_to_fp12(fp12_t *b ,fp12cv_t *a){
  int i;
  fp_t tmp1;
  fp_t ans[12];
  for(i=0;i<12;i++){
    fp_mul(&ans[i],&a->x0.x0,&matrix_of_fp12cv_to_fp12.x[0][i]);
    fp_mul(&tmp1,&a->x0.x1,&matrix_of_fp12cv_to_fp12.x[1][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x2,&matrix_of_fp12cv_to_fp12.x[2][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x3,&matrix_of_fp12cv_to_fp12.x[3][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x0,&matrix_of_fp12cv_to_fp12.x[4][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x1,&matrix_of_fp12cv_to_fp12.x[5][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x2,&matrix_of_fp12cv_to_fp12.x[6][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x3,&matrix_of_fp12cv_to_fp12.x[7][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x2.x0,&matrix_of_fp12cv_to_fp12.x[8][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x2.x1,&matrix_of_fp12cv_to_fp12.x[9][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x2.x2,&matrix_of_fp12cv_to_fp12.x[10][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x2.x3,&matrix_of_fp12cv_to_fp12.x[11][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
  }
  fp_set(&b->x0.x0.x0,&ans[0]);
  fp_set(&b->x0.x0.x1,&ans[1]);
  fp_set(&b->x0.x1.x0,&ans[2]);
  fp_set(&b->x0.x1.x1,&ans[3]);
  fp_set(&b->x0.x2.x0,&ans[4]);
  fp_set(&b->x0.x2.x1,&ans[5]);
  fp_set(&b->x1.x0.x0,&ans[6]);
  fp_set(&b->x1.x0.x1,&ans[7]);
  fp_set(&b->x1.x1.x0,&ans[8]);
  fp_set(&b->x1.x1.x1,&ans[9]);
  fp_set(&b->x1.x2.x0,&ans[10]);
  fp_set(&b->x1.x2.x1,&ans[11]);
}
void fp12_to_fpm2(fpm2_t *b, fp12_t *a){
  int i;
  fp_t tmp1;
  fp_t ans[12];

  for(i=0;i<12;i++){
    fp_mul(&ans[i],&a->x0.x0.x0,&matrix_of_fp12_to_fpm2.x[0][i]);
    fp_mul(&tmp1,&a->x0.x0.x1,&matrix_of_fp12_to_fpm2.x[1][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x1.x0,&matrix_of_fp12_to_fpm2.x[2][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x1.x1,&matrix_of_fp12_to_fpm2.x[3][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x2.x0,&matrix_of_fp12_to_fpm2.x[4][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x2.x1,&matrix_of_fp12_to_fpm2.x[5][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x0.x0,&matrix_of_fp12_to_fpm2.x[6][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x0.x1,&matrix_of_fp12_to_fpm2.x[7][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x1.x0,&matrix_of_fp12_to_fpm2.x[8][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x1.x1,&matrix_of_fp12_to_fpm2.x[9][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x2.x0,&matrix_of_fp12_to_fpm2.x[10][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x2.x1,&matrix_of_fp12_to_fpm2.x[11][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
  }
  fp_set(&b->x[0].x[0],&ans[0]);
  fp_set(&b->x[0].x[1],&ans[1]);
  fp_set(&b->x[0].x[2],&ans[2]);
  fp_set(&b->x[0].x[3],&ans[3]);
  fp_set(&b->x[1].x[0],&ans[4]);
  fp_set(&b->x[1].x[1],&ans[5]);
  fp_set(&b->x[1].x[2],&ans[6]);
  fp_set(&b->x[1].x[3],&ans[7]);
  fp_set(&b->x[2].x[0],&ans[8]);
  fp_set(&b->x[2].x[1],&ans[9]);
  fp_set(&b->x[2].x[2],&ans[10]);
  fp_set(&b->x[2].x[3],&ans[11]);
}
void fpm2_to_fp12(fp12_t *b, fpm2_t *a){
  int i;
  fp_t tmp1;
  fp_t ans[12];
  for(i=0;i<12;i++){
    fp_mul(&ans[i],&a->x[0].x[0],&matrix_of_fpm2_to_fp12.x[0][i]);
    fp_mul(&tmp1,&a->x[0].x[1],&matrix_of_fpm2_to_fp12.x[1][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[0].x[2],&matrix_of_fpm2_to_fp12.x[2][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[0].x[3],&matrix_of_fpm2_to_fp12.x[3][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[1].x[0],&matrix_of_fpm2_to_fp12.x[4][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[1].x[1],&matrix_of_fpm2_to_fp12.x[5][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[1].x[2],&matrix_of_fpm2_to_fp12.x[6][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[1].x[3],&matrix_of_fpm2_to_fp12.x[7][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[2].x[0],&matrix_of_fpm2_to_fp12.x[8][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[2].x[1],&matrix_of_fpm2_to_fp12.x[9][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[2].x[2],&matrix_of_fpm2_to_fp12.x[10][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[2].x[3],&matrix_of_fpm2_to_fp12.x[11][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
  }
  fp_set(&b->x0.x0.x0,&ans[0]);
  fp_set(&b->x0.x0.x1,&ans[1]);
  fp_set(&b->x0.x1.x0,&ans[2]);
  fp_set(&b->x0.x1.x1,&ans[3]);
  fp_set(&b->x0.x2.x0,&ans[4]);
  fp_set(&b->x0.x2.x1,&ans[5]);
  fp_set(&b->x1.x0.x0,&ans[6]);
  fp_set(&b->x1.x0.x1,&ans[7]);
  fp_set(&b->x1.x1.x0,&ans[8]);
  fp_set(&b->x1.x1.x1,&ans[9]);
  fp_set(&b->x1.x2.x0,&ans[10]);
  fp_set(&b->x1.x2.x1,&ans[11]);
}

/*----------test----------*/
void matrix_build_fp12_and_fpm(){
  int i,j;
  fp12_t a[12];
  fpm_t b[12];
  matrix_t ma,ma_inv;
  matrix_t mb,mb_inv;
  mpz_t exp;
  matrix_t x1,x2;

  //exp
  mpz_init(exp);
  mpz_pow_ui(exp,prime_z,12);
  mpz_sub_ui(exp,exp,1);
  mpz_tdiv_q_ui(exp,exp,13);//exp = (p^12-1)/13

  //karatsuba
  while(1){
    fp12_set_random(&a[0],state);
    fp12_pow(&a[0],&a[0],exp);
    if(fp12_cmp_one(&a[0])!=0){
      break;
    }
  }
  for(i=1;i<12;i++){
    fp12_frobenius_map_p1(&a[i],&a[i-1]);
  }
  for(i=0;i<12;i++){
    fp_set(&ma.x[i][0],&a[i].x0.x0.x0);
    fp_set(&ma.x[i][1],&a[i].x0.x0.x1);
    fp_set(&ma.x[i][2],&a[i].x0.x1.x0);
    fp_set(&ma.x[i][3],&a[i].x0.x1.x1);
    fp_set(&ma.x[i][4],&a[i].x0.x2.x0);
    fp_set(&ma.x[i][5],&a[i].x0.x2.x1);
    fp_set(&ma.x[i][6],&a[i].x1.x0.x0);
    fp_set(&ma.x[i][7],&a[i].x1.x0.x1);
    fp_set(&ma.x[i][8],&a[i].x1.x1.x0);
    fp_set(&ma.x[i][9],&a[i].x1.x1.x1);
    fp_set(&ma.x[i][10],&a[i].x1.x2.x0);
    fp_set(&ma.x[i][11],&a[i].x1.x2.x1);
  }
  matrix_inv(&ma_inv,&ma);
  matrix_mul(&x1,&ma,&ma_inv);
  matrix_mul(&x2,&ma_inv,&ma);
  printf("a * a^-1 cmp a^-1 * a = %d\n",matrix_cmp(&x1,&x2));

  //fpm
  while(1){
    fpm_set_random(&b[0],state);
    fpm_pow(&b[0],&b[0],exp);
    if(fpm_cmp_ui(&b[0],1)!=0){
      break;
    }
  }
  for(i=1;i<12;i++){
    fpm_frobenius(&b[i],&b[i-1]);
  }
  for(i=0;i<12;i++){
    for(j=0;j<12;j++){
      fp_set(&mb.x[i][j],&b[i].x[j]);
    }
  }
  matrix_inv(&mb_inv,&mb);
  //set
  matrix_mul(&matrix_of_fp12_to_fpm,&ma_inv,&mb);
  matrix_inv(&matrix_of_fpm_to_fp12,&matrix_of_fp12_to_fpm);
}
void fp12_to_fpm(fpm_t *b, fp12_t *a){
  int i,j;
  fp_t tmp1;
  fp_t ans[12];

  for(i=0;i<12;i++){
    fp_mul(&ans[i],&a->x0.x0.x0,&matrix_of_fp12_to_fpm.x[0][i]);
    fp_mul(&tmp1,&a->x0.x0.x1,&matrix_of_fp12_to_fpm.x[1][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x1.x0,&matrix_of_fp12_to_fpm.x[2][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x1.x1,&matrix_of_fp12_to_fpm.x[3][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x2.x0,&matrix_of_fp12_to_fpm.x[4][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x2.x1,&matrix_of_fp12_to_fpm.x[5][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x0.x0,&matrix_of_fp12_to_fpm.x[6][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x0.x1,&matrix_of_fp12_to_fpm.x[7][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x1.x0,&matrix_of_fp12_to_fpm.x[8][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x1.x1,&matrix_of_fp12_to_fpm.x[9][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x2.x0,&matrix_of_fp12_to_fpm.x[10][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x2.x1,&matrix_of_fp12_to_fpm.x[11][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
  }
  for(i=0;i<12;i++){
    fp_set(&b->x[i],&ans[i]);
  }
}
void fpm_to_fp12(fp12_t *b, fpm_t *a){
  int i;
  fp_t tmp1;
  fp_t ans[12];
  for(i=0;i<12;i++){
    fp_mul(&ans[i],&a->x[0],&matrix_of_fpm_to_fp12.x[0][i]);
    fp_mul(&tmp1,&a->x[1],&matrix_of_fpm_to_fp12.x[1][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[2],&matrix_of_fpm_to_fp12.x[2][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[3],&matrix_of_fpm_to_fp12.x[3][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[4],&matrix_of_fpm_to_fp12.x[4][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[5],&matrix_of_fpm_to_fp12.x[5][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[6],&matrix_of_fpm_to_fp12.x[6][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[7],&matrix_of_fpm_to_fp12.x[7][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[8],&matrix_of_fpm_to_fp12.x[8][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[9],&matrix_of_fpm_to_fp12.x[9][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[10],&matrix_of_fpm_to_fp12.x[10][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[11],&matrix_of_fpm_to_fp12.x[11][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
  }
  fp_set(&b->x0.x0.x0,&ans[0]);
  fp_set(&b->x0.x0.x1,&ans[1]);
  fp_set(&b->x0.x1.x0,&ans[2]);
  fp_set(&b->x0.x1.x1,&ans[3]);
  fp_set(&b->x0.x2.x0,&ans[4]);
  fp_set(&b->x0.x2.x1,&ans[5]);
  fp_set(&b->x1.x0.x0,&ans[6]);
  fp_set(&b->x1.x0.x1,&ans[7]);
  fp_set(&b->x1.x1.x0,&ans[8]);
  fp_set(&b->x1.x1.x1,&ans[9]);
  fp_set(&b->x1.x2.x0,&ans[10]);
  fp_set(&b->x1.x2.x1,&ans[11]);
}
/*----------test----------*/
void matrix_build_fpm_and_fp12cv(){
  int i,j;
  
  fpm_t a[12];
  fp12cv_t b[12];
  matrix_t ma,ma_inv;
  matrix_t mb,mb_inv;
  mpz_t exp;
  matrix_t x1,x2;

  //exp
  mpz_init(exp);
  mpz_pow_ui(exp,prime_z,12);
  mpz_sub_ui(exp,exp,1);
  mpz_tdiv_q_ui(exp,exp,13);//exp = (p^12-1)/13

  //fpm
  while(1){
    fpm_set_random(&a[0],state);
    fpm_pow(&a[0],&a[0],exp);
    if(fpm_cmp_ui(&a[0],1)!=0){
      break;
    }
  }
  for(i=1;i<12;i++){
    fpm_frobenius(&a[i],&a[i-1]);
  }
  for(i=0;i<12;i++){
    for(j=0;j<12;j++){
      fp_set(&ma.x[i][j],&a[i].x[j]);
    }
  }
  matrix_inv(&ma_inv,&ma);
  matrix_mul(&x1,&ma,&ma_inv);
  matrix_mul(&x2,&ma_inv,&ma);
  printf("a * a^-1 cmp a^-1 * a = %d\n",matrix_cmp(&x1,&x2));

  //fp12cv
  while(1){
    fp12cv_set_random(&b[0],state);
    fp12cv_pow(&b[0],&b[0],exp);
    if(fp12cv_cmp_one(&b[0])!=0){
      break;
    }
  }
  for(i=1;i<12;i++){
    fp12cv_frobenius(&b[i],&b[i-1]);
  }
  for(i=0;i<12;i++){
    fp_set(&mb.x[i][0],&b[i].x0.x0);
    fp_set(&mb.x[i][1],&b[i].x0.x1);
    fp_set(&mb.x[i][2],&b[i].x0.x2);
    fp_set(&mb.x[i][3],&b[i].x0.x3);
    fp_set(&mb.x[i][4],&b[i].x1.x0);
    fp_set(&mb.x[i][5],&b[i].x1.x1);
    fp_set(&mb.x[i][6],&b[i].x1.x2);
    fp_set(&mb.x[i][7],&b[i].x1.x3);
    fp_set(&mb.x[i][8],&b[i].x2.x0);
    fp_set(&mb.x[i][9],&b[i].x2.x1);
    fp_set(&mb.x[i][10],&b[i].x2.x2);
    fp_set(&mb.x[i][11],&b[i].x2.x3);
  }
  matrix_inv(&mb_inv,&mb);
  matrix_mul(&x1,&mb,&mb_inv);
  matrix_mul(&x2,&mb_inv,&mb);
  printf("b * b^-1 cmp b^-1 * b = %d\n",matrix_cmp(&x1,&x2));

  
  //set
  /*
  matrix_mul(&matrix_of_fpm_to_fp12cv,&ma_inv,&mb);
  matrix_inv(&matrix_of_fp12cv_to_fpm,&matrix_of_fpm_to_fp12cv);
  */
  matrix_mul(&matrix_of_fp12cv_to_fpm,&mb_inv,&ma);
  matrix_inv(&matrix_of_fpm_to_fp12cv,&matrix_of_fp12cv_to_fpm);
  matrix_mul(&x1,&matrix_of_fpm_to_fp12cv,&matrix_of_fp12cv_to_fpm);
  matrix_mul(&x2,&matrix_of_fp12cv_to_fpm,&matrix_of_fpm_to_fp12cv);
  printf("M2 * M2^-1 cmp M2^-1 * M2 = %d\n",matrix_cmp(&x1,&x2));
  
}
void fpm_to_fp12cv(fp12cv_t *b, fpm_t *a){
  int i;
  fp_t tmp1;
  fp_t ans[12];
  for(i=0;i<12;i++){
    fp_mul(&ans[i],&a->x[0],&matrix_of_fpm_to_fp12cv.x[0][i]);
    fp_mul(&tmp1,&a->x[1],&matrix_of_fpm_to_fp12cv.x[1][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[2],&matrix_of_fpm_to_fp12cv.x[2][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[3],&matrix_of_fpm_to_fp12cv.x[3][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[4],&matrix_of_fpm_to_fp12cv.x[4][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[5],&matrix_of_fpm_to_fp12cv.x[5][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[6],&matrix_of_fpm_to_fp12cv.x[6][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[7],&matrix_of_fpm_to_fp12cv.x[7][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[8],&matrix_of_fpm_to_fp12cv.x[8][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[9],&matrix_of_fpm_to_fp12cv.x[9][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[10],&matrix_of_fpm_to_fp12cv.x[10][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x[11],&matrix_of_fpm_to_fp12cv.x[11][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
  }
  fp_set(&b->x0.x0,&ans[0]);
  fp_set(&b->x0.x1,&ans[1]);
  fp_set(&b->x0.x2,&ans[2]);
  fp_set(&b->x0.x3,&ans[3]);
  fp_set(&b->x1.x0,&ans[4]);
  fp_set(&b->x1.x1,&ans[5]);
  fp_set(&b->x1.x2,&ans[6]);
  fp_set(&b->x1.x3,&ans[7]);
  fp_set(&b->x2.x0,&ans[8]);
  fp_set(&b->x2.x1,&ans[9]);
  fp_set(&b->x2.x2,&ans[10]);
  fp_set(&b->x2.x3,&ans[11]);
}
void fp12cv_to_fpm(fpm_t *b, fp12cv_t *a){
  int i,j;
  fp_t tmp1;
  fp_t ans[12];

  for(i=0;i<12;i++){
    fp_mul(&ans[i],&a->x0.x0,&matrix_of_fp12cv_to_fpm.x[0][i]);
    fp_mul(&tmp1,&a->x0.x1,&matrix_of_fp12cv_to_fpm.x[1][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x2,&matrix_of_fp12cv_to_fpm.x[2][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x0.x3,&matrix_of_fp12cv_to_fpm.x[3][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x0,&matrix_of_fp12cv_to_fpm.x[4][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x1,&matrix_of_fp12cv_to_fpm.x[5][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x2,&matrix_of_fp12cv_to_fpm.x[6][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x1.x3,&matrix_of_fp12cv_to_fpm.x[7][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x2.x0,&matrix_of_fp12cv_to_fpm.x[8][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x2.x1,&matrix_of_fp12cv_to_fpm.x[9][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x2.x2,&matrix_of_fp12cv_to_fpm.x[10][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
    fp_mul(&tmp1,&a->x2.x3,&matrix_of_fp12cv_to_fpm.x[11][i]);
    fp_add(&ans[i],&ans[i],&tmp1);
  }
  for(i=0;i<12;i++){
    fp_set(&b->x[i],&ans[i]);
  }
}
void matrix_build(){
  matrix_build_fp12_and_fpm();
  matrix_build_fpm_and_fp12cv();
  matrix_mul(&matrix_of_fp12_to_fp12cv,&matrix_of_fp12_to_fpm,&matrix_of_fpm_to_fp12cv);
  //matrix_mul(&matrix_of_fp12_to_fp12cv,&matrix_of_fpm_to_fp12cv,&matrix_of_fp12_to_fpm);
  matrix_mul(&matrix_of_fp12cv_to_fp12,&matrix_of_fp12cv_to_fpm,&matrix_of_fpm_to_fp12);
  //matrix_inv(&matrix_of_fp12cv_to_fp12,&matrix_of_fp12_to_fp12cv);
}