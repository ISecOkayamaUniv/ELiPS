#include <ELiPS/matrix.h>

void matrix_build_fp12_and_fp12cv(){
  fp_t m_fp12[12][12],m_fp12cv[12][12];
  fp12_t tmp[12],buf;
  fp12cv_t tmp_cv[12],buf_cv;
  mpz_t exp;
  mpz_init(exp);
  mpz_pow_ui(exp,prime_z,12);
  mpz_sub_ui(exp,exp,1);
  mpz_tdiv_q_ui(exp,exp,35);//exp = (p^12-1)/35

  //fp12
  fp12_set_random(&buf,state);
  fp12_pow(&tmp[0],&buf,exp);
  for(int i=1;i<12;i++){
    fp12_mul(&tmp[i],&tmp[i-1],&tmp[0]);
  }
  
  //set
  for(int i=0;i<12;i++){
    fp_set(&m_fp12[i][0],&tmp[i].x0.x0.x0);
    fp_set(&m_fp12[i][1],&tmp[i].x0.x0.x1);
    fp_set(&m_fp12[i][2],&tmp[i].x0.x1.x0);
    fp_set(&m_fp12[i][3],&tmp[i].x0.x1.x1);
    fp_set(&m_fp12[i][4],&tmp[i].x0.x2.x0);
    fp_set(&m_fp12[i][5],&tmp[i].x0.x2.x1);
    fp_set(&m_fp12[i][6],&tmp[i].x1.x0.x0);
    fp_set(&m_fp12[i][7],&tmp[i].x1.x0.x1);
    fp_set(&m_fp12[i][8],&tmp[i].x1.x1.x0);
    fp_set(&m_fp12[i][9],&tmp[i].x1.x1.x1);
    fp_set(&m_fp12[i][10],&tmp[i].x1.x2.x0);
    fp_set(&m_fp12[i][11],&tmp[i].x1.x2.x1);
  }

  //fp12cv
  fp12cv_set_random(&buf_cv,state);
  fp12cv_pow(&tmp_cv[0],&buf_cv,exp);
  for(int i=1;i<12;i++){
    fp12cv_mul(&tmp_cv[i],&tmp_cv[i-1],&tmp_cv[0]);
  }

  //set
  for(int i=0;i<12;i++){
    fp_set(&m_fp12cv[i][0],&tmp_cv[i].x0.x0);
    fp_set(&m_fp12cv[i][1],&tmp_cv[i].x0.x1);
    fp_set(&m_fp12cv[i][2],&tmp_cv[i].x0.x2);
    fp_set(&m_fp12cv[i][3],&tmp_cv[i].x0.x3);
    fp_set(&m_fp12cv[i][4],&tmp_cv[i].x1.x0);
    fp_set(&m_fp12cv[i][5],&tmp_cv[i].x1.x1);
    fp_set(&m_fp12cv[i][6],&tmp_cv[i].x1.x2);
    fp_set(&m_fp12cv[i][7],&tmp_cv[i].x1.x3);
    fp_set(&m_fp12cv[i][8],&tmp_cv[i].x2.x0);
    fp_set(&m_fp12cv[i][9],&tmp_cv[i].x2.x1);
    fp_set(&m_fp12cv[i][10],&tmp_cv[i].x2.x2);
    fp_set(&m_fp12cv[i][11],&tmp_cv[i].x2.x3);
  }




  mpz_clear(exp);
}