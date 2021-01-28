#include "ELiPS/bls12.h"

/**
 * @brief Set mp_limb_t to string(hex) and return size of string.
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
int mpn_to_hexstr(unsigned char *ANS, mp_limb_t *A){
  mp_size_t a_size=-1;
  int str_size = 0,i=0;
  mp_limb_t tmp[FPLIMB];
  mpn_copyd(tmp,A,FPLIMB);
  for (int i = FPLIMB - 1; i >= 0; i--){
    if(tmp[i]==0) continue;
    a_size=i;
    break;
  }
  a_size++;
  str_size = mpn_get_str(ANS, 16, tmp, a_size);
  return str_size;
}

/**
 * @brief Set string(hex) to mp_limb_t.
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]P --a integer size of string.
 * @param[in]P --a pointer to set.
 */
void hexstr_to_mpn(mp_limb_t *ANS,int str_size, unsigned char *A){
  mpn_set_str(ANS,A,(mp_size_t)str_size,16);
}

void test_str_mpn_cast(){
  fp_t A_x32,B_x64;
  unsigned char a_char[FPLIMB_BITS/4];//Up to 461/4
  int str_size;

  fp_init(&A_x32);
  fp_set_random(&A_x32,state);
  //fp_set_ui(&A_x32, 257);
  fp_println("A_x32=", &A_x32);
  str_size=mpn_to_hexstr(a_char,A_x32.x0);

  /************************/
  //please send a_char and str_size to other device
  /************************/

  fp_init(&B_x64);
  hexstr_to_mpn(B_x64.x0,str_size,a_char);
  fp_println("B_x64=",&B_x64);

}

// void mpn_x32_to_x64(mp_limb_t *ANS_x64,mp_limb_t *A_x32){
//   int mplimb_x64 = 8, mplimb_x32 = 15, i=0;
//   for (i = 0; i < mplimb_x64;i++){

//   }
// }

// void mpn_x64_to_x32(mp_limb_t *ANS_x32,mp_limb_t *A_x64){
//   int mplimb_x64 = 8, mplimb_x32 = 15, i=0;
//   for (i = 0; i < mplimb_x64;i++){
//     ANS[i] = ;
//     ANS[i + 1] = ;
//   }
// }

// void test_mpn_x64_x32(){
//   /*************************************/
//   //PLEASE EXECUTE ON RASPI(32 bit OS)!!!
//   fp_t A_x32;
//   fp_set_random(&A_x32, state);
//   fp_println("A_x32=", &A_x32);
//   /*************************************/

//   /*************************************/
//   //PLEASE EXECUTE ON 64bit karnel!!!
//   /*************************************/
//   fp_t B_x64;
  
// }

int main(){
  bls12_init();
  bls12_print_parameters();
  test_str_mpn_cast();
  return 0;
}

