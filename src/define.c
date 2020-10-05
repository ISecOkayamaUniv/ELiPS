#include <ELiPS/define.h>

int bn12_X_binary[bn12_X_length+1];
int bn12_X6_2_binary[bn12_X6_2_length+1];

int bls12_X_binary[bls12_X_length+1];
int bls12_X2_binary[bls12_X2_length+1];


//int cost_add,cost_add_ui,cost_sub,cost_sub_ui,cost_mul,cost_mul_ui,cost_sqr,cost_inv,cost_mod;
int cost_add,cost_add_ui,cost_sub,cost_sub_ui,cost_mul,cost_set_neg,cost_sqr,cost_inv,cost_mod;
int cost_add_nonmod,cost_add_nonmod_double,cost_sub_nonmod,cost_sub_nonmod_double,cost_div2,cost_mod_nomal;


mp_limb_t buf[FPLIMB],tmp_mul[FPLIMB2],tmp1[FPLIMB],tmp2[FPLIMB];

/*============================================================================*/
/* Pairing functions                                                          */
/*============================================================================*/

gmp_randstate_t state;
mpz_t X_z,prime_z,order_z,trace_z;
mp_limb_t X,prime[FPLIMB],order[FPLIMB],trace[FPLIMB];
mp_limb_t prime2[FPLIMB2];
fp2_t Alpha_1,Alpha_1_inv;
mp_limb_t epsilon1[FPLIMB],epsilon2[FPLIMB];
mp_limb_t Two_inv[FPLIMB];
mpz_t Two_inv_z;
mpz_t root_2,root_X;
mpz_t efp_total,efp12_total;
fp2_t frobenius_constant[12][6];
fp2_t skew_frobenius_constant[12][2];
mp_limb_t curve_b[FPLIMB];

//montgomery
mp_limb_t R[FPLIMB],Ri[FPLIMB],R1[FPLIMB],RR[FPLIMB],Ni[FPLIMB];
int m;
mp_limb_t u[FPLIMB+1];
mp_limb_t N[FPLIMB2],R2[FPLIMB],R3[FPLIMB],RmodP[FPLIMB];
/*============================================================================*/
/* Test functions                                                             */
/*============================================================================*/
struct timeval tv_start,tv_end;

float MILLER_OPT_AFFINE,FINALEXP_OPT_AFFINE;
float MILLER_OPT_PROJECTIVE,FINALEXP_OPT_PROJECTIVE;

cost MILLER_OPT_AFFINE_COST,FINALEXP_OPT_AFFINE_COST;
cost MILLER_OPT_PROJECTIVE_COST,FINALEXP_OPT_PROJECTIVE_COST;
