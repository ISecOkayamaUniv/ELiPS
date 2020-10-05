#include <ELiPS/time.h>
/*----------------------------------------------------------------------------*/
//time
float timedifference_msec(struct timeval tv_start, struct timeval tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) * 1000.0f + (tv_end.tv_usec - tv_start.tv_usec) / 1000.0f;
}

float timedifference_usec(struct timeval tv_start, struct timeval tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) + (tv_end.tv_usec - tv_start.tv_usec);
}

/*
float timedifference_nsec(struct timespec tv_start, struct timespec tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) + (tv_end.tv_nsec - tv_start.tv_nsec);
}
*/

// void cost_zero(){
//     cost_add=0;
//     cost_add_ui=0;
//     cost_sub=0;
//     cost_sub_ui=0;
//     cost_mul=0;
//     cost_mul_ui=0;
//     cost_sqr=0;
//     cost_inv=0;
//     cost_mod=0;
// }
// void cost_init(cost *A){
// 	A->add=0;
// 	A->add_ui=0;
// 	A->sub=0;
// 	A->sub_ui=0i;
// 	A->mul=0;
// 	A->mul_ui=0;
// 	A->sqr=0;
// 	A->inv=0;
// 	A->mod=0;
// }

// void cost_check(cost *A){
// 	A->add=cost_add;
// 	A->add_ui=cost_add_ui;
// 	A->sub=cost_sub;
// 	A->sub_ui=cost_sub_ui;
// 	A->mul=cost_mul;
// 	A->mul_ui=cost_mul_ui;
// 	A->sqr=cost_sqr;
// 	A->inv=cost_inv;
// 	A->mod=cost_mod;
// }
// void cost_addition(cost *A,cost *B){
// 	A->add+=B->add;
// 	A->add_ui+=B->add_ui;
// 	A->sub+=B->sub;
// 	A->sub_ui+=B->sub_ui;
// 	A->mul+=B->mul;
// 	A->mul_ui+=B->mul_ui;
// 	A->sqr+=B->sqr;
// 	A->inv+=B->inv;
// 	A->mod+=B->mod;
// }
// void cost_substruction(cost *ANS,cost *A,cost *B){
// 	ANS->add=A->add-B->add;
// 	ANS->add_ui=A->add_ui-B->add_ui;
// 	ANS->sub=A->sub-B->sub;
// 	ANS->sub_ui=A->sub_ui-B->sub_ui;
// 	ANS->mul=A->mul-B->mul;
// 	ANS->mul_ui=A->mul_ui-B->mul_ui;
// 	ANS->sqr=A->sqr-B->sqr;
// 	ANS->inv=A->inv-B->inv;
// 	ANS->mod=A->mod-B->mod;
// }
// void cost_printf(char *str,cost *A,int n){
//     printf("COST<%s>\n",str);
//     printf("ADD   :%d\n",A->add/n);
//     printf("ADD UI:%d\n",A->add_ui/n);
//     printf("SUB   :%d\n",A->sub/n);
//     printf("SUB UI:%d\n",A->sub_ui/n);
//     printf("MUL   :%d\n",A->mul/n);
//     printf("MUL UI:%d\n",A->mul_ui/n);
//     printf("SQR   :%d\n",A->sqr/n);
//     printf("INV   :%d\n",A->inv/n);
//     printf("MOD   :%d\n",A->mod/n);
// }
void cost_zero(){
	cost_add=0;
	cost_add_ui=0;
	cost_add_nonmod=0;
	cost_add_nonmod_double=0;
	cost_sub=0;
	cost_sub_ui=0;
	cost_sub_nonmod=0;
	cost_sub_nonmod_double=0;
	cost_mul=0;
	cost_set_neg=0;
	cost_div2=0;
	cost_sqr=0;
	cost_inv=0;
	cost_mod=0;
	cost_mod_nomal=0;
}
void cost_init(cost *A){
	A->add=0;
	A->add_ui=0;
	A->add_nonmod=0;
	A->add_nonmod_double=0;
	A->sub=0;
	A->sub_ui=0;
	A->sub_nonmod=0;
	A->sub_nonmod_double=0;
	A->mul=0;
	A->set_neg=0;
	A->div2=0;
	A->sqr=0;
	A->inv=0;
	A->mod=0;
	A->mod_nomal=0;
}

void cost_check(cost *A){
	A->add=cost_add;
	A->add_ui=cost_add_ui;
	A->add_nonmod=cost_add_nonmod;
	A->add_nonmod_double=cost_add_nonmod_double;
	A->sub=cost_sub;
	A->sub_ui=cost_sub_ui;
	A->sub_nonmod=cost_sub_nonmod;
	A->sub_nonmod_double=cost_sub_nonmod_double;
	A->mul=cost_mul;
	A->set_neg=cost_set_neg;
	A->div2=cost_div2;
	A->sqr=cost_sqr;
	A->inv=cost_inv;
	A->mod=cost_mod;
	A->mod_nomal=cost_mod_nomal;
}
void cost_addition(cost *A,cost *B){
	A->add+=B->add;
	A->add_ui+=B->add_ui;
	A->add_nonmod+=B->add_nonmod;
	A->add_nonmod_double+=B->add_nonmod_double;
	A->sub+=B->sub;
	A->sub_ui+=B->sub_ui;
	A->sub_nonmod+=B->sub_nonmod;
	A->sub_nonmod_double+=B->sub_nonmod_double;
	A->mul+=B->mul;
	A->set_neg+=B->set_neg;
	A->div2+=B->div2;
	A->sqr+=B->sqr;
	A->inv+=B->inv;
	A->mod+=B->mod;
	A->mod_nomal+=B->mod_nomal;
}
void cost_substruction(cost *ANS,cost *A,cost *B){
	ANS->add=A->add-B->add;
	ANS->add_ui=A->add_ui-B->add_ui;
	ANS->add_nonmod=A->add_nonmod-B->add_nonmod;
	ANS->add_nonmod_double=A->add_nonmod_double-B->add_nonmod_double;
	ANS->sub=A->sub-B->sub;
	ANS->sub_ui=A->sub_ui-B->sub_ui;
	ANS->sub_nonmod=A->sub_nonmod-B->sub_nonmod;
	ANS->sub_nonmod_double=A->sub_nonmod_double-B->sub_nonmod_double;
	ANS->mul=A->mul-B->mul;
	ANS->set_neg=A->set_neg-B->set_neg;
	ANS->div2=A->div2-B->div2;
	ANS->sqr=A->sqr-B->sqr;
	ANS->inv=A->inv-B->inv;
	ANS->mod=A->mod-B->mod;
	ANS->mod_nomal=A->mod_nomal-B->mod_nomal;
}
void cost_printf(char *str,cost *A,int n){
	printf("COST<%s>\n",str);
	printf("ADD   :%d\n",A->add/n);
	printf("ADD UI:%d\n",A->add_ui/n);
	printf("ADD NM:%d\n",A->add_nonmod/n);
	printf("ADD2NM:%d\n",A->add_nonmod_double/n);
	printf("SUB   :%d\n",A->sub/n);
	printf("SUB UI:%d\n",A->sub_ui/n);
	printf("SUB NM:%d\n",A->sub_nonmod/n);
	printf("SUB2NM:%d\n",A->sub_nonmod_double/n);
	printf("MUL   :%d\n",A->mul/n);
	printf("SET NEG:%d\n",A->set_neg/n);
	printf("DIV2  :%d\n",A->div2/n);
	printf("SQR   :%d\n",A->sqr/n);
	printf("INV   :%d\n",A->inv/n);
	printf("MOD   :%d\n",A->mod/n);
	printf("MOD nomal:%d\n",A->mod_nomal/n);
}
