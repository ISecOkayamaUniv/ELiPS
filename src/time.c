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

void cost_zero(){
    cost_add=0;
    cost_add_ui=0;
    cost_sub=0;
    cost_sub_ui=0;
    cost_mul=0;
    cost_mul_ui=0;
    cost_sqr=0;
    cost_inv=0;
    cost_mod=0;
}
void cost_init(cost *A){
	A->add=0;
	A->add_ui=0;
	A->sub=0;
	A->sub_ui=0i;
	A->mul=0;
	A->mul_ui=0;
	A->sqr=0;
	A->inv=0;
	A->mod=0;
}

void cost_check(cost *A){
	A->add=cost_add;
	A->add_ui=cost_add_ui;
	A->sub=cost_sub;
	A->sub_ui=cost_sub_ui;
	A->mul=cost_mul;
	A->mul_ui=cost_mul_ui;
	A->sqr=cost_sqr;
	A->inv=cost_inv;
	A->mod=cost_mod;
}
void cost_addition(cost *A,cost *B){
	A->add+=B->add;
	A->add_ui+=B->add_ui;
	A->sub+=B->sub;
	A->sub_ui+=B->sub_ui;
	A->mul+=B->mul;
	A->mul_ui+=B->mul_ui;
	A->sqr+=B->sqr;
	A->inv+=B->inv;
	A->mod+=B->mod;
}
void cost_substruction(cost *ANS,cost *A,cost *B){
	ANS->add=A->add-B->add;
	ANS->add_ui=A->add_ui-B->add_ui;
	ANS->sub=A->sub-B->sub;
	ANS->sub_ui=A->sub_ui-B->sub_ui;
	ANS->mul=A->mul-B->mul;
	ANS->mul_ui=A->mul_ui-B->mul_ui;
	ANS->sqr=A->sqr-B->sqr;
	ANS->inv=A->inv-B->inv;
	ANS->mod=A->mod-B->mod;
}
void cost_printf(char *str,cost *A,int n){
    printf("COST<%s>\n",str);
    printf("ADD   :%d\n",A->add/n);
    printf("ADD UI:%d\n",A->add_ui/n);
    printf("SUB   :%d\n",A->sub/n);
    printf("SUB UI:%d\n",A->sub_ui/n);
    printf("MUL   :%d\n",A->mul/n);
    printf("MUL UI:%d\n",A->mul_ui/n);
    printf("SQR   :%d\n",A->sqr/n);
    printf("INV   :%d\n",A->inv/n);
    printf("MOD   :%d\n",A->mod/n);
}