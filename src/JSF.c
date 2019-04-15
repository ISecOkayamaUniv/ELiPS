#include <ELiPS/JSF.h>
//JSF
void Joint_sparse_form(int **binary,mpz_t scalar[2],int *loop_length){
	int i,j;
	unsigned long int u;
	mpz_t mod_2,mod_4,mod_8;
	mpz_init(mod_2);
	mpz_init(mod_4);
	mpz_init(mod_8);
	
	mpz_t k[2];
	mpz_init(k[0]);
	mpz_init(k[1]);
	//set
	j=0;
	mpz_set(k[0],scalar[0]);
	mpz_set(k[1],scalar[1]);
	
	while(mpz_cmp_ui(k[0],0)>0 || mpz_cmp_ui(k[1],0)>0){
		for(i=0; i<2; i++){
			mpz_mod_ui(mod_2,k[i],2);
			if(mpz_cmp_ui(mod_2,0)==0){
				u=0;
			}else{
				mpz_mod_ui(mod_4,k[i],4);
				u=mpz_get_ui(mod_4);
				if(u==3){
					u=-1;
				}
				mpz_mod_ui(mod_8,k[i],8);
				mpz_mod_ui(mod_4,k[1-i],4);
				if((mpz_cmp_ui(mod_8,3)==0 || mpz_cmp_ui(mod_8,5)==0) && mpz_cmp_ui(mod_4,2)==0){
					u=-u;
				}
			}
			binary[i][j]=u;
		}
		for(i=0; i<2; i++){
			u=binary[i][j];
			switch (u){
				case 1:
					mpz_sub_ui(k[i],k[i],1);
					break;
				case -1:
					mpz_add_ui(k[i],k[i],1);
					break;
				default:
					break;
			}
			mpz_tdiv_q_ui(k[i],k[i],2);
		}
		j=j+1;
	}
	*loop_length=j-1;
	
	mpz_clear(mod_2);
	mpz_clear(mod_4);
	mpz_clear(mod_8);
	mpz_clear(k[0]);
	mpz_clear(k[1]);	
}
