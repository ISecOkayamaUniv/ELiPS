#include "ELiPS/bls12.h"
#include "ELiPS/fr.h"

int main(){
    bls12_init();
    fr_order_init();
    debug_pairing(100);
    
    return 0;
}