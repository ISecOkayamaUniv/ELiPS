#include <stdio.h>
#include "ELiPS/bls12.h"
#include "ELiPS/file.h"


//test function
void test_write_efp12_single(const char *file){
    printf("\nwrite : %s\n", file);
    FILE *file_p;
    efp12_t A;
    efp12_init(&A);
    //efp12_rational_point(&A);
    bls12_generate_g2(&A);
    efp12_println("A=", &A);

    file_p = fopen(file, "wb");
    if (file_p == NULL){
        printf("cannot open\n");
        return;
    }

    efp12_write(file_p, &A);
    if (fclose(file_p) == EOF) {
        return;
    }
    printf("success to save\n");
}

void test_read_efp12_single(const char *file){
    printf("\nread : %s\n", file);
    FILE *file_p;
    efp12_t A;
    efp12_init(&A);

    file_p = fopen(file, "rb");
    if (file_p == NULL){
        printf("cannot open\n");
        return;
    }

    efp12_read(&A, file_p);
    efp12_println("A=",&A);
    if (fclose(file_p) == EOF) {
        return;
    }
    printf("success to read\n");
    return;

}

void test_write_efp12_10(const char *file){
    printf("\nwrite : %s\n", file);
    FILE *file_p;
    efp12_t A[10];
    int i;
    for (i = 0; i < 10;i++){
        efp12_init(&A[i]);
        bls12_generate_g2(&A[i]);
        printf("A[%d]", i);
        efp12_println("=", &A[i]);
    }
    file_p = fopen(file, "wb");
    if (file_p == NULL){
        printf("cannot open\n");
        return;
    }
    for (i = 0; i < 10;i++){
        efp12_write(file_p, &A[i]);
    }
    if (fclose(file_p) == EOF) {
        return;
    }
    printf("success to save\n");
}

void test_read_efp12_10(const char *file){
    printf("\nread : %s\n", file);
    FILE *file_p;
    efp12_t A[10];
    int i;
    file_p = fopen(file, "rb");
    for (i = 0; i < 10;i++){
        efp12_init(&A[i]);
    }

    if (file_p == NULL){
        printf("cannot open\n");
        return;
    }
    for (i = 0; i < 10;i++){
        efp12_read(&A[i], file_p);
        printf("A[%d]", i);
        efp12_println("=", &A[i]);
    }
    if (fclose(file_p) == EOF) {
        return;
    }

    printf("show\n");
}

int main(){
    bls12_init();
    
    //file name
    const char *file = "testpc_txt.txt";
    //const char *file = "testraspi_txt.txt";

    test_write_efp12_single(file);
    test_read_efp12_single(file);
    //test_write_efp12_10(file);
    //test_read_efp12_10(file);
}