#include <ELiPS/file.h>
void mpn_write(FILE *file,mp_limb_t *A, mp_size_t size){
#if ARCBIT == 64
    fwrite(A, sizeof(mp_limb_t)*size, 1, file);
#else
    mp_limb_t dumy[size];
    mpn_zero(dumy,size);
    mpn_copyd(dumy,A,size-1);
    fwrite(dumy,sizeof(mp_limb_t)*size,1,file);
#endif
}
void mpn_read(mp_limb_t *A,FILE *file,mp_size_t size){
#if ARCBIT == 64
    fread(A, sizeof(mp_limb_t) * size, 1, file);
#else
    mp_limb_t dumy[size];
    mpn_zero(dumy,size);
    fread(dumy, sizeof(mp_limb_t) * size, 1, file);
    mpn_copyd(A,dumy,size-1);
#endif  
}
void fp_write(FILE *file,fp_t *A){
    mpn_write(file, A->x0,FPLIMB);
}
void fp_read(fp_t *A,FILE *file){
    mpn_read(A->x0, file,FPLIMB);
}

void fp2_write(FILE *file,fp2_t *A){
    fp_write(file, &A->x0);
    fp_write(file, &A->x1);
}
void fp2_read(fp2_t *A,FILE *file){
    fp_read(&A->x0, file);
    fp_read(&A->x1, file);
}

void fp6_write(FILE *file,fp6_t *A){
    fp2_write(file, &A->x0);
    fp2_write(file, &A->x1);
    fp2_write(file, &A->x2);
}
void fp6_read(fp6_t *A,FILE *file){
    fp2_read(&A->x0, file);
    fp2_read(&A->x1, file);
    fp2_read(&A->x2, file);
}

void fp12_write(FILE *file,fp12_t *A){
    fp6_write(file, &A->x0);
    fp6_write(file, &A->x1);
}
void fp12_read(fp12_t *A,FILE *file){
    fp6_read(&A->x0, file);
    fp6_read(&A->x1, file);
}

void efp_write(FILE *file,efp_t *A){
    fp_write(file, &A->x);
    fp_write(file, &A->y);
    fwrite(&A->infinity,sizeof(int),1,file);
}
void efp_read(efp_t *A,FILE *file){
    fp_read(&A->x, file);
    fp_read(&A->y, file);
    fread(&A->infinity,sizeof(int),1,file);
}

void efp2_write(FILE *file,efp2_t *A){
    fp2_write(file, &A->x);
    fp2_write(file, &A->y);
    fwrite(&A->infinity,sizeof(int),1,file);
}
void efp2_read(efp2_t *A,FILE *file){
    fp2_read(&A->x, file);
    fp2_read(&A->y, file);
    fread(&A->infinity,sizeof(int),1,file);
}

void efp6_write(FILE *file,efp6_t *A){
    fp6_write(file, &A->x);
    fp6_write(file, &A->y);
    fwrite(&A->infinity,sizeof(int),1,file);
}
void efp6_read(efp6_t *A,FILE *file){
    fp6_read(&A->x, file);
    fp6_read(&A->y, file);
    fread(&A->infinity,sizeof(int),1,file);
}

void efp12_write(FILE *file,efp12_t *A){
    fp12_write(file, &A->x);
    fp12_write(file, &A->y);
    fwrite(&A->infinity,sizeof(int),1,file);
}
void efp12_read(efp12_t *A,FILE *file){
    fp12_read(&A->x, file);
    fp12_read(&A->y, file);
    fread(&A->infinity,sizeof(int),1,file);
}

void g1_write(FILE *file,g1_t *A){
    efp_write(file,A);
}
void g1_read(g1_t *A,FILE *file){
    efp_read(A,file);
}

void g2_write(FILE *file,g2_t *A){
    efp2_write(file,A);
}
void g2_read(g2_t *A,FILE *file){
    efp2_read(A,file);
}

void g3_write(FILE *file,g3_t *A){
    g3_write(file,A);
}
void g3_read(g3_t *A,FILE *file){
    g3_read(A,file);
}

void fr_write(FILE *file,fr_t *A){
    mpn_write(file,A->x0,FRLIMB);
}
void fr_read(fr_t *A,FILE *file){
    mpn_read(A->x0,file,FRLIMB);
}
