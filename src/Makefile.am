
lib_LTLIBRARIES = libelips.la

libelips_la_LIBADD = bls12/bls12_generate_curve.lo\
bls12/bls12_generate_parameter.lo\
bls12/bls12_init.lo\
bls12/bls12_g1_scm.lo\
bls12/bls12_g2_scm.lo\
bls12/bls12_g3_exp.lo\
bls12/bls12_miller.lo\
bls12/bls12_final_exp.lo\
bls12/bls12_pairing.lo\
bls12/bls12_sym.lo\
bls12/bls12_test.lo\
bls/bls_final_exp.lo\
bls/bls_g1.lo\
bls/bls_g2.lo\
bls/bls_generate.lo\
bls/bls_line_calculate.lo\
bls/bls_miller.lo\
bls/bls_pairing.lo\
bls/bls_twist.lo

SUBDIRS = bls12 bls

AM_CPPFLAGS = -I../include

LDADD = bls12/libelipsbls12.la bls/libelipsbls.la

libelips_la_SOURCES = define.c\
define_cvma.c\
scalar.c\
mpn.c\
fp.c\
fp2.c\
fp2cv.c\
fp4cv.c\
fp6.c\
fp12.c\
fp12cv.c\
fpm.c\
fpm2.c\
matrix.c\
efp.c\
efp2.c\
efp2cv.c\
efp6.c\
efp12.c\
efpm.c\
efpm2.c\
twist.c\
line_calculate.c\
jsf_naf.c\
test.c\
field_test.c\
time.c\
fr.c\
file.c\
../include/ELiPS/define.h\
../include/ELiPS/define_cvma.h\
../include/ELiPS/scalar.h\
../include/ELiPS/mpn.h\
../include/ELiPS/fp.h\
../include/ELiPS/fp2.h\
../include/ELiPS/fp2cv.h\
../include/ELiPS/fp4cv.h\
../include/ELiPS/fp6.h\
../include/ELiPS/fp12.h\
../include/ELiPS/fp12cv.h\
../include/ELiPS/fpm.h\
../include/ELiPS/fpm2.h\
../include/ELiPS/matrix.h\
../include/ELiPS/efp.h\
../include/ELiPS/efp2.h\
../include/ELiPS/efp2cv.h\
../include/ELiPS/efp6.h\
../include/ELiPS/efp12.h\
../include/ELiPS/efpm.h\
../include/ELiPS/efpm2.h\
../include/ELiPS/twist.h\
../include/ELiPS/line_calculate.h\
../include/ELiPS/jsf_naf.h\
../include/ELiPS/test.h\
../include/ELiPS/field_test.h\
../include/ELiPS/time.h\
../include/ELiPS/bls.h\
../include/ELiPS/bls12.h\
../include/ELiPS/fr.h\
../include/ELiPS/file.h