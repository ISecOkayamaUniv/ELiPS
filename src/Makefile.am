
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
bls12/bls12_test.lo\
bn12/bn12_generate_curve.lo\
bn12/bn12_generate_parameter.lo\
bn12/bn12_init.lo\
bn12/bn12_g1_scm.lo\
bn12/bn12_g2_scm.lo\
bn12/bn12_g3_exp.lo\
bn12/bn12_miller.lo\
bn12/bn12_final_exp.lo\
bn12/bn12_pairing.lo\
bn12/bn12_test.lo

SUBDIRS = bls12 bn12

AM_CPPFLAGS = -I../include

LDADD = bls12/libelipsbls12.la bn12/libelipsbn12.la

libelips_la_SOURCES = define.c\
scalar.c\
mpn.c\
fp.c\
fp2.c\
fp6.c\
fp12.c\
efp.c\
efp2.c\
efp6.c\
efp12.c\
twist.c\
line_calculate.c\
jsf_naf.c\
test.c\
time.c\
../include/ELiPS/define.h\
../include/ELiPS/scalar.h\
../include/ELiPS/mpn.h\
../include/ELiPS/fp.h\
../include/ELiPS/fp2.h\
../include/ELiPS/fp6.h\
../include/ELiPS/fp12.h\
../include/ELiPS/efp.h\
../include/ELiPS/efp2.h\
../include/ELiPS/efp6.h\
../include/ELiPS/efp12.h\
../include/ELiPS/twist.h\
../include/ELiPS/line_calculate.h\
../include/ELiPS/jsf_naf.h\
../include/ELiPS/test.h\
../include/ELiPS/time.h\
../include/ELiPS/bn12.h\
../include/ELiPS/bls12.h
