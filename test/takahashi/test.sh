gcc test_g1_scm.c -o test_g1_scm -lelips -lgmp
./test_g1_scm > result_g1_test.txt

gcc test_g2_scm.c -o test_g2_scm -lelips -lgmp
./test_g2_scm > result_g2_test.txt

gcc test_g3_exp.c -o test_g3_exp -lelips -lgmp
./test_g3_exp > result_g3_test.txt

