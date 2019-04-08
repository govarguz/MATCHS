# CMake generated Testfile for 
# Source directory: /u/gvargas/code/e++SCv1/testsuite/ljcos
# Build directory: /u/gvargas/code/e++SCv1/buildV1/testsuite/ljcos
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(ljcos "/mpcdf/soft/SLE_12_SP3/packages/x86_64/anaconda/2.5.1.0/bin/python2" "/u/gvargas/code/e++SCv1/testsuite/ljcos/test_ljcos.py")
set_tests_properties(ljcos PROPERTIES  ENVIRONMENT "PYTHONPATH=/u/gvargas/code/e++SCv1/buildV1:/u/gvargas/code/e++SCv1/buildV1/contrib:")
