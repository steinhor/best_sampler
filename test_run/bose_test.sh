
#compile and run test for bose correction terms
cd ../run
make sampler

cd ../test_run
make bose_test
./bose_test

#graph output
python fcomp.py
