make viscoustest;
for i in -0.00005 -0.0001 -0.0002 -0.0003 -0.0004 -0.0005 -0.0006; do echo ${i} | caffeinate viscoustest >> logfiles/vt${i}.txt &; done
