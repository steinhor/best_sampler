#! /bin/bash
make viscoustest;
for i in 0.0001 0.0002;
	`echo "-"${i} | ./viscoustest >> logfiles/vt${i}.txt` &
done
