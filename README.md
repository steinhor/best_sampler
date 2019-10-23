# How to Set Up best_sampler

Once you have a copy of best_sampler on your device, enter the directory and open the file "makefile_defs.mk". You will need to change the pathways in this file to match the location of the best_sampler directory on your device and the locations of eigen3 and gsl.

Now enter `best_sampler/software` and type

	make install

If this is successful, exit the software directory and enter `best_sampler/run`. Then type

	make sampler
	./sampler

The output should be

	will read res info from ../local/resinfo/pdg-SMASH.dat
	NResonances:493
	will read decay info from ../local/resinfo/pdg-SMASH.dat
	opening ../local/include/surface_2D.dat
	Exiting ReadHyper2D() happily, TotalVolume=7170.271412, nelements=43461
	nparts=503589	totvol=999999.999911	nparts/totvol=0.503589	nhad=0.503943
	YIPPEE!!!!! We made it all the way through!


