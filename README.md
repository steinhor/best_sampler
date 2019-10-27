# How to Set Up best_sampler

mkdir build
cd build
Try cmake .., if Eigen or GSL not found, manual way is
cmake .. -DCMAKE_INSTALL_PREFIX=[...]/eigen-eigen-b3f3d4950030/ -DGSL_ROOT_DIR=[...]
make -j

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

## Changing Parameters

There is a file in the run directory called `parameters.dat`. You can change the file from which the sampler reads resonance, decay, and hyper-cell information by changing the macros `RESONANCES_INFO_FILE`, `RESONANCES_DECAYS_FILE`, and `HYPER_INFO_FILE` respectively within `parameters.dat` to different file pathways. Bear in mind that the file reading is set up specifically for the format of the current resonance info and hyper-cell files. The sampler also creates different sampler objects depending on the value of the sigma field and the temperature--you can change the values of both of these in this same file; the variables for these are fairly self-explanatory.
