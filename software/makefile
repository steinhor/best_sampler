include ../makefile_defs.mk

install : installdirs ${SAMPLER_INSTALLDIR}/lib/libsampler.a ${SAMPLER_INSTALLDIR}/include/pratt_sampler/classdefs.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/randy.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/sampler.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/constants.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/parametermap.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/resonances.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/sf.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/master.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/part.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/misc.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/hyper.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/meanfield.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/eos.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/test.h

installdirs :
	mkdir -p ${SAMPLER_INSTALLDIR}/include;\
	mkdir -p ${SAMPLER_INSTALLDIR}/lib;\

clean :
	rm -f build/*.o lib/*.a

INC=-I include -I ${EIGEN_HOME}

uninstall :
	rm -f ${SAMPLER_INSTALLDIR}/lib/libsampler.a ${SAMPLER_INSTALLDIR}/include/pratt_sampler/classdefs.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/parametermap.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/randy.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/sampler.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/master.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/part.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/misc.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/hyper.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/meanfield.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/eos.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/test.h

SAMPLER_HFILES = include/pratt_sampler/sampler.h include/pratt_sampler/classdefs.h include/pratt_sampler/parametermap.h include/pratt_sampler/randy.h include/pratt_sampler/constants.h include/pratt_sampler/misc.h include/pratt_sampler/part.h include/pratt_sampler/sf.h include/pratt_sampler/master.h include/pratt_sampler/meanfield.h include/pratt_sampler/eos.h include/pratt_sampler/test.h

SAMPLER_OBJFILES = build/sampler.o build/randy.o build/parametermap.o build/resonances.o build/misc.o build/part.o build/bess.o build/master.o build/makeparts.o build/meanfield.o build/hyper.o build/eos.o build/test.o build/resonances_spectral.o

####################

${SAMPLER_INSTALLDIR}/include/pratt_sampler/master.h : include/pratt_sampler/master.h
	cp -f include/pratt_sampler/master.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/sampler.h : include/pratt_sampler/sampler.h
	cp -f include/pratt_sampler/sampler.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/classdefs.h : include/pratt_sampler/classdefs.h
	cp -f include/pratt_sampler/classdefs.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/constants.h : include/pratt_sampler/constants.h
	cp -f include/pratt_sampler/constants.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/randy.h : include/pratt_sampler/randy.h
	cp -f include/pratt_sampler/randy.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/parametermap.h : include/pratt_sampler/parametermap.h
	cp -f include/pratt_sampler/parametermap.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/resonances.h : include/pratt_sampler/resonances.h
	cp -f include/pratt_sampler/resonances.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/misc.h : include/pratt_sampler/misc.h
	cp -f include/pratt_sampler/misc.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/part.h : include/pratt_sampler/part.h
	cp -f include/pratt_sampler/part.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/sf.h : include/pratt_sampler/sf.h
	cp -f include/pratt_sampler/sf.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/hyper.h : include/pratt_sampler/hyper.h
	cp -f include/pratt_sampler/hyper.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/meanfield.h : include/pratt_sampler/meanfield.h
	cp -f include/pratt_sampler/meanfield.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/eos.h : include/pratt_sampler/eos.h
	cp -f include/pratt_sampler/eos.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/include/pratt_sampler/test.h : include/pratt_sampler/test.h
	cp -f include/pratt_sampler/test.h ${SAMPLER_INSTALLDIR}/include/pratt_sampler/

${SAMPLER_INSTALLDIR}/lib/libsampler.a : lib/libsampler.a
	cp -r lib/libsampler.a ${SAMPLER_INSTALLDIR}/lib/

######################

lib/libsampler.a : ${SAMPLER_HFILES} ${SAMPLER_OBJFILES}
	rm -f lib/libsampler.a;\
	ar -ru lib/libsampler.a ${SAMPLER_OBJFILES};

build/master.o : src/master.cc ${SAMPLER_HFILES}
	${CPP} -c ${CFLAGS} ${INC} -o build/master.o src/master.cc

build/sampler.o : src/sampler.cc ${SAMPLER_HFILES}
	${CPP} -c ${CFLAGS} ${INC} -o build/sampler.o src/sampler.cc

build/makeparts.o : src/makeparts.cc ${SAMPLER_HFILES}
	${CPP} -c ${CFLAGS} ${INC} -o build/makeparts.o src/makeparts.cc

build/randy.o : src/randy.cc include/pratt_sampler/randy.h
	${CPP} -c ${CFLAGS} ${INC} -o build/randy.o src/randy.cc

build/parametermap.o : src/parametermap.cc include/pratt_sampler/parametermap.h
	${CPP} -c ${CFLAGS} ${INC} -o build/parametermap.o src/parametermap.cc

build/resonances.o : src/resonances.cc include/pratt_sampler/resonances.h
	${CPP} -c ${CFLAGS} ${INC} -o build/resonances.o src/resonances.cc

build/resonances_spectral.o : src/resonances_spectral.cc include/pratt_sampler/resonances.h
	${CPP} -c ${CFLAGS} ${INC} -o build/resonances_spectral.o src/resonances_spectral.cc

build/misc.o : src/misc.cc include/pratt_sampler/misc.h
	${CPP} -c ${CFLAGS} ${INC} -o build/misc.o src/misc.cc

build/part.o : src/part.cc ${SAMPLER_HFILES}
	${CPP} -c ${CFLAGS} ${INC} -o build/part.o src/part.cc

build/bess.o : src/bess.cc include/pratt_sampler/sf.h
	${CPP} -c ${CFLAGS} ${INC} -o build/bess.o src/bess.cc

build/meanfield.o : src/meanfield.cc include/pratt_sampler/meanfield.h
	${CPP} -c ${CFLAGS} ${INC} -o build/meanfield.o src/meanfield.cc

build/hyper.o : src/hyper.cc include/pratt_sampler/hyper.h
	${CPP} -c ${CFLAGS} ${INC} -o build/hyper.o src/hyper.cc

build/eos.o : src/eos.cc include/pratt_sampler/eos.h
	${CPP} -c ${CFLAGS} ${INC} -o build/eos.o src/eos.cc

build/test.o : src/test.cc include/pratt_sampler/test.h
	${CPP} -c ${CFLAGS} ${INC} -o build/test.o src/test.cc

############
