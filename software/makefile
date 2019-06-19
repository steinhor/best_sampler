include ../makefile_defs.mk

install : installdirs ${SAMPLER_INSTALLDIR}/lib/libsampler.a ${SAMPLER_INSTALLDIR}/include/classdefs.h ${SAMPLER_INSTALLDIR}/include/randy.h ${SAMPLER_INSTALLDIR}/include/sampler.h ${SAMPLER_INSTALLDIR}/include/constants.h ${SAMPLER_INSTALLDIR}/include/parametermap.h ${SAMPLER_INSTALLDIR}/include/resonances.h ${SAMPLER_INSTALLDIR}/include/sf.h ${SAMPLER_INSTALLDIR}/include/master.h ${SAMPLER_INSTALLDIR}/include/part.h ${SAMPLER_INSTALLDIR}/include/misc.h ${SAMPLER_INSTALLDIR}/include/hyper.h ${SAMPLER_INSTALLDIR}/include/meanfield.h ${SAMPLER_INSTALLDIR}/include/eos.h

installdirs : 
	mkdir -p ${SAMPLER_INSTALLDIR}/include;\
	mkdir -p ${SAMPLER_INSTALLDIR}/lib;\

clean :
	rm -f build/*.o lib/*.a

INC=-I include -I ${EIGEN_HOME}

uninstall :
	rm -f ${SAMPLER_INSTALLDIR}/lib/libsampler.a ${SAMPLER_INSTALLDIR}/include/classdefs.h ${SAMPLER_INSTALLDIR}/include/parametermap.h ${SAMPLER_INSTALLDIR}/include/randy.h ${SAMPLER_INSTALLDIR}/include/sampler.h ${SAMPLER_INSTALLDIR}/include/master.h ${SAMPLER_INSTALLDIR}/include/part.h ${SAMPLER_INSTALLDIR}/include/misc.h ${SAMPLER_INSTALLDIR}/include/hyper.h ${SAMPLER_INSTALLDIR}/include/meanfield.h ${SAMPLER_INSTALLDIR}/include/eos.h

SAMPLER_HFILES = include/sampler.h include/classdefs.h include/parametermap.h include/randy.h include/constants.h include/misc.h include/part.h include/sf.h include/master.h include/meanfield.h include/eos.h

SAMPLER_OBJFILES = build/sampler.o build/randy.o build/parametermap.o build/resonances.o build/misc.o build/part.o build/bess.o build/master.o build/makeparts.o build/meanfield.o build/hyper.o build/eos.o

####################

${SAMPLER_INSTALLDIR}/include/master.h : include/master.h
	cp -f include/master.h ${SAMPLER_INSTALLDIR}/include/

${SAMPLER_INSTALLDIR}/include/sampler.h : include/sampler.h
	cp -f include/sampler.h ${SAMPLER_INSTALLDIR}/include/
	
${SAMPLER_INSTALLDIR}/include/classdefs.h : include/classdefs.h
	cp -f include/classdefs.h ${SAMPLER_INSTALLDIR}/include/

${SAMPLER_INSTALLDIR}/include/constants.h : include/constants.h
	cp -f include/constants.h ${SAMPLER_INSTALLDIR}/include/

${SAMPLER_INSTALLDIR}/include/randy.h : include/randy.h
	cp -f include/randy.h ${SAMPLER_INSTALLDIR}/include/

${SAMPLER_INSTALLDIR}/include/parametermap.h : include/parametermap.h
	cp -f include/parametermap.h ${SAMPLER_INSTALLDIR}/include/

${SAMPLER_INSTALLDIR}/include/resonances.h : include/resonances.h
	cp -f include/resonances.h ${SAMPLER_INSTALLDIR}/include/

${SAMPLER_INSTALLDIR}/include/misc.h : include/misc.h
	cp -f include/misc.h ${SAMPLER_INSTALLDIR}/include/

${SAMPLER_INSTALLDIR}/include/part.h : include/part.h
	cp -f include/part.h ${SAMPLER_INSTALLDIR}/include/

${SAMPLER_INSTALLDIR}/include/sf.h : include/sf.h
	cp -f include/sf.h ${SAMPLER_INSTALLDIR}/include/

${SAMPLER_INSTALLDIR}/include/hyper.h : include/hyper.h
	cp -f include/hyper.h ${SAMPLER_INSTALLDIR}/include/

${SAMPLER_INSTALLDIR}/include/meanfield.h : include/meanfield.h
	cp -f include/meanfield.h ${SAMPLER_INSTALLDIR}/include/

${SAMPLER_INSTALLDIR}/include/eos.h : include/eos.h
	cp -f include/eos.h ${SAMPLER_INSTALLDIR}/include/
		
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

build/randy.o : src/randy.cc include/randy.h
	${CPP} -c ${CFLAGS} ${INC} -o build/randy.o src/randy.cc

build/parametermap.o : src/parametermap.cc include/parametermap.h
	${CPP} -c ${CFLAGS} ${INC} -o build/parametermap.o src/parametermap.cc

build/resonances.o : src/resonances.cc include/resonances.h
	${CPP} -c ${CFLAGS} ${INC} -o build/resonances.o src/resonances.cc

build/misc.o : src/misc.cc include/misc.h
	${CPP} -c ${CFLAGS} ${INC} -o build/misc.o src/misc.cc

build/part.o : src/part.cc ${SAMPLER_HFILES}
	${CPP} -c ${CFLAGS} ${INC} -o build/part.o src/part.cc

build/bess.o : src/bess.cc include/sf.h
	${CPP} -c ${CFLAGS} ${INC} -o build/bess.o src/bess.cc

build/meanfield.o : src/meanfield.cc include/meanfield.h
	${CPP} -c ${CFLAGS} ${INC} -o build/meanfield.o src/meanfield.cc	

build/hyper.o : src/hyper.cc include/hyper.h
	${CPP} -c ${CFLAGS} ${INC} -o build/hyper.o src/hyper.cc

build/eos.o : src/eos.cc include/eos.h
	${CPP} -c ${CFLAGS} ${INC} -o build/eos.o src/eos.cc
			
############
