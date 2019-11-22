include ../makefile_defs.mk

INCLUDELIBDIRS = -I${SAMPLER_INSTALLDIR}/include -I${EIGEN_HOME} -I${GSL_INCLUDE}/ -L${GSL_LIB} -L${SAMPLER_INSTALLDIR}/lib

samplerlibs : 
	make -C ${SAMPLER_HOME} install

sampler : samplermain.cc samplerlibs
	${CPP} -o sampler samplermain.cc ${CFLAGS} ${INCLUDELIBDIRS} -lsampler -lgsl

samplertest : samplertest.cc samplerlibs
	${CPP} -o samplertest samplertest.cc ${CFLAGS} ${INCLUDELIBDIRS} -lsampler -lgsl

spectraltest : spectraltest.cc samplerlibs
	${CPP} -o spectraltest spectraltest.cc ${CFLAGS} ${INCLUDELIBDIRS} -lsampler -lgsl
