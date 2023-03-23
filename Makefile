ifndef STAPL
  STAPL = $(shell echo "/home/seliny2/stapl-developer")
endif
include $(STAPL)/GNUmakefile.STAPLdefaults

OBJS=$(shell ls *.cpp | sed 's/.cpp//g')
LIB+=-l$(shell ls  $(Boost_LIBRARY_DIRS) | grep -e boost_program_options | head -n 1 | sed 's/lib//g' | sed 's/\.so//g' | sed 's/\.a//g')

default: compile

test: all

all: compile

compile: $(OBJS)

clean:
	rm -rf *.o core* a.out ii_files rii_files $(OBJS)
	rm serial
	rm parallelBLS
	rm debug
	rm ompBLS
	rm omp-serial


serial: serial.cpp
	mpic++ -D_STAPL -I/home/seliny2/stapl-developer/./tools/libstdc++/11.2.0 -I/home/seliny2/stapl-developer/tools -I/home/seliny2/stapl-developer -Wall -Werror -Wno-unused-local-typedefs -Wno-unknown-pragmas -Wno-misleading-indentation -Wno-deprecated-declarations -Wno-aligned-new -Wno-noexcept-type -Wno-parentheses -fno-new-ttp-matching  -I/usr/local/boost-1.78.0/include -DBOOST_RESULT_OF_USE_TR1_WITH_DECLTYPE_FALLBACK -o serial serial.cpp -L/home/seliny2/stapl-developer/lib -lstapl -lrt -L/usr/local/boost-1.78.0/lib -lboost_serialization -lboost_system -I/home/seliny2/exoplanet/CCfits-2.6/usr/local/include -I/home/seliny2/exoplanet/cfitsio-4.2.0/include -L/home/seliny2/exoplanet/CCfits-2.6/usr/local/lib -L/home/seliny2/exoplanet/cfitsio-4.2.0/lib -lCCfits -lcfitsio -lboost_program_options
	# -I/home/seliny2/exoplanet/CCfits-2.6/usr/local/include \
	#-I/home/seliny2/exoplanet/cfitsio-4.2.0/include \
	#-L/home/seliny2/exoplanet/CCfits-2.6/usr/local/lib \
	#-L/home/seliny2/exoplanet/cfitsio-4.2.0/lib \
	#-lCCfits \
	#-lcfitsio

parallelBLS: parallelBLS.cpp
	mpic++ -D_STAPL -I/home/seliny2/stapl-developer/./tools/libstdc++/11.2.0 -I/home/seliny2/stapl-developer/tools -I/home/seliny2/stapl-developer -Wall -Werror -Wno-unused-local-typedefs -Wno-unknown-pragmas -Wno-misleading-indentation -Wno-deprecated-declarations -Wno-aligned-new -Wno-noexcept-type -Wno-parentheses -fno-new-ttp-matching -I/usr/local/boost-1.78.0/include -DBOOST_RESULT_OF_USE_TR1_WITH_DECLTYPE_FALLBACK -o parallelBLS parallelBLS.cpp -L/home/seliny2/stapl-developer/lib -lstapl -lrt -L/usr/local/boost-1.78.0/lib -lboost_serialization -lboost_system -I/home/seliny2/exoplanet/CCfits-2.6/usr/local/include -I/home/seliny2/exoplanet/cfitsio-4.2.0/include -L/home/seliny2/exoplanet/CCfits-2.6/usr/local/lib -L/home/seliny2/exoplanet/cfitsio-4.2.0/lib -lCCfits -lcfitsio -lboost_program_options

debug: find_bug.cpp
	mpic++ -D_STAPL -I/home/seliny2/stapl-developer/./tools/libstdc++/11.2.0 -I/home/seliny2/stapl-developer/tools -I/home/seliny2/stapl-developer -Wall -Werror -Wno-unused-local-typedefs -Wno-unknown-pragmas -Wno-misleading-indentation -Wno-deprecated-declarations -Wno-aligned-new -Wno-noexcept-type -Wno-parentheses -fno-new-ttp-matching -I/usr/local/boost-1.78.0/include -DBOOST_RESULT_OF_USE_TR1_WITH_DECLTYPE_FALLBACK -o debug find_bug.cpp -L/home/seliny2/stapl-developer/lib -lstapl -lrt -L/usr/local/boost-1.78.0/lib -lboost_serialization -lboost_system -I/home/seliny2/exoplanet/CCfits-2.6/usr/local/include -I/home/seliny2/exoplanet/cfitsio-4.2.0/include -L/home/seliny2/exoplanet/CCfits-2.6/usr/local/lib -L/home/seliny2/exoplanet/cfitsio-4.2.0/lib -lCCfits -lcfitsio -lboost_program_options

omp-serial : omp-serialBLS.cpp
	g++  -Wall -g -o  omp-serial omp-serialBLS.cpp  -fopenmp -I/home/seliny2/exoplanet/CCfits-2.6/usr/local/include -I/home/seliny2/exoplanet/cfitsio-4.2.0/include -L/home/seliny2/exoplanet/CCfits-2.6/usr/local/lib -L/home/seliny2/exoplanet/cfitsio-4.2.0/lib -lCCfits -lcfitsio

ompBLS: ompBLS.cpp
	g++  -Wall -g -o  ompBLS ompBLS.cpp  -fopenmp -I/home/seliny2/exoplanet/CCfits-2.6/usr/local/include -I/home/seliny2/exoplanet/cfitsio-4.2.0/include -L/home/seliny2/exoplanet/CCfits-2.6/usr/local/lib -L/home/seliny2/exoplanet/cfitsio-4.2.0/lib -lCCfits -lcfitsio