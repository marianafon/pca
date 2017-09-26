##1) Sem vetorização
#CPPFLAGS=-g -Wall

#2) Com vetorização -> Escolher uma dessas opções #optimized/missed/all
#CPPFLAGS=-g -ftree-vectorize -Wall -O2 -fopt-info-vec-optimized

#3) Com vetorização
#CPPFLAGS=-g -O2 -ftree-vectorize -msse2 -ftree-vectorizer-verbose=10

#4) Com vetorização -> Escolher uma dessas opções #optimized/missed/all ffast-math
CPPFLAGS=-g -ftree-vectorize -Wall -O2 -fopt-info-vec-optimized -ffast-math

#5) Com vetorização -> Escolher uma dessas opções #optimized/missed/all ffast-math
#CPPFLAGS=-g -ftree-vectorize -Wall -O2 -fopt-info-vec-optimized -ffast-math -msse2

BOOST_HOME=/usr
TNT_HOME= ./lib/tnt
JAMA_HOME=./lib/jama

LIBS= -L /usr/include
INCLUDE = -I ${TNT_HOME} -I ${JAMA_HOME}

CPLUS_INCLUDE_PATH=${BOOST_HOME}/include/c++
export CPLUS_INCLUDE_PATH

.PHONY : clean all

all: $(subst .cpp,.o,$(SOURCES))  pca


%.O: %.cpp
	$(CXX) $(CPPFLAGS) ${INCLUDE} ${LIBS} $^ $@
pca: pca.cpp
	$(CXX) $(CPPFLAGS) ${INCLUDE} $^  ${LIBS} -o $@

clean:
	rm -rf  *.o  pca