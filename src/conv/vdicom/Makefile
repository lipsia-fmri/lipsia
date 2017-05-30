
PROG = vdicom

SRC  = vdicom.cpp i2o_dcm.cpp
HDRS = i2o_dcm.h

LDLIBS = -lgslcblas -lgsl -lviaio3.0 -lm -lz

LIPSIA ?= /usr/include
LIPSIA_INCLUDE ?= ${LIPSIA}/include
LIPSIA_LIB ?= ${LIPSIA}/lib

INCLUDE_PATH += ${LIPSIA_INCLUDE}
LIB_PATH += ${LIPSIA_LIB}

DBG_OPT=-O2
DBG_OPT=-g

#CXXFLAGS += ${DBG_OPT} $(addprefix -I,${INCLUDE_PATH}) -std=c++11
CXXFLAGS += ${DBG_OPT} $(addprefix -I,${INCLUDE_PATH}) -std=c++0x
LDFLAGS +=  $(addprefix -L,${LIB_PATH})

CC=${CXX}
OBJ=$(SRC:.cpp=.o)

.PHONY : clean tar

${PROG}: ${OBJ}

vdicom.o: ${HDRS}
i2o_dcm.o: ${HDRS}

clean:
	-rm -f ${PROG} ${OBJ} *~

tar:
	tar cfz vdicom.tgz ${SRC} $(HDRS) Makefile
