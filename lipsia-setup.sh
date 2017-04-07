#!/bin/bash
#
# /bk2015
#

[[ $0 != *bash* ]] && { echo $0 must be sourced not executed to take effect.; exit 1; }

testDir () { test -d "$1" || ( echo Mandatory \"$1\" not found ! Sorry.; return 1; ) }

command -v realpath >/dev/null || realpath () {
	[[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}";
    }
#

I_AM=$(realpath ${BASH_SOURCE[0]})
MY_PLACE=${I_AM%/*} # does 'dirname' exist on a mac ?

# LIPSIA DEVELOPMENT BASE DIRECTORY
export LIPSIA_DEV=${LIPSIA_DEV:=${MY_PLACE}}
#############################################
# LIPSIA INSTALLATION FOLDER
LIPSIA_DEFAULT_INST_DIR=/opt/lipsia3.0
#LIPSIA_DEFAULT_INST_DIR=/usr
#LIPSIA_DEFAULT_INST_DIR=/usr/local
#LIPSIA_DEFAULT_INST_DIR=$HOME/bin
#LIPSIA_DEFAULT_INST_DIR=you_name_it

export LIPSIA_INST=${LIPSIA_INST:=$LIPSIA_DEFAULT_INST_DIR}

testDir "${LIPSIA_DEV}/src"     || return $?
testDir "${LIPSIA_DEV}/include" || return $?
testDir "${LIPSIA_DEV}/lib"     || return $?
testDir "${LIPSIA_DEV}/bin"     || return $?

# adjust if necessary, i.e. these libraries are not at the usual places or you
# want/have to use your own versions ( will be silently ignored if nonexistent )
: ${LIPSIA_EXT:=$LIPSIA_DEV/ext}

export LIPSIA_GSL=${LIPSIA_EXT}/gsl
export LIPSIA_FFTW=${LIPSIA_EXT}/fftw
# it's not on git, so get it the old fashioned way: #####################
# cd ${LIPSIA_DEV:?"set LIPSIA_DEV first!"} && mkdir -p ext && cd $ _  && 
# wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.6-pl1.tar.gz && tar xfz $_ &&
# wget wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz && tar xfz $_  &&
# ln -s gsl-latest gsl; ln -s fftw-* fftw  #-- the last one may fail   ##
## ... and don't forget to adjust LD_LIBRARY_PATH #######################
#########################################################################

export LIPSIA_INCLUDE=${LIPSIA_DEV}/include
export LIPSIA_LIB=${LIPSIA_DEV}/lib
export LIPSIA_BIN=${LIPSIA_DEV}/bin

echo \$LIPSIA_DEV set to \"${LIPSIA_DEV}\"

#
# everything else is fixed or derived
#
############################################################################
export LIPSIA_CFLAGS="-O2 -ansi -Wall"
# if all binaries will exclusively run on homogenous hardware like a cluster
# of identical machines or maybe on your own PC only, try this instead:
#export LIPSIA_CFLAGS="-march=native -O3 -ansi -Wall"
############################################################################

if [ "${LIPSA_DEV}" != "/usr" ]; then
    LIPSIA_CPPFLAGS=-I${LIPSIA_DEV}/include
    LIPSIA_LDFLAGS=-L${LIPSIA_DEV}/lib
fi

if [ -d "${LIPSIA_GSL}" ]; then
    if [ ! -d ${LIPSIA_EXT}/include/gsl ]; then
	echo making GSL ...
	cd ${LIPSIA_GSL} && {
	    CFLAGS="$LIPSIA_CFLAGS" ./configure --prefix=${LIPSIA_EXT} && make -j 4 && make install
	}
	cd -
    fi
    if [ ${LIPSIA_EXT}/include/gsl ]; then
	echo using ${LIPSIA_GSL} ...
	LIPSIA_CPPFLAGS="${LIPSIA_CPPFLAGS} -I${LIPSIA_EXT}/include"
	#   LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_EXT}/lib -rpath=${LIPSIA_EXT}/lib"
	LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_EXT}/lib"
	if [ -z ${LD_LIBRARY_PATH} ]; then
	    export LD_LIBRARY_PATH=$LIPSIA_EXT/lib
	else
	    if [ "${LD_LIBRARY_PATH/$LIPSIA_EXT/}" = "${LD_LIBRARY_PATH}" ]; then
		LD_LIBRARY_PATH="${LIPSIA_EXT}/lib:${LD_LIBRARY_PATH}"
	    fi
	fi
	#echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
    fi
else
    :    # we might check for an installed/usable gsl here
fi

if [ -d "${LIPSIA_FFTW}" ]; then
    echo using ${LIPSIA_FFTW} ...
    if [ ! -f "${LIPSIA_EXT}/include/fftw3.h" ]; then
	echo makeing FFTW ...
	cd ${LIPSIA_FFTW} && {
	    CFLAGS="$LIPSIA_CFLAGS" ./configure --prefix=${LIPSIA_EXT} && make -j 4 && make install
	}
	cd -
    fi
    if [ ! -f "${LIPSIA_EXT}/include/fftw3.h" ]; then
	if [ ${LIPSIA_CPPFLAGS/${LIPSIA_EXT}/} != ${LIPSIA_CPPFLAGS} ]; then
	    LIPSIA_CPPFLAGS="${LIPSIA_CPPFLAGS} -I${LIPSIA_EXT}/include"
	    #     LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_EXT}/lib -rpath=${LIPSIA_FFTW}/lib"
	    LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_EXT}/lib"
	    if [ -z ${LD_LIBRARY_PATH} ]; then
		export LD_LIBRARY_PATH=$LIPSIA_FFTW/lib
	    else
		if [ "${LD_LIBRARY_PATH/$LIPSIA_EXT/}" = "${LD_LIBRARY_PATH}" ]; then
		    LD_LIBRARY_PATH="${LIPSIA_EXT}/lib:${LD_LIBRARY_PATH}"
		fi
	    fi
	    #echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
	fi
    fi
else
    :    # we might check for an installed/usable fftw here
fi

if [ "${LD_LIBRARY_PATH/$LIPSIA_EXT/}" != "${LD_LIBRARY_PATH}" ]; then
    echo
    echo "( \"${LIPSIA_EXT}/lib\" should be added to LD_LIBRARY_PATH permanently )"
    echo
fi

export CPPFLAGS="${LIPSIA_CPPFLAGS}"
export LDFLAGS="${LIPSIA_LDFLAGS}"
export LDLIBS=""
export CFLAGS="${LIPSIA_CFLAGS}"

# compiler setup
if [ `uname` == "Linux " ]; then
    # for Linux: ( why is this necessary ? )
    export CC="gcc"
else
    :
    # for MacOs, e.g.:
    # export CC="/usr/local/bin/gcc-4.8"
fi

####################################################################################
if [ ${PWD} != ${LIPSIA_DEV}/src ];
then
    echo "\"cd \${LIPSIA_DEV}/src\" and try \"make -j 4\"."
else
    echo Ok, try \"make\"
fi
####################################################################################
