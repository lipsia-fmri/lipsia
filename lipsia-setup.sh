#!/bin/bash
#
# /bk2015
#

[[ $0 != *bash* ]] && { echo $0 must be sourced not executed to take effect.; exit 1; }

command -v realpath >/dev/null || realpath () {
	[[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}";
    }
#
function testDir ()
{
    test -d "$1" || ( echo Mandatory \"$1\" not found ! Sorry.; return 1; )
}

function Prepend ()
{
    local arg1_value=${!1}
    if [ -z $arg1_value ]; then
	export eval $1=$2
    else
	if [ "${arg1_value/$2/}" = "$arg1_value" ]; then
	    eval $1="$2:$arg1_value"
	fi
    fi
}

I_AM=$(realpath ${BASH_SOURCE[0]})
MY_PLACE=${I_AM%/*} # does 'dirname' exist on a mac ?

#############################################
# LIPSIA DEVELOPMENT BASE DIRECTORY
export LIPSIA_DEV=${LIPSIA_DEV:=${MY_PLACE}}

#############################################
# LIPSIA INSTALLATION FOLDER
LIPSIA_DEFAULT_INST_DIR=/opt/lipsia3.0
#LIPSIA_DEFAULT_INST_DIR=/usr
#LIPSIA_DEFAULT_INST_DIR=/usr/local
#LIPSIA_DEFAULT_INST_DIR="$HOME/bin"
#LIPSIA_DEFAULT_INST_DIR=you_name_it

export LIPSIA_INST=${LIPSIA_INST:=$LIPSIA_DEFAULT_INST_DIR}

testDir "${LIPSIA_DEV}/src"     || return $?
testDir "${LIPSIA_DEV}/include" || return $?
testDir "${LIPSIA_DEV}/lib"     || return $?
testDir "${LIPSIA_DEV}/bin"     || return $?

#########################################################################
# adjust if necessary, e.g. if these libraries are not at the usual places or you
# want/have to use your own versions ( will be silently ignored if nonexistent )
: ${LIPSIA_EXT:=$LIPSIA_DEV/ext}

: ${LIPSIA_GSL:=${LIPSIA_EXT}/gsl}
: ${LIPSIA_FFTW:=${LIPSIA_EXT}/fftw}
# it's not on git, so get it the old fashioned way: #####################
# cd ${LIPSIA_DEV:?"set LIPSIA_DEV first!"} && mkdir -p ext && cd $_   && 
# wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.6-pl1.tar.gz && tar xfz $_ &&
# wget wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz && tar xfz $_  &&
# ln -s gsl-latest gsl; ln -s fftw-* fftw  #-- the last one may fail   ##
#########################################################################

export LIPSIA_INCLUDE=${LIPSIA_DEV}/include
export LIPSIA_LIB=${LIPSIA_DEV}/lib
export LIPSIA_BIN=${LIPSIA_DEV}/bin

echo LIPSIA_DEV set to \"${LIPSIA_DEV}\"

#
# everything else is fixed or derived
#
############################################################################
export LIPSIA_CFLAGS="-O2 -ansi -Wall"
# if all binaries will exclusively run on homogenous hardware like a cluster
# of identical machines or maybe on your own PC _only_, try this instead:
#export LIPSIA_CFLAGS="-march=native -O3 -ansi -Wall"
############################################################################

if [ "${LIPSA_DEV}" != "/usr" ]; then
    LIPSIA_CPPFLAGS=-I${LIPSIA_DEV}/include
    LIPSIA_LDFLAGS=-L${LIPSIA_DEV}/lib
fi

# try to adjust the paths to gsl and fftw

# gsl might already be installed here
if [ -d "${LIPSIA_INST}/gsl" ]; then
    LIPSIA_GSL=${LIPSIA_INST}/gsl
    LIPSIA_CPPFLAGS="${LIPSIA_CPPFLAGS} -I${LIPSIA_GSL}/include"
    LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_GSL}/lib"
    Prepend LD_LIBRARY_PATH ${LIPSIA_GSL}/lib
    echo "(using $LIPSIA_GSL)"
else
    if [ -d "${LIPSIA_GSL}" ]; then
	#echo LIPSIA_GSL=$LIPSIA_GSL
	if [ ! -d ${LIPSIA_EXT}/include/gsl ]; then
	    if [ -x ${LIPSIA_GSL}/configure ]; then
		echo making GSL ...
		cd ${LIPSIA_GSL} && {
		    CFLAGS="$LIPSIA_CFLAGS" ./configure --prefix=${LIPSIA_EXT} && make -j 4 && make install
		}
		cd -
	    fi
	fi
	# successfull installation should put it here
	if [ -d ${LIPSIA_EXT}/include/gsl ]; then
	    echo "(using ${LIPSIA_GSL} ...)"
	    LIPSIA_CPPFLAGS="${LIPSIA_CPPFLAGS} -I${LIPSIA_EXT}/include"
	    #   LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_EXT}/lib -rpath=${LIPSIA_EXT}/lib"
	    LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_EXT}/lib"
	    Prepend LD_LIBRARY_PATH $LIPSIA_EXT/lib
	    #echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
	fi
    fi
fi

# fftw might already be installed here
if [ -d "${LIPSIA_INST}/fftw" ]; then
    export LIPSIA_FFTW=${LIPSIA_INST}/fftw
    LIPSIA_CPPFLAGS="${LIPSIA_CPPFLAGS} -I${LIPSIA_FFTW}/include"
    LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_FFTW}/lib"
    Prepend LD_LIBRARY_PATH ${LIPSIA_FFTW}/lib
    echo "(using $LIPSIA_FFTW)"
else
    if [ -d "${LIPSIA_FFTW}" ]; then
	#echo LIPSIA_FFTW=${LIPSIA_FFTW}
	if [ ! -f "${LIPSIA_EXT}/include/fftw3.h" ]; then
	    if [ -x "${LIPSIA_FFTW}/configure" ]; then
		echo makeing FFTW ...
		cd ${LIPSIA_FFTW} && {
		    CFLAGS="$LIPSIA_CFLAGS" ./configure --prefix=${LIPSIA_EXT} && make -j 4 && make install
		}
		cd -
	    fi
	fi
	if [ -f "${LIPSIA_EXT}/include/fftw3.h" ]; then
	    if [ "${LIPSIA_CPPFLAGS/${LIPSIA_EXT}/}" != "${LIPSIA_CPPFLAGS}" ]; then
		echo "(using ${LIPSIA_FFTW} ...)"
		LIPSIA_CPPFLAGS="${LIPSIA_CPPFLAGS} -I${LIPSIA_EXT}/include"
		#     LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_EXT}/lib -rpath=${LIPSIA_FFTW}/lib"
		LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_EXT}/lib"
		Prepend LD_LIBRARY_PATH ${LIPSIA_EXT}/lib
		#echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
	    fi
	fi
    fi
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

# compiler/linker setup
if [ "$(uname)" == "Linux" ]; then
    # for Linux:
    export CC="gcc"
    
    LSB_RELEASE_D="$(lsb_release -d)"
    if [ "${LSB_RELEASE_D/Ubuntu/}" != "${LSB_RELEASE_D}" ]
    then
	# ubuntu's ld per default sets "--as-needed", switch it off:
	LDFLAGS="$LDFLAGS -Xlinker --no-as-needed"
    fi
else
    :
    # for MacOs, e.g.:
    # export CC="/usr/local/bin/gcc-4.8"
fi

####################################################################################
#if [ ${PWD} != ${LIPSIA_DEV}/src ];
#then
#    echo "\"cd \${LIPSIA_DEV}/src\" and try \"make -j 4\"."
#else
#    echo Ok, try \"make\"
#fi
####################################################################################
