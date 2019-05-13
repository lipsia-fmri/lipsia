#!/bin/bash
#
# /bk2015 - 2019
#

[[ $0 != *bash* ]] && {
    echo $0 must be sourced not executed to take effect.
    exit 1
}

command -v realpath >/dev/null || realpath () {
	[[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}";
    }

VERBOSE=0
SILENCE=0

function Usage () {
    echo -e "Usage:
    	 ${BASH_SOURCE[0]} [OPTION]
	 -v : be more verbose than necessary
	 -s : be silent if everything's fine
	 -h : this help text"
    return
}

function Summary () {
    if [ $SILENCE -eq 0 ];then
	ECHO LIPSIA_INST=\"$LIPSIA_INST\"
	ECHO LIPSIA_DEV=\"${LIPSIA_DEV}\"
	ECHO CC=$CC
	ECHO CXX=$CXX
	ECHO CPPFLAGS="${LIPSIA_CPPFLAGS}"
	ECHO LDFLAGS="${LIPSIA_LDFLAGS}"
	ECHO CFLAGS="${LIPSIA_CFLAGS}"
	ECHO "(using GSL from $LIPSIA_GSL)"
	ECHO
    fi
}

function testDir () {
    test -d "$1" || ( echo Mandatory \"$1\" not found ! Sorry.; return 1; )
}

function Prepend () {
    local arg1_value=${!1}
    if [ -z $arg1_value ]; then
	export eval $1=$2
    else
	if [ "${arg1_value/$2/}" = "$arg1_value" ]; then
	    eval $1="$2:$arg1_value"
	fi
    fi
}

function ECHO () {
    if [ $SILENCE -eq 0 ];then
	echo $*
    fi
}

function Find_DIR () {
    local foundIn=""
    local dirList=$1
    for a in $dirList;do
	foundIn=$(find $a -name $2 -print -quit 2>/dev/null)
	if [ -n $foundIn ];then
	    echo ${foundIn%/*}
	    break;
	fi
    done
    echo
}
function Find_H () {
    local dirList="${LIPSIA_INST} ${LIPSIA_EXT}/include /usr/include /usr/local/include /opt/include"
    echo $(Find_DIR "$dirList" $1)
}
function Find_LIB () {
    local dirList="${LIPSIA_INST} ${LIPSIA_EXT}/lib /usr/lib64 /usr/local/lib64 /opt/lib64"
    echo $(Find_DIR "$dirList" $1)
}

I_AM=$(realpath ${BASH_SOURCE[0]})
MY_PLACE=${I_AM%/*} # does 'dirname' exist on a mac ?

for((ARG=1;ARG<=$#;ARG++)); do
    OPTION=${!ARG}
    case $OPTION in
	-v ) : $((VERBOSE++)); if [ $VERBOSE -gt 1 ]; then ECHO "verbose"; fi ;;
	-s ) SILENCE=1 ;;
	-* ) Usage; return ;;	
    esac    
done

#############################################
# LIPSIA DEVELOPMENT BASE DIRECTORY
export LIPSIA_DEV=${LIPSIA_DEV:=${MY_PLACE}}

testDir "${LIPSIA_DEV}/src"     || return $?
testDir "${LIPSIA_DEV}/include" || return $?
testDir "${LIPSIA_DEV}/lib"     || return $?
testDir "${LIPSIA_DEV}/bin"     || return $?

#############################################
# LIPSIA INSTALLATION FOLDER
LIPSIA_DEFAULT_INST_DIR=/opt/lipsia3.0
#LIPSIA_DEFAULT_INST_DIR=/usr
#LIPSIA_DEFAULT_INST_DIR=/usr/local
#LIPSIA_DEFAULT_INST_DIR="$HOME/bin"
#LIPSIA_DEFAULT_INST_DIR=you_name_it

# might be set already, but if not we use our default
export LIPSIA_INST=${LIPSIA_INST:=$LIPSIA_DEFAULT_INST_DIR}

#########################################################################
# folder to hold packages like GSL to resolve external dependencies which
# are not or shall not or cannot be resolved by the system 
: ${LIPSIA_EXT:=$LIPSIA_DEV/ext}

: ${LIPSIA_GSL:=${LIPSIA_EXT}/gsl}
# it's not on git, so get it the old fashioned way: #####################
# cd ${LIPSIA_DEV:?"set LIPSIA_DEV first!"} && mkdir -p ext && cd $_   && 
# wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz && tar xfz ${_##*/}    &&
# ln -s ${_%.tar.gz} gsl                                               ##
#########################################################################

export LIPSIA_INCLUDE=${LIPSIA_DEV}/include
export LIPSIA_LIB=${LIPSIA_DEV}/lib
export LIPSIA_BIN=${LIPSIA_DEV}/bin

#
# everything else is fixed or derived
#
############################################################################
export LIPSIA_CFLAGS="-O2 -ansi -Wall"
# if all binaries will exclusively run on homogenous hardware like a cluster
# of identical machines or maybe on your own PC _only_, try this instead:
#export LIPSIA_CFLAGS="-march=native -fast -ansi -Wall"
############################################################################

if [ "${LIPSA_DEV}" != "/usr" ]; then
    LIPSIA_CPPFLAGS=-I${LIPSIA_DEV}/include
    LIPSIA_LDFLAGS=-L${LIPSIA_DEV}/lib
fi

# try to adjust the paths to gsl
D=$(Find_H 'gsl_*.h')
#assuming all gsl_*.h files share a common subdirectory (usually "gsl")
if [ -n "$D" ]; then
    GSL_INCLUDEDIR=${D%/*}
fi

GSL_LIBDIR=$(Find_LIB 'libgsl.*')
#echo $GSL_INCLUDEDIR $GSL_LIBDIR

# gsl might already be installed here
if [ -d "${LIPSIA_INST}/gsl" ]; then
    LIPSIA_GSL=${LIPSIA_INST}/gsl
    LIPSIA_CPPFLAGS="${LIPSIA_CPPFLAGS} -I${LIPSIA_GSL}/include"
    LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_GSL}/lib"
    Prepend LD_LIBRARY_PATH ${LIPSIA_GSL}/lib
else
    if [ -d "${LIPSIA_GSL}" ]; then
	#ECHO LIPSIA_GSL=$LIPSIA_GSL
	if [ ! -d ${LIPSIA_EXT}/include/gsl ]; then
	    if [ -x ${LIPSIA_GSL}/configure ]; then
		ECHO making GSL ...
		cd ${LIPSIA_GSL} && {
		    CFLAGS="$LIPSIA_CFLAGS" ./configure --prefix=${LIPSIA_EXT} && make -j 4 && make install
		}
		cd -
	    fi
	fi
	# successfull installation should put it here
	if [ -d ${LIPSIA_EXT}/include/gsl ]; then
	    ECHO "(using ${LIPSIA_GSL} ...)"
	    LIPSIA_CPPFLAGS="${LIPSIA_CPPFLAGS} -I${LIPSIA_EXT}/include"
	    #   LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_EXT}/lib -rpath=${LIPSIA_EXT}/lib"
	    LIPSIA_LDFLAGS="${LIPSIA_LDFLAGS} -L${LIPSIA_EXT}/lib"
	    Prepend LD_LIBRARY_PATH $LIPSIA_EXT/lib
	    #ECHO "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
	fi
    fi
fi

if [ "${LD_LIBRARY_PATH/$LIPSIA_EXT/}" != "${LD_LIBRARY_PATH}" ]; then
    ECHO
    ECHO "( \"${LIPSIA_EXT}/lib\" should be added to LD_LIBRARY_PATH permanently )"
    ECHO
fi

export CPPFLAGS="${LIPSIA_CPPFLAGS}"
export LDFLAGS="${LIPSIA_LDFLAGS}"
export CFLAGS="${LIPSIA_CFLAGS}"

# compiler/linker setup
if [ "$(uname)" == "Linux" ]; then
    # for Linux:
    export CC="gcc"
    export CXX="g++"
    
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

if [ $VERBOSE -ne 0 ];then
    Summary
fi
####################################################################################
if [ ${PWD} != ${LIPSIA_DEV}/src ];
then
    ECHO "\"cd ${LIPSIA_DEV}/src; make \"."
else
    ECHO Ok, try \"make\"
fi
####################################################################################
