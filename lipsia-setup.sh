#!/bin/bash
#
# /bk2015 - 2019
#
# TODO: handle boost dependancy
#

[[ $0 != *bash* ]] && {
    echo $0 must be sourced not executed to take effect.
    exit 1
}

command -v realpath >/dev/null || realpath () {
	[[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}";
    }

I_AM=$(realpath ${BASH_SOURCE[0]})
MY_PLACE=${I_AM%/*} # does 'dirname' exist on a mac ?

VERBOSE=0
SILENCE=0
ISMAC=0

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
	ECHO
	ECHO CC=$CC
	ECHO CXX=$CXX
	ECHO CPPFLAGS="${LIPSIA_CPPFLAGS}"
	ECHO LDFLAGS="${LIPSIA_LDFLAGS}"
	ECHO CFLAGS="${LIPSIA_CFLAGS}"
	ECHO LIPSIA_INST="$LIPSIA_INST"
	ECHO LIPSIA_DEV="${LIPSIA_DEV}"
	if [ $ISMAC -eq 0 ];then
	    ECHO LIPSIA_LD_LIBRARY_PATH="${LIPSIA_LD_LIBRARY_PATH}"
	fi
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
    if [ $SILENCE -eq 0 -a $VERBOSE -ne 0 ];then
	echo $*
    fi
}

function Find_DIR () {
    local foundIn=""
    local dirList=$1
    for a in $dirList;do
	foundIn=$(find -L $a -name $2 -print 2>/dev/null | head -1)
	if [ ${#foundIn} -ne 0 ];then
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
    local dirList="${LIPSIA_INST} ${LIPSIA_EXT}/lib /usr/lib64 /usr/local/lib64 /usr/local/lib /opt/lib"
    echo $(Find_DIR "$dirList" $1)
}

# here we go

for((ARG=1;ARG<=$#;ARG++)); do
    OPTION=${!ARG}
    case $OPTION in
	-v ) : $((VERBOSE++)); if [ $VERBOSE -gt 1 ]; then ECHO "verbose"; fi ;;
	-s ) SILENCE=1 ;;
	-* ) Usage; return ;;	
    esac    
done

if [ $SILENCE -eq 0 -a $VERBOSE -eq 0 ];then
    echo --- ${BASH_SOURCE[0]}: make me quiet with -s or chatty with -v
fi

# compiler/linker setup
if [ "$(uname)" == "Darwin" ]; then
    export ISMAC=1
    export CC=${CC:=$(ls /usr/local/bin/gcc* 2>/dev/null |head -1)}
    export CXX=${CXX:=$(ls /usr/local/bin/g++* 2>/dev/null |head -1)}
    if [ ${#CC} -eq 0 ]; then
	echo '*************************************************************************************************'	
	echo '********* No gcc found in \"/usr/local/bin\" ****************************************************'
	echo '************* Did you do the following ? ********************************************************'
	echo 'xcode-select --install'
	echo '/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"'
	echo 'brew install gcc'
	echo 'brew install gsl'
	echo 
	echo 'In case gcc ist installed somewhere else, please specify $CC and $CXX appropriately.'
	echo '( /usr/bin/gcc (clang) is currently insufficient )'
	echo '*************************************************************************************************'
	return
    else
	# MacOS marks a link with a '@' 
	CC=${CC/@/}
	CXX=${CXX/@/}
    fi
elif [ "$(uname)" == "Linux" ]; then
    # not much to do here
    export CC=${CC:=gcc}
    export CXX=${CXX:=g++}

    if [ "$(uname -a | grep -c Ubuntu)" -ne 0 ]; then
	# ubuntu's ld per default sets "--as-needed", switch it off:
	LDFLAGS="$LDFLAGS -Xlinker --no-as-needed"
    elif [ "$(uname -a | grep -c .fc[1-9]..)" -ne 0 ]; then
	# fedora, should just work
	:
    elif [ "$(uname -a | grep -c .el[6-9].)" -ne 0 ]; then
	# rh,centos,...
	:
    fi
fi

# lipsia development base directory
export LIPSIA_DEV=${LIPSIA_DEV:=${MY_PLACE}}

testDir "${LIPSIA_DEV}/src"     || return $?
testDir "${LIPSIA_DEV}/include" || return $?
testDir "${LIPSIA_DEV}/lib"     || return $?
testDir "${LIPSIA_DEV}/bin"     || return $?

# lipsia installation folder
LIPSIA_DEFAULT_INST_DIR=/opt/lipsia3.0
#LIPSIA_DEFAULT_INST_DIR=/usr
#LIPSIA_DEFAULT_INST_DIR=/usr/local
#LIPSIA_DEFAULT_INST_DIR="$HOME/bin"
#LIPSIA_DEFAULT_INST_DIR=you_name_it

export LIPSIA_INST=${LIPSIA_INST:=$LIPSIA_DEFAULT_INST_DIR}

#########################################################################
# folder to hold packages like GSL to resolve external dependencies which
# are not or shall not or cannot be resolved by the system 
: ${LIPSIA_EXT:=$LIPSIA_DEV/ext}

export LIPSIA_INCLUDE="${LIPSIA_DEV}/include"
export LIPSIA_LIB="${LIPSIA_DEV}/lib"
export LIPSIA_BIN="${LIPSIA_DEV}/bin"
export LIPSIA_LD_LIBRARY_PATH="${LIPSIA_LIB}"

Prepend LD_LIBRARY_PATH "${LIPSIA_LIB}"

############################################################################
export LIPSIA_CFLAGS="-O2 -ansi -Wall"
# if all binaries will exclusively run on homogenous hardware like a cluster
# of identical machines or maybe on your own PC _only_, try this instead:
#export LIPSIA_CFLAGS="-march=native -fast -ansi -Wall"
############################################################################

LIPSIA_CPPFLAGS=-I${LIPSIA_DEV}/include
LIPSIA_LDFLAGS=-L${LIPSIA_DEV}/lib

# try to adjust the paths to gsl
D=$(Find_H 'gsl_math.h')
# gsl_*.h files usually share a common subdirectory
if [ -n "$D" ]; then
    GSL_INCLUDE_BASE=${D%/*}
    if [ "${GSL_INCLUDE_BASE/\/usr\/include/}" = "${GSL_INCLUDE_BASE}" ]; then
	LIPSIA_CPPFLAGS+=" -I${GSL_INCLUDE_BASE}"
    fi
    ECHO GSL_INCLUDE_BASE=$GSL_INCLUDE_BASE
else
    echo '********* GSL header files not found ***********************************'
    echo '* on MacOS install it via \"brew install gsl\"'
    echo '* on Linux use the system package manager (yum, dnf, apt, ...) to install'
    echo '  the gsl development package (libgsl-dev, gsl-devel, ...).'
    echo '* Or get and compile it by yourself like:'
    echo '  cd ${LIPSIA_DEV:?"set LIPSIA_DEV first!"} && mkdir -p ext && cd $_ && '
    echo '  wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz && tar xfz ${_##*/}  && '
    echo '  ln -s ${_%.tar.gz} gsl && cd gsl                                   && '
    echo '  CFLAGS=\"$LIPSIA_CFLAGS\" ./configure --prefix=${LIPSIA_EXT}       && '
    echo '  make -j 4 && make install'
    echo '************************************************************************'
    echo
    return
fi

# try to adjust the paths to zlib
D=$(Find_H zlib.h)
if [ -n "$D" ]; then
    if [ "$D" != "/usr/include" -a "${LIPSIA_CPPFLAGS/$D/}" = "${LIPSIA_CPPFLAGS}" ];then
	LIPSIA_CPPFLAGS+=" -I${D}"
	ECHO added path to zlib.h \($D\) 
    fi
else
    echo '********* ZLIB header files not found **********************************'
#    echo 'on MacOS install it via \"brew install zlib\"' # "xcode-select --install" solves this
    echo 'on Linux use the system package manager (yum, dnf, apt, ...) to install'
    echo 'the zlib development package (zlib1g-dev, zlib-devel, ...).'
    echo
    return
fi

GSL_LIB_DIR=$(Find_LIB 'libgsl.*')
if [ -n "$GSL_LIB_DIR" ]; then
    ECHO GSL_LIB_DIR=$GSL_LIB_DIR
    D="${GSL_LIB_DIR}"
    if [ "${D/\/usr/\/lib/}" = "$D" ];then
	LIPSIA_LDFLAGS+=" -L$D"
    fi
    if [ $ISMAC -eq 0 -a "${D/\/usr\/lib/}" = "$D" -a "${D/\/usr\/local\/lib/}" = "${D}" ];then
	Prepend LIPSIA_LD_LIBRARY_PATH "${D}"
	Prepend LD_LIBRARY_PATH "${D}"
    fi
fi

export CPPFLAGS="${LIPSIA_CPPFLAGS}"
export LDFLAGS="${LIPSIA_LDFLAGS}"
export CFLAGS="${LIPSIA_CFLAGS}"

if [ $ISMAC -eq 0 -a $SILENCE -eq 0 ];then
    echo
    echo "( \"${LIPSIA_LD_LIBRARY_PATH}\" should be added to \$LD_LIBRARY_PATH permanently )"
fi

if [ $VERBOSE -ne 0 ];then
    Summary
fi
####################################################################################
if [ $SILENCE -eq 0 ];then
    if [ "${PWD}" != "${LIPSIA_DEV}/src" ]; then
	echo "Ok: \"cd ${LIPSIA_DEV}/src; make\"."
    else
	echo Ok: try \"make\"
    fi
fi
####################################################################################
