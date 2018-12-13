How to install LIPSIA
===========================


There are two alternatives for installing lipsia. The first possibility is to compile and build everything from scratch on your local machine, which requires the manual installation of several libraries (as described below). The second option is to install lipsia using docker, which may be a easier in some cases, especially for non-linux systems.



Install LIPSIA from scratch
===============================

1) Install the necessary compilers and libraries
`````````````````````````````````````````````````````
* gcc, g++
* fftw3 dev
* gsl dev
* boost dev
* zlib dev
* BLAS (e.g. openblas)
* python: numpy, matplotlib

Ubuntu:
 ::

    sudo apt-get install build-essential libfftw3-dev libgsl0-dev libboost-dev zlib1g-dev libopenblas-dev python-tk
    sudo pip install numpy matplotlib


Mac OSX:
 ::

     #open a terminal: Applications/Utilities/Terminal

     #install command line tools
     xcode-select --install

     #install homebrew package manager
     ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

     #install GCC
     brew tap homebrew/versions
     brew install [flags] gcc48

     #install gsl
     brew install gsl

     #install fftw
     brew install fftw

     Furthermore, make sure the last line in the file "lipsia-setup.sh" should be:
     export CC="/usr/local/bin/gcc-4.8"


2) Execute the script "lipsia-setup.sh" (in a shell)
``````````````````````````````````````````````````````
 ::

   cd <lipsia_dir>
   source ./lipsia-setup.sh


3) Compile lipsia
`````````````````````````
 ::

   cd <lipsia_dir>/src
   make
where <lipsia_dir> is the directory of the lipsia repository.
All executables can be found in <lipsia_dir>/bin.


4) Change bash profile
`````````````````````````
 ::

The following changes need to be performed to your local bash profile (for ubuntu ~/.bashrc, for OSX ~/.bash_profile) Furthermore, the library path needs to be set. Add the following two lines to

 ::

    PATH=<lipsia_dir>/bin:$PATH
    export LD_LIBRARY_PATH=<lipsa_dir>/lib:$LD_LIBRARY_PATH



Install LIPSIA using Docker
===============================

The first step is to install docker on your local machine, see: https://docs.docker.com/install/
The second step is to build the docker image as follows:

 ::

   docker build -t lipsia .

After the installation, you can run any lipsia program by prepending *docker run lipsia*, e.g.

 ::

   docker run -v ${dir}:${dir} lipsia vecm -in ${dir}/fmri.v -mask ${dir}/mask.v -out ${dir}/ecm.v

where $dir is the path to the local directory, which is needed by docker. (e.g. dir=/home/user/).
