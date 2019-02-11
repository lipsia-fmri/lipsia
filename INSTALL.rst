How to install LIPSIA
===========================


There are two alternatives for installing lipsia. The first possibility is to compile and build everything from scratch on your local machine, which requires the manual installation of several libraries (as described below). The second option is to install a virtual version of lipsia using docker. This is the preferable option for operatings systems outside the linux/unix world and delivers the same performance as the native build.



Install LIPSIA from source (Linux only)
========================================

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

2) Clone the git repository:
 ::
	
    git clone https://github.com/lipsia-fmri/lipsia.git


3) Execute the script "lipsia-setup.sh" (in a shell)
``````````````````````````````````````````````````````
 ::

   cd <lipsia_dir>
   source ./lipsia-setup.sh


4) Compile lipsia
`````````````````````````
 ::

   cd <lipsia_dir>/src
   make
where <lipsia_dir> is the directory of the lipsia repository.
All executables can be found in <lipsia_dir>/bin.


5) Change bash profile
`````````````````````````
 ::

The following changes need to be performed to your local bash profile (for ubuntu ~/.bashrc, for OSX ~/.bash_profile) Furthermore, the library path needs to be set. Add the following two lines to

 ::

    PATH=<lipsia_dir>/bin:$PATH
    export LD_LIBRARY_PATH=<lipsa_dir>/lib:$LD_LIBRARY_PATH



Install LIPSIA using Docker
===============================

1) Install docker on your local machine
`````````````````````````````````````````````````````

Please follow the instructions: https://docs.docker.com/install/


2) Clone the git repository:
 ::
	
    git clone https://github.com/lipsia-fmri/lipsia.git


3) Build the Dockerfile
`````````````````````````````````````````````````````

 ::
   cd <lipsia_dir>
   docker build -t lipsia .


After the installation, you can run any lipsia program by prepending *docker run lipsia*, e.g.

 ::

   docker run -v ${dir}:${dir} lipsia vecm -in ${dir}/fmri.v -mask ${dir}/mask.v -out ${dir}/ecm.v

where $dir is the path to the local directory containing your data file, e.g. ${dir}/fmri.v must exist as a file on your local system. 
