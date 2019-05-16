How to install LIPSIA
===========================


You can either build lipsia from scratch (installing all dependencies) or use the docker version (see below).



Install LIPSIA from source (Linux/Unix/MacOS)
========================================

1) Install the necessary compilers and libraries
`````````````````````````````````````````````````````
* gcc, g++
* gsl dev
* boost dev
* zlib dev
* BLAS (e.g. openblas)
* python: numpy, matplotlib
* git

Ubuntu:
------------
 ::

    sudo apt-get install build-essential libgsl0-dev libboost-dev zlib1g-dev libopenblas-dev python-tk git-core
    sudo pip install numpy matplotlib


MacOS:
-------------

You will need to install

* homebrew
* command-line-tools
* gsl
* git

To install these, open a terminal (Applications/Utilities/Terminal). Then run the following commands:
 ::

    xcode-select --install
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    brew install git
    brew install gsl



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
