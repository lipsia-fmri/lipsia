How to install LIPSIA using Singularity
=========================================

The easiest way to install Lipsia is via a singularity container. To build a singularity container for lipsia, execute the following command:

> singularity build lipsia.sif ./singularity-recipe.txt 

This creates a container called "lipsia.sif". Its size is around 220 MByte.
To execute programs within this container, do the following:

> ./lipsia.sif exec vsml

This will execute the lipsia-program called "vsml". Other lipsia-programs can be called in a similar way.




Install LIPSIA from source (Linux/Unix/MacOS)
========================================

1) Install the necessary compilers and libraries
`````````````````````````````````````````````````````
* gcc, g++
* gsl dev
* boost dev
* zlib dev
* python: numpy, matplotlib
* git

Ubuntu:
------------
 ::

    sudo apt-get install build-essential libgsl0-dev libboost-dev zlib1g-dev python-tk git-core
    sudo pip install numpy matplotlib

Fedora:
------------
 ::

    sudo dnf install g++ gcc git zlib-devel boost-devel gsl-devel zlib-devel
    sudo pip install numpy matplotlib

    ( "yum" instead of "dnf" for Red Hat, Centos, and the like )

MacOS:
-------------

You will need to install

* command-line-tools
* homebrew
* gcc
* gsl
* boost

To install these, open a terminal (Applications/Utilities/Terminal). Then run the following commands:
 ::

    xcode-select --install
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    brew install gcc
    brew install gsl
    brew install boost    


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
