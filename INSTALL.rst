How to install LIPSIA
===========================
1) Install the necessary compilers and libraries
`````````````````````````````````````````````````````
* gcc, g++
* fftw3 dev
* gsl dev
* boost dev
* zlib dev

Ubuntu:
 ::

    sudo apt-get install build-essential libfftw3-dev libgsl0-dev libboost-dev zlib1g-dev


OSX:
 ::
     at the moment, lipsia does not support Mac OSX. We are working on the issue...


2) Execute the script "lipsia-setup.sh" (in a shell)
``````````````````````````````````````````````````````
 ::

   cd <lipsa_dir>
   source ./lipsia-setup.sh

3) Compile lipsia
`````````````````````````
 ::

   cd <lipsa_dir>/src
   make
where <lipsa_dir> is the directory of the lipsia repository.
All executables can be found in <lipsa_dir>/bin.

4) Change bash profile
`````````````````````````
 ::

The following changes need to be performed to your local bash profile (for ubuntu ~/.bashrc, for OSX ~/.bash_profile) Furthermore, the library path needs to be set. Add the following two lines to 

 ::

    PATH=<lipsia_dir>/bin:$PATH
    export LD_LIBRARY_PATH=<lipsa_dir>/lib:$LD_LIBRARY_PATH
