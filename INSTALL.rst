How to install LIPSIA
===========================
1) Install the necessary compilers and libraries
* gcc, g++
* fftw3 dev
* gsl dev
* boost dev

Ubuntu:
``````````
 ::

    sudo apt-get install build-essential libfftw3-dev libgsl0-dev libboost-dev


OSX:
`````````

 ::

     #install homebrew package manager
     ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

     #install GCC
     brew tap homebrew/versions
     brew install [flags] gcc48

     #install gsl
     brew install gsl

     #install fftw
     brew install fftw

     Furtermore, make sure the last line in the file "lipsia-setup.sh" should be:
     export CC="/usr/local/bin/gcc-4.8"


2) Then execute the script "lipsia-setup.sh" using the command line

 ::

   source ./lipsia-setup.sh

3) Compile lipsia

 ::

   cd <lipsa_dir>/src
   make
where <lipsa_dir> is the directory of the lipsia repository.
All executables can be found in <lipsa_dir>/bin.
We recommend adding them to the .bashrc, so that you can start lipsia from anywhere in the terminal. Add this line to ~/.bashrc (for ubuntu, for OSX ~/.bash_profile)

 ::

    PATH=<lipsia_dir>/bin:$PATH
