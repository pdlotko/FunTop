README file for the code directory and its subdirectories.

AUTHORS: 
   Pawel Dlotko    
   Thomas Wanner
   Thomas Stephens tstephe3@gmu.edu
   2013,2014,2015,2016

THIS SOFTWARE IS LICENSED UNDER THE GPL

Required third-party dependencies:
----------------------------------

cxsc  http://www.xsc.de/
phat  https://code.google.com/p/phat/
-------
GNU g++ compiler 4.6 or later(?)
GNU make
(OPTIONAL) Python2.7.x with Matplotlib 1.3.1 (or later)
-------

Installation:
-------------

For UNIX-like platforms, including Cygwin:
 
Step 1. Run the install script int this directory
        
            ./install
         
        With your permissions, this will do the following:
        
        1. Unarchive packages/phat_1_2_1.tar.gz in the
           directory dependencies/phat_1_2_1 (or the directory of your choice)
        2. Unarchive packages/cxsc-2-5-4_RANVALINSTALL.tar.gz 
           into the directory dependencies/cxsc-2-5-4 (or the directory of your
           choice)
        3. Take you through a modified install script for the cxsc
           library.  All defaults *should be* okay.
        4. Autogenerate Makefiles in the application directories
           below ranval3d/code/applications/ with linkages based on directories
           established in Step 1 and Step 2
        
        * NOTE: sometimes linking in a Makefile is not enough to
          guarantee the linkage at runtime -- you may need to append
          the path to the cxsc lib directories to the LD_LIBRARY_PATH
          environment variable.  The best way to do this is to add 
          
            export LD_LIBRARY_PATH=/path/to/cxsc/lib:$LD_LIBRARY_PATH

          to your ~/.bashrc file.

Usage:
------

The library contain three main functionalities presented in src/persistent_homology
First one, presented in computePersistenceOfFunction.cpp allows to compute persistent
homology of a scalar value function on a given rectangle.
Second one, presented in compute_PC_approximation_scalar_valued_function.cpp allows
to compute partially constant approximation of a scalar valued function.
Third, presented in compute_PC_approximation_vector_valued_function.cpp allows
to compute partially constant approximation of a vector valued function.

Numerous examples of various function can be found in src/persistent_homology/configuration.hpp
Feel free to add your own.

For any questions, please contact pdlotko at impan dot pl
