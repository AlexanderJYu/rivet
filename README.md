# RIVET



Program to visualize multi-dimensional persistence.
by Michael Lesnick and Matthew Wright
created December 2013

## Requirements

Before starting to build RIVET, make sure you have the following installed:
 
* A C++ compiler (g++ or clang are what we use)
* CMake
* Qt 5
* Boost (including boost serialization)

All of these are generally available using your operating systems's package
manager. For example, on Mac OS X:
    
    brew install cmake qt5 boost
    
On Ubuntu:

    sudo apt-get install cmake qt5-default qt5-qmake qtbase5-dev-tools libboost-all-dev

Make sure to use the qmake from Qt 5, not Qt 4!  You can 
check which version of qmake is on your path, if any, with the command:

    qmake --version
    

## Building

After cloning to $RIVET_DIR:


    cd $RIVET_DIR
    mkdir build
    cd build
    cmake ..
    make
    cd .. 
    qmake  #make sure this is qmake from Qt 5, not Qt 4!
    make
    

You may see compiler warnings during either of the `make` executions. 
These can safely be ignored. 

After this, you will have two executables built: the viewer (RIVET.app, 
on a Mac, RIVET elsewhere), and the computation engine (rivet_console).

It is then necessary to move or symlink the console into the same folder
where the viewer was built. For example, on Mac OS X:

    cd RIVET.app/Contents/MacOS
    ln -s ../../../build/rivet_console
    
On Ubuntu and most other systems:

    ln -s build/rivet_console
    
In the future, all these steps will be automated so that a single cmake
build will create both executables, and put everything in the right place.

## Contributors

Michael Lesnick (Princeton)
Matthew Wright (St. Olaf College)
Bryn Keller (Intel Labs)

TODO: Other contributors please add yourselves!
     
## Contributing
    
We welcome your contribution! Code, documentation, unit tests, 
interesting sample data files are all welcome!

Before submitting your branch for review, please run the following from the
top level RIVET folder you cloned from Github:

```
clang-format -i **/*.cpp **/*.h
```

This will format the source code using the project's established source
code standards (these are captured in the `.clang-format` file in the
project root directory).