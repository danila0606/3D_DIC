# 3D DIC/PT

This is an implementation of a naive method for calculating correlations between images in 3D using DIC. Also, code provides possibility to do Particle Tracking using 3D DIC output as predicitions. 

## Setting up
The main dependency is `Ncorr`. Check required libraries on `http://www.ncorr.com`.
The code has been tested on MAC/Linux.

Build the project after installing all the dependencies:
```bash
cd build/
cmake . -DCMAKE_BUILD_TYPE=Release
cmake --build .
sudo make install
cd ../DIC_3D/build/
cmake . -DCMAKE_BUILD_TYPE=Release
cmake --build .
cd ../bin
python ../GUI.py
```

## References
```
Ncorr: open-source 2D digital image correlation matlab software
J Blaber, B Adair, A Antoniou
Experimental Mechanics 55 (6), 1105-1122
```