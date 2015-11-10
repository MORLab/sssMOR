# sssMOR
A sparse state-space, model order reduction toolbox developed at the Chair of 
Automatic Control, Technische Universität München

*Programmed and tested with: MATLAB R2015b*

## Testing

For the tests, a unitTest environment has been created. New functionalities, as well as new ways of calling the functions, should be included as test cases. You can find more information about this in the "test" directory of sssMOR.

**Important:** in oder to run the tests you need the add the benchmarks form the SLICOT library (available [here](http://www.icm.tu-bs.de/NICONET/benchmodred.html)) to the directory "benchmarks", since ".mat" files are not included by Git.

Please make sure to **avoid using following benchmarks** for testing since they are badly conditioned: *LF10, beam, random, SpiralInductorPeec*

Here is a table of different benchmarks and the worst condition number for (A - s0 E) using IRKA shifts:

benchmark       |  O(max(condest(A-s0 E)))
----------------| ------------------------
CDplayer        | 10^5
build           | 10^5
eady            | 10^3
fom             | 10^3
heat-cont       | 10^3
iss             | 10^5
rail_1357       | 10^3
----------------|------------------------ 
LF10            | 10^9
beam            | 10^8
random          | 10^7
SpiralInductorPeec | 10^6 