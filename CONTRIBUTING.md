# Contributing to sssMOR

Contributions to sssMOR are welcome, as we try to keep the software up to date with **newest MOR routines** and **numerical linear algebra** algorithms.

In the following are a few guidelines to respect when developing code for sssMOR.

If you have any questions on how to contribute, just get in touch with us at sssMOR@rt.mw.tum.de.

## Styleguide
The file ``programmingGuidelines.m`` contains some general programming and naming rules aimed at improving the quality of the code, its readability and the overall functionality. Here is a summary of the most important points.

### Naming conventions
- general variable and function names should be mixed-case, starting with a lower case, and be self explainatory. Acronyms should be lower case. Underscores should be avoided! Abbreviations should be avoided unless they are commonly used in the domain
```
krylovSubspace, shift, reducedOrderModel, html, tum
rom = reducedOrderModel %examples of acceptable abbreviations
fom = fullOrderModel
```

- amounts, number of objects: n*Objects*
```
nShifts, nPoints
```

- avoid plurals
```
shift % one
shiftArray, shiftVec % many
```

- indices: i*Object*
```
iShift, iFile, jPosition
for iPoint = 1:nPoints %example of nested loops
    for jFile = 1:nFiles
    end
end
```

- booleans: use positive bolean names
```
isFound %instead of isNotFound
```


- Structures begin with upper case letter. Fieldnames do not repeat the name of the structure
```
Segment.length, Options.MAX_ITERATIONS, GeneralOptions.fieldName
```

- Functions that return only one value are named by that value
```
mean, moment, maxError
```

### Definition of execution options

Whenever a function can be executed with different parameters or even in different modes, the execution options should be passed at the end of the inputs such as
```
  sysr = cure(sys,..,Opts)
```
Opts is a structure containing some of the parameters that the function
accepts.
Default paramters are defined at the beginning of the function in which
they are used with a default structure
```
  Def.<fieldname> = ...
```
Then, the input should be parsed and the current options structure
updated with the function parseOpts:
```
  if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
  else
        Opts = parseOpts(Opts,Def);
  end     
```

### REFERENCES:
1. Johnson (2014): MATLAB Style Guidelines 2.0
2. Hung (1999): The Pragmatic Programmer

## Testing

As the size of the project and interconnections between functions increases rapidly, it is essential to test the bais functionality of the toolbox every time we committ a change.
For this reason, we have implemented ``unittest`` routines for every function that can be run within the test suites
``test`` and ``testSssMOR``.

Whenever functions are changed and new functionalities added, the unitest checking for the correct usage **must** be adapted as well. You can run the unittest for a single function by executing

```run(testXXX)```

or even of a just one test function within ``testXXX`` by calling

```run(testXXX,'functionYYY')```.

In addition, **before pushing any changes to GIT** you must make sure that all test run clear by executing the function ``test`` (sss + sssMOR) or ``testSssMOR`` and that the header is formatted according to ``headerTemplate.m``.

### Unittest

For the tests, a unitTest environment has been created. New functionalities, as well as new ways of calling the functions, should be included as test cases. You can find more information about this in the "test" directory of sssMOR.

**Important:** in oder to run the tests you need the add the benchmarks from the SLICOT library (available [here](http://www.icm.tu-bs.de/NICONET/benchmodred.html)) to the directory "benchmarks", since ".mat" files are not included by Git.

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

### Testing on RT84 for new releases
Here is a simple workflow to respect when preparing the toolbox for a new release

1. Make sure everything is ready **on your machine**
  1. check that sssMOR (and all subfolders) are in the path
  1. check in the **Add-On Manager** that no other sssMOR or sss version is installed
  1. run ``test`` and make sure all tests run through
  2. run ``publishHelp('all','evalcode',true)`` and make sure that the documentation is generated correctly
  3. merge the **master** branch into the **realease** branch and create a new tag for the release **sssMORvX.XX**

1. Run the tests **on different MATLAB releases** on **Windows** and **Linux**
  1. log into RT84
  1. copy the current raw files for release onto *C:\_sssMORtests*
  2. open following MATLAB versions (Windows)
    * MATLAB2014b
    * MATLAB2015b
    * MATLAB2016a
  1. Open the VM Virtual Box (if empty, add the Virtual Machine under *C:\_Virtualbox_Linux\Ubuntu_64*)
  3. open following MATLAB versions (Linux Virtual Box)
    * MATLAB2014b
    * MATLAB2015b
    * MATLAB2016a
  4. Make sure the current release files are the only sssMOR files in the path
  5. Run ``test``
