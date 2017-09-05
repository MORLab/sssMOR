sssMOR - Changelog
==================

A list of (major) changes between releases, as well as information about MATLAB versions used and toolbox dependencies. Sometimes we add also changes to come to our **roadmap**.
***

Roadmap (changes to come)
-------------------------
- Second order models: direct reduction
***

v2.00 [06 September 2017]
-----------------------

|                 |                                       |
|:----------------|:--------------------------------------|
| Dependencies    | sss Toolbox                           |
| Programmed with | MATLAB R2016b, R2015b, R2015a         |
| Tested with     | MATLAB R2014a, R2015b, R2016b, R2017a |
| on              | Windows 7                             |

### New Features
- DOCUMENTATION
  * added documentation files to a folder doc/ inside sssMOR. In this way, everybody can profit from the doc documentation of sssMOR, not just those who download the release version.
  * added p-Functions publishDoc.p and publishFunction.p to be able to update the documentation whenever headers are changed.
- CPLXPAIRALL
  * added function to sort several arrays with respect to the sorting of the first input argument as of cplxpair
- STABSEP
  * added an own function to separate an ustable model into a stable and antistable part. In addition to the built-in function, it returns also the projection matrices resulting in the submodels.
- L2NORM
  * added function to compute the L2 norm of a dynamical system. This is especially interesting when dealing with unstable models where the H2-norm is undefined.

### Changes
- SSS
  * updated to v2.00
  * SSS is no more a submodule of sssMOR, hence it is not anymore in src/ folder. Instead, sss and sssMOR are now distributed independently and have to be both in the path for sssMOR to work.
- SSRED
  * added load and save functions to the ssRed class
  * allowing conversion from sss or ss to ssRed objectsf
  * made the definition of ssRed objects easier, allowing also the definition of "userDefined" ssRed objects and the definition of arbitrary reduction parameters to be stored in the ssRed object
- IRKA
  * allowing to call irka with a reduced order q instead of initial shifts and/or tangential directions
- ARNOLDI, RK, IRKA, CIRKA, ...
  * added nLU as output argument to count the number of LU decompositions required by the reduction algorithms
- CIRKA
  * added option to take only the stable subsystem of the model function
- RKICOP
  * allowing matrix valued shifts as input (for MIMO models)
- CURE
  * added as output argument a cell array of all reduced models at each step of the cumulative framework

***


v1.08 [02 Frebruary 2017]
-----------------------

|                 |                               |
|:----------------|:------------------------------|
| Dependencies    | sss Toolbox                   |
| Programmed with | MATLAB R2015b                 |
| Tested with     | MATLAB R2014b, R2015b, R2016b |
| on              | Windows 7, Ubuntu 16.04.1 LTS |

### Changes
- TBR
  * changed the definition of the Cholesky (like) factors of the Gramians to unify the notation with the low-rank method. They are now lower triangular matrices (A=LL') as opposed to MATLAB's built-in notation.
- CIRKA, MODELFCT, MODELFCTMOR
  * added functions for the H2-optimal reduction using the model function framework, which leads to significant speedup, especially for very large scale models.
- ADDED EXTRAS
  * moments, isH2opt, getSylvester moved to extras
  * shiftVec factorized as own function
  * ismemberf2, getDesiredOutput added
- SSS
  * updated to v1.03

### Bugfixes
- MODALMOR
  * invariant subspace is computed only using eigs due to its higher efficiency
- TBR
  * checking for empty E matrix if model is ssRed
- RK
  * checking for empty E matrix if model is ssRed
***

v1.07 [05 October 2016]
-----------------------

|                 |                               |
|:----------------|:------------------------------|
| Dependencies    | sss Toolbox                   |
| Programmed with | MATLAB R2015b                 |
| Tested with     | MATLAB R2014b, R2015b, R2016b |
| on              | Windows 7, Ubuntu 16.04.1 LTS |

### Changes
- minor improvements in MODALMOR, TBR, CURE, SPARK
- SSS toolbox: using v1.02

### Bugfixes
- SSRED
  - Adding hidden properties for system matrices to cope  with changed ss-class definition after R2016a
***

v1.06 [16 September 2016]
-------------------------

|                 |                               |
|:----------------|:------------------------------|
| Dependencies    | sss Toolbox                   |
| Programmed with | MATLAB R2015b                 |
| Tested with     | MATLAB R2014b, R2015b, R2016b |
| on              | Windows 7, Ubuntu 16.04.1 LTS |

### Changes
- APP
 	- new functionality: plotting of impulse- and step-response

- SSRED
	- **new class for reduced state-space-models**

- SPARK
	- **new output (reduced model, reduction with porkV or porkW)**

- ISRK
	- **new function: iterative SVD-Rational Krylov Algorithm**

- RKOP
	- **new function: determine optimal expansion point for Laguerre series**

- RKICOP
	- **new function: Rational Krylov with an iteratively calculated optimal point**
***

v1.05 [9 May 2016]
------------------
|                 |                               |
|:----------------|:------------------------------|
| Dependencies    | sss Toolbox                   |
| Programmed with | MATLAB R2015b                 |
| Tested with     | MATLAB R2014b, R2015b, R2016b |
| on              | Windows 7, Ubuntu 16.04.1 LTS |

### Changes
- APP
 - **New sssMOR app added**

- ARNOLDI
	- **new outputs (matrices of Sylvester equation)**
	- **new input (Opts structure)**
	- Opts.orth (1xMGS,2x MGS,DGKS)
	- Opts.dgksTol
	- Opts.reorth (1xMGS, qr)
	- Opts.real (0,'real')
	- Opts.krylov (cascaded krylov basis for SISO)
	- Opts.lse (sparse, full, hess) for solving linear systems of equations

- RK
	- **new outputs (matrices of Sylvester equation)**
	- **new input (Opts structure for arnoldi)**

- IRKA
	- **new outputs (matrices of Sylvester equation)**
	- **new output (kIter)**
	- Opts.suppressverbose for speed-up

- MODALMOR
	- **new output (dominance analysis values)**
	- Opts.dominance: dominance analysis, most dominant eigenvalue reduction

- TBR
	- **Alternating Direction Implicit (ADI) functionality added using the M.E.S.S. toolbox**
	- plot decay of hsv and ask user (if q not defined)
***

v1.00 - First Release [16 November 2015]
-----------------------------------------

|                 |                               |
|:----------------|:------------------------------|
| Dependencies    | sss Toolbox                   |
| Programmed with | MATLAB R2015b                 |
| Tested with     | MATLAB R2014b, R2015b, R2016b |
| on              | Windows 7, Ubuntu 16.04.1 LTS |
