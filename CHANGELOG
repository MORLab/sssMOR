sssMOR - Changelog
==================

A list of (major) changes between releases, as well as information about MATLAB versions used and toolbox dependencies. Sometimes we add also changes to come to our **roadmap**.
***

Roadmap (changes to come)
-------------------------
- Second order models: direct reduction
***

v1.09 [tbd]
-----------------------

|                 |    |
|:----------------|:---|
| Dependencies    |    |
| Programmed with |    |
| Tested with     |    |
| on              |    |


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
