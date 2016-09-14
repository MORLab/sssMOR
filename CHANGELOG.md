# sssMOR - Changelog

## Roadmap (changes to come)

## v1.05
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

## v1.06
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
