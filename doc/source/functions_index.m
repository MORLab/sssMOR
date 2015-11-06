%% Functions Reference
% UI-functions of the sssMOR toolbox sorted by categories

%% LTI MOR Classic
%
% * <.\arnoldihelp.html |arnoldi|>  --  Arnoldi algorithm using multiple expansion points
% * <.\irkahelp.html |irka|>  --  Iterative Rational Krylov Algorithm
% * <.\isH2opthelp.html |isH2opt|>  --  evaluate Maier-Luenberger conditions for H2-optimality
% * <.\modalMorhelp.html |modalMor|>  --  Modal model order reduction of LTI SISO systems
% * <.\momentshelp.html |moments|>  --  Returns the moments or Markov parameters of an LTI system
% * <.\rkhelp.html |rk|>  --  Model Order Reduction by Rational Krylov
% * <.\tbrhelp.html |tbr|>  --  Performs model order reduction by the Truncated Balanced Realization

%% LTI MOR State of the Art
%
% * <.\curehelp.html |cure|>  --  CUmulative REduction framework
% * <.\porkVhelp.html |porkV|>  --  Pseudo-Optimal Rational (Input) Krylov PORK [Wolf et al. 2013]
% * <.\porkWhelp.html |porkW|>  --  Pseudo-Optimal Rational (Output) Krylov PORK [Wolf et al. 2013]
% * <.\sparkhelp.html |spark|>  --  Stability Preserving Adaptive Rational Krylov

%% Sparse State Space
%
% * <.\appendhelp.html |append|>  --  Appends a set of sparse LTI system (sss)
% * <.\bodehelp.html |bode|>  --  Plots the bode diagram of an LTI system
% * <.\c2dhelp.html |c2d|>  --  Converts a sss object from continues to discrete
% * <.\clearhelp.html |clear|>  --  Deletes all state-space dimension related properties
% * <.\connecthelp.html |connect|>  --  Connects a set of sparse LTI system (sss) by evaluating the names of in- and outputs
% * <.\connectSsshelp.html |connectSss|>  --  Connects an appended sparse state space  LTI system (sss) with feedback matrix K
% * <.\decayTimehelp.html |decayTime|>  --  Computes the time period in which a sparse LTI system levels off
% * <.\diaghelp.html |diag|>  --  Transforms an LTI system to (block)-diagonal representation
% * <.\disphelp.html |disp|>  --  Displays information about a sparse state-space model
% * <.\eighelp.html |eig|>  --  Compute eigenvalues and eigenvectors of a sparse state-space model
% * <.\eigshelp.html |eigs|>  --  compute eigenvalues of the sparse state space system using sparse matrices.
% * <.\freqresphelp.html |freqresp|>  --  Evaluates complex transfer function of LTI systems
% * <.\impulsehelp.html |impulse|>  --  Computes and/or plots the impulse response of a sparse LTI system
% * <.\issdhelp.html |issd|>  --  Check strict dissipativity of LTI sss system
% * <.\isstablehelp.html |isstable|>  --  check stability of LTI sss system
% * <.\normhelp.html |norm|>  --  Computes the p-norm of an sss LTI system
% * <.\plushelp.html |plus|>  --  Computes sum of two LTI systems: u-->(sys1+sys2)-->y
% * <.\pzmaphelp.html |pzmap|>  --  Pole-zero plot of sparse state-space system
% * <.\residuehelp.html |residue|>  --  Computes residues, poles and feedthrough of an LTI system
% * <.\sigmahelp.html |sigma|>  --  Plots the amplitude of the frequency response of an LTI system
% * <.\simhelp.html |sim|>  --  Simulates a sss system using iddata input time series
% * <.\simBackwardEulerhelp.html |simBackwardEuler|>  --  Integrates sss model using backward Euler
% * <.\simDiscretehelp.html |simDiscrete|>  --  Integrates discrete time model
% * <.\simForwardEulerhelp.html |simForwardEuler|>  --  Integrates sss model using forward Euler
% * <.\simRK4help.html |simRK4|>  --  Integrates sss model using Runge-Kutta 4
% * <.\sizehelp.html |size|>  --  Computes the size of a sparse LTI system (sss)
% * <.\spyhelp.html |spy|>  --  plot sparsity pattern of sss system
% * <.\sshelp.html |ss|>  --  Converts sparse LTI system (sss) to Matlab\control\ss
% * <.\stephelp.html |step|>  --  Computes and/or plots the step response of a sparse LTI system
% * <.\truncatehelp.html |truncate|>  --  truncates a sparse LTI system (sss)
% * <.\minushelp.html |minus|>  --  Computes difference of two LTI systems: u-->(sys1-sys2)-->y
% * <.\mtimeshelp.html |mtimes|>  --  mtimes computes the product of two LTI systems: u-->sys2-->sys1-->y

%% Header Templates
%
% * <.\headerTemplatehelp.html |headerTemplate|>  --  sssMOR toolbox header template
% * <.\headerTemplateSsshelp.html |headerTemplateSss|>  --  sss toolbox header template
%
%%
% <html>
%   <hr>
%   <p class="copy">&copy; 2015 RT Technische Universität München
%        <tt class="minicdot">&#149;</tt>
%        <a href="http://www.rt.mw.tum.de">Website</a>
%        <tt class="minicdot">&#149;</tt>
%        <a href="file:./LICENSE.txt">Terms of Use</a>
%        <tt class="minicdot">&#149;</tt>
%        <a href="file:./README.txt">Read me</a>
%   </p>
% </html>
