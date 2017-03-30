function [s0out,varargout] = cplxpairAll(s0in,varargin)
% IRKA - Iterative Rational Krylov Algorithm
%
% Syntax:
%       sysr                            = IRKA(sys, q)
%       sysr                            = IRKA(sys, s0)
%       sysr                            = IRKA(sys, s0, Rt, Lt)
%       sysr                            = IRKA(sys, s0,..., Opts)
%       [sysr, V, W]                    = IRKA(sys, s0,... )
%       [sysr, V, W, s0]                = IRKA(sys, s0,... )
%       [sysr, V, W, s0, Rt, Lt]        = IRKA(sys, s0,... )
%       [sysr, V, W, s0, Rt, Lt, B_, Sv, Rv, C_, Sw, Lw, kIter] = IRKA(sys, s0,... )
%       [sysr, V, W, s0, Rt, Lt, B_, Sv, Rv, C_, Sw, Lw, kIter,  s0Traj, RtTraj, LtTraj] = IRKA(sys, s0,... )
%
% Description:
%       This function executes the Iterative Rational Krylov
%       Algorithm (IRKA) as proposed by Gugergin and Beattie in [1].
% 
%       The IRKA iteration is conducted to search for an optimal set of 
%       shifts in Krylov subspace-based model reduction. If IRKA converges,
%       then the reduced model is known to be a local optimum with respect
%       to the H2 norm of the error.
%
%       Convergence is determined by observing the shifts and norm of the
%       reduced model over the iterations. This behavior can be changed
%       with the optional Opts structure.
%
%       The algorithm can be called either with a reduced order q, in which
%       case all initial shifts are set at the origin (and tangential
%       directions to ones), or with a prescribed initial set of shifts s0
%       (and tangential directions Rt,Lt).
%
% Input Arguments:  
%       *Required Input Arguments:*
%       -sys:			full order model (sss)
%
%       *Optional Input Arguments:*
%       -q:             reduced order
%       -s0:			vector of initial shifts
%       -Rt/Lt:			initial right/left tangential directions for MIMO
%       -Opts:			structure with execution parameters
%			-.maxiter:	maximum number of iterations;
%						[{50} / positive integer]
%			-.tol:		convergence tolerance;
%						[{1e-3} / positive float]
%			-.type:		choose between different irka modifications;
%						[{''} / 'stab']
%			-.verbose:	show text output during iterations;
%						[{0} / 1]
%			-.stopCrit:	stopping criterion;
%						[{'combAny'} / 's0' / 'sysr' / 'combAll']
%           -.suppressverbose: suppress any type of verbose for speedup;
%                       [{0} / 1]
%           -.(refer to *arnoldi* or *rk* for other options)
%
% Output Arguments:      
%       -sysr:              reduced order model (sss)
%       -V,W:               resulting projection matrices
%       -s0:                final choice of shifts
%       -Rt,Lt:             final choice of tangential directions
%       -B_,Sv,Rv:          matrices of the input Sylvester equation
%       -C_,Sw,Lw:          matrices of the output Sylvester equation
%       -kIter:             number of iterations
%		-s0Traj,RtTraj,..:  trajectory of all shifst and tangential directions for all iterations
%
% Examples:
%       This code computes an H2-optimal approximation of order 8 to
%       the benchmark model 'building'. One can use the function isH2opt to
%       verify if the necessary conditions for optimality are satisfied.
%
%> sys = loadSss('building')
%> [sysr, ~, ~, s0opt] = irka(sys, -eigs(sys,8).');
%> bode(sys,'-',sysr,'--r');
%> isH2opt(sys, sysr, s0opt)
%
% See Also: 
%       arnoldi, rk, isH2opt
%
% References:
%       * *[1] Gugercin et al. (2008)*, H2 model reduction for large-scale linear dynamical systems
%       * *[2] Beattie et al. (2014)*, Model reduction by rational interpolation
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  20 Jan 2017
% Copyright (c) 2016,2017 Chair of Automatic Control, TU Muenchen


s0out = cplxpair(s0in);

if nargin > 1
    [~,cplxSorting] = ismember(s0out,s0in); %determine sorting order
    varargout = varargin; %preallocate
    for iArg = 1:length(varargin)
        varargout{iArg} = varargin{iArg}(:,cplxSorting);
    end
end