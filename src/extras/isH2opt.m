function isH2opt = isH2opt(sys,sysr,s0,Opts)
% ISH2OPT - Evaluate Meier-Luenberger necessary conditions for H2-optimality
%
% Syntax:
%       ISH2OPT(sys,sysr,s0)
%       ISH2OPT(sys,sysr,s0,Opts)
%       isH2opt = ISH2OPT(sys,sysr,s0)
%
% Description:
%       This function evaluates the Maier-Luenberger conditions for
%       H2-optimality according to the input data, i.e. 
%
%       * sys, the high fidelity model
%       * sysr, the reduced order model
%       * s0,  the shifts at which sysr interpolates sys
% 
%       Following conditions must be met:
%
%       * *a)* the reduced eigenvalues are the mirror images of the shifts
%       * *b)* two moments are matched at each shift
% 
%       Then sysr is said to be a locally H2-optimal approximation of sys
%
%		//Note: This function currently works only for SISO systems, but
%       the generalization to MIMO will be implemented in one of the next
%       releases.
%
% Input Arguments:
%		*Required Input Arguments:*
%		- sys/sysr:		original and reduced model to be tested
%		- s0:			vector of shifts for Hermite reduction
%
%		*Optional Input Arguments:*
%		- Opts:			Option structure with follwing fields
%			-.tol:		Relative tolerance 
%                       [{'1e-3'} / positive float]
%
% Output Arguments:
%       - isH2opt:      Boolean
%
% Examples:
%       This code computes an H2-optimal approximation of order 8 to
%       the benchmark model 'building' and uses the function isH2opt to
%       verify if the necessary conditions for optimality are satisfied.
%
%> sys = loadSss('building')
%> [sysr, ~, ~, s0opt] = irka(sys, -eigs(sys,8).');
%> bode(sys,'-',sysr,'--r');
%> isH2opt(sys, sysr, s0opt)
%
% See Also: 
%		moments, irka, spark, eig
%
% References:
%		* *[1] Gugercin et al. (2008)*, H2 model reduction for large-scale linear dynamical systems
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
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  31 Aug 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Parse input and initialize

%   Get execution parameters
Def.tol = 1e-3;
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%   Predefine variables
isH2opt = false; %initialize
%%  Computations 
%   Check that eigenvalues are mirror images of shifts
if any(abs(s0))<1e-3
    condition = (norm(setdiffVec(s0',-conj(eig(sysr))),1))/sysr.n;
else
    condition = (norm(setdiffVec(s0',-conj(eig(sysr))),1)/norm(s0))/sysr.n;
end
if condition <=Opts.tol
    %   check Meier-Luenberger conditions
    for iShift = 1:length(s0)
        m   = moments(sys,s0(iShift),2);
        mr  = moments(sysr,s0(iShift),2);
        if sys.isSiso
            if norm(squeeze(m-mr))/norm(squeeze(m))> Opts.tol
                % moments do not match
                if nargout == 0
                    fprintf('Reduced model is NOT locally H2-optimal to given accuracy\n');
                end
                return
            end
        else
            %MIMO
            error('MIMO implementation of this is missing')
        end
    end
    %   if you reached this point, all verifications passed!
    isH2opt = true;
end
%%  Print result if no ouput was defined
if nargout == 0
    if isH2opt
        fprintf('Reduced model IS locally H2-optimal\n');
    else
        fprintf('Reduced model is NOT locally H2-optimal\n');
    end
end