function M = moments(sys, s0, n, Opts)
% MOMENTS - Returns the moments or Markov parameters of an LTI system
%
% Syntax:
%       M = moments(sys, s0, n)
%
% Description:
%       This function computes n moments, i.e. the negative Taylor series
%       coefficients of the transfer function of an LTI system about a
%       complex frequency s0. If s0 is inf, Markov parameters are
%       computed.
%
%       The i-th moment Mi of a transfer function about s0 is defined as
%
%       $$ M_i := C \left[ (A - s_0E)^{-1}E\right]^i(A - s_0E)^{-1}B $$
% 
%       where i spans from 0 to infinity.
%
%       To compute multiple moments about different shifts, s0 and n can be
%       arrays.
%
%       The ouput is a 3D array with the moments ordered according to the
%       shifts.
%       
% Input Arguments:
%      -sys:	an sss-object containing the LTI system
%      -s0:     (array of) expansion point (inf -> Markov parameters)
%      -n:      (array of) number of moments to be computed
%
% Output Arguments: 
%      -M:		3D-array of moments / Markov parameters
%
% Examples:
%      This code computes 3 moments of the transfer function of the
%      benchmark model 'build' about the frequencies 1, 1+i and 1-i
%
%> sys = loadSss('building');
%> M = moments(sys,[1,1+i,1-i], 3)
%
%      To compute different order moments for different frequencies,
%      specify n as an array of te same length as s0
%
%> sys = loadSss('CDplayer');
%> M = moments(sys,[1,1+i,1-i],[1,2,4])
%
% See Also:
%       rk, lu
%
% References:
%       * *[1] Antoulas (2005)*, Approximation of large-scale dynamical systems
%       * *[2] Panzer (2014)*, Model Order Reduction by Krylov Subspace Methods
%              with Global Error Bounds and Automatic Choice of Parameters
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
% Last Change:  06 April 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Parse input
Def.refine.active   = false; %use iterative refinement for higher accuracy
Def.refine.tol      = 1e-6;  %default refinement tolerance
Def.refine.maxiter  = 1e2;
        
% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end       

if Opts.refine.active, I = speye(size(sys.A)); end %auxiliary variable

%%  Defining execution parameters
nS0 = length(s0);
if nS0 > 1
    if length(n) == 1
        % moment order defined once for all shifts
        nM = nS0 * n;
        n = repmat(n,1,nS0); %create a vector for the subsequent loop
    elseif length(n) == nS0
        % moment order defined individually
        nM = sum(n);
    else
        error('combination of s0 and n is incompatible');
    end
else %only one shift
    nM = n;
end
%   Preallocate
M=zeros(size(sys.C,1),size(sys.B,2),nM);
%%  Compute the moments
for iS0 = 1:nS0
    currIdx = sum(n(1:iS0-1));
    if isinf(s0(iS0))
        M(:,:,currIdx+1) = sys.D;
        if issparse(sys.E)
            [L,U,p,q,R]=lu(sys.E, 'vector');
        else
            [L,U,p,q,R]=lu(sparse(sys.E), 'vector');
        end
        tempR=sys.B; %temp right side
        for iO=2:n(iS0)
           tempS(q,:) = U\(L\(R(:,p)\tempR)); %temp solution
           
           % Iterative refinement
           if Opts.refine.active
               Opts.refine.L = R*I(p,:).'*L; Opts.refine.U = U*I(:,q).';
               tempS = iterativeRefinement(sys.E,tempR,tempS,Opts.refine);
           end
           
           % Compute moment
            M(:,:,currIdx + iO) = sys.C*tempS - sys.D;
            tempR=sys.A*tempS;
        end
    else %finite frequency
        if issparse(sys.A)
            [L,U,p,q,R]=lu(sys.A - s0(iS0)*sys.E, 'vector');
        else
            [L,U,p,q,R]=lu(sparse(sys.A - s0(iS0)*sys.E), 'vector');
        end
        tempR=sys.B;
        for iO=1:n(iS0)
            tempS(q,:) = U\(L\(R(:,p)\tempR));
            
           % Iterative refinement
           if Opts.refine.active
               Opts.refine.L = R*I(p,:).'*L; Opts.refine.U = U*I(:,q).';
               tempS = iterativeRefinement(sys.A - s0(iS0)*sys.E,tempR,tempS,Opts.refine);
           end
            
            M(:,:,currIdx+iO) = sys.C*tempS - sys.D;
            tempR=sys.E*tempS;
        end
    end
end
