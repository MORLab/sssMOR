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
% Last Change:  11 September 2016
% Copyright (c) 2015, 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Parse input
Def.refine      = 0; %use iterative refinement for higher accuracy (0, 'wilkinson','cgs')
Def.refTol      = 1e-6;  %default refinement tolerance
Def.refMaxiter  = 1e2;
Def.lse = 'sparse';
        
% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end  

Opts.krylov='standardKrylov';
index=1;

% Create n array if n is a scalar
if isscalar(n)
    n = repmat(n,1,length(s0));
end

%   Preallocate
M=zeros(size(sys.C,1),size(sys.B,2),sum(n));

%%  Compute the moments
for iO=1:length(s0)
    if isinf(s0(iO))
        % first Markov parameter=sys.D
        M(:,:,index) = sys.D;
        index=index+1;
        n(iO)=n(iO)-1;
    end
    
    Rt=[];
    for j=1:size(sys.B,2)
        Rt=blkdiag(Rt,ones(1,n(iO)));
    end
    
    tempS = solveLse(sys.A,sys.B,sys.E,s0(iO)*ones(1,n(iO)*size(sys.B,2)),Rt,Opts);
   
    
    for j=1:n(iO)
        temp=[];
        for i=1:size(sys.B,2)
            temp=[temp,tempS(:,j+(i-1)*n(iO))];
        end
        M(:,:,index) = sys.C*temp - sys.D;
        index=index+1;
    end
end
end
