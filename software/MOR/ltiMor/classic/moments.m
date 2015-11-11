function M = moments(sys, s0, n)
% MOMENTS - Returns the moments or Markov parameters of an LTI system
%
% Syntax:
%       M = moments(sys, s0, n)
%
% Description:
%       This function computes n moments, i.e. the negative Taylor series
%       coefficients of the transfer function of an LTI system about a
%       complex frequency s0. If s0 is inf, Markov parameters are
%       computers.
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
%> sys = loadSss('build');
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
% Last Change:  26 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

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
        [L,U,k,l,S]=lu(sys.E, 'vector');
        B=sys.B;
        for i=2:n(iS0)
            b=S\B; b=b(k,:);
            x=L\b;
            x(l,:)=U\x;
            M(:,:,currIdx + i) = sys.C*x - sys.D;
            B=sys.A*x;
        end
    else
        [L,U,k,l,S]=lu(sys.A - s0(iS0)*sys.E, 'vector');
        B=sys.B;
        for i=1:n
            b=S\B; b=b(k,:);
            x=L\b;
            x(l,:)=U\x;
            M(:,:,currIdx+i) = sys.C*x - sys.D;
            B=sys.E*x;
        end
    end
end
