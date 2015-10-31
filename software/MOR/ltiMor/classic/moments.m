function M = moments(sys, s0, n)
% MOMENTS - Returns the moments or Markov parameters of an LTI system
%
% Syntax:
%       M = moments(sys, s0, n)
%
%
% Inputs:
%      -sys:    an sss-object containing the LTI system
%      -s0:     (array of) expansion point (inf -> Markov parameters)
%      -n:      (array of) number of moments to be computed
%
%
% Output: 
%      -M:     3D-array of moments / Markov parameters
%
%
% Examples:
%       No examples
% 
% 
%------------------------------------------------------------------
%   This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State Space, Model Order 
%   Reduction and System Analysis Toolbox developed at the Chair of 
%   Automatic Control, Technische Universitaet Muenchen. For updates 
%   and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
%   More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
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
