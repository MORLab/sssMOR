function M = moments(sys, s0, n)
% Returns the moments or Markov parameters of an LTI system
% ------------------------------------------------------------------
% m = moments(sys, s0, n)
% Inputs:       * sys: an sss-object containing the LTI system
%               * s0: expansion point (inf -> Markov parameters)
%               * n: number of moments to be computed
% Output:       * M: array of moments / Markov parameters
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer, Alessandro Castagnotto
% Last Change:  26 Oct 2015
% ------------------------------------------------------------------

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
