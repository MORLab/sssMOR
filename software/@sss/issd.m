function [issd, numericalAbscissa] = issd(sys)
%ISSD  check strict dissipativity of LTI sss system
%
% Syntax:
%   ISSD(sys)
%   issd = ISSD(sys)
%   [issd,numericalAbscissa] = ISSD(sys)
%
% Description:
%   This function determines wether the LTI, sss system "sys" is 
%   given in a strictly dissipative realization. The computations are meant 
%   to avoid operations on full matrices. However, whenever this is not 
%   possible, a warning is issued.
%
%   If no output is defined, then the result is printed on the screen.
%   Depending on the number of ouputs defined the function can return
%   - issd: a boolean value (1=true, 0=false, NaN= dissipative but not strictly)
%   - numericalAbscissa: i.e. the largest real part of the eigenvalues of the.
%
%   NaN is returned either, when the computation was not possible, or 
%   when the numerical abscissa is zero. In the latter case, the system
%   might be stable (in the sense of Lyapunov) or unstable if the
%   multiplicity of the eigenvalues at the origin is greater than one.
%
% See also:
%   ISPD, EIGS, CHOL, SPARSE,
%
% References:
%   [1] Panzer 2014 Model order reduction by Krylov subspace methods
%       with global error bounds and automatic choice of parameters
%
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                   -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Heiko Panzer, Alessandro Castagnotto
% Last Change:  01 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%%  Parse input
if condest(sys.e)>1e16, error('issd does not support DAEs'),end

%%  Perform computations
% E >0?
isPosDef = ispd(sys.e);
if ~isPosDef
    if nargout == 0, warning('System is not strictly dissipative (E~>0).'); 
    else issd = 0; end
    return
end

% A + A' <0?  
isNegDef = ispd(-sys.a-sys.a');
if isNegDef
    if nargout == 0, fprintf('System is strictly dissipative.\n'); else issd = 1; end
else
    if nargout == 0, warning('System is not strictly dissipative (E>0, A+A''~<0)'); else issd = 0; end
end
    
if nargout==2 % computation of the numerical abscissa required
    p    = 20;		% number of Lanczos vectors
    tol  = 1e-10;	% convergence tolerance
    opts = struct('issym',true, 'p',p, 'tol',tol, 'v0',sum(sys.e,2));
    try
        numericalAbscissa  = eigs((sys.a+sys.a')/2, sys.e, 1, 'la', opts);
    catch err
        warning('Computation of the numerical abscissa failed with message:%s',err.message);
        numericalAbscissa = NaN;
    end
end

