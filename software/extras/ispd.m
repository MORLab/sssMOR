function ispd = ispd(A)
% ISPD - Determine if a matrix is positive definite 
%
% Syntax:
%   ISPD(A)
%   ispd = ISPD(A)
%
% Description:
%   Positive (or negative) definiteness is an important property of
%   matrices that can be exploited if given. Therefore, it is often
%   desirable to determine whether a matrix has this property or not.
%
%   This function tries a Cholesky factorization of A. If it succedes, then
%   the matrix is positive definite.
%
% See also:
%   CHOL, ISSD
%
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Last Change:  01 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%%  Computation
[~,p] = chol(A); 
ispd = ~p>0;

%%  Print result if no ouput was defined
if nargout == 0
    if p>0
        fprintf('Matrix is NOT positive definite \n');
    else
        fprintf('Matrix IS positive definite \n');
    end
end
