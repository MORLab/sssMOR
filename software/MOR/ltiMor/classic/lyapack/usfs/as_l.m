function Y = as_l(tr,X)
%
%  Solves linear systems with the real, symmetric, negative definite matrix A, 
%  i.e., Y = inv(A)*X.
%
%  The Cholesky factor of -A is provided as global data. This data must 
%  be generated by calling 'as_l_i' before calling this routine!
%
%  Calling sequence:
%
%    Y = as_l(tr,X)
%
%  Input:
%
%    tr        is not referenced;
%    X         matrix of proper size.
%
%  Output:
%
%    Y         the solution matrix. 
%
% 
%   LYAPACK 1.0 (Thilo Penzl, May 1999)

if nargin~=2
  error('Wrong number of input arguments.');
end

global LP_U

if isempty(LP_U)
  error('This routine needs global data which must be generated by calling ''as_l_i'' first.');
end 

Y = -LP_U\(LP_U'\X);      % Note the minus!


