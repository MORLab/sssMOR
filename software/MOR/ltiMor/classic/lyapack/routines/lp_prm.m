function [prm,iprm] = lp_prm(A,mth,E)
%
%  Reordering of SPARSE standard/generalized systems
%
%  This routine permutes the states of the system 
%    .
%    x  =  A x + B u
%    y  =  C x
%
%  or the generalized system
%      .
%    M x  =  N x + Btilde u
%      y  =  Ctilde x
%
%  such that the reordered matrix A or the reordered matrices M and N 
%  produce less fill in sparse factorization methods. This
%  preprocessing step is often essential for the good performance
%  of LYAPACK routines. 
%
%  Calling sequence for Sparse reverse Cuthill-McKee ordering:
%
%    [prm,iprm] = lp_prm(A,'rcm')
%    [prm,iprm] = lp_prm(M,'rcm',N)
%
%  Input:
%
%    A, M, N   n-x-n matrices A, M, or N (must be sparse);
%    mth       reordering method, possible values:
%                 'rcm'  -> Sparse reverse Cuthill-McKee ordering
%                 'samd' -> Symmetric approximate minimum degree
%                           permutation
%                 'camd' -> Column approximate minimum degree
%                           permutation (use with care. This is not
%                           a similarity transformation since
%                           it is only applied to the columns. 
%                           We explicitly only include it for 
%                           completeness reasons)
%
%  Output:
%
%    prm       a vector containing a permutation of the numbers 1,...,n
%              such that the matrices A(prm,prm) or M(prm,prm) and
%              N(prm,prm) have a small bandwidth;
%    iprm      a vector containing the inverse permutation to prm.
%
%  Remarks:
%
%    This routine is useful when implementing preprocessing routines
%    'usfs.pre'.
%
%    It might be necessary to re-reorder the data computed by certain
%    Lyapack routines (postprocessing).
% 
%
%  Reference:
%
%  [1] T.A. Davis 
%      Direct Methods for Sparse Linear Systems (Chapter 7)
%      SIAM 2006
%
%  LYAPACK 1.1 (Jens Saak, July 2007)

% Input data not completely checked!

if ~issparse(A), error('A must be sparse'); end

n = size(A,1);

if nargin==2, E = []; end
if length(E)>=1
  if ~issparse(E), error('E must be sparse'); end
else
  E = sparse(n,n);
end
switch mth
 case 'rcm',
  prm = symrcm(spones(A)+spones(E));  
 case 'samd',
  prm = symamd(spones(A)+spones(E));  
 case 'camd',
  prm = colamd(spones(A)+spones(E));    
 otherwise
  error('Requested reordering method unknown!');
end

iprm = zeros(size(prm)); 

for i = 1:n 
  iprm(prm(i)) = i; 
end 



