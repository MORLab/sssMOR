%
%rrqr     Rank-Revealing QR Factorization
%
%   [Q,R,p,r] = rrqr(A), where A is m-by-n, produces an m-by-n upper
%   triangular matrix R and an m-by-m unitary matrix Q so that
%        A*P = Q*[ R11 R12 ]
%                [  0  R22 ],
%   p is a permutation vector and r is the rank of A.
%
%   [Q,R,p,r] = rrqr(A,'s') produces the "economy size" decomposition.
%   if m <= n, R is m-by-n and Q is m-by-m, otherwise R is n-by-n and Q
%   is m-by-n.
%
%   [Q,R,p,r] = rrqr(A,tol). 1/tol specifies an upper bound on the
%   condition number of R11. If tol == 0 or tol is unset, tol = "machine
%   precision" is chosen as default. tol must be >= 0.
%   The tol parameter can be combined with the 's' parameter.
%
%   [B,R,p,r] = rrqr(A,C) returns a matrix B so that B=C*Q.
%   The tol parameter is accepted as well.
%
%   [R,p,r] = rrqr(A) is identical to the upper cases but does not compute
%   Q. The tol and 's' parameters are accepted as well.
%
function [varargout]=rrqr2(A,varargin)

if(nargout == 4)
      %Fall [Q,R,p,r] = rrqrGate(A,tol,str) oder
      %Fall [Q,R,p,r] = rrqrGate(A,str,tol) oder
      %Fall [Q,R,p,r] = rrqrGate(A,C,tol)
      if (nargin == 3)
            a = varargin{1};
            b = varargin{2};
            if(a == 's' & isscalar(b))%Fall [Q,R,p,r] = rrqrGate(A,tol,str)
                  [m,n] = size(A);
                  if(m <= n)
                        Q=zeros(n,m);
                        for i=1:1:m
                              Q(i,i) = 1;
                        end
                        [Q,R,p,r]=rrqrGate(A,Q,b);
                        Q = Q(1:m,1:m);
                  else
                        Q=diag(ones(1,m));
                        [Q,R,p,r]=rrqrGate(A,Q,b);
                        Q = Q(:,1:n);
                        R = R(1:n,1:n);
                  end
                  varargout{1} = Q;
                  varargout{2} = R;
                  varargout{3} = p;
                  varargout{4} = r;
            elseif(isscalar(a) & b == 's')%Fall [Q,R,p,r] = rrqrGate(A,str,tol)
                  [m,n] = size(A);
                  if(m <= n)
                        Q=zeros(n,m);
                        for i=1:1:m
                              Q(i,i) = 1;
                        end
                        [Q,R,p,r]=rrqrGate(A,Q,a);
                        Q = Q(1:m,1:m);
                  else
                        Q=diag(ones(1,m));
                        [Q,R,p,r]=rrqrGate(A,Q,a);
                        Q = Q(:,1:n);
                        R = R(1:n,1:n);
                  end
                  varargout{1} = Q;
                  varargout{2} = R;
                  varargout{3} = p;
                  varargout{4} = r;
            elseif(not(isscalar(a)) & isscalar(b))%Fall [Q,R,p,r] = rrqrGate(A,C,tol)
                  [m,n] = size(A);
                  [mC,nC] = size(a);
                  if(not(nC==m))
                        error('Please check your input arguments');
                  end
                  [Q,R,p,r]=rrqrGate(A,a,b);
                  varargout{1} = Q;
                  varargout{2} = R;
                  varargout{3} = p;
                  varargout{4} = r;
            else
                  error('Please check your input arguments');
            end
      elseif (nargin == 2)
            %Fall [Q,R,p,r] = rrqrGate(A,tol) oder
            %Fall [Q,R,p,r] = rrqrGate(A,str) oder
            %Fall [Q,R,p,r] = rrqrGate(A,C)
            a = varargin{1};
            if(a == 's') %duenne Besetzung
                  [m,n] = size(A);
                  if(m <= n)
                        Q=zeros(n,m);
                        for i=1:1:m
                              Q(i,i) = 1;
                        end
                        [Q,R,p,r]=rrqrGate(A,Q,0);
                        Q = Q(1:m,1:m);
                  else
                        Q=diag(ones(1,m));
                        [Q,R,p,r]=rrqrGate(A,Q,0);
                        Q = Q(:,1:n);
                        R = R(1:n,1:n);
                  end
                  varargout{1} = Q;
                  varargout{2} = R;
                  varargout{3} = p;
                  varargout{4} = r;
            elseif(isscalar(a)) %tol wurde uebergeben, dann volle Besetzung.
                  [m,n] = size(A);
                  Q = diag(ones(1,m));

                  [Q,R,p,r]=rrqrGate(A,Q,a);

                  varargout{1} = Q;
                  varargout{2} = R;
                  varargout{3} = p;
                  varargout{4} = r;
            elseif(not(isscalar(a))) %C wurde uebergeben
                  [m,n] = size(A);
                  [mC,nC] = size(a);
                  if(not(nC==m))
                        error('Please check your input arguments');
                  end
                  [Q,R,p,r]=rrqrGate(A,a,0);
                  varargout{1} = Q;
                  varargout{2} = R;
                  varargout{3} = p;
                  varargout{4} = r;
            else
                  error('Please check your input arguments')
            end
      elseif (nargin == 1)
      %Fall [Q,R,p,r] = rrqrGate(A)
            [m,n] = size(A);
            Q = diag(ones(1,m));

            [Q,R,p,r] = rrqrGate(A,Q,0);

            varargout{1} = Q;
            varargout{2} = R;
            varargout{3} = p;
            varargout{4} = r;
      else
            error('Invalid number of input arguments')
      end
elseif (nargout == 3)
      %Fall [R,p,r] = rrqrGate(A,tol)
      if(nargin == 1)
            [R,p,r] = rrqrGate(A);
            varargout{1} = R;
            varargout{2} = p;
            varargout{3} = r;
      elseif(nargin == 2)
            a = varargin{1};
            if(a == 's') %duenne Besetzung
                  [m,n] = size(A);
                  if(m <= n)
                        [R,p,r]=rrqrGate(A,0);
                  else
                        [R,p,r]=rrqrGate(A,0);
                        R = R(1:n,1:n);
                  end
                  varargout{1} = R;
                  varargout{2} = p;
                  varargout{3} = r;
            elseif(a >= 0) %tol wurde uebergeben, dann volle Besetzung.
                  [R,p,r]=rrqrGate(A,a);

                  varargout{1} = R;
                  varargout{2} = p;
                  varargout{3} = r;
            else
                  error('Please check your input arguments')
            end
      elseif(nargin == 3)
            a = varargin{1};
            b = varargin{2};
            if(a == 's' & b >= 0)
                  [m,n] = size(A);
                  if(m <= n)
                        [R,p,r]=rrqrGate(A,b);
                  else
                        [R,p,r]=rrqrGate(A,b);
                        R = R(1:n,1:n);
                  end
                  varargout{1} = R;
                  varargout{2} = p;
                  varargout{3} = r;
            elseif(a >= 0 & b == 's')
                  [m,n] = size(A);
                  if(m <= n)
                        [R,p,r]=rrqrGate(A,a);
                  else
                        [R,p,r]=rrqrGate(A,a);
                        R = R(1:n,1:n);
                  end
                  varargout{1} = R;
                  varargout{2} = p;
                  varargout{3} = r;
            else
                  error('Please check your input arguments')
            end
      else
            error('Invalid number of input arguments')
      end
else
        error('Invalid number of output arguments')
end
