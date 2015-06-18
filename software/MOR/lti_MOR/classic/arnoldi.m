function [V,AV,EV,Ct] = arnoldi(E,A,b,s0,IP)
% Arnoldi Algorithm for order reduction of LTI SISO systems using
% multiple expansion points
% ------------------------------------------------------------------
% V = arnoldi(E,A,b,s0)
% Inputs:       * E,A,b: System matrices
%               * s0:    Vector of expansion points
% Outputs:      * V: Orthogonal matrix spanning the Krylov subspace
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de)
% Last Change:  12 Oct 2012
% ------------------------------------------------------------------

%fprintf('Multipoint Arnoldi algorithm by RT-TUM (www.rt.mw.tum.de) panzer@tum.de\n');

if ~exist('IP', 'var'), IP=@(x,y) (x'*E*y); end

if size(s0,1)>1
    error('s0 must be a vector containing the expansion points.')
end
q=length(s0); % order of the reduced model

% remove one element of complex pairs
k=find(imag(s0));
if ~isempty(k)
    s0c = cplxpair(s0(k));
    s0(k) = [];
    s0 = [s0 s0c(1:2:end)];
end

% preallocate memory
V=zeros(length(b),q);

EV=zeros(length(b),q);
AV=zeros(length(b),q);
Ct=eye(1,q); %**


for j=1:length(s0)
%     fprintf([repmat(8,1,4) '%3.0f%%'], 100*(j-1)/q);

    % new basis vector
    temp=b; newlu=1; 
    Ct(j)=1; %**
    if j>1
        if s0(j)==s0(j-1)
            temp=V(:,j-1);
            newlu=0;
            Ct(j)=0; %**
        end
    end
    
    if isinf(s0(j))
        % s0=inf, match Markov parameter instead of moment
        if newlu==0
            temp=A*temp;
        end
        if newlu==1
            try
                % compute Cholesky factors of E
                clear L U o p S
                [R,~,S] = chol(sparse(E));
%                 R = chol(sparse(E));
            catch err
                if (strcmp(err.identifier,'MATLAB:posdef'))
                    % E is not pos. def -> use LU instead
                    [L,U,p,o,S]=lu(sparse(E),'vector');
                else
                    rethrow(err);
                end
            end
        end
        if exist('U', 'var')
            temp = invsolve(temp,L,U,p,o,S);
        else
            temp = S*(R\(R'\(S'*temp)));
        end
    else
        if newlu==0
            temp=E*temp;
        end
        if newlu==1
            % vector LU for sparse matrices
            [L,U,p,o,S]=lu(sparse(A-s0(j)*E),'vector');
        end
        temp = invsolve(temp,L,U,p,o,S);
    end
    
    temp_EV = E*temp;
    temp_AV = A*temp;

    % split complex conjugate columns into real (->j) and imag (->j+length(s0c)/2
    if ~isreal(s0(j))
        V(:,j+length(s0c)/2)=imag(temp);
        temp=real(temp);

        EV(:,j+length(s0c)/2)=imag(temp_EV);
        temp_EV=real(temp_EV);

        AV(:,j+length(s0c)/2)=imag(temp_AV);
        temp_AV=real(temp_AV);
    end

%     temp2=temp;
    % orthogonalize vectors
    for i=1:j-1
      h=IP(temp,V(:,i));
      temp=temp-V(:,i)*h;
      temp_EV=temp_EV-EV(:,i)*h;
      temp_AV=temp_AV-AV(:,i)*h;
      Ct(j)=Ct(j)-h*Ct(i);
    end

    % normalize new basis vector
    h = sqrt(IP(temp,temp));
    V(:,j)=temp/h;
    EV(:,j)=temp_EV/h;
    AV(:,j)=temp_AV/h;
    Ct(j) = Ct(j)/h;
   
end

%orthogonalize columns from imaginary components
for j=length(s0)+1:q
    temp=V(:,j);
    temp_EV=EV(:,j);
    temp_AV=AV(:,j);
    for i=1:j-1
      h=IP(temp, V(:,i));
      temp=temp-h*V(:,i);
      temp_AV=temp_AV-AV(:,i)*h;
      temp_EV=temp_EV-EV(:,i)*h;
      Ct(j) = Ct(j)-h*Ct(i);
    end
%     V(:,j)=temp/sqrt(temp'*IP*temp);
    h = sqrt(IP(temp,temp));
    V(:,j)=temp/h;
    AV(:,j)=temp_AV/h;
    EV(:,j)=temp_EV/h;
    Ct(j) = Ct(j)/h;
end

%QR-orthogonalization in case V is badly scaled
s = svd(IP(V,V)); %determinant computation would also work
if max(s)>1.02
    warning('Modified Gram-Schmidt could not manage to fully keep numerics well behaved. Performing a QR decomposition...');
    [Q,~] = qr(V);
    V = Q(:,1:q);
    AV = A*V;
    EV = E*V;
end

% assignin('base', 'Ct', Ct);

% fprintf([repmat(8,1,4) 'done.\n']);
end


function newcol = invsolve(newcol,L,U,p,o,S)
    b=S\newcol;
%     opts1.LT = true;
%     newcol=linsolve(L,b(p,:),opts1);
    newcol=L\b(p,:);
%     opts2.UT = true;
%     newcol=linsolve(U,newcol,opts2);
    newcol(o,:)=U\newcol;
end
