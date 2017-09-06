function [rk_stand,rk_cascad,Vstand,Vcascad,alpha,aa,bb] = test1(Opts)


sys = loadSss('iss');
A = sys.A;
B = sys.B;
C = sys.C;
E = sys.E;

m = size(B,2);      n = size(A,2);
p = size(C,1);


% compute shifts using irka
shifts = zeros(1,10);
Rt = ones(m,size(shifts,10));
Lt = ones(p,size(shifts,10));
[~, ~, ~, s, T, Lt] = irka(sys, shifts,Rt,Lt);
%T = ones(m,size(s,2));


% build standard basis
for ii= 1:1:2*size(s,2)
    if ii > size(s,2)
        jj = ii-size(s,2);
    else
        jj = ii;
    end
    Anew = (A-s(1,jj)*E);
    Ainv = (Anew\ones(n));
    vnew = Ainv*B;%*T(:,jj);
%     rhs = B;
%     AinvB = (Anew,rhs);
%     vnew = AinvB*T(:,jj);
   jstart = ii*size(vnew,2)-(size(vnew,2)-1);  jende = ii*size(vnew,2);
    Vstand(:,jstart:jende) = vnew;
    Vstand = real(Vstand); 
    tic
    Vstand = gramSchmidt(ii,Vstand,Opts);
    aa=toc
end

% build cascaded basis
jj = 1;
for ii= 1:1:2*size(s,2)
   if ii == 1
       Anew = (A-s(1,jj)*E);
       Ainv = (Anew\ones(n));
       Ainv_old = Ainv;
       vnew = Ainv*B;%*T(:,jj);
%       rhs = B;
%        Ainv = solveLse(Anew,rhs);
%        Ainv_old = Ainv;
%        vnew = Ainv*T(:,jj);
       jstart = ii*size(vnew,2)-(size(vnew,2)-1);  jende = ii*size(vnew,2);
       Vcascad(:,jstart:jende) = vnew;
       Vcascad = real(Vcascad);
       jj = jj+1;
   else
       if ii > size(s,2)
           jj = ii-size(s,2);
       else
           jj = ii;
       end
       Anew = (A-s(1,jj)*E);
       Ainv = (Anew\ones(n))*E;
       Adyn = Ainv*E*Ainv_old;
       Ainv_old = Adyn;
       vnew = Adyn*B;%*T(:,jj);
%       rhs = B*T(:,jj);
%        rhs = Ainv_old;
%        %Ainv = solveLse(Anew,rhs);
%        Ainv = Anew\rhs;
%        Ainv_old = Ainv;
%        vnew = Ainv*T(:,jj);
       jstart = ii*size(vnew,2)-(size(vnew,2)-1);  jende = ii*size(vnew,2);
       Vcascad(:,jstart:jende) = vnew;
       Vcascad = real(Vcascad);
   end
   tic
   Vcascad = gramSchmidt(ii,Vcascad,Opts);
   bb=toc
end


rk_stand = rank(Vstand);
rk_cascad = rank(Vcascad);

alpha = subspace(Vstand,Vcascad);





function [V, TRv, W, TLw] = gramSchmidt(jCol, V, varargin)
%   Gram-Schmidt orthonormalization
%   Input:  jCol:  Column to be treated
%           V, W:  Krylov-Subspaces
%   Output: V, W:  orthonormal basis of Krylov-Subspaces
%           TRv, TLw: Transformation matrices

% global variable hermite
global hermite_gram_sch 
%if hermite_gram_sch == 1,    hermite = 1;    else    hermite = 0;    end

% input
if isa(varargin{end},'struct')
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end

if length(varargin) == 1
    W = varargin{1};
else
    W = eye(size(V));
end
    
% hier muss ich jetzt noch das IP standardm??ig definieren
IP = @(x,y) (x.'*y);

TRv=eye(size(V,2));
TLw=eye(size(V,2));
if jCol>1
    switch Opts.orth
        case 'dgks'
            % iterates standard gram-schmidt
            orthError=1;
            count=0;
            while(orthError>Opts.dgksTol)
                h=IP(V(:,1:jCol-1),V(:,jCol));
                V(:,jCol)=V(:,jCol)-V(:,1:jCol-1)*h;
                TRv(:,jCol)=TRv(:,jCol)-TRv(:,1:jCol-1)*h;
                if hermite_gram_sch
                    h=IP(W(:,1:jCol-1),W(:,jCol));
                    W(:,jCol)=W(:,jCol)-W(:,1:jCol-1)*h;
                    TLw(:,jCol)=TLw(:,jCol)-TLw(:,1:jCol-1)*h;
                end
                orthError=norm(IP([V(:,1:jCol-1),V(:,jCol)/sqrt(IP(V(:,jCol),V(:,jCol)))],...
                    [V(:,1:jCol-1),V(:,jCol)/sqrt(IP(V(:,jCol),V(:,jCol)))])-speye(jCol),'fro');
                if count>50 % if dgksTol is too small, Matlab can get caught in the while-loop
                    error('Orthogonalization of the Krylov basis failed due to the given accuracy.');
                end
                count=count+1;
            end
        case 'mgs'
            for iCol=1:jCol-1
              h=IP(V(:,jCol),V(:,iCol));
              V(:,jCol)=V(:,jCol)-V(:,iCol)*h;
              TRv(:,jCol)=TRv(:,jCol)-h*TRv(:,iCol);
              if hermite_gram_sch
                h=IP(W(:,jCol),W(:,iCol));
                W(:,jCol)=W(:,jCol)-W(:,iCol)*h;
                TLw(:,jCol)=TLw(:,jCol)-h*TLw(:,iCol);
              end 
            end
       case '2mgs'
            for k=0:1
                for iCol=1:jCol-1
                  h=IP(V(:,jCol),V(:,iCol));
                  V(:,jCol)=V(:,jCol)-V(:,iCol)*h;
                  TRv(:,jCol)=TRv(:,jCol)-h*TRv(:,iCol);
                  if hermite_gram_sch
                    h=IP(W(:,jCol),W(:,iCol));
                    W(:,jCol)=W(:,jCol)-W(:,iCol)*h;
                    TLw(:,jCol)=TLw(:,jCol)-h*TLw(:,iCol);
                  end 
                end
            end
        otherwise
            error('Opts.orth is invalid.');
    end  
end

% normalize new basis vector
h = sqrt(IP(V(:,jCol),V(:,jCol)));
V(:,jCol)=V(:,jCol)/h;
TRv(:,jCol) = TRv(:,jCol)/h;
if hermite_gram_sch
    h = sqrt(IP(W(:,jCol),W(:,jCol)));
    W(:,jCol)=W(:,jCol)/h;
    TLw(:,jCol) = TLw(:,jCol)/h;
end
end



end

