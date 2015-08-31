function [V, W] = lanczos(E,A,b,c,s0_inp,s0_out)
% Lanczos Algorithm for order reduction of LTI SISO systems using
% multiple expansion points
% ------------------------------------------------------------------

warning('Obacht. Noch nicht 100% getestet!')

q=length(s0_inp);

% fprintf('Arnoldi   0%%');


s0 = unique([s0_inp s0_out]);

% remove one element of complex pairs
k=find(imag(s0)<0);
if ~isempty(k)
    s0c = conj(s0(k));
    s0(k) = [];
end

% preallocate memory
V=zeros(length(b),q);
W=V;
kV=0;
kW=0;

for j=1:length(s0)
%     fprintf([repmat(8,1,4) '%3.0f%%'], 100*(j-1)/q);
   
    if isinf(s0(j))
        % s0=inf, match Markov parameter instead of moment
        
        [L,U,p,o,S]=lu(E,'vector');
        temp = b;
        for k=1:length(find(s0_inp==inf))
            temp = invsolve(temp,L,U,p,o,S);
            [V,W] = GS(temp, V,W,kV,kW, E);
            kV=kV+1;
            temp = A*V(:,kV);
        end
        temp = c';
        for k=1:length(find(s0_out==inf))
            temp = invsolve(temp,U',L',o,p,S);
            [W,V] = GS(temp, W, V, kW, kV, E);
            kW=kW+1;
            temp = A*W(:,kW);
        end
        
    else
        
        [L,U,p,o,S]=lu(sparse(A-s0(j)*E),'vector');
        temp = b;
        for k=1:length(find(s0_inp==s0(j)))
            temp = invsolve(temp,L,U,p,o,S);
            [V,W] = GS(temp, V, W, kV, kW, E);
            kV=kV+1;
            temp = E*V(:,kV);
        end
        temp = c';
        for k=1:length(find(s0_out==s0(j)))
            temp = invsolve(temp,U',L',o,p,S);
            [W,V] = GS(temp, W, V, kW, kV, E);
            kW=kW+1;
            temp = E*W(:,kW);
        end
        
    end

%     % split complex conjugate columns into real (->j) and imag (->j+length(s0c)/2
%     if ~isreal(s0(j))
%         V(:,j+length(s0c)/2)=imag(temp);
%         temp=real(temp);
%     end

    %orthogonalize columns from imaginary components

end

% fprintf([repmat(8,1,4) 'done.\n']);
end

function [V,W] = GS(newcol, V, W, kV, kW, E)
    if ~exist('E', 'var')
        E=speye(size(V,1));
    end
    if kV<=kW
        for i=1:kV
            h=W(:,i)'*E*newcol;
            newcol=newcol-V(:,i)*h;
        end
    else
        for i=1:kW
            h=W(:,i)'*E*newcol;
            newcol=newcol-V(:,i)*h;
        end
    end
    
%     if kV<size(W,2)
%         newcol=-newcol/(W(:,size(V,2)+1)'*E*newcol);
%     else
%         x = newcol'*E*newcol;
%         if x<0
%             newcol=-newcol/sqrt(-x);
%         else
%             newcol=newcol/sqrt(x);
%         end
%     end
    
    if kV<kW
        xV = newcol'*E*newcol;
        xW = W(:,kV+1)'*E*W(:,kV+1);
        x = newcol'*E*W(:,kV+1) / sqrt(xV*xW);
%        x=1;
        if x<0
            V(:,kV+1)=-newcol/sqrt(-x*xV);
            W(:,kV+1)=W(:,kV+1)/sqrt(-x*xW);
        else
            V(:,kV+1)=newcol/sqrt(x*xV);
            W(:,kV+1)=W(:,kV+1)/sqrt(x*xW);
        end
    else
        V(:,kV+1)=newcol;
    end

    for i=(kV+2):kW
        h=W(:,i)'*E*V(:,kV+1);
        W(:,i)=W(:,i)-W(:,kV+1)*h;
    end


end

function newcol = invsolve(newcol,L,U,p,o,S)
    b=S\newcol;
    newcol=L\b(p,:);
    newcol(o,:)=U\newcol;
end



