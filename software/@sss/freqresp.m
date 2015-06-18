
function G = freqresp(sys, s, varargin)
% Evaluates complex transfer function of LTI systems
% ------------------------------------------------------------------
% G = freqresp(sys, s, varargin)
% Inputs:       * sys: an sss-object containing the LTI system
%               * s: vector of complex frequencies
%    [optional] * number of input, number of output
% Outputs:      * G: vector of complex frequency response values
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer, Rudy Eid
% Last Change:  02 Feb 2011
% ------------------------------------------------------------------

[A,B,C,D,E] = ABCDE(sys);
m=sys.m; p=sys.p;
if nargin >= 4
    % I/O-pair selected
    if ~isempty(varargin{1}) && ~isnan(varargin{1})
        B=B(:,varargin{1});
        D=D(:,varargin{1});
        m=size(B,2);
    end
    if ~isempty(varargin{2}) && ~isnan(varargin{2})
        C=C(varargin{2},:);
        D=D(varargin{2},:);
        p=size(C,1);
    end
end

% calculate value of transfer function for given vector of frequencies
G=zeros(p,m,length(s));
for i=1:length(s)
    if isinf(s(i)) || sys.n==0
        G(:,:,i) = D;
        continue
    end
    %   [L,U]=lu(A - s(i)*E, 'vector');
    [L,U,k,l,S]=lu(A - s(i)*E, 'vector');
    warning('OFF', 'MATLAB:nearlySingularMatrix')
    b=S\B; b=b(k,:);
    x=L\b;
    x(l,:)=U\x;
    warning('ON', 'MATLAB:nearlySingularMatrix')
    G(:,:,i) = -C*x + D;
end

warning('off', 'MATLAB:mat2cell:TrailingUnityVectorArgRemoved');
G = mat2cell(G,ones(p,1),ones(m,1),length(s));
G = cellfun(@(x) x(:,:),G, 'UniformOutput', false);
    
%     else
% %         if find(E-eye(size(A)))
% %             % [AA,BB,Q,Z] = HESS(A,B) => Q*A*Z = AA,  Q*B*Z = BB.
% %             [A, E, Q, Z] = hess(A, E);
% %             C = C*Q;
% %             B = Z*B;
% %             for i=1:length(s)
% %                 [L,U] = lu(full(A - s(i)*E));
% %                 G(i) = -C/U*(L\B) + D;
% %             end
% %         else
% %             [T, A] = balance(A);
% %             B = T \ B;
% %             C = C * T;
% %             [P, A] = hess(A);
% %             C = C*P;
% %             B = P'*B;
% %             sys = struct('A', A, 'B', B, 'C', C, 'D', D, 'E', E);
% %             g = ssfresp(A,B,[],[],[],s);
% %             g = reshape(g,[length(A) numel(s)]);
% %             z = C*g  + D;
%             if ~all(isfield(sys,{'residuen','poles'}))
%                 [temp1, temp2,temp3,sys] = residuen(sys);
%                 clear('temp1','temp2','temp3')
%             end
%             res=sys.residuen;
%             p=sys.poles;
%             if ~isempty(in)
%                 res=transpose({res{:,in}});
%             end
%             if ~isempty(out)
%                 res={res{out,:}};
%             end
%             G=cellfun(@(x) transpose(sum(transpose(ones(size(s)))*x./(transpose(s)*ones(size(p))-ones(size(s))'*p),2)+D),res,'UniformOutput',false);
% 
% %            res=cell2mat(res);
% %            for i=1:length(s)
% %                 G(i) = sum( res./(s(i)-p));
% % %                 if isreal(s(i)) % remove numerical noise in imaginary part
% % %                     G(i) = real(G(i)) + D;
% % %                 end
% %             end
% %         end
% 
%     end
end
