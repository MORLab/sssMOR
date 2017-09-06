function [S,R,output_data] = crksm(varargin)
% CRKSM - Solve Laypunov equations with a cummulative rational Krylov subspace method
% Info: Funktionen fuer neue shifts muessen ab Zeile 470 unter der 'mess'-Option eingebunden werden 
%
%              APE' + EPA' + BB' = 0 (I)
%              AQE' + EQA' + C'C = 0 (II)
%
% Synthax MY_RKSM
%       [S,R,output_data]         = MY_RKSM(sys)
%       [S,R,output_data]         = MY_RKSM(sys, s0_inp)
%       [S,R,output_data]         = MY_RKSM(sys, s0_inp, Rt)
%       [S,R,output_data]         = MY_RKSM(sys, [], s0_out)
%       [S,R,output_data]         = MY_RKSM(sys, [], s0_out, [], Lt)
%       [S,R,output_data]         = MY_RKSM(sys, s0_inp, s0_out)
%       [S,R,output_data]         = MY_RKSM(sys, s0_inp, s0_out, Rt, Lt)
%       [S,R,output_data]         = MY_RKSM(sys,...,Opts_rksm)
%
%       [S,output_data]           = MY_RKSM(A,B,[],[],s0_inp)
%       [S,output_data]           = MY_RKSM(A,B,[],[],s0_inp,Rt)
%       [S,R,output_data]         = MY_RKSM(A,B,C, [],s0_inp)
%       [S,R,output_data]         = MY_RKSM(A,B,C, [],s0_inp,s0_out)
%       [S,R,output_data]         = MY_RKSM(A,B,C, [],s0_inp,Rt)
%       [S,R,output_data]         = MY_RKSM(A,B,C, [],s_inp,s_out,Rt,Lt)
%       [S,output_data]           = MY_RKSM(A,B,[] ,E,s0_inp)
%       [S,output_data]           = MY_RKSM(A,B,[] ,E,s0_inp,Rt)
%       [S,R,output_data]         = MY_RKSM(A,B,C,E,s0_inp)
%       [S,R,output_data]         = MY_RKSM(A,B,C,E,s0_inp,s0_out)
%       [S,R,output_data]         = MY_RKSM(A,B,C,E,s0_inp,s0_out,Rt,Lt)
%       [S,R,output_data]         = MY_RKSM(A,B,...,s_inp,...,Opts_rksm) 
%
% Input:
%       - system matrices A, B, C, and E or a sys-object
%       - intial shift vector s_inp (must have at least two entries)
%       - tangential directions Rt and Lt in the MIMO case
%       - Opts-struct
%
% Output:
%       - Cholseky factors S, R of the solutions P, Q of the linear
%         Lyapunov equations
%       - data-Struct output_data containing:
%         V_basis, W_basis, norm values of the iterations norm_val,
%         additional shifts, last reduced system Ar, Br, ...
%         last rhs, reference residual norm res0, last norm value
%
%
% possible Options:
% - Opts.shifts:      choose which shifts and how they should be used
%                     - cyclic: use the same shifts for the whole iteration
%                     - adaptive: a new shift is computed online for every
%                       single iteration (only for SISO-systems)
%                     - mess: use a function to generate shifts out of the
%                       mess-toolbox
%                       ['cyclic' / 'adaptive' / 'mess']
%
% - Opts.reduction:   specify the sort of reduction 
%                     - onesided: get reduced system by only using the V-matrix
%                     - twosided: get reduced system by using V an W matrices
%                       ['onesided' / 'twosided']
%
% - Opts.residual:    specify determination criteria
%                     - residual_lyap: compute residual of Lyapunov equation
%                     - norm_chol: compare  the norm of the two last Cholesky
%                       factors
%                       ['residual_lyap' / 'norm_chol']
%
% - Opts.rksmnorm:    specify norm
%                     - 'H2': use 2-norm (Euclidian Norm)
%                     - 'fro': use Frobenius Norm
%
% - Opts.lowrank:     compute the low rank factor of the final solution
%                     [0 / 1]
%
% - Opts.orth:        specify orhtogonalization method in Gram-Schmidt
%                     ['2mgs' / 'dgks' / 'mgs']
%
% - Opts.tol:         specify tolerance
%
% - Opts.maxiter:     specify maximal number of iterations
%
% - Opts.info_rksm:   shows status information during the function
%                     call of rksm [0 / 1]

%% Still To Come
% - extended case nach moeglichkeit aufloesen
% - extended input f?r matrizen machen

%% Create Def-struct containing default values and make global settings
% Note: there may be no point named Def.reuseLU, otherwise one gets a conflict with solveLse/lyapchol/bilyapchol

Def.crksmMethod          = 'rational';         % ['rational' / 'extended']
Def.shifts               = 'cyclic';           % ['cyclic' / 'adaptive' / 'mess']
Def.reduction            = 'onesided';         % ['onesided' / 'twosided']
Def.residual             = 'residual_lyap';    % ['residual_lyap' / 'norm_chol']
Def.rksmnorm             = 'H2';               % ['H2' / 'fro']
Def.lowrank              =  0;                 % [0 / 1]
Def.orth                 = '2mgs';             % for Gramschmidt
Def.rctol                =  1e-12;             % default tolerance
Def.maxiter_rksm         =  200;               % default number of iterations
Def.info_rksm            =  0;                 % show programme status information, [0 / 1]
Def.lse                  = 'sparse';           % for solveLse function
Def.real                 = true;

% global variables
global hermite_gram_sch  

%% Parsing of Inputs

if isa(varargin{end},'struct') 
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

output_data = struct('V_basis',[], 'W_basis',[], 'norm_val',[],'shifts',[],...
                     'Ar',[], 'Br',[], 'Er',[], 'Cr',[], 'rhsb',[], 'res0',[], 'lastnorm',[]);

if length(varargin) == 1 || (ismatrix(varargin{1}) && length(varargin) <= 4)
    
    s0_inp = [0 inf 0 inf];              s0_out = []; 
    Rt     = [];                         Lt     = [];
    if strcmp(Opts.reduction,'twosided'), s0_out = s0_inp; end
else
    % input of sys-object
    if isa(varargin{1},'ss') || isa(varargin{1},'sss') || isa(varargin{1},'ssRed')

        % read in sys object
        sys = varargin{1};

        % use hess if sys is ssRed object
        if isa(sys,'ssRed')
            Opts.lse='hess'; 
            if isempty(sys.E), sys.E = eye(size(sys.A)); end %ssRed robust compatibility
        end

        % set system matrices
        A = sys.A;       n = size(A,1);      B = sys.B;       m = size(B,2);
        C = sys.C;       p = size(C,1);      E = sys.E; 

        % check usage and inputs
        if length(varargin) == 2                            % usage: CRKSM(sys, s0_inp)
            s0_inp = varargin{2};                s0_out = [];
            Rt     = [];                         Lt     = [];
        elseif length(varargin) == 3       
            if size(varargin{3},1) == m &&...
               size(varargin{3},2) == size(varargin{2},2)   % usage: CRKSM(sys, s0_inp, Rt)
                s0_inp = varargin{2};            s0_out = [];
                Rt     = varargin{3};            Lt     = [];
            elseif isempty(varargin{2})                     % usage: CRKSM(sys, [], s0_out)
                s0_inp = [];                     s0_out = varargin{3};
                Rt     = [];                     Lt     = [];
            else                                            % usage: CRKSM(sys, s0_inp, s0_out)
                s0_inp = varargin{2};            s0_out = varargin{3};
                Rt     = [];                     Lt     = [];
            end 
        elseif length(varargin) == 5
            if isempty(varargin{2})                         % usage: CRKSM(sys, [], s0_out, [], Lt)
                s0_inp = [];                     s0_out = varargin{3};
                Rt     = [];                     Lt     = varargin{5};         
            else                                            % usage: CRKSM(sys, s0_inp, s0_out, Rt, Lt)
                s0_inp = varargin{2};            s0_out = varargin{3};
                Rt     = varargin{4};            Lt     = varargin{5};
            end  
        else
            error('Input not compatible with current crksm implementation');
        end

    % input of matrices    
    elseif length(varargin) > 1

        % set A, B matrix
        A = varargin{1};         n = size(A,1);
        B = varargin{2};         m = size(B,2);

        % set other inputs
        if nargin == 5
            if isempty(varargin{3}) && isempty(varargin{4})
                C      = [];               E      = speye(size(A));
                s0_inp = varargin{5};      s0_out = [];
                Rt     = [];               Lt     = [];
            elseif isempty(varargin{4})
                C      = varargin{3};      E      = speye(size(A));
                s0_inp = varargin{5};      s0_out = [];
                Rt     = [];               Lt     = [];
            elseif isempty(varargin{3})
                C      = [];               E      = speye(size(A));
                s0_inp = varargin{5};      s0_out = [];
                Rt     = [];               Lt     = [];   
            else
                C      = varargin{3};      E      = speye(size(A));
                s0_inp = varargin{5};      s0_out = [];
                Rt     = [];               Lt     = [];   
            end 
        elseif nargin == 6
            if isempty(varargin{3}) && isempty(varargin{4})
                C      = [];               E      = speye(size(A));
                s0_inp = varargin{5};      s0_out = [];
                Rt     = varargin{6};      Lt     = []; 
            elseif isempty(varargin{3})
                C      = [];               E      = varargin{4};
                s0_inp = varargin{5};      s0_out = [];
                Rt     = varargin{6};      Lt     = []; 
            elseif isempty(varargin{4}) && size(varargin{6,1}) == m
                C      = varargin{3};      E      = speye(size(A));
                s0_inp = varargin{5};      s0_out = [];
                Rt     = varargin{6};      Lt     = [];  
            elseif isempty(varargin{4}) && size(varargin{6,1}) ~= m
                C      = varargin{3};      E      = speye(size(A));
                s0_inp = varargin{5};      s0_out = varargin{6};
                Rt     = [];               Lt     = []; 
            else
                C      = varargin{3};      E      = varargin{4};
                s0_inp = varargin{5};      s0_out = varargin{6};
                Rt     = [];               Lt     = []; 
            end   
        elseif nargin == 8
            if isempty(varargin{4})
                C      = varargin{3};      E      = speye(size(A));
                s0_inp = varargin{5};      s0_out = varargin{6};
                Rt     = varargin{7};      Lt     = varargin{8}; 
            else
                C      = varargin{3};      E      = varargin{4};
                s0_inp = varargin{5};      s0_out = varargin{6};
                Rt     = varargin{7};      Lt     = varargin{8};
            end    
        else
            if length(varargin) > 8 || length(varargin) < 3
                error('Wrong input');
            end   
        end
    end

    % check settings 
    if ~isempty(s0_out) && strcmp(Opts.reduction,'onesided'),  Opts.reduction = 'twosided';  end
    if isempty(s0_out) && strcmp(Opts.reduction,'twosided'),           s0_out = s0_inp;      end
    if ~isempty(s0_out) && ~exist('p','var'),                               p = size(C,1);   end

    % check shifts and tangential directions, make column vectors, matrices
    [s0_inp,s0_out,Rt,Lt] = parseShifts(s0_inp,s0_out,Rt,Lt,n,m,p,Opts);

    % preallocate memory for V,W
    % V = zeros(n,Opts.maxiter_rksm);
    % if ~isempty(s0_out),     W = zeros(n,Opts.maxiter_rksm);    end

    % clear arguements
    clear Def sys varargin
end
   
%% RKSM Method
switch Opts.rksm_method
    case 'rational'
        % build the first basis with arnoldi, containing two directions or three if second shift is complex
        if strcmp(Opts.reduction,'onesided') == true
            % build first input subspace (onesided), set shift variable jCol_inp, for SISO
            if (s0_inp(1,1) == conj(s0_inp(1,2)))
                if isempty(Rt)
                    [V] = arnoldi(E,A,B,s0_inp(1,1:2));
                    vnew = zeros(n,nnz(V(1,:))/2);
                else
                    [V] = arnoldi(E,A,B,s0_inp(1,1:2),Rt(:,1:2));
                    vnew = zeros(n,nnz(V(1,:))/2);
                end       
            else
                if ~isempty(Rt)
                    [V] = arnoldi(E,A,B,s0_inp(1,1),Rt(:,1));
                    vnew = zeros(n,nnz(V(1,:))/2);
                else
                    [V] = arnoldi(E,A,B,s0_inp(1,1));
                    vnew = zeros(n,nnz(V(1,:))/2);
                end
            end    
        % build first input and output subspace (twosided, hermite)
        elseif strcmp(Opts.reduction,'twosided') == true && isempty(s0_out)
            if (s0_inp(1,1) == conj(s0_inp(1,2)))
                if isempty(Rt) && isempty(Lt)
                    [V,~,~,W] = arnoldi(E,A,B,C,s0_inp(1,1:2));
                    vnew = zeros(n,nnz(V(1,:))/2);
                else
                    [V,~,~,W] = arnoldi(E,A,B,C,s0_inp(1:2,1),Rt(:,1:2),Lt(:,1:2));
                    vnew = zeros(n,nnz(V(1,:))/2);
                end
            else
                if isempty(Rt) && isempty(Lt)
                    [V,~,~,W] = arnoldi(E,A,B,C,s0_inp(1,1));
                    vnew = zeros(n,nnz(V(1,:))/2);
                else
                    [V,~,~,W] = arnoldi(E,A,B,C,s0_inp(1,1),Rt(:,1),Lt(:,1));
                    vnew = zeros(n,nnz(V(1,:))/2);
                end
            end           
        % build first input and output subspace (twosided, not hermite) 
        elseif strcmp(Opts.reduction,'twosided') == true && ~isempty(s0_out)
            if (s0_inp(1,1) == conj(s0_inp(1,2)))
                if isempty(Rt) && isempty(Lt)
                    [V] = arnoldi(E,A,B,s0_inp(1,1:2));
                    vnew = zeros(n,size(V,2)/2);
                    [~,~,~,W] = arnoldi(E,A,B,C,s0_out(1,1:2));
                else
                    [V] = arnoldi(E,A,B,s0_inp(1:2,1),Rt(:,1:2));
                    vnew = zeros(n,size(V,2)/2);
                    [~,~,~,W] = arnoldi(E,A,B,C,s0_out(1,1:2),Rt(:,1:2),Lt(:,1:2));
                end
            else
                if isempty(Rt) && isempty(Lt)
                    [V] = arnoldi(E,A,B,s0_inp(1,1));
                    vnew = zeros(n,size(V,2));
                    [~,~,~,W] = arnoldi(E,A,B,C,s0_out(1,1));
                else
                    [V] = arnoldi(E,A,B,s0_inp(1,1),Rt(:,1),Lt(:,1));
                    vnew = zeros(n,size(V,2));
                    [~,~,~,W] = arnoldi(E,A,B,C,s0_out(1,1),Rt(:,1),Lt(:,1));
                end
            end
        end
        
        % calculating the first reduced system % der teil arbeitet schon
        % mit real subsoace von arnoldi
        if strcmp(Opts.reduction,'onesided')
            Ar = V'*A*V;
            Br = V'*B;
            Er = V'*E*V;
        else
            Ar = W'*A*V;
            Br = W'*B;
            Er = W'*E*V;
            Cr = C*V;
        end
                        
        % start iteration
        for ii = (size(V,2)/size(vnew,2)):1:Opts.maxiter_rksm
            if exist('S','var')
                if size(V,2) == (ii-1)*size(vnew,2)
                    
                    % start again with the first shift in s0_inp/s0_out
                    if ii > size(s0_inp,2)
                        % get new shifts
                        if strcmp(Opts.shifts,'cyclic')
                            % enlarge shifts, tangential directions etc. 
                            s0_inp = repmat(s0_inp,1,2);
                            if ~isempty(s0_out),    s0_out = repmat(s0_out,1,2); end
                            if ~isempty(Rt),        Rt     = repmat(Rt,1,2);     end
                            if ~isempty(Lt),        Lt     = repmat(Lt,1,2);     end
                        elseif strcmp(Opts.shifts,'adaptive')
                        % hier adaptive shifts
                        elseif strcmp(Opts.shifts,'mess')
                        % hier irgendeine andere Funktion um shifts zu berechnen
                        end
                    end

                    % set the rhs for the diffrent input cases
                    % case: onesided, block
                    if isempty(s0_out) && isempty(Rt)
                        rhsB = E*V(:,size(V,2)-(size(vnew,2)-1):size(V,2));
                    % case: twosided, hermite, block
                    elseif isempty(Rt) && isequal(s0_inp,s0_out)
                        rhsB = E*V(:,size(V,2)-(size(vnew,2)-1):size(V,2));
                        rhsC = E*W(:,size(W,2)-(size(wnew,2)-1):size(W,2));
                    % case: onesided, tangential
                    elseif ~isempty(Rt) && isempty(s0_out)
                        rhsB = B*Rt(:,ii);
                        %rhsB = E*V(:,size(V,2)-(size(vnew,2)-1):size(V,2))*s_inpRt;
                    % case: twosided, tangential, not hermite
                    elseif ~isempty(Rt) && ~isempty(Lt)
                        rhsB = B * Rt(:,ii);
                        rhsC = C' * Lt(:,ii);
                        %rhsB = E*V(:,size(V,2)-(size(vnew,2)-1):size(V,2))*Rt(:,ii);
                        %rhsC = E*W(:,size(W,2)-(size(wnew,2)-1):size(W,2))*Lt(:,ii); 
                    end
                    
                    % A matrix for solveLse function: (A-s_inp*E), (A-s_out*E)'; get new directions for V,W
                    if s0_inp(1,ii) ~= s0_inp(1,ii-1) && s0_inp(1,ii) ~= conj(s0_inp(1,ii-1)) 
                        Anew_V = (A-s0_inp(1,ii)*E);
                        [vnew] = solveLse(Anew_V,rhsB);
                    else
                        Opts.reuseLU = 1;
                        [vnew] = solveLse(Anew_V,E*rhsB,Opts);
                    end 
                    if strcmp(Opts.reduction,'twosided') && s0_out(1,ii) ~= s0_out(1,ii-1) && s0_out(1,ii) ~= conj(s0_out(1,ii-1))   
                        Anew_W = (A-s0_out(ii,1)*E)'; 
                        [wnew] = solveLse(Anew_W,rhsC);
                    elseif strcmp(Opts.reduction,'twosided') && s0_out(1,ii) == s0_out(1,ii-1) && s0_out(1,ii) == conj(s0_out(1,ii-1)) 
                         Opts.reuseLU = 1;
                         [wnew] = solveLse(Anew_W,E'*rhsC,Opts);
                    end
                    
                    % build and orthogonalize new basis
                   if strcmp(Opts.reduction,'twosided')
                       V = [V real(vnew)];
                       W = [W real(wnew)];
                       %V(:,nnz(V(1,:))+1:nnz(V(1,:))+size(vnew,2)) = real(vnew);
                       %W(:,nnz(W(1,:))+1:nnz(W(1,:))+size(wnew,2)) = real(wnew);
                   else % onesided case
                       V = [V real(vnew)];
                       %V(:,nnz(V(1,:))+1:nnz(V(1,:))+size(vnew,2)) = real(vnew);
                   end
               end % end of complex wait-sequence (if isreal(vnew))

               % orthogonalize new basis; note: for loop is neccessary for block
               for jj=size(V,2)-(size(vnew,2)-1):1:size(V,2)
                   if strcmp(Opts.reduction,'twosided')
                       hermite_gram_sch = 1;
                       [V,~,W] = gramSchmidt(jj,V,W,Opts);
                   else
                        V = gramSchmidt(jj,V,Opts);
                   end
               end
               
               % reduction step
               if strcmp(Opts.reduction,'twosided')
                   [Ar,Er,Br,Cr] = reduction(A,B,E,Ar,Br,Er,size(V,2)-(size(vnew,2)-1),size(V,2),V,W,C,Cr);
               else
                   [Ar,Er,Br] = reduction(A,B,E,Ar,Br,Er,size(V,2)-(size(vnew,2)-1),size(V,2),V);
               end
           end % end of if ii> size(V,2)
           
           
           
           % save last Cholesky factor of solution for norm_chol
           if strcmp(Opts.residual,'norm_chol') && ii > 2,  S_last = S;  end
                   
           % try to solve Lyapunov equation first with lyapchol or second with lyap
           try
               S = lyapchol(Ar,Br,Er);
               if nargout >= 2 && strcmp(Opts.reduction,'twosided')
                   R = lyapchol(Ar',Cr',Er');
               end
           catch
               S = lyap(Ar,Br*Br',[],Er);
               if nargout >= 2 && strcmp(Opts.reduction,'twosided')
                   R = lyap(Ar',Cr'*Cr,[],Er');
               end
           end
                      
           % choose computation of norm/stopping criteria
           if strcmp(Opts.residual,'residual_lyap') || ii <= 2
               % test determination (computation of residual after Panzer/Wolff)
               % compute Er^-1*Br and Er^-1*Ar and other factors
               Er_inv_Br = solveLse(Er,Br);
               Opts.reuseLU = 1;
               Er_inv_Ar = solveLse(Er,Ar,Opts);
               AV = A*V;        EV = E*V;
               %Pr = S'*S;

               % compute factors from residual
               Borth = B-EV*Er_inv_Br;
               %B_s     = B-EV*(Er\Br);

               Borth2 = Borth'*Borth;
               %BstBs   = B_s'*B_s;
               
               Cr_hat_rhs = Borth'*(AV-EV*Er_inv_Ar);
               Cr_hat = solveLse(Borth2,Cr_hat_rhs);
               %F = E*V*(Er_inv_Br+(S*S')*Cr_hat');
               F = E*V*(Er_inv_Br+(S'*S)*Cr_hat');
 
               % compute residual norm (Euclidean Norm)
               if strcmp(Opts.rksmnorm, 'H2')
                   res0  = norm(Br' * Br,2);
                   output_data.res0(1,ii) = res0;
                   Rnorm = max(abs(eig(full([Borth2+Borth'*F, Borth2; F'*Borth+F'*F, F'*Borth])))) / res0; 
               else
                   % Frobenius Norm
                   res0  = norm(Br' * Br,'fro');
                   output_data.res0(1,ii) = res0;
                   Rnorm = sqrt(sum(eig(full([Borth2+Borth'*F, Borth2; F'*Borth+F'*F, F'*Borth])))^2) / res0;
               end
           else
               if strcmp(Opts.rksmnorm, 'H2')
                   %X_lastnorm = norm(S_last);
                   X_lastnorm = NormFrobEfficient(1,S_last);
                   %X_norm = norm(S);
                   X_norm = NormFrobEfficient(1,S);
               else
                   X_lastnorm = NormFrobEfficient(0,S_last);
                   X_norm = NormFrobEfficient(0,S);
               end
               output_data.lastnorm(1,ii)=X_norm;  
               output_data.lastnorm(2,ii)=X_lastnorm ;
               Rnorm = abs(X_norm-X_lastnorm);
           end
           output_data.norm_val(1,ii) = Rnorm;
           
           % show status information
           if mod(ii,2) == 0 && Opts.info_rksm == 1
               fprintf('\nRKSM step: \t %d',ii);
               if strcmp(Opts.residual,'residual_lyap') 
                   fprintf('\t current residual norm: \t %d',Rnorm);
               else
                   fprintf('\t relative change in S: \t %d',Rnorm);
               end
           end
          
           % stop program
           if Rnorm < Opts.rctol
               chol_rank_S = rank(S);
               if ~exist('R','var')
                   fprintf('\n RKSM step: %d \t Convergence\n last residual norm: %d, tolerance: %d,\n rank of Cholesky factorS: %d',ii,Rnorm,Opts.rctol,chol_rank_S);
               else
                   chol_rank_R = rank(R);
                   fprintf('\n RKSM step \t %d \t Convergence, last residual: \t %d, tolerance: \t %d,\n rank of Cholesky factor S: \t %d\n rank of Cholesky factor R: \t %d',ii,Rnorm,Opts.rctol,chol_rank_S,chol_rank_R);
               end
                    break;
           elseif size(S,2) == n
               disp('\n V has reached the dimension of the original System without converging!');
               disp('\n If you want to continue with the iterations until the maximal number of iterations is reached press 1!');
               disp('\n If you want to stop the program and see the current output arguements press 2!');
               decision = input('\n\n press one or two\n');
               if decision == 2
                   break;
               end
           end
           
           % enlarge subspace with the imaginary part of the direction and eventually calculate new real directions
           if ~isreal(vnew)
               if strcmp(Opts.reduction,'twosided')
                   if ~isreal(vnew) && ~isreal(wnew)
                       V = [V imag(vnew)];
                       W = [W imag(wnew)];
                       %V(:,nnz(V(1,:))+1:nnz(V(1,:))+size(vnew,2)) = imag(vnew);
                       %W(:,nnz(W(1,:))+1:nnz(W(1,:))+size(wnew,2)) = imag(wnew);
                   elseif ~isreal(vnew) && isreal(wnew)
                       V = [V imag(vnew)];
                       %V(:,nnz(V(1,:))+1:nnz(V(1,:))+size(vnew,2)) = imag(vnew);
                       
                       % calculate new real direction, if last vnew is real and last wnew is complex       
                       if isempty(Lt)
                           rhsC = E*W(:,size(W,2)-size(wnew,2)+1:size(W,2));
                       else
                           rhsC = B*Lt(:,ii);
                           %rhsC = E*W(:,size(W,2)-size(wnew,2)+1:size(W,2))*Lt(:,ii);
                       end
                       
                       if s0_out(1,ii) ~= s0_out(1,ii-1)
                            Anew_W = (A-s0_out(1,ii)*E)'; 
                            [wnew] = solveLse(Anew_W,rhsC);
                       else
                           Opts.reuseLU = 1;
                           [wnew] = solveLse(Anew_W,E'*rhsC,Opts);
                       end
                       W = [W real(wnew)];
                       %W(:,nnz(W(1,:))+1:nnz(W(1,:))+size(wnew,2)) = real(wnew);
           
                   elseif isreal(vnew) && ~isreal(wnew)
                       W = [W imag(wnew)];
                       %W(:,nnz(W(1,:))+1:nnz(W(1,:))+size(wnew,2)) = imag(wnew);
                       % calculate new real direction, if last vnew is real and last wnew is complex
                       if isempty(Rt)
                           rhsB = E*V(:,size(V,2)-size(vnew,2)+1:size(V,2));
                       else
                           rhsB = B*Rt(:,ii);
                           %rhsB = E*V(:,size(V,2)-size(vnew,2)+1:size(V,2))*Rt(:,ii);
                       end
                       if s0_inp(1,ii) ~= s0_inp(1,ii-1)
                           Anew_V = (A-s0_inp(1,ii)*E);
                           [vnew] = solveLse(Anew_V,rhsB);
                       else
                           Opts.reuseLU = 1;
                           [vnew] = solveLse(Anew_V,E*rhsB,Opts);
                       end
                       V(:,nnz(V(1,:))+1:nnz(V(1,:))+size(vnew,2)) = real(vnew);
                   end
               else % onesided case
                   V = [V imag(vnew)];
                   %V(:,nnz(V(1,:))+1:nnz(V(1,:))+size(vnew,2)) = imag(vnew);
               end
               % make vnew real
               vnew = zeros(n,size(vnew,2));
           end
        end  % end of for loop
        if ii == Opts.maxiter_rksm  && Rnorm > Opts.rctol
            warning('\n maximum number of iterations is reached without converging!' )
        end
             
               
    case 'extended' 
        % EKSM method 
        % the shifts for building the Krylov subsubspace are fixed at 0 and
        % inf (Markov parameter). The subspace is built in the following
        % way:
        % for V we have
        %       1. subspace: span{E^-1*B}
        %       2. subspace: span{E^-1*B, A^-1*B}
        %       3. subspace: span{E^-1*B, A^-1*B, E^-1*A*E^-1*B}
        %       4. subspace: span{E^-1*B, A^-1*B, E^-1*A*E^-1*B, A^-1*E*A^-1*B} and so on
        % for W we have
        %       1. subspace: span{E^-T*C_T}
        %       2. subspace: span{E^-T*C_T, A^-T*C_T}
        %       3. subspace: span{E^-T*C_T, A^-T*C_T, E^-T*A_T*E^-T*C_T}
        %       4. subspace: span{E^-T*C_T, A^-T*C_T, E^-T*A_T*E^-T*C_T, A^-T*E_T*A^-T*C_T} and so on
        
        % start iteration 
        [vnew2] = solveLse(A,B);             % A^-1*B
        [vnew1] = solveLse(E,B);             % E^-1*B
        
        if withoutC == 0
            [wnew2] = solveLse(A',C');       % A^-T*C_T
            [wnew1] = solveLse(E',C');       % E^-T*C_T
        end
        
        for ii = 1:1:maxiter
            % build first/second subspace
            if ii == 1 
                V = vnew1;
                jbasis    = (ii-m);                  % index for the first column of the directions calculated before
                jnew      = jbasis+m;                % index for the first new column
                jnew_last = jnew+m-1;                % index of the last new column
                for jj=jnew:1:jnew_last
                    V = gramSchmidt(jj,V,Opts); 
                end
                if withoutC == 0
                    W = wnew1;
                    for jj=jnew:1:jnew_last
                        V = gramSchmidt(jj,W,Opts); 
                    end
                end
            elseif ii == 2
                V = [V vnew2];
                jbasis    = (ii-1)*m+-(m-1);
                jnew      = jbasis+m;    
                jnew_last = jnew+m-1;
                for jj=jnew:1:jnew_last
                    V = gramSchmidt(jj,V,Opts); 
                end
                if withoutC == 0
                    W = [W wnew2];
                    for jj=jnew:1:jnew_last
                        V = gramSchmidt(jj,W,Opts); 
                    end
                end
            end
                
            % build higher subspaces
            if ii == 3
                % calculate E^1*A and A^1*E
                if withoutE == 1
                    E_invA = A;
                    A_invE = solveLse(A,eye(size(A)));
                    if withoutC == 0
                        E_invTA_T = A';
                        A_invTE_T = solveLse(A',eye(size(A)));
                    end
                else
                    E_invA = solveLse(E,A);
                    A_invE = solveLse(A,E);
                    if withoutC == 0
                        E_invTA_T = solveLse(E',A');
                        A_invTE_T = solveLse(A',E');
                    end
                end
            end
            
            if mod(ii,2) ~= 0 && ii > 2
                vnew1 = E_invA*vnew1;
                V = [V vnew1];
                jbasis    = (ii-1)*m-(m-1);
                jnew      = jbasis+m;    
                jnew_last = jnew+m-1;
                for jj=jnew:1:jnew_last
                    V = gramSchmidt(jj,V,Opts);
                end
                if withoutC == 0
                    wnew1 = E_invTA_T*wnew1;
                    W = [W wnew1];
                    for jj=jnew:1:jnew_last
                        V = gramSchmidt(jj,W,Opts);
                    end
                end
            elseif mod(ii,2) == 0 && ii > 2
                vnew2 = A_invE*vnew2;
                V = [V vnew2];
                jbasis    = (ii-1)*m-(m-1);
                jnew      = jbasis+m;    
                jnew_last = jnew+m-1;
                for jj=jnew:1:jnew_last
                    V = gramSchmidt(jj,V,Opts);
                end
                if withoutC == 0
                    wnew2 = A_invTE_T*wnew2;
                    W = [W wnew2];
                    for jj=jnew:1:jnew_last
                        V = gramSchmidt(jj,W,Opts);
                    end
                end
            end
            
            % reduce system
            if ii == 1
                if withoutC == 1
                    Ar = V'*A*V;
                    Br = V'*B;
                    Er = V'*E*V;
                else
                    Ar = W'*A*V;
                    Br = W'*B;
                    Er = W'*E*V;
                    Cr = C*V;
                end
            else   
                if withoutC == 1
                   [Ar,Er,Br] = reduction(A,B,E,Ar,Br,Er,jnew,jnew_last,V);
               else
                   [Ar,Er,Br,Cr] = reduction(A,B,E,Ar,Br,Er,jnew,jnew_last,V,W,C,Cr);
               end
            end
            
            % save last Cholesky factor of solution for norm_chol
            if strcmp(Opts.residual,'norm_chol')
                if ii > 2,  S_last = S; end
            end
            
            % try to solve Lyapunov equation first with lyapchol or second with lyap
            try
                S = lyapchol(Ar,Br,Er);
                if nargout >= 2 && withoutC == 0
                   R = lyapchol(Ar',Cr',Er');
               end
            catch
                S = lyap(Ar,Br*Br',[],Er);
               if nargout >= 2 && withoutC == 0 
                   R = lyap(Ar',Cr'*Cr,[],Er');
               end
            end
            
            % test determination (computation of residual after
            % Panzer/Wolff) or norm of Cholesky factor
            if strcmp(Opts.residual,'residual_lyap') || ii <= 2
               % test determination (computation of residual after Panzer/Wolff)
               % compute Er^-1*Br and Er^-1*Ar and other factors
               Er_inv_Br = solveLse(Er,Br);
               Opts.reuseLU = 1;
               Er_inv_Ar = solveLse(Er,Ar,Opts);
               AV = A*V;        EV = E*V;

               % compute factors from residual
               Borth = B-EV*Er_inv_Br;
               Borth2 = Borth'*Borth;
               Cr_hat_rhs = Borth'*(AV-EV*Er_inv_Ar);
               Cr_hat = solveLse(Borth2,Cr_hat_rhs);
               F = E*V*(Er_inv_Br+(S'*S)*Cr_hat');
 
               % compute residual norm (Euclidean Norm)
               if strcmp(Opts.rksmnorm, 'H2')
                   res0  = norm(Br' * Br,2);
                   output_data.res0(1,ii) = res0;
                   Rnorm = max(abs(eig(full([Borth2+Borth'*F, Borth2; F'*Borth+F'*F, F'*Borth])))) / res0; 
               else
                   % Frobenius Norm
                   res0  = norm(Br' * Br,'fro');
                   output_data.res0(1,ii) = res0;
                   Rnorm = sqrt(sum(eig(full([Borth2+Borth'*F, Borth2; F'*Borth+F'*F, F'*Borth])))^2) / res0;
               end
           else
               if strcmp(Opts.rksmnorm, 'H2')
                   X_lastnorm = NormFrobEfficient(1,S_last);
                   X_norm = NormFrobEfficient(1,S);
               else
                   X_lastnorm = NormFrobEfficient(0,S_last);
                   X_norm = NormFrobEfficient(0,S);
               end
                Rnorm = abs(X_norm-X_lastnorm);
           end
           output_data.norm_val(1,ii) = Rnorm;
           
           % show status information
           if mod(ii,2) == 0 && Opts.info_rksm == 1 
               fprintf('\nRKSM step: \t %d',ii);
               if strcmp(Opts.residual,'residual_lyap') 
                   fprintf('\t current residual norm: \t %d',Rnorm);
               else
                   fprintf('\t relative change in S: \t %d',Rnorm);
               end
           end
            
            if Rnorm < Opts.rctol
              break; 
            end 
            
        end
end

% create output
% compute low rank factor of solution
if Opts.lowrank == 1
    S = V*S;    
    if exist('R','var'),  R = V*R;  end
end

% fill output struct
output_data.V_basis = V;
output_data.Ar = Ar;
output_data.Br = Br;
output_data.Er = Er;
if exist('W','var')
    output_data.W_basis = W;
    output_data.Cr = Cr;
end

if ~exist('R','var')
    R = [];
end

% clear global variables
clearvars hermite_gram_sch 
end


% *************************************************************************
%% ***************************** AUXILIARY ********************************
% *************************************************************************
function [s0_inp,s0_out,Rt,Lt] = parseShifts(s0_inp,s0_out,Rt,Lt,n,m,p,Opts)

if  (~exist('s0_inp', 'var') || isempty(s0_inp)) && ...
    (~exist('s0_out', 'var') || isempty(s0_out))
    error('sssMOR:rk:NoExpansionPoints','No expansion points assigned.');
end

% sort expansion points & tangential directions
s0_inp = shiftVec(s0_inp);
s0old = s0_inp;
if Opts.real, 
    s0_inp = cplxpair(s0_inp);  %make sure shifts can be paired 
else
    s0_inp = sort(s0_inp);
end

if ~isempty(Rt)
    if size(Rt,2) ~= length(s0_inp), error('Inconsistent size of Rt');end
    if size(Rt,1) ~= m,              error('Inconsistent size of Lt');end
    [~,cplxSorting] = ismember(s0_inp,s0old); 
    Rt = Rt(:,cplxSorting);

    % hiermit wird gepueft, ob tang. richtungen konjugiert komplex sind, unabh?ngig von den shifts  
    %if sum(sum(imag(Rt),2)) ~= 0, error('wrong input'); end

    % hier wird gepr?ft, ob ein kompl. konjugierter shift auch eine
    % komplex conjugierte tang. Richtung hat
    if size(find(imag(s0_inp)),2) ~= size(find(imag(Rt(1,:))),2)
        error('wrong input');
    end
end
clear s0old

if ~isempty(s0_out)
    % sort expansion points & tangential directions
    s0_out = shiftVec(s0_out);
    s0old = s0_out;
    if Opts.real, 
        s0_out = cplxpair(s0_out); %make sure shifts can be paired 
    else
        s0_out = sort(s0_out);
    end

    if ~isempty(Lt)
        if size(Lt,2) ~= length(s0_out), error('Inconsistent size of Lt');end
        if size(Lt,1) ~= p,              error('Inconsistent size of Lt');end

        [~,cplxSorting] = ismember(s0_out,s0old); 
        Lt = Lt(:,cplxSorting);
        % hiermit wird gepueft, ob tang. richtungen konjugiert komplex sind, unabh?ngig von den shifts  
        %if sum(sum(imag(Lt),2)) ~= 0, error('wrong input'); end

        % hier wird gepr?ft, ob ein kompl. konjugierter shift auch eine
        % komplex conjugierte tang. Richtung hat
        if size(find(imag(s0_out)),2) ~= size(find(imag(Lt(1,:))),2)
            error('wrong input');
        end
    end
end
if length(s0_inp) > n || length(s0_out) > n
    error('sssMOR:arnoldi:reducedOrderExceedsOriginal','The desired reduced order exceeds the original order');
end

if ~isempty(s0_inp) && ~isempty(s0_out)
    % check if number of input/output expansion points matches
    if length(s0_inp) ~= length(s0_inp)
        error('Inconsistent length of expansion point vectors.');
    end
end
end

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

function [Ar,Er,Br,varargout] = reduction(A,B,E,Ar,Br,Er,jnew,jnew_last,V,varargin)
% Description
% function for cheaply computung the reduced system matrices
% Input:  - V,W (Krylov subspace bases, W optional if twosided)
%         - A,B,C,E system matrices
%         - Ar,Br,Cr,Er reduced system matrices from the iterations before
%         - jnew, jnew_last (indices for calculating the reduced system)
%
% Synthax REDUCTION:
%       [Ar,Er,Br]               = REDUCTION(A,B,E,Ar,Br,Er,jnew,jnew_last,V) 
%       [Ar,Er,Br,Cr]            = REDUCTION(A,B,E,Ar,Br,Er,jnew,jnew_last,V,W,C,Cr)
%
% Output: - Ar,Br,Er,Cr (reduced system)

if length(varargin) > 1
    W = varargin{1};    C = varargin{2};    Cr = varargin{3};
    Ar(1:jnew-1,jnew:jnew_last) = W(:,1:jnew-1)'*A*V(:,jnew:jnew_last);              % 1. step: new columns to the existing matrix Ar      
    Ar(jnew:jnew_last,1:jnew_last) = W(:,jnew:jnew_last)'*A*V(:,1:jnew_last);        % 2. step: rows are filled up so that a reduced square matrix results
    Er(1:jnew-1,jnew:jnew_last) = W(:,1:jnew-1)'*E*V(:,jnew:jnew_last);
    Er(jnew:jnew_last,1:jnew_last) = W(:,jnew:jnew_last)'*E*V(:,1:jnew_last);
    Br(jnew:jnew_last,:) = W(:,jnew:jnew_last)'*B;
    Cr(:,jnew:jnew_last) = C*V(:,jnew:jnew_last);
else
    Ar(1:jnew-1,jnew:jnew_last) = V(:,1:jnew-1)'*A*V(:,jnew:jnew_last);        
    Ar(jnew:jnew_last,1:jnew_last) = V(:,jnew:jnew_last)'*A*V(:,1:jnew_last);   
    Er(1:jnew-1,jnew:jnew_last) = V(:,1:jnew-1)'*E*V(:,jnew:jnew_last);
    Er(jnew:jnew_last,1:jnew_last) = V(:,jnew:jnew_last)'*E*V(:,1:jnew_last);
    Br(jnew:jnew_last,:) = V(:,jnew:jnew_last)'*B;
end
if nargout == 4,     varargout{1} = Cr;      end
end


function [ Xnorm ] = NormFrobEfficient( Y,varargin )
% NormFrobEfficient: an efficient way of calculating the Frobenius Norm of
% a symmetric matrix X

% Detailed explanation
% Input: factors G, D of the low rank factorization of X=G*D*G'
% Output: Frobenius Norm x (scalar)
%
% Remark:
% If the Frobenius Norm of e.g. X=B*B' is needed then set G=B and D=I
% because X=G*D*G' = B*I*B' = B*B'.
% For detailed explanation on low rank factorization see paper "Efficient
% low rank solution of generalized Lyapunov equations"

% read in input
if nargin == 2
    X = varargin{1};
    if istril(X) || istriu(X)
        X_L = X;
        X_D = speye(size(X_L,1));
    end
else
    X_L = varargin{2};
    X_D = varargin{3};
end

% if no matrix decomposition is available
if ~exist('X_L','var')
    try    
        X_L = chol(X);
        X_D = speye(size(X));
    catch
        try 
            [X_L,X_D] = ldl(X);
        catch
            % if all other matrix decompositions failed use onlx the
            % symmetric part of X
            [X_L,X_D] = ldl((X+X')/2);
        end
    end
end
% QR-decomposition of X_L
[~,X_R] = qr(full(X_L),0);

% determine Frobenius norm: || X ||_F
if Y == 1
    Xnorm = norm(X_R*X_D*X_R');
else
    Xnorm = norm(X_R*X_D*X_R','fro');
end
end




