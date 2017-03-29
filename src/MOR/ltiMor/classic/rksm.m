function [S,R,output_data] = rksm(varargin)
% My_rksm - Solve Laypunov equations with a cummulative rational Krylov
% subspace method
%
%              APE' + EPA' + BB' = 0 (I)
%              AQE' + EQA' + C'C = 0 (II)
%
% Synthax MY_RKSM
%       [S,output_data]           = MY_RKSM(A,B,s_inp) 
%       [S,output_data]           = MY_RKSM(A,B,s_inp,Rt) 
%       [S,R,output_data]         = MY_RKSM(A,B,C,s_inp) 
%       [S,R,output_data]         = MY_RKSM(A,B,C,s_inp,s_out)
%       [S,R,output_data]         = MY_RKSM(A,B,C,s_inp,Rt,Lt) 
%       [S,R,output_data]         = MY_RKSM(A,B,C,s_inp,s_out,Rt,Lt) 
%       [S,output_data]           = MY_RKSM(A,B,E,s_inp) 
%       [S,output_data]           = MY_RKSM(A,B,E,s_inp,Rt) 
%       [S,R,output_data]         = MY_RKSM(A,B,C,E,s_inp)
%       [S,R,output_data]         = MY_RKSM(A,B,C,E,s_inp,s_out)
%       [S,R,output_data]         = MY_RKSM(A,B,C,E,s_inp,Rt,Lt)
%       [S,R,output_data]         = MY_RKSM(A,B,C,E,s_inp,s_out,Rt,Lt) 
%       [S,R,output_data]         = MY_RKSM(A,B,C,s_inp,...,Opts_rksm) 
%       [S,R,output_data]         = MY_RKSM(sys) 
%       [S,R,output_data]         = MY_RKSM(sys,s_inp) 
%       [S,R,output_data]         = MY_RKSM(sys,s_inp,s_out)
%       [S,R,output_data]         = MY_RKSM(sys,s_inp,Rt)
%       [S,R,output_data]         = MY_RKSM(sys,s_inp,Rt,Lt) 
%       [S,R,output_data]         = MY_RKSM(sys,s_inp,s_out,Rt,Lt)
%       [S,R,output_data]         = MY_RKSM(sys,s_inp,...,Opts)
%
% Info: Funktionen f?r neue shifts muessen ab Zeile 470 unter der 'mess'-Option eingebunden werden 
%
%
% possible Options:
% - Opts.rksm_method: choose the method to be executed
%                     - standard: "real" RKSM-method using special shifts to expand the subspace
%                     - extended: extended Krylov subspaces are used
%                       ['rational' / 'extended']
%
% - Opts.shifts:      choose which shifts and how they should be used
%                     - cyclic: use the same shifts for the whole iteration
%                     - adaptive: a new shift is computed online for every
%                       single iteration (only for SISO-systems)
%                     - mess: use a function to generate shifts out of the
%                       mess-toolbox
%                       ['cyclic' / 'adaptive' / 'mess']
%
% - Opts.qrksm:       if no shifts are known, specify how many irka-shifts
%                     should be used
%                     default value:10
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
%
%% still to come/immer noch zu erledigen
% - adaptive SISO testen --> ganz nach hinten geschoben
%% Create Def-struct containing default values and make global settings
% Note: there may be no point named Def.reuseLU, otherwise one gets a conflict with solveLse/lyapchol/bilyapchol

Def.rksm_method          = 'rational';         % ['rational' / 'extended']
Def.shifts               = 'cyclic';           % ['cyclic' / 'adaptive' / 'mess']
Def.qrksm                =  10;                % default value for number of irka-shifts
Def.reduction            = 'onesided';         % ['onesided' / 'twosided']
Def.residual             = 'residual_lyap';    % ['residual_lyap' / 'norm_chol']
Def.rksmnorm             = 'H2';               % ['H2' / 'fro']
Def.lowrank              =  1;                 % [0 / 1]
Def.orth                 = '2mgs';             % for Gramschmidt
Def.rctol                =  1e-12;             % default tolerance
Def.maxiter_rksm         =  200;               % default number of iterations
Def.info_rksm            = 0;                  % show programme status information, [0 / 1]

% global variables
global hermite_gram_sch hermite withoutC Rt Lt s_out

% struct containing output information
output_data = struct('V_basis',[], 'W_basis',[], 'norm_val',[],'shifts',[],...
                     'Ar',[], 'Br',[], 'Er',[], 'Cr',[], 'rhsb',[], 'res0',[], 'lastnorm',[]);


%% parsing of inputs
if isa(varargin{end},'struct')  % wenn opts gesetzt ist, dann wird bei varargin opts nicht mitgezaehlt
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

% input of sys-object
if isa(varargin{1},'ss') || isa(varargin{1},'sss') || isa(varargin{1},'ssRed')
    
    sys = varargin{1};
    % set A, B matrix
    A = sys.A;      n = size(A,2);     
    B = sys.B;      m = size(B,2);    
    
    % check if there is an E matrix in the sys 
    try
        E = sys.E;
    catch
        E = speye(size(A,1));
    end
    
    % check if there is a C matrix in the sys and wether a one- or twosided reduction is performed 
    try
        C = sys.C;          p = size(C,1);
    catch
        disp('no C matrix available, onesided Krylov subspaces are used!');
    end
    if strcmp(Opts.reduction,'onesided') && exist('C','var')
        clear C p
    end
    
    switch length(varargin)
        
        case 2
            s_inp = varargin{2};
            
        case 3
            s_inp = varargin{2};
            if exist('C','var')
                s_out = varargin{3};
                if size(s_inp,1) > 2 || size(s_inp,2) ~= size(s_out,2) || size(s_out,1) > 2
                    error('shift vector s_inp (input space) and s_out (output space) must have the same dimensions');
                end
            else
                Rt = varargin{3};
                if size(Rt,1) ~= m || size(Rt,2) ~= size(s_inp,2)
                    error('shifts vector and tangential directions vector for input space must have the same dimensions');
                end
            end
            
        case 4
            s_inp = varargin{2};
            Rt = varargin{3};
            Lt = varargin{4};
            if (size(Rt,1) ~= m || size(Rt,2) ~= size(s_inp,2)) && (size(Lt,1) ~= p || size(Rt,2) ~= size(Lt,2))
                error('shifts vector and tangential directions vector for input space must have the same dimensions');
            end
            
        case 5
            s_inp = varargin{2};
            s_out = varargin{3};
            if size(s_inp,1) ~= size(s_out,1) || size(s_inp,2) ~= size(s_out,2)
                error('shift vector s0 (input space) and s1 (output space) must have the same dimensions');
            end
            Rt = varargin{4};
            Lt = varargin{5};
            if size(Rt,1) ~= size(Lt,1) || size(Rt,2) ~= size(Lt,2)
                error('tangential directions t0 (input space) and t1 (output space) must have the same dimensions');
            end
            if size(Rt,1) ~= m || size(Lt,1) ~= p || size(Rt,2) ~= size(s_inp,2) || size(Lt,2) ~= size(s_out,2)
                error('shifts vector and tangential directions vector for input/output space must have the same dimensions');
            end
    end % end of switch
    tol = Opts.rctol;
    maxiter = Opts.maxiter_rksm;

    if size(s_inp,2) == 1 && exist('s_inp','var')
        error('At least two starting shifts are necessary');
    end
              
% input of matrices    
elseif length(varargin) > 1
    
    % set A, B matrix
    A = varargin{1};         n = size(A,1);
    B = varargin{2};         m = size(B,2);
    
    % set C, E, shifts, tangential directions, tolerance and maxiter
    switch length(varargin)
        case 3
            s_inp = varargin{3};
            
        case 4
            if size(varargin{3},2) == n && size(varargin{3},1) ~= n
                C = varargin{3};
                p = size(C,1);
                s_inp = varargin{4};
            elseif size(varargin{3}) == size(A)
                E = varargin{3};
                s_inp = varargin{4};
            else
                s_inp = varargin{3};
                Rt = varargin{4};
                if size(Rt,1) ~= m || size(Rt,2) ~= size(s_inp,2)
                    error('shifts vector and tangential directions vector for input space must have the same dimensions');
                end
            end
            
        case 5
            if (size(varargin{1},1) == size(varargin{3},1) && size(varargin{1},2) == size(varargin{3},2)) ||...
               (size(varargin{1},1) == size(varargin{4},1) && size(varargin{1},2) == size(varargin{4},2))
                if size(varargin{1},1) == size(varargin{3},1)
                    E = varargin{3};
                    s_inp = varargin{4};
                    Rt = varargin{5};
                    if size(Rt,1) ~= m || size(Rt,2) ~= size(s_inp,2)
                        error('shifts vector and tangential directions vector for input space must have the same dimensions');
                    end
                else
                    C = varargin{3};
                    p = size(C,1);
                    E = varargin{4};
                    s_inp = varargin{5};
                end
            else
                C = varargin{3};
                p = size(C,1);
                s_inp = varargin{4};
                s_out = varargin{5};
                if size(s_inp,1) > 2 || size(s_inp,2) ~= size(s_out,2) || size(s_out,1) > 2
                    error('shift vector s0 (input space) and s1 (output space) must have the same dimensions');
                end
            end  
            
        case 6
            C = varargin{3};
            p = size(C,1);
            if size(varargin{1},1) == size(varargin{4},1) && size(varargin{1},2) == size(varargin{4},2)
                E = varargin{4};
                s_inp = varargin{5};
                s_out = varargin{6};
                if size(s_inp,1) > 2 || size(s_inp,2) ~= size(s_out,2) || size(s_out,1) > 2
                    error('shift vector s0 (input space) and s1 (output space) must have the same dimensions');
                end
            else
                s_inp = varargin{4};
                Rt = varargin{5};
                Lt = varargin{6};
                if (size(Rt,1) ~= m || size(Rt,2) ~= size(s_inp,2)) && (size(Lt,1) ~= p || size(Rt,2) ~= size(Lt,2))
                    error('shifts vector and tangential directions vector for input space must have the same dimensions');
                end
            end
            
        case 7
            C = varargin{3};
            p = size(C,1);
            if size(varargin{1},1) == size(varargin{4},1) && size(varargin{1},2) == size(varargin{4},2)
                E = varargin{4};
                s_inp = varargin{5};
                Rt = varargin{6};
                Lt = varargin{7};
                if (size(Rt,1) ~= m || size(Rt,2) ~= size(s_inp,2)) || (size(Lt,1) ~= p || size(Rt,2) ~= size(Lt,2))
                    error('shifts vector and tangential directions vector for input space must have the same dimensions');
                end
            else
                s_inp = varargin{4};
                s_out = varargin{5};
                if size(s_inp,1) ~= size(s_out,1) || size(s_inp,2) ~= size(s_out,2)
                    error('shift vector s0 (input space) and s1 (output spce) must have the same dimensions');
                end
                Rt = varargin{6};
                Lt = varargin{7};
                if size(Rt,1) ~= m || size(Rt,2) ~= size(s_inp,2) || size(Lt,1) ~= p
                    error('tangential directions t0 (input space) and t1 (output spce) must have the same dimensions');
                end
                if size(Rt,1) ~= m || size(Rt,2) ~= size(s_inp,2) || size(Lt,1) ~= p || size(Lt,2) ~= size(s_out,2)
                    error('shifts vector and tangential directions vector for input/output space must have the same dimensions');
                end
            end
            
        case 8
            C = varargin{3};
            p = size(C,1);
            E = varargin{4};
            s_inp = varargin{5};
            s_out = varargin{6};
            if size(s_inp,1) ~= size(s_out,1) || size(s_inp,2) ~= size(s_out,2)
                error('shift vector s0 (input space) and s1 (output spce) must have the same dimensions');
            end
            Rt = varargin{7};
            Lt = varargin{8};
            if size(Rt,1) ~= size(Lt,1) || size(Rt,2) ~= size(Lt,2)
                error('tangential directions t0 (input space) and t1 (output space) must have the same dimensions');
            end
            if size(Rt,1) ~= m || size(Lt,1) ~= p || size(Rt,2) ~= size(s_inp,2) || size(Lt,2) ~= size(s_out,2)
                error('shifts vector and tangential directions vector for input/output space must have the same dimensions');
            end
      
    end
    tol = Opts.rctol;
    maxiter = Opts.maxiter_rksm;
    
    
    if length(varargin) > 8 || length(varargin) < 3
        error('Wrong input');
    elseif (size(s_inp,1) > 2 && size(s_inp,2) > 2) ||  (size(s_inp,2) > 2 && size(s_inp,1) > 2)
        error('Wrong input, s0 must be either a vector or a 2 cross 1 matrix or a 1 cross 2 matrix ');
    elseif size(s_inp,1) == 1 && exist('s_inp','var')
        error('At least two starting shifts are necessary');
    end
end
    
%% Check and prepare inputs

% check if E matrix is available
if ~exist('E','var') || isempty(E)
    E = speye(size(A));
    withoutE = true;
else
    withoutE = false;
end

% check C matrix (twosided case), hermite case 
if ~exist('C','var') || isempty(C) || strcmp(Opts.reduction,'onesided')
    withoutC = true;
    p = 1;
else
    withoutC = false;
end

% check if hermite
if ~exist('s_out','var') || isempty(s_out)
    hermite = true;
else 
    hermite = false;
end

% check if there are shifts, if not use irka-shifts
if ~exist('s_inp','var') || isempty(s_inp)
    s_inp = zeros(1,Opts.qrksm);
    [~, ~, ~, s_inp] = irka(sys, s_inp);
end

%% RKSM Methods

% choose method
switch Opts.rksm_method
    
    case 'rational'
        % create shift matrix
        if withoutC == 1
            [s_ma,tangential] = make_shiftmatrix(s_inp,m);
        else
            [s_ma,tangential] = make_shiftmatrix(s_inp,m,p,A);
        end
        
        % sort shifts in an increasing way
        s_ma = mysort(s_ma);
        output_data.shifts = s_ma;
        clear s_inp s_out Rt Lt 
 
        % build the first basis with arnoldi, containing two directions or three if second shift is complex
        if withoutC == 1 && hermite == 1
            % build first input subspace (onesided), set shift variable jCol_inp, for SISO
            if (s_ma(1,1) == conj(s_ma(1,2)))
                [V] = arnoldi(E,A,B,s_ma(1,1:2));
                % initialize index j for current shift and tangential direction
                j = 2;
            else
                [V] = arnoldi(E,A,B,s_ma(1,1));
                % initialize index j for current shift and tangential direction
                j = 1;
            end    
        % build first input and output subspace (twosided, hermite)
        elseif withoutC == 0 && hermite == 1
            if (s_ma(1,1) == conj(s_ma(1,2)))
                [V,~,~,W] = arnoldi(E,A,B,C,s_ma(1,1:2));
                j = 2;
            else
                [V,~,~,W] = arnoldi(E,A,B,C,s_ma(1,1));
                j = 1;
            end           
        % build first input and output subspace (twosided, not hermite) 
        elseif withoutC == 0 && hermite == 0
            if (s_ma(1,1) == conj(s_ma(1,2)))
                [V] = arnoldi(E,A,B,s_ma(1,1:2));
                [~,~,~,W] = arnoldi(E,A,B,C,s_ma(2+m+p,1:2));
                j = 2;
            else
                [V] = arnoldi(E,A,B,s_ma(1,1:2));
                [~,~,~,W] = arnoldi(E,A,B,C,s_ma(2+m+p,1:2));
                j = 1;
            end
        end
        
        % calculating the first reduced system, ensure, that V is real
        V = real(V);
        if withoutC == 1
            Ar = V'*A*V;
            Br = V'*B;
            Er = V'*E*V;
        else
            W = real(W);
            Ar = W'*A*V;
            Br = W'*B;
            Er = W'*E*V;
            Cr = C*V;
        end
                
        % initialize counter for higher moments
        if j == 1
            counter_inp = 0;
        elseif (s_ma(1,j) == s_ma(1,j-1))
            counter_inp = 1;
        elseif (s_ma(1,j) ~= s_ma(1,j-1))
            counter_inp = 0;
        end
        if hermite == 0 
            if j == 1
                counter_out = 0;
            elseif (s_ma(2+m+p,j) == s_ma(2+m+p,j-1))
                counter_out = 1;
            elseif (s_ma(1,j) ~= s_ma(1,j-1))
                counter_out = 0;
            end
        end
        
        % start iteration, in the first iteration one tries to find a solution
        % in the subspace calculated above with the arnoldi algorithm
        % Initializations
        if  size(V,2) == 2*m, k = 2; else k = 1; end
        %end
        wait1 = false;   index = 0;
        
        for ii = k:1:maxiter
            if ii > k
                if wait1 == 0
                    % increase index
                    j = j+1;

                    % start again with the first shift in s_inp/s_out and reset
                    % counter and complex flag or use new shifts
                    if mod(j,length(s_ma(1,:))+1) == 0
                        % get new shifts
                        if strcmp(Opts.shifts,'cyclic')
                            j = 1;
                        elseif strcmp(Opts.shifts,'adaptive')
                            [new_s] = adaptive_shift(Ar,output_data.shifts(1,size(output_data.shifts,2)-1),output_data.shifts(1,size(output_data.shifts,2)));
                            
                            % create new shift matrix and save shift
                            [s_ma] = make_shiftmatrix(new_s,m,p,A);
                            output_data.shifts = [output_data.shifts s_ma];
                            j = 1;
                        elseif strcmp(Opts.shifts,'mess')
                            % hier irgendeine andere Funktion um shifts zu
                            % berechnen

                            % hier einstellen, ab welcher iteration neue
                            % shifts gerufen werden sollen
                            if ii < 25 
                                j=1;
                                kk = 1;
                            else
                                if ~exist('kk','var')
                                    kk = 1;
                                end
                                % start-Werte f?r irka, wenn mess dann
                                % auskommentieren
                                shifts = zeros(1,4+kk);
                                Rt = ones(m,4+kk);        
                                Lt = ones(p,4+kk);
            
                                % sys-Objekt erstellen, es ist m?glich auch
                                % mit dem rteduzierten Modell zu arbeiten
                                C = ones(1,n);
                                D = zeros(size(C,1),m);
                                sys = sss(A,B,C,D,E);

                                %Cr = ones(1,size(Ar,1));
                                %Dr = zeros(size(Cr,1),size(Br,2));
                                %Ar = real(Ar); Br = real(Br); Cr = real(Cr);  Er = real(Er);
                                %sys = sss(Ar,Br,Cr,Dr,Er);
                                if ~exist('C','var')
                                    C = Opts.Cma(1,:);
                                end
%                                 [~, ~, ~, s0, Rt, Lt] = irka(sys,shifts,Rt,Lt);
% 
%                                 [s_ma] = make_shiftmatrix(s0,m,p,A);
%                                 kk = kk + 1;

                                    % mess-shifts zum einkommentieren


                                  eqn=struct('A_',sys.A,'E_',sys.E,'B',sys.B,'C',sys.C,'prm',speye(size(sys.A)),'type','N','haveE',sys.isDescriptor);
                                  % opts struct: MESS options
                                  messOpts.adi=struct('shifts',struct('l0',20,'kp',50,'km',25,'b0',ones(sys.n,1),...
                                  'info',1),'maxiter',Opts.maxiter,'restol',0,'rctol',Opts.rctol,...
                                  'info',1,'norm','fro');
                                  Opts.lse         = 'gauss'; 
                                  lseType='solveLse';
                                  oper = operatormanager(lseType);
                                  messOpts.solveLse.lse=Opts.lse;
                                  messOpts.solveLse.krylov=0;
                                  % get adi shifts
                                  [s0,~,~,~,~,~,~,eqn]=mess_para(eqn,messOpts,oper); s0=s0';
                                  Rt = ones(m,size(s0,2));
                                  Lt = ones(p,size(s0,2));
                                  [~, ~, ~, s0, Rt, Lt] = irka(sys,s0,Rt,Lt);
                                  [s_ma] = make_shiftmatrix(s0,m,p,A);

                                  % hier das j+1 muss immer gesetzt werden
                                  % fuer einen guten Programmablauf
                                  j = 1;
                            end
                        end
                        counter_inp = 0;       
                        high_inp = false;
                        
                        if hermite == 0
                            counter_out = 0;       high_out = false; 
                        end
                    end

                    % save last shift,read in next shift, check order for input and output shifts
                    % save last shift
                    if ii <= length(s_ma(1,:))
                        jCol_inp = s_ma(1,j);     jCol_inp_old = s_ma(1,j-1); 
                    elseif ii > length(s_ma(1,:)) 
                        jCol_inp_old = jCol_inp;  jCol_inp = s_ma(1,j); 
                    elseif ~exist('jCol_inp','var')
                        jCol_inp = s_ma(1,j);  jCol_inp_old = s_ma(1,j-1);  
                    end
                        jCol_Rt = s_ma(2,j);       jCol_Lt = s_ma(3,j);  

                    % check order of current shift
                    if isreal(jCol_inp)
                        if jCol_inp == jCol_inp_old
                            high_inp = true;
                            counter_inp = counter_inp+1;
                        else
                            counter_inp = 0;
                            high_inp = 0;
                        end
                    end
                    if hermite == 0
                        if ii < length(s_ma(1,:))
                            jCol_out = s_ma(4,j);     jCol_out_old = s_ma(4,j-1); 
                        else
                        jCol_out_old = jCol_out;       jCol_out = s_ma(4,j);
                        end
                        if isreal(jCol_out)
                            if jCol_out == jCol_out_old
                                high_out = true;
                                counter_out = counter_out+1;
                            else   counter_out = 0;       high_out = false;
                            end
                        end
                    end

                    % check if shift is complex, check order of shift
                    if ~isreal(jCol_inp) && j>2
                        % check order of shift
                        if (jCol_inp == s_ma(1,j-2)) && imag(jCol_inp) < 0
                            high_inp = true;
                            counter_inp = counter_inp+1;
                        elseif jCol_inp ~= conj(s_ma(1,j-1))
                            counter_inp = 0;    high_inp = false;
                        end
                    elseif ~isreal(jCol_inp) && j<=2
                        high_inp = 0;
                    end

                    if hermite == 0 % hier hermite komplex muss ich machen
                        if ~isreal(jCol_out) && j>2
                            % check order of shift
                            if (jCol_out == s_ma(4,j-2)) && imag(jCol_out) < 0
                                high_out = true;
                                counter_out = counter_out+1;
                            elseif jCol_out ~= conj(s_ma(4,j-1))
                                counter_out = 0;    high_out = false;
                            elseif ~isreal(jCol_out) && j<=2
                               high_out = 0;
                            end
                        end
                    end

                    % Build cummulative basis for SISO and MIMO
                    if tangential == 1
                        if size(V,2) == m
                            jbasis    = 1;                   % index for the first column of the directions calculated before
                            jnew      = jbasis+m;            % index for the first new column
                        elseif size(V,2) == 2*m
                            jbasis    = m+1;
                            jnew      = jbasis+m;
                        else
                            jbasis    = size(V,2);  
                            jnew      = jbasis+1;       
                        end
                    else
                        if size(V,2) == m || size(V,2) > 2*m
                            jbasis    = (ii-2)*m+1;       
                            jnew      = jbasis+m;        
                        else
                            jbasis    = m+1;              
                            jnew      = jbasis+m;         
                        end
                    end

                    % calculating the new directions input space
                    % Opts.getLU = 1;       comment in if Opts.reuseLU does not work
                    rhsB = V(:,jbasis:jnew-1);     % new rhs
                    %output_data.rhsb(:,ii) = rhsB;
                    Anew_V = (A-jCol_inp*E);       % A matrix for solveLse function for V direction
                    if withoutC == 0 && hermite == 0
                        Anew_W = (A-jCol_out*E)';
                    elseif withoutC == 0 && hermite == 1
                        Anew_W = (A-jCol_inp*E)'; 
                    end
                    
                    if tangential == 1,    rhsB = rhsB*jCol_Rt;    end
                    
                    if tangential ==1 && withoutC == 0
                        rhsC = W(:,jbasis:jnew-1);
                        rhsC = rhsC*jCol_Lt;
                    elseif tangential ==0 && withoutC == 0
                        rhsC = W(:,jbasis:jnew-1);
                    end
                    
                    % calculation of new V directions
                    if high_inp == 0
                        [vnew] = solveLse(Anew_V,rhsB);
                    else
                        Opts.reuseLU = 1;
                        [vnew] = solveLse(Anew_V,E*rhsB,Opts);
                    end
                    
                    % calculation of new W directions
                    if withoutC == 0 && hermite == 1
                        if high_inp == 0
                            [wnew] = solveLse(Anew_W,rhsC);
                        else
                            Opts.reuseLU = 1;
                            [wnew] = solveLse(Anew_W,E'*rhsC,Opts);
                        end
                    end
                    
                    if withoutC == 0 && hermite == 0
                        if high_out == 0
                            [wnew] = solveLse(Anew_W,rhsC);
                        else
                            Opts.reuseLU = 1;
                            [wnew] = solveLse(Anew_W,E'*rhsC,Opts);
                        end
                    end
                    
                    if tangential == 1
                        vnew = vnew(:,1);
                        if withoutC == 0
                            wnew = wnew(:,1);
                        end
                    else
                        vnew = vnew(:,1:m);
                        if withoutC == 0
                            wnew = wnew(:,1:m);
                        end
                    end
               end % end of wait

               % build and orthogonalize new basis
               if ~isreal(vnew)
                   % if new direction vnew is complex, now after expanding
                   % the basis with the real part it follows the complex part
                   % note: until now no further shift is read out!
                   if index == ii-1
                       vnew_imag = imag(vnew);
                       V = [V vnew_imag];
                       
                       if withoutC == 0 && ~isreal(wnew)
                           wnew_imag = imag(wnew);
                           W = [W wnew_imag];
                       
                       % here a higher order is calculated with the current s_out because s_out is real and not complex
                       elseif withoutC == 0 && isreal(wnew)
                           % calculate additional direction wnew by determining a higher order of s_out_j
                           Opts.reuseLU = 1;
                           rhsC = W(:,jbasis+m:jnew-1+m);   % new rhs achte hier auf +m
                           if hermite == 1 && tangential == 1
                               rhsC = rhsC*jCol_Lt;
                           end
                           if hermite == 1
                               Anew_W = (A-jCol_inp*E)';
                           else
                               Anew_W = (A-jCol_out*E)';
                           end
                            
                           [wnew] = solveLse(Anew_W,E'*rhsC,Opts);
                           
                           % if the higher order of wnew is not real use the real part
                           if ~isreal(wnew)
                               wnew = wnew(:,m);   wnew = real(wnew);
                               W = [W wnew];
                           else
                               W = [W wnew];
                           end  
                       end
                       % reset variables
                       wait1 = false;
                       index = 0;
                   else
                       % if new direction vnew is complex, first expand basis with the real part of vnew
                       % hier gehe ich immer zuerst rein
                       vnew_real = real(vnew);
                       V = [V vnew_real];
                       
                       if withoutC == 0 && ~isreal(wnew)
                           wnew_real = real(wnew);
                           W = [W wnew_real];
                       elseif withoutC == 0 && isreal(wnew)
                           W = [W wnew];
                       end
                       
                       
                       % some settings for the code
                       wait1 = true;
                       index = ii;
                   end
               else
                   % if vnew is real
                   if wait1 == 0
                   V = [V vnew];
                   end
                   
                   % if wnew is complex
                   if withoutC == 0 && ~isreal(wnew)

                       if index == ii-1
                           % if new direction wnew is complex, now after expanding
                           % the basis with the real part it follows the complex part of wnew
                           % note: until now no further shift is read out!
                           wnew_imag = imag(wnew);
                           W = [W wnew_imag];
                           
                           % calculate additional direction vnew with the same shift but of higher order
                           Opts.reuseLU = 1;
                           rhsB = V(:,jbasis+m:jnew-1+m);    % new rhs for B
                           Anew_V = (A-jCol_inp*E);
                           if tangential == 1
                               rhsB = rhsB*jCol_Rt;
                           end
                           
                           [vnew] = solveLse(Anew_V,E*rhsB,Opts);

                           if ~isreal(vnew)
                               vnew = vnew(:,m);   vnew = real(vnew);
                               V = [V Vnew];
                           else
                               V = [V vnew];
                           end 
                           
                           % settings
                           wait1 = false;
                           index = 0;
                           
                       else
                           % if new direction wnew is complex, first expand basis W with the real part of wnew
                           % hier gehe ich immer zuerst rein
                           wnew_real = real(wnew);
                           W = [W wnew_real];
                       
                           % some settings for the code
                           wait1 = true;
                           index = ii;   
                       end
                       
                   elseif withoutC == 0 && isreal(wnew)
                      W = [W wnew];    
                   end
               end
              
               % set indices for calculation of reduced system
               if tangential == 1
                   jbasis    = size(V,2)-1;        % index for the first column of the directions calculated before   
                   jnew      = jbasis+1;        
                   jnew_last = jnew;             % index of the last new column     
               else
                   jbasis    = (ii-2)*m+1;     
                   jnew      = jbasis+m;
                   jnew_last = jnew+m-1;         % index of the last new column
                
               end
               
               % ensure, that V,W is real
               V = real(V);
               if withoutC == 0
                   W = real(W);
               end
                
               for jj=jnew:1:jnew_last
                   if withoutC == 1
                       V = gramSchmidt(jj,V,Opts);
                   else
                       hermite_gram_sch = 1;
                       [V,~,W] = gramSchmidt(jj,V,W,Opts);
                   end
               end

               % calculating the ruduced system in a cheap manner
               if withoutC == 1
                   [Ar,Er,Br] = reduction(A,B,E,Ar,Br,Er,jnew,jnew_last,V);
               else
                   [Ar,Er,Br,Cr] = reduction(A,B,E,Ar,Br,Er,jnew,jnew_last,V,W,C,Cr);
               end
           end % end of ii>k
           
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
           
           % choose computation of norm/stopping criteria
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
           if Rnorm < tol
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
        end  % end of for loop
        if ii == maxiter && Rnorm > tol
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
               F = E*V*(Er_inv_Br+(S*S')*Cr_hat');
 
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
            
            if Rnorm < tol
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

% clear clobal variables
clearvars hermite_gram_sch hermite withoutC Rt Lt s_out
end


% *************************************************************************
%% ***************************** AUXILIARY ********************************
% *************************************************************************


function x = make_vect(x,y)
% function which transforms a matrix with not more than 2 rows and 1 column or 2 columns and 1 row in a vector
% or just transposes a vector
% function also sorts the vector in an increasing way
% Input: - shift vector s_inp / s_out
%        - tangential directions Rt an Lt
% Output: sorted vectors, see example
%
% Example: x = [0 1 2 3 4; ---> x = [0 1 1 2 2 2 3 4 4]
%               1 2 3 1 2]


if size(x,1) == 2
    temp = zeros(1,sum(x(2,:)));
    for ii=1:1:size(x,2)
        k = sum(x(2,1:(ii-1)));
        k = (k+1):(k+x(2,ii));
        temp(k) = x(1,ii)*ones(1,x(2,ii));
    end
    x = temp;
else
    temp = zeros(y,sum(x(y+1,:)));
    for ii=1:1:size(x,2)
        k = sum(x(y+1,1:(ii-1)));
        k = (k+1):(k+x(y+1,ii));
        temp(1:y,k) = x(1:y,ii)*ones(1,x(y+1,ii));
    end
     x = temp;
end

end


function [s_ma,tangential] = make_shiftmatrix(s_inp,m,varargin)
% This function ensures that for every single shift the same number of
% moments is matched especially if there are diffrent shift vectors for the
% input and output Krylov space or if there are complex shifts.
%
% 1. Example: only real shifts: s_inp = [0 1 2 3 4;   s_out = [3 5 6 7 9;
%                                        2 1 2 2 3]            1 1 1 1 1]
%
% As you can see the input space V consists of more directions as the
% outputspace W which leads to an error. To avoid this the number of higher
% moments of the shift with the lower number of higher moments of one pair 
% is increased so that both numbers are equal, see the following:
%
%       First pair of shifts:   s_inp(1,1) = 0                   s_out(1,1) = 3
%                               --> number of moments: 2         --> number of moments = 1 
%                               --> unequal number of directions in V and W                            
%                               --> adding one higher moment to shift s_out = 3                             
% 
% Filling up the shift vectors in this way is possible because one just
% works with a Taylor approximation of higher order
%
% If there are complex shifts note that one gets for evry complex shift and
% for every higher order two directions. Just imagine the first shift in
% the above example is complex, then there are there are 4 directions for
% the first shift in s_inp and only one for the first shift in s_out
% --> 3 higher order moments have to be added in s_out(2,1).

global hermite withoutC Rt Lt s_out 

if length(varargin) == 1
    p = varargin{1};
elseif length(varargin) == 2
    p = varargin{1};
    A = varargin{2};
else
    p = 1;
end

 % transpose shift vectors and tangential directions
if size(s_inp,1) >= 2 && size(s_inp,2) == 1 
    s_inp = s_inp';
elseif exist('s_out','var') && ~isempty(s_out) && size(s_out,1) >= 2 && size(s_out,2) == 1
    s_out = s_out';
end

% build matrix containing all the information of shifts, tangential directions
% preallocate memory
if hermite == 1 && withoutC ==1
    s_info = zeros(3+m,size(s_inp,2));
elseif hermite == 1 && withoutC == 0
    s_info = zeros(2+m+p,size(s_inp,2));
else
    s_info = zeros(4+m+p,size(s_inp,2));
end

% read in shifts and order of shifts
if size(s_inp,1) == 2
    s_info(1,:) = s_inp(1,:);
    s_info(2,:) = s_inp(2,:);
else
    s_info(1,:) = s_inp;
    s_info(2,:) = ones(1,size(s_inp,2));
end

% read in tangential directions
if (~exist('Rt','var') || isempty(Rt)) && (~exist('Lt','var') || isempty(Lt))
    s_info(3:2+m,:) = ones(m,size(s_inp,2));
    s_info(3+m:2+p+m,:) = ones(p,size(s_inp,2));
    tangential = false;
elseif exist('Rt','var') && (~exist('Lt','var') || isempty(Lt))
    s_info(3:2+m,:) = Rt;
    s_info(3+m:2+p+m,:) = ones(p,size(s_inp,2));
    tangential = true;
elseif exist('Rt','var') && exist('Lt','var') && (~isempty(Rt) && ~isempty(Lt))
    s_info(3:2+m,:) = Rt;
    s_info(3+m:2+p+m,:) = Lt;
    tangential = true;
end

% read in s_out
if hermite == 0
    s_info(3+m+p,:) = s_out(1,:);
    if size(s_out,1) == 2
        s_info(4+p+m,:) = s_out(2,:);
    else
        s_info(4+p+m,:) = ones(1,size(s_inp,2));
    end
end

%neu
s_info = single(s_info);
%neu

% set variables
s_in = s_info(1,:);
if size(s_info,1) > (2+m+p)
    s_out = s_info(3+m+p,:);
end

% if there are complex shifts
if ~isreal(s_in(1,:))
    % look for complex pairs
    try
        s_in = cplxpair(s_in);
        % neu
        %s_in_new = cplxpair(s_in);
        % neu ende
    catch
        % add conjugate complex elements 
        for ii=1:1:size(s_in,2)
            if ~isreal(s_in(1,ii))
                counter = 1;
                while counter ~= 0
                    if s_in(1,ii) == conj(s_in(1,counter))
                        counter = 0;
                    else
                        counter = counter+1;
                    end
                    if counter > size(s_in,2)
                        counter = 0;
                        s_in(1,:) = [s_in(1,:) conj(s_in(1,ii))];
                    end
                end   
            end   
        end
        s_in(1,:) = cplxpair(s_in(1,:));
        clear counter
    end
    
    %if output == 1 && no_complex_out_s == 1
    if exist('s_out','var') && ~isreal(s_out) && ~isempty(s_out)
        try
            s_out = cplxpair(s_out);
            if size(s_in,2) ~= size(s_out,2)
                lastwarn('try to use a shift vectors containing complex pairs of shifts if there are complex shifts')
            end
        catch
            % add conjugate complex elements
            for ii=1:1:size(s_out,2)
                if ~isreal(s_out(1,ii))
                    counter = 1;
                    while counter ~= 0
                        if s_out(1,ii) == conj(s_out(1,counter))
                            counter = 0;
                        else
                            counter = counter+1;
                        end
                        if counter > size(s_out,2)
                            counter = 0;
                            s_out(1,:) = [s_out(1,:) conj(s_out(1,ii))];
                        end
                    end   
                end   
            end
            %s_out = cplxpair(s_out);
            if size(s_in,2) ~= size(s_out,2)
                lastwarn('try to add the following shifts to your vector:')
            end
            clear counter   
       end
    end
    
    % build expanded shift information matrix
    % preallocate memory
    if ~exist('s_out','var') || isempty(s_out)
        s_info_new = zeros(2+m+p,size(s_in,2));
    else
        s_info_new = zeros(4+m+p,size(s_in,2));
    end
    
    for ii=1:1:size(s_in,2)
        counter = 1;
        while counter ~= 0
            if s_in(1,ii) == s_info(1,counter)
                s_info_new(1,ii) = s_in(1,ii);
                if ~exist('s_out','var') || isempty(s_out)
                    s_info_new(2:2+p+m,ii) = s_info(2:2+p+m,counter);
                else
                    s_info_new(2:4+p+m,ii) = s_info(2:4+p+m,counter);
                end
                counter = 0;
            elseif s_in(1,ii) ~= s_info(1,counter) && counter == size(s_info,2)
                s_info_new(1,ii) = s_in(1,ii);
                s_info_new(2,ii) = s_info_new(2,ii-1);       % fill up order
                s_info_new(3:2+m,ii) = ones(m,1);            % fill up Rt 
                s_info_new(3+m:2+m+p,ii) = ones(m,1);        % fill up Lt
                if ~exist('s_out','var') || isempty(s_out)
                    %add_shifts = -eigs(A); % spaltenvektor
                    %[~, ~, ~, s0, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = irka(sys, add_shifts);
                    s_info_new(3+m+p:4+m+p,ii) = ones(2,1);  % fill up output
                end
                counter = 0;
            elseif s_in(1,ii) ~= s_info(1,counter) || counter <= size(s_info,2)
                counter = counter+1;
            end  
        end  
    end
    s_info = s_info_new;
    clear counter s_info_new
        
    for ii=1:1:size(s_in,2)
        if ~isreal(s_info(1,ii)) && imag(s_info(1,ii)) < 0 && s_info(2,ii) >= s_info(2,ii+1)
            s_info(2,ii+1) = s_info(2,ii);
        elseif ~isreal(s_info(1,ii)) && imag(s_info(1,ii)) < 0 && s_info(2,ii) < s_info(2,ii+1)
            s_info(2,ii) = s_info(2,ii+1);          
        end      
    end
    
    % compare oder of input and output shifts
    if exist('s_out','var') && ~isempty(s_out)
        for ii=1:1:size(s_info,2)
            if s_info(2,ii) >= s_info(4+m+p,ii)
                s_info(4+m+p,ii) = s_info(2,ii);
            else
                s_info(2,ii) = s_info(4+m+p,ii);
            end        
        end
    end
   
% only real shifts
else
    if exist('s_out','var') && ~isempty(s_out)
        for ii=1:1:size(s_info,2)
            if s_info(2,ii) >= s_info(4+m+p,ii)
                s_info(4+m+p,ii) = s_info(2,ii);
            else
                s_info(2,ii) = s_info(4+m+p,ii);
            end        
        end
    end
end

% create input shift matrix (contains shifts and order of each shift) and
% build shift vector
s_inp = s_info(1:2,1:size(s_info,2));
s_inp = make_vect(s_inp);
s_inp = cplxpair(s_inp);

% create right tangential directions matrix (contains directions and order of each shift) and
% build tangetial vector/matrix
Rt(1:m,1:size(s_info,2)) = s_info(3:2+m,1:size(s_info,2));
Rt(m+1,1:size(s_info,2)) = s_info(2,1:size(s_info,2));
Rt = make_vect(Rt,m);

% create left tangential directions matrix (contains directions and order of each shift) and
% build tangetial vector/matrix
Lt(1:p,1:size(s_info,2)) = s_info(3+m:2+m+p,1:size(s_info,2));
Lt(p+1,1:size(s_info,2)) = s_info(2,1:size(s_info,2));
Lt = make_vect(Lt,p);

% eventually create output shift vector
if hermite == 0
    s_out = s_info(3+m+p:4+m+p,1:size(s_info,2));
    s_out = make_vect(s_out);
end

% collect all shift, direction vectors in one matrix
s_ma(1,:) = s_inp;       
s_ma(2:m+1,:) = Rt;
s_ma(2+m:1+m+p,:) = Lt;
if hermite == 0
    s_ma(2+m+p,:) = s_out;
end
clear s_in s_info s_inp
clearvars -global Rt Lt s_out
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




function x = mysort(x)
n = size(x,2);
while n>1
    newn = 1;
    for ii=1:1:n-1
        if real(x(1,ii)) > real(x(1,ii+1)) 
            temp = x(:,ii);
            x(:,ii) = x(:,ii+1);
            x(:,ii+1) = temp;
            newn = ii+1;
        end
        if isnan(real(x(1,ii))) && ~isnan(real(x(1,ii+1)))
            temp = x(:,ii);
            x(:,ii) = x(:,ii+1);
            x(:,ii+1) = temp;
            newn = ii+1;
        end
    end
    n = newn;
end
end


function [snew] = adaptive_shift(Ar,s1,s2)

s_vect(1,1) = s1;       s_vect(1,2) = s2;

% compute eigenvalues of Ar
eigAr = eig(Ar);

% build spectral set of subspace (spec_set is a vector) / eHpoints ist spec
% set
if isreal(s1) && isreal(s2)
    spec_set = sort([s1 -s1 s2 -s2 (-eigAr)']);
elseif isreal(s1) && ~isreal(s2)
    spec_set = sort([s1 -s1 s2 conj(s2) (-eigAr)']);
elseif ~isreal(s1) && isreal(s2)
    spec_set = sort([s1 conj(s1) s2 -s2 (-eigAr)']);
else
    spec_set = sort([s1 conj(s1) s2 conj(s2) (-eigAr)']);
end

% build convex hull or real poles
if any(abs(imag(spec_set))) == 1
    hull_conv = convhull(real(spec_set),imag(spec_set));
    spec_set = hull_conv;
end

% build residual; product das mit den z siehe formel 2.2
for ii=1:1:length(spec_set)-1
    v = linspace(spec_set(ii),spec_set(ii+1),500);
    for jj=1:1:length(v)
        res(jj) = prod((v(jj)-s1)./(v(jj)-eigAr));
    end
    [~,temp] = max(abs(res));
    snew(ii) = v(temp);
end
for ii=1:1:length(snew)
    res(ii) = prod((snew(ii)-s1)./(snew(ii)-eigAr));
end
[~,temp] = max(abs(res));
snew = snew(temp);  
s_vect(1,ii) = snew; 


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




