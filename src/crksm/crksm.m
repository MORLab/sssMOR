function [sysr,data] = crksm(varargin)
% CRKSM - Cumulative Rational Krylov Subspace Method for cumulative reduction and/or approximately solving Lyapunov equations 
%
% Syntax:
%       [sysr,data]         = CRKSM(sys, s0_inp)
%       [sysr,data]         = CRKSM(sys, s0_inp, Rt)
%       [sysr,data]         = CRKSM(sys, [], s0_out)
%       [sysr,data]         = CRKSM(sys, [], s0_out, [], Lt)
%       [sysr,data]         = CRKSM(sys, s0_inp, s0_out)
%       [sysr,data]         = CRKSM(sys, s0_inp, s0_out, Rt, Lt)
%       [sysr,data]         = CRKSM(sys,...,Opts_rksm)
%       [sysr,data]         = CRKSM(A,B,[],[],s0_inp) 
%       [sysr,data]         = CRKSM(A,B,[],[],s0_inp,Rt) 
%       [sysr,data]         = CRKSM(A,B,C, [],s0_inp)
%       [sysr,data]         = CRKSM(A,B,C, [],s0_inp,s0_out) 
%       [sysr,data]         = CRKSM(A,B,C, [],[],s0_out) 
%       [sysr,data]         = CRKSM(A,[],C,[],[],s0_out,[],Lt) 
%       [sysr,data]         = CRKSM(A,B,C, [],s_inp,s_out,Rt,Lt) 
%       [sysr,data]         = CRKSM(A,B,[] ,E,s0_inp) 
%       [sysr,data]         = CRKSM(A,B,[] ,E,s0_inp,Rt) 
%       [sysr,data]         = CRKSM(A,B,C,E,s0_inp) 
%       [sysr,data]         = CRKSM(A,B,C,E,s0_inp,s0_out) 
%       [sysr,data]         = CRKSM(A,B,C,E,[],s0_out) 
%       [sysr,data]         = CRKSM(A,B,C,E,s0_inp,s0_out,Rt,Lt) 
%       [sysr,data]         = CRKSM(A,B,...,s_inp,...,Opts_rksm) 
%
% Description:
%
%              $$    A X E^T + E X A^T + B B^T = 0 \quad   (1)    $$
%
%              $$    A^T Y E + E^T Y A + C^T C = 0 \quad   (2)    $$
%
% Input Arguments:
%		*Required Input Arguments:*
%       -sys:                   sss-object containing LTI system 
%       -A/B/C/E:               system matrices
%       -s0_inp:                initial expansion points for input Krylov subspace (must have at least two entries)
%
%		*Optional Input Arguments:*
%       -s0_out:                expansion points for output Krylov subspace
%       -Rt/Lt:                 right/left tangential directions (MIMO case)
%       -Opts:                  a structure containing following options
%           -.purpose:          purpose of using CRKSM [{'lyapunov'} / 'MOR']
%                -'lyapunov':       use CRKSM for approximately solving Lyapunov equations
%                -'MOR':            use CRKSM for cumulative and adaptive model order reduction
%           -.shifts:           choose which shifts and how they should be used
%                -'cyclic':         use the same shifts for the whole iteration
%                -'adaptive':       a new shift is computed online for every single iteration
%                -'mess':           use a function to generate shifts out of the mess-toolbox
%
%           -.residual:         specify determination criteria [{'residual_lyap'} / 'norm_chol']
%                -'residual_lyap':  compute residual of Lyapunov equation
%                -'norm_chol':      compare  the norm of the two last Cholesky factors
%
%           -.rksmnorm:         specify norm [{'H2'} / 'fro']
%                -'H2':             use 2-norm (Euclidian Norm)
%                -'fro':            use Frobenius Norm
%
%           -.lowrank:          compute the low-rank factor of the final solution 
%                               [{0} / 1]
%
%           -.orth:             specify orhtogonalization method in Gram-Schmidt
%                               [{'2mgs'} / 'dgks' / 'mgs']
%
%           -.tol:              specify tolerance
%
%           -.maxiter:          specify maximal number of iterations
%
% Output Arguments:
%       -sysr:                  reduced system
%       -S:                     Cholesky factor X=S*S' of Lyapunov equation A*X*E'+E*X*A'+B*B'=0
%       -R:                     Cholesky factor Y=R*R' of lyapunov equation A'*Y*E+E'*Y*A+C'*C=0
%       -data:                  struct containing the following data
%           -V/W:               V_basis, W_basis
%           -norm_val:          norm values of the iterations
%           -TODO:              additional shifts, 
%           -Ar/Br/Cr:          reduced system matrices of last iteration
%           -rhs:               last rhs
%           -res0:              reference residual norm
%           -TODO:              last norm value
%
% Examples:
%       This code computes ....
%
%> sys = sss('fom');
%> Opts.purpose = 'MOR';
%> sysr = crksm(sys, -eigs(sys,8).');
%> bode(sys,'-',sysr,'--r');
%
% See Also: 
%       rk, arnoldi, solveLse, sss/lyapchol, mess_lradi, getShifts
%
% References:
%       * *[1] Druskin, Simoncini (2011)*, Adaptive Rational Krylov Subspaces
%       for large-scale dynamical systems
%       * *[2] Druskin, Simoncini, Zaslavsky (2014)*, Adaptive Tangential
%       Interpolation in Rational Krylov Subspaces for MIMO Dynamical Systems
%       * *[3] Kürschner (2016)*, Efficient Low-Rank Solution of Large-Scale Matrix Equations
%       * *[4] Wolf (2014)*, H2 Pseudo-Optimal Moder Order Reduction
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Paul Heidenreich, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Oct 2017
% Copyright (c) 2016-2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Create Def-struct containing default values
% Note: there may be no point named Def.reuseLU, otherwise one gets a conflict with solveLse/lyapchol/bilyapchol

% general default option settings for MOR and Lyapunov
Def.purpose              = 'lyapunov';         % [{'lyapunov'} / 'MOR']
Def.shifts               = 'fixedCyclic';      % [{'fixedCyclic'} / 'dynamical']
Def.shiftTol             =  0.1;               % default value for new shifts
Def.strategy             = 'ADI';              % [{'ADI'} / 'adaptive' / '' / '' / '']
Def.orth                 = '2mgs';             % for Gram-Schmidt
Def.lse                  = 'sparse';           % for solving a LTI-system
Def.maxiter_rksm         =  200;               % default number of iterations
Def.shiftsTol            =  0.1;               % default value for new shifts

% default option settings for MOR
Def.real                 = true;               % [{true} / false], true means to keep the subspace real

% default option settings for Lyapunov
Def.residual             = 'residual_lyap';    % [{'residual_lyap'} / 'norm_chol']
Def.rctol                =  1e-12;             % default tolerance
Def.rksmnorm             = 'H2';               % [{'H2'} / 'fro']
Def.lowrank              =  0;                 % [{0} / 1] 

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

% input of sys-object
if isa(varargin{1},'ss') || isa(varargin{1},'sss') || isa(varargin{1},'ssRed')

    % read in sys object
    sys = varargin{1};

    % use hess if sys is ssRed object
    if isa(sys,'ssRed')
        Opts.lse='hess'; 
        if isempty(sys.E), sys.E = eye(size(sys.A)); end %ssRed robust compatibility
    end

    % check usage and inputs
    if length(varargin) == 2                            % usage: CRKSM(sys, s0_inp)
        s0_inp = varargin{2};                s0_out  = [];        
        Rt     = [];                         Lt      = [];
        input  = 1;                          pointer = @blockV;
    elseif length(varargin) == 3       
        if size(varargin{3},1) == size(sys.B,2) && sys.isSiso == 0 &&...
           size(varargin{3},2) == size(varargin{2},2)   % usage: CRKSM(sys, s0_inp, Rt)
            s0_inp = varargin{2};            s0_out  = [];           
            Rt     = varargin{3};            Lt      = [];
            input  = 2;                      pointer = @tangentialV;
        elseif isempty(varargin{2})                     % usage: CRKSM(sys, [], s0_out)
            s0_inp = [];                     s0_out  = varargin{3};  
            Rt     = [];                     Lt      = [];
            input  = 3;                      pointer = @blockW;
        else                                            % usage: CRKSM(sys, s0_inp, s0_out)
            s0_inp = varargin{2};            s0_out  = varargin{3};      
            Rt     = [];                     Lt      = [];           
            input  = 4;                      pointerV = @blockV;
                                             pointerW = @blockW;
        end 
    elseif length(varargin) == 5
        if isempty(varargin{2})                         % usage: CRKSM(sys, [], s0_out, [], Lt)
            s0_inp = [];                     s0_out  = varargin{3};  
            Rt     = [];                     Lt      = varargin{5};
            input  = 5;                      pointer = @tangentialW;
        else                                            % usage: CRKSM(sys, s0_inp, s0_out, Rt, Lt)
            s0_inp = varargin{2};            s0_out  = varargin{3};  
            Rt     = varargin{4};            Lt      = varargin{5};   
            input  = 6;                      pointerV = @tangentialV;
                                             pointerW = @tangentialW;
        end  
    else
        error('Input not compatible with current crksm implementation');
    end

% input of matrices    
elseif length(varargin) > 1

    % set A, B matrix
    A = varargin{1};       D = [];  
    B = varargin{2};   

    % set other inputs
    if nargin == 5
        C      = varargin{3};              E       = varargin{4};  
        s0_inp = varargin{5};              s0_out  = [];
        Rt     = [];                       Lt      = [];
        input  = 1;                        pointer = @blockV;
        if isempty(E),  E = speye(size(A)); end
    elseif nargin == 6
        if isempty(varargin{3})                                   
            C      = zeros(1,size(A,1));   E       = varargin{4};
            s0_inp = varargin{5};          s0_out  = [];
            Rt     = varargin{6};          Lt      = []; 
            input  = 2;                    pointer = @tangentialV;
            if isempty(E),  E = speye(size(A)); end
        elseif isempty(varargin{4}) && size(varargin{6,1}) ~= size(B,2)    
            C      = varargin{3};          E      = varargin{4};
            s0_inp = varargin{5};          s0_out = varargin{6};
            Rt     = [];                   Lt     = []; 
            input  = 4;                    pointerV = @blockV;
                                           pointerW = @blockW; 
            if isempty(E),  E = speye(size(A)); end
            if isempty(s0_inp)
                input  = 3;  
                pointer = @blockW;     
                clear pointerV pointerW
            end
        end
    elseif nargin == 8
            C      = varargin{3};          E      = varargin{4};
            s0_inp = varargin{5};          s0_out = varargin{6};
            Rt     = varargin{7};          Lt     = varargin{8}; 
            input  = 6;                    pointerV = @tangentialV;
                                           pointerW = @tangentialW;
            if isempty(E),       E = speye(size(A));        end
            if isempty(s0_inp)
                input  = 5;  
                pointer = @tangentialW;     
                clear pointerV pointerW
            end                      
    else
        if length(varargin) > 8 || length(varargin) < 3
            error('Wrong input');
        end   
    end
    sys = sss(A,B,C,D,E);       % build sys-object
end

% check shifts and tangential directions, check Opts-field OPts.crksmUsage, check if extended or rational krylov
if  (~exist('s0_inp', 'var') || isempty(s0_inp)) && ...
    (~exist('s0_out', 'var') || isempty(s0_out))
    error('sssMOR:rk:NoExpansionPoints','No expansion points assigned.');
end

% check extended case
if (input == 1 && s0_inp(1,1) == 0 && s0_inp(1,2) == inf) || ...
   (input == 3 && s0_out(1,1) == 0 && s0_out(1,2) == inf) 
    Opts.shifts = 'fixedCyclic';
else
    % sort expansion points & tangential directions
    s0_inp = shiftVec(s0_inp);
    s0old = s0_inp;
    if Opts.real
        s0_inp = cplxpair(s0_inp);  %make sure shifts can be paired 
    else
        s0_inp = sort(s0_inp);
    end

    if ~isempty(Rt)
        if size(Rt,2) ~= length(s0_inp),         error('Inconsistent size of Rt');end
        if size(Rt,1) ~= size(sys.B,2),          error('Inconsistent size of Lt');end
        [~,cplxSorting] = ismember(s0_inp,s0old); 
        Rt = Rt(:,cplxSorting);
        % hiermit wird gepueft, ob tang. richtungen konjugiert komplex sind, unabh?ngig von den shifts  
        %if sum(sum(imag(Rt),2)) ~= 0, error('wrong input'); end 
        %hier wird gepr?ft, ob ein kompl. konjugierter shift auch eine komplex conjugierte tang. Richtung hat
        if size(find(imag(s0_inp)),2) ~= size(find(imag(Rt(1,:))),2)
            error('wrong input');
        end
    end

    if ~isempty(s0_out)
        % sort expansion points & tangential directions
        s0_out = shiftVec(s0_out);
        s0old = s0_out;
        if Opts.real 
            s0_out = cplxpair(s0_out); %make sure shifts can be paired 
        else
            s0_out = sort(s0_out);
        end
        if ~isempty(Lt)
            if size(Lt,2) ~= length(s0_out),   error('Inconsistent size of Lt');end
            if size(Lt,1) ~= size(sys.C,1),    error('Inconsistent size of Lt');end
            [~,cplxSorting] = ismember(s0_out,s0old); 
            Lt = Lt(:,cplxSorting);
            % hiermit wird gepueft, ob tang. richtungen konjugiert komplex sind, unabh?ngig von den shifts  
            %if sum(sum(imag(Lt),2)) ~= 0, error('wrong input'); end
            % hier wird gepr?ft, ob ein kompl. konjugierter shift auch eine komplex conjugierte tang. Richtung hat
            if size(find(imag(s0_out)),2) ~= size(find(imag(Lt(1,:))),2)
                error('wrong input');
            end
        end
    end
    if length(s0_inp) > size(sys.A,1) || length(s0_out) > size(sys.A,1)
        error('sssMOR:arnoldi:reducedOrderExceedsOriginal','The desired reduced order exceeds the original order');
    end

    if ~isempty(s0_inp) && ~isempty(s0_out)
        % check if number of input/output expansion points matches
        if length(s0_inp) ~= length(s0_inp)
            error('Inconsistent length of expansion point vectors.');
        end
    end
end
   
%% RKSM Method
% built first subspace (two columns in case of cplx. conj. shifts) with arnoldi
switch input
    case 1
        [basis1] = arnoldi(sys.E,sys.A,sys.B,s0_inp(1,1:2),Opts);
    case 2
        [basis1] = arnoldi(sys.E,sys.A,sys.B,s0_inp(1,1:2),Rt(:,1:2),Opts);
    case 3
        [basis1] = arnoldi(sys.E',sys.A',sys.C',s0_out(1,1:2),Opts);
    case 4
        [basis1] = arnoldi(sys.E,sys.A,sys.B,s0_inp(1,1:2),Opts);
        [basis2] = arnoldi(sys.E',sys.A',sys.C',s0_out(1,1:2),Opts);
    case 5
        [basis1] = arnoldi(sys.E',sys.A',sys.C',s0_out(1,1:2),Lt(:,1:2),Opts);
    case 6
        [basis1] = arnoldi(sys.E,sys.A,sys.B,s0_inp(1,1:2),Rt(:,1:2),Opts);
        [basis2] = arnoldi(sys.E',sys.A',sys.C',s0_out(1,1:2),Lt(:,1:2),Opts);
end
newdir1 = zeros(size(sys.A,1),size(basis1,2)/2); % declare newdir-variable
newdir2 = newdir1;
basis1_no = basis1;

% check usage of crksm-function
if strcmp(Opts.purpose,'lyapunov')
    usage = @crksmLyap;  
else
    usage = @crksmSysr;
end

clear input Def varargin s0old A B C D E 

% reduction step
if  ~exist('basis2','var')  
    Ar = basis1'*sys.A*basis1;   Br = basis1'*sys.B;   Er = basis1'*sys.E*basis1;   Cr = sys.C*basis1;
    basis2 = [];
else
    Ar = basis2'*sys.A*basis1;   Br = basis2'*sys.B;   Er = basis2'*sys.E*basis1;   Cr = sys.C*basis1;
end
% build ssRed object
sysr = ssRed(Ar,Br,Cr,sys.D,Er);

% call usage handle function for the first solving step
[sysr,data,Opts] = usage(sys,sysr,basis1,1,Opts);

% start iteration
if ~exist('data.out2','var')
    for ii = (size(basis1,2)/size(newdir1,2))+1:1:Opts.maxiter_rksm
        if size(basis1,2) == (ii-1)*size(newdir1,2)

            % get new shifts and tangetial directions
            if (~isempty(s0_inp) && ii > size(s0_inp,2)) || (~isempty(s0_out) && ii > size(s0_out,2))
                if strcmp(Opts.shifts,'cyclic') ||  strcmp(Opts.shifts,'fixedCyclic')
                    s0_inp = repmat(s0_inp,1,2);        Rt = repmat(Rt,1,2); 
                    s0_out = repmat(s0_out,1,2);        Lt = repmat(Lt,1,2); 
                else
                    [s0_inp,s0_out,Rt,Lt] = getShifts(sys,sysr,s0_inp,s0_out,Rt,Lt,basis1,basis2,basis1_no,Opts);
                end
            end

            % get new direction, enlarge basis
            if isempty(basis2) || isempty(basis1)
                newdir1 = pointer(sys,basis1,basis2,s0_inp,s0_out,Rt,Lt,ii,size(newdir1,2));
                if Opts.real == 1
                    basis1 = [basis1 real(newdir1)];       % basis1 is either V or W 
                else
                    basis1 = [basis1 newdir1];
                end
            else
                newdir1 = pointerV(sys,basis1,basis2,s0_inp,s0_out,Rt,Lt,ii,size(newdir1,2));
                newdir2 = pointerW(sys,basis2,basis1,s0_inp,s0_out,Rt,Lt,ii,size(newdir1,2));
                if Opts.real == 1
                    basis1 = [basis1 real(newdir1)];    % basis1 is V  
                    basis2 = [basis2 real(newdir2)];    % basis2 is W
                else
                    basis1 = [basis1 newdir1];    % basis1 is V  
                    basis2 = [basis2 newdir2];    % basis2 is W
                end       
            end                   
       end % end of complex wait-sequence (if isreal(vnew))

       % orthogonalize new basis; note: for loop is neccessary for block
       if isempty(basis2) || isempty(basis1)
           hermite_gram_sch = 0;
           for jj=size(basis1,2)-(size(newdir1,2)-1):1:size(basis1,2)
                basis1 = gramSchmidt(jj,basis1,hermite_gram_sch,Opts);
           end
       else
           if s0_inp == s0_out, hermite_gram_sch = 1; end
           for jj=size(basis1,2)-(size(newdir1,2)-1):1:size(basis1,2)
               [basis1,~,basis2] = gramSchmidt(jj,basis1,basis2,hermite_gram_sch,Opts);
           end
       end
       % reduction step
       jnew = size(basis1,2)-(size(newdir1,2)-1);    jnew_last = size(basis1,2);  
       if ~isempty(basis2)
           Ar(1:jnew-1,jnew:jnew_last)     = basis2(:,1:jnew-1)'*sys.A*basis1(:,jnew:jnew_last);              % 1. step: new columns to the existing matrix Ar      
           Ar(jnew:jnew_last,1:jnew_last)  = basis2(:,jnew:jnew_last)'*sys.A*basis1(:,1:jnew_last);           % 2. step: rows are filled up so that a reduced square matrix results
           Er(1:jnew-1,jnew:jnew_last)     = basis2(:,1:jnew-1)'*sys.E*basis1(:,jnew:jnew_last);
           Er(jnew:jnew_last,1:jnew_last)  = basis2(:,jnew:jnew_last)'*sys.E*basis1(:,1:jnew_last);
           Br(jnew:jnew_last,:)            = basis2(:,jnew:jnew_last)'*sys.B;
           Cr(:,jnew:jnew_last)            = sys.C*basis1(:,jnew:jnew_last);
       else
           Ar(1:jnew-1,jnew:jnew_last)     = basis1(:,1:jnew-1)'*sys.A*basis1(:,jnew:jnew_last);        
           Ar(jnew:jnew_last,1:jnew_last)  = basis1(:,jnew:jnew_last)'*sys.A*basis1(:,1:jnew_last);   
           Er(1:jnew-1,jnew:jnew_last)     = basis1(:,1:jnew-1)'*sys.E*basis1(:,jnew:jnew_last);
           Er(jnew:jnew_last,1:jnew_last)  = basis1(:,jnew:jnew_last)'*sys.E*basis1(:,1:jnew_last);
           Br(jnew:jnew_last,:)            = basis1(:,jnew:jnew_last)'*sys.B;
           Cr(:,jnew:jnew_last)            = sys.C*basis1(:,jnew:jnew_last);
       end
       sysr = ssRed(Ar,Br,Cr,sys.D,Er);     % ssRed-object

       % call usage handle function for Lyapunov/sysr
       [sysr,data,Opts] = usage(sys,sysr,basis1,ii,Opts);
       
       % check shift capacity
       if mod(ii,10) == 1 && strcmp(Opts.shifts,'dynamical') && (data.out1(ii,1)-data.out1(ii-10,1)>Opts.shiftTol)...
          && (data.out1(ii,1)) ~= inf % reuse old shifts
           Opts.shifts = 'cyclic';
       elseif mod(ii,10) == 1 && ~strcmp(Opts.shifts,'fixedCyclic') && (data.out1(ii,1)-data.out1(ii-10,1)<Opts.shiftTol)... 
              && (data.out1(ii,1)) ~= inf % get new shifts
           Opts.shifts = 'dynamical';   
       end

       % quit for loop, show information
       if ~isempty(data.out2), break;   end

       % enlarge subspace with the imaginary part of the direction and eventually calculate new real directions
       if Opts.real == 1 && (~isreal(newdir1) || ~isreal(newdir2))
           if isempty(basis2) || isempty(basis1)
               basis1 = [basis1 imag(newdir1)];
           else
               if ~isreal(newdir1) && ~isreal(newdir2)
                   basis1 = [basis1 imag(newdir1)];
                   basis2 = [basis2 imag(newdir2)];
               elseif ~isreal(newdir1) && isreal(newdir2)
                   basis1 = [basis1 imag(newdir1)];

                   % calculate new real direction, if last vnew is real and last wnew is complex 
                   newdir2 = pointerW(sys,basis2,basis1,s0_inp,s0_out,Rt,Lt,ii,size(newdir1,2));
                   basis2 = [basis2 real(newdir2)];

               elseif isreal(newdir1) && ~isreal(newdir2)
                   basis2 = [basis2 imag(newdir2)];

                   % calculate new real direction, if last vnew is real and last wnew is complex
                   newdir1 = pointerV(sys,basis1,basis2,s0_inp,s0_out,Rt,Lt,ii,size(newdir1,2));
                   basis1 = [basis1 real(newdir1)];
               end
           end
           % empty newdir of basis1 and/or basis2 
           newdir1 = zeros(size(sys.A,1),size(newdir1,2));
           newdir2 = newdir1;
       end
    end  % end of for loop

    % create output
    data.out4 = basis1;
    data.out5 = basis2;
    data.out6 = s0_inp;
    if ii == Opts.maxiter_rksm
        warning('\n maximum number of iterations is reached without converging!' )
    end
end % end of: "if ~exist('data.out2','var')"

% clear global variables
clear global hermite_gram_sch
end


%% ***************************** AUXILIARY ********************************
function [V, TRv, W, TLw] = gramSchmidt(jCol,V,varargin)
%   Gram-Schmidt orthonormalization
%   Input:  jCol:  Column to be treated
%           V, W:  Krylov-Subspaces
%   Output: V, W:  orthonormal basis of Krylov-Subspaces
%           TRv, TLw: Transformation matrices

% input
if isa(varargin{end},'struct')
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end

if length(varargin) == 1
    hermite_gram_sch = varargin{1};
    W = eye(size(V));
elseif length(varargin) == 2
    W = varargin{1};
    hermite_gram_sch = varargin{2};
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

function [vnew] = blockV(sys,V,~,s0_inp,~,~,~,iter,colIndex)
    rhsB = V(:,size(V,2)-(colIndex-1):size(V,2));
    if s0_inp(1,iter) ~= s0_inp(1,iter-1) && s0_inp(1,iter) ~= conj(s0_inp(1,iter-1))
       vnew = solveLse(sys.A,rhsB,sys.E,s0_inp(1,iter));
    else
       Opts.reuseLU = 1;
       vnew = solveLse(sys.A,rhsB,sys.E,s0_inp(1,iter),Opts);
    end
    vnew = vnew(:,size(rhsB,2)); % because solution of solveLse comes sometimes with one zero-column  
end

function [wnew] = blockW(sys,W,~,~,s0_out,~,~,iter,colIndex)
    rhsC = W(:,size(W,2)-(colIndex-1):size(W,2));
    if s0_out(1,iter) ~= s0_out(1,iter-1) && s0_out(1,iter) ~= conj(s0_out(1,iter-1))
        wnew(:,size(W,2)) = solveLse(sys.A,rhsC,sys.E,s0_out(1,iter)); 
    else
        Opts.reuseLU = 1;
        wnew(:,size(W,2)) = solveLse(sys.A,rhsC,sys.E,s0_out(1,iter),Opts);
    end
end

function [vnew] = tangentialV(sys,V,~,s0_inp,~,Rt,~,iter,colIndex)
    %rhsB = E*V(:,size(V,2)-(colIndex-1):size(V,2))*rt;
    if s0_inp(1,iter) ~= s0_inp(1,iter-1) && s0_inp(1,iter) ~= conj(s0_inp(1,iter-1))
        vnew(:,size(V,2)) = solveLse(sys.A,sys.B*Rt(:,iter),sys.E,s0_inp(1,iter));
    else
        Opts.reuseLU = 1;
        vnew(:,size(V,2)) = solveLse(sys.A,sys.B*Rt(:,iter),sys.E,s0_inp(1,iter),Opts);
    end
end

function [wnew] = tangentialW(sys,~,W,~,s0_out,~,Lt,iter,colIndex)
    %rhsC = E'*W(:,size(W,2)-(colIndex-1):size(W,2))*lt; 
    if s0_out(1,iter) ~= s0_out(1,iter-1) && s0_out(1,iter) ~= conj(s0_out(1,iter-1))
        wnew(:,size(W,2)) = solveLse(sys.A,sys.C'*Lt(:,iter),sys.E,s0_out(1,iter));
    else
        Opts.reuseLU = 1;
        wnew(:,size(W,2)) = solveLse(sys.A,sys.C'*Lt(:,iter),sys.E,s0_out(1,iter),Opts);
    end
end

function [sysr,data,Opts] = crksmLyap(sys,sysr,basis1,iter,Opts)
    % set persistent variables 
    persistent S Rnorm
    if exist('S','var'), S_last = S; else,  S_last = []; end
    % try to solve Lyapunov equation first with lyapchol or second with lyap
    try
       S = lyapchol(sysr.A,sysr.B,sysr.E);
       if nnz(sysr.C) ~= 0
           R = lyapchol(sysr.A',sysr.C',sysr.E');
       else
           R = [];
       end
    catch
       S = lyap(sysr.A,sysr.B*sysr.B',[],sysr.E);
       if nnz(sysr.C) ~= 0
           R = lyap(sysr.A',sysr.C'*sysr.C,[],sysr.E');
       else
           R = [];
       end
    end
        
    % choose computation of norm/stopping criteria
    if strcmp(Opts.residual,'residual_lyap')
       % test determination (computation of residual after Panzer/Wolff), compute Er^-1*Br and Er^-1*Ar and other factors
       Er_inv_Br = solveLse(sysr.E,sysr.B);
       Opts.reuseLU = 1;
       Er_inv_Ar = solveLse(sysr.E,sysr.A,Opts);   

       % compute factors from residual
       Opts.Bbot = sys.B-(sys.E*basis1)*Er_inv_Br;
       Cr_hat = solveLse(Opts.Bbot'*Opts.Bbot,Opts.Bbot'*((sys.A*basis1)-(sys.E*basis1)*Er_inv_Ar));
       F = sys.E*basis1*(Er_inv_Br+(S'*S)*Cr_hat');

       % compute residual norm (Euclidean Norm)
       if strcmp(Opts.rksmnorm, 'H2')
           res0  = norm(sysr.B' * sysr.B,2);
           Rnorm(iter,1) = max(abs(eig(full([Opts.Bbot'*Opts.Bbot+Opts.Bbot'*F, Opts.Bbot'*Opts.Bbot; F'*Opts.Bbot+F'*F, F'*Opts.Bbot])))) / res0; 
       else
           % Frobenius Norm
           res0  = norm(sysr.B' * sysr.B,'fro');
           Rnorm(iter,1) = sqrt(sum(eig(full([Opts.Bbot'*Opts.Bbot+Opts.Bbot'*F, Opts.Bbot'*Opts.Bbot; F'*Opts.Bbot+F'*F, F'*Opts.Bbot])))^2) / res0;
       end
    else
        
       if strcmp(Opts.rksmnorm, 'H2')
           X_lastnorm = NormFrobEfficient(1,S_last);
           X_norm = NormFrobEfficient(1,S);
       else
           X_lastnorm = NormFrobEfficient(0,S_last);
           X_norm = NormFrobEfficient(0,S);
       end
       Rnorm(iter,1) = abs(X_norm-X_lastnorm);
     end
    data.out1 = Rnorm; data.out2 = []; data.out3 = [];
    
    % stop program
   if Rnorm(iter,1) < Opts.rctol
       if Opts.lowrank == 1            % compute low rank factor of solution
           S = basis1*S;    
           if ~isempty('R','var'),  R = basis1*R;  end
       end
       data.out2 = S;  data.out3 = R;  % wirte output in info-struct
       % show information of programme
       if ~exist('R','var')
           fprintf('\n RKSM, usage Lyapunov, step: %d \t Convergence\n last residual norm: %d, tolerance: %d,\n' ,iter,Rnorm(iter,1),Opts.rctol);
       else
           fprintf('\n RKSM, usage Lyapunov, step \t %d \t Convergence, last residual: \t %d, tolerance: \t %d,\n',iter,Rnorm(iter,1),Opts.rctol);
       end
   elseif size(S,2) == size(sys.A,2)
       disp('\n V has reached the dimension of the original System without converging!');
   end
end

function [sysr,data,Opts] = crksmSysr(~,sysr,~,iter,Opts)
    persistent stopCrit sysr_last
    if iter == 1, sysr_last = sss([],[],[]);  end
    % stopping criteria
    if all(real(eig(sysr))<0) && all(real(eig(sysr_last))<0)
       stopCrit(iter,1)=norm(sysr-sysr_last)/norm(sysr);
    else
       stopCrit(iter,1) = inf; %initialize in case the reduced model is unstable 
    end
    
    data.out1 = stopCrit; data.out2 = [];
    
    if stopCrit(iter,1) < Opts.rctol
       data.out1(2,1) = data.out1(1,1);
       data.out2 = 1;
       fprintf('\n RKSM, usage MOR, step: %d \t Convergence\n last system norm: %d, tolerance: %d,\n' ,iter,stopCrit(iter,1),Opts.rctol);
    end
    sysr_last = sysr;
end


