classdef ssRed < ss
% SSRED - Reduced state-space LTI system (ssRed) class
%
% Syntax:
%       sysr = ssRed(method,params,A,B,C)
%       sysr = ssRed(method,params,A,B,C,paramsList)
%       sysr = ssRed(method,params,A,B,C,D)
%       sysr = ssRed(method,params,A,B,C,D,paramsList)
%       sysr = ssRed(method,params,A,B,C,D,E)
%       sysr = ssRed(method,params,A,B,C,D,E,paramsList)
%       sysr = ssRed(method,params,sys_sss)
%       sysr = ssRed(method,params,sys_sss,paramsList)
%       sysr = ssRed(method,params,sys_ss)
%       sysr = ssRed(method,params,sys_ss,paramsList)
%       sysr = ssRed(method,params,sys_ssRed)
%
% Description:
%       This class is derived from the ss-class. It is used to represent
%       models which are the result of reducing the order of large models
%       with specific algorithms. 
%
% Input Arguments:
%       -method: name of the used reduction algorithm;
%                ['tbr' / 'modalMor' / 'irka' / 'rk' / 'projectiveMor' / 'porkV' / 'porkW' / 'cure_spark' / 'cure_irka' / 'cure_rk+pork']
%       -params (tbr):          structure with the parameters for the tbr-algorithm;
%           -.originalOrder:    Model order before reduction;
%           -.type:             select amongst different tbr algorithms
%                               ['tbr' / 'adi' / 'matchDcGain' ]
%           -.redErr:           upper bound of reduction error
%                               ['0' / positive float]
%           -.hsvTol:           tolerance for Hankel-Singular values
%                               ['1e-15' / positive float]
%           -.warnOrError:      display warnings or errors
%                               ['warn' / 'error' / '0']
%           -.lse:              solve linear system of equations (only for adi)
%                               ['gauss' / 'luChol']
%           -.hsv:              Hankel singular values
%       -params (modalMor):     structure with the parameters for the
%                               modalMor-algorithm
%           -.originalOrder:    Model order before reduction
%           -.type:             option to eigs command;
%                               ['SM' / 'LM' / 'SA' / 'LA' / 'SR' / 'LR' / real or complex scalar]
%           -.orth:             orhtogonalization;
%                               ['0' / 'qr']
%           -.real:             real reduced system;
%                               ['real' / '0']
%           -.tol:              tolerance for the eigenspace computation;
%                               [positive float]
%           -.dominance:        perform dominance analysis;
%                               ['0' / 'analyze' / '2q' / '3q' / '4q']
%       -params (irka):         structure with the parameters for the
%                               irka-algorithm
%           -.originalOrder:    Model order before reduction
%           -.maxiter:          maximum number of iterations;
%                               [50 / positive integer]
%           -.tol:              convergence tolerance;
%                               [1e-3 / positive float]
%           -.type:             choose between different irka modifications;
%                               ['' / 'stab']
%           -.verbose:          show text output during iterations;
%                               [0 / 1]
%           -.stopCrit:	stopping criterion;
%                               ['combAny' / 's0' / 'sysr' / 'combAll']
%           -.suppressverbose: suppress any type of verbose for speedup;
%                               [0 / 1]
%           -.orth:             orthogonalization of new projection direction
%                               ['2mgs' / 0 / 'dgks' / 'mgs']
%           -.reorth:           reorthogonalization
%                               ['gs' / 0 / 'qr']
%           -.lse:              use LU or hessenberg decomposition
%                               ['sparse' / 'full' / 'hess']
%           -.dgksTol:          tolerance for dgks orthogonalization
%                               [1e-12 / positive float]
%           -.krylov:           standard or cascaded krylov basis
%                               [0 / 'cascade]
%           -.kIter:            number of iterations
%           -.s0:               final choice of shifts
%           -.Rt:               final choice of right tangential directions
%                               for MIMO
%           -.Lt:               final choice of left tangential directions
%                               for MIMO
%           -.s0Traj:           trajectory of all shift
%           -.RtTraj:           trajectory of right tangential directions
%           -.LtTraj:           trajectory of left tangential directions
%       -params (rk):           structure with the parameters for the
%                               rk-algorithm
%           -.originalOrder:    Model order before reduction
%           -.real:             keep the projection matrices real
%                               [true / false]
%           -.orth:             orthogonalization of new projection direction
%                               ['2mgs' / 0 / 'dgks' / 'mgs']
%           -.reorth:           reorthogonalization
%                               ['gs' / 0 / 'qr']
%           -.lse:              use LU or hessenberg decomposition
%                               ['sparse' / 'full' / 'hess']
%           -.dgksTol:          tolerance for dgks orthogonalization
%                               [1e-12 / positive float]
%           -.krylov:           standard or cascaded krylov basis
%                               [0 / 'cascade]
%           -.IP:               inner product
%           -.s0_inp:           expansion points for input Krylov subspaces
%           -.s0_out:           expansion points for output Krylov
%                               subspaces
%           -.Rt:               right tangential directions for MIMO
%           -.Lt:               left tangential directions for MIMO
%       -params (projectiveMor):structure with the parameters for the
%                               projectiveMor-algorithm
%           -.originalOrder:    Model order before reduction
%           -.trans:            choose how W should be transposed
%                               [ {T} / H ]
%       -params (porkV):        structure with the parameters for the
%                               porkV-algorithm
%           -.originalOrder:    Model order before reduction
%       -params (porkW):        structure with the parameters for the
%                               porkW-algorithm
%           -.originalOrder:    Model order before reduction
%       -params (cure):         structure with the parameters for the
%                               cure-algorithm
%           -.originalOrder:    Model order before reduction
%           -.cure.redfun:      reduction algorithm
%                               [{'spark'} / 'irka' / 'rk+pork']
%           -.cure.nk:          reduced order at each iteration 
%                               [{'2'} / positive integer]
%           -.cure.fact:        factorization mode 
%                               [{'V'} / 'W']
%           -.cure.init:        shift initialization mode 
%                               [{'sm'} / 'zero' / 'lm' / 'slm']
%           -.cure.stop:        stopping criterion
%                               [{'nmax'} / 'h2Error']
%           -.cure.stopval:     value according to which the stopping criterion is evaluated
%                               [{'round(sqrt(sys.n))'} / positive integer]
%           -.cure.verbose:     display text during cure 
%                               [{'0'} / '1']
%           -.cure.SE_DAE:      reduction of index 1 semiexplicit DAE 
%                               [{'0'} / '1']
%           -.cure.test:        execute analysis code 
%                               [{'0'} / '1']
%           -.cure.gif:         produce a .gif file of the CURE iteration
%                               [{'0'} / '1']
%           -.cure.maxIter:     maximum number of CURE iterations
%                               [{'20'} / positive integer]
%           -.warn:             show warnings
%                               [{'0'} / '1']
%           -.w:                frequencies for analysis plots
%                               [{''} / '{wmin,wmax}' / vector of frequencies]
%           -.zeroThers:        value that can be used to replace 0 
%                               [{'1e-4'} / postivie float]
%       -A: system matrix
%       -B: input matrix
%       -C: output matrix
%       -D: static gain matrix
%       -E: descriptor matrix
%       -sys_ss:    control system toolbox state-space (ss)-object
%       -sys_sss:   sparse state-space (sss)-object
%       -paramsList{i}: Cell-Array of structs. Each field of the cell-array 
%                       represents one reduction (without the current reduction).
%           -.method:           name of the used reduction algorithm (see above)                               
%           -.params:           structure with the parameters of the used
%                               algorithm (see above)
%
% Output Arguments:
%       -sys: reduced state-space (ssRed)-object
%
% Examples:
%		This code creates an instance of the ssRed-class
%
%> A = rand(10);
%> B = rand(10,1);
%> C = rand(2,10);
%> params.originalOrder = 20;
%> params.real = true;
%> sysr = ssRed('rk',params,A,B,C);
%
% See Also: 
%        ss, dss, sss
%
% References:
%		* *[1] Documentation of the Control System Toolbox from MATLAB*
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Niklas Kochdumper
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  15 Jun 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------
    
    properties(SetAccess = private)
        reductionParameters
    end
    properties(Dependent, Hidden)
        n,p,m
        isSiso, isSimo, isMiso, isMimo, isBig
        isDae, isDescriptor
    end
    properties(Hidden)
        HankelSingularValues
        TBal, TBalInv
        ConGram, ConGramChol
        ObsGram, ObsGramChol
        
        isSym
    end
    
    methods
        function obj = ssRed(varargin)
            % parse and check the arguments
            if nargin < 3 || nargin > 8 || nargin == 4
                error('Invalid syntax for the "ssRed" command. Type "help ssRed" for more information.');
            end
            paramsList = [];
            nargin_new = nargin;
            if ~isnumeric(varargin{nargin}) && ~isa(varargin{nargin},'ss') && ...
               ~isa(varargin{nargin},'sss')         %paramsList specified
                paramsList = varargin{nargin};
                nargin_new = nargin_new-1;
            end
            if nargin_new == 3
                sys = varargin{3};
                if isa(sys,'ssRed')
                    if ~isempty(paramsList)
                        error('Invalid syntax for the "ssRed" command. Type "help ssRed" for more information.');
                    end
                    paramsList = sys.reductionParameters;
                end
                if isa(sys,'ss')    %ss-objects or ssRed-objects
                   A = sys.a;
                   B = sys.b;
                   C = sys.c;
                   D = sys.d;
                   E = sys.e;
                elseif isa(sys,'sss')   %sss-objects
                   A = full(sys.a);
                   B = full(sys.b);
                   C = full(sys.c);
                   D = full(sys.d);
                   E = full(sys.e);
                else
                   error('The third argument has to be an object of type "ss" or "sss"'); 
                end
            else
                A = full(varargin{3});
                B = full(varargin{4});
                C = full(varargin{5});
                
                if nargin_new == 6
                    D = full(varargin{6});
                    E = [];
                elseif nargin_new == 7
                    D = full(varargin{6});
                    E = full(varargin{7});
                else                
                    D = [];
                    E = [];
                end             
            end
            
            if ~isa(varargin{1},'char') || ismember(varargin{1},{'tbr','modalMor','rk','irka','projectiveMor','porkV','porkW','cure'}) == 0
                error('The first argument has a wrong format. Type "help ssRed" for more information.');
            end
            
            % call the construktor of the superclass ss
            obj@ss(A,B,C,D);
            obj.e = E;
            
            % check the params struct
            try
               obj.checkParamsStruct(varargin{2},varargin{1});
            catch ex
               error(ex.message); 
            end
            
            % check all fields of paramsList
            if ~isempty(paramsList)
                try
                   for i = 1:size(paramsList,1)
                      obj.checkParamsStruct(paramsList{i}.params,paramsList{i}.method);
                      if length(fieldnames(paramsList{i})) ~= 2
                          error('The argument "paramsList" has the wrong format. Type "help ssRed" for more information.');
                      end
                      paramsList{i}.params = obj.parseParamsStruct(paramsList{i}.params,paramsList{i}.method);
                   end
                catch ex
                   error('The argument "paramsList" has the wrong format. Type "help ssRed" for more information.'); 
                end
            end
            
            % update the reductionParameters list
            if isempty(paramsList)
               obj.reductionParameters = cell(1,1);
               obj.reductionParameters{1,1}.method = varargin{1};
               obj.reductionParameters{1,1}.params = obj.parseParamsStruct(varargin{2},varargin{1});
            else
               len = size(paramsList,1);
               obj.reductionParameters = paramsList;
               obj.reductionParameters{len+1,1}.method = varargin{1};
               obj.reductionParameters{len+1,1}.params = obj.parseParamsStruct(varargin{2},varargin{1});
            end
            
            % remove "projectiveMor" from the reduction history
            obj.reductionParameters = obj.removeProjectiveMor(obj.reductionParameters);
        end
        
        
        %% Get Basic Properties
        function m = get.m(sys) % number of inputs
            m = size(sys.b,2);
        end
        function n = get.n(sys) % system order
            n = size(sys.a,1);
        end
        function p = get.p(sys) % number of outputs
            p = size(sys.c,1);
        end
        
        %% Get helper functions
        function isSiso = get.isSiso(sys); isSiso=(sys.p==1)&&(sys.m==1); end
        function isSimo = get.isSimo(sys); isSimo=(sys.p>1)&&(sys.m==1); end
        function isMiso = get.isMiso(sys); isMiso=(sys.p==1)&&(sys.m>1); end
        function isMimo = get.isMimo(sys); isMimo=(sys.p>1)||(sys.m>1); end
        function isBig = get.isBig(sys); isBig=(sys.n>5000);end
        
        function isDae = get.isDae(sys)
            if condest(sys.e)==Inf
                isDae = 1;
            else
                isDae = 0;
            end
        end
        
        function isDescriptor = get.isDescriptor(sys)
            isDescriptor = logical(full(any(any(sys.e-speye(size(sys.e))))));
        end
        
        function sys = resolveDescriptor(sys)
            sys.a = sys.e\sys.a;
            sys.b = sys.e\sys.b;
            sys.e = [];
        end
        
        function isSym = get.isSym(sys) %A=A', E=E'
            if isequal(sys.isSym,0) || isequal(sys.isSym,1)
                isSym = sys.isSym;
            else
                if max(max(sys.a-sys.a.'))<1e-6 && max(max(sys.e-sys.e.'))<1e-6
                    isSym = 1;
                else
                    isSym = 0;
                end
            end
        end        
        
        %% Override operators and build-in-functions
        function varargout = eig(sys, varargin)
            if sys.isBig
                warning(['System order is very large: ',num2str(sys.n),'. You may want to try eigs(sys) instead.'])
            end

            if sys.isDescriptor
                if nargout==1||nargout==0
                    [varargout{1}] = eig(full(sys.a), full(sys.e),varargin{:});
                elseif nargout == 2
                    [varargout{1}, varargout{2}] = eig(full(sys.a), full(sys.e),varargin{:});
                elseif nargout == 3
                    [varargout{1}, varargout{2}, varargout{3}]  = eig(full(sys.a), full(sys.e),varargin{:});
                end
            else
                if nargout==1||nargout==0
                    [varargout{1}] = eig(full(sys.a),varargin{:});
                elseif nargout == 2
                    [varargout{1}, varargout{2}] = eig(full(sys.a),varargin{:});
                elseif nargout == 3
                    [varargout{1}, varargout{2}, varargout{3}]  = eig(full(sys.a),varargin{:});
                end
            end
        end
        
        function diff = minus(sys1, sys2)     
            if size(sys1.a,1) == 0
                if isa(sys2,'sss')
                    diff = sss(sys2.a, sys2.b, -sys2.c, sys2.d, sys2.e);
                else
                    diff = dss(sys2.a, sys2.b, -sys2.c, sys2.d, sys2.e);
                end
                return
            end
            if size(sys2.a,1) == 0
                diff = dss(sys1.a, sys1.b, sys1.c, sys1.d, sys1.e);
                return
            end

            if size(sys1.b,2) ~= size(sys2.b,2)
                error('sys1 and sys2 must have same number of inputs.')
            end
            if size(sys1.c,1) ~= size(sys2.c,1)
                error('sys1 and sys2 must have same number of outputs.')
            end
            if isa(sys1,'ss') && isempty(sys1.e)
                sys1.e = eye(size(sys1.a,1));
            end
            if isa(sys2,'ss') && isempty(sys2.e)
                sys2.e = eye(size(sys2.a,1));
            end

            if isa(sys2,'sss')           
                diff = sss([sys1.a sparse(size(sys1.a,1),size(sys2.a,1)); sparse(size(sys2.a,1),size(sys1.a,1)) sys2.a], ...
                    [sys1.b; sys2.b], ...
                    [sys1.c, -sys2.c], ...
                    sys1.d - sys2.d, ...
                    [sys1.e sparse(size(sys1.a,1),size(sys2.a,1)); sparse(size(sys2.a,1),size(sys1.a,1)) sys2.e]);
            else
                diff = dss([sys1.a zeros(size(sys1.a,1),size(sys2.a,1)); zeros(size(sys2.a,1),size(sys1.a,1)) sys2.a], ...
                    [sys1.b; sys2.b], ...
                    [sys1.c, -sys2.c], ...
                    sys1.d - sys2.d, ...
                    [sys1.e zeros(size(sys1.a,1),size(sys2.a,1)); zeros(size(sys2.a,1),size(sys1.a,1)) sys2.e]);
            end
        end
        
        function [r,p,d] = residue(sys, Opts)       
            Def.rType = 'res';

            if ~exist('Opts','var') || isempty(Opts)
                Opts = Def;
            else
                Opts = parseOpts(Opts,Def);
            end

            %perform eigen-decomposition of system
            try
                [T,J] = eig(sys);
            catch err
                error('Computation of the eigenvalues and eigenvectors failed with message:%s',err.message);
            end

            % transform system to diagonal form
            p=diag(J).';
            if issparse(T)
                rcondNumber = 1/condest(T);
            else
                rcondNumber=rcond(T);
            end
            if rcondNumber<eps
                warning(['Matrix of eigenvectors is close to singular or badly scaled. Results may be inaccurate. RCOND =',num2str(rcondNumber)]);
            end
            B=(sys.e*T)\sys.b;
            C=sys.c*T;
            d=sys.d;

            % calculate residues
            if strcmp(Opts.rType,'dir')
                % return the residual directions instead of the residuals
                r = {C, B};  
            else
                % return the residuals
                r = cell(1,sys.n);
                for i=1:sys.n
                    r{i} = full(C(:,i)*B(i,:));
                end
            end
        end
        
        function tmax = decayTime(sys)
            [res,p]=residue(sys);

            % is system stable?
            if any(real(p)>0 & real(p)<1e6) % larger than 0 but, smaller than infinity-threshold
                % no -> tmax=NaN
                tmax=NaN;
                warning('sss:decayTime:UnstableSys','The system is not stable. The decay time is set to tmax=NaN.');
                return
            end

            tmax=0; temp = cat(3,res{:}); 
            for i=1:sys.p
                for j=1:sys.m
                    % how much does each pole contribute to energy flow?
                    h2=zeros(size(p));
                    for k=1:length(p)
                        %we need the siso residual for all poles into on vectors
                        resIJvec = squeeze(temp(i,j,:)).';
                        h2(k)=res{k}(i,j)*sum(resIJvec./(-p(k)-p));
                    end

                    [h2_sorted, I] = sort(real(h2));
                    % which pole contributes more than 1% of total energy?
                    I_dom = I( h2_sorted > 0.01*sum(h2_sorted) );
                    if isempty(I_dom)
                        % no poles are dominant, use slowest
                        I_dom = 1:length(p);
                    end
                    % use slowest among dominant poles
                    [h2_dom, I2] = sort(abs(real(p(I_dom))));

                    % when has slowest pole decayed to 1% of its maximum amplitude?
                    tmax=max([tmax, log(100)/abs(real(p(I_dom(I2(1)))))]);
                end
            end
        end
        
        function [issd, numericalAbscissa] = issd(sys)
            %  Parse input
            if condest(sys.e)>1e16, error('issd does not support DAEs'),end

            %  Perform computations
            % E >0?
            isPosDef = ispd(sys.e);
            if ~isPosDef
                if nargout == 0, warning('System is not strictly dissipative (E~>0).'); 
                else issd = 0; end
                return
            end

            % A + A' <0?  
            isNegDef = ispd(-sys.a-sys.a');
            if isNegDef
                if nargout == 0, fprintf('System is strictly dissipative.\n'); else issd = 1; end
            else
                if nargout == 0, warning('System is not strictly dissipative (E>0, A+A''~<0)'); else issd = 0; end
            end

            if nargout==2 % computation of the numerical abscissa required
                p    = 20;		% number of Lanczos vectors
                tol  = 1e-10;	% convergence tolerance
                opts = struct('issym',true, 'p',p, 'tol',tol, 'v0',sum(sys.e,2));
                try
                    numericalAbscissa = eigs((sys.a+sys.a')/2, sys.e, 1, 'la', opts);
                catch err
                    warning('Computation of the numerical abscissa failed with message:%s',err.message);
                    numericalAbscissa = NaN;
                end
            end
        end

    end
    
    %%Private and static helper methods
    methods(Hidden, Access = private, Static)
        
        function correct = checkParamsStruct(params,method)
        % Checks if the struct with the parameters has the correct
        % structure
            correct = 1;
            
            % check original Order
            if ~isfield(params,'originalOrder')
                  error('Struct params does not contain the field "originalOrder"!')
            end           
            if mod(params.originalOrder,1)~=0
                  error('params.originalOrder: Wrong value. Integer value required!'); 
            end
            
            % check algorithm-specific parameters
            if strcmp(method,'tbr')                     %tbr
               if ~isfield(params,'type')
                   error('Struct params does not contain the field "type"!');
               elseif ~isfield(params,'redErr')
                   error('Struct params does not contain the field "redErr"!');
               elseif ~isfield(params,'hsvTol')
                   error('Struct params does not contain the field "hsvTol"!');
               elseif ~isfield(params,'warnOrError')
                   error('Struct params does not contain the field "warnOrError"!');
               elseif ~isfield(params,'lse')
                   error('Struct params does not contain the field "lse"!');
               elseif ~isfield(params,'hsv')
                   error('Struct params does not contain the field "hsv"!');
               end                    
            elseif strcmp(method,'modalMor')            %modalMor          
               if ~isfield(params,'type')
                   error('Struct params does not contain the field "type"!');
               elseif ~isfield(params,'orth')
                   error('Struct params does not contain the field "orth"!');
               elseif ~isfield(params,'real')
                   error('Struct params does not contain the field "real"!');
               elseif ~isfield(params,'tol')
                   error('Struct params does not contain the field "tol"!');
               elseif ~isfield(params,'dominance')
                   error('Struct params does not contain the field "dominance"!');
               end
            elseif strcmp(method,'irka')                %irka
               if ~isfield(params,'maxiter')
                   error('Struct params does not contain the field "maxiter"!');
               elseif ~isfield(params,'tol')
                   error('Struct params does not contain the field "tol"!');
               elseif ~isfield(params,'type')
                   error('Struct params does not contain the field "type"!');
               elseif ~isfield(params,'verbose')
                   error('Struct params does not contain the field "verbose"!');
               elseif ~isfield(params,'stopCrit')
                   error('Struct params does not contain the field "stopCrit"!');
               elseif ~isfield(params,'suppressverbose')
                   error('Struct params does not contain the field "suppressverbose"!');
               elseif ~isfield(params,'orth')
                   error('Struct params does not contain the field "orth"!');
               elseif ~isfield(params,'lse')
                   error('Struct params does not contain the field "lse"!');
               elseif ~isfield(params,'dgksTol')
                   error('Struct params does not contain the field "dgksTol"!');
               elseif ~isfield(params,'krylov')
                   error('Struct params does not contain the field "krylov"!');
               elseif ~isfield(params,'s0')
                   error('Struct params does not contain the field "s0"!');
               elseif ~isfield(params,'Rt')
                   error('Struct params does not contain the field "Rt"!');
               elseif ~isfield(params,'Lt')
                   error('Struct params does not contain the field "Lt"!');
               elseif ~isfield(params,'kIter')
                   error('Struct params does not contain the field "kIter"!');
               elseif ~isfield(params,'s0Traj')
                   error('Struct params does not contain the field "s0Traj"!');
               elseif ~isfield(params,'RtTraj')
                   error('Struct params does not contain the field "RtTraj"!');
               elseif ~isfield(params,'LtTraj')
                   error('Struct params does not contain the field "LtTraj"!');
               end    
            elseif strcmp(method,'rk')                  %rk
               if ~isfield(params,'real')
                   error('Struct params does not contain the field "real"!');
               elseif ~isfield(params,'orth')
                   error('Struct params does not contain the field "orth"!');
               elseif ~isfield(params,'reorth')
                   error('Struct params does not contain the field "reorth"!');
               elseif ~isfield(params,'lse')
                   error('Struct params does not contain the field "lse"!');
               elseif ~isfield(params,'dgksTol')
                   error('Struct params does not contain the field "dgksTol"!');
               elseif ~isfield(params,'krylov')
                   error('Struct params does not contain the field "krylov"!');
               elseif ~isfield(params,'IP')
                   error('Struct params does not contain the field "IP"!');
               elseif ~isfield(params,'Rt')
                   error('Struct params does not contain the field "Rt"!');
               elseif ~isfield(params,'Lt')
                   error('Struct params does not contain the field "Lt"!');
               elseif ~isfield(params,'s0_inp')
                   error('Struct params does not contain the field "s0_inp"!');
               elseif ~isfield(params,'s0_out')
                   error('Struct params does not contain the field "s0_out"!');
               end
            elseif strcmp(method,'projectiveMor')       %projectiveMor
               if ~isfield(params,'trans')
                   error('Struct params does not contain the field "trans"!');
               end
            elseif strcmp(method,'porkV')               %porkV
            elseif strcmp(method,'porkW')               %porkW
            elseif strcmp(method,'cure')
               if ~isfield(params,'cure')
                   error('Struct params does not contain the field "cure"!');
               elseif ~isfield(params,'warn')
                   error('Struct params does not contain the field "warn"!');
               elseif ~isfield(params,'w')
                   error('Struct params does not contain the field "w"!');
               elseif ~isfield(params,'zeroThres')
                   error('Struct params does not contain the field "zeroThres"!');
               end
               
               paramCure = params.cure;               
               
               if ~isfield(paramCure,'redfun')
                   error('Struct params does not contain the field "cure.redfun"!');
               elseif ~isfield(paramCure,'nk')
                   error('Struct params does not contain the field "cure.nk"!');
               elseif ~isfield(paramCure,'fact')
                   error('Struct params does not contain the field "cure.fact"!');
               elseif ~isfield(paramCure,'init')
                   error('Struct params does not contain the field "cure.init"!');
               elseif ~isfield(paramCure,'stop')
                   error('Struct params does not contain the field "cure.stop"!');
               elseif ~isfield(paramCure,'stopval')
                   error('Struct params does not contain the field "cure.stopval"!');
               elseif ~isfield(paramCure,'verbose')
                   error('Struct params does not contain the field "cure.verbose"!');
               elseif ~isfield(paramCure,'SE_DAE')
                   error('Struct params does not contain the field "cure.SE_DAE"!');
               elseif ~isfield(paramCure,'test')
                   error('Struct params does not contain the field "cure.test"!');
               elseif ~isfield(paramCure,'gif')
                   error('Struct params does not contain the field "cure.gif"!');
               elseif ~isfield(paramCure,'maxIter')
                   error('Struct params does not contain the field "cure.maxiter"!');
               end           
            end
        end
        
        function parsedParams = parseParamsStruct(params,method)
        %Selects the values of the fields of the params-struct that belong 
        %to the given reduction method and cuts of all other fields
        
            % paramters that are identical for all algorithms
            parsedParams.originalOrder = params.originalOrder;
            
            % algorithm-specific parameters
            if strcmp(method,'tbr')                     %tbr
               parsedParams.type = params.type;
               parsedParams.redErr = params.redErr;
               parsedParams.hsvTol = params.hsvTol;
               parsedParams.warnOrError = params.warnOrError;
               parsedParams.lse = params.lse;
               parsedParams.hsv = params.hsv;                  
            elseif strcmp(method,'modalMor')            %modalMor 
               parsedParams.type = params.type;
               parsedParams.orth = params.orth;
               parsedParams.real = params.real;
               parsedParams.tol = params.tol;
               parsedParams.dominance = params.dominance;
            elseif strcmp(method,'irka')                %irka
               parsedParams.maxiter = params.maxiter;
               parsedParams.tol = params.tol;
               parsedParams.type = params.type;
               parsedParams.verbose = params.verbose;
               parsedParams.stopCrit = params.stopCrit;
               parsedParams.suppressverbose = params.suppressverbose;
               parsedParams.orth = params.orth;
               parsedParams.lse = params.lse;
               parsedParams.dgksTol = params.dgksTol;
               parsedParams.krylov = params.krylov;
               parsedParams.s0 = params.s0;
               parsedParams.Rt = params.Rt;
               parsedParams.Lt = params.Lt;
               parsedParams.kIter = params.kIter;
               parsedParams.s0Traj = params.s0Traj;
               parsedParams.RtTraj = params.RtTraj;
               parsedParams.LtTraj = params.LtTraj;  
            elseif strcmp(method,'rk')                  %rk
               parsedParams.real = params.real;
               parsedParams.orth = params.orth;
               parsedParams.reorth = params.reorth;
               parsedParams.lse = params.lse;
               parsedParams.dgksTol = params.dgksTol;
               parsedParams.krylov = params.krylov;
               parsedParams.IP = params.IP;
               parsedParams.Rt = params.Rt;
               parsedParams.Lt = params.Lt;
               parsedParams.s0_inp = params.s0_inp;
               parsedParams.s0_out = params.s0_out;
            elseif strcmp(method,'projectiveMor')       %projectiveMor
               parsedParams.trans = params.trans;
            elseif strcmp(method,'porkV')               %porkV
            elseif strcmp(method,'porkW')               %porkW
            elseif strcmp(method,'cure')
               
            end
            
        end
        
        function parsedParamsList = removeProjectiveMor(paramsList)
        %Removes the reductionParameters for "projectiveMor" from the
        %reduction history, because "projectiveMor" performs just the 
        %projection, but is not a standalone reduction algorithm
        
            if size(paramsList,1) > 1
                counter = 1;
                parsedParamsList = cell(1,1);
                for i = 1:size(paramsList,1)
                   if ~strcmp(paramsList{i,1}.method,'projectiveMor')
                       parsedParamsList{counter,1} = paramsList{i,1};
                       counter = counter + 1;
                   end
                end
            else
                parsedParamsList = paramsList;
            end           
        end
    end    
end

