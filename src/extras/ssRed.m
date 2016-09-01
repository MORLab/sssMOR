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
%           -.stopCrit:	        stopping criterion;
%                               ['combAny' / 's0' / 'sysr' / 'combAll']
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
%       -params (cure_spark):   structure with the parameters for the
%                               cure-algorithm (reduction algorithm spark)
%           -.originalOrder:    Model order before reduction
%           -.currentReducedOrder: Model order at the current iteration
%                                  step of the cure algorithm
%           -.cure.fact:        factorization mode 
%                               [{'V'} / 'W']
%           -.cure.init:        shift initialization mode 
%                               [{'sm'} / 'zero' / 'lm' / 'slm']
%           -.cure.stop:        stopping criterion
%                               [{'nmax'} / 'h2Error']
%           -.cure.stopval:     value according to which the stopping criterion is evaluated
%                               [{'round(sqrt(sys.n))'} / positive integer]
%           -.cure.maxIter:     maximum number of CURE iterations
%                               [{'20'} / positive integer]
%           -.spark.type:       chooses between standard SPARK, where the original 
%                               model is reduced directly, or MESPARK, where a 
%                               model function is created and updated after convergence.
%                               [{'model'} / 'standard']
%           -.spark.test:       specifies weather the user desires to get insight 
%                               in what is happening. This is realized by  
%                               plotting intermediate results during optimization.
%                               [{'0'} / '1']
%           -.spark.mfe:        maximum functions evaluations
%                               [{'5e3'} / positive integer]
%           -.spark.mi:         maximum iterations in solver
%                               [{'150'} / positive integer]
%           -.spark.xTol:       step tolerance in solver
%                               [{'1e-10'} / positive float]
%           -.spark.fTol:       function value tolerance
%                               [{'1e-10'} / positive float]
%           -.spark.modelTol:   convergence tolerance for model funciton
%                               [{'1e-5'} / positive float]
%           -.mespark.ritz:     use eigenvalues of model function to initialize
%                               the shifts 
%                               [{'1'} / '0']
%           -.mespark.pertIter: number of iterations after which a
%                               pertubation of the shifts starts to avoid
%                               stagnation of the model function
%                               [{'5'} / positive integer]
%           -.mespark.maxIter:  maximum number of model function updates
%                               [{'20'} / positive integer]
%       -params (cure_irka):    structure with the parameters for the
%                               cure-algorithm (reduction algorithm irka)
%           -.originalOrder:    Model order before reduction
%           -.currentReducedOrder: Model order at the current iteration
%                                  step of the cure algorithm
%           -.cure.fact:        factorization mode 
%                               [{'V'} / 'W']
%           -.cure.init:        shift initialization mode 
%                               [{'sm'} / 'zero' / 'lm' / 'slm']
%           -.cure.stop:        stopping criterion
%                               [{'nmax'} / 'h2Error']
%           -.cure.stopval:     value according to which the stopping criterion is evaluated
%                               [{'round(sqrt(sys.n))'} / positive integer]
%           -.cure.SE_DAE:      reduction of index 1 semiexplicit DAE 
%                               [{'0'} / '1']
%           -.cure.maxIter:     maximum number of CURE iterations
%                               [{'20'} / positive integer]
%           -.maxiter:          maximum number of iterations;
%                               [50 / positive integer]
%           -.tol:              convergence tolerance;
%                               [1e-3 / positive float]
%           -.type:             choose between different irka modifications;
%                               ['' / 'stab']
%           -.verbose:          show text output during iterations;
%                               [0 / 1]
%           -.stopCrit:	        stopping criterion;
%                               ['combAny' / 's0' / 'sysr' / 'combAll']
%           -.suppressverbose:  suppress any type of verbose for speedup;
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
%       -params (cure_rk+pork): structure with the parameters for the
%                               cure-algorithm (reduction algorithm rk+pork)
%           -.originalOrder:    Model order before reduction
%           -.currentReducedOrder: Model order at the current iteration
%                                  step of the cure algorithm
%           -.cure.fact:        factorization mode 
%                               [{'V'} / 'W']
%           -.cure.init:        shift initialization mode 
%                               [{'sm'} / 'zero' / 'lm' / 'slm']
%           -.cure.stop:        stopping criterion
%                               [{'nmax'} / 'h2Error']
%           -.cure.stopval:     value according to which the stopping criterion is evaluated
%                               [{'round(sqrt(sys.n))'} / positive integer]
%           -.cure.SE_DAE:      reduction of index 1 semiexplicit DAE 
%                               [{'0'} / '1']
%           -.cure.maxIter:     maximum number of CURE iterations
%                               [{'20'} / positive integer]
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
            
            if ~isa(varargin{1},'char') || ismember(varargin{1},{'tbr', ...
                    'modalMor','rk','irka','projectiveMor','porkV','porkW', ...
                    'cure_spark','cure_irka','cure_rk+pork'}) == 0
                error('The first argument has a wrong format. Type "help ssRed" for more information.');
            end
            
            % call the construktor of the superclass ss
            obj@ss(A,B,C,D);
            obj.e = E;
            
            % check all fields of paramsList
            if ~isempty(paramsList)
                try
                   for i = 1:size(paramsList,1)
                      if length(fieldnames(paramsList{i})) ~= 2
                          error('The argument "paramsList" has the wrong format. Type "help ssRed" for more information.');
                      end
                      if i > 1 && ismember(paramsList{i-1}.method,{'cure_spark','cure_irka','cure_rk+pork'})
                          paramsList{i}.params = obj.parseParamsStruct(paramsList{i}.params,paramsList{i}.method,0);
                      else
                          paramsList{i}.params = obj.parseParamsStruct(paramsList{i}.params,paramsList{i}.method,1);
                      end
                   end
                catch ex
                   error('The argument "paramsList" has the wrong format. Type "help ssRed" for more information.'); 
                end
            end
            
            % update the reductionParameters list
            if isempty(paramsList)
               obj.reductionParameters = cell(1,1);
               obj.reductionParameters{1,1}.method = varargin{1};
               obj.reductionParameters{1,1}.params = obj.parseParamsStruct(varargin{2},varargin{1},1);
            else
               len = size(paramsList,1);
               obj.reductionParameters = paramsList;
               obj.reductionParameters{len+1,1}.method = varargin{1};
               obj.reductionParameters{len+1,1}.params = obj.parseParamsStruct(varargin{2},varargin{1},1);
            end
            
            % remove "projectiveMor" from the reduction history
            obj.reductionParameters = obj.removeProjectiveMor(obj.reductionParameters);
            obj.reductionParameters = obj.removeCureParameters(obj.reductionParameters); 
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
        
        function parsedStruct = parseParamsStruct(params,method,cureRequired)
        % Checks if the struct with the parameters  "params" has the correct
        % structure for the specified reduction method "method". If this is
        % the case, then the values of the required fields of the struct 
        % "params" are copied to the struct "parsedStruct". If ths struct
        % "params" contains additional fields, these fields are not copied
        % to the struct "parsedStruct". Because the field "cure" is only
        % stored for the first cure iteration, the parameter "cureRequired"
        % indicates whether the field "cure" is required in the struct
        % "params" or not
            
            % check algorithm-specific parameters
            if strcmp(method,'tbr')                     %tbr
               list = {'originalOrder','type','redErr','hsvTol','lse','hsv'};
               parsedStruct = ssRed.parseStructFields(params,list,'params');                 
            elseif strcmp(method,'modalMor')            %modalMor  
               list = {'originalOrder','type','orth','real','tol','dominance'};
               parsedStruct = ssRed.parseStructFields(params,list,'params'); 
            elseif strcmp(method,'irka')                %irka
               list = {'originalOrder','maxiter','tol','type','stopCrit', ...
                       'orth','lse','dgksTol','krylov', ...
                       's0','Rt','Lt','kIter','s0Traj','RtTraj','LtTraj'};
               parsedStruct = ssRed.parseStructFields(params,list,'params');   
            elseif strcmp(method,'rk')                  %rk
               list = {'originalOrder','real','orth','reorth','lse','dgksTol','krylov', ...
                       'IP','Rt','Lt','s0_inp','s0_out'};
               parsedStruct = ssRed.parseStructFields(params,list,'params');
            elseif strcmp(method,'projectiveMor')       %projectiveMor
               list = {'originalOrder','trans'};
               parsedStruct = ssRed.parseStructFields(params,list,'params');
            elseif strcmp(method,'porkV')               %porkV
               list = {'originalOrder'};
               parsedStruct = ssRed.parseStructFields(params,list,'params');
            elseif strcmp(method,'porkW')               %porkW
               list = {'originalOrder'};
               parsedStruct = ssRed.parseStructFields(params,list,'params');
            elseif strcmp(method,'cure_spark')          %cure_spark  
               list = {'originalOrder','currentReducedOrder'};
               parsedStruct = ssRed.parseStructFields(params,list,'params'); 
                
               if cureRequired
                    list = {'cure','spark','mespark'};
                    ssRed.parseStructFields(params,list,'params');
               
                    list = {'fact','init','stop','stopval','maxIter'};
                    parsedStruct.cure = ssRed.parseStructFields(params.cure,list,'params.cure');
               else
                    list = {'spark','mespark'};
                    ssRed.parseStructFields(params,list,'params');
               end

               list = {'type','test','mfe','mi','xTol', ...
                       'fTol','modelTol'};
               parsedStruct.spark = ssRed.parseStructFields(params.spark,list,'params.spark');
               
               list = {'ritz','pertIter','maxIter'};
               parsedStruct.mespark = ssRed.parseStructFields(params.mespark,list,'params.mespark');
            elseif strcmp(method,'cure_irka')           %cure_irka
               list = {'originalOrder','currentReducedOrder','maxiter','tol', ...
                       'type','stopCrit','orth','lse','dgksTol','krylov', ...
                       's0','Rt','Lt','kIter','s0Traj','RtTraj','LtTraj'};
               parsedStruct = ssRed.parseStructFields(params,list,'params');
               
               if cureRequired
                    list = {'cure'};
                    ssRed.parseStructFields(params,list,'params');
               
                    list = {'fact','init','stop','stopval','maxIter'};
                    parsedStruct.cure = ssRed.parseStructFields(params.cure,list,'params.cure');  
               end
            elseif strcmp(method,'cure_rk+pork')        %cure_rk+pork
               list = {'originalOrder','currentReducedOrder','real', ...
                       'orth','reorth','lse','dgksTol','krylov', ...
                       'IP','Rt','Lt','s0_inp','s0_out'};
               parsedStruct = ssRed.parseStructFields(params,list,'params');
               
               if cureRequired
                    list = {'cure'};
                    ssRed.parseStructFields(params,list,'params');
               
                    list = {'fact','init','stop','stopval','maxIter'};
                    parsedStruct.cure = ssRed.parseStructFields(params.cure,list,'params.cure');
               end
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
        
        function parsedParamsList = removeCureParameters(paramsList)
        %The parameters under the field "cure" stay the same for all cure
        %iterations. Because of this, the field "cure" is only stored for
        %the first iteration step. Therefore, this function removes all the
        %redundant cure paramters from the reduction history
            
            parsedParamsList = paramsList;
        
            for i = 2:length(paramsList)
                if ismember(paramsList{i-1}.method,{'cure_spark','cure_irka','cure_rk+pork'}) && ...
                   isfield(paramsList{i}.params,'cure')
                    parsedParamsList{i}.params = rmfield(parsedParamsList{i}.params,'cure');
                end
            end
        end
        
        function outputStruct = parseStructFields(structure,fields,structName)
        %This function checks if the the struct "structure" contains all fields
        %specified in the cell-array "fields". If the case, all the values 
        %of the fields of the struct "structure" are copied to the struct 
        %"outputStruct". The variable "structName" is used in the error
        %message to inform the user which field caused the error.
             
            for i = 1:length(fields)
               if ~isfield(structure,fields{i})
                  error(strcat('Struct "',structName, ...
                      '" does not contain the field "',fields{i},'"!'));
               else
                  outputStruct.(fields{i}) = structure.(fields{i});
               end
            end            
            t = 1;        
        end
        
    end    
end

