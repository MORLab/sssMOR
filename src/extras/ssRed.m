classdef ssRed < ss
% SSRED - Reduced state-space LTI system (ssRed) class
%
% Syntax:
%       sysr = ssRed(method,params,A,B,C)
%       sysr = ssRed(method,params,A,B,C,D)
%       sysr = ssRed(method,params,A,B,C,D,E)
%       sysr = ssRed(method,params,sys_sss)
%       sysr = ssRed(method,params,sys_ss)
%
% Description:
%       This class is derived from the ss-class. It is used to represent
%       models which are the result of reducing the order of large models
%       with specific algorithms. 
%
% Input Arguments:
%       -method: name of the used reduction algorithm;
%                ['tbr' / 'modalMor' / 'irka' / 'rk' / 'projectiveMor' / 'cure']
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
%       -params (rk):           structure with the parameters for the
%                               rk-algorithm
%           -.originalOrder:    Model order before reduction
%           -.real:             keep the projection matrices real
%                               [true / false]
%       -params (projectiveMor):structure with the parameters for the
%                               projectiveMor-algorithm
%           -.originalOrder:    Model order before reduction
%           -.trans:            choose how W should be transposed
%                               [ {T} / H ]
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
        reductionMethod
        reductionParameters
    end
    properties(Dependent, Hidden)
        n,p,m
        isSiso, isSimo, isMiso, isMimo, isBig
        isDae, isDescriptor
    end
    
    methods
        function obj = ssRed(varargin)
            % parse and check the arguments
            if nargin < 3 || nargin > 7 || nargin == 4
                error('Invalid syntax for the "ssRed" command. Type "help ssRed" for more information.');
            elseif nargin == 3
                sys = varargin{3};
                if isa(sys,'ss')
                   A = sys.a;
                   B = sys.b;
                   C = sys.c;
                   D = sys.d;
                   E = sys.e;
                elseif isa(sys,'sss')
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
                
                if nargin == 6
                    D = full(varargin{6});
                    E = [];
                elseif nargin == 7
                    D = full(varargin{6});
                    E = full(varargin{7});
                else                
                    D = [];
                    E = [];
                end             
            end
            
            if ~isa(varargin{1},'char') || ismember(varargin{1},{'tbr','modalMor','rk','irka','projectiveMor','cure'}) == 0
                error('The first argument has a wrong format. Type "help ssRed" for more information.');
            end
            
            % call the construktor of the superclass ss
            obj@ss(A,B,C,D);
            obj.e = E;
            obj.reductionMethod = varargin{1};
            obj.reductionParameters = varargin{2};
            
            % check the params struct
            try
               obj.checkParamsStruct(varargin{2},varargin{1});
            catch ex
               error(ex.message); 
            end
                
        end
        
        function value = get.reductionMethod(obj)
            value = obj.reductionMethod;
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
               end
               
%                if ~(ischar(params.type) && ismember(params.type,{'tbr','adi','matchDcGain'}))
%                    error('params.type: Wrong value. Type "help ssRed" for more information.');
%                elseif ~(ischar(params.redErr) && strcmp(params.redErr,'0')) && ...
%                       ~(isscalar(params.redErr) && params.redErr >= 0)
%                    error('params.redErr: Wrong value. Type "help ssRed" for more information.');
%                elseif ~isscalar(params.hsvTol) || params.hsvTol <= 0
%                    error('params.redErr: Wrong value. Positive float value required!');
%                elseif ~(ischar(params.warnOrError) && ismember(params.warnOrError,{'warn','error','0'}))
%                    error('params.warnOrError: Wrong value. Type "help ssRed" for more information.');
%                elseif ~(ischar(params.lse) && ismember(params.lse,{'gauss','luChol'}))
%                    error('params.lse: Wrong value. Type "help ssRed" for more information.');
%                end               
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
               
%                if ~isscalar(params.type) && ~(ischar(params.type) && ...
%                   ismember(params.type,{'SM','LM','SA','LA','SR','LR'}))
%                    error('params.type: Wrong value. Type "help ssRed" for more information.');
%                elseif ~(ischar(params.orth) && ismember(params.orth,{'0','qr'}))
%                    error('params.orth: Wrong value. Type "help ssRed" for more information.');
%                elseif ~(ischar(params.real) && ismember(params.real,{'real','0'}))
%                    error('params.real: Wrong value. Type "help ssRed" for more information.');
%                elseif ~isscalar(params.tol) || params.tol <= 0
%                    error('params.tol: Wrong value. Positive float value required');
%                elseif ~(ischar(params.dominance) && ismember(params.dominance,{'0','analyze','2q','3q','4q'}))
%                    error('params.dominance: Wrong value. Type "help ssRed" for more information.');
%                end
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
               end
               
%                if ~isscalar(params.maxiter) || mod(params.maxiter,1) ~= 0 || ...
%                   params.maxiter <= 0
%                    error('params.maxiter: Wrong value. Positive integer value required');
%                elseif ~isscalar(params.tol) || params.tol <= 0
%                    error('params.tol: Wrong value. Positive float value required');
%                elseif ~(ischar(params.type) && ismember(params.type,{'','stab'}))
%                    error('params.type: Wrong value. Type "help ssRed" for more information.');
%                elseif params.verbose ~= 0 && params.verbose ~= 1
%                    error('params.verbose: Wrong value. Value 0 or 1 required');
%                elseif ~(ischar(params.stopCrit) && ismember(params.stopCrit,{'combAny','s0','sysr','combAll'}))
%                    error('params.stopCrit: Wrong value. Type "help ssRed" for more information.');
%                elseif params.suppressverbose ~= 0 && params.suppressverbose ~= 1
%                    error('params.suppressverbose: Wrong value. Value 0 or 1 required');
%                end               
            elseif strcmp(method,'rk')                  %rk
               if ~isfield(params,'real')
                   error('Struct params does not contain the field "real"!');
               end
               
%                if params.real ~= 1 && params.real ~= 0
%                    error('params.verbose: Wrong value. Value 0 or 1 required');
%                end
            elseif strcmp(method,'projectiveMor')       %projectiveMor
               if ~isfield(params,'trans')
                   error('Struct params does not contain the field "trans"!');
               end
%                if ~(ischar(params.trans) && ismember(params.trans,{'T','H'})) 
%                    error('params.trans: Wrong value. Type "help ssRed" for more information.');
%                end
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
    end    
end

