classdef ssRed < ss
% SSRED - Reduced state-space LTI system (ssRed) class
%
% Syntax:
%       sysr = ssRed(A,B,C)
%       sysr = ssRed(A,B,C,D)
%       sysr = ssRed(A,B,C,D,E)
%       sysr = ssRed(A,B,C,D,E,method,params)
%       sysr = ssRed(A,B,C,D,E,method,params,sys)
%       sysr = ssRed(A,B,C,D,E,method,params,paramsList)
%       
%
% Description:
%       This class is derived from the ss-class. It is used to represent
%       models which are the result of reducing the order of large models
%       with specific algorithms. ssRed objects have all attributes that ss
%       objects normally possess
%
% Input Arguments:
%       -A: system matrix
%       -B: input matrix
%       -C: output matrix
%       -D: static gain matrix
%       -E: descriptor matrix
%       -method: name of the used reduction algorithm;
%                ['tbr' / 'modalMor' / 'irka' / 'rk' / 'projectiveMor' / 'porkV' / 'porkW' / 'spark' / 'cure_spark' / 'cure_irka' / 'cure_rk+pork' / 'stabsep' / 'rkOp' / 'rkIcop' / 'modelFct' / 'cirka' / 'userDefined']
%       -sys:   original state space system before reduction (class sss or ssRed)
%       -paramsList:    Structure array. Each entry represents one 
%                       reduction (without the current reduction).
%           -.method:           name of the applied reduction algorithm (see above)                               
%           -.params:           structure with the parameters of the
%                               applied algorithm (see below)
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
%       -params (spark):        structure with the parameters for the
%                               spark-algorithm
%           -.originalOrder:    Model order before reduction
%           -.spark.type:       chooses between standard SPARK, where the original 
%                               model is reduced directly, or MESPARK, where a 
%                               model function is created and updated after convergence.
%                               [{'model'} / 'standard']
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
%           -.spark.pork:       projection with porkV (input krylov subspace)
%                               or porkW (output krylov subspace)
%                               [{'V'} / 'W']
%           -.mespark.ritz:     use eigenvalues of model function to initialize
%                               the shifts 
%                               [{'1'} / '0']
%           -.mespark.pertIter: number of iterations after which a
%                               pertubation of the shifts starts to avoid
%                               stagnation of the model function
%                               [{'5'} / positive integer]
%           -.mespark.maxIter:  maximum number of model function updates
%                               [{'20'} / positive integer]
%       -params (cure_spark):   structure with the parameters for the
%                               cure-algorithm (reduction algorithm spark)
%           -.originalOrder:    Model order before reduction
%           -.currentReducedOrder: Model order at the current iteration
%                                  step of the cure algorithm
%           -.shifts:           Shifts used for the current iteration step 
%                               of the cure algorihm
%           -.cure.fact:        factorization mode 
%                               [{'V'} / 'W']
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
%           -.shifts:           Shifts used for the current iteration step 
%                               of the cure algorihm
%           -.cure.fact:        factorization mode 
%                               [{'V'} / 'W']
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
%           -.shifts:           Shifts used for the current iteration step 
%                               of the cure algorihm
%           -.cure.fact:        factorization mode 
%                               [{'V'} / 'W']
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
%       -params (stabsep):      structure with the parameters for the
%                               model order reduction resulting from 
%                               applying the stabsep-function
%           -.originalOrder:    Model order before reduction
%           -.reducedOrder:     Model order after reduction
%       -params (rkOp):         structure with the parameters for the
%                               rkOp-algorithm
%           -.originalOrder:    Model order before reduction
%           -.sOpt:             optimal expansion point
%           -.rk:               reduction type
%                               [{'twoSided'} / 'input' / 'output']
%           -.lse:              solve linear system of equations
%                               [{'sparse'} / 'full' / 'gauss' / 'hess' / 'iterative' ]
%       -params (rkIcop):       structure with the parameters for the
%                               rkIcop-algorithm
%           -.originalOrder:    Model order before reduction
%           -.sOpt:             optimal expansion point
%           -.s0:               inital expansion point (scalar)
%           -.rk:               reduction type
%                               [{'twoSided'} / 'input' / 'output']
%           -.maxIter:          maximum number of iterations;
%                               [{20} / positive integer]
%           -.tol:              convergence tolerance;
%                               [{1e-2} / positive float]
%           -.lse:              solve linear system of equations
%                               [{'sparse'} / 'full' / 'gauss' / 'hess' / 'iterative' ]
%       -params (modelFct):     structure with the parameters for the
%                               modelFct-algorithm
%           -.originalOrder:    Model order before reduction
%           -.s0mTot:           vector of all shifts used
%           -.updateModel:      chooses which shifts of s0mr are included in update;
%                               [{'new'} / 'all' / 'lean' ]
%           -.modelTol:         tolerance for identifying new shifts;
%                               [{1e-2} / positive float ]
%       -params (cirka):        structure with the parameters for the
%                               cirka-algorithm
%           -.originalOrder:    Model order before reduction
%           -.modelFctOrder:    Final order of the model function
%           -.kIrka:            Vector of irka iterations
%           -.s0:               final shifts for reduction
%           -.qm0:              initial size of model function;
%                               [{2*length(s0)} / positive integer]
%           -.s0m:              initial shifts for surrogate;
%                               [{[s0,s0]} / vector ]
%           -.maxiter:          maximum number of iterations;
%                               [{15} / positive integer]
%           -.tol:              convergence tolerance;
%                               [{1e-3} / positive float]
%           -.stopCrit:         convergence criterion for CIRKA;
%                               ['s0' / 'sysr' / 'sysm' / {'combAny'} / 'combAll']
%           -.updateModel:      type of model function update;
%                               [{'new'},'all']
%           -.clearInit:        reset the model function after first iteration;
%                               [{true}, false]
%           -.irka:             irka options (cmp irka)
%       -params (userDefined):  []
%
% Output Arguments:
%       -sys: reduced state-space (ssRed)-object
%
% Examples:
%       This code creates an instance of the ssRed-class
%
%> A = rand(10);
%> B = rand(10,1);
%> C = rand(1,10);
%> params.originalOrder = 20;
%> sysr = ssRed(A,B,C,'porkV',params);
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
% Authors:      Niklas Kochdumper, Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  20 Jan 2017
% Copyright (c) 2016,2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------
    
    properties
        x0
        redParam        
        issymmetric
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
        residues
        
        reductionParameters %deprecated; backward compatibility; use redParam
        isSym       %deprecated; backward compatibility; use issymmetric
    end
    properties(Hidden,Access = private)
        a_,b_,c_,d_,e_
    end
    
    methods
        function obj = ssRed(varargin)
            
            % parse input arguments
            if nargin~=1 && ~isempty(varargin{1})   % not an empty model            
                if nargin < 3 || nargin > 8
                    error('Invalid syntax for the "ssRed" command. Type "help ssRed" for more information.');
                end
                
                % default values
                D = []; E = []; method = 'userDefined'; paramsList = [];
                params = []; name = [];
                
                % system matrices
                A = full(varargin{1});
                B = full(varargin{2});
                C = full(varargin{3});
                iTemp = 4;          % next index in varargin
                if nargin >= iTemp && isnumeric(varargin{4})
                   D = full(varargin{4});
                   iTemp = 5;
                   if nargin >= iTemp && isnumeric(varargin{5})
                      E = full(varargin{5});
                      iTemp = 6;
                   end
                end
                
                % additional input arguments
                switch nargin
                    case iTemp
                        error('Invalid syntax for the "ssRed" command. Type "help ssRed" for more information.');
                    case iTemp+1
                        method = varargin{iTemp};
                        params = varargin{iTemp+1};
                    case iTemp+2
                        method = varargin{iTemp};
                        params = varargin{iTemp+1};
                        if isa(varargin{iTemp+2},'sss') || isa(varargin{iTemp+2},'ss')
                            name = varargin{iTemp+2}.Name;
                        end
                        if isa(varargin{iTemp+2},'ssRed')
                            paramsList = varargin{iTemp+2}.reductionParameters;
                        elseif isstruct(varargin{iTemp+2})
                            paramsList = varargin{iTemp+2};
                        end
                end
            else
               % ensures that syntax "ssRed([])" gives back an empty model
               A=[];B=[];C=[];D=[];E=[];name=[];
            end
            
            % call the constructor of the superclass ss
            obj@ss(A,B,C,D,'e',E);
            
            % set the name property of the model
            obj.Name = name;
            
            % store the correct names for the properties containing the
            % system matrices. This is a compatibility fix because the
            % property names changed from lower to upper case in Matlab
            % R2016a. To access the system matrices, use "obj.(obj.a_)"
            obj.a_ = ltipack.matchProperty('a',...
                            ltipack.allprops(obj),class(obj));
            obj.b_ = ltipack.matchProperty('b',...
                            ltipack.allprops(obj),class(obj));
            obj.c_ = ltipack.matchProperty('c',...
                            ltipack.allprops(obj),class(obj));
            obj.d_ = ltipack.matchProperty('d',...
                            ltipack.allprops(obj),class(obj));
            obj.e_ = ltipack.matchProperty('e',...
                            ltipack.allprops(obj),class(obj));
            
            % verify correctness of input parameters
            if nargin~=1 && ~isempty(varargin{1})   % not an empty model
                            
                % parameter "method"
                if ~isa(method,'char') 
                    error('Argument "method" has to be a string. Type "help ssRed" for more information.');
                end
                
                % check all fields of paramsList
                obj.checkParamsList(paramsList)

                % update the reductionParameters list
                if isempty(paramsList)
                   paramsTemp(1).method = method;
                   try
                       paramsTemp(1).params = obj.parseParamsStruct(params,method,1);
                   catch ex
                       error('Argument "params" has the wrong format. Type "help ssRed" for more information.');
                   end
                   obj.reductionParameters = paramsTemp;
                else
                   len = size(paramsList,2);
                   paramsList(len+1).method = method;
                   try
                       paramsList(len+1).params = obj.parseParamsStruct(params,method,1);
                   catch ex
                       error('Argument "params" has the wrong format. Type "help ssRed" for more information.');
                   end
                   obj.reductionParameters = paramsList;
                end
                
                % postprocess cure parameters
                obj.reductionParameters = obj.removeCureParameters(obj.reductionParameters);
                
                % set initial values
                obj.x0 = [];
            end
        end
                
        %% Get Basic Properties
        function x0 = get.x0(sys)
            x0 = sys.x0;
            if isempty(x0)
                x0 = zeros(sys.n,1);
            end
        end
        
        function m = get.m(sys) % number of inputs
            m = size(sys.(sys.b_),2);
        end
        function n = get.n(sys) % system order
            n = size(sys.(sys.a_),1);
        end
        function p = get.p(sys) % number of outputs
            p = size(sys.(sys.c_),1);
        end
        
        %% Set Basic Properties
        function sys = set.x0(sys, x0)
            if (~isempty(x0)) && (any(size(x0) ~= [sys.n,1]))
                error('A and x0 must have the same number of rows.')
            end
            sys.x0 = x0;
        end
        
        function sys = set.redParam(sys,redParam)
            try
                sys.checkParamsList(redParam)
            catch ex
                error('Invalid value for the property "reductionParameters". Type "help ssRed" for more information.');
            end
            sys.redParam = redParam;            
        end
        
        %% Get helper functions
        function isSiso = get.isSiso(sys); isSiso=(sys.p==1)&&(sys.m==1); end
        function isSimo = get.isSimo(sys); isSimo=(sys.p>1)&&(sys.m==1); end
        function isMiso = get.isMiso(sys); isMiso=(sys.p==1)&&(sys.m>1); end
        function isMimo = get.isMimo(sys); isMimo=(sys.p>1)||(sys.m>1); end
        function isBig = get.isBig(sys); isBig=(sys.n>5000);end
        
        function isDae = get.isDae(sys)
            if condest(sys.(sys.e_))==Inf
                isDae = 1;
            else
                isDae = 0;
            end
        end
        
        function isDescriptor = get.isDescriptor(sys)
            isDescriptor = logical(full(any(any(sys.(sys.e_)-eye(size(sys.(sys.e_)))))));
        end
        
        function sys = resolveDescriptor(sys)
            sys.(sys.a_) = sys.(sys.e_)\sys.(sys.a_);
            sys.(sys.b_) = sys.(sys.e_)\sys.(sys.b_);
            sys.(sys.e_) = eye(size(sys.(sys.a_))); 
            %makes the usage of sys.e in computations more robust than []
        end
        
        function reductionParameters = get.reductionParameters(sys)
           reductionParameters = sys.redParam;
        end
        
        function sys = set.reductionParameters(sys,redParam)
           sys.redParam = redParam; 
        end
        
        function isSym = get.isSym(sys)
            isSym = sys.issymmetric;
        end   
        
        function issymmetric = get.issymmetric(sys) %A=A', E=E'
            if isequal(sys.issymmetric,0) || isequal(sys.issymmetric,1)
                issymmetric = sys.issymmetric;
            else
                if full(max(max(sys.(sys.a_)-sys.(sys.a_).')))<1e-6 && full(max(max(sys.(sys.e_)-sys.(sys.e_).')))<1e-6
                    issymmetric = true;
                else
                    issymmetric = false;
                end
            end
        end
        
        %% Overload subsref to cope with compatibility issues in MATLAB
        function result = subsref(sys, arg)
            %   Parts are taken from built-in subsref
            %   Copyright 1986-2010 The MathWorks, Inc.
            
            ni = nargin;
            if ni==1,
               result = M;  return
            end
            switch arg(1).type
              case '.'
                  %COMPATIBILITY ISSUE BETWEEN 2015b and 2016a
                  % make sure system matrices have the right (lower/upper)
                  % case
                  if any(strcmp(arg(1).subs,{'a','b','c','d','e','A','B','C','D','E'}))
                        arg(1).subs = ltipack.matchProperty(arg(1).subs,...
                            ltipack.allprops(sys),class(sys));
                  end
                 result = builtin('subsref',sys,arg(1));
              case '()'
                 result = subparen(sys,arg(1).subs);
              case '{}'
                 ctrlMsgUtils.error('Control:ltiobject:subsref3')
            end
           if length(arg)>1,
              % SUBSREF for InputOutputModel objects can't be invoked again
              % inside this method so jump out to make sure that downstream
              % references to InputOutputModel properties are handled correctly,
              result = ltipack.dotref(result,arg(2:end));
           end
    end
        %% Override operators and build-in-functions
        
        function infostr = disp(sys)
        % Displays information about a reduced state-space model (Similar
        % to sss/disp, but with additional information)
            if isempty(sys)
                fprintf(1,'  Empty reduced state-space model.\n\n');
            else
                mc = metaclass(sys);
                str = [];
                if ~isempty(mc.Name) && ~isempty(sys.Name)
                    str = [mc.Name ' Model ' sys.Name, ' '];
                end

                if sys.isDae;            str = [str '(DAE)'];
                elseif sys.isDescriptor; str = [str '(DssRed)'];
                else                     str = [str '(ssRed)'];
                end

                if sys.isSiso;       str = [str '(SISO)'];
                elseif sys.isSimo;   str = [str '(SIMO)'];
                elseif sys.isMiso;   str = [str '(MISO)'];
                elseif sys.isMimo;   str = [str '(MIMO)'];
                end

                str = [str  char(10), num2str(sys.n) ' states, ' num2str(sys.m) ...
                    ' inputs, ' num2str(sys.p) ' outputs'];

                if sys.Ts==0
                    str = [str  char(10) 'Continuous-time state-space model.'];
                else
                    str = [str  char(10) 'Sample time: ' num2str(sys.Ts) ' seconds'];
                    str = [str  char(10) 'Discrete-time state-space model.'];
                end
                
                % display all reduction algorithms used, but each reduction
                % algorithm only once, even if it was used multiple times
                usedMethods = [];
                methodString = '';
                for i = 1:size(sys.reductionParameters,2)
                   if ~ismember(sys.reductionParameters(1,i).method,usedMethods)
                       usedMethods{end+1} = sys.reductionParameters(1,i).method;
                       methodString = strcat(methodString,sys.reductionParameters(1,i).method,',',char(1));
                   end
                end
                methodString = methodString(1:size(methodString,2)-2);
                
                if strcmp(sys.reductionParameters(end,1).method,'userDefined')
                    str = [str char(10) 'Reduction Method(s): ' methodString char(10)];
                else
                    str = [str char(10) 'Reduction Method(s): ' methodString char(10) ...
                            'Original order: ' num2str(sys.reductionParameters(1,1).params.originalOrder)];
                end

                if nargout>0
                    infostr = {str};
                else
                    str = strrep(str, char(10), [char(10) '  ']);
                    disp(['  ' str char(10)]);
                end
            end
        end
        
        function sys = clear(sys)
            sys = ssRed([]);
        end
        
        function varargout = diag(varargin)
            [varargout{1:nargout}] = sssFunc.diag(varargin{:});
        end
        
        function varargout = eig(varargin)
            [varargout{1:nargout}] = sssFunc.eig(varargin{:});
        end
        
        function varargout = minus(varargin) 
            [varargout{1:nargout}] = sssFunc.minus(varargin{:});
        end
        
        function varargout = plus(varargin) 
            [varargout{1:nargout}] = sssFunc.plus(varargin{:});
        end
        
        function varargout = mtimes(varargin) 
            [varargout{1:nargout}] = sssFunc.mtimes(varargin{:});
        end
        
        function varargout = residue(varargin)  
            [varargout{1:nargout}] = sssFunc.residue(varargin{:});
        end
        
        function varargout = spy(varargin)  
            [varargout{1:nargout}] = sssFunc.spy(varargin{:});
        end
        
        function varargout = sss(varargin)
            error(['In sssMOR convertion from class "ssRed" to class ', ...
                   '"sss" is prohibited because this would result in a loss ', ...
                   'of information. If conversion is nevertheless desired, ', ...
                   'use syntax "sysSss = sss(sys.A,sys.B,sys.C,sys.D,sys.E)."']);
        end
        
        function varargout = decayTime(varargin)
            [varargout{1:nargout}] = sssFunc.decayTime(varargin{:});
        end
        
        function varargout = issd(varargin)
            [varargout{1:nargout}] = sssFunc.issd(varargin{:});
        end
        
        function varargout = eigs(varargin)
            [varargout{1:nargout}] = sssFunc.eigs(varargin{:});
        end
        
        function syst = truncate(sys, idxOut, idxIn)
            args.type = '()';
            args.subs{1} = idxOut;
            args.subs{2} = idxIn;
            syst = subsref(sys,args);
        end
        
        function varargout = freqresp(varargin)
            % check if Options are specified
            if ~isempty(varargin) && isstruct(varargin{end})
                Opts = varargin{end};
                varargin = varargin(1:end-1);
            end
            % call the correct build in functions
            if exist('Opts','var')
                if isfield(Opts,'maxPoints')
                    warning('Value for option "maxPoints" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'lse')
                    warning('Value for option "lse" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'frd') && Opts.frd == 1 && nargout == 1
                    [G,w] = freqresp@ss(varargin{:});
                    varargout{1} = frd(G,w);
                else
                    [varargout{1:nargout}] = freqresp@ss(varargin{:});
                end
            else
                [varargout{1:nargout}] = freqresp@ss(varargin{:});
            end
        end
        
        function  varargout = impulse(varargin)
            % check if Options are specified
            if ~isempty(varargin) && isstruct(varargin{end})
                Opts = varargin{end};
                varargin = varargin(1:end-1);
            end
            % call the correct function depending on the value of Opts.tf
            if exist('Opts','var')
                if isfield(Opts,'odeset')
                    warning('Value for option "odeset" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'tolOutput')
                    warning('Value for option "tolOutput" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'tolState')
                    warning('Value for option "tolState" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'ode')
                    warning('Value for option "ode" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'tsMin')
                    warning('Value for option "tsMin" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'tLin')
                    warning('Value for option "tLin" remains ineffective for ssRed-objects'); 
                end
                
                if isfield(Opts,'tf') && Opts.tf == 1 
                    if nargout > 3
                       error('Maximal three output arguments if Opts.tf == 1'); 
                    end
                    varargout{1} = tf(varargin{1});
                    if nargout > 1
                       [varargout{2},varargout{3}] = impulse@ss(varargin{:}); 
                    end
                else
                    [varargout{1:nargout}] = impulse@ss(varargin{:});
                end
            else
                [varargout{1:nargout}] = impulse@ss(varargin{:});
            end
        end
        
        function  varargout = step(varargin)
            % check if Options are specified
            if ~isempty(varargin) && isstruct(varargin{end})
                Opts = varargin{end};
                varargin = varargin(1:end-1);
            end
            % call the correct function depending on the value of Opts.tf
            if exist('Opts','var')
                if isfield(Opts,'odeset')
                    warning('Value for option "odeset" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'tolOutput')
                    warning('Value for option "tolOutput" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'tolState')
                    warning('Value for option "tolState" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'ode')
                    warning('Value for option "ode" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'tsMin')
                    warning('Value for option "tsMin" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'htCell')
                    warning('Value for option "htCell" remains ineffective for ssRed-objects'); 
                end
                if isfield(Opts,'tLin')
                    warning('Value for option "tLin" remains ineffective for ssRed-objects'); 
                end
                
                if isfield(Opts,'tf') && Opts.tf == 1 
                    if nargout > 3
                       error('Maximal three output arguments if Opts.tf == 1'); 
                    end
                    varargout{1} = tf(varargin{1});
                    if nargout > 1
                       [varargout{2},varargout{3}] = step@ss(varargin{:}); 
                    end
                else
                    [varargout{1:nargout}] = step@ss(varargin{:});
                end
            else
                [varargout{1:nargout}] = step@ss(varargin{:});
            end
        end
        
        function  varargout = lyapchol(varargin)
            [varargout{1:nargout}] = sssFunc.lyapchol(varargin{:});
        end
        
        function varargout = stabsep(varargin)
            [varargout{1:nargout}] = stabsep@ss(varargin{:});
            % add an entry to the reduction history if the model order was
            % changed
            if nargout >= 1
                if varargout{1}.n < varargin{1}.n
                    sys = varargout{1};
                    params.originalOrder = varargin{1}.n;
                    params.reducedOrder = sys.n;
                    varargout{1} = ssRed(sys.(sys.a_),sys.(sys.b_), ...
                                         sys.(sys.c_),sys.(sys.d_), ...
                                         sys.(sys.e_),'stabsep', ...
                                         params, varargin{1});
                    
                    % make robust to computations with sys.e_
                    if isempty(varargout{1}.(sys.e_))
                        varargout{1}.(sys.e_) = eye(sys.n);
                    end
                        
                    if nargout >= 2
                        sys = varargout{2};
                        params.originalOrder = varargin{1}.n;
                        params.reducedOrder = sys.n;
                        varargout{2} = ssRed(sys.(sys.a_),sys.(sys.b_), ...
                                             sys.(sys.c_),sys.(sys.d_), ...
                                             sys.(sys.e_),'stabsep', ...
                                             params, varargin{1});
                    end
                    % make robust to computations with sys.e_
                    if isempty(varargout{1}.(sys.e_))
                        varargout{1}.(sys.e_) = eye(sys.n);
                    end
                end
            end
        end
    end
    
    %%Private and static helper methods
    methods(Hidden, Access = private, Static)
        
        function checkParamsList(paramsList)
        % Checks if the list "paramsList" with reduction parameters 
        % resulting from multiple reduction steps has the correct format.
        % If this is not the case, an error is produced
        
           if ~isempty(paramsList)
                try
                   for i = 1:size(paramsList,1)
                      if length(fieldnames(paramsList(i))) ~= 2
                          error('Argument "paramsList" has the wrong format. Type "help ssRed" for more information.');
                      end
                      if i > 1 && ismember(paramsList(i-1).method,{'cure_spark','cure_irka','cure_rk+pork'})
                          ssRed.parseParamsStruct(paramsList(i).params,paramsList(i).method,0);
                      else
                          ssRed.parseParamsStruct(paramsList(i).params,paramsList(i).method,1);
                      end
                   end
                catch ex
                   error('Argument "paramsList" has the wrong format. Type "help ssRed" for more information.'); 
                end
            end 
        end
        
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
        
            % check if the reduction method belongs to the predefined 
            % method
            if ~ismember(method,{'tbr', ...
                    'modalMor','rk','irka','projectiveMor','porkV','porkW', ...
                    'spark','cure_spark','cure_irka','cure_rk+pork', ...
                    'stabsep','rkOp','rkIcop','modelFct','cirka', ...
                    'userDefined'})
                % for user defined custom reduction methods, the struct can
                % be arbitrary. Therefore there are no checks neccesary
                parsedStruct = params;
            else
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
                elseif strcmp(method,'modelFct')            %modelFct
                   list = {'originalOrder','s0mTot','updateModel','modelTol'};
                   parsedStruct = ssRed.parseStructFields(params,list,'params');
                elseif strcmp(method,'cirka')               %cirka
                   list = {'originalOrder','modelFctOrder','kIrka','s0',...
                            'qm0','s0m','maxiter','tol','stopCrit','updateModel',...
                            'clearInit','irka'};
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
                elseif strcmp(method,'spark')               %spark
                   list = {'originalOrder'};
                   parsedStruct = ssRed.parseStructFields(params,list,'params');

                   list = {'spark','mespark'};
                   ssRed.parseStructFields(params,list,'params');

                   list = {'type','mfe','mi','xTol', ...
                           'pork','fTol','modelTol'};
                   parsedStruct.spark = ssRed.parseStructFields(params.spark,list,'params.spark');

                   list = {'ritz','pertIter','maxIter'};
                   parsedStruct.mespark = ssRed.parseStructFields(params.mespark,list,'params.mespark');
                elseif strcmp(method,'cure_spark')          %cure_spark  
                   list = {'originalOrder','currentReducedOrder','shifts'};
                   parsedStruct = ssRed.parseStructFields(params,list,'params'); 

                   if cureRequired
                        list = {'cure','spark','mespark'};
                        ssRed.parseStructFields(params,list,'params');

                        list = {'fact','stop','stopval','maxIter'};
                        parsedStruct.cure = ssRed.parseStructFields(params.cure,list,'params.cure');
                   else
                        list = {'spark','mespark'};
                        ssRed.parseStructFields(params,list,'params');
                   end

                   list = {'type','mfe','mi','xTol', ...
                           'fTol','modelTol'};
                   parsedStruct.spark = ssRed.parseStructFields(params.spark,list,'params.spark');

                   list = {'ritz','pertIter','maxIter'};
                   parsedStruct.mespark = ssRed.parseStructFields(params.mespark,list,'params.mespark');
                elseif strcmp(method,'cure_irka')           %cure_irka
                   list = {'originalOrder','currentReducedOrder','shifts','maxiter','tol', ...
                           'type','stopCrit','orth','lse','dgksTol','krylov', ...
                           's0','Rt','Lt','kIter','s0Traj','RtTraj','LtTraj'};
                   parsedStruct = ssRed.parseStructFields(params,list,'params');

                   if cureRequired
                        list = {'cure'};
                        ssRed.parseStructFields(params,list,'params');

                        list = {'fact','stop','stopval','maxIter'};
                        parsedStruct.cure = ssRed.parseStructFields(params.cure,list,'params.cure');  
                   end
                elseif strcmp(method,'cure_rk+pork')        %cure_rk+pork
                   list = {'originalOrder','currentReducedOrder','shifts','real', ...
                           'orth','reorth','lse','dgksTol','krylov', ...
                           'IP','Rt','Lt','s0_inp','s0_out'};
                   parsedStruct = ssRed.parseStructFields(params,list,'params');

                   if cureRequired
                        list = {'cure'};
                        ssRed.parseStructFields(params,list,'params');

                        list = {'fact','stop','stopval','maxIter'};
                        parsedStruct.cure = ssRed.parseStructFields(params.cure,list,'params.cure');
                   end
                elseif strcmp(method,'stabsep')             %stabsep
                   list = {'originalOrder','reducedOrder'};
                   parsedStruct = ssRed.parseStructFields(params,list,'params');
                elseif strcmp(method,'rkOp')                %rkOp
                   list = {'originalOrder','sOpt','rk','lse'};
                   parsedStruct = ssRed.parseStructFields(params,list,'params');  
                elseif strcmp(method,'rkIcop')              %rkIcop
                   list = {'originalOrder','sOpt','s0','rk','maxIter','tol','lse'};
                   parsedStruct = ssRed.parseStructFields(params,list,'params'); 
                elseif strcmp(method,'userDefined')         %userDefined
                   parsedStruct = params;
                end
            end
        end
        
        function parsedParamsList = removeCureParameters(paramsList)
        %The parameters under the field "cure" stay the same for all cure
        %iterations. Because of this, the field "cure" is only stored for
        %the first iteration step. Therefore, this function removes all the
        %redundant cure paramters from the reduction history
            
            parsedParamsList = paramsList;
        
            for i = 2:length(paramsList)
                if ismember(paramsList(i-1).method,{'cure_spark','cure_irka','cure_rk+pork'}) && ...
                   isfield(paramsList(i).params,'cure')
                    parsedParamsList(i).params = rmfield(parsedParamsList(i).params,'cure');
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
        end
        
    end    
end

