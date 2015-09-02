classdef sss
    % sparse state space model of LTI dynamical system
    %   Detailed explanation goes here
    
    properties
        A,B,C,D,E
    end
    properties(Dependent)
    end
    properties(Dependent, Hidden)
        a,b,c,d,e       
        n,p,m
        is_mimo
    end
    properties(Hidden)
        poles, residues, invariant_zeros

        H_inf_norm = []
        H_inf_peakfreq = []
        H_2_norm = []
        
        HankelSingularValues
        T_bal, T_bal_inv
        ConGram, ConGramChol
        ObsGram, ObsGramChol

        isdescriptor=0
        
        SimulationTime
        
        decay_time
        
        mor_info
        % mor_info must be a struct containing the fields 'time', 'method', 'orgsys'
    end
    
    methods
        function sys = sss(varargin)
            if nargin==1
                if strcmp(class(varargin{1}), 'ss')
                    sys_ss = varargin{1};
                    % convert ss to sss
                    sys = sss(sys_ss.A, sys_ss.B, sys_ss.C, sys_ss.D, sys_ss.E);
                    return
                end
            end
            if nargin<3
                error('Please specify matrices A, B, C.')
            end
            sys.A = sparse(varargin{1});
            if size(sys.A,1) ~= size(sys.A,2)
                error('A must be square.')
            end
            sys.B = sparse(varargin{2});
            if size(sys.B,1) ~= size(sys.A,1)
                error('A and B must have the same number of rows.')
            end
            sys.C = sparse(varargin{3});
            if size(sys.C,2) ~= size(sys.A,2)
                error('A and C must have the same number of columns.')
            end
            
            if nargin>=4
                if (~isscalar(varargin{4}) && ~isempty(varargin{4})) || any(any(varargin{4}~=0))
                    % non-zero D is given
                    if size(sparse(varargin{4}),2) ~= size(sys.B,2)
                        error('B and D must have the same number of columns.')
                    end
                    if size(sparse(varargin{4}),1) ~= size(sys.C,1)
                        error('C and D must have the same number of rows.')
                    end
                    sys.D = sparse(varargin{4});
                else
                    % zero or a zero matrix is given
                    sys.D = sparse(sys.p, sys.m);
                end
            else
                sys.D = sparse(sys.p, sys.m);
            end
            if nargin>=5 && ~isempty(varargin{5})
                if size(varargin{5},1) ~= size(varargin{5},2)
                    error('E must be square.')
                end
                if size(sys.A,1) ~= size(varargin{5},1)
                    error('A and E must have the same size.')
                end
                if any(any(sparse(varargin{5})-speye(sys.n)))
                    sys.E = sparse(varargin{5});
                    sys.isdescriptor = 1;
                else
                    sys.E = [];
                end
            end
        end
        
        function n = get.n(sys)
            % system order
            n = size(sys.A,1);
        end
        function m = get.m(sys)
            % number of inputs
            m = size(sys.B,2);
        end
        function p = get.p(sys)
            % number of outputs
            p = size(sys.C,1);
        end
        
        function is_mimo = get.is_mimo(sys)
            is_mimo = (sys.p>1)||(sys.m>1);
        end        

        function D = get.D(sys)
            D = sys.D;
            if isempty(D)
                D = sparse(sys.m,sys.p);
            end
        end
        function E = get.E(sys)
            E = sys.E;
            if isempty(E)
                E = speye(sys.n);
            end
        end
        
        function a = get.a(sys)
            a = sys.A;
        end
        function b = get.b(sys)
            b = sys.B;
        end
        function c = get.c(sys)
            c = sys.C;
        end
        function d = get.d(sys)
            d = sys.D;
        end
        function e = get.e(sys)
            e = sys.E;
        end
        
        function [A,B,C,D,E] = ABCDE(sys)
            % returns system matrices
            A=sys.A; B=sys.B; C=sys.C; D=sys.D; E=sys.E;
        end
        
        function str = disp_mor_info(sys)
            if isempty(sys.mor_info)
                str=[];
            else
                str = ['Created by MOR on ' datestr(sys.mor_info.time) '.' char(10) ...
                    'Reduction Method: ' sys.mor_info.method  '.' char(10) ...
                    'Original system: ' sys.mor_info.orgsys '.'];
            end
        end
        
        function varargout = disp(sys)
            % display system information
            str = 'Sparse State Space Model';
            str = [str '.' char(10)];
            if sys.n==1
                str = [str '1 state, '];
            else
                str = [str num2str(sys.n) ' states, '];
            end
            if sys.m==1
                str = [str '1 input, '];
            else
                str = [str num2str(sys.m) ' inputs, '];
            end
            if sys.p==1
                str = [str '1 output.'];
            else
                str = [str num2str(sys.p) ' outputs.'];
            end
            if ~isempty(sys.mor_info)
                str = [str char(10) sys.disp_mor_info];
            end
            if nargout>0
                varargout = {str};
            else
                str = strrep(str, char(10), [char(10) '  ']);
                disp(['  ' str char(10)]);
            end
        end
        
        function [varargout] = subsref(sys, arg)
            % Returns selected I/O-channel of a sparse LTI MIMO system
            if strcmp(arg(1).type, '()')
                if length(arg.subs)==2
                    varargout = {sss(sys.A, sys.B(:,arg.subs{2}), sys.C(arg.subs{1},:), sys.D(arg.subs{1},arg.subs{2}), sys.E)};
                end
            else
                [varargout{1:nargout}] = builtin('subsref',sys,arg);
            end
        end
        
        function p = eig(sys)
            if isempty(sys.poles)
                if sys.isdescriptor
                    sys.poles = eig(full(sys.a), full(sys.e));
                else
                    sys.poles = eig(full(sys.a));
                end
            end
            p = sys.poles;
        end
        
        function p = eigs(sys,varargin)
            if sys.isdescriptor
                p = eigs(sys.A,sys.E,varargin{:});
            else
                p = eigs(sys.A,varargin{:});
            end
        end
        
        function sys = resolve_dae(sys, varargin)
            if nargin==1
                sys.A = sys.E\sys.A;
                sys.B = sys.E\sys.B;
                sys.E = [];
                sys.isdescriptor = 0;
            elseif varargin{1}==1
                %***
            elseif varargin{1}==2
            end                    
        end
        
    end
    
    methods(Access = private)
    end
    
end

