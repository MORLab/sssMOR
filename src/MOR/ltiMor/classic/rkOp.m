function varargout = rkOp(varargin)
% RKOP - Determination of optimal expansion point for Laguerre series
%
% Syntax:
%       sOpt = RKOP(sys)
%       sOpt = RKOP(h,t)
%       [sysr, V, W, sOpt] = RKOP(sys,q)
%       [sysr, V, W, sOpt] = RKOP(sys,q,Opts)
%
% Description:
%       This function deteremines an optimal expansion point for 
%       Krylov-based model order reduction. 
%
%       The computed point is the optimal Laguerre parameter that 
%       guarantees an optimal approximation of the impulse response of the
%       original system.
%
%       If a reduction order q is specified, a reduced model is computed at
%       the optimal expansion point using Rational Krylov.
%
% Input Arguments:  
%       *Required Input Arguments:*
%       -sys:			full oder model (ss or sss)
%       -h:             impulse response vector of the system
%       -t:             time vector corresponding to h
%       *Optional Input Arguments:*
%       -q:             reduction order
%       -Opts:			structure with execution parameters
%			-.rk:       reduction type
%						[{'twoSided'} / 'input' / 'output']
%
% Output Arguments:      
%       -sOpt:          optimal expansion point
%       -sysr:          reduced system
%       -V, W:          projection matrices spanning Krylov subspaces
%
% Examples:
%       This code computes the optimal expansion point and a reduced system
%       of order q=10 for the benchmark model 'building'.
%
%> sys = loadSss('building')
%> [sysr,V,W,sOpt] = rkOp(sys,10);
%> bode(sys,'-',sysr,'--r');
%
% See Also: 
%       rk, rkIcop, irka, arnoldi
%
% References:
%       * *[1] R. Eid: Time domain Model Reduction by Moment Matching, 
%               Ph.D thesis, Institute of Automatic Control, Technische 
%               Universitaet Muenchen, 2009.
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
% Authors:      Heiko Panzer (heiko@mytum.de), Rudy Eid
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  28 Jun 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if length(varargin)==3 && isa(varargin{3},'struct')
    Opts=varargin{3};
    varargin=varargin(1:2);
end

%% Parse the inputs
%   Default execution parameters
Def.rk = 'twoSided'; % 'twoSided','input','output'

% create the options structure
if ~exist('Opts','var') || isempty(fieldnames(Opts))
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

% compute sOpt
if isa(varargin{1},'double') && length(varargin)==2 && length(varargin{1})==length(varargin{2}) % from impulse response
    h = varargin{1};
    th = varargin{2};
    sOpt=zeros(size(h,2),size(h,3));
    for j=1:size(h,2)
        for k=1:size(h,3)
            m1=0;
            m2=0;
            for i=1:length(th)
                m1=m1+th(i)*h(i,j,k)^2;
                if i<length(th)
                    m2=m2+th(i)*((h(i+1,j,k)-h(i,j,k))/(th(2)-th(1)))^2;
                end
            end
            sOpt(j,k) = sqrt(m2/m1);
        end
    end
    
elseif isa(varargin{1},'sss') || isa(varargin{1},'ssRed') || isa(varargin{1},'ss') % from lyapunov equation
    if isa(varargin{1},'ss')
        sys = sss(varargin{1});
    else
        sys = varargin{1};
    end
    sOpt=zeros(sys.p,sys.m);
    for i=1:sys.p
        for j=1:sys.m
            try
                if sys.isDescriptor
                    P=lyap(sys.A,sys.B(:,j)*sys.B(:,j)',[],sys.E);
                    Y=lyap(sys.A,P,[],sys.E);
                else
                    P=lyap(sys.A,sys.B(:,j)*sys.B(:,j)');
                    Y=lyap(sys.A,P);
                end
            catch ex
                error(['Error during calculation of rkOp: ' ex.message]);
            end

            sOpt(i,j)=sqrt((sys.C(i,:)*sys.A*Y*sys.A'*sys.C(i,:)')/(sys.C(i,:)*Y*sys.C(i,:)'));
            if nnz(sOpt(i,j)<0) || ~isreal(sOpt(i,j))
                warning(['The optimal expansion point is negative or complex.'...
                    'It has been replaced by its absolute value.']);
                sOpt(i,j) = abs(sOpt(i,j));
            end
        end
    end
else
    error('Wrong input.');
end

% return reduced system
if length(varargin)==2 && ~isa(varargin{1},'double') && isscalar(varargin{2})
    q=varargin{2};
    if sys.isSiso
        switch(Opts.rk)
            case 'twoSided'
                [varargout{1},varargout{2},varargout{3}] = rk(sys,[sOpt;q],[sOpt;q]);
            case 'input'
                [varargout{1},varargout{2},varargout{3}] = rk(sys,[sOpt;q]);
            case 'output'
                [varargout{1},varargout{2},varargout{3}] = rk(sys,[],[sOpt;q]);
            otherwise
                error('Wrong Opts.');
        end
    else        
        switch(Opts.rk)
            case 'twoSided'
                sOpt=sOpt';
                sOpt=sOpt';
                tempLt=[];
                Rt=[];
                Lt=[];

                for j=1:sys.m
                   Rt=blkdiag(Rt,ones(1,q*sys.p));
                end
                
                for j=1:sys.p
                    tempLt=blkdiag(tempLt,ones(1,q));
                end                    

                for j=1:sys.m
                    Lt=[Lt,tempLt];
                end
                
                [varargout{1},varargout{2},varargout{3}] = rk(sys,[sOpt(:)';ones(1,sys.m*sys.p)*varargin{2}],[sOpt(:)';ones(1,sys.m*sys.p)*varargin{2}],Lt,Lt);
            case 'input'
                sOpt=sOpt';
                Rt=[];

                for j=1:sys.m
                   Rt=blkdiag(Rt,ones(1,q*sys.p));
                end
                
                [varargout{1},varargout{2},varargout{3}] = rk(sys,[sOpt(:)';ones(1,sys.m*sys.p)*varargin{2}],Rt);
            case 'output'
                sOpt=sOpt';
                tempLt=[];
                Lt=[];
                
                for j=1:sys.p
                    tempLt=blkdiag(tempLt,ones(1,q));
                end                    

                for j=1:sys.m
                    Lt=[Lt,tempLt];
                end
                
                [varargout{1},varargout{2},varargout{3}] = rk(sys,[],[sOpt(:)';ones(1,sys.m*sys.p)*varargin{2}],[],Lt);
            otherwise
                error('Wrong Opts.');
        end
    end
    varargout{4}=sOpt;
    
elseif length(varargin)==1 || (length(varargin)==2 && length(varargin{1})==length(varargin{2}))
    varargout{1} = sOpt;
else
    error('Wrong input.');
end
    
    
