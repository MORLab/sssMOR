function [sysr, varargout] = tbr(sys, varargin)
% TBR - Performs model order reduction by Truncated Balanced Realization
%
% Syntax:
%       sys				= TBR(sys)
%       sysr			= TBR(sys,q)
%       sysr            = TBR(sys,q,Opts)
%       [sysr,V,W]		= TBR(sys,q)
%       [sysr,V,W,hsv]	= TBR(sys,q)
%
% Description:
%       Computes a reduced model of order q by balancing and truncation,
%       i.e. by transforming the system to a balanced realization where all
%       states are equally controllable and observable and selecting only
%       the first q modes responsible for the highest energy transfer in
%       system [1]. 
%
%       If q is not specified, then TBR computes only a balanced
%       realization of the system without truncation.
%
%       Hankel singular values and the matrices for transformation to
%       balanced realization are stored in the sss object sys.
%
%
% Input Arguments:
%		*Required Input Arguments:*
%       -sys:   an sss-object containing the LTI system
%		*Optional Input Arguments:*
%       -q:     order of reduced system
%       -Opts:			a structure containing following options
%			-.type:		option for reduction method;
%						[{'tbr'} / 'matchDcGain']
%
% Output Arguments:
%       -sysr:  reduced system
%       -V,W:   (opt.) projection matrices (only if q is given!)
%       -hsv:   Hankel singular values
%
%//Note: If no q is given, the balancing transformation and calculation of
%		the Hankel singular values is performed without subsequent model
%		reduction.
%
% Examples:
%       To compute a balanced realization, use
%
%> sys = loadSss('building');
%> sysBal = tbr(sys)
%
%       To performe balanced reduction, specify a reduced order q
%
%> sysr = tbr(sys,10);
%> bode(sys,'-b',sysr,'--r')
%
% See Also:
%       rk, modalMor, gram, balancmr
%
% References:
%       * *[1] Moore (1981)*, Principal component analysis in linear systems: controllability,
%       observability and model reduction
%       * *[2] Antoulas (2005)*, Approximation of large-scale dynamical systems
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
% Authors:      Heiko Panzer, Sylvia Cremer, Rudy Eid, 
%               Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  02 Dec 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

% Default execution parameters
Def.type = 'tbr';

orderGiven = 0;

%Parse the optional input parameters
if nargin == 3
    
    if isempty(varargin{2})
        Opts = Def;
    else
        Opts = parseOpts(varargin{2},Def);
    end
    
    if isnumeric(varargin{1})
       orderGiven = 1;
    end
    
elseif nargin == 2
    
    if isnumeric(varargin{1})
       Opts = Def;
       orderGiven = 1;
    else
       Opts = parseOpts(varargin{1},Def);
    end
    
elseif nargin == 1
    
    Opts = Def;
    
end



%% Is Controllability Gramian available?
if isempty(sys.ConGramChol)
    if isempty(sys.ConGram)
        % No, it is not. Solve Lyapunov equation.
        try
            if sys.isDescriptor
                sys.ConGramChol = lyapchol(sys.A,sys.B,sys.E);
            else
                sys.ConGramChol = lyapchol(sys.A,sys.B);
            end
            R = sys.ConGramChol;
        catch ex
            warning(ex.message, 'Error in lyapchol. Trying without Cholesky factorization...')
            if sys.isDescriptor
                try
                    sys.ConGram = lyap(sys.A, sys.B*sys.B', [], sys.E);
                catch ex2
                    warning(ex2.message, 'Error solving Lyapunov equation. Premultiplying by E^(-1)...')
                    tmp = sys.E\sys.B;
                    sys.ConGram = lyap(sys.E\sys.A, tmp*tmp');
                end
            else
                sys.ConGram = lyap(full(sys.A), full(sys.B*sys.B'));
            end
            try
                R = chol(sys.ConGram);
            catch ex2
                myex = MException(ex2.identifier, ['System seems to be unstable. ' ex2.message]);
                throw(myex)
            end
        end
    else
        R = chol(sys.ConGram);
    end
else
    R = sys.ConGramChol;
end


%% Is Observability Gramian available?
if isempty(sys.ObsGramChol)
    if isempty(sys.ObsGram)
        % No, it is not. Solve Lyapunov equation. 
       try
            if sys.isDescriptor
                L = lyapchol(sys.A'/sys.E', sys.C');
            else
                L = lyapchol(sys.A',sys.C');
            end
            sys.ObsGramChol = sparse(L);
        catch ex
            warning(ex.message, 'Error in lyapchol. Trying without Cholesky factorization...')
            if sys.isDescriptor
                sys.ObsGram = lyap(sys.A'/sys.E', sys.C'*sys.C);
            else
                sys.ObsGram = lyap(full(sys.A'), full(sys.C'*sys.C));
            end
            try
                L = chol(sys.ObsGram);
            catch ex
                myex = MException(ex2.identifier, ['System seems to be unstable. ' ex2.message]);
                throw(myex)
            end
       end
    else
        L = chol(sys.ObsGram);
    end
else
    L = sys.ObsGramChol;
end


% calculate balancing transformation and Hankel Singular Values
[K,S,M]=svd(R*L');
hsv = diag(S);
sys.HankelSingularValues = real(hsv);
sys.TBalInv = R'*K/diag(sqrt(hsv));
sys.TBal = diag(sqrt(hsv))\M'*L/sys.E;


% store system
if inputname(1)
    assignin('caller', inputname(1), sys);
end

%% if MOR is to be performed, calculate V, W and reduced system
if orderGiven == 0
    h=figure(1);
    bar(1:sys.n,abs(hsv./hsv(1)),'r');
    title('Hankel Singular Values');
    xlabel('Order');
    ylabel({'Relative hsv decay';sprintf('abs(hsv/hsv(1)) with hsv(1)=%.4d',hsv(1))});
    set(gca,'YScale','log');
    set(gca, 'YLim', [-Inf;1.5]);
    prompt='Please enter the desired order: (>=0) ';
    q=input(prompt);
    if ishandle(h)
        close Figure 1;
    end
    if q<0 || round(q)~=q
        error('Invalid reduction order.');
    end
else
    q=varargin{1};
end
if q>sys.n
    warning('Reduction order exceeds system order. It is replaced by the system order.');
    q=sys.n;
end

%% Perform the reduction with the selected method

if strcmp(Opts.type,'tbr')
    
    V = sys.TBalInv(:,1:q);
    W = sys.TBal(1:q,:)';

    sysr = sss(W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);
    
elseif strcmp(Opts.type,'matchDcGain')
    
    W=sys.TBalInv;
    V=sys.TBal;
    ABal=W*sys.A*V;
    BBal=W*sys.B;
    CBal=sys.C*V;

    A11=ABal(1:q,1:q);
    A12=ABal(1:q,q+1:end);
    A21=ABal(q+1:end,1:q);
    A22=ABal((q+1):end,(q+1):end);
    B1=BBal(1:q);B2=BBal(q+1:end);
    C1=CBal(1:q);C2=CBal(q+1:end);

    ARed=A11-A12/A22*A21;
    BRed=B1-A12/A22*B2; 
    
    if sys.isDescriptor
        EBal=W'*sys.E*V;
        E11=EBal(1:q,1:q); % E12=E_bal(1:q,1+q:end);
        E21=EBal(1+q:end,1:q); % E22=E_bal(q+1:end,1+q:end);
        ERed=E11-A12/A22*E21;
        CRed=C1-C2/A22*A21+C2*A22i*E21/ERed*ARed;
        DRed=sys.D-C2/A22*B2+C2/A22*E21/ERed*BRed;
        sysr = sss(ARed, BRed, CRed, DRed, ERed);
    else % Er=I
        CRed=C1-C2/A22*A21;
        DRed=sys.D-C2/A22*B2;
        sysr = sss(ARed, BRed, CRed, DRed);
    end
    
else
   error('Wrong option. Only options "tbr" and "matchDcGain" are available'); 
end


if nargout>1
    varargout{1} = V;
    varargout{2} = W;
    varargout{3} = real(hsv);
end
