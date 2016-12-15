function [sysr, varargout] = tbr(sys, varargin)
% TBR - Performs model order reduction by Truncated Balanced Realization
%
% Syntax:
%       sys                 = TBR(sys)
%       sysr                = TBR(sys,q)
%       [sysr,V,W]          = TBR(sys,q)
%       [sysr,V,W,hsv]      = TBR(sys,q)
%       [sysr,V,W,hsv,R,L]	= TBR(sys,q)
%       [sysr,...]          = TBR(sys,Opts)
%       [sysr,...]          = TBR(sys,q,Opts)
%       [sysr,...]          = TBR(sys,q,R,L,Opts)
%
% Description:
%       Computes a reduced model by balanced truncation,
%       i.e. by transforming the system to a balanced realization where all
%       states are equally controllable and observable and selecting only
%       the modes responsible for the highest energy transfer in
%       system [1]. 
%
%       If the value of q equals the system order, then TBR computes a balanced 
%       realization of the system.
%
%       Hankel singular values and the matrices for transformation to
%       balanced realization are stored in the sss object sys.
%
%       If a reduction order q is passed to the function, the reduced
%       system will be of this size (type 'tbr' or 'matchDcGain') or smaller depending on
%       'rctol' (type 'adi') with the options 'hsvTol' and 'redErr' ignored.
%       Use 'forceOrder' to keep the desired order with adi. If not, the 
%       option reduction error 'redErr' is crucial. This error is defined 
%       as two times the sum of all Hankel-Singular values truncated after 
%       reduction. To avoid this option it can be set to zero 
%       ('redErr'=0). If so, the Hankel-Singular values (satisfying the 
%       option 'hsvTol') will be plotted for the user to enter a desired 
%       reduction order.
%
%       If the option 'type' is set to 'adi', a low rank approximation of the
%       Cholesky factor is performed. If the option 'type' is not 
%       defined, ADI is applied to systems with sys.n>500.
%
%       If the option 'type' is set to 'matchDcGain', then a
%       residualization is computed to match the DC gain of the original
%       model. Note that this is only possible using direct methods (tbr)
%       and not with adi.
%
%//Note: When ADI is used, only a small number of Hankel-Singular values
%       are computed. To determine the reduction error, the unknown
%       Hankel-Singular values are assumed to have the same value as the 
%       last one computed (worst-case scenario).
%
% Input Arguments:
%		*Required Input Arguments:*
%       -sys:   an sss-object containing the LTI system
%		*Optional Input Arguments:*
%       -q:     order of reduced system
%       -R,L:   low rank Cholesky factors of the Gramian matrices
%       -Opts:              a structure containing following options
%           -.type:         select amongst different tbr algorithms
%                           [{'tbr'} / 'adi' / 'matchDcGain' ]
%           -.redErr:       upper bound of reduction error
%                           [{'0'} / positive float]
%           -.hsvTol:       tolerance for Hankel-Singular values
%                           [{'1e-15'} / positive float]
%           -.warnOrError:  display warnings or errors
%                           [{'warn'} / 'error' / '0']
%           -.lse:          solve linear system of equations (only for adi)
%                           [{'gauss'} / 'luChol']
%           -.rctol:        tolerance for relative change (only for adi)
%                           [{'1e-9'} / positive float]
%           -.forceOrder    return desired order (only for adi)
%                           [{'false'} / 'true']
%
% Output Arguments:
%       -sysr:  reduced system
%       -V,W:   projection matrices
%       -hsv:   Hankel singular values
%       -R,L:   low rank Cholesky factors of the Gramian matrices
%
% Examples:
%       To compute a balanced realization, use
%
%> sys = loadSss('building');
%> sysBal = tbr(sys,sys.n)
%
%       To performe balanced reduction, specify a reduced order q
%
%> sysr = tbr(sys,10);
%> bode(sys,'-b',sysr,'--r')
%
% See Also:
%       rk, modalMor, gram, balancmr, lyapchol
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
%               Alessandro Castagnotto, Lisa Jeschek
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  11 Dec 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

% Default execution parameters
Def.type        = 'tbr'; % select tbr method (tbr, adi, matchDcGain)
Def.redErr      = 0; % reduction error (redErr>2*sum(hsv(q+1:end)))
Def.hsvTol      = 1e-15; % hsv tolerance (hsv(q)<hsvTol)
Def.warnOrError = 'warn'; % display warnings or errors (0,'warn','error')
Def.lse         = 'gauss'; % usfs for adi ('gauss', 'luChol')
Def.rctol       = 1e-9; % ADI stopping criterion (relative change)
Def.forceOrder  = false; % ADI force order q

% check input for q and Opts
if nargin>1
    if nargin==2 && ~isa(varargin{1},'double')
        Opts=varargin{1};
    elseif nargin>=4
        q=varargin{1};
        R=varargin{2};
        L=varargin{3};
        if nargin==5
            Opts=varargin{4};
        end
    else
        q=varargin{1};
        if nargin==3
            Opts=varargin{2};
        end
    end
end

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
    if sys.n>500
        if isempty(sys.ConGramChol) && isempty(sys.ObsGramChol) && isempty(sys.ConGram) && isempty(sys.ObsGram)
            Opts.type='adi';
        end
    end
else
    if ~isfield(Opts,'type') && sys.n>500
        if isempty(sys.ConGramChol) && isempty(sys.ObsGramChol) && isempty(sys.ConGram) && isempty(sys.ObsGram)
            Opts.type='adi';
        else
            Opts.type='tbr';
        end
    end
    Opts = parseOpts(Opts,Def);
end

if sys.isDae
    error('tbr does not work with DAE systems.');
end

if strcmp(Opts.type,'adi')
    if sys.n<100 && strcmp(Opts.warnOrError,'error')
        error('System is too small for ADI (sys.n >= 100 required).');
    elseif sys.n<100 && strcmp(Opts.warnOrError,'warn')
        warning('System is too small for ADI. Trying without ADI...');
        Opts.type='tbr';
    end
end

if strcmp(Opts.type,'adi')
    lyapOpts.method='adi';
    lyapOpts.lse=Opts.lse;
    lyapOpts.rctol=Opts.rctol;
    if exist('q','var')
        lyapOpts.q=q;
    end
    if ~exist('R','var') || ~exist('L','var')
        [R,L]=lyapchol(sys,lyapOpts);
    end
    
    qAdi=min([size(L,1),size(R,1)]);

else
    lyapOpts.method='hammarling';
    if ~exist('R','var') || ~exist('L','var')
        [R,L]=lyapchol(sys,lyapOpts);
    end
end
    
% calculate balancing transformation and Hankel Singular Values
[K,S,M]=svd(full(R*L'),0);
hsv = diag(S);
sys.HankelSingularValues = real(hsv);
sys.TBalInv = R'*K*diag(ones(size(hsv))./sqrt(hsv));
sys.TBal = diag(ones(size(hsv))./sqrt(hsv))*M'*solveLse(sys.E',L',struct('lse','gauss'))';

qmax=size(sys.TBal,1);

% determine reduction order
if exist('q','var') || Opts.redErr>0
    if ~exist('q','var')
        if strcmp(Opts.type,'adi') && qmax<sys.n
            % worst case for unknown hsv
            hsvSum=2*real(hsv(qmax))/real(hsv(1))*(sys.n-qmax+1);
        else
            hsvSum=0;
        end
        for i=qmax:-1:0
            if hsvSum>Opts.redErr || i==0
                q=i+1;
                if q>qmax
                    q=qmax;
                end
                break;
            else
                hsvSum=hsvSum+2*real(hsv(i))/real(hsv(1));
            end
        end
        
        if strcmp(Opts.type,'adi') && qmax==qAdi
            warning(['Reduction order was set to q = ', num2str(q,'%d'),...
            ' to satisfy the upper bound for the reduction error. ',10,...
            'The upper bound can be unprecise due to the use of ADI.']);
        else
            warning(['Reduction order was set to q = ', num2str(q,'%d'),...
                ' to satisfy the upper bound for the reduction error. ']);
        end
    end
    if q>sys.n && sys.n==qmax
        if strcmp(Opts.warnOrError,'error')
            error('Reduction order exceeds system order.');
        elseif strcmp(Opts.warnOrError,'warn')
            warning(['Reduction order exceeds system order. q has been changed to ',...
                'the system order qmax = ', num2str(qmax,'%d'), '.']);
        end
        q=sys.n;
    end
    if q>qmax && strcmp(Opts.type,'adi')
        warning(['The reduction order was set to the size of the ADI ',...
            'iterates q = ',num2str(qmax,'%d'),'.',...
            ' To force the desired order (q = ',num2str(q,'%d'),'), use the option',...
            ' Opts.forceOrder.']);
        q=qmax;
    end
    if sum(hsv>=Opts.hsvTol*hsv(1))<q
        if strcmp(Opts.warnOrError,'error')
            error(['The reduction order of q = ', num2str(q,'%d'),' includes ',...
                'Hankel-Singular values smaller than the chosen tolerance (see Opts.hsv).']);
        elseif strcmp(Opts.warnOrError,'warn')
            warning(['The reduced system of desired order may not be a minimal ',...
                'realization and it may not be stable. The recommended reduced',...
                ' order is q = ',num2str(sum(hsv>=Opts.hsvTol*hsv(1)),'%d'),...
                ' (see Opts.hsvTol).']);
        end
    end
else
    qmax = min([sum(hsv>=Opts.hsvTol*hsv(1)), qmax]);
    h=figure(1);
    bar(1:qmax,abs(hsv(1:qmax)./hsv(1)),'r');
    title('Hankel Singular Values');
    xlabel('Order');
    ylabel({'Relative hsv decay';sprintf('abs(hsv/hsv(1)) with hsv(1)=%.4d',hsv(1))});
    set(gca,'YScale','log');
    set(gca, 'YLim', [-Inf;1.5]);
    set(gca, 'XLim', [0; qmax]);
    prompt=['Please enter the desired order: (0<= q <=', num2str(qmax,'%d)'),' '];
    q=input(prompt);
    if ishandle(h)
        close Figure 1;
    end
    if q<0 || round(q)~=q
        error('Invalid reduction order.');
    end
    if q>sys.n && qmax==sys.n
        if strcmp(Opts.warnOrError,'error')
            error('Reduction order exceeds system order.');
        elseif strcmp(Opts.warnOrError,'warn')
            warning(['Reduction order exceeds system order. q has been changed to ',...
                'the system order qmax = ', num2str(qmax,'%d'), '.']);
        end
        q=qmax;
    elseif q>qmax
        if strcmp(Opts.type,'adi') && qmax==qAdi
            if strcmp(Opts.warnOrError,'error')
                error(['Reduction order must be smaller than q = ', num2str(qmax,'%d'),...
                    ' due to ADI.']);
            elseif strcmp(Opts.warnOrError,'warn')
                warning(['A reduction is only possible to qmax = ', num2str(qmax,'%d'),...
                    ' due to ADI. q has been changed accordingly.']);
            end
        else
            if strcmp(Opts.warnOrError,'error')
                error(['Reduction order must be smaller than q = ', num2str(qmax,'%d'),...
                    ' due to Hankel-Singular values smaller than the given tolerance',...
                    ' (see Opts.hsvTol).']);
            elseif strcmp(Opts.warnOrError,'warn')
                warning(['q has been changed to qmax = ', num2str(qmax,'%d'),...
                        ' due to Hankel-Singular values smaller than the given '...
                        'tolerance (see Opts.hsvTol).']);
            end
        end
        q=qmax;
    end
end

V = sys.TBalInv(:,1:q);
W = sys.TBal(1:q,:)';

% Storing additional parameters
%Storing additional information about tbr reduction in the object 
%containing the reduced model:
%   1. Define a new field for the Opts struct and write the information
%      that should be stored to this field
%   2. Adapt the method "parseParamsStruct" of the class "ssRed" in such a
%      way that the new defined field passes the check
Opts.originalOrder = sys.n;
Opts.hsv = hsv;

switch Opts.type
    case 'tbr'
        if isa(sys,'ssRed')
            sysr = ssRed('tbr',Opts,W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V, sys.reductionParameters);
        else
            sysr = ssRed('tbr',Opts,W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);
        end
    case 'matchDcGain'
        W=sys.TBalInv;
        V=sys.TBal;
        ABal=W*sys.A*V;
        BBal=W*sys.B;
        CBal=sys.C*V;

        [A11,A12,A21,A22] = partition(ABal,q);
        B1=BBal(1:q,:);B2=BBal(q+1:end,:);
        C1=CBal(:,1:q);C2=CBal(:,q+1:end);
        
        if rcond(A22)<eps
            if strcmp(Opts.warnOrError,'warn')
                % don't display Matlab's warning several times, but display 
                % only 1 warning that informs user of the consequences
                warning('tbr:rcond','MatchDcGain might be inaccurate because of a nearly singular matrix.');
                warning('off','MATLAB:nearlySingularMatrix');
            elseif strcmp(Opts.warnOrError,'error')
                error('tbr:rcond','Nearly singular matrix in matchDcGain.');
            end
        end

        lse0=solveLse(A22',A12')';
        ARed=A11-lse0*A21;
        BRed=B1-lse0*B2; 

        if sys.isDescriptor
            EBal=W'*sys.E*V;
            E11=EBal(1:q,1:q); % E12=E_bal(1:q,1+q:end);
            E21=EBal(1+q:end,1:q); % E22=E_bal(q+1:end,1+q:end);
            ERed=E11-lse0*E21;
            lse1=solveLse(A22',C2')';
            lse2=solveLse(ERed',E21')';
            CRed=C1-lse1*A21+C2*A22*lse2*ARed;
            DRed=sys.D-lse1*B2+lse1*lse2*BRed;
            if isa(sys,'ssRed')
                sysr = ssRed('tbr',Opts,ARed, BRed, CRed, DRed, ERed,sys.reductionParameters);
            else
                sysr = ssRed('tbr',Opts,ARed, BRed, CRed, DRed, ERed);
            end
        else % Er=I
            CRed=C1-solveLse(A22',C2')'*A21;
            DRed=sys.D-solveLse(A22',C2')'*B2;
            if isa(sys,'ssRed')
                sysr = ssRed('tbr',Opts,ARed, BRed, CRed, DRed,sys.reductionParameters);
            else
                sysr = ssRed('tbr',Opts,ARed, BRed, CRed, DRed); 
            end
        end
        
        warning('on','MATLAB:nearlySingularMatrix');    
    case 'adi'
        if isa(sys,'ssRed')
            sysr = ssRed('tbr',Opts,W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V,sys.reductionParameters);
        else
            sysr = ssRed('tbr',Opts,W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);
        end
end

%   Rename ROM
sysr.Name = sprintf('%s_%i_tbr',sys.Name,sysr.n);

if nargout>1
    varargout{1} = V;
    varargout{2} = W;
    varargout{3} = real(hsv);
    if nargout>3
        varargout{4} = R;
        varargout{5} = L;
    end
end
