function [sysr, varargout] = tbr(sys, varargin)
% TBR - Performs model order reduction by Truncated Balanced Realization
%
% Syntax:
%       sys				= TBR(sys)
%       sysr			= TBR(sys,q)
%       [sysr,V,W]		= TBR(sys,q)
%       [sysr,V,W,hsv]	= TBR(sys,q)
%       [sysr,...]      = TBR(sys,Opts)
%       [sysr,...]      = TBR(sys,q,Opts)
%
% Description:
%       Computes a reduced model of order q by balancing and truncation,
%       i.e. by transforming the system to a balanced realization where all
%       states are equally controllable and observable and selecting only
%       the first q modes responsible for the highest energy transfer in
%       system [1]. 
%
%       If q is specified as system order, then TBR computes a balanced 
%       realization of the system.
%
%       Hankel singular values and the matrices for transformation to
%       balanced realization are stored in the sss object sys.
%
%       If a reduction order q is passed to the function, the reduced
%       system will be of this size with the options 'hsvTol' and
%       'redErr' ignored. If not, the option reduction error 'redErr' is 
%       crucial. This error is defined as two times the sum of all 
%       Hankel-Singular values bigger than the reduction order. To avoid 
%       this option it can be set to zero ('redErr'=0). If so, the
%       Hankel-Singular values (satisfying the option 'hsvTol') will be
%       plotted for the user to enter a desired reduction order.
%
%       If the option 'type' is set to 'adi', a low rank approximation of the
%       cholseky factor is performed. If the option 'adi' is not 
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
%       -Opts:              a structure containing following options
%           -.real:         return real reduced system
%                           [{'real'} / '0']
%           -.type:         select amongst different tbr algorithms
%                           [{'tbr'} / 'adi' / 'matchDcGain' ]
%           -.redErr:       upper bound of reduction error
%                           [{'0'} / positive float]
%           -.hsvTol:       tolerance for Hankel-Singular values
%                           [{'1e-15'} / positive float]
%           -.warnOrError:  display warnings or errors
%                           [{'warn'} / 'error' / '0']
%
% Output Arguments:
%       -sysr:  reduced system
%       -V,W:   projection matrices
%       -hsv:   Hankel singular values
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
%               Alessandro Castagnotto, Lisa Jeschek
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  16 Feb 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

% Default execution parameters
Def.type = 'tbr'; % select tbr method (tbr, adi, matchDcGain)
Def.redErr = 0; % reduction error (redErr>2*sum(hsv(q+1:end)))
Def.hsvTol = 1e-15; % hsv tolerance (hsv(q)<hsvTol)
Def.real = 'real'; % real reduced system ('0', 'real')
Def.warnOrError = 'warn'; % display warnings or errors (0,'warn','error')

% check input for q and Opts
if nargin>1
    if nargin==2 && ~isa(varargin{1},'double')
        Opts=varargin{1};
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
    if ~isfield(Opts,'adi') && sys.n>500
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
    % lyapack options (lyaOpts):
    %   method='heur': find ADI parameters with a heuristic method
    %   (l0,kp,km)=(20,50,25): shifts p consist of l0 or l0+1 values of kp('LM')+km('SM') Ritz values
    %   zk='Z': return cholesky factor of solution X=Z*Z' of lyapunov equation
    %   rc='R'('C'): return real (complex) Z
    %   adi.type='B'('C'): layponov equation type ('B','C')
    %   adi.max_it=100: maximum number of iterations for ADI iteration (stopping criteria)
    %   adi.min_res=0: minimum residual for ADI iteration - expensive(stopping criteria)
    %   adi.with_rs='N': (S/N) stop ADI iteration if stagnation of error - very expensive (stopping criteria)
    %   adi.min_in=1e-12: tolerance for difference between ADI iterates - inexpensive(stopping criteria)
    %   adi.cc_upd=0: column compression parameter (0=never)
    %   adi.cc_tol=sqrt(eps): column compression tolerance (default=sqrt(eps))
    %   adi.info=0; information level
    %   usfs.au: no descriptor (E=I), sparse, possibly unsymmetric
    %   usfs.as: no descriptor (E=I), sparse, sys.A symmetric
    %   usfs.munu: descriptor (E~=I), sparse, possibly unsymmetric
    %   usfs.msns: descriptor (E~=I), sparse, sys.A and sys.E symmetric
    
    % required changes of lyapack functions (functions will error without changes):
    %   routine\lp_lradi.m - line 427: insert 'full': svd(full(...)) - else 
    %       error if Opts.real='real' because svd(sparse)
    %   usfs\as_s - line 37: replace 'LP_U' with 'LP_UC' - else error if
    %       Opts.sym='sym' and system is not descriptor because LP_U does not exist
    
    % optional changes of lyapack functions (function work without changes):
    %   usfs\au_l_i, usfs\au_l, usfs\au_l_d: [L,U]=lu(A) replaced with [L,U,a,o,S] = lu(A)
    %   usfs\munu_l_i, usfs\munu_l, usfs\munu_l_d: [L,U]=lu(A) replaced with [L,U,a,o,S] = lu(A)
    %   usfs\au_s_i, usfs\au_s, usfs\au_s_d: [L,U]=lu(A) replaced with [L,U,a,o,S] = lu(A)
    %   usfs\munu_s_i, usfs\munu_s, usfs\munu_s_d: [L,U]=lu(A) replaced with [L,U,a,o,S] = lu(A)
 
    lyaOpts.l0=20;
    lyaOpts.kp=50;
    lyaOpts.km=25;
    lyaOpts.method='heur';
    lyaOpts.zk='Z';
    
    if strcmp(Opts.real,'real')
        lyaOpts.rc='R';
    else
        lyaOpts.rc='C';
    end

    lyaOpts.adi=struct('type','B','max_it', 300,'min_res',0,'with_rs','N',...
        'min_in',1e-12,'min_end',0,'info',0,'cc_upd',0,'cc_tol',0);
    
    if exist('q','var') %size of cholesky factor [sys.n x q] -> qmax=q
        lyaOpts.adi.max_it=q;
        lyaOpts.adi.min_end=1;
    end

    if sys.isSym
        if ~sys.isDescriptor
            lyaOpts.usfs=struct('s','as_s','m','as_m');
            [A0,B0,C0]=as_pre(sys.A,sys.B,sys.C); %preprocessing: reduce bandwith of A
            as_m_i(A0);
            as_l_i;
            p=lp_para(as,[],[],lyaOpts, ones(length(B0),1)); %determine ADI parameters p (Ritz values of A)
            lyaOpts.p=p.p;
            as_s_i(lyaOpts.p);
        else
            lyaOpts.usfs=struct('s','msns_s','m','msns_m');
            [M0,MU0,N0,B0,C0]=msns_pre(sys.E,sys.A,sys.B,sys.C);
            msns_m_i(M0,MU0,N0); 
            msns_l_i;
            p=lp_para(msns,[],[],lyaOpts,ones(length(B0),1));
            lyaOpts.p=p.p;
            msns_s_i(lyaOpts.p);
        end
    else
        if ~sys.isDescriptor
            lyaOpts.usfs=struct('s','au_s','m','au_m');
            [A0,B0,C0]=au_pre(sys.A,sys.B,sys.C);
            au_m_i(A0);
            au_l_i;
            p=lp_para(au,[],[],lyaOpts, ones(length(B0),1));
            lyaOpts.p=p.p;
            au_s_i(lyaOpts.p);
        else
            lyaOpts.usfs=struct('s','munu_s','m','munu_m');
            [M0,ML0,MU0,N0,B0,C0]=munu_pre(sys.E,sys.A,sys.B,sys.C);
            munu_m_i(M0,ML0,MU0,N0); 
            munu_l_i;
            p=lp_para(munu,[],[],lyaOpts,ones(length(B0),1));
            lyaOpts.p=p.p;
            munu_s_i(lyaOpts.p);
        end
    end 
    
    % ADI solution of lyapunov equation
    [R,Ropts]=lp_lradi([],[],B0,lyaOpts);
    if sys.isSym && norm(full(sys.B-sys.C'))==0
        L=R;
        Lopts=Ropts;
    else
        lyaOpts.adi.type='C';
        [L,Lopts]=lp_lradi([],[],C0,lyaOpts);
    end
    
    if lyaOpts.adi.min_end==1
        q_min_in=max([Ropts.adi.min_iter,Lopts.adi.min_iter]);
        if q_min_in>0 && q_min_in < q
        warning(['After q=', num2str(q_min_in,'%d'),...
            ' the contribution of the ADI iterates was very small.']);
        end
    end

    qmax=min([size(R,2),size(L,2)]); %make sure R and L have the same size
    R=R(:,1:qmax);
    L=L(:,1:qmax);
    qmaxAdi=qmax;
    
    % calculate balancing transformation and Hankel Singular Values
    [K,S,M]=svd(full(L'*R));
    hsv = diag(S);
    sys.HankelSingularValues = real(hsv);
    sys.TBalInv = R*M/diag(sqrt(hsv));
    sys.TBal = diag(sqrt(hsv))\K'*L';
    
else
    qmax=sys.n;
    % Is Controllability Gramian available?
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

    % Is Observability Gramian available?
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
    [K,S,M]=svd(full(R*L'));
    hsv = diag(S);
    sys.HankelSingularValues = real(hsv);
    sys.TBalInv = R'*K/diag(sqrt(hsv));
    sys.TBal = diag(sqrt(hsv))\M'*L/sys.E;
end

% store system
if inputname(1)
    assignin('caller', inputname(1), sys);
end

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
        
        if strcmp(Opts.type,'adi')
            warning(['Reduction order was set to q = ', num2str(q,'%d'),...
            ' to satisfy the upper bound for the reduction error. ',10,...
            'The upper bound can be unprecise due to the use of ADI.']);
        else
            warning(['Reduction order was set to q = ', num2str(q,'%d'),...
                ' to satisfy the upper bound for the reduction error. ']);
        end
    end
    if q>sys.n
        if strcmp(Opts.warnOrError,'error')
            error('Reduction order exceeds system order.');
        elseif strcmp(Opts.warnOrError,'warn')
            warning(['Reduction order exceeds system order. q has been changed to ',...
                'the system order qmax = ', num2str(qmax,'%d'), '.']);
        end
        q=sys.n;
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
        if strcmp(Opts.type,'adi') && qmax==qmaxAdi
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

switch Opts.type
    case 'tbr'
        sysr = sss(W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);
    case 'matchDcGain'
        W=sys.TBalInv;
        V=sys.TBal;
        ABal=W*sys.A*V;
        BBal=W*sys.B;
        CBal=sys.C*V;

        [A11,A12,A21,A22] = partition(ABal,q);
        B1=BBal(1:q,:);B2=BBal(q+1:end,:);
        C1=CBal(:,1:q);C2=CBal(:,q+1:end);

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
    
    case 'adi'
        sysr=sss(W'*feval(lyaOpts.usfs.m,'N',V),W'*B0,C0*V,sys.D);

        %delete global data
        as_m_d; as_l_d; as_s_d(p.p);
        msns_m_d; msns_l_d; msns_s_d(p.p)
        au_m_d; au_l_d; au_s_d(p.p);
        munu_m_d; munu_l_d; munu_s_d(p.p)
end

if nargout>1
    varargout{1} = V;
    varargout{2} = W;
    varargout{3} = real(hsv);
end
