function [sysr, polesvf3]  = vectorFitting(sys,n,s0,Opts)
% vectorFitting - rational approximation of frequency response
%
% Syntax:
%       sysr                            = VECTORFITTING(sys, n, s0)
%       [sysr,polesvf3]                 = VECTORFITTING(sys, n, s0)
%       [sysr,polesvf3]                 = VECTORFITTING(sys, n, s0,Opts)
%
% Description:
%       
%
% Input Arguments:  
%       *Required Input Arguments:*
%       -sys:			full oder model (sss)
%       -s0:			vector of initial shifts
%
%       *Optional Input Arguments:*
%       -Opts:			structure with execution parameters
%			-.tol:		convergence tolerance;
%						[{1e-3} / positive float]
%
% Output Arguments:      
%       -sysr:              reduced order model (sss)
%     
% Examples:
%
% See Also: 
%       arnoldi, rk, isH2opt
%
% References:
%       * *[1] Gustavsen and Semlyen (1999)*, Rational approximation of frequency domain responses by vector fitting
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
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  31 Aug 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


    %%  Initialize by plotting the frequencies for which the data is available
    if Opts.plot
        fh = figure('Name','Sampling frequencies for VF');
        plot(complex(s0),'o');xlabel('real');ylabel('imag')
    end

    %% Preprocessing

    % take only one complex conjugate partner (the one with imag(s0)>0)
    s0 = cplxpair(s0); idx = find(imag(s0)); s0(idx(1:2:end)) = [];

    % adapt the ID size to the imput size
    m = size(sys.b,2); p = size(sys.c,1);
    if m>1, n = round(n/m);end %avoid blowing-up for MIMO
    % if mod(n,2) ~= 0, n = n-1; end   %make even
    nSample = length(s0);

    %%  Collect frequency data
    fprintf(1,'Generating frequency data...');
    f = zeros(p,m,nSample);
    for iW = 1:nSample;
        f(:,:,iW) = sys.C*((s0(iW)*sys.E-sys.A)\sys.B)+sys.D;
    end
    fprintf(1,'...done!\n');

    %%  Initialize poles
    polesvf3 = initializePoles(sys,Opts.vf.poles,n,Opts.wLims);
    if Opts.plot
        figure(fh); hold on; plot(complex(polesvf3),'rx');
    end

    %%  Adpatively choose reduced order

    %   VF solves iteratively a problem of the form Ax=b
    %   An indicator for a correct choice of input space dimension for x is
    %   given by the rank of A. 
    %   We compute the A matrix only for the first column of G(s) and assume
    %   similar information content for the other columns
    %
    %   The A matrix is built up of
    %      | A1          Atilde1 |
    %   A =|     A2      Atilde2 |
    %      |         Ap  Atilde3 |
    %
    %%%
    %  WARNING: since A depends on the a priori choice of n and not just on
    %  the data in f, this method might not be very well suited for
    %  determining a good approximation order
    %%%

    if Opts.vf.adaptiveOrder
        A = zeros(p*nSample,n*(p+1));
        Ak      = zeros(nSample,n,p);
        Atildek = zeros(nSample,n,p);

        for kOut = 1:p
            for iRow = 1:nSample
                for jCol = 1:n
                    Ak(iRow,jCol,kOut) = 1/(s0(iRow)-polesvf3(jCol));
                    Atildek(iRow,jCol,kOut) = - f(kOut,1,iRow)/(s0(iRow)-polesvf3(jCol));
                end
            end
        %     A = blkdiag(A,Ak(:,:,kOut));
            cIdxX = (kOut-1)*nSample+1 : kOut*nSample;
            cIdxY = (kOut-1)*n+1 : kOut*n;
            A(cIdxX,cIdxY)       = Ak(:,:,kOut);
            A(cIdxX,n*p+1:end)   = Atildek(:,:,kOut);  
        end

        s = svd(A);
        figure; semilogy(s/s(1),'o-');
        title('Normalized singular values of vector fitting matrix');
        ylabel('$\sigma/\sigma(1)$','interpreter','latex');
        hold on; plot([1,length(s)],Opts.vf.svdTol*[1,1],'r--')

        keyboard
    end
    %%  Run VF
    switch Opts.vf.method
        case 1 %old code by Alessandro
            %MIMO systems have to be fitted columnwise
            f = reshape(f,m*p,nSample);
            AA = []; BB = []; CC = []; DD = [];
            rho=ones(1,nSample);
            for iCol = 1:m
                fm = f((iCol-1)*p+1 : iCol*p,:);
                for iter = 1:Opts.vf.maxiter
                    [SER,polesvf3,rmserr] =vectfit3(fm,s0,polesvf3,rho);
                    fprintf(1,'VF iteration %i, error %e \n',iter,rmserr);
                    if rmserr <= Opts.vf.tol, break, end
                end
                AA = blkdiag(AA,SER.A);
                BB = blkdiag(BB,SER.B);
                CC = [CC, SER.C];
                DD = [DD, SER.D];
            end
            sysr = ss(full(AA),BB,CC,DD);
        case 2 %using code from Serkan
            f = reshape(f,m*p,nSample);
            rho=ones(1,nSample);
            for iter = 1 : Opts.vf.maxiter  % number of vecfit3 steps
                [SS_vectfit3,polesvf3,rmserr,fit,opts]=vectfit3(f,s0,polesvf3,rho);
                fprintf(1,'VF iteration %i, error %e \n',iter,rmserr);
                if rmserr <= Opts.vf.tol, break, end
            end
            sysr = vecfit3_to_ss(SS_vectfit3,polesvf3,n,m,p);    
        case 3 %using Matrix Fitting Toolbox
            %%%
            % WARNING! It assumes symmetry of G(s): Gij(s)=Gji(s)
            %%%
            for iter = 1 : Opts.vf.maxiter  % number of vecfit3 steps
                [SER,rmserr] = VFdriver(f,s0,polesvf3);
                fprintf(1,'VF iteration %i, error %e \n',iter,rmserr);
                if rmserr <= Opts.vf.tol, break, end
            end

            sysr = ss(SER.A,SER.B,SER.C,SER.D);
         case 4 %new code by Alessandro
            %MIMO systems have to be fitted columnwise

            AA = []; BB = []; CC = []; DD = [];
            rho=ones(1,nSample);
            for iCol = 1:m
                fm = f(:,iCol,:);
                fm = reshape(fm,p,nSample);
                for iter = 1:Opts.vf.maxiter
                    [SER,polesvf3,rmserr] =vectfit3(fm,s0,polesvf3,rho);
                    fprintf(1,'VF iteration %i, error %e \n',iter,rmserr);
                    if rmserr <= Opts.vf.tol, break, end
                end
                AA = blkdiag(AA,SER.A);
                BB = blkdiag(BB,SER.B);
                CC = [CC, SER.C];
                DD = [DD, SER.D];
            end
            sysr = ss(full(AA),BB,CC,DD);
    end
    if Opts.plot, 
        figure(fh); plot(complex(polesvf3),'+g'); 
        legend('s0 data','init poles','final poles')
    end
end

%%  Auxiliary
function poles = initializePoles(sys,type,nm,wLim,wReLim)
switch type
    case 'eigs'
        onemore = nm - 2*fix(nm/2) ;
        poles = eigs(sss(sys),nm-onemore,'lm').';
        if onemore, poles = [-abs(imag(poles(1))), poles]; end
    case 'vectfit3'
       %generate initial poles
        if nm > 1
            bet=logspace(-2,log10(wLim(2)),nm/2);
            poles=[];
            for k=1:length(bet)
                alf=-bet(k)*1e-2;
                poles=[poles (alf-1i*bet(k)) (alf+1i*bet(k)) ];
            end
        else %initialize at least one pole
            poles = wLim(2);
        end
    case 'serkan'
        onemore = nm - 2*fix(nm/2) ;
        bet=logspace(max([log10(wLim(1)),-6]),max([log10(wLim(2)),-6]),...
            fix(nm/2)+onemore);
        poles=[];
        if onemore, poles(1) = -bet(1); end ;
        for iIter = 1 + onemore : length(bet)
            alf=-bet(iIter);
            poles=[poles ; (alf-1i*bet(iIter)) ; (alf+1i*bet(iIter)) ];
        end
        poles = poles.';
    case 'serkan+ale'
        onemore = nm - 2*fix(nm/2) ;
        bet=logspace(max([log10(wLim(1)),-2]),max([log10(wLim(2)),-2]),...
            fix(nm/2)+onemore);
        if abs(diff(wReLim))< 1e-3 %all shifts have the same real part
            alf = bet;
        else
            alf=logspace(max([log10(wReLim(1)),-2]),log10(wReLim(2)),fix(nm/2)+onemore);
        end
        poles=[];
        if onemore, poles(1) = -alf(1); end ;
        for iIter = 1 + onemore : length(bet)
%             alf=-bet(iIter);
            poles=[poles ; (-alf(iIter)-1i*bet(iIter)) ; (-alf(iIter)+1i*bet(iIter)) ];
        end
        poles = poles.';
    case 'gershgorin'
end
end

%%  Unused
function sysrvf3 = vecfit3_to_ss(SS_vectfit3,polesvf3,r,m,p)

% INPUTS
% the output of SS_vectfit3 with 1 input m*p outputs

% OUTPUT
% "COMPLEX" state-space form with m inputs, p outputs

Avf3 = SS_vectfit3.A;
Cvec = SS_vectfit3.C;
Phi = {};
for i=1:r
    cc = Cvec(:,i);
    Phi{i} = reshape(cc,p,m);
end

Avf3 = []; Bvf3 = []; Cvf3 = [];

for i=1:r
    Avf3 = blkdiag(Avf3,polesvf3(i)*eye(m));
    Bvf3 = [Bvf3; eye(m)];
    Cvf3 = [Cvf3 Phi{i}];
end

sysrvf3 = ss(Avf3,Bvf3,Cvf3,0*zeros(p,m));
end