function [Z,V,Y,xi,err,timings] = RKSMat_lyap(A,E,B,sfreq, tol,maxit,sh,npts,D,opts)
% Adaptive tangential rational Arnoldi method for solving
%
% L(X)=AXE'+EXA'+BB'=0
%
%  CALLING SEQUENCE:
%     [y,U,conv,T] = RKSMat_lyap(A,E,B,sfreq, tol,maxit,sh,f,npts,D,opts)
%
% INPUT:
%    E,A         n-x-n matrices
%    B         n-x-m matrix (m << n);
%    sfreq   reduced cale is solved every sfreq iteration
%    tol      tolerances: tol(1) for rel.res (obligatory), 
%               tol(2) for rank.trunc., tol(2)=0 for no rank trunc 
%    maxit    maximal number iterations
%    sh      shifts (poles) for RKSM (if empty-> adaptive shifts)
%             can also be a string:
%           'im' -- choose along imaginary axis
%           'real' -- choose along real interval
%           'conv' -- choose from convex hull of Ritz values
%           'convR' -- real version of 'conv'
%    npts     numbers of test points for adaptive shift selection
%   D     - middle factor of rhs (for indefinite problems)
%   opts    = structure containing info on mass matrix handling 
%           opts.ML,opts.MU --cholesky or LU factors of E=ML*MU, if E=E^T,
%           MU=ML'
%           opts.facM - compute factors of E if empty ML,MU, 
%           compute Chol. / LU factors of E (facM=1)
%           compute Chol. / LU factors of E with pivoting (facM=2,default)
%           or directly work with E (facM=0)
%           opts.trueGres:  estimate the true residual norm of GCALE
%           (0,1), i.e. ||L(X)|| if trueGres=1, otherwise ||E\L(X)/E'||.
% OUTPUT:
%    Z,Y       low-rank solution factors of GCALE solution ZYZ'~X
%               (Y=I if tol(2)~=0)
%    V         orthonormal basis of the tangential rational Krylov subspace
%    xi         generated RKSM poles
%    err      convergence history containing the the relative change 
%              1st row: current subspace dim
%              2nd row: scaled residual norm 
%    timings   comp. times for different parts
%              1st row: linear solves
%              2nd row: orthogonalization
%              3rd row: small space solution
%              4th row: shift and direction generation
%              5th row: residual norm computation
%
% [1] Druskin/Simoncini/Zaslavsky, Adaptive Tangential Interpolation 
% in Rational {K}rylov Subspaces for {MIMO} Dynamical Systems,
% SIAMMatrix, 35(2):476–498,2014
%
% [2] Druskin/Simoncini, Adaptive rational {K}rylov subspaces for 
% large-scale dynamical systems, SCL, 60(8):546–560,2011
% 
%-----------------------------------------------------------------------
% P. Kuerschner, 2015

if nargin < 10 || isempty(opts),
    opts.facM=1; opts.ML=[]; opts.MU=[]; trueGres=0;
else
    if isfield(opts,'trueGres'), trueGres=opts.trueGres;
    else
       trueGres=0; 
    end
end

if nargin < 8 || isempty(npts), npts = 1601; end
% if nargin < 8 || isempty(exact), exact = []; end;
if nargin < 7 || isempty(sh),
    sh='conv'; adaptive = 1;
end

if length(tol)>1
    tolsol=tol(2); tol=tol(1);
else
    tolsol=1e-12;
end

if (norm(A-A',1)<1e-14), symA=1; else symA=0;end

if isnumeric(sh)
    xi=sh;adaptive = 0;
else
    
    switch(sh)
        case 'conv'
            adaptive = 1; xi=[];
        case 'convR'
            adaptive = 1; xi=[];
        otherwise sh='conv';
    end; end
[N,s] = size(B);
if nargin < 9 || isempty(D), D=eye(s); end
timings=zeros(5,maxit+1);
timings(1,1)=0;
% handling of mass matrix
pM=1:N; qM=pM;
if isempty(E)
    ML=speye(N); MU=ML;E=ML; noE=1;trueGres=0;
else
    if trueGres,  nBo=max(abs(eig(B'*B*D))); end
    noE=0;
    if (norm(E-E',1)<1e-14), symM=1; else symM=0;end
    if ~isempty(opts.ML)
        ML=opts.ML;
        if ~isempty(opts.MU), MU=opts.MU;
        else
            if symM, MU=ML'; else MU=speye(N); end
        end
    else
        switch opts.facM
            case 1 % use LU / Cholesky factorization
                ls=tic;
                if symM
                    ML=chol(E,'lower'); MU=ML';
                else
                    [ML,MU]=lu(E);
                end
                timings(1,1)=toc(ls);
            case 2 % use LU / Cholesky factorization with pivoting
                ls=tic;
                if symM
                    [ML,~,pM]=chol(E,'lower','vector'); MU=ML';
                    qM(pM)=1:N;
                else
                    [ML,MU,pM,qM]=lu(E,'vector');
                end
                timings(1,1)=toc(ls);
            otherwise %no factor.
                %         else
                ML=E;
                MU=speye(N);
        end
    end
end
Is=speye(s);nrmrestot=[];
ls=tic;
switch opts.facM,
    case 1
        vE=ML\B;
    case 2
        vE=ML\B(pM,:);
    otherwise
        vE=ML\B;
end
timings(1,1)= timings(1,1)+toc(ls);
if nargin < 4 || isempty(maxit), maxit = min(N,50); end;
if nargin < 3 || isempty(tol), tol = 1e-8; end;
solver = @(A,b) A\b;    % can use your favorite solver here
zeta = 1*zeros(maxit,1);   % may be optimized
nu1 = s;sdir=nu1;
ort=tic;
vEo=vE;
[vE,beta]=qr(vE,0);
roff=beta;
timings(2,1)=toc(ort);
if adaptive && s==1 && (strcmp(sh,'conv') || strcmp(sh,'convR') ) 
    % initial shift(s) for tangential selection
    tsh=tic;
    eops.tol=1e-1;
%           try
%             s1=eigs(A,E,1,'LR');
%           catch
              try
                   s1=real(eigs(A,E,1,'SM',eops));
              catch
                   s1=-0.1;
              end
%           end
   try
        nA=-abs(eigs(A,E,1,'LM',eops));
   catch
        nA=-normest(A);
   end
    
    %         s1=min(abs(ss));
    %         nA=max(abs(ss));
    %         normest(A);
    timings(4,1)= timings(4,1)+toc(tsh);
end

V = zeros(N,s*2*maxit+1);
AV=V;
V(:,1:nu1) = vE;
ls=tic;
% AV(:,1:nu1)=ML\(A*(MU\vE));
switch opts.facM  
    case 1
         AV(:,1:nu1) = ML\(A*(MU\vE)); 
%          AB = ML\(A*(MU\vEo)); 
    case 2
         w=  MU\vE;
         w = A*w(qM,:);
         AV(:,1:nu1) = ML\w(pM,:);
%          w=  MU\vEo;
%          w = A*w(qM,:);
%          AB = ML\w(pM,:);
         
    otherwise
         AV(:,1:nu1) = ML\(A*vE);
%          AB = ML\(A*vEo);
         
end
timings(1,1)= timings(1,1)+toc(ls);
% initialize data for  tangential direction and residual computations
[~,RC]=qr([AV(:,1:nu1)*beta,vE],0);
tsh=tic;
R=vE'*AV(:,1:nu1);
[UR,theta]=eig(R); theta=diag(theta);
timings(4,1)= timings(4,1)+toc(tsh);
H = zeros(2*s*maxit+1,2*s*maxit);
WR=RC*[eye(s);-R];
Xi=[];Dt=[];
Xi=zeros(s);
I = speye(N);O=0*Is;
xi = xi(:).';
fmR=zeros(N,1);
n=1;jc=1; 
proj_end=s;
blk(1)=1;blk(2)=nu1; % save block sizes
curr_start = 1; n = 1; sigma = -1;
while n <= maxit, %n
    % polynomial Arnoldi step for computing R = V'*A*V
%     Av=AV(:,prev_start:prev_end);
%     w = ML*Av; % - V(:,prev_start:prev_end)*diag(zeta(n)*ones(s,1));
if adaptive, %shift + direction computation
    sigma_prev = sigma;
    % look for minimum of nodal rational function sn
    tsh=tic;
    theta(real(theta)>0)=-theta(real(theta)>0);
    if strcmp(sh,'convR'), theta=real(theta); end
    thetao=theta;
    %         theta_aug=[theta;-normest(A)];
    
    if s==1 % SISO case ... only shifts via argmax(rat.fun)
        if strcmp(sh,'conv') || strcmp(sh,'convR') %&& length(theta)>2 %
            %           convex hull test set ala
            %           Simoncini/Druskin
            if any(imag(theta)) && length(theta)>2
                theta=[theta;nA];
                ch=convhull(real(theta),imag(theta));
                eH=-theta(ch);
                ieH=length(eH); missing=n*s-ieH;
                while missing>0, % include enough points from the border
                    neweH=(eH(1:ieH-1)+eH(2:ieH))/2;
                    eH=[eH;neweH];
                    missing=n*s-length(eH);
                end
                
            else
                eH=sort([-real(theta);-s1;-nA]);
            end
            if n==1, eH=sort(-[s1;nA]);
                if strcmp(sh,'conv') || strcmp(sh,'convR'),
                    gs=-s1*ones(s,1);
                else, gs=inf*ones(s,1); end
            else
                gs=kron([xi(2:end)],ones(1,s))';
                
            end
            xi_cand=[]; % candidate values
            if strcmp(sh,'conv')
                for j=1:length(eH)-1
                    xi_cand=[xi_cand,linspace(eH(j),eH(j+1),npts/(length(eH)))];
                end
            elseif strcmp(sh,'convR')
                xi_cand=linspace(min(abs(eH)),max(abs(eH)),npts);
            end
            sn=ratfun(xi_cand,thetao,gs);
            %         xi(n+1)=newpolei(eH',thetao,gs);
            [~,jx]=max(abs(sn));
            if real(xi_cand(jx))<0, xi_cand(jx)=-xi_cand(jx); end
            if abs(imag(xi_cand(jx))) / abs( xi_cand(jx))<1e-8
                xi(n+1)=real(xi_cand(jx));
            else
                xi(n+1)=xi_cand(jx);
            end
        end
        dir=1;
        %             xi_cand=sval;
    else % adaptive shift and tangential directions generation
        xi_cand=unique(-thetao); % candidate values
%         funshift=@(x)norm(WR*((R(1:proj_end,1:proj_end)-x*eye(proj_end))\roff));
        
         WRU=WR*UR;
        roffU=UR\roff;
        SR=spdiags(theta,0,proj_end,proj_end);
        
        funshift=@(x)norm(WRU*((SR-x*speye(proj_end))\roffU));
%           svdopts.tol=1e-4;
%           svdopts.maxit=10;
%           svdopts.v0=ones(s,1);
%           svdopts.issym=1;
%         funshift=@(x)sqrt(abs(eigs(@(v) resfun(v,x,R(1:proj_end,1:proj_end),roff,WR),s,1,'LM',svdopts)));
        rms=[];
        % find shift
        for hh=1:length(xi_cand),
            ff=funshift(xi_cand(hh));
            rms(hh)=ff(1);
        end
        [~,jx]=max(rms);
        if real(xi_cand(jx))<0, xi_cand(jx)=-xi_cand(jx); end
        if abs(imag(xi_cand(jx))) / abs( xi_cand(jx))<1e-8
            xi(n+1)=real(xi_cand(jx));
        else
            xi(n+1)=xi_cand(jx);
        end
        % tangential direction(s)
        [~,sings,Vs]=svd(WR*((R(1:proj_end,1:proj_end)-xi(n+1)*eye(proj_end))\roff),0);
        [~,idx]=find(sings/sings(1)>.1);  
        dir=Vs(:,idx);
        if isreal(xi(n+1))
            Dt=[Dt,dir];
        else
            Dt=[Dt,real(dir),imag(dir)];
        end
    end
    timings(4,jc)= timings(4,jc)+toc(tsh);
end    
    sigma=xi(n+1);
    if( ~strcmp(sh,'convR') && (sigma_prev == conj(sigma) ))
        % .. Complex conjugate shift, we did this one already, skip
        continue;
    end
    nu_old=nu1;
%     
    prev_start = curr_start;      % .. The first index where columns were added in the previous step.
        prev_end   = prev_start+nu_old-1; % .. The  last index where columns were added in the previous step. [curr_start, curr_end]
        real_start = prev_end+1;      % .. The first index where columns (for the real part) will be added in the current step.
        nu1=size(dir,2);
        real_end   = real_start+nu1-1; % .. The  last index where columns (for the real part) will be added in the current step. [real_start, real_end].
        imag_start = real_end+1;      % .. The first index where columns for the imaginary part will be added in the current step if the shift is complex.
        imag_end   = imag_start+nu1-1; 
    ls=tic;
    w1=vEo;
    switch opts.facM  %% explicit for first steps (????)
        case 1
            w1 = -MU*(solver(xi(n+1)*E-A,ML*w1*dir));
        case 2
            w1=w1*dir;
            w1(pM,:)=ML*w1;
            w1=-solver(xi(n+1)*E-A,w1);
            w1=MU*w1(pM,:);
        otherwise
            w1 = -solver(xi(n+1)*E-A,ML*w1*dir);
    end
    timings(1,jc)= timings(1,jc)+toc(ls);
    cplx=0;
    ort=tic;
    wR = real(w1);
    k1=1;k2=proj_end;
    blk(n+2)=blk(n+1)+nu1;
    for it=1:2,
        for kk=1:n+1
            k1=blk(kk); k2=blk(kk+1);
            gamma=V(:,k1:k2)'*wR;
            H(k1:k2,real_start-s:real_end-s) = H(k1:k2,real_start-s:real_end-s)+ gamma;
            wR = wR - V(:,k1:k2)*gamma;
        end
    end
    [V(:,real_start:real_end),H(real_start:real_end, real_start-s:real_end-s)]=qr(wR,0); 
%     [H(real_start:real_end, real_start-s:real_end-s),V(:,real_start:real_end)]=gram_sh_new(wR); 
    
    if ~isreal(xi(n+1))
        cplx=1;
%         Dt=[Dt,conj(dir)];
        wR=imag(w1);
        blk(n+3)=blk(n+2)+nu1;
        for it=1:2,
            for kk=1:n+2
                k1=blk(kk); k2=blk(kk+1);
                gamma=V(:,k1:k2)'*wR;
                H(k1:k2,imag_start-s:imag_end-s) = H(k1:k2,imag_start-s:imag_end-s)+ gamma;
                wR = wR - V(:,k1:k2)*gamma;
            end
        end
        [V(:,imag_start:imag_end),H(imag_start:imag_end, imag_start-s:imag_end-s),]=qr(wR,0); 
%         [H(imag_start:imag_end, imag_start-s:imag_end-s),V(:,imag_start:imag_end)]=gram_sh_new(wR); 
        
        curr_start = real_start;  % .. Where the new columns in Z start
        curr_end   = imag_end;    % .. Where the new columns in Z end
        proj_end   = imag_end;    % .. Where the columns used for projection end
        Xi=blkdiag(Xi,kron([real(xi(n+1)),imag(xi(n+1));-imag(xi(n+1)),real(xi(n+1))],eye(nu1))); 
    else
        cplx=0;
        % .. The shift is real
       curr_start = real_start;
        curr_end   = real_end;
        proj_end   = real_end;
       Xi=blkdiag(Xi,xi(n+1)*eye(nu1));
    end
    timings(2,jc)= timings(2,jc)+toc(ort);

    ssol=tic;
%     AV(:,curr_start:curr_end)=ML\(A*(MU\V(:,curr_start:curr_end)));
    %TODO: make this more efficient
    switch opts.facM  
        case 1
            AV(:,curr_start:curr_end)=ML\(A*(MU\V(:,curr_start:curr_end)));
        case 2
            w=  MU\V(:,curr_start:curr_end);
            w = A*w(qM,:);
            AV(:,curr_start:curr_end) = ML\w(pM,:);
        otherwise
            AV(:,curr_start:curr_end)=ML\(A*(V(:,curr_start:curr_end)));
    end
%     newAv=AV(:,j1s-s*(1+cplx)+1:j1s);
    g  = V(:, 1:prev_end)' * AV(:, curr_start:curr_end);
    g3 = V(:, curr_start:curr_end)' * AV(:, 1:curr_end);
    R=[R,g; g3];
%     AV(:,j1s-s*(1+cplx)+1:j1s)=ML\(A*(MU\V(:,j1s-s*(1+cplx)+1:j1s)));
%     newAv=AV(:,j1s-s*(1+cplx)+1:j1s);
%     g=V(:,1:j1s-s*(1+cplx))'*newAv;%(:,end-s+1:end);
%     %     g=V(:,1:j1s-s*(1+cplx))'*AV(:,j1s-s*(1+cplx)+1:j1s-s);
%     %     g2 = V(:,js1-s:j1s-s)'*(E\(A*V(:,1:js-s)));
%     g3 = V(:,js1-s*cplx:j1s)'*AV(:,1:j1s);
%     %     g3=  V(:,js1-s:j1s-s)'*AV(:,1:j1s-s);
%     %     R=V(:,1:j1s-s)'*(E\(A*V(:,1:j1s-s)));
%     R=[R,g; g3];
    %     R=V(:,1:j1s-s)'*AV(:,1:j1s-s);
    %        R=V(:,1:j1s)'*AV(:,1:j1s);
    %     R=V(:,1:j1s-s)'*(ML\(A*(ML'\V(:,1:j1s-s))));
    %     R=V(:,1:n*s-s*cplx)'*(E\(A*V(:,1:n*s-s*cplx)));%     invXi=inv(Xi);
    %     K = H(1:n*s-s*cplx,1:n*s-s*cplx)*blkdiag(zeros(s),invXi(1:n*s-s*cplx,1:n*s-s*cplx)) + eye(n*s-s*cplx);
    %     R = (H(1:n*s-s*cplx,1:n*s-s*cplx) + diag(kron(zeta(1:n-cplx),ones(s,1))))/K;
    timings(3,jc)=timings(3,jc)+toc(ssol);
    if (mod(n,sfreq)==0)
        ssol=tic;
        roff=[beta;zeros(proj_end-s,s)];
        rhs2=roff*D*roff';
        if trueGres, nBf=nBo; else nBf=norm(rhs2); end
        %            rhs2=speye(ih*p,p)*beta2*speye(ih*p,p)';
%         if n==34,
%             disp('now')
%         end
        Y = lyap(R(1:proj_end, 1:proj_end),rhs2);
        timings(3,jc)=timings(3,jc)+toc(ssol);
        
        % computed residual   (exact, in exact arithmetic)
        tres=tic;
        [~,RC]=qr(ML*[AV(:,1:s)*beta,V(:,1:proj_end)],0);
        
        VB=V(:,1:proj_end)'*B;
        
        Rb=[roff,H(1:proj_end,1:proj_end-s)];
        sdir=sdir+nu1+(~isreal(xi(n+1)))*nu1;
        W=RC*([(blkdiag(eye(s),roff*Dt)+[zeros(s,size(Rb,2));Rb*Xi])/Rb,[zeros(s,sdir);eye(sdir)]]);
        WR=W*[eye(sdir);-R(1:proj_end,1:proj_end)];
        nrmres=norm(W*[zeros(size(Y)),Y;Y,rhs2]*W')/nBf;
        nrmrestot=[nrmrestot,nrmres];
        err(1,jc) = proj_end;
        err(2,jc) = nrmres;
        timings(5,jc)=timings(5,jc)+toc(tres);
        fprintf('nres: \t %d\t %e\t %d\t%d\n',n,nrmres,proj_end,nu1)
%                 abs(eigs(@(x) ML\(A*(ML'\((V(:,1:j1s-s)*Y)*(V(:,1:j1s-s)'*(x)))))+((V(:,1:j1s-s)*Y)*(V(:,1:j1s-s)'*(ML\(A'*(ML'\x)))))+B*(B'*x),N,1,'LM'))/(norm(B'*B,2))
%         abs(eigs(@(x) ML\(A*(MU'\((V(:,1:proj_end)*Y)*(V(:,1:proj_end)'*(x)))))+((V(:,1:proj_end)*Y)*(V(:,1:proj_end)'*(MU\(A'*(ML'\x)))))+vEo*(vEo'*x),N,1,'LM'))/nBf
%         switch opts.facM  
%             case 1
%                         Z=MU\V(:,1:size(Y,1));
%             case 2
%                         Z(pM,1:size(Y,1)) = MU\V(:,1:size(Y,1));
%             otherwise
%                         Z=V(:,1:size(Y,1));
%         end
%         resN_gen=abs(eigs(@(x) (A*(Z*Y*(Z'*(E'*x))))+(E*(Z*Y*(Z'*(A'*x))))+B*(B'*x),N,1,'LM'))/norm(B'*B,2)
        
        
%         [ nrmres*nBf; abs(eigs(@(x) (A*(MU\((V(:,1:proj_end)*Y)*(V(:,1:proj_end)'*(ML'*x)))))+ML*((V(:,1:proj_end)*Y)*(V(:,1:proj_end)'*(MU\(A'*(x)))))+B*(B'*x),N,1,'LM'))]/(norm(B'*B,2))
        if( ~isreal( sigma ) )
            % Two step were done at once.
            n = n + 1;
            curr_start      = imag_start;
            xi(n+1)=conj(xi(n));
        end
        if (nrmres<tol),
            % factored solution
             if tolsol % do rank truncation
                [uY,sY,~]=svd(Y); [sY,id]=sort(diag(sY),'descend');
                %                 sY=flipud(sY);
                %                 sY=sY(id);
                uY=uY(:,id);
                is=find(sY>eval(tolsol)*sY(1));
                Y0 = uY(:,is)*diag(sqrt(sY(is)));
                Z = V(:,1:size(Y0,1))*Y0;
                ls=tic;
                switch opts.facM  %% explicit for first steps (????)
                    case 1
                        Z=MU\Z;
                    case 2
                        Z(pM,:) = MU\Z;
                end
                timings(1,jc)= timings(1,jc)+toc(ls);
                Y=eye(size(Z,2));
            else
                ls=tic;
                switch opts.facM  %% explicit for first steps (????)
                    case 1
                        Z=MU\V(:,1:size(Y,1));
                    case 2
                        Z(pM,:) = MU\V(:,1:size(Y,1));
                        %                             Z = Z(qM,:);
                    otherwise
                        Z=V(:,1:size(Y,1));
                end
                timings(1,jc)= timings(1,jc)+toc(ls);
            end
                V=V(:,1:proj_end);
            break,
            
        end
        jc=jc+1;
    end;
    if n >= maxit
        V=V(:,1:proj_end);
        switch opts.facM  %% explicit for first steps (????)
            case 1
                Z=MU\V(:,1:proj_end);
            case 2
                Z(pM,:) = MU\V(:,1:proj_end);
                %                             Z = Z(qM,:);
            otherwise
                Z=V(:,1:proj_end);
        end
        V=V(:,1:proj_end);
        disp(['Max. iteration number of ' num2str(maxit) ' reached.']);
    end;
    tsh=tic;
    % compute Ritz values of order n
   [UR,theta]=eig(R(1:proj_end,1:proj_end)); theta=diag(theta);
    timings(4,jc)=timings(4,jc)+toc(tsh);
    n=n+1; 
end
timings=timings(:,1:jc);

function r=ratfun(x,eH,s)

r=zeros(1,length(x));
for j=1:length(x)
r(j)=abs(prod( (x(j)-s)./(x(j)-eH) ));
end    



function z=resfun(v,x,R,b,W) % function for res.norm estimation for the use in eigs
    z=W*((R-x*eye(size(R)))\b)*v;
    z=z'*z;
