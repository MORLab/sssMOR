function [Z,V,Y,xi,err,timings] = RKSMa_lyap(A,E,B,sfreq, tol,maxit,sh,npts,D,opts)
% Rational block Arnoldi method for solving
%
% L(X)=APE'+EPA'+BB'=0
%
%  CALLING SEQUENCE:
%     [y,U,conv,T] = RKSMa_lyap(A,E,B,sfreq, tol,maxit,sh,f,npts,D,opts)
%
% INPUT:
%    E,A         n-x-n matrices
%    B         n-x-m matrix (m << n);
%    sfreq   reduced cale is solved every sfreq iteration
%    tol        stopping tolerance
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
%           opts.trueGres:  estimate the true residual norm of GCARE
%           (0,1), i.e. ||L(X)|| if trueGres=1, otherwise ||E\L(X)/E'||.
% OUTPUT:
%    Z,Y       low-rank solution factors of GCALE solution ZYZ'~X
%    V         orthonormal basis of the rational matrix Krylov subspace
%    xi         generated RKSM poles
%    err      convergence history containing the the relative change 
%              1st row: number of iteration
%              2nd row: current subspace dim
%              3rd row: scaled residual norm 
%    timings   comp. times for different parts
%              1st row: linear solves
%              2nd row: orthogonalization
%              3rd row: small space solution
%              4th row: shift generation
%              5th row: residual norm computation
%
%-----------------------------------------------------------------------
% P. Kuerschner, 2014

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
            xi_cand = logspace(-10,10,npts);%+1i*logspace(-10,10,npts);
    end; end
[N,s] = size(B);
if nargin < 9 || isempty(D), D=eye(s); end
timings=zeros(6,maxit+1);
timings(1,1)=0;
% handling of mass matrix
if isempty(E)
    ML=speye(N); MU=ML;E=ML;noE=1;trueGres=0;
   
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
fAb_conv=0;
ls=tic;

% vE=ML\B;

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
if N == 1, fmR = fm(A)*v; err = 0; return; end;
solver = @(A,b) A\b;    % can use your favorite solver here
zeta = 1*zeros(maxit,1);   % may be optimized
nu1 = s;
ort=tic;
[vE,beta]=qr(vE,0);


% old normalization
% [vE,k01,k02,beta]=defl(vE,nu1,0,1e-12);
% nu1 = length(k01); nu2 = length(k02);
% beta= full(beta(k01,k01));

timings(2,1)=toc(ort);
if adaptive && (strcmp(sh,'conv') || strcmp(sh,'convR') )
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
    case 2
         w=  MU\vE;
         w = A*w(qM,:);
         AV(:,1:nu1) = ML\w(pM,:);
    otherwise
         AV(:,1:nu1) = ML\(A*vE);
end
timings(1,1)= timings(1,1)+toc(ls);


tsh=tic;
R=vE'*AV(:,1:nu1);
theta=eig(R);
timings(4,1)= timings(4,1)+toc(tsh);
H = zeros(2*s*maxit+1,2*s*maxit);
Xi=[];
I = speye(2*maxit+1);O=0*Is;
xi = xi(:).';
fmR=zeros(N,1);
n=1;jc=1; 
curr_start = 1; n = 1; sigma = -1;
while n <= maxit, %n
    prev_start = curr_start;      % .. The first index where columns were added in the previous step.
        prev_end   = prev_start+nu1-1; % .. The  last index where columns were added in the previous step. [curr_start, curr_end]
        real_start = prev_end+1;      % .. The first index where columns (for the real part) will be added in the current step.
        real_end   = real_start+nu1-1; % .. The  last index where columns (for the real part) will be added in the current step. [real_start, real_end].
        imag_start = real_end+1;      % .. The first index where columns for the imaginary part will be added in the current step if the shift is complex.
        imag_end   = imag_start+nu1-1; % .. The  last index where columns for the imaginary part will be added in the current step if the shift is complex. [imag_start, imag_end].

    if adaptive, %shift computation
        sigma_prev = sigma;
        % look for minimum of nodal rational function sn
        tsh=tic;
        theta(real(theta)>0)=-theta(real(theta)>0);
        if strcmp(sh,'convR'), theta=real(theta); end
        thetao=theta;
        %         theta_aug=[theta;-normest(A)];
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
            if n==1, eH=sort(-[s1;nA]); end
            xi_cand=[];
            for j=1:length(eH)-1
                xi_cand=[xi_cand,linspace(eH(j),eH(j+1),500/s)];
            end
            %             xi_cand=sval;
        end
        if n==1
            if strcmp(sh,'conv') || strcmp(sh,'convR'), 
                gs=-s1*ones(s,1); 
            else, gs=inf*ones(s,1); 
            end
        else
            gs=kron([xi(2:end)],ones(1,s))';
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
        timings(4,jc)= timings(4,jc)+toc(tsh);
    end;
    sigma=xi(n+1);
    if( sigma_prev == conj(sigma) )
        % .. Complex conjugate shift, we did this one already, skip
        continue;
    end
    ls=tic;
    
    %     w = MU*(solver(A - xi(n+1)*E,ML*V(:,prev_start:prev_end)));
    switch opts.facM  
%         case 1
%             w1 = MU*(solver(A - xi(n+1)*E,ML*V(:, prev_start:prev_end)));
%         case 2
%             w1(pM,:)=ML*V(:, prev_start:prev_end);
%             w1=solver(A - xi(n+1)*E,w1);
%             w1=MU*w1(pM,:);
%         otherwise
%             w1 = solver(A - xi(n+1)*E,ML*V(:, prev_start:prev_end));
        case 1
        w1 = -MU*(solver(xi(n+1)*E-A,ML*V(:, prev_start:prev_end)));
    case 2
        w1(pM,:)=ML*V(:, prev_start:prev_end);
        w1=-solver(xi(n+1)*E-A,w1);
        w1=MU*w1(pM,:);
    otherwise
        w1 = -solver(xi(n+1)*E-A,ML*V(:, prev_start:prev_end));    
            
            
    end
    timings(1,jc)= timings(1,jc)+toc(ls);
    cplx=0;
    ort=tic;
    wR = real(w1);
    for it=1:2,
        for kk=1:n
            k1=(kk-1)*nu1+1; k2=kk*nu1;
                gamma=V(:,k1:k2)'*wR;
                H(k1:k2,real_start-nu1:real_end-nu1) = H(k1:k2,real_start-nu1:real_end-nu1)+ gamma;
                wR = wR - V(:,k1:k2)*gamma;
            end
    end
%     [H(real_start:real_end, real_start-nu1:real_end-nu1),V(:,real_start:real_end)]=gram_sh_new(wR); 
    [V(:,real_start:real_end),H(real_start:real_end, real_start-nu1:real_end-nu1)]=qr(wR,0); 
    
    if ~isreal(xi(n+1))
        cplx=1;
        wR=imag(w1);
        for it=1:2,
            for kk=1:n+1
                k1=(kk-1)*s+1; k2=kk*s;
                gamma=V(:,k1:k2)'*wR;
                H(k1:k2,imag_start-nu1:imag_end-nu1) = H(k1:k2,imag_start-nu1:imag_end-nu1)+ gamma;
                wR = wR - V(:,k1:k2)*gamma;
            end
        end
%         [H(imag_start:imag_end, imag_start-nu1:imag_end-nu1),V(:,imag_start:imag_end)]=gram_sh_new(wR); 
        [V(:,imag_start:imag_end),H(imag_start:imag_end, imag_start-nu1:imag_end-nu1),]=qr(wR,0); 
        
        curr_start = real_start;  % .. Where the new columns in Z start
        curr_end   = imag_end;    % .. Where the new columns in Z end
        proj_end   = real_end;    % .. Where the columns used for projection end
        
        I(n+1,n+1)=0;
        Xi=blkdiag(Xi,[real(xi(n+1)),imag(xi(n+1));-imag(xi(n+1)),real(xi(n+1))]);
    else
        cplx=0;
        % .. The shift is real
        curr_start = real_start;
        curr_end   = real_end;
        proj_end   = prev_end;
        Xi=blkdiag(Xi,xi(n+1));
       
    end
    timings(2,jc)= timings(2,jc)+toc(ort);
%     invXi=inv(Xi);
    %     V(:,1:n)'*V(:,1:n)
    ssol=tic;
%     AV(:,curr_start:curr_end)=ML\(A*(MU\V(:,curr_start:curr_end)));
    %TODO: make this more efficient
    switch opts.facM  
        case 1 %     newAv
            AV(:,curr_start:curr_end)=ML\(A*(MU\V(:,curr_start:curr_end)));
        case 2
            w=  MU\V(:,curr_start:curr_end);
            w = A*w(qM,:);
            AV(:,curr_start:curr_end) = ML\w(pM,:);
        otherwise
            AV(:,curr_start:curr_end)=ML\(A*(V(:,curr_start:curr_end)));
    end

    g  = V(:, 1:prev_end)' * AV(:, curr_start:curr_end);
    g3 = V(:, curr_start:curr_end)' * AV(:, 1:curr_end);
    R=[R,g; g3];
    
    % analytic formula
%     HT = H(curr_end-nu1+1:curr_end, proj_end-nu1+1:proj_end); % p x p matrix
%     R = (kron(I(1:proj_end,1:proj_end),eye(s))+H(1:proj_end,1:proj_end)*kron(Xi,eye(s)))/H(1:proj_end,1:proj_end)-V(:,1:proj_end)'*AV(:, curr_start+cplx*s:curr_end)*([zeros(nu1, proj_end-nu1) HT] / H(1:proj_end, 1:proj_end));


    
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
    justsolved=0;
    cplxssolve=n>1 && ((xi(n-1)==conj(xi(n))) && mod(n-1,sfreq)==0);
    if (mod(n,sfreq)==0) || cplxssolve%% small scale solution
        justsolved=1;
        ssol=tic;
        roff=[beta;zeros(proj_end-nu1,s)];
        rhs2=roff*D*roff';
        nBf=norm(rhs2);
        %            rhs2=speye(ih*p,p)*beta2*speye(ih*p,p)';
        Y = lyap(R(1:proj_end, 1:proj_end),rhs2);
        timings(3,jc)=timings(3,jc)+toc(ssol);
        
        % computed residual   (exact, in exact arithmetic)
%         tres=tic;
%         g1 = V(:, 1:curr_end-nu1)' * AV(:, curr_end-nu1+1:curr_end);
%         u1 = AV(:, curr_end-nu1+1:curr_end) - V(:, 1:curr_end-nu1) * g1;
%         
%         g2 = V(:, curr_end-nu1+1:curr_end)' * u1; % p x p matrix
%         u2 = u1 - V(:, curr_end-nu1+1:curr_end) * g2;
%         
%         g3 = qr( u2, 0 ); g3 = triu( g3(1:size(g3, 2), :) ); % p x p matrix; note Y is a proj_end x proj_end matrix, proj_end = k*nu1
%         
%         

%         if( proj_end >= 2*nu1 )
%             temp = [zeros(nu1, proj_end-nu1) HT] / H(1:proj_end, 1:proj_end); % nu1 x proj_end matrix
%             R1 = ( [zeros(nu1, proj_end-2*nu1) -imag(xi(n+1))*HT real(xi(n+1))*HT] / H(1:proj_end, 1:proj_end) - g2*temp ) * Y;
%             R2 = -g3*temp*Y;
%         else
%             temp = [zeros(nu1, proj_end-nu1) HT] / H(1:proj_end, 1:proj_end);
%             R1 = ( (real(xi(n+1))*HT) / H(1:proj_end, 1:proj_end) - g2*temp ) * Y;
%             R2 = -g3*temp*Y;
%         end
%         nrmres=norm( [R1; R2] ) / nBf;
        
% %                 u1=newAv(:,end-s+1:end)-V(:,1:j1s-s*(1+cplx))*g(:,1:s);
%         g = V(:,1:j1s-s*(1+cplx))'*newAv(:,end-s+1:end);
%         u1=newAv(:,end-s+1:end)-V(:,1:j1s-s*(1+cplx))*g(:,1:s);
%         %         u1=newAv(:,end-s+1:end)-V(:,1:j1s-s)*g;
%         d=-V(:,1:j1s-s)*(Y*(H(1:js,1:js)'\[sparse(s*(n-1),s);Is])*H(s*n+1:s*(n+1),s*n-s+1:s*n)');
%         U=[-V(:,js1-s:j1s-s)*xi(n+1),  d u1 ];
%         %         d=-V(:,1:j1s)*(Y*(H(1:n*s,1:n*s)'\[sparse(s*(n-1),s);Is])*H(s*n+1:s*(n+1),s*n-s+1:s*n)');
%         %         U=[-V(:,js1:j1s)*xi(n+1),  d u1 ];
%         rr=qr(full(U),0); rr=triu(rr(1:size(rr,2),:));
        
%       nrmres=  norm(rr*sparse([O Is O; Is O Is; O Is O ])*rr',2)/(nBf);%+singE*nrma*nrmx);
%         nrmrestot=[nrmrestot,nrmres];
%         err(1,jc) = proj_end;
%         err(2,jc) = nrmres;
%         timings(5,jc)=timings(5,jc)+toc(tres);
        
%                 abs(eigs(@(x) ML\(A*(ML'\((V(:,1:j1s-s)*Y)*(V(:,1:j1s-s)'*(x)))))+((V(:,1:j1s-s)*Y)*(V(:,1:j1s-s)'*(ML\(A'*(ML'\x)))))+B*(B'*x),N,1,'LM'))/(norm(B'*B,2))
%         abs(eigs(@(x) (A*(MU\((V(:,1:proj_end)*Y)*(V(:,1:proj_end)'*(ML'*x)))))+ML*((V(:,1:proj_end)*Y)*(V(:,1:proj_end)'*(MU\(A'*(x)))))+B*(B'*x),N,1,'LM'))/(norm(B'*B,2))
%         abs(eigs(@(x)(A*(((V(:,1:proj_end)*Y)*(V(:,1:proj_end)'*(E'*x))))+E*((V(:,1:proj_end)*Y)*(V(:,1:proj_end)'*((A'*x))))+B*(B'*x)),N,1,'LM'))/(norm(B'*B,2))
%         fq=uu*([zeros(nu1, proj_end-nu1) HT] / H(1:proj_end, 1:proj_end))*Y;
%       an alternative?? 
        trr=tic;
        HT = H(curr_end-nu1+1:curr_end, proj_end-nu1+1:proj_end); % p x p matrix
        g1 = V(:, 1:curr_end-nu1)' * AV(:, curr_end-nu1+1:curr_end);
        u1 = AV(:, curr_end-nu1+1:curr_end) - V(:, 1:curr_end-nu1) * g1;
        if noE || ~trueGres % standard CALE res.norm (||E\L(X)/E'|| for GCALEs)
            if isreal(xi(n+1)) % real shift
                [~,uu]=qr((V(:, curr_end-nu1+1:curr_end)*xi(n+1)-u1),0);
                qf=(uu*([zeros(nu1, proj_end-nu1) HT] / H(1:proj_end, 1:proj_end)))*Y;
            else % pair of complex shifts
                y1 = [zeros(nu1, proj_end-nu1) HT];
                y2 = (([zeros(nu1, proj_end-2*nu1) -imag(xi(n+1))*HT real(xi(n+1))*HT]));
                [~,uu]=qr([V(:, curr_end-nu1+1:curr_end),u1],0);
                qf=uu*([y2;-y1]/ H(1:proj_end, 1:proj_end))*Y;
            end
            nrmres=norm(qf)/nBf;
        elseif trueGres % true res norm of GCALE (only for E~=I)
            if isreal(xi(n+1)) % real shift
                uf=V(:, curr_end-nu1+1:curr_end)*xi(n+1)-u1;
                hw=(V(:,1:proj_end)*(([zeros(nu1, proj_end-nu1) HT] / H(1:proj_end, 1:proj_end))*Y)');
                
               [~,g3]=qr(ML*[uf,hw],0);
                nrmres=norm(g3*[O,Is;Is,O]*g3')/nBo;
            else % pair of complex shifts
                y1 = [zeros(nu1, proj_end-nu1) HT];
                y2 = (([zeros(nu1, proj_end-2*nu1) -imag(xi(n+1))*HT real(xi(n+1))*HT]));
                hw =V(:,1:proj_end)*(([y2;-y1]/ H(1:proj_end, 1:proj_end))*Y)';
                [~,g3]=qr(ML*[V(:, curr_end-nu1+1:curr_end),u1,hw],0);
                nrmres=norm(g3*kron([O,Is;Is,O],eye(2))*g3')/nBo;
            end
        end
        timings(5,jc)=timings(5,jc)+toc(trr);
        nrmrestot=[nrmrestot,nrmres];
        err(1,jc) = proj_end;
        err(2,jc) = nrmres;
        fprintf('Lyap: \t %d\t %e\t %d\n',n,nrmres,n*s)
%         uf=V(:, curr_end-nu1+1:curr_end)*xi(n+1)-u1;
%         uf=uf;
%         hw=(V(:,1:proj_end)*(([zeros(nu1, proj_end-nu1) HT] / H(1:proj_end, 1:proj_end))*Y)');
%         g3=triu(qr(ML*[uf,hw],0));
%         nrmE=norm(g3*[O,Is;Is,O]*g3')/norm(B'*B,2)
        
%         f=[V(:, curr_end-nu1+1:curr_end),u1]*([y2;-y1]/ H(1:proj_end, 1:proj_end))*Y*V(:,1:proj_end)';
%         
%         norm(f+f')/nBf

        if (nrmres<tol),
            % factored solution
            
             if tolsol % do rank truncation
                [uY,sY,~]=svd(Y); [sY,id]=sort(diag(sY),'descend');
                %                 sY=flipud(sY);
                %                 sY=sY(id);
                uY=uY(:,id);
                is=find(sY>tolsol*sY(1));
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
    if( ~isreal( sigma ) )
            % Two step were done at once.
            n = n + 1;
            curr_start      = imag_start;
            xi(n+1)=conj(xi(n));
    end
    if n >= maxit
        if ~justsolved,
                ssol=tic;
                roff=[beta;zeros(proj_end-nu1,s)];
                rhs2=roff*D*roff';
                nBf=norm(rhs2);
                Y = lyap(R(1:proj_end, 1:proj_end),rhs2);
                timings(3,jc)=timings(3,jc)+toc(ssol);
        end
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
        
        disp(['Max. iteration number of ' num2str(maxit) ' reached.']);
    end;
    tsh=tic;
    % compute Ritz values of order n
    theta=eig(R(1:proj_end,1:proj_end));
    timings(4,jc)=timings(4,jc)+toc(tsh);
    n=n+1; 
end

function sn = evalrat( theta,xi,x )
%EVALRAT Evaluate rational nodal function
%           (x^n + ...)
%   sn(x) = -----------
%            1 + ...
% with zeros theta and poles xi at points x
% in a stable way.
 
n = length(theta);
sn = 0*x + 1;
for j = 1:n,
    % find xi closest to min(abs(sn))
    [dummy,ind] = min(abs(sn));
    xm = x(ind);
    [dummy,ind] = min(abs(xi-xm));
    sn = sn ./ (1-x/xi(ind)+eps);
    xi(ind) = NaN;
    
     % find theta closest to max(abs(sn))
    [dummy,ind] = max(abs(sn));
    tm = x(ind);
    [dummy,ind] = min(abs(theta-tm));
    sn = sn .* (x-theta(ind));
    theta(ind) = NaN;
end;
        
function r=ratfun(x,eH,s)

r=zeros(1,length(x));
for j=1:length(x)
r(j)=abs(prod( (x(j)-s)./(x(j)-eH) ));
end    


function snew=newpolei(eHpoints,eH,s)
snew=zeros(1,length(eHpoints)-1);
for j=1:length(eHpoints)-1
%   snew(j) = fminbnd( @(x) ratfuni(x,eH,s), eHpoints(j),eHpoints(j+1),optimset('TolX',1e-3));
     sval=linspace(eHpoints(j),eHpoints(j+1),200);
     [~,jx] = max (abs(ratfun(sval,eH,s)));
     snew(j)=sval(jx);
end
[~,jx]=max(abs(ratfun(snew,eH,s)));
snew=snew(jx);
return


% function r=ratfun(x,eH,s)
% 
% r=zeros(1,length(x));
% for j=1:length(x)
% r(j)=abs(prod( (x(j)-eH)./(x(j)-s) ));
%