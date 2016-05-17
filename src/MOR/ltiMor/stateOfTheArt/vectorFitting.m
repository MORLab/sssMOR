function sysr  = vectorFitting(sys,n,s0,Opts)

figure('Name','Sampling frequencies for VF');
plot(complex(s0),'o');xlabel('real');ylabel('imag')

% take only one complex conjugate partner
s0 = cplxpair(s0); idx = find(imag(s0)); s0(idx(1:2:end)) = [];

m = size(sys.b,2); p = size(sys.c,1);
%                 if m>1, n = round(n/m);end %avoid blowing-up for MIMO
if mod(n,2) ~= 0, n = n-1; end   %make even

nSample = length(s0);

f = freqresp(sss(sys),s0); 
f = reshape(f,m*p,nSample);

polesvf3 = initializePoles(sys,Opts.vf.poles,n,Opts.vf.wLims);
hold on; plot(complex(polesvf3),'rx');


if 1 %old code by Alessandro
    %MIMO systems have to be fitted columnwise
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
else %using code from Serkan
    rho=ones(1,nSample);
    for iter = 1 : Opts.vf.maxiter  % number of vecfit3 steps
        [SS_vectfit3,polesvf3,rmserr,fit,opts]=vectfit3(f,s0,polesvf3,rho);
        fprintf(1,'VF iteration %i, error %e \n',iter,rmserr);
        if rmserr <= Opts.vf.tol, break, end
    end
    sysr = vecfit3_to_ss(SS_vectfit3,polesvf3,n,m,p);
end
end

function poles = initializePoles(sys,type,nm,wLim)
switch type
    case 'eigs'
        poles = eigs(sss(sys),nm,'sm').';
    case 'vectfit3'
        try
            wMax = abs(imag(eigs(sss(sys),1,'li')));
        catch
            wMax = 1e3;
        end
        
        %generate initial poles
        if nm > 1
            bet=logspace(-2,log10(wMax),nm/2);
            poles=[];
            for k=1:length(bet)
                alf=-bet(k)*1e-2;
                poles=[poles (alf-1i*bet(k)) (alf+1i*bet(k)) ];
            end
        else %initialize at least one pole
            poles = wMax;
        end
    case 'serkan'
        onemore = nm - 2*fix(nm/2) ;
        bet=logspace(log10(wLim(1)),log10(wLim(2)),fix(nm/2)+onemore);
        poles=[];
        if onemore, poles0(1) = -bet(1); end ;
        for iIter = 1 + onemore : length(bet)
            alf=-bet(iIter);
            poles=[poles ; (alf-1i*bet(iIter)) ; (alf+1i*bet(iIter)) ];
        end
        poles = poles.';
    case 'gershgorin'
end
end

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