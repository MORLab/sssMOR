function [SER,rmserr,bigHfit,opts]=VFdriver(bigH,s,poles,opts); 

%   [SER,rmserr,bigHfit]=poleresiduefit(bigH,s,poles); 
%   [SER,rmserr,bigHfit]=poleresiduefit(bigH,s,poles,opts); 
%     =============================================================================
%     =   Routine: poleresiduefit.m                                               =
%     =   Version 1.0                                                             =
%     =   Last revised: 27.08.2009                                                = 
%     =   Programmed by: Bjorn Gustavsen,                                         =
%     =   SINTEF Energy Research, N-7465 Trondheim, NORWAY                        =
%     =============================================================================
%
% PURPOSE : Calculate a rational model of a symmetric matrix of frequency domain 
%           data (s,H(s)), on pole-residue form, and on state space form
%                     
%              H(s)=SUM(Rm/(s-am)) +D +s*E    %pole-residue  
%                        m 
%
%              H(s)=C*(s*I-A)^(-1)*B +D +s*E  %state-space  
%
% APPROACH: The elements of the upper triangle of H are stacked into a single vector 
%           that is fitted using a Fast implementation of the Relaxed version of Vector Fitting
%           (FRVF) as implemented in vectfit3.m. All matrix elements are fitted with a common 
%           pole set.
% ==========================================================================            
% INPUT :
%
% bigH(Nc,Nc,Ns) : matrix function to be fitted. 
%                 Nc : dimension of H(s)
%                 Ns : number of frequency samples 
%
% 
% s(1,Ns)    : vector of frequency points [rad/sec] 
%
% poles(1,N) :  vector of Initial poles. Specify "poles=[]" if automated 
%               pole generation is to be used.
% ==========================================================================     
% opts. : structure for 
%           -automated generation of initial poles
%           -plotting of results
%           -overriding default (Def) settings.
%
% Automated generation of initial poles:
%   opts.N         : Fitting order 
%   opts.polestype : valid options: 
%                   'lincmplx' : linearly spaced, complex conjugate pairs
%                   'logcmplx' : logarithmically spaced, complex conjugate pairs
%   opts.nu        : ratio between imag.part and real part. (Def=0.001). 
%
% opts.Niter1 :  N.o. VF iterations: Fitting sum of matrix elements (Def=4)
%               (obtains improved initial poles with little computational effort)
% opts.Niter2 :  N.o. VF iterations: Fitting entire matrix (Def=4)
%
% opts.weight(Nc,Nc,Ns): The rows in the system matrix are weighted using this array       
%              Specify "opts.weight=[]" if you wish to use the opts.weightparam option 
%
% opts.weightparam (Def=1)
%            =1 --> weight=1 for all elements in Least Sq. problem, at all freq. 
%            =2 --> weight(s)=1/abs(Hij(s))      ; indvidual element weight
%            =3 --> weight(s)=1/sqrt(abs(Hij(s))); indvidual element weight
%            =4 --> weight(s)=1/norm(H(s))       ; common weight for all matrix elements
%            =5 --> weight(s)=1/sqrt(norm(H(s)) ; common weight for all matrix elements
%
% opts.asymp=1 --> D=0,  E=0  ('Strictly proper') 
% opts.asymp=2 --> D~=0, E=0  ('Proper')                           (Def=2)
% opts.asymp=3 --> D~=0, E~=0 ('Improper') 
% 
% opts.stable=1 --> poles are 'flipped' into the left half plane   (Def=1)
% opts.relaxed=1 --> Use VF with relaxed non-triviality constraint (Def=1)
%
% opts.plot=1 -->  Plotting results to screen:  (Def=1)
%                  figure(1): magnitude functions
%                    cyan trace    : H
%                    magenta trace : (H)fit 
%                    green trace   : H - (H)fit
%
% opts.logx=0 --> Plotting using linear abscissa axis      (Def=0)
% opts.logx=1 --> Plotting using logarithmic absissa axis             
%
% opts.logy=0 --> Plotting using linear ordinate axis (Def=1) 
% opts.logy=1 --> Plotting using logarithmic ordinate axis
% opts.errplot=1--> Will include fitting error in the plot (Def=1)
%
% opts.phaseplot=1--> Will make separeate plot of phase angles (figure(2)). (Def=1)
%
% opts.screen         %=1--> Will write info to screen during fitting (Def=1) 
%
% opts.cmplx_ss=1 --> Will generate complex state space model with diagonal A  (Def=1)
% opts.cmplx_ss=0 --> Will generate real-only state space model with block-diagonal A
%
%
% opts.remove_HFpoles=1 --> Will throw out high frequency out-of-band poles   (Def=0)
% opts.factor_HF          %Will throw out poles above opts.factor_HF*s(end). (Def=1.1)
%                         (Used when opts.remove_HFpoles=1).  
% opts.passive_DE    =1 --> Enforce positive realness for D,E during fitting (Def=0)
%
% ==========================================================================     
% OUTPUT :
% 
%       H(s)=C*(s*I-A)^(-1)*B +D +s*E    %state space model
%                     
%           =SUM(Rm/(s-am)) +D +s*E      %pole-residue model 
%             m 
%
% SER.A (N,N)    : matrix A (diagonal)
% SER.B(Nc*N,Nc) : matrix B
% SER.C(Nc,Nc*N) : matrix C
% SER.D(Nc,Nc)   : matrix D (is produced if asymp=2 or 3)
% SER.E(Nc,Nc)   : matrix E (is produced if asymp=3)
%
% SER.R(Nc,Nc,N) : Residue matrices
% SER.poles (N,1): Poles for pole-reside model
%
% rmserr        : root-mean-square error of approximation for H(s)

%***************************************************************************** 
% NOTE: The use of this program is limited to NON-COMMERCIAL usage only.
% If the program code (or a modified version) is used in a scientific work, 
% then reference should be made to the following:  
%
% [1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency       
%     domain responses by Vector Fitting", IEEE Trans. Power Delivery,        
%     vol. 14, no. 3, pp. 1052-1061, July 1999.                                
%
% [2] B. Gustavsen, "Improving the pole relocating properties of vector
%     fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592,
%     July 2006.   
%
% [3] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter,
%     "Macromodeling of Multiport Systems Using a Fast Implementation of
%     the Vector Fitting Method", IEEE Microwave and Wireless Components 
%     Letters, vol. 18, no. 6, pp. 383-385, June 2008.
%**************************************************************************
%******
%
% The relaxation (opts.relaxed) is described in                                  
%       B. Gustavsen, "Improving the pole relocating properties of vector       
%       fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592,    
%       July 2006.                                                                   
%
%********************************************************************************
% Bug-fixes :
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Default settings:
def.N=[];
def.poletype='lincmplx';
def.nu=1e-3;
def.Niter1=4;
def.Niter2=4;
def.weight=[];
def.weightparam=1;
def.asymp=1; %changed from original
def.stable=1;
def.relaxed=1;
def.plot=0;%changed from original
def.logx=1;
def.logy=1;
def.errplot=0; %changed from original
def.phaseplot=1;
def.screen=0; %changed from original
def.cmplx_ss=0; %changed from original
def.remove_HFpoles=0;
def.factor_HF=1.1;
def.passive_DE=0;
def.passive_DE_TOLD=1e-6;
def.passive_DE_TOLE=1e-16;


if nargin<4
  opts=def;
else 
  %Merge default values into opts  
  A=fieldnames(def);    
  for m=1:length(A)
    if ~isfield(opts,A(m))
      dum=char(A(m)); dum2=getfield(def,dum); opts=setfield(opts,dum,dum2);
    end
  end  
end 

VF.asymp=opts.asymp; 
if opts.stable==1,   VF.stable=1; else, VF.stable=0; end;  
if opts.relaxed==1,  VF.relax=1; else, VF.relax=0; end   
if opts.plot==1,     VF.spy2=1; else, VF.spy2=0; end  
if opts.logx==1,     VF.logx=1; else, VF.logx=0; end  
if opts.logy==1,     VF.logy=1; else, VF.logy=0; end  
if opts.errplot==1,  VF.errplot=1; else, VF.errplot=0; end  
if opts.phaseplot==1,VF.phaseplot=1; else, VF.phaseplot=0; end 
if opts.cmplx_ss==1, VF.cmplx_ss=1; else,VF.cmplx_ss=0; end  


Niter1=opts.Niter1;  
Niter2=opts.Niter2; 
weightparam=opts.weightparam;        
remove_HFpoles=opts.remove_HFpoles; 
factor_HF=opts.factor_HF;  
passive_DE=opts.passive_DE;
passive_DE_TOLD=opts.passive_DE_TOLD;
passive_DE_TOLE=opts.passive_DE_TOLE;


fit1=[];fit2=[];fit3=[];

Ns=length(s); 

if isempty(poles)
  if isempty(opts.N)  
    disp('===> ERROR in poleresiduefit.m: You did not specify a value for opts.N (fitting order). Must stop.')   
    disp('     Solution: either specify value for opts.N, or provide initial poles in array poles.') 
    return
  end
  N=opts.N;
  oldpoletype=opts.poletype;
  if N<6;
    if strcmp(opts.poletype,'linlogcmplx') 
      opts.poletype='logcmplx';
    end  
  end   
  nu=opts.nu;
  if length(opts.poletype)==8
      
    if opts.poletype=='lincmplx' %Complex, linearly spaced starting poles 
      bet=linspace(s(1)/i,s(Ns)/i,floor(N/2));
      poles=[];
      for n=1:length(bet)
        alf=-nu*bet(n);
        poles=[poles (alf-i*bet(n)) (alf+i*bet(n)) ]; 
      end
    elseif opts.poletype=='logcmplx' %Complex, logaritmically spaced starting poles
      bet=logspace(log10(s(1)/i),log10(s(Ns)/i),floor(N/2));
      poles=[];
      for n=1:length(bet)
        alf=-nu*bet(n);
        poles=[poles (alf-i*bet(n)) (alf+i*bet(n)) ]; 
      end
    else
      disp('-->ERROR in poleresiduefit.m: Illegal value for opts.poletype')
      disp('   Valid input: ''lincmplex'' and ''logcmplx''')  
      disp('   Given input:')
      disp(opts.poletype)
      return    
    end
    
  elseif length(opts.poletype)==11
    if opts.poletype=='linlogcmplx'    
      bet=linspace(s(1)/i,s(Ns)/i,ceil((N-1)/4));
      poles1=[];
      for n=1:length(bet)
        alf=-nu*bet(n);
        poles1=[poles1 (alf-i*bet(n)) (alf+i*bet(n)) ]; 
      end
      bet=logspace(log10(s(1)/i),log10(s(Ns)/i),2+floor(N/4)); 
      bet(1)=[];bet(end)=[];
      poles2=[];
      for n=1:length(bet)
        alf=-nu*bet(n);
        poles2=[poles2 (alf-i*bet(n)) (alf+i*bet(n)) ]; 
      end  
      poles=[poles1 poles2];
   
    else
      disp('-->ERROR in poleresiduefit.m: Illegal value for opts.poletype')
      disp('   Valid input: ''lincmplex'' and ''logcmplx''')  
      disp('   Given input:')
      disp(opts.poletype)
      return          
    end    
  end    

  %if 2*length(bet)<N %An odd number of poles was prescribed
  if length(poles)<N %An odd number of poles was prescribed      
    if strcmp(opts.poletype,'lincmplx')  
      pole_extra=-(s(1)/i+s(end)/i)/2; %Placing surplus pole in midpoint    
    elseif strcmp(opts.poletype,'logcmplx') | strcmp(opts.poletype,'linlogcmplx') 
      pole_extra=-10^((log10(s(1)/i)+log10(s(end)/i))/2); %Placing surplus pole in midpoint          
    end
    poles=[poles pole_extra];
  end
opts.poletype=oldpoletype;
end %if isempty(poles)


 

Nc=length(bigH(:,1,1));
Ns=length(s);

warning('off', 'MATLAB:nearlySingularMatrix')
warning('off', 'MATLAB:rankDeficientMatrix')




if opts.screen==1
tic    
disp('-----------------S T A R T--------------------------')
end

if opts.screen==1
  disp('****Stacking matrix elements (lower triangle) into single column ...')
end  
tell=0;
for col=1:Nc
  for row=col:Nc
    tell=tell+1;
    f(tell,:)=squeeze(bigH(row,col,:)).'; %stacking elements into a single vector
  end
end
nnn=tell;



%Fitting options
VF.spy1=0;
VF.skip_pole=0; 
%VF.skip_pole=opts.skippole, 'VFdriver on line 342'
VF.skip_res=1;
VF.legend=1;

oldspy2=VF.spy2;
VF.spy2=0;

if Nc==1
  f_sum=f;
end  
if Nc>1 %Will do only for multi-terminal case
  %Forming columns sum and associated LS weight:
  f_sum=0;
  tell=0;
  for row=1:Nc
    for col=row:Nc
      tell=tell+1;  
      %if row==col  
        if weightparam==1 | weightparam==4  | weightparam==5 
          f_sum=f_sum+f(tell,:); %unweighted sum 
        elseif weightparam==2  
          f_sum=f_sum+f(tell,:)/norm(f(tell,:));
        elseif weightparam==3            
          f_sum=f_sum+f(tell,:)/sqrt(norm(f(n,:)));
        end  
      %end
    end  
  end
end  

%Creating LS weight
if isempty(opts.weight) %Automatic specification of weight
  if weightparam==1 %weight=1 for all elements in LS problem, at all freq.
    weight=ones(1,Ns); weight_sum=ones(1,Ns); 
  elseif weightparam==2 %Indvidual element weighting   
    weight=1./abs(f);       weight_sum=1./abs(f_sum);            
  elseif weightparam==3 %Indvidual element weighting   
    weight=1./sqrt(abs(f)); weight_sum=1./sqrt(abs(f_sum));
  elseif weightparam==4 %Common weighting for all matrix elements
    for k=1:Ns 
      weight(k)=1./norm(f(:,k));  
    end
    weight_sum=weight;
  elseif weightparam==5 %Common weighting for all matrix elements
    for k=1:Ns 
      weight(k)=1./sqrt(norm(f(:,k)));  
    end
    weight_sum=weight;
  else
    disp('-->ERROR in mtrxVectfit: Illegal value for opts.weight'),return  
  end
else
  weight=zeros(nnn,Ns);
  tell=0;
  for row=1:Nc
    for col=row:Nc  
      tell=tell+1;  
      weight(tell,:)=squeeze(opts.weight(row,col,:));
    end
  end  
  weight_sum=ones(1,Ns);
end  

if Nc>1 %Will do only for multi-terminal case
  if opts.screen==1
    disp('****Calculating improved initial poles by fitting column sum ...')
  end
  for iter=1:Niter1
     if opts.screen==1 
       disp(['   Iter ' num2str(iter)])
     end  
    [SER,poles,rmserr,fit]=vectfit3(f_sum,s,poles,weight_sum,VF);     
  end
end

if opts.screen==1
  disp('****Fitting column ...')
end  
VF.skip_res=1;
for iter=1:Niter2
  if opts.screen==1  
    disp(['   Iter ' num2str(iter)])
  end  
  if iter==Niter2, VF.skip_res=0; end
  %[SER,poles,rmserr,fit1]=vectfit2(f,s,poles,weight,VF);  
  [SER,poles,rmserr,fit1]=vectfit3(f,s,poles,weight,VF);    
end
if Niter2==0
   VF.skip_res=0; VF.skip_pole=1;
   %!![SER,poles,rmserr,fit1]=vectfit2(f,s,poles,weight,VF);     
   [SER,poles,rmserr,fit1]=vectfit3(f,s,poles,weight,VF);       
end

%===========================================
% Throwing out high-frequency poles:
fit2=fit1;
if remove_HFpoles==1 
  if opts.screen==1  
    disp('****Throwing out high-frequency poles: ...')
  end  
  ind=find( abs(poles)>factor_HF*abs(s(end)) ); 
  poles(ind)=[]; %Deleting poles above upper frequency limit
  N=length(poles);
  if opts.screen==1
    disp('****Refitting residues: ...')
  end  
  VF.skip_pole=1;   
  %!![SER,poles,rmserr,fit2]=vectfit2(fit1,s,poles,weight,VF); 
  [SER,poles,rmserr,fit2]=vectfit3(fit1,s,poles,weight,VF);   
end
%===========================================


%===========================================
if passive_DE==1 & VF.asymp>1
  if opts.screen==1  
    if VF.asymp==2  
      disp('****Enforcing positive realness for D...')
    elseif VF.asymp==3  
      disp('****Enforcing positive realness for D, E...')  
    end  
  end  
  tell=0;
  DD=zeros(Nc);EE=zeros(Nc);
  for col=1:Nc
    for row=col:Nc
      tell=tell+1;
      DD(row,col)=SER.D(tell);  EE(row,col)=SER.E(tell);
    end
  end
  DD=DD+(tril(DD,-1)).';
  EE=EE+(tril(EE,-1)).';
%passive_DE_TOLD,passive_DE_TOLE
  %Calculating Dmod, Emod:
  [V,L]=eig(DD); for n=1:Nc, if L(n,n)<0, L(n,n)=passive_DE_TOLD; end, end; DD=V*L*V^(-1);
  [V,L]=eig(EE); for n=1:Nc, if L(n,n)<0, L(n,n)=passive_DE_TOLE; end, end; EE=V*L*V^(-1);
  tell=0;
  %Calculating fmod:
  Emod=zeros(Nc,1);Dmod=zeros(Nc,1);
  for col=1:Nc
    for row=col:Nc
      tell=tell+1;
      Dmod(tell)=DD(row,col); Emod(tell)=EE(row,col);
      fmod(tell,:)=fit2(tell,:) -Dmod(tell) -s.*Emod(tell);
    end
  end
  if opts.screen==1  
    if VF.asymp==2  
      disp('****Refitting C while enforcing D=0...')
    elseif VF.asymp==3  
      disp('****Refitting C while enforcing D=0, E=0...')  
    end 
  end  
  VF.skip_pole=0;
  VF.asymp=1;
  for iter=1:1
    %[SER,poles,rmserr,fit3]=vectfit2(fmod,s,poles,weight,VF);
    [SER,poles,rmserr,fit3]=vectfit3(fmod,s,poles,weight,VF);   
  end  
  SER.D=Dmod;
  SER.E=Emod;
  for tell=1:length(fit3(:,1))
    fit3(tell,:)=fit3(tell,:) +SER.D(tell) +s.*SER.E(tell);
  end  
end %if opts.passive_DE==1  


if Nc>1 
  if opts.screen==1
    toc  
    disp('****Transforming model of lower matrix triangle into state-space model of full matrix ...')
  end
  [SER]=tri2full(SER); 
end

if opts.screen==1
  disp('****Generating pole-residue model ...')
end  
[R,a]=ss2pr(SER.A,SER.B,SER.C);
SER.R=R;
SER.poles=a;


%rmserror of fitting:
if isempty(fit3)==0 
  fit=fit3;
elseif isempty(fit2)==0 
  fit=fit2; 
elseif isempty(fit1)==0 
  fit=fit1; 
end  
diff=fit-f; rmserr=sqrt(sum(sum(abs(diff.^2))))/sqrt(nnn*Ns);


VF.spy2=oldspy2;
if VF.spy2==1
  if opts.screen==1  
    disp('****Plotting of results ...')  
  end  

  freq=s./(2*pi*i);
  if VF.logx==1     
    if VF.logy==1
      figure(1),
      h1=loglog(freq,abs(f),'b'); xlim([freq(1) freq(Ns)]);hold on
      h2=loglog(freq,abs(fit),'r--'); hold off
      if VF.errplot==1
        hold on,h3=loglog(freq,abs(f-fit),'g');hold off; 
      end
    else %logy=0 
      figure(1),
      h1=semilogx(freq,abs(f),'b'); xlim([freq(1) freq(Ns)]);hold on
      h2=semilogx(freq,abs(fit),'r--'); hold off
      if VF.errplot==1
        hold on,h3=semilogx(freq,abs(f-fit),'g');hold off; 
      end
    end
    if VF.phaseplot==1
      figure(2),
      h4=semilogx(freq,180*unwrap(angle(f))/pi,'b'); xlim([freq(1) freq(Ns)]);hold on
      h5=semilogx(freq,180*unwrap(angle(fit))/pi,'r--');hold off
    end  
  else %logx=0
    if VF.logy==1
      figure(1),
      h1=semilogy(freq,abs(f),'b'); xlim([freq(1) freq(Ns)]);hold on
      h2=semilogy(freq,abs(fit),'r--'); hold off
      if VF.errplot== 1
        hold on,h3=semilogy(freq,abs(f-fit),'g');hold off; 
      end
    else %logy=0 
      figure(1), 
      h1=plot(freq,abs(f),'b'); xlim([freq(1) freq(Ns)]);hold on
      h2=plot(freq,abs(fit),'r--'); hold off
      if VF.errplot== 1
        hold on,h3=plot(freq,abs(f-fit),'g');hold off; 
      end
    end
    if VF.phaseplot==1    
      figure(2),
      h4=plot(freq,180*unwrap(angle(f))/pi,'b'); xlim([freq(1) freq(Ns)]);hold on
      h5=plot(freq,180*unwrap(angle(fit))/pi,'r--');hold off
    end  
  end %logy=0
  figure(1),
  xlabel('Frequency [Hz]'); ylabel('Magnitude [p.u.]');
  title('Approximation of f');
  if VF.legend==1
    if VF.errplot==1
      legend([h1(1) h2(1) h3(1)],'Original','FRVF','Deviation');
    else
      legend([h1(1) h2(1)],'Original','FRVF');       
    end  
  end
  if VF.phaseplot==1
    figure(2),
    xlabel('Frequency [Hz]'); ylabel('Phase angle [deg]');
    title('Approximation of f');
    if VF.legend==1
      legend([h4(1) h5(1)],'Original','FRVF'); 
    end
  end
  drawnow;
end

if opts.screen==1
disp('-------------------E N D----------------------------')
end

bigHfit=zeros(Nc,Nc,Ns);
tell=0;
for row=1:Nc
  for col=row:Nc
    tell=tell+1; 
    bigHfit(row,col,:)=fit(tell,:);  
    if row~=col
      bigHfit(col,row,:)=fit(tell,:);  
    end  
  end
end  

warning('on', 'MATLAB:nearlySingularMatrix')
warning('on', 'MATLAB:rankDeficientMatrix')
end



function [SER2]=tri2full(SER);
%
% Convert rational model of lower matrix triangle into state-space model of
% full matrix.
%
% Input:
% SER: Output structure from vectfit2 when fitting lower triangle of a square matrix
% using a common pole set.
% Both formats determined by parameter VF.cmplx_ss are valid
% 
% Output:
% SER: State-space model of full matrix (with common pole set)
%
% This routine is part of the vector fitting package (v2.1) 
% Last revised: 27.10.2005. 
% Created by:   Bjorn Gustavsen.
%

A=SER.A; B=SER.B; C=SER.C; D=SER.D; E=SER.E;

tell=0;
for k=1:1e4
  tell=tell+k;
  if tell==length(D), Nc=k; break, end
end

N=length(A);
tell=0;
CC=zeros(Nc,Nc*N);
AA=[]; BB=[];
for col=1:Nc
  AA=blkdiag(AA,A);
  BB=blkdiag(BB,B);  
  for row=col:Nc
    tell=tell+1;
    DD(row,col)=D(tell);
    EE(row,col)=E(tell);
    CC(row,(col-1)*N+1:col*N)=C(tell,:); 
    CC(col,(row-1)*N+1:row*N)=C(tell,:);     
  end
end
DD=DD+(DD-diag(diag(DD))).';
EE=EE+(EE-diag(diag(EE))).';

SER2.A=AA; SER2.B=BB; SER2.C=CC; SER2.D=DD; SER2.E=EE;

end


function [R,a]=ss2pr(A,B,C)
%
% Convert state-space model having COMMON POLE SET into pole-residue model.
%
% Input:
% A,B,C: must have the format produced by vectfit2.m. Both
% formats determined by parameter VF.cmplx_ss are valid
% 
% Output:
% R(Nc,Nc,N) %Residues
% a(N)       %poles 
%
% This routine is part of the vector fitting package (v2.1) 
% Last revised: 27.10.2005. 
% Created by:   Bjorn Gustavsen.
%

%Converting real-only state-space model into complex model, if necessary:  
if max(max(abs(A-diag(diag(A)))))~=0
  errflag=0;
  for m=1:length(A)-1
    if A(m,m+1)~=0
      A(m,m)    =A(m,m)+i*A(m,m+1); 
      A(m+1,m+1)=A(m+1,m+1)-i*A(m,m+1);  
      
      B(m,:)  =(B(m,:)+B(m+1,:))/2;
      B(m+1,:)=B(m,:); 
      
      C(:,m)  =C(:,m)+i*C(:,m+1);
      C(:,m+1)=conj(C(:,m));  
    end
  end
end  

%Converting complex state-space model into pole-residue model
Nc=length(C(:,1));
N=length(A)/Nc;
R=zeros(Nc,Nc,N);
for m=1:N
  Rdum=zeros(Nc);
  for n=1:Nc
    ind=(n-1)*N+m;
    Rdum=Rdum +C(:,ind)*B(ind,:);
  end  
  R(:,:,m)=Rdum;
end
a=full(diag(A(1:N,1:N)));

end
