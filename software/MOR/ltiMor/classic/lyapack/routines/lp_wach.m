function p=lp_wach(a,b,alpha,TOL)
%
% function p=myadipars(a,b,alpha,TOL) 
%  
% calculates the optimal ADI shiftparameters (for equations where
% Matrix A is stable and symmetric) according to Jing-Rebecca Li
% and Jakob Whites "Low Rank Solution of Lyapunov equation" (which
% gives an overview of Wachspress`s method form e.g. "The ADI model Problem" 
%
% p   is the array of shift parameters
%  
% a   is assumed to be the absolute value of the smallest magnitude
%     eigenvalue of the Systemmatrix A 
%
% b   is assumed to be the absolute value of the largest magnitude eigenvalue
%
% alpha is the arctan of the maximum of abs(imag(lamba))/abs(real(lambda))
%       for all lambda in the spectrum of A
%
% TOL is the epsilon1 of the above paper. The smaller the better
%     the approximation but the more parameters are calculated 
%
  if (alpha==0)
    kprime=a/b;
  else
    c2 = 2/(1+(a/b+b/a)/2);
    m = 2*cos(alpha)*cos(alpha)/c2 -1;
    if (m<1)
      error(['Shift parameters would be complex, function not aplicable, ' ...
             'aborting!']); 
      
      %
      % FIX ME: if m<1 parameter become complex! switch back to the
      % heuristics by Thilo or complex parameters suggested by
      % Wachspress.
      % 
      % ALA2006 -> also test switching to Krylov projection based
      % method KPIK by V.Simoncini
      %
    end
    kprime = 1/(m+sqrt(m^2-1));
  end
  k=min(1-eps,sqrt(1-kprime^2));%this is a workaround for the case
				%k=1 that works for Peters Model
				%reduction problems not really
				%nice but it works great.
  
  %TODO: check the computation of k, kprime to avoid roundoff errors
  %and probably replace the hack above. 
  
  [K,E]=ellip(k,pi/2);
  if (alpha==0)
    [v,E]=ellip(kprime,pi/2);
  else
    [v,E]=ellip(kprime,asin(sqrt(a/(b*kprime))));
  end
  J=ceil(K/(2*v*pi)*log(4/TOL));
  
  p=ones(J,1);
  for i=1:J
    p(i)=-sqrt(a*b/kprime)*dn((i-0.5)*K/J,k); 
    %here we have the choice to take the
    %matlab function ellipj or the local
    %one dn. the later can be ported to
    %FORTRAN or C Code very easily
  end
  
function [F,E]=ellip(hk,phi)
%  function [F,E]=ellip(hk,phi);
%  Computes complete and incomplete elliptic integrals F(k,phi) and E(k,phi)
%       Input  : hk  --- Modulus k ( 0 < k < 1 )
%                phi --- Argument 
%       Output : F   --- F(k,phi)
%                E   --- E(k,phi)
%       ==================================================
  
  g=0.0;
  a0=1.0;
  b0=min(1-eps,sqrt(1.0-hk.*hk));
  d0=phi;
  r=hk.*hk;
  if (hk == 1.0 && phi == pi/2) ;
    F=1.0e+300;
    E=1.0;
  elseif (hk == 1.0);
    F=log((1.0+sin(d0))./cos(d0));
    E=sin(d0);
  else
    fac=1.0;
    for  n=1:40;
      a=(a0+b0)./2.0;
      b=sqrt(a0.*b0);
      c=(a0-b0)./2.0;
      fac=2.0.*fac;
      r=r+fac.*c.*c;
      if (phi ~= pi/2) ;
        d=d0+atan((b0./a0).*tan(d0));
        g=g+c.*sin(d);
        d0=d+pi.*fix(d./pi+.5);
      end;
      a0=a;
      b0=b;
      if (c < 1.0e-15) break; end;
    end;
    ck=pi./(2.0.*a);
    ce=pi.*(2.0-r)./(4.0.*a);
    if (phi == pi/2) ;
      F=ck;
      E=ce;
    else
      F=d0./(fac.*a);
      E=F.*ce./ck+g;
    end;
  end;
  return;
  
function out=dn(u,k)
%function out=dn(u,k) calculates the value of the elliptic
%function dn (see Abramowitz/Stegun Handbook of mathematical
%functions '65)
%
  
  a(1)=1;
  c(1)=k;
  b(1)=min(1-eps,sqrt(1-k*k));
  i=1;
  while abs(c(i))>eps
    i=i+1;
    a(i)=(a(i-1)+b(i-1))/2;
    b(i)=sqrt(a(i-1)*b(i-1));
    c(i)=(a(i-1)-b(i-1))/2;
  end
  phi1=(2.^(i-1)).*a(i).*u;%here 2^(i-1) and not 2^i as in
                           %Abramowitz/Stegun because counting
                           %starts at 1 not at 0 like in the book
  phi0=0;
  for j=i:-1:2 
    if (j<i) 
      phi1=phi0;
    end;
    phi0=(phi1+asin(c(j)*sin(rem(phi1,2*pi))/a(j)))/2;
  end
  arg=1-k*k*sin(rem(phi0,2*pi))^2;
  if (arg<.1)
    out=sqrt(arg);
  else
    out=cos(rem(phi0,2*pi))/cos(phi1-phi0);
  end
  
  %the last two are both representations found in the
  %Abramowitz/Stegun book. if arg is close to zero the cosine
  %version should be better to avoid numerical inexactness resulting
  %from the substraction.  