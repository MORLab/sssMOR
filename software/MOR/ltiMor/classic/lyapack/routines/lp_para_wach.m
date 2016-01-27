function [p,err_code,rw,Hp,Hm] = lp_para_wach(usfs,Bf,Kf,opts,b0)
%
%  Estimation of suboptimal ADI shift parameters for the matrix 
%
%    F = A-Bf*Kf'.  
%
%  Calling sequence:
%
%    [p,err_code,rw,Hp,Hm] = lp_para(usfs,Bf,Kf,opts,b0)
%
%  Input:
%
%    usfs      structure of function_handles for the user supplied functions;
%    Bf        matrix Bf;
%              Set Bf = [] if not existing or zero!
%    Kf        matrix Kf;
%              Set Kf = [] if not existing or zero!
%    opts      parameter structure with members:
%        l0        desired number of shift parameters (kp+km > 2*l0)
%                  (The algorithm delivers either l0 or l0+1 parameters!);
%        kp, km    numbers of Arnoldi steps w.r.t. F and inv(F), 
%                  respectively (kp, km << order of A);
%
%    b0         Arnoldi start vector (optional; chosen at random if not
%              provided).
%
%  Output:
%
%    p         an l0- or l0+1-vector of suboptimal ADI parameters;
%    err_code  Error code; = 1, if Ritz values with positive real parts
%              have been encountered; otherwise, err_code = 0;
%    rw        vector containing the Ritz values;
%    Hp        Hessenberg matrix in Arnoldi process w.r.t. F;
%    Hm        Hessenberg matrix in Arnoldi process w.r.t. inv(F);
%
%  User-supplied functions called by this function:
%
%    'usfs.m', 'usfs.l'    
%
%  Remarks:
%
%    Typical values are l0 = 10..40, kp = 20..80, km = 10..40.
%    The harder the problem is the large values are necessary.
%    Larger values mostly result in a faster convergence, but also in a
%    larger memory requirement.
%    However, for "well-conditioned" problems small values of l0 can 
%    lead to the optimal performance.
%
%  References:
%
%  [1] T. Penzl.
%      LYAPACK (Users' Guide - Version 1.0).
%      1999.
%
%   
%  LYAPACK 1.8 (Jens Saak, February 2008)

% Input data not completely checked!

ni = nargin;

err_code = 0;

n = feval(usfs.m);                    % Get system order.
if (ni < 4)||(isempty(opts))||(isempty(opts.l0))||(isempty(opts.kp)) ...
      ||(isempty(opts.km)) 
  error('missing process parameters kp, km or l0'); 
end
if opts.kp >= n, error('kp must be smaller than n!'); end
if opts.km >= n, error('km must be smaller than n!'); end
if 2*opts.l0 >= opts.kp+opts.km, error('2*l0 must be smaller than kp+km!'); end

if ni < 5
  b0 = randn(n,1);
end
b0 = (1/norm(b0))*b0;

rwp = [];
rwm = [];
rw = [];
Hp = [];
Hm = [];

if opts.kp > 0
  [Hp,V] = lp_arn_p(usfs,Bf,Kf,opts.kp,b0);
  rwp = eig(Hp(1:opts.kp,1:opts.kp));                  % =: R_+
  rw = [rw; rwp];
end

if opts.km > 0
  [Hm,V] = lp_arn_m(usfs,Bf,Kf,opts.km,b0);
  rwm = ones(opts.km,1)./eig(Hm(1:opts.km,1:opts.km));      % =: 1 / R_- 
  rw = [rw; rwm];                           % =: R 
end

if any(real(rw) >= zeros(size(rw)))
  err_code = 1;
  disp('These are the Ritz values computed by the Arnoldi process w.r.t. F:')
  disp(rwp)
  disp('These are the Ritz values computed by the Arnoldi process w.r.t. inv(F):')
  disp(rwm)
  disp(' ');
  disp('####################################################################');
  disp('WARNING in ''lp_para'': NON-STABLE RITZ VALUES DETECTED!!!')
  disp(' ');
  disp('This is quite a serious problem, that can be caused by  ');
  disp('(i)   non-stable matrices F (Be sure that F is stable. ADI like');
  disp('      methods only work for stable or antistable problems. If your');
  disp('      Lyapunov equation is antistable, multiply it by -1.)');
  disp('(ii)  matrices F that are stable but have an indefinite symmetric')
  disp('      part (This is THE weak point of this algorithm. Try to work')
  disp('      with the "reduced" Ritz values, i.e., the unstable values are')
  disp('      simply removed. This is not an elegant measure but it may work.')
  disp('      However, the convergence of ADI can be poor. This measure is')
  disp('      taken automatically. Another measure might be to enlarge the')
  disp('      values of kp or km, and run the program again.')
  disp('(iii) matrices F with a negative definite, but ill-conditioned')
  disp('      symmetric part (This is quite unlikely. The problem is')
  disp('      caused by round-off errors).')
  disp(' ')
  disp('#####################################################################')
  disp(' ');
  disp(' ');
  disp('NOTE: The unstable Ritz values will be ignored in the further computation!!! ');
  disp(' ')
  pause(3);
  rw0 = rw; rw = [];
  for j = 1:length(rw0)
    if real(rw0(j))<0
      rw = [rw; rw0(j)];
    end

  end
end  
a=min(-real(rw));
b=max(-real(rw));
alpha=atan(max(imag(rw)./real(rw)));
p = lp_wach(a,b,alpha,eps);






