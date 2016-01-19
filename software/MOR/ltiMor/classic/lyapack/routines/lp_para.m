function [opts,err_code,rw,Hp,Hm] = lp_para(usfs,Bf,Kf,opts,b0)
%
%  Estimation of suboptimal ADI shift parameters for the matrix 
%
%    F = A-Bf*Kf'.  
%
%  Calling sequence:
%
%    [opts,err_code,rw,Hp,Hm] = lp_para(usfs,Bf,Kf,opts,b0)
%
%  Input:
%
%    usfs      structure of user supplied functions for matrix operations
%    Bf        matrix Bf;
%              Set Bf = [] if not existing or zero!
%    Kf        matrix Kf;
%              Set Kf = [] if not existing or zero!
%    b0         Arnoldi start vector (optional; chosen at random if not
%              provided).
%    opts      structure of input parameters with members:
%        l0        desired number of shift parameters (kp+km > 2*l0)
%                  (The algorithm delivers either l0 or l0+1 parameters!);
%        kp, km    numbers of Arnoldi steps w.r.t. F and inv(F), 
%                  respectively (kp, km << order of A);
%        meth      'heur'  heuristic choice of parameters see [1]
%                  'h_ro'  as above but takes only the real parameters
%                          computed to avoid complex arithmetics and
%                          storage requirement 
%                  'h_rp'  as 'heur' but uses the real parts of the
%                          possibly complex parameters to avoid complex
%                          arithmetics and storage requirement
%                  'wach'  (sub)optimal parameters following the theory
%                          of Wachspress et.al. see [2], [3] for
%                          details using Ritz values for spectral approximation
%                  'wach2' (sub)optimal parameters following the theory
%                          of Wachspress et.al. see [2], [3] for
%                          details using eigs for the spectral approximation
%
%  Output:
%
%    opts      is extended by
%              p    an l0- or l0+1-vector of suboptimal ADI parameters;
%    err_code  Error code; = 1, if Ritz values with positive real parts
%              have been encountered; otherwise, err_code = 0;
%    rw        vector containing the Ritz values;
%    Hp        Hessenberg matrix in Arnoldi process w.r.t. F;
%    Hm        Hessenberg matrix in Arnoldi process w.r.t. inv(F);
%
%  User-supplied functions called by this function:
%
%    'usfs.m', 'usfs.l' or 'usfs.e' depending on the method used.
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
%  [2] P. Benner, H. Mena, J. Saak
%      On the Parameter Selection Problem in the Newton-ADI
%      Iteration for Large-Scale Riccati Equations; 
%      Chemnitz Scientific Computing Preprints 06-03, TU Chemnitz;
%      2006. ISBN/ISSN: 1864-0087 
%  
%  [3] E. Wachspress
%      The ADI Model Problem
%      1995 self published
% 
%  LYAPACK 1.8 (Jens Saak, February 2008)

% Input data not completely checked!

err_code = 0;
rw=[];
Hp=[];
Hm=[];
p=[];

switch opts.method
 case 'heur'
  [p,err_code,rw,Hp,Hm] = lp_para_heur(usfs,Bf,Kf,opts,b0);
 case 'h_ro'
  [p,err_code,rw,Hp,Hm] = lp_para_heur(usfs,Bf,Kf,opts,b0);
  for cnt = 1:length(p0)
    if isreal(p0(cnt))
      p=[p p0(cnt)];
    end
  end
 case 'h_rp'
  [p,err_code,rw,Hp,Hm] = lp_para_heur(usfs,Bf,Kf,opts,b0);
  p=real(p);
 case 'wach' 
  p = lp_para_wach(usfs,Bf,Kf,opts,b0);
 case 'wach2' 
  p = lp_para_wach2(usfs,Bf,Kf);
 otherwise
  error('Unknown ADI shift parameter computation method requested!');
end
opts.p = p;




