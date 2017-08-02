%% Functions Reference
% UI-functions of the sssMOR toolbox sorted by categories

%% Reduced Order Models (ssRed)
%
% * <.\ssRedhelp.html |ssRed|>  --  Reduced state-space LTI system (ssRed) class
% * <.\l2normhelp.html |l2norm|>  --  L2-norm of a dynamical system (ssRed)
% * <.\stabsephelp.html |stabsep|>  --  Stable-unstable decomposition

%% LTI MOR Classic
%
% * <.\arnoldihelp.html |arnoldi|>  --  Arnoldi algorithm for Krylov subspaces with multiple shifts
% * <.\modalMorhelp.html |modalMor|>  --  Modal truncation order reduction of LTI systems
% * <.\projectiveMorhelp.html |projectiveMor|>  --  Reduces dynamic model by projection
% * <.\rkhelp.html |rk|>  --  Model Order Reduction by Rational Krylov
% * <.\tbrhelp.html |tbr|>  --  Performs model order reduction by Truncated Balanced Realization

%% LTI MOR State-of-the-Art
%
% * <.\cirkahelp.html |cirka|>  --  Confined Iterative Rational Krylov Algorithm
% * <.\curehelp.html |cure|>  --  CUmulative REduction framework
% * <.\irkahelp.html |irka|>  --  Iterative Rational Krylov Algorithm
% * <.\isrkhelp.html |isrk|>  --  Iterative SVD-Rational Krylov Algorithm
% * <.\modelFcthelp.html |modelFct|>  --  computes or updates the model function of an sss object
% * <.\modelFctMorhelp.html |modelFctMor|>  --  model function-based model order reduction
% * <.\porkVhelp.html |porkV|>  --  Pseudo-Optimal Rational Krylov (Input)
% * <.\porkWhelp.html |porkW|>  --  Pseudo-Optimal Rational Krylov (Output)
% * <.\rkIcophelp.html |rkIcop|>  --  Rational Krylov with an iteratively calculated optimal point
% * <.\rkOphelp.html |rkOp|>  --  Determination of optimal expansion point for Laguerre series
% * <.\sparkhelp.html |spark|>  --  Stability Preserving Adaptive Rational Krylov

%% Extras
%
% * <.\getDesiredOutputhelp.html |getDesiredOutput|>  --  get only the desired output from a function
% * <.\getSylvesterhelp.html |getSylvester|>  --  Get matrices of Sylvester's equation for Krylov subspaces
% * <.\isH2opthelp.html |isH2opt|>  --  Evaluate Meier-Luenberger necessary conditions for H2-optimality
% * <.\ismemberf2help.html |ismemberf2|>  --  determine which elements are within a radius
% * <.\momentshelp.html |moments|>  --  Returns the moments or Markov parameters of an LTI system
% * <.\setdiffVechelp.html |setdiffVec|>  --  Computes the difference between two unsorted vectors
% * <.\shiftVechelp.html |shiftVec|>  --  convert shift definition to single row notation
% * <.\cplxpairAllhelp.html |cplxpairAll|>  --  Sort arrays of complex numbers into complex conjugate pairs

%% Demos
%
% * <.\sssMOR_gettingStartedhelp.html |sssMOR_gettingStarted|>  --  Introductory demo to sssMOR toolbox
%
%%
% <html>
%   <hr>
%   <p class="copy">&copy; 2015-2017 RT Technische Universit&auml;t M&uuml;nchen
%        <tt class="minicdot">&#149;</tt>
%        <a href="https://www.rt.mw.tum.de/?sssMOR">Website</a>
%        <tt class="minicdot">&#149;</tt>
%        <a href="file:txts/LICENSE.txt">License</a>
%        <tt class="minicdot">&#149;</tt>
%        <a href="file:txts/RELEASE.txt">Release Notes</a>
%   </p>
% <div>
% <table>
%  <tr>
%   <td style="background-color:#ffffff; border:0; width:25%; vertical-align:middle; text-align:center">
%             <img src="img/logo_sssMOR_long.png" alt="sssMOR_Logo" height="40px">
%      </td>
%   <td style="background-color:#ffffff; border:0; width:25%; vertical-align:middle; text-align:center">
%      <img src="img/MORLAB_Logo.jpg" alt="MORLAB_Logo" height="40px"></td>
%   <td style="background-color:#ffffff; border:0; width:25%; vertical-align:middle; text-align:center">
%      <img src="img/Logo_Textzusatz_rechts_engl_Chair.png" alt="RT_Logo" height="40px"></td>
%   <td style="background-color:#ffffff; border:0; width:25%; vertical-align:middle; text-align:center">
%      <img src="img/TUM-logo.png" alt="TUM_Logo" height="40px"></td>
%  </tr>
% </table>
% </div>
% </html>
