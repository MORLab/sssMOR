%% sssMOR - Model Order Reduction
%
% Even when using |sss| functions to exploit the sparsity, computations such
% as simulations, optimization and control design algorithms based on the
% full order models (FOM) will require a substantial amount of time,
% provided they can be carried through.
%
% For this reason, in the large-scale setting we often seek reduced order
% models (ROM) of much smaller dimension that capture the relevant dynamics.
%
% For the ISS benchmark example above, the figure below illustrates how
% appropriate model order reduction (MOR) techniques can yield good
% approximations of the dynamics and drastically reduced the model size.
%

%%
%
% <html>
% <div style="text-align: center; position: relative; height: 400px; width:100%;">
%     <img src="img/iss_MOR.png" alt="iss" width="500" style="vertical-align: bottom;">
% </div>
% </html>
%

%%
% For linear systems, this is generally done by applying Petrov-Galerkin 
% projections of the form:
% 
% $W^T E V \dot{x}_r = W^T A V x_r +  W^T B u$
% 
% $y_r = C V x_r + D u$
%
% There are several ways to compute the projection matrices $V, W$
% depending on what properties of the FOM should be preserved.
%
% *sssMOR* contains some classic reduction methods such as _modal
% truncation_, _truncated balanced realizations_ and _rational Krylov
% subspace methods_. At the same time, it implements well-known
% state-of-the-art algorithms like the _iterative rational Krylov
% algorithm_ (IRKA) as well as some more recent algorithms such as
% _CUmulative REduction framework_ (CURE) and the _Stability-Preserving,
% Adaptive Rational Krylov algorithm_ (SPARK).
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
