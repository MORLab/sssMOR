%% sssMOR - model order reduction
%
% Even when using sss functions to exploit the sparsity, 
% computations such as simulations, optimization and control design algorithms based on the full order models (FOM) 
% will require a substantial amount of time, provided they can be carried through.
%
% For this reason, in the large-scale setting we often seek reduced 
% order models (ROM) of much smaller size that capture the relevant dynamics. 
%
% For the benchmark examples above, the figures below illustrate how 
% appropriate model order reduction (MOR) techniques can yield good approximations of the dynamics and drastically reduced the model size.
%	
%   << iss_MOR.png >>
%
%   <<gyro_MOR.png >>
%	
%	For linear systems as in \eqref{eq:FOM}, this is generally done by applying Petrov-Galerkin projections of the form
% 
% $W^T E V \dot{x}_r = W^T A V x_r +  W^T B u$
% 
% $ y_r = C V x_r + D u $
%
% There are several ways to compute the projection matrices $V, W$ depending on what properties of the FOM should be preserved.
%
% *sssMOR* contains some classic reduction methods such as _modal truncation_, _truncated balanced realizations_ and _rational Krylov subspace methods_. 
% At the same time, it implements so well known state-of-the-art algorithms 
% like the _iterative rational Krylov algorithm_ (IRKA) as well as some more 
% recent algorithms such as _cumulative redction framework_ (CURE) and the _stability-preserving, rational Krylov algorithm_ (SPARK)
% 

%%
% <html>
%   <hr>
%   <p class="copy">&copy; 2015 RT Technische Universit&auml;t M&uuml;nchen
%        <tt class="minicdot">&#149;</tt>
%        <a href="http://www.rt.mw.tum.de">Website</a>
%        <tt class="minicdot">&#149;</tt>
%        <a href="file:./LICENSE.txt">Terms of Use</a>
%        <tt class="minicdot">&#149;</tt>
%        <a href="file:./README.txt">Read me</a>
%   </p>
%   <div style="position: relative; width:100%; height: 100px;">
%        <div style="position: absolute; left:0%;">
%             <img src="img/logo_sssMOR_long.png" alt="sssMOR_Logo" height="40px">
%        </div>
%        <div style="position: absolute; left:20%;">
%             <img src="img/MORLAB_Logo.jpg" alt="MORLAB_Logo" height="40px">
%        </div>
%        <div style="position: absolute; right:25%;">
%             <img src="img/Logo_Textzusatz_rechts_engl_Chair.png" alt="RT_Logo" height="40px">
%        </div>
%        <div style="position: absolute; right:0%;">
%             <img src="img/TUM-logo.png" alt="TUM_Logo" height="40px">
%        </div>
%   </div>
% </html>