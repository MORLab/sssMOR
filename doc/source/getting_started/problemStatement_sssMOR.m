%% Why sssMOR? - The curse of dimensionality
%
% *So why is there a need for sssMOR?*
%
% The accurate modeling of dynamical systems often results in a large
% number of state variables and differential equations describing the
% system behavior in time. This is often the case, for instance, when
% discretizing partial differential equations in space over a fine grid or
% in the large-scale integration of electrical circuits.
%
% If the dynamics are linear, which is often the case at least in a region
% around the operating point, then the system can by modeled using
% state-space equations of the form
%
% $E \dot{x} = A x + B u$
%
% $y = C x + D u$
%
% where $E \in R^{N  \times  N}$ is the _descriptor matrix_, $A \in R^{N
% \times  N}$ is the system matrix and $x \in R^N$, $u \in R^m$, $y \in
% R^p$ ($p,m \ll  N$) represent the state, input and output variables of 
% the system, respectively.
%	
% Note that in a large-scale setting, $N$ can grow quite large, easily
% above an order of magnitude of $\mathcal{O}(10^4)$. This poses a big
% challenge on the numerical treatment of such models, in first place due
% to storage limitations. In fact, following code would cause MATLAB to
% quit with an error
try
    A = eye(1e5);
catch err
    fprintf('%s\n',err.message);
end
%%
% The probles is that MATLAB is trying to store all $N^2 = 10^{10}$ entries
% of the matrix. It is easy to estimate that in IEEE Standard 754
% double-precision floating point computation, this requires roughly 8
% Byte $\cdot10^{10}$ = 80 GB of memory.
% 
% Fortunately in this case, the information of the identity matrix eye(N)
% is only on the diagonal, i.e. there are only $N=10^5$ nonzero entries. By
% recognizing and exploting the sparsity of the problem, for example
% calling,
A = speye(1e5); whos A
%%
% MATLAB stores only a triple of indices and values for each nonzero entry
% of A, reducing the storage requirements from 80GB to 2.4 MB! Matrices
% that have large dimensions but only a significantly smaller amount of
% nonzero entries are called _sparse_ and can be defined in MATLAB by using
% appropriate functions like speye(N) in the example before.
%
% In fact, in this fashion we are able to double the order of magnitude of
% storable identity matrices to $\mathcal{O}(10^8)$:
N = 10^8; A = speye(N); whos A
%%	
% A similar problem arises when considering *large-scale dynamical
% systems*. In fact, the number $N$ of state variables in $x$ can easily
% grow larger than $\mathcal{O}(10^4)$, making it impossible to define and
% analyze state-space objects ( |ss| or |dss| ) with MATLAB's Control 
% System Toolbox. Luckily also for large-scale systems, high dimensionality
% goes along with a high degree of sparsity, as it is
% shown in figures below for a selection of benchmark problems.
%

%%
% <html>
% <div style="text-align: center; position: relative; height: 400px; width:100%;">
%     <div style="text-align: left; height: 150px; width:100%; display: inline-block">
%         <div style="text-align: center; width:33%; display: inline-block">
%             <img src="img/iss.jpg" alt="iss" width="200" style="vertical-align: bottom;">
%         </div>
%         <div style="text-align: center; width:32%; display: inline-block">
%             <img src="img/gyro.jpg" alt="gyro" width="200" style="vertical-align: bottom;">
%         </div>
%         <div style="text-align: center; width:33%; display: inline-block;">
%             <img src="img/rail.jpg" alt="rail" width="200" style="vertical-align: bottom;">
%         </div>
%     </div>
%     <div style="text-align: left; height: 150px; width:100%; display: inline-block">
%         <div style="text-align: center; width:33%; display: inline-block">
%             <img src="img/iss_sparsity.png" alt="iss" width="250" style="vertical-align: top;">
%         </div>
%         <div style="text-align: center; width:32%; display: inline-block">
%             <img src="img/gyro_sparsity.png" alt="gyro" width="250" style="vertical-align: top;">
%         </div>
%         <div style="text-align: center; width:33%; display: inline-block;">
%             <img src="img/rail_5177_sparsity.png" alt="rail" width="250" style="vertical-align: top;">
%         </div>
%     </div>
% </div>
% </html>
% 

%%	
% Unfortunaltely, the built-in functions |ss| and |dss| convert all
% matrices to full, so that sparsity of the system matrices cannot be
% exploited and the size limitation stated above hold, as demonstrated by
% this short code:
%
A = speye(1e5); b = rand(1e5,1); c = b';
try
    sys = ss(A,b,c,[]);
catch err
    fprintf('%s\n',err.message);
end
%%
% *Is there a way to deal with such large-scale models in MATLAB
% nonetheless?*
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
