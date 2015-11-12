%% Why sssMOR? - The curse of dimensionality
%
% So why is there a need for sssMOR?
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
% R^p$ ($p,m \ll  N$) represent the state, input and output of the system
% respectively.
%	
% Note that in a large-scale setting, $N$ can grow quite large, easily
% above an order of magnitude of $\mathcal{O}(10^4)$. This poses a big
% challenge on the numerical treatment of such models, in first place due
% to storage limitations. In fact, following code would cause MATLAB to
% quit with an error
N = 10^5; A = eye(N);
%%
% The probles is that MATLAB is trying to store all $N^2 = 10^10$ entries
% of the matrix. It is easy to estimate that in IEEE Standard 754
% double-precision floating point computation, this requires roughly 8
% Byte$\cdot10^{10}$ = 80 GB of memory.
% 
% Fortunately in this case, the information of the identity matrix eye(N)
% is only on the diagonal, i.e. there are only $N=10^5$ nonzero entries. By
% recognizing and exploting the sparsity of the problem, for example
% calling,
N = 10^5; A = speye(N); whos A
%%
% MATLAB stores only a triple of indices and values for each nonzero entry
% of A, reducing the storage requirements from 80GB to 2.4 MB! Matrices
% that have large dimensions but only a significantly smaller amount of
% nonzero entries are called sparse and can be defined in MATLAB by using
% appropriate function like speye(N) in the example before.
%
% In fact, in this fashion we are able to double the order of magnitude of
% storable identity matrices to $\mathcal{O}(10^8)$
N = 10^8; A = speye(N); whos A
%%	
% A similar problem arises when considering *large-scale dynamical
% systems*. In fact, the number $N$ of state variables in $x$ can easily
% grow larger than $\mathcal{O}(10^4)$, making it impossible to define and
% analyze state-space objects (\mcode{ss} or \mcode{dss}) with
% MATLAB's Control System Toolbox. Luckily also for large-scale systems,
% high dimensionality goes along with a high degree of sparsity, as it is
% shown in figures below for a selection of benchmark problems.
%
% <<img/iss.jpg>>
%
% <<img/iss_sparsity.png>>
%
% <<img/gyro.jpg>>
% 
% <<img/gyro_sparsity.png>>   
%
% <<img/rail.jpg>>
% 
% <<img/rail_5177_sparsity.png>>  
%
%%	
% Unfortunaltely, the built-in functions _ss_ and _dss_ convert all
% matrices to full, so that sparsity of the system matrices cannot be
% exploited and the size limitation stated above hold, as demonstrated by
% this short code
%
N = 10^8; A = speye(N); b = rand(N,1);
sys = ss(A,b,b',0);	
%%
% *Is there a way* to deal with such *large-scale models* in MATLAB
% nonetheless?
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