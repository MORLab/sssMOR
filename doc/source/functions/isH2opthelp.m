%% isH2opt
% evaluate Maier-Luenberger conditions for H2-optimality
%
%% Syntax
%
% <html>
%        <div class="syntax">
% isH2opt(sys,sysr,s0) <br>
% isH2opt(sys,sysr,s0,Opts) <br>
% isH2opt = isH2opt(sys,sysr,s0) <br>
%        </div>
% </html>
%
%% Examples
% No examples
%
%
%% Description
% This function evaluates the Maier-Luenberger conditions for
% H2-optimality according to the input data, i.e.:
%
% * sys, the high fidelity model
% * sysr, the reduced order model
% * s0,  the shifts at which sysr interpolates sys
%
% Following conditions must be met:
%
% # the reduced eigenvalues are the mirror images of the shifts
% # two moments are matched at each shift
%
% Then sysr is said to be a locally H2-optimal approximation of sys
%
%
%% See also
% <momentshelp.html |moments|>, <irkahelp.html |irka|>, <matlab:doc('spark') spark>, <matlab:doc('eig') eig>
%
%% References
% * *[1] Gugercin et al. (2008)*, H2 model reduction for large-scale linear dynamical systems
%
%

%%
% <html>
%   <hr>
%   <p class="copy">&copy; 2015 RT Universit&auml;t M&uuml;nchen
%        <tt class="minicdot">&#149;</tt>
%        <a href="http://www.rt.mw.tum.de">Website</a>
%        <tt class="minicdot">&#149;</tt>
%        <a href="file:./LICENSE.txt">Terms of Use</a>
%        <tt class="minicdot">&#149;</tt>
%        <a href="file:./README.txt">Read me</a>
%   </p>
% </html>
