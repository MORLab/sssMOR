%% What is sssMOR?
%
% *sssMOR* is a sparse state-space and model order reduction toolbox for
% MATLAB designed to work with large-scale dynamical systems in
% state-space.
%
% The *sssMOR* toolbox is composed of two parts: *sss* and *sssMOR*.
%
% *sss* extends the capabilities of the Control System Toolbox by defining
% sparse state-space (sss) objects and implementing large-scale and
% sparsity exploiting versions of common system and control functions (such
% as bode, step, pzmap, eig,...).
%
% *sssMOR* contains classic and state-of-the-art model order reduction
% (MOR) techniques to capture the dynamics of large-scale systems in
% reduced order models of significantly smaller size.
%
% By using *sss* and *sssMOR*, it is possible to define and analyze dynamic
% system objects with state-space dimensions higher than
% $\mathcal{O}(10^4)$, which is generally the limit for standard built-in
% ss and dss objects.
%
% *sssMOR* is a MATLAB toolbox developed at the Model Order Reduction Lab
% (MORLAB) of the Chair of Automatic Control of TU M&uuml;nchen.
%
% *sss* is a MATLAB toolbox developed at MORLAB in collaboration with the
% Chair of Thermofluid Dynamics, TU M&uuml;nchen.
%
% 

%%
% <html>
%   <hr>
%   <p class="copy">&copy; 2016 RT Technische Universit&auml;t M&uuml;nchen
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
