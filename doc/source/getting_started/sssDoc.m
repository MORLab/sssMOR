%% sss - Sparse State-Space Objects
%
% To overcome the size limitations of |ss| and |dss| objects, we have introduced
% a toolbox called *sss* which stands for sparse state-space. 
%
% This toolbox contains an *sss* class that allows the definition of *sss* objects, 
% i.e. state-space objects defined by sparse matrices. 
%
% In fact, the large-scale dynamical system we could not define in the example above 
% can now find place in your workspace
%
A = speye(10^8); b(10^8,1)=1;
sys = sss(A,b,b'); whos sys
clear
%%	
% Further, many functions control engineers use on a daily basis to analyze 
% and manipulate dynamic system objects, such as |bode|, |step|, |pzmap|, 
% |eig|, |isstable| etc. are included in the *sss* toolbox to work with *sss* 
% objects and exploit sparsity, avoiding dense computations, whenever this 
% is possible.

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
