%% References
%
% * *[1] Antoulas (2005)*, Approximation of large-scale dynamical systems
% * *[2] Antoulas (2010)*, Interpolatory model reduction of large-scale...
% * *[3] Beattie et al. (2014)*, Model reduction by rational interpolation
% * *[4] Castagnotto et al. (2015)*, Stability-preserving, adaptive
%      model order reduction of DAEs by Krylov subspace methods
% * *[5] Documentation of the Control System Toolbox from MATLAB*
% * *[6] Eid, Panzer and Lohmann (2009)*, How to choose a single 
%      expansion point in Krylov-based model reduction? Technical 
%      reports on Automatic Control, vol. TRAC-4(2), Institute of 
%      Automatic Control, Technische Universitaet Muenchen.
% * *[7] Eid (2009)*, Time domain Model Reduction by Moment Matching, Ph.D
%      thesis, Institute of Automatic Control, Technische 
%      Universitaet Muenchen.
% * *[8] Foellinger (2013)*, Regelungstechnik (pp. 305-319)
% * *[9] Gallivan et al. (2002)*, Sylvester equations and projection
%      based model reduction
% * *[10] Gear (1971)*, Numerical Initial Value Problems in 
% Ordinary Differential Equations.
% * *[11] Giraud (2005)*, The loss of orthogonality in the Gram-Schmidt... 
% * *[12] Grimme (1997)*, Krylov projection methods for model reduction
% * *[13] Gugercin et al. (2008)*, H2 model reduction for large-scale linear dynamical systems
% * *[14] Gugercin (2008)*, An iterative SVD-Krylov based method for
%      model reduction of large-scale dynamical systems
% * *[15] Lehoucq and Sorensen (1996)*, Deflation Techniques for an 
%      Implicitly Re-Started Arnoldi Iteration.
% * *[16] Moore (1981)*, Principal component analysis in linear systems: controllability,
% observability and model reduction
% * *[17] Panzer (2014)*, Model Order Reduction by Krylov Subspace Methods
%      with Global Error Bounds and Automatic Choice of Parameters
% * *[18] Penzl (2000)*, LYAPACK - A MATLAB Toolbox for Large Lyapunov
% and Riccati Equations, Model Reduction Problems, and
% Linear-Quadratic Optimal Control Problems.
% * *[19] Saak, J. and Koehler, M. and Benner, P. (2016)*, M-M.E.S.S.-1.0 -- 
% The Matrix Equations Sparse Solvers library
% * *[20] Shampine and Gordon (1975)*, Computer Solution of Ordinary Differential 
% Equations: the Initial Value Problem, W. H. Freeman, San Francisco.
% * *[21] Shampine (1994)*, Numerical Solution of Ordinary Differential Equations, 
% Chapman & Hall, New York.
% * *[22] Sorensen (1992)*, Implicit Application of Polynomial Filters 
%      in a k-Step Arnoldi Method.
% * *[23] Wolf (2014)*, H2 Pseudo-Optimal Moder Order Reduction
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
