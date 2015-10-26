%% irka
% Iterative Rational Krylov Algorithm
%
%% Syntax
%
% <html>
%        <div class="syntax">
% [sysr, V, W, s0, s0_traj] = irka(sys, s0, Opts) <br>
%        </div>
% </html>
%
%% Inputs
%
% <html>
% <table cellspacing="0" cellpadding="4" width="" border="1" frame="box" rules="none" class="">
% <tbody valign="top">
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            sys
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            full oder model (sss)
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            s0
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            vector of initial shifts
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            Opts
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            (opt.) structure with execution parameters
%        </p>
%    </td>
%    </tr>
% </tbody>
% </table>
% </html>
%
%% Outputs
%
% <html>
% <table cellspacing="0" cellpadding="4" width="" border="1" frame="box" rules="none" class="">
% <tbody valign="top">
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            sysr
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            reduced order model (sss)
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            V<br>W
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            resulting projection matrices
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            s0
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            final choice of shifts
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            s0_traj
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            trajectory of all shifst for all iterations
%        </p>
%    </td>
%    </tr>
% </tbody>
% </table>
% </html>
%
%% Examples
% No examples
%
%
%% Description
% This function executes the Iterative Rational Krylov
% Algorithm (|irka|) as proposed by Gugergin and Beattie in [1].
%
% The |irka| iteration is conducted to search for an optimal set of
% shifts in Krylov subspace-based model reduction. If |irka| converges,
% then the reduced model is known to be a local optimum with respect
% to the H2 norm of the error.
%
%
%% See also
% <arnoldihelp.html |arnoldi|>, <rkhelp.html |rk|>
%
%
%% References
% * *[1] Gugercin (2008)*, H2 model reduction for large-scale linear dynamical systems
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
