%% rk
% Model Order Reduction by Rational Krylov
%
%% Syntax
%
% <html>
%        <div class="syntax">
% [sysr, V, W] = rk(sys, s0_inp, s0_out, IP) <br>
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
%            an sss-object containing the LTI system
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            s0_inp
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            Expansion points for Input Krylov Subspace
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            s0_out
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            Expansion points for Output Krylov Subspace
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            IP
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            Inner product (optional)
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
%            reduced system
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
%            Projection matrices spanning Krylov subspaces
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
% s0 may either be horizontal vectors containing the desired
% expansion points, e.g. [1 2 3] matches one moment about 1, 2 and 3,
% respectively. [1+1j 1-1j 5 5 5 5 inf inf] matches one moment about 1+1j,
% 1-1j, 4 moments about 5 and 2 Ma|rk|ov parameters.
%
% An alternative notation for s0 is a two-row matrix, containing the
% expansion points in the first row and their multiplicity in the second,
% e.g. [4 pi inf; 1 20 10] matches one moment about 4, 20 moments about pi
% and 10 Ma|rk|ov parameters.
%
% To perform one-sided |rk|, set s0_inp or s0_out to [], respectively.
%
%
%% See also
% <arnoldihelp.html |arnoldi|>, <rkhelp.html |rk|>
%
%
%% References
% * *[1] Grimme (1997)*, Krylov projection methods for model reduction
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
