%% moments
% Returns the moments or Markov parameters of an LTI system
%
%% Syntax
%
% <html>
%        <div class="syntax">
% m = moments(sys, s0, n) <br>
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
%            s0
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            expansion point (inf -> Markov parameters)
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            n
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            number of <tt>moments</tt> to be computed
%        </p>
%    </td>
%    </tr>
% </tbody>
% </table>
% </html>
%
%% Output
%
% <html>
% <table cellspacing="0" cellpadding="4" width="" border="1" frame="box" rules="none" class="">
% <tbody valign="top">
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            m
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            vector of <tt>moments</tt> / Markov parameters
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
