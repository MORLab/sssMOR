%% tbr
% Performs model order reduction by the Truncated Balanced Realization
%
%% Syntax
%
% <html>
%        <div class="syntax">
% [sysr, varargout] = tbr(sys, varargin) <br>
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
%            q
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            (opt.) order of reduced system
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
%            (opt.) projection matrices (only if q is given!)
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            hsv
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            Hankel singular values
%        </p>
%    </td>
%    </tr>
% </tbody>
% </table>
% </html>
%
%
% <html>
%        <div class="info">
%        <b>Note</b>   If no q is given, the balancing transformation and calculation of the
%        Hankel Singular Values is performed without subsequent model reduction .
%        </div>
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
