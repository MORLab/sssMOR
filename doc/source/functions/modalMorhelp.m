%% modalMor
% Modal model order reduction of LTI SISO systems
%
%% Syntax
%
% <html>
%        <div class="syntax">
% [sysr, V, W] = modalMor(sys, q, Opts) <br>
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
%            order of reduced system
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            opts
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            a structure containing following options
%        </p>
%        <p class="table">
% <table cellspacing="0" cellpadding="4" width="" border="1" frame="box" rules="none" class="inner">
% <tbody valign="top">
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            type
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            options to eigs command - &#123;<tt>'SM'</tt>,'LM',...&#125;
%        </p>
%        <p class="table">
% <table cellspacing="0" cellpadding="4" width="" border="1" frame="box" rules="none" class="inner">
% <tbody valign="top">
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            <tt>'SM'</tt>
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            for eigenvalues of Smallest Magnitude
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            <tt>'LM'</tt>
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            for eigenvalues of Largest Magnitude
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            <tt>'SA'</tt>
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            for smallest algebraic
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            <tt>'LA'</tt>
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            for largest algebraic
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            <tt>'SR'</tt>
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            for smallest real part
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            <tt>'LR'</tt>
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            for largest real part
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            scalar number
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            for eigenvalues next to it
%        </p>
%    </td>
%    </tr>
% </tbody>
% </table>
%        </p>
%    </td>
%    </tr>
% </tbody>
% </table>
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
%            projection matrices
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
% This function computes the reduced order system sysr and the
% projection matrices V and W by the modal reduction technique [1].
%
% Only a few eigenvalues and left and right eigenvectors of the pair (A,E)
% are computed with the eigs command. The right eigenvectors build the
% columns of V, while the left eigenvectors build the columns of W.
%
%
%% See also
% <tbrhelp.html |tbr|>, <rkhelp.html |rk|>, <irkahelp.html |irka|>
%
%
%% References
% * *[1] Antoulas (2005)*, Approximation of large-scale dynamical systems
% * *[2] Lehoucq and Sorensen (1996)*, Deflation Techniques for an Implicitly Re-Started Arnoldi Iteration.
% * *[3] Sorensen (1992)*, Implicit Application of Polynomial Filters in a k-Step Arnoldi Method
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
