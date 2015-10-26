%% arnoldi
% Arnoldi algorithm using multiple expansion points
%
%% Syntax
%
% <html>
%        <div class="syntax">
% [V,Ct]        = arnoldi(E,A,b,s0,IP) <br>
% [V,Ct,W,Bt]   = arnoldi(E,A,b,c,s0,IP) <br>
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
%            E/A/b/c
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            System matrices
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
%            Vector of expansion points
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
%            (opt.) function handle for inner product
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
%            V
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            Orthonormal basis spanning the input Krylov subsp.
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            Ct
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            Right tangential directions of Sylvester Eq.
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            W
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            Orthonormal basis spanning the output Krylov subsp.
%        </p>
%    </td>
%    </tr>
%    <tr bgcolor="#F2F2F2">
%    <td>
%        <p class="table">
%            Bt
%        </p>
%    </td>
%    <td>
%        <p class="table">
%            Left tangential directions of Sylvester Eq.
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
% This function is used to compute the matrix V spanning the
% input Krylov subspace corresponding to E, A, b and s0 [1,2].
%
% The columns of V build an orthonormal basis of the input Krylov
% subspace. The orthogonalization is conducted using a reorthogonalized
% modified Gram-Schmidt procedure [3] with respect to the inner product
% defined in IP (optional). If no inner product is specified, then the
% elliptic product corresponding to E is chosen by default:
%
% IP=@(x,y) (x'*E*y)
%
% which requires E to be a positive definite matrix.
%
%
%% See also
% <rkhelp.html |rk|>
%
%
%% References
% * *[1] Grimme (1997)*, Krylov projection methods for model reduction
% * *[2] Antoulas (2005)*, Approximation of large-scale dynamical systems
% * *[3] Giraud (2005)*, The loss of orthogonality in the Gram-Schmidt...
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
