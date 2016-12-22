function varargout = getDesiredOutput(func,outputNo,varargin)
% GETDESIREDOUTPUT - get only the desired output from a function
%
% Syntax:
%       GETDESIREDOUTPUT(func,outputNo)
%       GETDESIREDOUTPUT(func,outputNo,varargin)
%       varargout = GETDESIREDOUTPUT(func,outputNo,varargin)
%
% Description:
%       This function executes a function defined by its handle func and
%       returns only the outputs of said function defined by their indices in
%       the vector outputNo. If the function in func admits inputs, these can
%       be defined as varagins.
%
%       This function is especially useful when defining anonymous functions
%       where only some particular inputs of the base function are required.
%
% See also:
%   VARARGOUT, VARARGIN
%
% References:
%   [1] http://stackoverflow.com/questions/3096281/skipping-outputs-with-anonymous-function-in-matlab
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  02 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

    varargout = cell(max(outputNo),1);
    [varargout{:}] = func(varargin{:});
    varargout = varargout(outputNo);
end