function A = readMatrixMarket(filename)
% Imports a Matrix stored in the Matrix Market format
% ------------------------------------------------------------------
% A = readMatrixMarket(filename)
% Input:        * filename: name of file containing the data
% Output:       * A: matrix
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer
% Last Change:  08 Feb 2011
% ------------------------------------------------------------------

% open file
file = fopen(filename,'r');

while 1
    n=fgetl(file);
    % ignore comment lines containing '%'
    if ~isempty(strfind(n,'%'))
        continue
    end
    n=str2num(n); %#ok<ST2NM>
    if isempty(n)
        continue
    end
    if length(n)>1
        n1=n(1); % number of rows
        n2=n(2); % number of columns
        % number of elements
        if length(n)>2
            n3=n(3);
        else
            n3=n1;
        end
    end
	break
end

% read content from file
M = fscanf(file, '%f ',[3 n3]);

% close file
fclose(file);

% see if matrix is empty
if size(M,1)==0
    if isempty(n1)
        A=[];
        return
    end
    A = sparse([],[],[],n1,n2,0);
    return
end

% merge elements to sparse matrix
A = sparse(M(1,:), M(2,:), M(3,:), n1, n2);

% if matrix is triangular, mirror it
if (all(M(1,:)>=M(2,:)) || all(M(1,:)<=M(2,:))) && n1==n2
    x=stqd('title','Mirror matrix?','string','The matrix is triangular. Do you want to mirror all off-diagonal elements?');
    if ~strcmp(x, 'No')
        A = A + transpose(A) - diag(diag(A));
    end
end
