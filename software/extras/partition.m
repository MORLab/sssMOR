function [A11,A12,A21,A22] = partition(A,row,col)

%   function [A11,A12,A21,A22] = partition(A,row,col)
%
%   Partitions A into 4 submatrices. A11 has the dimension (row x col), all
%   other submatrices accordingly.
%
%   The argument "col" is optional

if ~exist('col','var') || isempty(col)
    col = row;
end

A11 = A(1:row,1:col);
if nargout > 1
    A12 = A(1:row,col+1:end);
    A21 = A(row+1:end,1:col);
    A22 = A(row+1:end,col+1:end);
end