function fh = spysys(varargin)

% spysys(E,A)
% spysys(sys)

%   parse input
if nargin == 2
    E = varagin{1};
    A = varargin{2};
else
    sys = varargin{1};
    A = sys.a;
    E = sys.e;
end

fh = nicefigure('Sparsity pattern of system matrices');

subplot(1,2,1); spy(E); title('E');
subplot(1,2,2); spy(A); title('A');

