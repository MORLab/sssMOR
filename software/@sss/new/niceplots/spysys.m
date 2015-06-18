function fh = spysys(A,E)

% spysys(A,E)
% spysys(sys)

%   parse input
if strcmp(class(A),'sss')
    sys = A;
    A = sys.a;
    E = sys.e;
end

fh = nicefigure('Sparsity pattern of system matrices');

% subplot(1,2,1); spy(A,'Color',TUM_Blau); title('A');
% subplot(1,2,2); spy(E,'Color',TUM_Blau); title('E');

subplot(1,2,1); spy(A); title('A');
subplot(1,2,2); spy(E); title('E');

