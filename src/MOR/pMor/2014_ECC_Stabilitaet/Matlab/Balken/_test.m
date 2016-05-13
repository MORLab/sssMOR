A = rand(3,3);
A=-A'*A;

Q = rand(3,3); Q = orth(Q);
%Q = 2*eye(3,3); Q = orth(Q);


% P = lyap(A',eye(3,3));
% L = chol(P);
% A/L
% 
% An = Q'*A*Q;
% Pn = lyap(An',eye(3,3));
% Ln = chol(Pn);
% An/(L*Q)

P = lyap(A',eye(3,3));
L = chol(P);
A/L

An = Q'*A*Q;
Pn = lyap(An',eye(3,3));
Ln = chol(Pn);
An/(Ln)

% P = lyap(A',eye(3,3));
% L = sqrt(P);
% A/L
% 
% An = Q'*A*Q;
% Pn = lyap(An',eye(3,3));
% Ln = sqrt(Pn);
% An/Ln