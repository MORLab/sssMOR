function sysd = diag(sys)
% Transforms an LTI system to (block)diagonal representation
% ------------------------------------------------------------------
% [mag, phase, omega] = bode(sys, omega, in, out, options)
% Input:        * sys: an sss-object containing the LTI system
% Output:       * sysd: system in diagonal form
% During the diagonalization, the C-vector is normalized to contain ones
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de)
% Last Change:  03 Feb 2011
% ------------------------------------------------------------------

%perform eigen-decomposition of system
if sys.is_dae
    [T,A] = eig(full(sys.A),full(sys.E));
    if max(max(A))==inf
        error('System contains algebraic states.')
    end
else
    [T,A] = eig(full(sys.A));
end

% transform system to diagonal form
B=(sys.E*T)\sys.B;
C=sys.C*T;

% save poles and residues
sys.poles = transpose(diag(A));
if sys.is_mimo
    %mimo
    r = cell(sys.p, sys.m);
    for j=1:sys.m
        for i=1:sys.p
            r{i,j} = C(i,:) .* conj(B(:,j)');
        end
    end
else
    %siso
    r = {C .* conj(B')};
end
sys.residues = r;

% find real system representation
i=0;
while i<length(B)
    i=i+1;
    if i<length(B)
        if abs(real(A(i,i)+A(i+1,i+1))) >= abs(imag(A(i,i)+A(i+1,i+1)))*10^3 && ...
           abs(imag(A(i,i)-A(i+1,i+1)))/abs(real(A(i,i)+A(i+1,i+1))) >= 10^(-3)
            delta = real(A(i,i) + A(i+1,i+1))/2;
            omega = imag(A(i,i) - A(i+1,i+1))/2;
            A(i:i+1,i:i+1) = [delta, omega; -omega, delta];
            
            % residues are shifted to input vector
            pc = real(B(i)*C(i)+B(i+1)*C(i+1))/2;
            pd = imag(B(i)*C(i)-B(i+1)*C(i+1))/2;
            B(i) = pc+pd;
            B(i+1) = pc-pd;
            i=i+1;
            continue
        end
    end
    B(i) = B(i) * C(i);
end

% remove remaining imaginary components (resulting from numerical noise)
B=real(B);

sysd = sss(A, B, ones(size(C)), sys.D);

if inputname(1)
    assignin('caller', inputname(1), sys);
end
