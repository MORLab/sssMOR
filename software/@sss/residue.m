function [r,p,d] = residue(sys)
% Computes residues, poles and feedthrough of an LTI system
% ------------------------------------------------------------------
% [r,p,d] = residue(sys)
% Inputs:       * sys: an sss-object containing the LTI system
% Outputs:      * residues r, poles p and feedthrough d, such that
%                   G(s) = r_i/(p_i+s) + d
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer
% Last Change:  10 Nov 2011
% ------------------------------------------------------------------

% are poles and residues already available?
if ~isempty(sys.poles) && ~isempty(sys.residues)
    p=sys.poles;
    r=sys.residues;
    d=sys.d;
    return
end


tic
%perform eigen-decomposition of system
if sys.is_dae
    [T,J] = eig(full(sys.A),full(sys.E));
    if max(max(J))==inf
    error('System contains algebraic states.')
    end
else
    [T,J] = eig(full(sys.A));
end

% transform system to diagonal form
p=conj(diag(J)');
B=(sys.E*T)\sys.B;
C=sys.C*T;
d=sys.D;

% calculate residues
if sys.is_mimo
    warning('Code edited 2013/10/16 HP');
    %mimo
%     r = cell(sys.p, sys.m);
    r = cell(sys.n, 1);
    for i=1:sys.n
%         for i=1:sys.p
%             r{i,j} = C(i,:) .* conj(B(:,j)');
%         end
        r{i} = C(:,i)*B(i,:);
    end
else
    %siso
    r = {C .* conj(B')};
end

return

% measure time required for diagonalization and computation of 
% 1000 values of the impulse response
sys.SimulationTime=toc;
tic
cellfun(@(x) sum((diag(x)*conj(exp(3'*p))'),1),r,'UniformOutput',false);
simtime=toc;

% store results to caller workspace
sys.residues=r;
sys.poles=p;
sys.SimulationTime=sys.SimulationTime+1000*simtime;
if inputname(1)
    assignin('caller', inputname(1), sys);
end
