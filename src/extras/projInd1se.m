function [SPC,SPB] = projInd1se(sys)
    
% structure:
% | E11 |  0  |   | x1p |   | A11 | A12 |   | x1 |   | B1 |
% |-----+-----| * |-----| = |-----+-----| * |----| + |----| * u
% |  0  |  0  |   | x2p |   | A21 | A22 |   | x2 |   | B2 |
%
%                   | x1 |
% y = | C1 | C2 | * |----|
%                   | x2 |
%
% Caution: E11 and A22 have to be regular!


    %% Determine the dynamic order nDyn
    nDyn = find(any(sys.E)==0,1)-1;

    %% generate function handles for (spectral) projection
    SPC = @(C) SpecProjFinC(C,sys,nDyn);
    SPB = @(B) SpecProjFinB(B,sys,nDyn);

end


%% define spectral projections
% multiplication of B from the left with the left spectral projector onto
% the deflating subspace of lambda*E-A corresponding to the finite
% eigenvalues
function PB = SpecProjFinB(B,sys,nDyn)
    %               | I | -A12*(A22)^(-1) |
    % Specleftfin = |---+-----------------| 
    %               | 0 |        0        | 

    % get blocks of A
    [~,A12,~,A22] = partition(sys.A,nDyn);
    
    % partition B
    B1 = B(1:nDyn,:);
    B2 = B(nDyn+1:end,:);
    
    % compute product of B and projector (multiplication with P from the left)
    PB1 = B1 - A12*(A22\B2);
    PB2 = zeros(size(B2));
    
    % return projected B
    PB = [PB1; PB2];
end

% multiplication of C from the right with the right spectral projector onto
% the deflating subspace of lambda*E-A corresponding to the finite
% eigenvalues
function CP = SpecProjFinC(C,sys,nDyn)
    %                |       I       | 0 |
    % Specrightfin = |---------------+---| 
    %                | -A22^(-1)*A21 | 0 | 
    
    % get blocks of A
    [~,~,A21,A22] = partition(sys.A,nDyn);
    
    % partition C
    C1 = C(:,1:nDyn);
    C2 = C(:,nDyn+1:end);
    
    % compute product of C and projector (multiplication with P from the right)
    CP1 = C1 - (C2/A22)*A21;
    CP2 = zeros(size(C2));
    
    % return projected C
    CP = [CP1, CP2];
end