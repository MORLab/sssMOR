function [SPC,SPB] = projInd2stokes(sys)
  
% structure:
% | E12 |  0  |   | vp |   | A11 | A12 |   | v |   | B1 |
% |-----+-----| * |----| = |-----+-----| * |---| + |----| * u
% |  0  |  0  |   | pp |   | A21 |  0  |   | p |   | B2 |
%
%                   | v |   v ... velocity
% y = | C1 | C2 | * |---|
%                   | p |   p ... pressure
%
% Caution: A21*E11^(-1)*A12 has to be regular!


    %% Determine the dynamic order nDyn
    nDyn = find(any(sys.E)==0,1)-1;
    
    % works currently only with E11 as identity!
    E11 = partition(sys.E,nDyn);
    if norm(E11-speye(nDyn),'fro')>eps
        error('This works only for E11 = I so far')
    end

    %% generate function handles for (spectral) projection
    SPC = @(C) SpecProjFinC(C,sys,nDyn);
    SPB = @(B) SpecProjFinB(B,sys,nDyn);

end


%% define spectral projections
% multiplication of B from the left with the left spectral projector onto
% the deflating subspace of lambda*E-A corresponding to the finite
% eigenvalues
function PB = SpecProjFinB(B,sys,nDyn)
    % SPECIAL CASE: E11 = I
    % M = I-A12*(A21*A12)^(-1)*A21
    %               | M | -M*A11*A12*(A21*A12)^(-1) |
    % Specleftfin = |---+---------------------------| 
    %               | 0 |             0             | 

    % get blocks of A
    [A11,A12,A21,~] = partition(sys.A,nDyn);

    % partition B
    B1 = B(1:nDyn,:);
    B2 = B(nDyn+1:end,:);

    % compute product of B and projector (multiplication with P from the left)
    PB1 = B1 - A12*((A21*A12)\(A21*B1)) ...
          - A11*A12*((A21*A12)\B2) ...
          + A12*((A21*A12)\(A21*A11*A12*((A21*A12)\B2)));
    PB2 = zeros(size(B2));

    % return projected B
    PB = [PB1; PB2];
end

% multiplication of C from the right with the right spectral projector onto
% the deflating subspace of lambda*E-A corresponding to the finite
% eigenvalues
function CP = SpecProjFinC(C,sys,nDyn)
    % SPECIAL CASE: E11 = I
    % M = I-A12*(A21*A12)^(-1)*A21
    %                |             M             | 0 |
    % Specrightfin = |---------------------------+---| 
    %                | -(A21*A12)^(-1)*A21*A11*M | 0 | 
    
    % get blocks of A
    [A11,A12,A21,~] = partition(sys.A,nDyn);

    % partition C
    C1 = C(:,1:nDyn);
    C2 = C(:,nDyn+1:end);

    % compute product of C and projector (multiplication with P from the right)
    CP1 = C1 - ((C1*A12)/(A21*A12))*A21 ...
          - (C2/(A21*A12))*A21*A11 ...
          + (((C2/(A21*A12))*A21*A11*A12)/(A21*A12))*A21;
    CP2 = zeros(size(C2));

    % return projected C
    CP = [CP1, CP2];
end