function [sysm, s0mTot, V, W,nLU] = modelFct(sys,s0m,s0mTot,V,W,Opts)
% MODELFCT - Computes or updates the model function of an sss object
%
% Syntax:
%       sysm = MODELFCT(sys,s0m)
%       [sysm, s0mTot, V, W]        = MODELFCT(sys,s0m)
%       [sysm, s0mTot, V, W]        = MODELFCT(sys,s0m,s0mTot,V,W)
%       [sysm, s0mTot, V, W]        = MODELFCT(sys,s0m,s0mTot,V,W,Opts)
%       [sysm, s0mTot, V, W, nLU]   = MODELFCT(sys,s0m,...)
%
% Description:
%       This function generates a surrogate model, called "model function" |sysm|,
%       of a high-dimensional model |sys| by Hermite interpolation at the complex
%       frequencies defined in the vector |s0m|.
%
%       If additional inputs |s0mTot|, |V|, |W| are passed, then the model function
%       is updated combining the already existing information, defined by
%       |s0mTot|, to |s0m|. |V| and |W| are updated respectively.
%
%       The optionl structure |Opts| allows the definition of additional
%       execution parameters.
%
%       This function is used in the model function based MOR approach
%       introduced in [1] and [2].
%
%       //Note: In its current form, MODELFCT supports only SISO models.
%       An extension to MIMO will be given in a later release.
% 
% Input Arguments:  
%       *Required Input Arguments:*
%       -sys:			full oder model (sss)
%       -s0m:			vector of shifts for model function (update)
%
%       *Optional Input Arguments:*
%       -s0mTot:        vector of all shifts already used
%       -V,W:           corresponding projection matrices V,W
%       -Opts:			structure with execution parameters
%		  -.updateModel:chooses which shifts of s0m are included in update;
%                       [{'new'} / 'all' / 'lean' ]
%            -.modelTol:tolerance for identifying new shifts;
%                       [{1e-2} / positive float ]
%            -.plot     : generate analysis plots;
%						[{false} / true]
%            -.tol:		convergence tolerance;
%						[{1e-3} / positive float]
%
% Output Arguments:      
%       -sysm:          (updated) model function;
%       -s0mTot:        cumulated vector of approximation frequencies
%       -V,W:           (updated) projection matrices
%       -nLU:           number of LU decompositions
%
% See also:
%       cirka, modelFctMor, solveLse, rk, spark
%
% References:
%       * *[1] Panzer (2014)*, Model Order Reduction by Krylov Subspace Methods
%              with Global Error Bounds and Automatic Choice of Parameters
%       * *[2] Castagnotto et al. (2016)*, Fast H2 optimal model order
%              reduction exploiting the local nature of Krylov-Subspace...
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  09 Aug 2017
% Copyright (c) 2016-2017 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

if ~sys.isSiso, error('sssMOR:modelFct:notSiso','This function currently works only for SISO models');end

    %%  Define default execution parameters
    Def.updateModel = 'new'; % 'all','new','lean'
    Def.modelTol    = 1e-2;
    Def.plot        = false;

    if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
    else
        Opts = parseOpts(Opts,Def);
    end  

    %% Parse input
    if nargin < 3,  %new model function
        s0mTot = s0m;
        V = []; W = []; 
    else %model function update
        s0m = updateModelFctShifts(s0mTot,s0m,Opts);           
        s0mTot = [s0mTot, s0m];
    end
    
    %% Initialize variables in nested functions
    N   = size(sys.A,1);
    L1  = sparse(N,N);U1=L1;P1=L1;Q1=L1; 
    L2  = sparse(N,N);U2=L2;P2=L2;Q2=L2; 
    
    %%  Compute the model function
    % Check model function is not larger than original
    if length(s0mTot)<size(sys.a,1)
        %   Update model function
        [sysm,V,W,nLU] = updateModelFct(s0m,V,W);
    else
        warning('sssMOR:modelFct:sizeLimit',...
            ['Model function is already as big as the original.',...
            ' Returning the original model.']);
        sysm    = sys;
        V       = speye(sys.n);
        W       = V;
    end
    
%%  Auxiliary functions --------------------------------------------------
    %%  Shift and model function update
    function s0m = updateModelFctShifts(s0mTot,s0new,Opts)
        switch Opts.updateModel
            case 'all'
                s0m = s0new;
            case 'new'
                idx = ismemberf2(s0new,s0mTot,Opts.modelTol); 
                s0m = s0new(~idx);
                if Opts.plot
                    fh = figure; lh(1) = plot(complex(s0mTot),'xb'); hold on
                    lh(2) = plot(complex(s0new),'or');
                    axis equal
                    for iS = 1:length(s0mTot);
                        [xp,yp] = circle(real(s0mTot(iS)),imag(s0mTot(iS)),Opts.modelTol*abs(s0mTot(iS)));
                        lh(3) = plot(xp,yp,'g');
                    end
                    legEntries = {'old','new','tolerance'};
                    if ~isempty(s0m)
                        lh(4) = plot(complex(s0m),'k+');
                        legEntries = [legEntries,'added'];
                    end
                    legend(lh,legEntries);
%                     set(gca,'xscale','log');
                    xlabel('Re');ylabel('Im'); title('New shifts for model function')
                    pause
                    close(fh)
                end
            otherwise
                error('selected model function update is not valid');
        end
    end
    function [sysm,V,W,nLU] = updateModelFct(s0,V,W)
        if isempty(V)
            %first run
            [sysm,V,W,~,~,~,~,~,~,nLU] = rk(sys,s0,s0);
        else
            idxComplex = find(imag(s0));
            if ~isempty(idxComplex)
                s0c = cplxpair(s0(idxComplex));
                s0(idxComplex) = []; %real shifts
                s0 = [s0 s0c(1:2:end)]; %add 1 complex shift per complex partner
            end
            nLU = 0;
            for iShift = 1:length(s0)
                if iShift > 1 && s0(iShift)==s0(iShift-1)
                        %do nothing: no new LU decomposition needed
                else %new LU needed
                    computeLU(s0(iShift));
                    nLU = nLU+1;
                end
                V = newColV(V);  W = newColW(W);
            end
            [V,~] = qr(V,0); [W,~] = qr(W,0);
%             sysm = sss(W'*sys.A*V,W'*sys.B,sys.C*V,...
%                             zeros(size(sys.C,1),size(sys.B,2)),W'*sys.E*V);
            sysm = projectiveMor(sys,V,W);
        end   
        %%  Storing additional parameters
        %Stroring additional information about thr reduction in the object 
        %containing the reduced model:
        %   1. Define a new field for the Opts struct and write the information
        %      that should be stored to this field
        %   2. Adapt the method "parseParamsStruct" of the class "ssRed" in such a
        %      way that the new defined field passes the check
        Opts.originalOrder  = sys.n;
        Opts.s0mTot         = s0mTot;
        sysm = ssRed(sysm.A,sysm.B,sysm.C,sysm.D,sysm.E,'modelFct',Opts,sys);
    end
    %%  Functions copied from "spark"
    function computeLU(s0)
        % compute new LU decompositions
        if imag(s0)  % complex conjugated 
            [L1,U1,P1,Q1] = lu(sparse(sys.A-s0*sys.E));  L2=conj(L1);U2=conj(U1);P2=P1;Q2=Q1;
        else % real shift
            [L1,U1,P1,Q1] = lu(sparse(sys.A-s0*sys.E)); L2 =[];%make empty
        end
    end
    function V = newColV(V)
        % add columns to input Krylov subspace
        iCol=size(V,2)+1;
        if iCol==1, x = sys.B; else x=sys.E*V(:,iCol-1);end
        r1  = Q1*(U1\(L1\(P1*x)));   
        if isempty(L2) %only one shift 
            v1 = r1;
%             V = GramSchmidt([V,v1],[],[],iCol*[1 1]);
            V = [V, v1];
%             [V,~] = qr([V,v1],0);
        else %complex conjugated pair
            tmp = Q2*(U2\(L2\(P2*x)));
            v1 = real(0.5*r1 + 0.5*tmp); v2  = real(Q2*(U2\(L2\(P2*(sys.E*r1)))));
%             V = GramSchmidt([V,v1,v2],[],[],[iCol,iCol+1]);
%             [V,~] = qr([V,v1,v2],0);
            V = [V, v1, v2];
        end
    end
    function W = newColW(W)
        % add columns to output Krylov subspace
        iCol=size(W,2)+1;
        if iCol==1, x=sys.C; else x=W(:,iCol-1)'*sys.E; end
        l1  = x*Q1/U1/L1*P1;          
        if isempty(L2) %only one shift
            w1 = l1;
%             W = GramSchmidt([W,w1'],[],[],iCol*[1 1]);
%             [W,~] = qr([W,w1'],0);
            W = [W, w1'];
        else %complex conjugated pair
            tmp = x*Q2/U2/L2*P2;
            w1 = real(0.5*l1 + 0.5*tmp);  w2  = real(l1*sys.E*Q2/U2/L2*P2);
%             W = GramSchmidt([W,w1',w2'],[],[],[iCol,iCol+1]);
%             [W,~] = qr([W,w1',w2'],0);
            W = [W, w1', w2'];
        end
    end
end
