function [sysr, sysr0, Hinf, HinfRatio, tOpt , bound] = HinfMor(sys, n, varargin) 
    % HINFMOR - H-infinity reduction by tangential interpolation
    % ------------------------------------------------------------------
    %
    % [sysr, Hinf] = HINFMOR(sys, varargin) 
    % Inputs:       * sys: full oder model (sss)
    %               * Opts (opt.) structure with execution parameters
    % Outputs:      * sysr: reduced order model (sss)
    %               * Hinf: H-infinity error (or indicator)
    % ------------------------------------------------------------------
    % USAGE:  
    %
    % See also IRKA
    %
    % ------------------------------------------------------------------
    % REFERENCES:
    % [1] Gugercin (2008), H2 model reduction for large-scale linear
    %     dynamical systems
    % [2] Flagg (2013), Interpolatory Hinf model reduction
    % [3] Beattie (2014), Model reduction by rational interpolation
    % ------------------------------------------------------------------
    % Authors:      Alessandro Castagnotto
    % Last Change:  23 Nov 2015
    % Copyright (c) ?
    % ------------------------------------------------------------------

    %%  Run IRKA
    if sys.isSiso
        % initialize
        try s0 = -eigs(sys,n,'sm').'; catch , s0 = zeros(1,n); end
        % run IRKA
        [sysr0, ~, ~, ~, ~, ~, ~, B_, Rt, C_, Lt] = irka(sys,s0);
    else %MIMO
        % initialize
        %   compute one step of tangential Krylov at 0 to get initial tangent 
        %   directions
        s0 = -eigs(sys,n,'sm').'; Rt = ones(sys.m,n); Lt = ones(sys.p,n);
        sysr = rk(sys,s0,s0,Rt,Lt);  [X,D,Y] = eig(sysr);
        Rt = full((Y.'*sysr.B).'); Lt = full(sysr.C*X); s0 = -diag(D).';
        %run IRKA
        [sysr0, ~, ~, ~, ~, ~, ~, B_, Rt, C_, Lt] = irka(sys,s0,Rt,Lt);
    end

    %   Transform (A- s0*E) to (s0*E- A)
    sysr0.C = -sysr0.C; sysr0.B = -sysr0.B;
    Rt = -Rt; Lt = -Lt;

    % % Check that the generalized tangential directions are correct
    % R = getSylvester(sys,sysr0,-V); L = getSylvester(sys,sysr0,-W,'W'); 
    % if norm(Rt-R)> 1e-5 || norm(Lt-L)> 1e-5
    %     warning('Residuals could be wrong')
    %     keyboard
    % end


    %%  Make Hinf correction
    %
    % steadyState: just take the current steady state error (can yield worse
    %               results)
    % steadyStateOpt: optimize over the difference in between the magnitude of
    %                 the steady state error and Dr
    % findGe0match: finds the feedthrough that best matches the steady-state
    %               error amplitude
    % normOpt:      Optimizes the actual inf norm of the error system!!
    % steadyState+normOpt: Initializes optimization at Ge(0)
    % DRange:       computes the actual error norm for a series of feedtrhoughs
    %               and takes the one with minimum value

    if nargin>2
        corrType = varargin{1};
    else
        corrType = 'normOpt';
    end

    %   Parametrize reduced order models by Dr
    sysrfun = @(Dr) sss(sysr0.A+Lt.'*Dr*Rt, sysr0.B+Lt.'*Dr, ...
                                     sysr0.C+Dr*Rt, Dr, sysr0.E);

    switch corrType
        case 'steadyState'
            warning(['This approach fails due to the same reasons stated for ',...
                'the forAll case']);
            %plus, Ge(0) changes depending on Dr

            Dr = freqresp(sys,0)-freqresp(sysr0,0);
            sysr = sysrfun(Dr);
        case 'steadyStateOpt'
            warning(['This approach fails due to the same reasons stated for ',...
                'the forAll case']);
            % however, is seems to work fine for build...

            G0 = freqresp(sys,0); %the only costly part
            Dr0 = G0-freqresp(sysr0,0);

            cost = @(Dr) abs(...
                        abs(G0-...
                        freqresp(sysrfun(Dr),0)) - abs(Dr));
            DrOpt = fmincon(cost,Dr0)
            sysr = sss(sysr0.A+Lt.'*DrOpt*Rt, sysr0.B-Lt.'*DrOpt, ...
                       sysr0.C-DrOpt*Rt, Dr0, sysr0.E);
        case 'findGe0match'
            warning(['This approach fails since there does not seem to be a ',...
                'reduced order model that yields a completely flat magnitude',...
                'response in the error']);

            G0 = freqresp(sys,0); %the only costly part

            Dr0 = G0 - freqresp(sysr0,0) %get an initial feedthrough
            deltaDr = 100*abs(Dr0), nStep = 5000; DrSet = []; dSet = [];
            syse0 = sys-sysr0;
            Hinf0 = norm(syse0,Inf); d0 = abs(Dr0); dMin = d0; %initial error
            drawnow
            for Dr = Dr0-deltaDr:(2*deltaDr)/nStep:Dr0+deltaDr;
                sysr = sysrfun(Dr);

                d = abs(abs(G0-freqresp(sysr,0)) - abs(Dr));
                if d < dMin
                    dMin = d;
                    dSet = [dSet, d];
                    DrSet = [DrSet,Dr];
    %                 syse = sys-sysr; sigma(syse,'Color',rand(1,3));
    %                 drawnow, keyboard
                end
            end
            %   best feedthrough in term of minimizing the error between the
            %   response at 0 and Inf
            Dr = DrSet(end), dMin
            sysr = sysrfun(Dr);
        case 'normOpt'
            Dr0 = zeros(sys.p,sys.m);
            [sysr, ~, Hinf,tOpt] = normOpt(Dr0);
        case 'steadyState+normOpt'
            % execution params
            plotCostOverDr = 1;
            
            % initialization at steady-state error amplitude response
            G0 = freqresp(sys,0); Dr0 = G0-freqresp(sysr0,0);
            [sysr, DrOpt, Hinf,tOpt] = normOpt(Dr0);

            %   See how the cost behaves around the chosen minimum?
            if plotCostOverDr
                nStep = 20; kRange = 6;
                if isscalar(DrOpt) %SISO
                    DrRange = linspace(DrOpt*(1-kRange),DrOpt*(1+kRange), nStep); 
                elseif any (size(DrOpt) == 1) % SIMO or MISO
                    DrRange = cell(size(DrOpt)); %get the right size
                    DrRange1 = linspace(DrOpt(1)*(1-kRange),DrOpt(1)*(1+kRange), nStep);
                    DrRange2 = linspace(DrOpt(2)*(1-kRange),DrOpt(2)*(1+kRange), nStep);
                    [DrRange1, DrRange2] = meshgrid(DrRange1,DrRange2);
                    DrRange{1} = DrRange1; DrRange{2} = DrRange2;
                end
                plotOverDrRange(DrOpt,Hinf);
            end
        case 'DrRange'
            % Get steady state response of the error system
            G0 = freqresp(sys,0); 
            Dr0 = abs(G0 - freqresp(sysr0,0));

            % Define a range for the feedthrough
            nStep = 100;
            DrRange = linspace(-Dr0,Dr0, nStep); 

            % Run the actual function
            plotOverDrRange;
        otherwise
            error('Specified Hinf optimization type not valid');
    end

    %% Additional outputs
    if nargout > 3
        Hinf = norm(sys-sysr,inf)/norm(sys,inf);
        if nargout > 4
            bound = HinfBound(sys,B_,C_);
        end
    end

    %% Auxiliary
    function [minDr, minVal] = plotOverDrRange(varargin)
        if isnumeric(DrRange) %vector of Dr values
            %initializing
            normVec = zeros(1,length(DrRange));
            normO = norm(ss(sys),inf); 
            normVec(1) = norm(ss(sys-sysr0),inf)/normO; 
            minVal = normVec(1); minDr = 0; DrRange = [0,DrRange];
            for iDr = 2:length(DrRange)
                sysrTest = sysrfun(DrRange(iDr));
                normVec(iDr) = norm(ss(sys-sysrTest),inf)/normO;
                if normVec(iDr) < minVal
                    minDr = DrRange(iDr); minVal = normVec(iDr);
                end
            end
            % sort Dr to be sure you cane use lines in plots
            [DrRange,idx] = sort(DrRange,'ascend'); normVec = normVec(idx);
            figure; h(1) = plot(DrRange,normVec,'-b'); hold on; 
            h(2) = plot(minDr,minVal,'*g');plot([DrRange(1),DrRange(end)],minVal*[1,1],'--g');
            legNames = {'Error', 'Min on plot'};
            %optimizer passed to function
            if ~isempty(varargin) 
                h(3) = plot(varargin{1},varargin{2}/normO,'or'); 
                plot([DrRange(1),DrRange(end)],varargin{2}/normO*[1,1],'--r');
                legNames = [legNames, {'Min from optimization'}];
            end
            
            % Add zero and Dr0 to the plot
            currDr = 0; currVal = norm(ss(sys-sysrfun(currDr)),inf)/normO;
            h(4) = plot(currDr, currVal,'ok');          
            currDr = freqresp(sys-sysrfun(0),0); currVal = norm(ss(sys-sysrfun(currDr)),inf)/normO;
            h(5) = plot(currDr,currVal,'sm');
            legNames = [legNames, {'val@0','Ge0(0)'}];
            
            % labeling
            xlabel('Dr'); ylabel('Relative Hinf error over Dr')
            title(sprintf('%s, n=%i',sys.Name,sysr.n),'Interpreter','none');
            
            % rescale y-axis (often you see high peaks)
            ylim(minVal*[.25,4]);
            
            legend(h,legNames,'Location','SouthOutside');
        else % Dr is a grid (MISO or SIMO)
            currDr = zeros(size(DrRange)); %get the right shape
            minDr = zeros(size(DrRange));
            
            Dr1Range = DrRange{1}; Dr2Range = DrRange{2};
            normMat = zeros(size(Dr1Range));
            normO = norm(ss(sys),inf); 
            % in this case, we don't add 0 as a value to the grid
            minVal = Inf;
            for iDr = 1:size(Dr1Range,1)
                for jDr = 1:size(Dr1Range,2)
                    currDr(1) = Dr1Range(iDr,jDr);
                    currDr(2) = Dr2Range(iDr,jDr);
                    sysrTest = sysrfun(currDr);
                    normMat(iDr,jDr) = norm(ss(sys-sysrTest),inf)/normO;
                    if normMat(iDr,jDr) < minVal
                        minDr(1) = Dr1Range(iDr,jDr);
                        minDr(2) = Dr2Range(iDr,jDr);
                        minVal = normMat(iDr,jDr);
                    end
                end
            end
            
            figure; h(1) = surf(Dr1Range,Dr2Range,normMat); hold on; 
            h(2) = plot3(minDr(1),minDr(2),minVal,'pg','MarkerFaceColor','g');
            plot3(Dr1Range,Dr2Range,minVal*ones(size(Dr1Range)),'--g');
            legNames = {'Error surface', 'Min on surface'};
            %optimizer passed to function
            if ~isempty(varargin) 
                DrOpti = varargin{1}; optimVal = varargin{2}/normO;
                h(3) = plot3(DrOpti(1),DrOpti(2),optimVal,'or','MarkerFaceColor','r'); 
                plot3(Dr1Range,Dr2Range,optimVal*ones(size(Dr1Range)),'--r');
                legNames = [legNames, {'Min from optimization'}];
            end
            
            % Add zero and Dr0 to the plot
            currDr = zeros(size(currDr)); currVal = norm(ss(sys-sysrfun(currDr)),inf)/normO;
            h(4) = plot3(currDr(1),currDr(2),currVal,'ok','MarkerFaceColor','k');            
            currDr = abs(freqresp(sys,0) - freqresp(sysr0,0));
            currVal = norm(ss(sys-sysrfun(currDr)),inf)/normO;
            h(5) = plot3(currDr(1),currDr(2),currVal,'sm','MarkerFaceColor','m');
            legNames = [legNames, {'val@0','Ge0(0)'}];
            
            % labeling
            xlabel('Dr'); zlabel('Relative Hinf error over Dr')
            title(sprintf('%s, n=%i',sys.Name,sysr.n),'Interpreter','none');
            
            % rescale y-axis (often you see high peaks)
            ZLims = minVal*[.25,4];
            zlim(ZLims);
            caxis(ZLims)
            
            % legend
            legend(h,legNames,'Location','SouthOutside')
        end
                
            drawnow
            
            if 1 %saving
                cdir = pwd;
                cd('..\res')
                sysSaveNm = sys.Name;
                k = strfind(sysSaveNm,':');
                if ~isempty(k)
                    if strcmp(sysSaveNm(k-1),'_')
                        sysSaveNm = strrep(sysSaveNm,':','simo');
                    else
                        sysSaveNm = strrep(sysSaveNm,':','miso');
                    end
                end
                saveas(gcf,sprintf('ErrorOverDr_%s_n%i',sysSaveNm,sysr.n));
                cd(cdir);
            end
    end

    function [sysr, DrOpt, Hinf,tOpt] = normOpt(Dr0)
        
            cost = @(Dr) norm(sys-sysrfun(Dr),Inf);
            solver = 'fminsearch';
            warning('optimizing over the actual error norm');
            warning off
            switch solver
                case 'fminsearch'
                %zero initialization
                tic, [DrOpt, Hinf] = fminsearch(cost,Dr0), tOpt = toc

                case 'ga'
                    options = gaoptimset('Display','iter','TimeLimit',5*60,...
                        'UseParallel',true, 'PopInitRange',[-1;1]);
                    tic, [DrOpt, Hinf] = ga(cost,1,[],[],[],[],[],[],[],options); tOpt = toc
            end
            warning on
            sysr = sysrfun(DrOpt); 
    end
    function bound = HinfBound(sys,B_,C_)
        % Panzer 2014 
%         function bndHinf = BoundHinf(L_S,P_S,B,C)
        % Upper bound on H-infinity norm of strictly dissipative system
        % Input: L_S,P_S: Cholesky factor of S=-A-A', and permutation matrix;
        % B_,C_ : Perp input and output matrices
        % Output: bndHinf: Upper bound
        % (c) 2014 Heiko K.F. Panzer, Tech. Univ. Muenchen.
        % This file is published under the BSD 3-Clause License. All rights reserved.
        
        S = -(sys.A + sys.A.'); [L_S,p,P_S] = chol(S); 
            if p
                warning('System is not in strictly dissipative form. HinfBound set to Inf');
                bound = inf;
            else
                B_S = L_S'\(P_S'*B_);
                C_S = (L_S'\(P_S'*C_'))';
                bound = norm(full(C_S*B_S)) + norm(full(B_S))*norm(full(C_S));
            end
    end

end




