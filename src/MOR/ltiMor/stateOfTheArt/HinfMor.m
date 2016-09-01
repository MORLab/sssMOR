function [sysr, HinfErr, sysr0, HinfRatio, tOpt, bound, surrogate, Virka, Rt] = HinfMor(sys, n, varargin) 
    % HINFMOR - H-infinity reduction by tangential interpolation
    % ------------------------------------------------------------------
    % TODO
    % sysr = HINFMOR(sys, varargin) 
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
   
    %%  Input parsing and execution parameters
    if ~isempty(varargin) 
        if isstruct(varargin{1}) %Opts passsed
            Opts = varargin{1};
        else %Opts.corrType passed
            Opts.corrType = varargin{1};        
        end
    else
        Opts = struct();
    end
    
    Def.plotCostOverDr = 0;
    Def.irka        = struct('stopCrit','s0','tol',1e-5,'type','stab');
    Def.corrType    = 'normOptCycle';
    Def.solver      = 'fmincon';    %optimization solver
    Def.algorithm   = 'sqp'; %'interior-point'
    Def.display     = 'off';
    Def.DrInit      = '0';          %0, '0', Ge0, matchGe0, maxGe
    Def.plot        = false;            % generate analysis plot
    Def.sampling    = 'random'; %sampling for sweepDr
    Def.sweepPoints = 5e3;
    
    Def.surrogateError  = true;     %create a surrogate of the error, not of sys
    Def.surrogate       = 'vf';     %original, 'model', 'vf', 'loewner'
    Def.whatData        = 'new';    %'all','new'
    Def.deflate         = 1;
    Def.tol             = 1e-1;     %ismemberf/getModelData
    Def.rankTol         = eps;      %rank tolerance Loewner
    Def.surrTol         = 1e-4;
    
    Def.vf.poles        = 'vectfit3'; %vectfit,eigs, serkan
    Def.vf.maxiter      = 10;
    Def.vf.tol          = 1e-5;
    Def.vf.adaptiveOrder = false; %automatically detect approx. order
    Def.vf.method       = 4;    % 4 is the only method working well    
    Def.wLim            = [1e-1,1e4];
    
    Def.debug           = false;
    
    % create the options structure
    if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
    else
        Opts = parseOpts(Opts,Def);
    end    
    %%  Run IRKA
    if sys.isSiso
        % initialize
        try s0 = -eigs(sys,n,'sm').'; catch , s0 = zeros(1,n); end
%         s0 = zeros(1,n);
        % run IRKA
        [sysr0, Virka, ~, s0opt, rt, lt, ~, ~, Rt, ~, ~, Lt,~, s0Traj,RtTraj, LtTraj] = irka(sys,s0,Opts.irka);
    else %MIMO
        % initialize
        %   compute one step of tangential Krylov at 0 to get initial tangent 
        %   directions
        
        %   Initialize trivially
        s0 = ones(1,n); Rt = ones(sys.m,n); Lt = ones(sys.p,n);
%         s0 = rand(1,n); Rt = rand(sys.m,n); Lt = rand(sys.p,n);

%         s0 = linspace(Opts.wLims(1),Opts.wLims(2),n/2); s0 = [1i*s0,-1i*s0];
%         Rt = ones(sys.m,n/2); Rt = [1i*Rt,-1i*Rt];
%         Lt = ones(sys.p,n/2); Lt = [1i*Lt,-1i*Lt];
% 
%         s0 = logspace(log10(Opts.wLims(1)),log10(Opts.wLims(2)),n);
%         Rt = ones(sys.m,n);
%         Lt = ones(sys.p,n); 
        
%         Rt = rand(sys.m,n); Lt = rand(sys.p,n);
%         s0 = -eigs(sys,n,'sm').'; Rt = ones(sys.m,n); Lt = ones(sys.p,n)
%         sysr = rk(sys,s0,s0,Rt,Lt);  
%         sysr = tbr(sys,n);
%         [X,D,Y] = eig(sysr);
%         Rt = full((Y.'*sysr.B).'); Lt = full(sysr.C*X); s0 = -diag(D).';

        %run IRKA
        [sysr0, Virka, Wirka, s0opt, rt, lt, ~, Sv, Rt, ~, Sw, Lt,~, s0Traj,RtTraj, LtTraj] = irka(sys,s0,Rt,Lt,Opts.irka);
        if Opts.plot; 
            figure; plot(complex(reshape(s0Traj,1,numel(s0Traj))),'x');
            plotName = sprintf('%s_%s_n%i_IRKAshifts',sys.Name,Opts.surrogate,n);
            title(plotName,'Interpreter','none');
            saveas(gcf,fullfile('..','res',sprintf('%s.fig',plotName)));
       end
    end
    
    sysr0 = ss(sysr0);
    
    %   Check Rt and Lt
%     normest(sys.A*Virka-sys.E*Virka*Sv-sys.B*Rt)
%     normest(sys.A.'*Wirka-sys.E.'*Wirka*Sw.'-sys.C.'*Lt)
%     
    % Check interpolation of the D-shifted model
%     keyboard
%  
%     Dr = rand(sys.p,sys.m);
%     sysrD = sysrfun(Dr);   
% %     sysrD = dss(sysr0.A+Lt.'*Dr*Rt, sysr0.B+Lt.'*Dr, sysr0.C+Dr*Rt, Dr, sysr0.E);
%     sysrD = sysr0 + sysrDelta(Dr); 
% 
%     H = @(sys,s) sys.C*((sys.E*s-sys.A)\sys.B) + sys.D;
%     
%     for iS0 = 1:length(s0)
%         (H(sysr0,s0opt(iS0))-H(sysrD,s0opt(iS0)))*rt(:,iS0)
%         lt(:,iS0).'*(H(sysr0,s0opt(iS0))-H(sysrD,s0opt(iS0)))
%     end  
    %%  Create Surrogate Model
    %   To reduce the cost of Hinf optimization, create a surrogate model
    %   from the data collected during irka
    if Opts.surrogateError
        syse0  = ss(sys)-sysr0; 
        syse0m = createSurrogate(syse0, Opts.surrogate);
        fprintf('Size of the surrogate model: %i \n',size(syse0m.a,1))
        fprintf('Stability of surrogate model: %i \n',isstable(syse0m))
        surrogate = syse0m;
    else
        sysm = createSurrogate(sys, Opts.surrogate);
        fprintf('Size of the surrogate model: %i \n',size(sysm.a,1))
        fprintf('Stability of surrogate model: %i \n',isstable(sysm))
        surrogate = sysm;
    end
%     keyboard
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

    switch Opts.corrType
        %{
        case 'steadyState'
            warning('this option is obsolete and shows bad performance');
            %plus, Ge(0) changes depending on Dr

            DrOpt = DrInit('Ge0');
            sysr = sysrfun(DrOpt);
        case 'steadyStateOpt'
            warning('this option is obsolete and shows bad performance');
            % however, is seems to work fine for build...

            G0 = freqresp(sysem,0); %the only costly part
            Dr0 = G0-freqresp(sysr0,0);

            cost = @(Dr) abs(...
                        abs(G0-...
                        freqresp(sysrfun(Dr),0)) - abs(Dr));
            DrOpt = fmincon(cost,Dr0);
            sysr = sss(sysr0.A+Lt.'*DrOpt*Rt, sysr0.B-Lt.'*DrOpt, ...
                       sysr0.C-DrOpt*Rt, Dr0, sysr0.E);
        case 'findGe0match'
            warning('this option is obsolete and shows bad performance');

            DrOpt = DrInit('matchGe0');
            sysr = sysrfun(Dr);
        case 'normOpt_siso'
            warning('this option is obsolete and shows bad performance');
             % running optimization wr to each entry of D individually
            % the cost function takes into account only the SISO tf
            DrOpt = zeros(sys.p,sys.m); tOpt = 0;
            for iOut = 1:sys.p
                for jIn = 1:sys.m
                    Dr0 = 0; 
                    sysrCurr = @(Dr) sss(....
                                    sysr0.A+Lt(iOut,:)'*Dr*Rt(jIn,:),...
                                    sysr0.B(:,jIn)+Lt(iOut,:).'*Dr, ...
                                    sysr0.C(iOut,:)+Dr*Rt(jIn,:), Dr, sysr0.E);
                    cost = @(Dr) norm(sysem(iOut,jIn)-sysrCurr(Dr),Inf);
                    [DrOptCurr, ~, tOptCurr] = normOpt(Dr0,cost);
                    tOpt = tOpt + tOptCurr;
                    DrOpt(iOut,jIn) = DrOptCurr;
                end
            end
            sysr = sysrfun(DrOpt);
        case 'DrRange'
            % Get steady state response of the error system
            Dr0 = DrInit('Ge0');

            % Define a range for the feedthrough
            nStep = 100;
            DrRange = linspace(-Dr0,Dr0, nStep); 

            % Run the actual function
            plotOverDrRange;
        %}
        case 'normOpt'
            Dr0 = DrInit(Opts.DrInit);
            cost = @(Dr) norm(syse0m - sysrDelta(Dr),Inf);
            [DrOpt, Hinf,tOpt] = normOpt(Dr0,cost);
            sysr = sysrfun(DrOpt); 
        case 'normOptCycle'
            % running optimization wr to each entry of D individually
            % the cost function takes into account the whole MIMO system
%             keyboard
            DrOpt = DrInit(Opts.DrInit); HinfVec = norm(syse0m,Inf); tOpt = 0;
            nCycles = 3; cycleCount = 0; stop = 0; %max number of cycles defined
            while cycleCount < nCycles && ~stop;
            cycleCount = cycleCount+1;
            for iOut = 1:sys.p
                for jIn = 1:sys.m
                    Dr0 = DrOpt(iOut,jIn);
                    cost = @(Dr) norm(syse0m - sysrDelta(Dr,iOut,jIn,DrOpt),Inf);
                    constr = @(Dr) stabilityConstraintCycle(Dr,iOut,jIn,DrOpt);
                    [DrOptCurr, Hinf,tOptCurr] = normOpt(Dr0,cost,constr);
                    tOpt = tOpt+tOptCurr;
                    DrOpt(iOut,jIn) = DrOptCurr;
                    HinfVec = [HinfVec, Hinf];
                    if length(HinfVec) > sys.p*sys.m +1 %cycled already once
                        crit = mean(abs(diff(HinfVec(end-2:end))./HinfVec(end-2:end-1)));
                        %mean improvement less that 1%
                        if crit < 1e-2; stop = 1; break, end
                    end
                end
                if stop, break, end %breaking out of nested for loop
            end
            end
            sysr = sysrfun(DrOpt); 
            if Opts.plot
                figure; plot(0:1:length(HinfVec)-1,HinfVec/HinfVec(1)); 
                ylabel('relative error decrease');
            end     
        case 'normOptCycleCombo'
            % running optimization wr to each entry of D individually
            % the cost function takes into account the whole MIMO system
            % after one cycle is run, we use the result to initialize a
            % multivariate optimization
            
            % 1) cycle optimization
%             Opts.solver = 'fmincon';
            DrOpt = DrInit(Opts.DrInit); 
%             if Opts.surrogateError
%                 HinfVec = norm(syse0m,Inf); 
%             else
%                 HinfVec = norm(sysm,Inf); 
%             end
            tOpt = 0;
                
            for iOut = 1:sys.p
                for jIn = 1:sys.m
                    Dr0 = DrOpt(iOut,jIn);
                    if Opts.surrogateError
                        cost = @(Dr) norm(syse0m-sysrDelta(Dr,iOut,jIn,DrOpt),Inf);
                    else
                        cost = @(Dr) norm(sysm-sysrfun(Dr,iOut,jIn,DrOpt),Inf);
                    end
                    constr = @(Dr) stabilityConstraintCycle(Dr,iOut,jIn,DrOpt);
                    [DrOptCurr, ~,tOptCurr] = normOpt(Dr0,cost,constr,'fmincon');
                    tOpt = tOpt+tOptCurr;
                    DrOpt(iOut,jIn) = DrOptCurr;
                end
            end
            
            % 2) multivariate optimization
            if Opts.surrogateError
                cost = @(Dr) norm(syse0m - sysrDelta(Dr),Inf);
            else
                cost = @(Dr) norm(sysm - sysrfun(Dr),Inf);
            end
            constr = @(Dr) stabilityConstraint(Dr);
            [DrOpt, Hinf,tOptCurr] = normOpt(DrOpt,cost,constr,Opts.solver);
            tOpt = tOpt + tOptCurr;
            
            sysr = sysrfun(DrOpt); 
%             sysr = ss(sysr0 + sysrDelta(DrOpt),'minimal');
%             figure; sigma(sysrfun(DrOpt),sysr0 + sysrDelta(DrOpt),'--r')
        case 'sweepDr'
            cost = @(Dr) norm(ss(syse0m-sysrDelta(Dr)),Inf);
            [DrOpt, tOpt, DrArray,costArray] = sweepDr(cost);
            assignin('caller','DrArray',DrArray);
            assignin('caller','costArray',costArray);
            sysr = sysrfun(DrOpt); 
        otherwise
            error('Specified Hinf optimization type not valid');
    end
    
    sysr = sss(sysr);

    %   See how the cost behaves around the chosen minimum?
    if Opts.plotCostOverDr
        nStep = 20; kRange = 5;
        if isscalar(DrOpt) %SISO
            DrRange = linspace(DrOpt*(1-kRange),DrOpt*(1+kRange), nStep); 
        elseif any (size(DrOpt) == 1) % SIMO or MISO
            DrRange = cell(size(DrOpt)); %get the right size
%             DrRange1 = linspace(DrOpt(1)*(1-kRange),DrOpt(1)*(1+kRange), nStep);
%             DrRange2 = linspace(DrOpt(2)*(1-kRange),DrOpt(2)*(1+kRange), nStep);
            DrRange1 = linspace(-6,6, nStep);
            DrRange2 = linspace(-8,8, nStep);
            
            [DrRange1, DrRange2] = meshgrid(DrRange1,DrRange2);
            DrRange{1} = DrRange1; DrRange{2} = DrRange2;
        end
        plotOverDrRange(DrOpt,Hinf);
    end
    %% Additional outputs
    if nargout > 1
        % Error indicators
%         if ~exist('Hinf','var')
%             Hinf = norm(sysm-sysr,inf);%optimized
%         end
%         HinfO = norm(sysm,inf); %original
%         HinfRel = Hinf/HinfO;
%         HinfRatio = Hinf/norm(sysm-sysr0,inf); %ratio to irka ROM

        % Real error
        HinfErr = norm(ss(sys-sysr),inf);%optimized
        Hinf0 = norm(ss(sys)-sysr0,inf); %before optimization
        HinfRatio = Hinf/Hinf0; %ratio to irka ROM
        
%         if Opts.debug, keyboard, end
       
        if nargout > 5
%             bound = HinfBound(sys,B_,C_);
            bound  =[];
        end
    end
    %% ======= Auxiliary ========
    function Dr0 = DrInit(type)
        %   Initialize Dr for optimization
        switch type
            case '0'
                Dr0 = zeros(sys.p,sys.m);
            case 0
                Dr0 = DrInit('0');    
                type = '0'; %overwrite for stability check at the end
            case 'Ge0'   
                Dr0 = freqresp(sysem,0)-freqresp(sysr0,0);      
            case 'Ge0half'
                Dr0 = DrInit('Ge0')/2;
            case 'matchGe0_old'        
                % this computation does not work well for MIMO systems
                % since is sweeping Dr0 only along a hyperline with
                % direction of the absolute value of Dr0
                G0 = freqresp(sysem,0); %the only costly part

                Dr0 = DrInit('Ge0'); %get an initial feedthrough
                deltaDr = 20*abs(Dr0); nStep = 100; 
                DrSet(:,:,1) = DrInit('0'); dSet(:,:,1) = Dr0; dMin = norm(Dr0); %initial error
                if Opts.plot, figure; sigma(sysem-sysr0,'b', 'LineWidth',2); end
                
                for k = 0:nStep
                    Dr = Dr0-deltaDr + k*(2*deltaDr)/nStep;
                    sysr = sysrfun(Dr);

%                     d = abs(abs(G0-freqresp(sysr,0)) - abs(Dr));
                      d = norm(G0 - (freqresp(sysr,0) + Dr));
                    if d < dMin && isstable(sysr)
                        dMin = d;
                        dSet(:,:,end+1) = d;
                        DrSet(:,:,end+1)= Dr;
                        if Opts.plot
                            syse = sysem-sysr; sigma(syse,'Color',rand(1,3)); 
                            drawnow
                        end
                    end
                end
                %   best feedthrough in term of minimizing the error between the
                %   response at 0 and Inf
                Dr0 = DrSet(:,:,end); 
                if Opts.plot 
                    syse = sysem-sysrfun(Dr0); sigma(syse,'Color',rand(1,3),...
                                                   'LineWidth',2); 
                    drawnow 
                end
            case 'matchGe0'      
                %   Run a sweep in an sys.m x sys.p dimensinal space
                %   The claim is that all computations are cheap since we
                %   compute G0 once.
                
                G0 = freqresp(sysem,0); %the only costly part
                sweepcost = @(Dr0) norm(G0 - (freqresp(sysrfun(Dr0),0) + Dr0));
                
                Dr0 = sweepDr(sweepcost);
            case 'maxGe'
                % finde the frequency w at which the maximum singular value
                % of the transfer function matrix is obtained. Conduct a
                % sweep to minimize the error at w
                
                %expensive computations
                [~, w] = norm(ss(sysem-sysr0),inf);
                Gew = freqresp(sysem,w);
                
                sweepcost = @(Dr0) norm(Gew - freqresp(sysrfun(Dr0),w));
                
                Dr0 = sweepDr(sweepcost);              
            otherwise
                error('Initialization option for Dr not valid');
        end
        % make sure the initialization yields a stable system 
        %(if not 0 and if not called by DrInit itself)
        caller = dbstack;
        if ~strcmp(caller(2).name,'HinfMor/DrInit')&& ~strcmp(type,'0') &&...
                ~isstable(sysrfun(Dr0))
            warning('Selected initialization for Dr0 would have yielded an unstable system. Changing it to 0');
            Dr0 = DrInit('0');
        end
    end
    function [DrOpt, tOpt, DrArray, costArray] = sweepDr(sweepcost)
        % This is used in parallel (parfor) 
        
        %   Determine the relevant frequency range for Dr based on the Hinf
        %   error after IRKA
        [~,wmax] = norm(ss(syse0m),inf); deltaDr = freqresp(syse0m,wmax);
        
        probSize = sys.m*sys.p; %dimension of the search space

        %CURSE of DIMENSIONALITY!!
        nPointsMax = Opts.sweepPoints;
        nStep = floor(nthroot(nPointsMax,probSize));
        if mod(nStep,2) == 0
          %number is even
          nStep = nStep-1;
        end 
        % take at least two points
        nStep = max([nStep,2]);

        % Create the grid
        x = {}; sampling = 'random';
        for iOut = 1:sys.p
            for jIn = 1:sys.m
                switch sampling
                    case 'grid'
                        x = [x, ...
                        {linspace(-deltaDr(iOut,jIn),deltaDr(iOut,jIn),nStep).'}];
                    case 'random'
                        x = [x, ...
                        {[0;randn(nStep,1)*3*deltaDr(iOut,jIn)]}];
                end
            end
        end
        Xl = cell(1,probSize); [Xl{:}] = ndgrid(x{:});   
        
        % Run the sweep
        nPoints = nStep^probSize;     
        DrArray = zeros(sys.p,sys.m,nPoints); 
        costArray = zeros(1,1,nPoints);  
        p = sys.p; m = sys.m;
        tic
        parfor k=1:nPoints  
            Dr0l = zeros(p,m);
            %   generate the current feedthrough
            for iEl = 1:probSize, Dr0l(iEl) = Xl{iEl}(k); end
            DrArray(:,:,k)= Dr0l;
            costArray(:,:,k) = sweepcost(Dr0l);
        end    
        tOpt = toc;
        [~, idxMin] = min(squeeze(costArray));
        DrOpt = DrArray(:,:,idxMin);    
        if Opts.plot, 
            figure; plot(squeeze(costArray)); hold on
            plot([1,length(costArray)],[1,1]*sweepcost(zeros(sys.p,sys.m)),'r--');
            title('Cost for Dr sweep'); ylabel('cost'); xlabel('index');
        end
    end
    function [DrOpt, Hinf,tOpt] = normOpt(Dr0,cost,constr,solver)       
            switch solver
                case 'fminsearch'
                tic, [DrOpt, Hinf] = fminsearch(cost,Dr0); tOpt = toc;               
                case 'fminunc'
                optOpts = optimoptions(@fminunc, 'algorithm','quasi-newton');
                tic, [DrOpt, Hinf] = fminunc(cost,Dr0,optOpts); tOpt = toc;               
                case 'fmincon'
                optOpts = optimoptions('fmincon','UseParallel',1,...
                                        'algorithm',Opts.algorithm,...
                                        'MaxFunEvals',5e2,...
                                        'Display',Opts.display);
                if ~exist('constr','var')
                    constr = @stabilityConstraint;
                end
                tic, [DrOpt, Hinf] = fmincon(cost,Dr0,[],[],[],...
                    [],[],[],constr,optOpts); tOpt = toc;
               
                case 'gs'                    
                % Restrict the search space to improve execution
                [lb,ub] = searchSpaceLimits(syse0m);    
                
                % Define optimization parameters
                optOpts = optimoptions('fmincon','UseParallel',1,...
                                        'algorithm',Opts.algorithm,...
                                        'MaxFunEvals',2e2);
                problem = createOptimProblem('fmincon',...
                            'objective',cost,'x0',Dr0,'options',optOpts,...
                            'nonlcon',@stabilityConstraint,...
                            'lb',lb,'ub',ub);
                gs = GlobalSearch('NumStageOnePoints',20,...%start points
                                  'NumTrialPoints',400,... %set of all potential start points
                                  'StartPointsToRun','bounds-ineqs',... %exclude certain points?
                                  'Display','iter',...
                                  'MaxTime',300); %stop after 5 min
                                 
                tic,  [DrOpt,Hinf] = run(gs,problem); tOpt = toc;                              
                case 'ms'
                optOpts = optimoptions('fminunc', 'algorithm','quasi-newton');
                problem = createOptimProblem('fminunc',....
                            'objective',cost, 'x0',Dr0,'options',optOpts);
                ms = MultiStart('UseParallel',true);
                tic, [DrOpt,Hinf,~,~,allmins] = run(ms,problem,25); tOpt = toc;
                
                %   Get the best stable result
                k = 1; nOpt = length(allmins);
                while ~isstable(sysrfun(DrOpt))
                    k = k+1;
                    DrOpt = allmins(k).X;
                    Hinf  = allmins(k).Fval;
                    if k >= nOpt
                        warning('none of the local optima from MS was stable')
                        DrOpt = allmins(k).X;
                        Hinf  = allmins(k).Fval;
                        break
                    end
                end                
                case 'ga'
                    options = gaoptimset('UseParallel',true,...
                                         'PopulationSize',25,...
                                         'Generations',20,...
                                         'InitialPopulation',Dr0(:).',...
                                         'Display','iter');
                    
                    hybridoptions = optimoptions('fmincon','UseParallel',1,...
                                    'Display','iter','algorithm','sqp');
                    options = gaoptimset(options,...
                    'HybridFcn',{@fmincon, hybridoptions});

                    gacost = @(x) cost(reshape(x,size(Dr0,1),size(Dr0,2)));

                    tic, [xOpt, Hinf] = ga(gacost,numel(Dr0),[],[],[],[],...
                            [],[],@stabilityConstraintGA,options); tOpt = toc;
                    DrOpt = reshape(xOpt,size(Dr0,1),size(Dr0,2));
            end
    end
    function sysr = sysrfun(Dr,iOut,jIn,DrMIMO)
        if nargin == 1
            %Dr is the full feedthrough
            DrMIMO = Dr;
        else
            DrMIMO(iOut,jIn) = Dr;
        end
        sysr = dss(full(sysr0.A)+Lt.'*DrMIMO*Rt, sysr0.B+Lt.'*DrMIMO, ...
                                     sysr0.C+DrMIMO*Rt, DrMIMO, full(sysr0.E));
    end
    function sysrDelta = sysrDelta(Dr,iOut,jIn,DrMIMO)
        if nargin == 1
            %Dr is the full feedthrough
            DrMIMO = Dr;
        else
            DrMIMO(iOut,jIn) = Dr;
        end
        sysrDelta1 = sysr0; sysrDelta1.B = Lt.'*DrMIMO; %sysrDelta1.D = Dr;
        sysrDelta2 = sysr0; sysrDelta2.C = DrMIMO*Rt; sysrDelta2.B = sysr0.B+Lt.'*DrMIMO;
        sysrDelta3 = sysr0; sysrDelta3.C = sysr0.C+DrMIMO*Rt; sysrDelta3.B = Lt.';
        sysrDelta4 = sysr0; sysrDelta4.C = -DrMIMO*Rt; sysrDelta4.B = Lt.'; sysrDelta4.D = eye(sys.p);

        sysrDelta = sysrDelta1 + sysrDelta2 + sysrDelta3*inv(sysrDelta4)*sysrDelta2 + Dr;    end
    function [minDr, minVal] = plotOverDrRange(varargin)
        if isnumeric(DrRange) %vector of Dr values
            %initializing
            normVec = zeros(1,length(DrRange));
            normO = norm(syse0m,inf); 
            normVec(1) = norm(syse0m-ss(sysr0),inf)/normO; 
            minVal = normVec(1); minDr = 0; DrRange = [0,DrRange];
            for iDr = 2:length(DrRange)
                sysrTest = sysrfun(DrRange(iDr));
                normVec(iDr) = norm(syse0m-ss(sysrTest),inf)/normO;
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
            currDr = 0; currVal = norm(syse0m-sysrfun(currDr),inf)/normO;
            h(4) = plot(currDr, currVal,'ok');          
            currDr = freqresp(syse0m-sysrfun(0),0); currVal = norm(syse0m-sysrfun(currDr),inf)/normO;
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
            normO = norm(syse0m,inf); 
            % in this case, we don't add 0 as a value to the grid
            minVal = Inf;
            for iDr = 1:size(Dr1Range,1)
                for jDr = 1:size(Dr1Range,2)
                    currDr(1) = Dr1Range(iDr,jDr);
                    currDr(2) = Dr2Range(iDr,jDr);
                    sysrTest = sysrfun(currDr);
                    normMat(iDr,jDr) = norm(syse0m-ss(sysrTest),inf)/normO;
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
            currDr = zeros(size(currDr)); currVal = norm(syse0m-sysrfun(currDr),inf)/normO;
            h(4) = plot3(currDr(1),currDr(2),currVal,'ok','MarkerFaceColor','k');            
            currDr = abs(freqresp(syse0m,0) - freqresp(sysr0,0));
            currVal = norm(syse0m-sysrfun(currDr),inf)/normO;
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
            
            % Add contour
            fh = figure;
            cValMin = min(min(normMat));
            cValMax = max(max(normMat));
            nMin = 20; nMax = 10;
            cVal1 = logspace(log10(cValMin),log10(cValMax/4),nMin);
            cVal2 = logspace(log10(cValMax/4),log10(cValMax),nMax);
            cVal = [cVal1(1:end-1),cVal2];
%             cVal = linspace(min(min(normMat)),max(max(normMat)),100);
            contour(Dr1Range,Dr2Range,normMat,cVal);
%             lh = contour(Dr1Range,Dr2Range,log10(normMat),50);
            colormap jet
            caxis([cValMin,cValMax])
            xlabel('$D_{r,1}$','Interpreter','latex'); 
            ylabel('$D_{r_2}$','Interpreter','latex');
            title('Countour plot of cost function for SIMO model','Interpreter','none');
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
    function [c,ceq]=stabilityConstraint(x)
        % define a nonlinear constraint to impose stability
        ceq = 1-isstable(sysrfun(x));
        c = [];
    end
    function [c,ceq]=stabilityConstraintGA(x)
        % define a nonlinear constraint to impose stability
        [c,ceq]=stabilityConstraint(reshape(x,size(Dr0,1),size(Dr0,2)));
    end
    function [c,ceq]= stabilityConstraintCycle(Dr,iOut,jIn,DrOpt)
        ceq = 1-isstable(sysrfun(Dr,iOut,jIn,DrOpt));
        c = [];
    end
    function [lb,ub] = searchSpaceLimits(syse)
        %syse: error system
        %[lb,ub] matrix values lower and upper bounds on Dr
        
        fac = 3; %scaling/robustness factor
        %initialize
        
        p = size(syse.c,1); m = size(syse.b,2);
        lb = zeros(p,m); ub = zeros(p,m);
        
        % define values
        for iO = 1:p
            parfor jI = 1:m
                val = norm(ss(syse(iO,jI)),inf);
                lb(iO,jI) = -fac*val; ub(iO,jI) = fac*val;
            end
        end
        
    end
    function syse0m = createSurrogate(syse0, type)
%         if Opts.debug, keyboard, end
        switch type  
            case 'original'
                syse0m = ss(syse0);
            case 'model' %is equivalent to Loewner
                [s0m,Rtm,Ltm] = getModelData(s0Traj,RtTraj,LtTraj);
                arnoldiOpts.makeOrth = 1;
                [~,V,W] = rk(syse0,s0m,s0m,Rtm,Ltm,arnoldiOpts);
                [V,~,~] = svd(V,0); [W,~,~] = svd(W,0);
                syse0m = dss(W'*syse0.A*V, W'*syse0.B, syse0.C*V,syse0.D,W'*syse0.E*V);               
            case 'loewner'
                %   Get the data
                ismemberTol = Opts.tol;
                [s0m,Rtm,Ltm] = getModelData(s0Traj,RtTraj,LtTraj,ismemberTol);
%                 while length(s0m)> size(syse0.A,1)
%                     ismemberTol = ismemberTol*10
%                     [s0m,Rtm,Ltm] = getModelData(s0Traj,RtTraj,LtTraj,ismemberTol);
%                 end
                
                rkOpts = struct('real',false,'orth',false);
%                 rkOpts = struct();
                
                [~,V,W] = rk(sss(syse0),s0m,s0m,Rtm,Ltm,rkOpts);
                                
                %   Create Loewner matrices
                L   = - W.'*syse0.E*V; %Loewner matrix
                sL  = - W.'*syse0.A*V; %shifted Loewner matrix

                %   Deflate
%                 s = svd(L - sL); figure; semilogy(s/s(1),'o-'); title('Normalized singular values of L - sL')
%                 s = svd([L, sL]); figure; semilogy(s/s(1)); title('Normalized singular values')
                s = svd([L, sL]); 
                if Opts.plot, 
                    figure; 
                    semilogy(s/s(1),'o-'); title('Normalized singular values of [L, sL]');
                    hold on; plot([1,length(s)],[Opts.surrTol,Opts.surrTol],'r--')
                    plotName = sprintf('%s_%s_n%i_LoewnerDecay',sys.Name,Opts.surrogate,n);
                    title(plotName,'Interpreter','none'); 
                    saveas(gcf,fullfile('..','res',[plotName,'.fig']));
                end

%                 r = find(s/s(1)<Opts.rankTol,1); if isempty(r), r = length(s); end 
                r = rank([L,sL],Opts.rankTol);
                if r == rank([L; sL], Opts.rankTol)
                    for iS = 1:length(s0m)
                        if rank(s0m(iS)*L-sL,Opts.rankTol) == r
                            break
                        end
                    end
                else
                    warning('Loewner conditions not satisfied');
                    iS = 1;
%                     syse0m = syse0; return
                end
                
                %   Compute the projection using a real, orthonormal basis
%                 [~,V,W] = rk(sss(syse0),s0m,s0m,Rtm,Ltm);
%                 L   = - W'*syse0.E*V; %Loewner matrix
%                 sL  = - W'*syse0.A*V; %shifted Loewner matrix

                [Ws, ~, Vs] = svd(s0m(iS)*L-sL,'econ');
                if Opts.debug
                    rs = input('Choose surrogate order: ');
                else
                    rs = find(s/s(1)< Opts.surrTol,1);
                    if isempty(rs), rs = r; end
                end
                V= V*Vs(:,1:rs); W= W*Ws(:,1:rs);
                
                %   Build the model function
                syse0m = dss(W.'*syse0.A*V, W.'*syse0.B, syse0.C*V,syse0.D,W.'*syse0.E*V);
            case 'vf'  
                
                happy = 'n';
                
                % Run Loewner (to determine reduced order)
%                 syse0mLoew = createSurrogate(syse0,'loewner');               
%                 nm = size(syse0mLoew.A,1);
                
                while ~strcmp(happy,'y')
                close all
                
                nm = input('Choose surrogate order: ');
                s0m = reshape(s0Traj,1,numel(s0Traj));
%                 [s0m] = getModelData(s0Traj,RtTraj,LtTraj,Opts.tol);
%                 s0m = reshape(s0Traj,1,numel(s0Traj));
%                 s0m = cplxpair(s0m); idx = find(imag(s0m)); s0m(idx(1:2:end)) = [];
                % remove real shifts
%                 s0m(imag(s0m)==0) = [];
                
                %avoid blowing-up for MIMO
%                 m = size(syse0.b,2);
%                 if m>1, nm = ceil(nm/m);end 
                %resize nm according to the data available
%                 nm = min([nm, ceil(2*length(s0m)/m)]);
%                 nm = min([nm, length(s0m)-2]);
                %                 
                % if mod(n,2) ~= 0, n = n-1; end   %make even
                Opts.forceReal = false;
                syse0m = vectorFitting(syse0,nm,s0m,Opts);
                
                if Opts.plot
                    plot(complex(eig(syse0m)),'xg');
                end
                
                % Compare to loewner
                if Opts.plot
                    figure('Name','Compare VF to Loewner');
                    bodemag(syse0,'b-',syse0m,'-.g'); 
                    legend(sprintf('original (n=%i)',size(syse0.A,1)),...
                           sprintf('vf (n=%i)',size(syse0m.A,1)))
                   drawnow
                end
                
                
                happy = input('Are you happy with the VF result? (y/n)  ','s');
                end
%                 if Opts.debug, keyboard; end
        end
        %                 isstable(sysm)
%         if Opts.plot
%             figure('Name','Original Vs surrogate models');
%             bodemag(syse0,'b-',syse0m,'--r'); 
%             legend(sprintf('original (n=%i)',size(syse0.A,1)),...
%                    sprintf('%s (n=%i)',type, size(syse0m.A,1)));  
%         end
%         if Opts.debug, keyboard, end
    end
    function [s0m,Rtm,Ltm] = getModelData(s0Traj,RtTraj,LtTraj,tol)
        % get interpolation data out of the trajectories
        kIrka = size(s0Traj,3); nRed = size(s0Traj,2);
        % initialize
        s0m = s0Traj(:,:,1); Rtm = RtTraj(:,:,1); Ltm = LtTraj(:,:,1);
        for iStep = 2:kIrka
            s0new = s0Traj(:,:,iStep);
            Rtnew = RtTraj(:,:,iStep); Ltnew = LtTraj(:,:,iStep);
            
            switch Opts.whatData
                case 'new'
                    idxS = ismemberf(s0new,s0m,'tol',tol); %available in the matlab central
                    idxR = ismemberf(Rtnew.',Rtm.','rows','tol',tol).';
                    idxL = ismemberf(Ltnew.',Ltm.','rows','tol',tol).';
                    idxNew = or(or(~idxS, ~idxR),~idxL);
                case 'all'
                    idxNew = 1:nRed; %take all
            end
            
            s0m = [s0m, s0new(idxNew)];
            Rtm = [Rtm, Rtnew(:,idxNew)]; Ltm = [Ltm, Ltnew(:,idxNew)];
        end
        
%         %   Do complexpair
%         s0mUnsrt = s0m; s0m = cplxpair(s0mUnsrt);
%         % get permutation indices, since cplxpair does not do it for you
%         [~,cplxSorting] = ismember(s0m,s0mUnsrt);
%         Rtm = Rtm(:,cplxSorting); Ltm = Ltm(:,cplxSorting);        
    end
    %% Trash
    function Dr0 = sweepDr0(sweepcost)
        %   Compute Dr0 from a sweep that tries to minimize "cost"
        
        if Opts.plot, figure; sigma(syse0m-sysr0,'b', 'LineWidth',2); end

        Dr0 = DrInit('0'); % Initialize with the data at the origin
        DrSet(:,:,1)= Dr0; dMin = sweepcost(Dr0); dSet(:,:,1) = dMin;

        deltaDr = 2*abs(DrInit('Ge0'));
        probSize = sys.m*sys.p; %dimension of the search space

        %CURSE of DIMENSIONALITY!!
        if probSize < 9 , nStep = 9; else nStep = 5; end

        % Create the grid
        x = {}; 
        for iOut = 1:sys.p
            for jIn = 1:sys.m
                x = [x, ...
                    {linspace(-deltaDr(iOut,jIn),deltaDr(iOut,jIn),nStep).'}];
            end
        end
        X = cell(1,probSize); [X{:}] = ndgrid(x{:});   
        
        % Run the sweep
        nPoints = nStep^probSize;                
        for k=1:nPoints    
            %   generate the current feedthrough
            for iEl = 1:probSize, Dr0(iEl) = X{iEl}(k); end

            sysr = sysrfun(Dr0);

            d = sweepcost(Dr0);
            if d < dMin && isstable(sysr)
                dMin = d;
%                 dSet(:,:,end+1) = d; %only for analysis
%                 DrSet(:,:,end+1)= Dr0;
                DrSet = Dr0;
                if Opts.plot
                    syse = syse0m-sysr; sigma(syse,'Color',rand(1,3)); 
                    drawnow
                end
            end
        end
        
        %   best feedthrough
        Dr0 = DrSet(:,:,end); 
        if Opts.plot 
            syse = syse0m-sysrfun(Dr0); sigma(syse,'Color',rand(1,3),...
                                           'LineWidth',2); drawnow 
        end          
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




