function [V,S,R,k] = spark(sys,s0,Opts)
% SPARK - Stability Preserving Adaptive Rational Krylov
% 
% Syntax:
%       [V,Sv,Rv,k] = SPARK(sys,s0)
%       [V,Sv,Rv,k] = SPARK(sys,s0,Opts)
%
% Description:
%       This function reduces a state-space, LTI model specified 
%       by the matrices A,B,C,E to a LTI model of order 2 using the trust
%       region optimization algorithm known as stability-preserving 
%       adaptive rational Krylov (SPARK).
%
%       Selecting the option 'model', the algorithm known as MESPARK is
%       computed, where a model function approximating the original model
%       is computed and reduced iteratively.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:       original full order model; 
%       -s0:        Initial shifts
%       *Optional Input Arguments:*
%       -Opts:      Structure containing computation options
%           -.spark.type:   chooses between standard SPARK, where the original 
%                           model is reduced directly, or MESPARK, where a 
%                           model function is created and updated after convergence.
%                           [{'model'} / 'standard']
%           -.spark.test:   specifies weather the user desires to get insight 
%                           in what is happening. This is realized by  
%                           plotting intermediate results during optimization.
%                           [{'0'} / '1']
%           -.spark.verbose: text output during the optimization
%                           [{'0'} / '1']
%           -.spark.mfe:    maximum functions evaluations
%                           [{'5e3'} / positive integer]
%           -.spark.mi:     maximum iterations in solver
%                           [{'150'} / positive integer]
%           -.spark.xTol:   step tolerance in solver
%                           [{'1e-10'} / positive float]
%           -.spark.fTol:   function value tolerance
%                           [{'1e-10'} / positive float]
%           -.spark.modelTol: convergence tolerance for model funciton
%                           [{'1e-5'} / positive float]
%           -.mespark.ritz: use eigenvalues of model function to initialize
%                           the shifts 
%                           [{'1'} / '0']
%           -.mespark.pertIter: number of iterations after which a
%                           pertubation of the shifts starts to avoid
%                           stagnation of the model function
%                           [{'5'} / positive integer]
%           -.mespark.maxIter: maximum number of model function updates
%                           [{'20'} / positive integer]
%
% Output Arguments:      
%       -V,S,R:     Krylov subspace [V or W], Sylvester matrices [Sv,Rv] or
%                   [Sw^T,Lw]
%       -k:         Number of iterations of MESPARK
%
% 
% See Also: 
%       cure, porkV, porkW, rk, getSylvester
%
% References:
%       * *[1] Panzer (2014)*, Model Order Reduction by Krylov Subspace Methods
%              with Global Error Bounds and Automatic Choice of Parameters
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko K.F. Panzer, Alessandro Castagnotto 
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  13 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Parse input and load default parameters
% default values
Def.spark.type = 'model'; %SPARK type, 'model' or 'standard'
Def.spark.test = 0; %execute analysis code
Def.spark.verbose = 0; %show text?
Def.spark.mfe = 5e3;
Def.spark.mi = 150; %5e3
Def.spark.xTol = 1e-10;
Def.spark.fTol = 1e-10;
Def.spark.modelTol = 1e-5;
    Def.mespark.ritz = 1;
    Def.mespark.pertIter = 5; % # iteration at which perturbation begins
    Def.mespark.maxIter = 20; %maximum number of model function updates

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
else
    if ~isfield(Opts,'spark') || isempty(Opts.spark)
        Opts.spark = Def.spark;
    else
        Opts.spark = parseOpts(Opts.spark,Def.spark);
    end  
    if ~isfield(Opts,'mespark') || isempty(Opts.mespark)
        Opts.mespark = Def.mespark;
    else
        Opts.mespark = parseOpts(Opts.mespark,Def.mespark);
    end
end
    
%   Check input size 
if size(sys.B,2)>1 || size(sys.C,1)>1, error('System must be SISO.'), end

%let all functions have access to the figure window created in ESPARK
global fh  

%% ---------------------------- CODE -------------------------------
warning('off','MATLAB:nearlySingularMatrix')
    % 
    p0 = [(s0(1)+s0(2))/2, s0(1)*s0(2)];   
    t = tic; precond = eye(2); 
    
    Opts.fmincon=optimset('TolFun',Opts.spark.fTol,'TolX',Opts.spark.xTol, ...
        'Display','none', 'Algorithm','trust-region-reflective', ...
        'GradObj','on','Hessian','on','MaxFunEvals',Opts.spark.mfe,'MaxIter',...
        Opts.spark.mi);
    if Opts.spark.test
        Opts.fmincon = optimset(Opts.fmincon,...
        'OutputFcn', @OuputFcn,...
        'PlotFcns',{@optimplotx, @optimplotfval, @optimplotfirstorderopt});
    end
    
    switch Opts.spark.type
        case 'standard'
            %definition of "mock model function" for consistency in CostFunction
            Am=sys.A; Bm=sys.B; Cm=sys.C; Em=sys.E;
            
            p_opt = ESPARK(p0);
            k = [];
        case 'model'
            [p_opt,k] = mespark(p0,s0);
    end
    
    % supply output variables
    v1  = Q1*(U1\(L1\(P1*sys.B))); v12= Q2*(U2\(L2\(P2*sys.B))); 
    v2 = Q2*(U2\(L2\(P2*(sys.E*v1))));
    V   = full(real([v1/2 + (v12/2+p_opt(1)*v2), v2*sqrt(p_opt(2))]));
    S = [2*p_opt(1), sqrt(p_opt(2)); -sqrt(p_opt(2)), 0]; R = [1 0];
    if Opts.spark.verbose
    disp(['spark required ca. ' num2str(2*(k+1)) ' LUs ', ...
        ' and converged in ' num2str(toc(t),'%.1f') 'sec.'])
    end
  
    warning('on','MATLAB:nearlySingularMatrix')
    
    %% ------------------ AUXILIARY FUNCTIONS --------------------------
    % a) PRIMARY
    function p_opt = ESPARK(p0)
        
        if Opts.spark.test, fh = plotCost(p0); end
        [~,~,H] = CostFunction(p0); precond = diag(1./abs(diag(H).^0.25));
        
        % run trust region algorithm to find minimum of model function
        [p_opt,~,flag,output] = fmincon(@CostFunction,p0/precond,[],[],[],[],[0;0],[inf;inf],[],Opts.fmincon);
        if ~flag==1 % premature abortion
            warning('%s\n',output.message);
        end
        p_opt = p_opt*precond; precond = eye(2);

        if Opts.spark.test
            figure(fh);
            plot(p_opt(1),p_opt(2),'p','LineWidth',2,'MarkerSize',10,...
                'MarkerFaceColor',TUM_Gruen,'MarkerEdgeColor','k');
            pause(2),close(fh);
        end

        % convert parameter to shifts and perform LU decompositions
        s_opt=p_opt(1)+[1,-1]*sqrt(p_opt(1)^2-p_opt(2)); computeLU(s_opt);
    end
    function [p_opt,k] = mespark(p0,s0)
        k = 0;
        % compute initial model function and cost function at p0
        computeLU(s0);  V = newColV([],3);  W = newColW([],3);
        Am=W'*sys.A*V; Bm=W'*sys.B; Cm=sys.C*V; Em=W'*sys.E*V;
        J_old = CostFunction(p0);
        
        if Opts.mespark.ritz
            % Model-function-based initialization
            if Opts.spark.verbose,fprintf('User initialization: p0 =[%e,%e]',p0(1),p0(2));end
            p0 = ritz_initial;  
        end
        
        count = 0; %counter for perturbation if not improving
        if Opts.spark.verbose, fprintf('Starting MESPARK...\n'),end
        while(1)
            k = k + 1;
            if Opts.spark.verbose,fprintf('\tIteration %i: q=%i\n',k,size(V,2)),end
            
            p_opt = ESPARK(p0);
    
            % update model function by two-sided (Hermite) projection
            V = newColV(V, 2);  W = newColW(W, 2);  
            Am=W'*sys.A*V; Bm=W'*sys.B; Cm=sys.C*V; Em=W'*sys.E*V;
            % evaluate cost functional at new parameter point
            J = CostFunction(p_opt);

            if Opts.spark.verbose
            disp(['  relative change:      ' num2str(norm((p0-p_opt)./p0), '%1.2e')]);
            disp(['  relative improvement: ' num2str((J-J_old)/J, '%1.2e')]);
            disp(['  absolute J = ' num2str(J, '%1.12e')]);
            end
 
            % decide how to proceed
            if abs((J-J_old)/J) < Opts.spark.modelTol || ...
                    norm((p0-p_opt)./p0) < Opts.spark.modelTol %|| size(Am,1)>=20
                if Opts.spark.verbose,fprintf('Tolerance reached! Quitting MESPARK...\n'),end
                break;                      % convergence in J or in p  => stop
            elseif J<J_old
                 J_old = J;  p0 = p_opt;	% improvement: continue with p_opt
                 if Opts.spark.verbose
                     fprintf('\t\t updating the model function...\n')
                     fprintf('\t\t restarting MESPARK where it converged...\n')
                 end
                 count = 0; %reset stagnation counter
            else %no improvement
                if Opts.spark.verbose,fprintf('\t\t no improvement!\n'),end
                count = count+1; %add one to the stagnation counter
                
                % maximum iterations reached
                if k >= Opts.mespark.maxIter 
                    if Opts.spark.verbose
                        warning('Maximum number of iterations in MESPARK reached! Aborting...')
                    end
                    break;
                end
                
                % going on with MESPARK
                if Opts.spark.verbose,fprintf('\t\t updating the model function...\n'),end
                if count <Opts.mespark.pertIter
                    if Opts.spark.verbose,fprintf('\t\t restarting MESPARK where it began...\n'),end
                else
                    if Opts.spark.verbose,warning('long-term stagnation: perturbing p0...'),end
                    p0 = perturb(p0,count);
                end
            end
        end
    end
    % b) SECONDARY
    function [J, g, H] = CostFunction(p)
        % H2 cost functional, gradient and Hessian
        p = p*precond;  a = p(1); b = p(2);  s1 = a+sqrt(a^2-b); s2 = a-sqrt(a^2-b);
        r1 = (Am-s1*Em)\Bm;  r2 = (Am-s1*Em)\(Em*r1); r3 = (Am-s1*Em)\(Em*r2); 
        l1 = Cm/(Am-s2*Em);  l2 = l1*Em/(Am-s2*Em);   l3 = l2*Em/(Am-s2*Em);
        [J, g, H] = CostFunctionH2(Am, Bm, Cm, Em, p, [r1,r2,r3], [l1;l2;l3]);
        g = g * precond;  H = precond * H * precond;
    end
    function [J, g, H] = CostFunctionH2(A, B, C, E, p, r, l)
    % Enhanced spark Cost Functional
    %   Input:  A,B,C,E: HFM matrices; 
    %           p:       parameter vector [a,b]; 
    %           r,l:     left an right rational Krylov sequence;
    %   Output: cost functional J; gradient g; Hessian H
    % $\MatlabCopyright$

        a = p(1); b = p(2);  l1 = l(1,:); r1 = r(:,1);
        Pr = [4*a, 4*a^2; 4*a^2 4*a*(a^2+b)];  Cr = [0.5*(C*r1 + l1*B), l1*E*r1];
        J = real(-Cr*Pr*Cr');
        if nargout==1, return, end
        l2 = l(2,:); r2 = r(:,2);  l3 = l(3,:); r3 = r(:,3); 

        dPrda =     [4, 8*a; 8*a, 12*a^2 + 4*b];    dPrdb =     [0, 0; 0,  4*a];
        ddPrdada =  [0, 8; 8, 24*a];                ddPrdadb =  [0, 0; 0, 4];
        dcrda =     [0.5*(C*r2 + l2*B) + a*(l2*E*r1 + l1*E*r2), 2*l2*A*r2];
        dcrdb =     [-0.5*(l2*E*r1 + l1*E*r2), -l2*E*r2];
        ddcrdada =  [C*r3 + l3*B + 4*a*l1*E*r3 + 4*a*l3*E*r1 + 2*a*l2*E*r2 + 2*l2*A*r2 ...
                          + 4*a^2*l2*E*r3 + 4*a^2*l3*E*r2, ...
                     4*l3*A*r2 + 4*l2*A*r3 + 8*a*l3*A*r3];
        ddcrdadb =  [-l3*E*r1 - l1*E*r3 - l2*E*r2 - 2*a*l2*E*r3 - 2*a*l3*E*r2, -4*l3*A*r3];
        ddcrdbdb =  [ l3*E*r2 + l2*E*r3,  2*l3*E*r3];

        g = real([-Cr*dPrda*Cr', -Cr*dPrdb*Cr'] - 2*Cr*Pr*[dcrda; dcrdb]');
        H = [-2*ddcrdada*Pr*Cr'-4*dcrda*dPrda*Cr'-2*dcrda*Pr*dcrda'-Cr*ddPrdada*Cr', ...
             -2*ddcrdadb*Pr*Cr'-2*dcrdb*dPrda*Cr'-2*dcrda*dPrdb*Cr'-2*dcrda*Pr*dcrdb'-...
             Cr*ddPrdadb*Cr'; 0,-2*ddcrdbdb*Pr*Cr'-4*dcrdb*dPrdb*Cr'-2*dcrdb*Pr*dcrdb'];
        H(2,1)=H(1,2); H=real(H);
    end
    function [X,Y,Z] = GramSchmidt(X,Y,Z,cols)
        % Gram-Schmidt orthonormalization
        %   Input:  X,[Y,[Z]]:  matrices in Sylvester eq.: V,Sv,Rv or W.',Sw.',Lw.'
        %           cols:       2-dim. vector: number of first and last column to be treated
        %   Output: X,[Y,[Z]]:  solution of Sylvester eq. with X.'*X = I
        % $\MatlabCopyright$

        if nargin<4, cols=[1 size(X,2)]; end
        for k=cols(1):cols(2)
            for j=1:(k-1)                       % orthogonalization
                T = eye(size(X,2)); T(j,k)=-X(:,k)'*X(:,j);
                X = X*T;
                if nargout>=2, Y=T\Y*T; end
                if nargout>=3, Z=Z*T; end
            end
            h = norm(X(:,k));  X(:,k)=X(:,k)/h; % normalization
            if nargout>=2, Y(:,k) = Y(:,k)/h; Y(k,:) = Y(k,:)*h; end
            if nargout==3, Z(:,k) = Z(:,k)/h; end
        end
    end
    function computeLU(s0)
        % compute new LU decompositions
        if real(s0(1))==real(s0(2))  % complex conjugated or double shift
            [L1,U1,P1,Q1] = lu(sparse(sys.A-s0(1)*sys.E));  L2=conj(L1);U2=conj(U1);P2=P1;Q2=Q1;
        else                         % two real shifts
            [L1,U1,P1,Q1] = lu(sparse(sys.A-s0(1)*sys.E));  
            [L2,U2,P2,Q2] = lu(sparse(sys.A-s0(2)*sys.E));
        end
    end
    function V = newColV(V, k)
        % add columns to input Krylov subspace
        for i=(size(V,2)+1):2:(size(V,2)+2*k)
            if i==1, x=sys.B; else x=sys.E*V(:,i-1); end
            r1  = Q1*(U1\(L1\(P1*x)));   tmp = Q2*(U2\(L2\(P2*x)));
            v1 = real(0.5*r1 + 0.5*tmp); v2  = real(Q2*(U2\(L2\(P2*(sys.E*r1)))));
            V = GramSchmidt([V,v1,v2],[],[],[i,i+1]);
        end
    end
    function W = newColW(W, k)
        % add columns to output Krylov subspace
        for i=(size(W,2)+1):2:(size(W,2)+2*k)
            if i==1, x=sys.C; else x=W(:,i-1)'*sys.E; end
            l1  = x*Q1/U1/L1*P1;          tmp = x*Q2/U2/L2*P2;
            w1 = real(0.5*l1 + 0.5*tmp);  w2  = real(l1*sys.E*Q2/U2/L2*P2);
            W = GramSchmidt([W,w1',w2'],[],[],[i,i+1]);
        end
    end
    function p0 = ritz_initial
        l = eig(Am,Em); %compute Ritz values
        % remove unstable eigenvalues
        l = l(real(l)<=0);
        % make sure they are sorted in ascending magnitude of real part
        [~,idx] = sort(abs(real(l)),'ascend');
        l = l(idx)';
        % make sure there are at least two valid eigenvalues
        if length(l)>1
            l_ritz = l(1:2);
            % make sure they are complex conjugate
            if sum(imag(l_ritz))
                l_ritz(2) = conj(l_ritz(1));
            end
        else
            l_ritz = -ones(2,1);
        end
        
        if Opts.spark.test
            bla=figure('Name','Ritz Values for Initialization');plot(real(l),imag(l),'b*');hold on
            plot(real(l_ritz),imag(l_ritz),'ro');
            legend('Ritz Values','initialization')
        end
        % redefine the optimization starting values according to the
        % Ritz values.
        % Take the mirror images to get positive real parts
        s0ritz = l_ritz-2*real(l_ritz);
        p0 = s2p(s0ritz); 
        if Opts.spark.verbose,fprintf('Initialization according to Ritz values: p0=[%e,%e]',p0(1),p0(2));end
        if Opts.spark.test,pause,close(bla);end
    end
    function p0 = perturb(p0,count)
        % function used to perturb p0 in case of long-term stagnation
        %
        % Generate a normally distributed random variable with mean in p0
        % and an iteration-step-dependent standard deviation
        sd = p0*((count+1-Opts.mespark.pertIter)/Opts.mespark.pertIter);
        p0 = p0 + randn(size(p0)).*sd;
        
        % replace negative values by 0
        p0(p0<=0) = Opts.zeroThres;
    end
    function varargout = s2p(varargin)
        % s2p: Shifts to optimization parameters for spark
        % ------------------------------------------------------------------
        % USAGE:  This function computes the two paramters a,b that are used within
        % the optimization in spark. 
        %
        % p = s2p(s)
        % [a,b] = s2p(s1,s2)
        % 
        % Computations:
        % a = (s1+s2)/2
        % b = s1*s2
        %
        % s1 = a + sqrt(a^2-b)
        % s2 = a - sqrt(a^2-b)
        %
        % See also CURE, SPARK.
        %
        % ------------------------------------------------------------------
        % REFERENCES:
        % [1] Panzer (2014), Model Order Reduction by Krylov Subspace Methods
        %     with Global Error Bounds and Automatic Choice of Parameters
        % ------------------------------------------------------------------
        % This file is part of MORLab, a Sparse State Space, Model Order
        % Reduction and System Analysis Toolbox developed at the Institute 
        % of Automatic Control, Technische Universitaet Muenchen.
        % For updates and further information please visit www.rt.mw.tum.de
        % For any suggestions, submission and/or bug reports, mail us at
        %                      -> MORLab@tum.de <-
        % ------------------------------------------------------------------
        % Authors:      Alessandro Castagnotto
        % Last Change:  27 April 2015
        % ------------------------------------------------------------------

        % parse input
        if nargin==1
            s1 = varargin{1}(1);
            s2 = varargin{1}(2);
        else
            s1 = varargin{1};
            s2 = varargin{2};
        end

        % compute
        a = (s1+s2)/2;
        b = s1*s2;

        % generate output
        if nargout<=1
            varargout{1} = [a,b];
        else
            varargout{1} = a;
            varargout{2} = b;
        end
    end
    % c) ANALYSIS
    function stop = OuputFcn(x,optimValues,state)
    
    stop = false;
    if ~strcmp(state,'init')
        if optimValues.firstorderopt < 1e-30 && optimValues.stepsize < 1e-30
            stop = true;
        end 
   
        figure(fh);
        x = x*precond;
        
        %   plot new parameter pair
        plot(x(1),x(2),'*k');

        %   plot the gradient
        g = optimValues.gradient;
        g = (g/norm(g,2)); %normalize, we care only about direction
        x_next = x - g';
        plot([x(1),x_next(1)],[x(2),x_next(2)],'-k');
          
    end
    end
    function fh = plotCost(p0)
        
    fh = figure('Name','Cost function for optimization');
    
    npoints = 25;
    nlines = 100;

    %     [DeltaA, DeltaB] = initialization_region(p0); %TBD
    Delta = [10,10]; %[DeltaA, DeltaB]
    pLog = round(log10(p0));
    
    x = logspace(pLog(1)-Delta(1)/2,pLog(1)+Delta(1)/2,npoints);
    y = logspace(pLog(2)-Delta(2)/2,pLog(2)+Delta(2)/2,npoints);
    
    [X,Y] = meshgrid(x,y);
    Z = zeros(size(X));
    for ii= 1:length(x)
        for jj = 1:length(x)
            Z(ii,jj) = CostFunction([X(ii,jj),Y(ii,jj)]);
        end
    end
    contour(X,Y,Z,nlines);colorbar;
    set(gca,'XScale','log','YScale','log');
    hold on
    plot(p0(1),p0(2),'p','LineWidth',2,'MarkerSize',10,...
    'MarkerFaceColor',TUM_Rot,'MarkerEdgeColor','k');
    
    % plot separatrix between complex and real shifts
    plot(x,x.^2,'--k');
    xlim([x(1),x(end)]); ylim([y(1),y(end)]);
    end
    function y = TUM_Gruen()
        y = [162 173 0]/255;
    end
    function y = TUM_Rot()
        y = [196 7 27]/255;
    end
end