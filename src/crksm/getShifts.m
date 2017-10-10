function [s0_inp,s0_out,Rt,Lt] = getShifts(sys,sysr,s0_inp,s0_out,Rt,Lt,V,W,Opts)

%% Still  to come
% - input output der funktionen machen
% - alles für c-sided machen
% - multiple directions mit option machen

%% Create Def-struct containing default values/options
Def.purpose         = 'lyapunov';       % choose purpose: ['lyapunov' / 'MOR'] 
Def.strategy        = 'eigs';           % strategy for shift generation: ['ADI' / 'const' / 'ROM' / 'eigs'] 
Def.shiftNotation   = 'row';            % select output of shifts as row- or column-vector
Def.Bbot            = [];               % recycle Bbot if it is already available
Def.Er_inv_Ar       = [];               % recycle Er^-1*Ar if it is already available

%% Parsing of Inputs

% create the options structure and check tangential directions
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end
clear Def

%% Compute New Shift/Shifts

% get new shifts: 'ADI' / 'const' / 'ROM' / 'eigs' 
if strcmp(Opts.strategy,'ADI') && strcmp(Opts.purpose,'MOR')
    error('This method is only implemented for solving a Lyapunov equations!');
elseif strcmp(Opts.strategy,'adaptive') 
    % choose case
    if isempty(s0_out)
        if isempty(Rt),     Rt = 1;   end
        [s0_inp,Rt] = newParaInp(sys,sysr,V,s0_inp,Rt,Opts);
    elseif isempty(s0_inp)
        if isempty(Lt),     Lt = 1;   end
        [s0_out,Lt] = newParaOut(sys,sysr,W,s0_out,Lt,Opts);
    elseif s0_inp == s0_out
        if isempty(Rt),     Rt = eye(size(s0_inp,2));   end
        [s0_inp,Rt] = newParaInp(sys,sysr,V,s0_inp,Rt,Opts);
        s0_out = s0_inp;
        %******************************* hier kommt noch was fuer die tangential direction Lt
    else
        if isempty(Rt),     Rt = eye(size(s0_inp,2));   end
        if isempty(Lt),     Lt = eye(size(s0_out,2));   end
        [s0_inp,Rt] = newParaInp(sys,sysr,V,s0_inp,Rt,Opts);
        [s0_out,Lt] = newParaOut(sys,sysr,W,s0_out,Lt,Opts);
    end
else 
    snew =  My_initializeShifts(sys,Opts);
    snew = reshape(snew,[1,20]);
    s0_inp = [s0_inp snew];
end

% change shift vector to column notation
if strcmp(Opts.shiftNotation,'column')  && ~iscolumn(s0_inp)
    s0_inp = reshape(s0_inp,[size(s0_inp,2), 1]);
    s0_out = reshape(s0_out,[size(s0_out,2), 1]);
end

end % end of getShifts

%% ***************************** AUXILIARY ********************************

function [s0_inp,Rt] = newParaInp(sys,sysr,V,s0_inp,Rt,Opts)
    ritzVal = eig(sysr);  % compute Ritz-Values of reduced system 
    resNorm_last = 0;     % Initialize and set variables

    % build spectral Set
    if ~isreal(s0_inp)
        specSet = sort([s0_inp'; -ritzVal]);               % complex sectrum -> Mayer-Luenberger conditions
        chull = convhull(real(specSet),imag(specSet));     % bulid convex hull of spectral set
        specSet = specSet(chull);
    else
        specSet = sort([s0_inp'; -ritzVal]);  
    end

    % qr-decompositons, residual matrices
    if isempty(Opts.Bbot)
        Er_inv_Br = solveLse(sysr.E,sysr.B);
        Opts.Bbot = sys.B-(sys.E*V)*Er_inv_Br;
        Opts.Er_inv_Ar = solveLse(sysr.E,sysr.A,Opts);   
    end
    res0 = norm(Opts.Bbot);

    % solve max-problem for new shift
    for ii = 1:1:size(specSet,1)
        Y = solveLse((sysr.A-specSet(ii,1)*sysr.E)*sysr.E,sysr.B);
        resNorm = norm(((sys.A*V-sys.E*V*Opts.Er_inv_Ar)*Y-Opts.Bbot))/res0;
        if resNorm > resNorm_last && ~any(ismember(s0_inp,specSet(ii,1)))
            resNorm_last = resNorm;
            snew = specSet(ii,1);
            Ypic = Y;
        end
    end
    
    % compute (multiple) tangential directions
    if size(Rt,1) == size(sys.B,2)
        Rinv = solveLse(V,Vhat);
        res = R_Bbot*(Rt*Rinv*Ypic-eye(size(sys.B,2)));
        [~,S,rSingVec] = svd(res); 
        for ii = 1:1:size(S,2)
            if Opts.multDir == false
                resNorm(ii,1) = norm((R_Bbot*(rRt*Rinv*Y-eye(size(sys.B,2))))*rSingVec(:,ii));
                [~,index] = max(resNorm);
                rt_new = rSingVec(:,index);
            else              
                if S(ii,ii) > 0.1*S(1,1) && ii ~= index
                    rt_new = [rt_new rSingVec(:,ii)];
                end    
            end
        end
    else
        rt_new = [];
    end

    % enlarge shift/direction vector
    if ~isreal(snew)
        s0_inp = [s0_inp cplxpair([snew conj(snew)])];
        Rt     = [Rt cplxpair([rt_new conj(rt_new)])];
    else
        s0_inp = [s0_inp snew];
        Rt     = [Rt rt_new];
    end
end
 
function [s0_out,Lt] = newParaOut(sys,sysr,W,s0_out,Lt,Opts)
    ritzVal = eig(sysr);  % compute Ritz-Values of reduced system 
    resNorm_last = 0;       % Initialize and set variables

    % build spectral Set
    if ~isreal(s0_out)
        specSet = sort([s0_out'; -ritzVal]);               % complex sectrum -> Mayer-Luenberger conditions
        chull = convhull(real(specSet),imag(specSet));     % bulid convex hull of spectral set
        specSet = specSet(chull);
    else
        specSet = sort([s0_out'; -ritzVal]);  
    end

    % qr-decompositons, residual matrices
    if isempty(Opts.Cbot) || ~exist('Opts.Cbot','var')
        ErT_inv_CrT = solveLse(sysr.E',sysr.C');
        Opts.Cbot = sys.C'-(sys.E'*W)*ErT_inv_CrT;
    end
    [~,R_Cbot] = qr(Opts.Cbot);
    [~,RW] = qr(W);
    [~,SRW,~] = svd(RW);

    % solve max-problem for new shift
    for ii = 1:1:size(specSet,1)
        Y = solveLse((sysr.A-specSet(ii,1)+sysr.E)',sysr.C');
        resNorm = norm(R_Cbot)*(norm(Lt*Y+SRW(1,1))-1);
        %resNorm = norm(R_Bbot*(Rt*Rinv*Y-eye(size(sys.B,2))));
        if resNorm > resNorm_last && all(ismember(s0_out,specSet(ii,1)))
            resNorm_last = resNorm;
            snew = specSet(ii,1);
            Ypic = Y;
        end
    end

    % compute (multiple) tangential directions
    %  res = R_Bbot*(rRt*Rinv*Ypic-eye(size(sys.B,2)));
    %  [~,S,rSingVec] = svd(res); 
    %  for ii = 1:1:size(rSingVec,2)
    %      resNorm(ii,1) = norm((R_Bbot*(rRt*Rinv*Y-eye(size(sys.B,2))))*rSingVec(:,ii));
    %  end
    %  [~,index] = max(resNorm);
    %  rt_new = rSingVec(:,index);
    %  
    %  if Opts.multDir == true
    %      % hier brauche ich das S
    %      % ******************************************************
    %  end

    % enlarge shift/direction vector
    if ~isreal(snew)
        s0_out = [s0_out cplxpair([snew conj(snew)])];
        % Lt     = [Lt cplxpair([lt_new conj(lt_new)])];
    else
        s0_out = [s0_out snew];
        %Lt     = [Lt lt_new];
    end
end



















