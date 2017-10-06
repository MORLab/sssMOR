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
    
    if isempty(s0_out)
        if isempty(Rt),     Rt = 1;   end
        [s0_inp,Rt] = newParaInp(sys,sysr,V,s0_inp,Rt,Opts);
    elseif isempty(s0_inp)
        if isempty(Lt),     Lt = 1;   end
        [s0_out,Lt] = newParaOut(sys,sysr,W,s0_out,Lt,Opts);
    elseif s0_inp == s0_out
        if isempty(Rt),     Rt = eye(size(s0_inp,2));   end
        [s0_inp,Rt] = newParaInp(sys,sysr,V,s0_inp,Rt,Opts);
        [s0_out,Lt] = newParaOut(sys,sysr,W,s0_inp,Lt,Opts);
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
 
 % compute Ritz-Values of reduced system 
 ritzVal = eig(sysr.A);
 
 % build spectral Set
 if ~isreal(s0_inp)
     specSet = sort([s0_inp'; -ritzVal]);               % complex sectrum -> Mayer-Luenberger conditions
     chull = convhull(real(specSet),imag(specSet));    % bulid convex hull of spectral set
     specSet = specSet(chull);
 else
     specSet = sort([s0_inp'; -ritzVal]);  
 end
 
 % get new shift/tangential direction, solve max-problem
 resNorm_last = 0;         % Initialize and set variables
 
 % qr-decompositons, residual matrices
 if isempty(Opts.Bbot)
     Er_inv_Br = solveLse(sysr.E,sysr.B);
     Opts.Bbot = sys.B-(sys.E*V)*Er_inv_Br;
 end
 [~,R_Bbot] = qr(Opts.Bbot);
 [~,RV] = qr(V);
 [~,SRV,~] = svd(RV);
 %snew = specSet';
 % solve max-problem for new shift
 for ii = 1:1:size(specSet,1)
     Y = solveLse((sysr.A-specSet(ii,1)+sysr.E),sysr.B);
     resNorm = norm(R_Bbot)*(norm(Rt*Y+SRV(1,1))-1);
     %resNorm = norm(R_Bbot*(Rt*Rinv*Y-eye(size(sys.B,2))));
     if resNorm > resNorm_last && specSet(ii,1) ~= s0_inp(end)
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
     s0_inp = [s0_inp cplxpair([snew conj(snew)])];
    % Rt     = [Rt cplxpair([rt_new conj(rt_new)])];
 else
     s0_inp = [s0_inp snew];
     %Rt     = [Rt rt_new];
 end
 end

function [s0_out,Lt] = adaptiveShiftOut(sys,sysr,W,s0_out,Lt,Opts)

end



















