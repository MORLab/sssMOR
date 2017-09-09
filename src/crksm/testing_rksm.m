clear
clc
clearvars -global

%% testdaten My_rksm

% load any benchmark system
sys = loadSss('fom');
% other models: rail_1357, rail_5177, iss, gyro, building, CDplayer, beam,
%               eady, heat-cont, LF10, random, bips98_606,
%               SpiralInductorPeec, fom

% declare system matices
A = sys.A;
B = sys.B;
C = sys.C;
E = sys.E;

m = size(B,2);
p = size(C,1);



%% testing function for several cases
% initialize method
method = 1; % 1 for standard rksm cyclic, reusing shifts, onesided, hermite
            % 2 for standard rksm cyclic, reusing shifts, onesided, hermite, tangential directions (MIMO)
            % 3 for standard rksm cyclic, reusing shifts, twosided, hermite
            % 4 for standard rksm cyclic, reusing shifts, twosided, not hermite
            % 5 for standard rksm cyclic, reusing shifts, twosided, not hermite, tangential
            % 6 for eksm method
            % 7 for standard rksm cyclic, adaptive shifts, onesided, hermite, only for SISO
            % 8 for standard rksm cyclic, adaptive shifts, twosided, hermite, only for SISO
            
% choose computation of residual norm and other options
Opts.residual  = 'residual_lyap';  
%Opts.residual  = 'norm_chol';
Opts.maxiter_rksm = 60;

Opts.maxiter = Opts.maxiter_rksm;
Opts.rctol = 1e-7; 

            
switch method
    case 1
        % set Opts for rksm 
        Opts.method = 'rksm';
        Opts.shifts = 'cyclic';
        %Opts.rksmnorm = 'fro';
        Opts.reduction = 'onesided';
        
        % compute shifts using irka
        shifts = zeros(1,10);
        Rt = ones(m,size(shifts,10));
        %Rt = rand(m,size(shifts,10));
        Lt = ones(p,size(shifts,10));
        %[~, ~, ~, s0, Rt, Lt] = irka(sys, shifts,Rt,Lt);
        sout =ones(3,size(shifts,10));

       % s0 = [0 inf 0 inf 0 inf 0 inf 0 inf];

        
    % sysr = rk(sys,s0);
        
        
        % call function crksm
        tic
        [S,R,data_out] = crksm(sys,Opts);
        toc
        %[S,R,data_out] = lyapchol(sys,s0,Opts);
        
    case 2
        % set Opts for rksm 
        Opts.method = 'rksm';
        Opts.shifts = 'cyclic';
        %Opts.rksmnorm = 'fro';
        Opts.reduction = 'twosided';
        
        % compute shifts using irka
        shifts = zeros(1,10);
        Rt = ones(m,size(shifts,2));
        Lt = ones(p,size(shifts,2));
        [~, ~, ~, s0, Rt, Lt] = irka(sys, shifts,Rt,Lt);
        
         % call function crksm
        tic
        [S,R,data_out] = crksm(sys,s0,Rt,Opts);
        %[S,R,data_out] = crksm(A,B,E,s0,Rt,Opts);
        toc
        
    case 3
        % set Opts for rksm 
        Opts.method = 'rksm';
        Opts.shifts = 'cyclic';
        %Opts.rksmnorm = 'fro';
        Opts.reduction = 'onesided';
        
        % compute shifts using irka
        shifts = zeros(1,10);
        Rt = ones(m,size(shifts,2));
        Lt = ones(p,size(shifts,2));
        [~, ~, ~, s0, Rt, Lt] = irka(sys, shifts,Rt,Lt);
        
        % call funtion
        tic
        [S,R,data_out] = rksm(A,B,C,E,s0,Opts);
        toc
        %[S,R,data_out] = lyapchol(sys,s0,Opts);
    case 4
        % compute shifts using irka
        shifts = zeros(1,10);
        Rt = ones(m,size(shifts,2));
        Lt = ones(p,size(shifts,2));
        [~, ~, ~, s0, Rt, Lt] = irka(sys, shifts,Rt,Lt);
        
        
        s0_inp = s0(1,1:6);
        s0_out = s0(1,7:12);
        tic
        [S,R,data_out] = crksm(A,B,C,E,s0_inp,s0_out,Opts);
        toc
        %[S,R,data_out] = My_rksm(sys,s0_inp,s0_out,Opts);
        
    case 5
        % compute shifts using irka
        shifts = zeros(1,10);
        Rt = ones(m,size(shifts,2));
        Lt = ones(p,size(shifts,2));
        [~, ~, ~, s0, Rt, Lt] = irka(sys, shifts,Rt,Lt);
        
        
        s0_inp = s0(1,1:6);
        s0_out = s0(1,7:12);
        Rt = Rt(1,1:6);
        Lt = Lt(1,1:6);
        tic
        [S,R,data_out] = crksm(A,B,C,E,s0_inp,s0_out,Rt,Lt,Opts);
        toc
        %[S,R,data_out] = My_rksm(sys,s0_inp,s0_out,Rt,Lt,Opts);
        
    case 6
        Opts.method = 'extended';
        s0 = ones(1,2);
        [S,data_out] = rksm(A,B,E,s0,Opts);
        
    case 7
        % compute shifts using irka
        shifts = zeros(1,10);
        Rt = ones(m,size(shifts,2));
        Lt = ones(p,size(shifts,2));
        [~, ~, ~, s0, Rt, Lt] = irka(sys, shifts,Rt,Lt);
        
        
        % Alternative f?r Shifts
        %s0(1,1) = norm(A, 'fro')/condest(A);
        %eigenwerte = -eigs(A,E);
        %s0(1,2) = eigenwerte(1,1);
        
        Opts.shifts = 'adaptive';
        tic
        [S,data_out] = crksm(A,B,E,s0,Opts);
        toc
    case 8
        % compute shifts using irka
        shifts = zeros(1,10);
        Rt = ones(m,size(shifts,2));
        Lt = ones(p,size(shifts,2));
        [~, ~, ~, s0, Rt, Lt] = irka(sys, shifts,Rt,Lt);
        
        % Alternative f?r Shifts
        %s0(1,1) = norm(A, 'fro')/condest(A);
        %eigenwerte = -eigs(A,E);
        %s0(1,2) = eigenwerte(1,1);

        Opts.Shifts = 'adaptive';
        tic
        [S,R,data_out] = crksm(A,B,C,E,s0,Opts);
        toc
    case 9
        disp('hallo');
        
     
end


% testing quality of solution (when dimension is high, comment out n>1500)
V = data_out.V_basis;
%W = data_out.W_basis;
alpha=subspace(V,Vrk);
if isempty(W)
    Pr = S'*S;
    P = V*Pr*V';
   % P  = Pr'*Pr;
else
    Pr = W*S*V';
    P  = Pr*Pr';
end
 
Y = A*P*E' + E*P*A' + B*B';
Y_norm = norm(Y);


% comparing with other hammarling-method
% Opts.method = 'hammarling';
% [Shammarling,Rhammarling] = lyapchol(sys,Opts);
% Phammarling = Shammarling*Shammarling';
% Yhammarling = A*Phammarling*E' + E*Phammarling*A' + B*B';
% Yhammarling_norm = norm(Yhammarling);



% mess-adi options
Opts.method = 'adi';
Opts.messPara = 'heur';    % only for MESS
Opts.rctol = 0;
Opts.restol = 1e-12;
Opts.norm = 2;

% call 
tic
[Sadi,Radi] = lyapchol(sys,Opts);
toc

% testing quality of solution (when dimension is high, comment out n>1500)
Padi = Sadi*Sadi';
Yadi = A*Padi*E' + E*Padi*A' + B*B';
Yadi_norm = norm(Yadi);

% save data from workspace in your path, comment the following in
% file_name = 'eady';
% save(file_name)

% genrate plot (if you want a plot of the residual norm comment the following in ;))

% plot(1:size(S,2)/m,data_out.norm_val, '-ob','LineWidth',2);
% hold on
% plot(1:size(Sadi,2)/m,Radi.res,'-+r','LineWidth',2);
% title('Convergence Rate "fom" model CRKSM ');
% xlabel('Number of Iterations');
% ylabel('Convergence Rate');
% legend('CRKSM','Location','southwest');
% legend('CRKSM','ADI','Location','southwest');


% generate a filename.tikz file for using matlab2tikz in latex,
% matlab2tikz-directory must be in yout path (I removed this directory, so you have to add it again) 
% note: cleanfigure-command is helpful!

% cleanfigure
% matlab2tikz( 'filename.tikz' );

