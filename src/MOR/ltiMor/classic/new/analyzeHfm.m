function HFM_data = analyze_HFM(path2model)

%   This function is used to analyze a given HFM before conducting MOR
%
%   "path2model" is the absolute path to the .mat file containing the model

%%  initialize
fprintf(1,'Function "analyze_HFM" started...\n\n')

%if no input was chosen, load the building model
if ~exist('path2model','var') || isempty(path2model)
%     path2model = 'P:\MOR\Matlab\Benchmark\SLICOT\build.mat';
      path2model = 'P:\MOR\Matlab\Benchmark\SLICOT\CDplayer.mat';
end

%%  load the model
fprintf(1,'Loading the model from "%s"... \n',path2model);
load(path2model)

%%  make all sorts of computations required to identify the model type
% 1st order, 2nd order, explicit/implicit form, DAE...

if exist('A','var')
    fprintf(1,'The model is in 1st order form. \n');
    N = size(A,1);
    fprintf(1,'The model has the state space dimension of %i. \n',N);
    fprintf(1,'The model has %i input(s) and %i output(s). \n',size(B,2),size(C,1));
end
% ...tbd

cond_tol = 1e10;
if exist('E','var')
    cond_E = condest(E);
    if cond_E > cond_tol;
        warning('Matrix E has a condition number of at least %4.1e and might be therefore singular. \n',cond_E);
    end
else
    E = speye(size(A));
end

if ~exist('D','var') || D == zeros(size(C,1),size(B,2)); %no feedthrough
    D = zeros(size(C,1),size(B,2));
    fprintf(1,'The model has no feedthrough. \n',size(A,1))
end

%%  Plot sparsity pattern
fh.sparse = figure('Name','Sparsity pattern of HFM','NumberTitle','off');
subplot(1,2,1)
    spy(A);title('spy(A)');
subplot(1,2,2)
    spy(E);title('spy(E)');
    

%%  Create an sss-object of the model
sys = sss(A,B,C,D,E);
HFM_data.sys = sys;

%%  Create a bode plot of the model

[mag,phase,w] = bode(sys); %compute the data for the bode plots

HFM_data.bode.mag = mag{1}; %store the data for later computation
HFM_data.bode.phase = phase{1};
HFM_data.bode.omega = w;

fh.bode = figure('Name','Bode Plot of HFM','NumberTitle','off');
options ={'Color',TUM_Blau};
bode_plot(mag,phase,w,options)

% %%  Plot the eigenvalues of the system
% 
% %   Since the computation of generalized eigenvalues for descriptor systems
% %   (= singularly implicit) is tricky, do this step only for explicit or
% %   regularly implicit systems.
% %
% %   -> FIND AN ALGORITHM SUITED FOR DAEs
% EigsLegend = {};
% 
% if cond_E < cond_tol 
%     if N < 1e4 %HFM might be small enough for direct pole computation
%         fh.eig = figure('Name','Pole Plot of HFM','NumberTitle','off');
%         poles = eig(sys);
% 
%         HFM_data.eig = poles;
% 
%         lh.eig = plot(real(poles),imag(poles),'*','Color',TUM_Blau);
%         EigsLegend = {EigsLegend{:},'eig'};
%         xlabel('Re(s)');ylabel('Im(s)');
%     else    %  Compute (some) eigenvalues through indirect methods
%         n_eigs = 100;
%         lm_eigs = eigs(A,E,n_eigs,'lm');
% 
%         HFM_data.eigs.lm = lm_eigs;
% 
%         fh.eigs = figure('Name','Pole Plot (eigs) of HFM','NumberTitle','off');
%         lh.lm_eigs = plot(real(lm_eigs),imag(lm_eigs),'*','Color',TUM_Blau);
%         EigsLegend = {EigsLegend{:},'eigs(lm)'};
%         hold on
% 
%         try
%             sm_eigs = eigs(A,E,n_eigs,'sm');
% 
%             HFM_data.eigs.sm = sm_eigs;
%             lh.sm_eigs = plot(real(sm_eigs),imag(sm_eigs),'*','Color',TUM_Orange);
%             EigsLegend = {EigsLegend{:},'eigs(sm)'};
%         catch err
%             fprintf(1,'Eigs computation for "smallest magnitude" option failed. \n')
%         end
%         try
%             lr_eigs = eigs(A,E,n_eigs,'lr');
% 
%             HFM_data.eigs.lr = lr_eigs;
%             lh.lr_eigs = plot(real(lr_eigs),imag(lr_eigs),'*','Color',TUM_Gruen);
%             EigsLegend = {EigsLegend{:},'eigs(lr)'};
%         catch err
%             fprintf(1,'Eigs computation for "largest real" option failed. \n')
%         end
%     end
%     legend(cell2mat(struct2cell(lh)),EigsLegend','Location','EastOutside');
% end




%%  Plot the hankel singular values

if exist('hsv','var')
    hankelSV = load(path2model,'hsv'); %overwrite name to avoid conflict with colormap
    hankelSV = hankelSV.hsv;
    
    fh.hsv = figure('Name','Hankel singular values of HFM','NumberTitle','off');
    lh.hsv = plot(hankelSV,'*','Color',TUM_Blau);
%     ylim([0,3e-3]);
    ylabel('hsv');
% else
%     fprintf(1,'[!] Hankel singular values were not stored in the ".mat" file. \n')
%     decision = input('Would you like to compute the hsv? (y/n)','s');
%     k = 0;
%     while k
%         if strcmp(decision,'y')
%             fprintf(1,'Starting computation of hsv... \n')
%             tic;
%             hsv = hankelsv(tf(sys));
%             time = toc;
%             fprintf(1,'Hsv computed in %4.1f s. \n',time)
%             
%             fh.hsv = figure('Name','Hankel singular values of HFM','NumberTitle','off');
%             lh.hsv = plot(hsv,'*','Color',TUM_Blau);
%             ylabel('hsv');
%     
%             k = 1;
%         elseif strcmp(decision,'n')
%             k = 1;
%         else
%             decision = input('INVALID ANSWER. Would you like to compute the hsv? (y/n)','s');
%         end
%     end
end


%%  Test, if the bases of the Krylov subspace can be computed

HFM_data.LSE_test = 0;
try
    s0 = 100*rand; %real, positive shift between 0 and 100
    tic
    v1 = (A-s0*E)\B;
    t_lse = toc;
    fprintf('Computation of the 1st Krylov direction with real shift in %4.1f s.\n',t_lse);
    
    try
        tic
        v2 = (A-s0*E)\(E*v1);
        t_lse = toc;
        fprintf('Computation of the 2nd Krylov direction with real shift in %4.1f s.\n',t_lse);
        
        try
            s0 = 100*rand + 100i*rand; %complex
            tic
            v3 = (A-s0*E)\B;
            t_lse = toc;
            fprintf('Computation of the 1st Krylov direction with complex shift in %4.1f s.\n',t_lse);
            
            try
                tic
                v4 = (A-s0*E)\(E*v3);
                t_lse = toc;
                fprintf('Computation of the 2st Krylov direction with complex shift in %4.1f s.\n',t_lse);
                
                
                HFM_data.LSE_test = 1;
                fprintf('The LSE solution test is passed! Proced with Krylov-MOR! \n',t_lse);
            catch err
                fprintf('Computation of the 2st Krylov direction with complex shift failed. \n');
                rethrow(err);
            end
            
        catch err
            fprintf('Computation of the 1st Krylov direction with complex shift failed. \n');
            rethrow(err);
        end
        
    catch err
        fprintf('Computation of the 2nd Krylov direction with real shift failed. \n');
        rethrow(err);
    end
    
catch err
    fprintf('Computation of the 1st Krylov direction with real shift failed. \n');
    rethrow(err);
end



function bode_plot(mag,phase,omega,options)

%   this function executes the exact same code that the bode command would,
%   if called without arguments

mag=cellfun(@(x) 20*log10(x), mag, 'UniformOutput',false);
[m, p] = size(mag); % inputs, outputs
plot_handles=zeros(2*m,p);

for i=1:m %for each input
    for j=1:p %for each output
        plot_handles(2*i-1,j)=subplot(2*m,p,2*((i-1)*p)+j);
        hold on        
%         if j==1
%             y_lab=sprintf('To Out(%i)',ceil(i/2));
%             ylabel(y_lab)
%         end
        if i==1 && m>1
            x_lab=sprintf('From In(%i)',j);
            title(x_lab)
        end
      
        % amplitude
        y_plot=mag{i,j};
        plot(omega,y_plot,options{:})
        set(gca, 'XScale', 'log');
        set(gca, 'XLim', [min(omega) max(omega)]);
        box on
        mx=max(mag{i,j}); mn=min(mag{i,j});
        set(gca, 'YLim', [mn-(mx-mn)/20, mx+(mx-mn)/20]);
        ylabel('Magnitude [dB]') 

        plot_handles(2*i,j)=subplot(2*m,p,(2*i-1)*p+j);
        hold on
        box on
        % phase
        y_plot=phase{i,j};
        plot(omega,y_plot,options{:})
        set(gca, 'XScale', 'log');
        set(gca, 'XLim', [min(omega) max(omega)]);
        ylabel('Phase [deg]');
        xlabel('Frequency [rad/sec]');
    end
end

% figure and axes properties
set(gcf,'UserData',plot_handles)
set(gcf,'Color', [1 1 1])

% link the axes for simultaneous zooming
linkaxes(reshape(plot_handles,2*m*p,1),'x');

