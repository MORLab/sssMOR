function redo_bodeplot(mag,phase,omega,options)

%   this function executes the exact same code that the bode command would,
%   if called without arguments

%%  Parse the input
%   The bode function for sss yields cell arrays "mag" and "phase". 
%   The bode function (built-in) returns 3D arrays...
%
%   However,
%   analyze_HFM saves the data into normal vectors. Parse the input to make
%   the function robust

%   this code transforms the data to 1D cell arrays!
if iscell(mag) %sss-bode
    mag=cellfun(@(x) 20*log10(x), mag, 'UniformOutput',false);
elseif length(size(mag))==3 %ss-bode
    mag = {reshape(mag,size(mag,3),1)};
    mag = cellfun(@(x) 20*log10(x), mag, 'UniformOutput',false);
    
    phase = {reshape(phase,size(phase,3),1)};
    phase = cellfun(@(x) unwrap(x), phase, 'UniformOutput',false);
else %if mag is a vector (e.g. coming from some other function)
    mag = {mag2db(mag)};
end

% if ~iscell(phase)
%     phase = {phase};
% end

%   if no options were passed, create an empty cell
if ~exist('options','var')
    options = {};
end


%%   Normal computations
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