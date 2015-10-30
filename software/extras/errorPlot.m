function errorPlot(sys,sysr,w)

%% Parse input
if ~exist('w','var')
    w = [];
end

nSysr = length(sysr);
if nSysr == 1 && ~iscell(sysr)
    sysr = {sysr};
end
magErr = cell(nSysr,1);

%% Compute frequency response of the original system sys
[mag, ~, w] = bode(sys,w);

%% Compute frequency responses of the error
%   Compute frequency data
for iSysr = 1:nSysr
    [magErr{iSysr},~] = bode(sys-sysr{iSysr},w);
end

%% Plotting
nicefigure('compareMor - results (error)');

%   Create axis
ax1 = axes;
ax1_pos = get(ax1,'Position'); % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none');
set(ax2,'YColor','b');

% Bode plot of original on right y axis
lih(1) = line(w,squeeze(mag),'Parent',ax2);
set(lih(1),'Color','b');
set(ax2,'XTick',[],'XScale','log')

% Error plots on left y axis
for iSysr = 1:nSysr
   lih(1+iSysr) = line(w,squeeze(magErr{iSysr}),'Parent',ax1);hold on
   randColor = rand(1,3);
   set(lih(1+iSysr),'Color',randColor,'LineStyle','--');
end
set(ax1,'XScale','log','YScale','log')

ylabel(ax1,'Error');
ylabel(ax2,'Magnitude /dB');

xlabel(ax1,'Frequency /rad/sec');
legend(lih,'Location','SouthEast');