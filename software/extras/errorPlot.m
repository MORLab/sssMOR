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
[mag,w] = sigma(sys,w);

%% Compute frequency responses of the error
%   Compute frequency data
for iSysr = 1:nSysr
    magErr{iSysr} = sigma(sys-sysr{iSysr},w);
end

%% Plotting
nicefigure('compareMor - results (error)');

%   Create axis
for iOut = 1:sys.p
    for jIn = 1:sys.m
        ax1 = subplot(sys.p,sys.m,(iOut-1)*sys.m + jIn);
        ax1_pos = get(ax1,'Position'); % position of first axes
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','bottom',...
            'YAxisLocation','right',...
            'Color','none');
        set(ax2,'YColor','b');

        % Bode plot of original on right y axis
        lih(1) = line(w,mag2db(squeeze(mag(iOut,jIn,:))),'Parent',ax2);
        set(lih(1),'Color','b');
        set(ax2,'Xlim',[w(1),w(end)],'XTick',[],'XScale','log')
        
        % Error plots on left y axis
        for iSysr = 1:nSysr
           currErr = magErr{iSysr};
           lih(1+iSysr) = line(w,squeeze(currErr(iOut,jIn,:)),'Parent',ax1);hold on
           randColor = rand(1,3);
           set(lih(1+iSysr),'Color',randColor,'LineStyle','--');
        end
        set(ax1,'Xlim',[w(1),w(end)],'XScale','log','YScale','log')

        if iOut == 1
            if jIn == 1
                ylabel(ax1,'Error');
            elseif jIn == sys.m
                ylabel(ax2,'Magnitude /dB');
            end
        elseif iOut == sys.p && jIn == 1
            xlabel(ax1,'Frequency /rad/sec');
        end
    end
end

legNames = {'FOM'};
for iSysr = 1:nSysr, legNames = [legNames, {sprintf('ROM%i',iSysr)}]; end
legend(lih,legNames,'Location','SouthEast');