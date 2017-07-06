function [AllOpts,nCases] = generateAllOpts(AllOpts)


%%  Extract fields and values from AllOpts
fields = fieldnames(AllOpts).';
vargin = cell(1,length(fields));

for iField = 1:length(fields)
    vargin{iField} = AllOpts.(fields{iField});
end
    
%% Generate all combinations using "allcomb"
% MATLAB central file exchange
% https://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin-
% downloaded 23 Nov 2016

P = allcomb(vargin{:});
nCases = size(P,1);

%% Construct AllOpts, a cell of Opts structures
AllOpts = cell(1,nCases);
for iCase = 1:nCases
    AllOpts{iCase} = cell2struct(P(iCase,:), fields, 2);
end
