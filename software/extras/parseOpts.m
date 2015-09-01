function Opts = parseOpts(Opts,Def)
%PARSEOPTS Summary of this function goes here
%   Detailed explanation goes here

fnames = fieldnames(Def); %target field names

for k = 1:length(fnames)
    if isstruct(Def.(fnames{k})) %nested structure
        % call parseOpts with the nested structure!
          % make sure the nested structure is defined 
            if ~isfield(Opts,fnames{k})
                Opts.(fnames{k}) = struct;
            end
          % call the function again
            Opts.(fnames{k}) = parseOpts(Opts.(fnames{k}),Def.(fnames{k}));
    else
        if ~isfield(Opts,(fnames{k})) || isempty(Opts.(fnames{k}))
            Opts.(fnames{k}) = Def.(fnames{k});
        end
    end
end