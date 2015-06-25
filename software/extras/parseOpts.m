function opts = parseOpts(opts,def)
%PARSEOPTS Summary of this function goes here
%   Detailed explanation goes here

fnames = fieldnames(def); %target field names

for k = 1:length(fnames)
    if isstruct(def.(fnames{k})) %nested structure
        % call parseOpts with the nested structure!
          % make sure the nested structure is defined 
            if ~isfield(opts,fnames{k})
                opts.(fnames{k}) = struct;
            end
          % call the function again
            opts.(fnames{k}) = parseOpts(opts.(fnames{k}),def.(fnames{k}));
    else
        if ~isfield(opts,(fnames{k})) || isempty(opts.(fnames{k}))
            opts.(fnames{k}) = def.(fnames{k});
        end
    end
end