function [newShift] = getShifts(sys,varargin)


%% Create Def-struct containing default values/options
Def.purpose         = 'lyapunov';       % ['lyapunov' / 'MOR']      

%% Parsing of Inputs
if isa(varargin{end},'struct') 
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%% Compute New Shift/Shifts

switch Opts.purpose
    case 'lyapunov'
        
        
    case 'MOR' 
        if strcmp(Opts.strategy,'ADI')
            error('This method id not availabele for MOR, use ADI only when solving a Lyapunov equation');
        end
end


end

