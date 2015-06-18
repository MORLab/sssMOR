function [sysr, V, W, alpha_opt] = RK_ICOP(sys, s0, varargin)
% ICOP Algorithm for Model Order Reduction by Krylov Subspace Methods
% ------------------------------------------------------------------
% [sysr, V, W, aopt] = RK_ICOP(sys, [s0;q], varargin)
% Inputs:       * sys: an sss-object containing the LTI system
%               * [s0;q]:   s0: initial expansion point
%                            q: order of reduced system
%    [optional] * 'projection': 'in' or 'out' 
%               * 'epsilon': maximum step size for convergence
%               * 'maxIter': maximum number of iterations
% Outputs:      * sysr: reduced system
%               * V, W: Projection matrices spanning Krylov subspaces
%               * alpha_opt: optimal expansion point
% ------------------------------------------------------------------
% Example of usage:
% [sysr,V,W,alpha] = RK_ICOP(sys,[0;10],'projection','in','epsilon',1e-3);
% ------------------------------------------------------------------
% For further explanation of the notation of [s0;q], please refer to
% <a href="matlab:help arnoldi">help arnoldi</a>.
%
% For details and extensibilities see the references:
% [1] R. Eid, H. Panzer and B. Lohmann: How to choose a single expansion
%   point in Krylov-based model reduction? Technical reports on
%   Automatic Control, vol. TRAC-4(2), Institute of Automatic Control, 
%   Technische Universitaet Muenchen, Nov. 2009. 
% [2] R. Eid: Time domain Model Reduction by Moment Matching, Ph.D thesis,
%   Institute of Automatic Control, Technische Universitaet Muenchen, 2009.
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Rudy Eid
% Last Change:  03 Feb 2011
% ------------------------------------------------------------------

projMode = 'Two-Sided Arnoldi';

epsilon = 1e-2;
maxIter = 20;

alpha_opt=s0(1);

s0_out=[alpha_opt;s0(2)];
s0_inp=[alpha_opt;s0(2)];

if nargin>3     % Eingangsparameter verarbeiten
    c=1;
    while c<=length(varargin)
        try
            switch lower(varargin{c})
            case 'projection',
                switch  varargin{c+1}
                case 'out',
                    projMode = 'One-Sided Arnoldi: Output Krylov Subspace';
                    s0_inp=[];
                case 'in',
                    projMode = 'One-Sided Arnoldi: Input Krylov Subspace';
                    s0_out=[];
                end
                c=c+2;
            case 'epsilon',
                epsilon = varargin{c+1};
                c=c+2;
            case 'maxiter',
                maxIter = varargin{c+1};
                c=c+2;
            otherwise
                c=c+1;
            end
        catch ex
             disp(ex.message);
             disp(varargin{c});
             c=c+1;
        end
    end
end

disp('--- Starting RK-ICOP ---')
disp(['  ' projMode])

i=0;
while 1
    i=i+1;

    disp(['  Step ' num2str(i) ': ']);
  
    alpha_opt_alt=alpha_opt;
    if ~isempty(s0_inp)
        %input subspace
        s0_inp=[alpha_opt;s0(2)];
    end
    if ~isempty(s0_out)
        %output subspace
        s0_out=[alpha_opt;s0(2)];
    end
    
    % calculate reduced system
    [sysr, V, W] = RK(sys, s0_inp, s0_out);
    
    try
        % calculate aopt
        alpha_opt = RK_OP(sysr);

        if isnan(alpha_opt)   % abort on error
            break
        elseif alpha_opt<0 || imag(alpha_opt)~=0
            disp('    Warning: RK_OP returns complex or negative result.')
            alpha_opt = abs(alpha_opt);   %using abs value instead
        end
    catch ex
        alpha_opt = NaN;  % Schleife verlassen bei Fehler
        disp(['Error: ' ex.message '. Aborting.'])
        break
    end
    disp(['  aopt = ' num2str(alpha_opt)]);
    
    if abs(alpha_opt-alpha_opt_alt)/alpha_opt <= epsilon        % abort on convergence
        disp(['Converged in step ' num2str(i)] )
        break
    end
    if i>=maxIter
        disp(['Not converged after ' num2str(i) ' steps.'])
        break
    end
end

disp('--- Leaving RK-ICOP ---')
