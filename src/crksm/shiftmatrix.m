function [s_ma,tangential] = shiftmatrix(s_inp,m,varargin)
% This function ensures that for every single shift the same number of
% moments is matched especially if there are diffrent shift vectors for the
% input and output Krylov space or if there are complex shifts.
% 
% 1. Example: only real shifts: s_inp = [0 1 2 3 4;   s_out = [3 5 6 7 9;
%                                        2 1 2 2 3]            1 1 1 1 1]
% 
% As you can see the input space V consists of more directions as the
% outputspace W which leads to an error. To avoid this the number of higher
% moments of the shift with the lower number of higher moments of one pair 
% is increased so that both numbers are equal, see the following:
% 
%       First pair of shifts:   s_inp(1,1) = 0                   s_out(1,1) = 3
%                               --> number of moments: 2         --> number of moments = 1 
%                               --> unequal number of directions in V and W                            
%                               --> adding one higher moment to shift s_out = 3                             
% 
% Filling up the shift vectors in this way is possible because one just
% works with a Taylor approximation of higher order
% 
% If there are complex shifts note that one gets for evry complex shift and
% for every higher order two directions. Just imagine the first shift in
% the above example is complex, then there are there are 4 directions for
% the first shift in s_inp and only one for the first shift in s_out
% --> 3 higher order moments have to be added in s_out(2,1).

global hermite withoutC Rt Lt s_out 

if length(varargin) == 1
    p = varargin{1};
elseif length(varargin) == 2
    p = varargin{1};
    A = varargin{2};
else
    p = 1;
end

% transpose shift vectors and tangential directions
if size(s_inp,1) >= 2 && size(s_inp,2) == 1 
    s_inp = s_inp';
elseif exist('s_out','var') && ~isempty(s_out) && size(s_out,1) >= 2 && size(s_out,2) == 1
    s_out = s_out';
end

% build matrix containing all the information of shifts, tangential directions
% preallocate memory
if hermite == 1 && withoutC ==1
    s_info = zeros(3+m,size(s_inp,2));
    clear Lt
elseif hermite == 1 && withoutC == 0
    s_info = zeros(2+m+p,size(s_inp,2));
else
    s_info = zeros(4+m+p,size(s_inp,2));
end

% read in input shifts and order of input shifts
if size(s_inp,1) == 2
    s_info(1,:) = s_inp(1,:);
    s_info(2,:) = s_inp(2,:);
else
    s_info(1,:) = s_inp;
    s_info(2,:) = ones(1,size(s_inp,2));
end

% read in tangential directions
if (~exist('Rt','var') || isempty(Rt)) && (~exist('Lt','var') || isempty(Lt))
    s_info(3:2+m,:) = ones(m,size(s_inp,2));
    s_info(3+m:2+p+m,:) = ones(p,size(s_inp,2));
    tangential = false;
elseif exist('Rt','var') && (~exist('Lt','var') || isempty(Lt))
    s_info(3:2+m,:) = Rt;
    s_info(3+m:2+p+m,:) = ones(p,size(s_inp,2));
    tangential = true;
elseif exist('Rt','var') && exist('Lt','var') && (~isempty(Rt) && ~isempty(Lt))
    s_info(3:2+m,:) = Rt;
    s_info(3+m:2+p+m,:) = Lt;
    tangential = true;
end

% read in output shifts s_out and order of them
if hermite == 1
    s_info(3+m+p,:) = s_inp(1,:);
    s_info(4+p+m,:) = s_inp(2,:);
elseif hermite == 0
    s_info(3+m+p,:) = s_out(1,:);
    if size(s_out,1) == 2
        s_info(4+p+m,:) = s_out(2,:);
    else
        s_info(4+p+m,:) = ones(1,size(s_inp,2));
    end
else
    s_info(3+p+m,:) = ones(1,size(s_inp,2));
    s_info(4+p+m,:) = s_info(3+p+m,:);
end

% neccessary in order to avoid numerical errors
s_info = single(s_info);

% set variables
s_in = s_info;
% if size(s_info,1) > (2+m+p)
%     s_out = s_info(3+m+p,:);
% end

% if there are complex shifts, try to pair them, if this does not work, add
% the conjugate complex part
if ~isreal(s_in(1,:))
    % look for complex pairs
    try
        %s_in = cplxpair(s_in);
        s_in = cplxpair(s_in(1,:));
    catch
        % add conjugate complex elements 
        for ii=1:1:size(s_in,2)
            if ~isreal(s_in(1,ii))
                counter = 1;
                while counter ~= 0
                    if s_in(1,ii) == conj(s_in(1,counter))
                        counter = 0;
                    else
                        counter = counter+1;
                    end
                    if counter > size(s_in,2)
                        counter = 0;
                        s_in = [s_in conj(s_in(:,ii))];
                    end
                end   
            end   
        end
        clear counter
    end
    
    %if output == 1 && no_complex_out_s == 1
    if exist('s_out','var') && ~isreal(s_out) && ~isempty(s_out)
        try
            s_out = cplxpair(s_out);
            if size(s_in,2) ~= size(s_out,2)
                lastwarn('try to use a shift vectors containing complex pairs of shifts if there are complex shifts')
            end
        catch
            % add conjugate complex elements
            for ii=1:1:size(s_out,2)
                if ~isreal(s_out(1,ii))
                    counter = 1;
                    while counter ~= 0
                        if s_out(1,ii) == conj(s_out(1,counter))
                            counter = 0;
                        else
                            counter = counter+1;
                        end
                        if counter > size(s_out,2)
                            counter = 0;
                            s_out = [s_out conj(s_out(1,ii))];
                        end
                    end   
                end   
            end
            %s_out = cplxpair(s_out);
            if size(s_in,2) ~= size(s_out,2)
                lastwarn('try to add the following shifts to your vector:')
            end
            clear counter   
       end
    end
    
    % build expanded shift information matrix
    % preallocate memory
    if ~exist('s_out','var') || isempty(s_out)
        s_info_new = zeros(2+m+p,size(s_in,2));
    else
        s_info_new = zeros(4+m+p,size(s_in,2));
    end
    
    for ii=1:1:size(s_in,2)
        counter = 1;
        while counter ~= 0
            if s_in(1,ii) == s_info(1,counter)
                s_info_new(1,ii) = s_in(1,ii);
                if ~exist('s_out','var') || isempty(s_out)
                    s_info_new(2:2+p+m,ii) = s_info(2:2+p+m,counter);
                else
                    s_info_new(2:4+p+m,ii) = s_info(2:4+p+m,counter);
                end
                counter = 0;
            elseif s_in(1,ii) ~= s_info(1,counter) && counter == size(s_info,2)
                s_info_new(1,ii) = s_in(1,ii);
                s_info_new(2,ii) = s_info_new(2,ii-1);       % fill up order
                for jj=1:1:size(s_info,2)
                    if s_in(1,ii) == conj(s_info(1,jj))
                        s_info_new(3:2+m,ii) = conj(s_info(3:2+m,jj));   
                    end
                end
%                 s_info_new(3:2+m,ii) = ones(m,1);            % fill up Rt 
                 s_info_new(3+m:2+m+p,ii) = ones(p,1);        % fill up Lt
                if ~exist('s_out','var') || isempty(s_out)
                    %add_shifts = -eigs(A); % spaltenvektor
                    %[~, ~, ~, s0, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = irka(sys, add_shifts);
                    s_info_new(3+m+p:4+m+p,ii) = ones(2,1);  % fill up output
                end
                counter = 0;
            elseif s_in(1,ii) ~= s_info(1,counter) || counter <= size(s_info,2)
                counter = counter+1;
            end  
        end  
    end
    s_info = s_info_new;
    clear counter s_info_new
    % bis hier her geht es
    
    % das ist daf?r da, um die Ordnung anzugleichen
    for ii=1:1:size(s_in,2)
        if ~isreal(s_info(1,ii)) && imag(s_info(1,ii)) < 0 && ii < size(s_in,2) && s_info(2,ii) >= s_info(2,ii+1)
            s_info(2,ii+1) = s_info(2,ii);
        elseif ~isreal(s_info(1,ii)) && imag(s_info(1,ii)) < 0 && ii < size(s_in,2) && s_info(2,ii) < s_info(2,ii+1)
            s_info(2,ii) = s_info(2,ii+1);          
        end      
    end
    
    % compare oder of input and output shifts
    if exist('s_out','var') && ~isempty(s_out)
        for ii=1:1:size(s_info,2)
            if s_info(2,ii) >= s_info(4+m+p,ii)
                s_info(4+m+p,ii) = s_info(2,ii);
            else
                s_info(2,ii) = s_info(4+m+p,ii);
            end        
        end
    end
   
% only real shifts
else
    if exist('s_out','var') && ~isempty(s_out)
        for ii=1:1:size(s_info,2)
            if s_info(2,ii) >= s_info(4+m+p,ii)
                s_info(4+m+p,ii) = s_info(2,ii);
            else
                s_info(2,ii) = s_info(4+m+p,ii);
            end        
        end
    end
end

% create input shift matrix (contains shifts and order of each shift) and
% build shift vector
s_inp = s_info(1:2,1:size(s_info,2));
s_inp = make_vect(s_inp);

% create right tangential directions matrix (contains directions and order of each shift) and
% build tangetial vector/matrix
Rt(1:m,1:size(s_info,2)) = s_info(3:2+m,1:size(s_info,2));
Rt(m+1,1:size(s_info,2)) = s_info(2,1:size(s_info,2));
Rt = make_vect(Rt,m);

% create left tangential directions matrix (contains directions and order of each shift) and
% build tangetial vector/matrix
Lt(1:p,1:size(s_info,2)) = s_info(3+m:2+m+p,1:size(s_info,2));
Lt(p+1,1:size(s_info,2)) = s_info(2,1:size(s_info,2));
Lt = make_vect(Lt,p);

% eventually create output shift vector
if hermite == 0
    s_out = s_info(3+m+p:4+m+p,1:size(s_info,2));
    s_out = make_vect(s_out);
end

% collect all shift, direction vectors in one matrix
s_ma(1,:) = s_inp;       
s_ma(2:m+1,:) = Rt;
s_ma(2+m:1+m+p,:) = Lt;
if hermite == 0
    s_ma(2+m+p,:) = s_out;
end
clear s_in s_info s_inp
clearvars -global Rt Lt s_out
end



function x = make_vect(x,y)
% function which transforms a matrix with not more than 2 rows and 1 column or 2 columns and 1 row in a vector
% or just transposes a vector
% function also sorts the vector in an increasing way
% Input: - shift vector s_inp / s_out
%        - tangential directions Rt an Lt
% Output: sorted vectors, see example
%
% Example: x = [0 1 2 3 4; ---> x = [0 1 1 2 2 2 3 4 4]
%               1 2 3 1 2]


if size(x,1) == 2
    temp = zeros(1,sum(x(2,:)));
    for ii=1:1:size(x,2)
        k = sum(x(2,1:(ii-1)));
        k = (k+1):(k+x(2,ii));
        temp(k) = x(1,ii)*ones(1,x(2,ii));
    end
    x = temp;
else
    temp = zeros(y,sum(x(y+1,:)));
    for ii=1:1:size(x,2)
        k = sum(x(y+1,1:(ii-1)));
        k = (k+1):(k+x(y+1,ii));
        temp(1:y,k) = x(1:y,ii)*ones(1,x(y+1,ii));
    end
     x = temp;
end

end


function x = mysort(x)
n = size(x,2);
while n>1
    newn = 1;
    for ii=1:1:n-1
        if real(x(1,ii)) > real(x(1,ii+1)) 
            temp = x(:,ii);
            x(:,ii) = x(:,ii+1);
            x(:,ii+1) = temp;
            newn = ii+1;
        end
        if x(1,ii) == conj(x(1,ii+1)) && imag(x(1,ii)) > 0
            temp = x(:,ii);
            x(:,ii) = x(:,ii+1);
            x(:,ii+1) = temp;
        end
        if isnan(real(x(1,ii))) && ~isnan(real(x(1,ii+1)))
            temp = x(:,ii);
            x(:,ii) = x(:,ii+1);
            x(:,ii+1) = temp;
            newn = ii+1;
        end
    end
    n = newn;
end
end



%                c_rs    = BstBs\(B_s'*(AV-EV*(Er\Ar)));
               
               
%                FHeiko = EV*((Er\Br) + (S'*S)*c_rs');
               
%                BstF = B_s'*FHeiko;
%                FtF  = FHeiko'*FHeiko;
               
%                R_normHeiko = max(abs(eig(full([BstBs, BstF; BstF', FtF] * [eye(m), eye(m);eye(m), zeros(m,m)]))));
%                 norm_APE = sqrt(max(abs(eig( (AV'*AV) * Pr * (EV'*EV) * Pr ))))
%                 norm_Pr  = max(abs(eig( full(V'*V*Pr) )));









% build and orthogonalize new basis
               if ~isreal(vnew) 
                   % if new direction vnew is complex, now after expanding
                   % the basis with the real part it follows the complex part
                   % note: until now no further shift is read out!
                   if index == ii-1
                       %vnew_imag = imag(vnew);
                       V = [V imag(vnew)];
                       
                       if withoutC == 0 && ~isreal(wnew)
                           wnew_imag = imag(wnew);
                           W = [W wnew_imag];
                       
                       % here a higher order is calculated with the current s_out because s_out is real and not complex
                       elseif withoutC == 0 && isreal(wnew)
                           % calculate additional direction wnew by determining a higher order of s_out_j

                           Opts.reuseLU = 1;
                           rhsC = W(:,jbasis+m:jnew-1+m);   % new rhs achte hier auf +m
                           if hermite == 1 && tangential == 1
                               rhsC = rhsC*jCol_Lt;
                           end
                           if hermite == 1
                               Anew_W = (A-jCol_inp*E)';
                           else
                               Anew_W = (A-jCol_out*E)';
                           end
                            
                           [wnew] = solveLse(Anew_W,E'*rhsC,Opts);
                           
                           % if the higher order of wnew is not real use the real part
                           if ~isreal(wnew)
                               wnew = wnew(:,m);   wnew = real(wnew);
                               W = [W wnew];
                           else
                               W = [W wnew];
                           end  
                       end
                       % reset variables
                       wait = false;
                       index = 0;
                   else
                       % if new direction vnew is complex, first expand basis with the real part of vnew
                       % hier gehe ich immer zuerst rein
                       %vnew_real = real(vnew);
                       V = [V real(vnew)];
                       
                       if withoutC == 0 && ~isreal(wnew)
                           wnew_real = real(wnew);
                           W = [W wnew_real];
                       elseif withoutC == 0 && isreal(wnew)
                           W = [W wnew];
                       end
                       
                       
                       % some settings for the code
                       wait = true;
                       index = ii;
                   end
               else
                   % if vnew is real
                   if wait == 0
                   V = [V vnew];
                   end
                   
                   % if wnew is complex
                   if withoutC == 0 && ~isreal(wnew)

                       if index == ii-1
                           % if new direction wnew is complex, now after expanding
                           % the basis with the real part it follows the complex part of wnew
                           % note: until now no further shift is read out!
                           wnew_imag = imag(wnew);
                           W = [W wnew_imag];
                           
                           % calculate additional direction vnew with the same shift but of higher order
                           Opts.reuseLU = 1;
                           rhsB = V(:,jbasis+m:jnew-1+m);    % new rhs for B
                           Anew_V = (A-jCol_inp*E);
                           if tangential == 1
                               rhsB = rhsB*jCol_Rt;
                           end
                           
                           [vnew] = solveLse(Anew_V,E*rhsB,Opts);

                           if ~isreal(vnew)
                               vnew = vnew(:,m);   vnew = real(vnew);
                               V = [V Vnew];
                           else
                               V = [V vnew];
                           end 
                           
                           % settings
                           index = 0;
                           
                       else
                           % if new direction wnew is complex, first expand basis W with the real part of wnew
                           % hier gehe ich immer zuerst rein
                           wnew_real = real(wnew);
                           W = [W wnew_real];
                       
                           % some settings for the code
                           wait = true;
                           index = ii;   
                       end
                       
                   elseif withoutC == 0 && isreal(wnew)
                      W = [W wnew];    
                   end
               end 
               
               
               
               
               % set indices for calculation of reduced system
               if tangential == 1
                   jbasis    = size(V,2)-1;        % index for the first column of the directions calculated before   
                   jnew      = jbasis+1;        
                   jnew_last = jnew;             % index of the last new column     
               else
                   jbasis    = (ii-2)*m+1;     
                   jnew      = jbasis+m;
                   jnew_last = jnew+m-1;         % index of the last new column
                
               end