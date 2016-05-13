
function w = gewichte(refs, p)

w=zeros(size(refs));
% for i = 1:1:length(refs) %[1 length(refs)]
%     if refs{i} == p
%         w=zeros(size(refs));
%         w(i)=1;
%         return
%     end
%     w(i) = sqrt(norm(refs{i}-p)^(-2));
%     %w(i) = abs(refs{i}-p)^(-1);
% end
%w = w/norm(w,1);

% Lineare Interpolation zwischen 2 Stützstellen
temp = 0;
for i = 1:1:length(refs) 
    if refs(i) == p
        w=zeros(size(refs));
        w(i)=1;
        return
    elseif refs(i) < p
        temp = i;
    end 
end
w(temp) = 1 - 1/(refs(temp+1)-refs(temp))*(p-refs(temp));
w(temp+1) = 1/(refs(temp+1)-refs(temp))*(p-refs(temp));



