

function spy3(A)

A=full(A);
figure

n = size(A);

k = norm(A);

for i=1:size(A,1)
    for j=1:size(A,2)
        if A(i,j)>0
            plot(j, n(1)-i, '.', 'Markersize', abs(A(i,j))/k*100 )
        elseif A(i,j)<0
            plot(j, n(1)-i, 'r.', 'Markersize', abs(A(i,j))/k*100 )
        end
        hold on
    end
end


