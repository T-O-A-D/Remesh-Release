function c = addmat(a,b,c,n,m)
% matrix A + matrix B ==> matrix C
for j = 1:m % execution loop label 1000
    for i = 1:n % execution loop label 2000
        c(i,j) = a(i,j) + b(i,j);
    end
end
return
end