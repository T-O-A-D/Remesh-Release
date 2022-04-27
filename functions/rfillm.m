function a = rfillm(a,n,m,value)
    for i = 1:n
        for j = 1:m
            a(i,j) = value;
        end
    end
    
    return
end