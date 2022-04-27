function fxba = fxba(xq,yq,ax,ay,alph)

    fxba = (ax*(ax*xq+ay*yq)/alph)+ay*(ay*xq-ax*yq);
    
    return
end