function fxtr = fxtr(xq,yq,ax,ay,alph)

    fxtr = ax*alph*(ax*xq+ay*yq)+ay*(ay*xq-ax*yq);
    
    return
end