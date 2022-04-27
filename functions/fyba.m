function fyba = fyba(xq,yq,ax,ay,alph)

    fyba = (ay*(ax*xq+ay*yq)/alph)-ax*(ay*xq-ax*yq);
    
    return
end