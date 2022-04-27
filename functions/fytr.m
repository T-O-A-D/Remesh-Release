function fytr = fytr(xq,yq,ax,ay,alph)

    fytr = ay*alph*(ax*xq+ay*yq)-ax*(ay*xq-ax*yq);

    return
end