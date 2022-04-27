function [coord, i3] = areach(npoin,nelem,intmat,coord)
 
    [r,c] = size(coord);
    temp = zeros(2, npoin);
    temp(:,1:c) = coord;
    coord = temp;

    % determinant function
    isucc = 0;
    for ie = 1:nelem % execution loop label 1000
        
        i1 = intmat(1,ie);
        i2 = intmat(2,ie);
        i3 = intmat(3,ie);
        
        x1 = coord(1,i1);
        x2 = coord(1,i2);
        x3 = coord(1,i3);
        y1 = coord(2,i1);
        y2 = coord(2,i2);
        y3 = coord(2,i3);

        dtmp = deter(x1,y1,x2,y2,x3,y3);

        if dtmp > 0.0 % jump gate label 1000
            continue
        else
            fprintf('Area error in element %f \n',ie)
            fprintf('node 1 = %f, x1 = %f, y1 = %f \n',i1,x1,y1)
            fprintf('node 2 = %f, x1 = %f, y1 = %f \n',i2,x2,y3)
            fprintf('node 3 = %f, x1 = %f, y1 = %f \n',i2,x3,y3)
            isucc = 1;
        end
        %jump gate exit label 1000
    end

    if isucc == 0
        fprintf('The checking was successful \n')
    end

    % output number of nodes and elements 
    fprintf('Total number of generated points : %f \n',npoin)
    fprintf('Total number of generated elements : %f \n',nelem)
    
    return
end