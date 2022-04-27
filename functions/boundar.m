function lcoor = boundar(npoin,nelem,lcoor,intmat)
    %      THIS SUBROUTINE FINDS OUT THE BOUNDARY POINTS
    %      LCOOR(IPOIN).EQ.0 --> INTERIOR
    %      LCOOR(IPOIN).NE.0 --> BOUNDARY
    % 

    %      INITIALIZE
    for ip = 1:npoin % execution loop label 1000

        lcoor(ip) = 0;

        % loop over the elements
        for ie = 1: nelem % execution loop label 2000
            for in = 1:3 % execution loop label 2001

                in1 = in + 1;

                if (in1 > 3) 
                    in1 = in1 - 3;
                end

                in2 = in + 2;

                if in2 > 3
                    in2 = in2 - 3;
                end

                ip = intmat(in,ie);

                lcoor(ip) = lcoor(ip) + intmat(in2,ie) - intmat(in1,ie);
                
            end
        end
        return
    end
end