function coord = smooth(ndimn,nnode,nelem,npoin,nsmoo, lcoor,lcore,intmat,coord,coor0)
    ndivn = 30;
    x = zeros(6,2);
    for idivn = 1:ndivn % execution loop label 500
        rdivn(idivn) = 1/idivn;
    end
    % smooth out the grid in nsmoo steps
    if nsmoo == 0 % jump gate label 10001
    else
        for ismoo = 1:nsmoo % execution loop label 10000
            % set coor0 = 0
            for ip = 1:npoin % execution loop label 1200
                for id = 1:2 % execution loop label 1201
                    coor0(id,ip) = 0.0;
                end
            end
            % loop over the elements
            for ielem = 1:nelem % execution loop label 2000
                for ic = 1:nnode % execution loop label 2100
                    in = intmat(ic,ielem);
                    for id = 1:ndimn % execution loop label 2101
                        x(ic,id) = coord(id,in);
                        x(nnode+ic,id) = coord(id,in);
                    end
                    node(ic) = in;
                end
                
                for ic = 1:nnode % exection loop label 2200
                    in = node(ic);
                    if(lcoor(in) ~= 0) % jump gate label 2199
                        continue
                    else
                        for jc = 1:nnode-1 % execution loop 2300
                            for id = 1:ndimn % execution loop 2301
                                coor0(id,in) = coor0(id,in) + x(ic+jc,id);
                            end
                        end
                    end
                    % jump gate exit label 2199
                end
            end
            for icoor = 1:npoin % execution loop 3000
                if(lcoor(icoor) ~= 0) % jump gate label 3000
                    continue
                else
                    is = 2*lcore(icoor);
                    cn = rdivn(is);
                    for id = 1:ndimn % execution loop 3100
                        coord(id,icoor) = cn*coor0(id,icoor);
                    end
                end
                % jump gate exit label 300
            end
        end
    end
    % jump gate exit label 10001
    
    return
end