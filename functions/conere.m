function lcore = conere(npoin,nelem,intmat,lcore)

    for ip = 1:npoin % execution loop label 1000
        lcore(ip) = 0;
    end

    for ie = 1:nelem
        for in = 1:3
            ip = intmat(in,ie);
            lcore(ip) = lcore(ip) + 1;
        end
    end

    return
end