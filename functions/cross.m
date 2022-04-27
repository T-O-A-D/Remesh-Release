function [ind, t] = cross(x0,xc,x1,xp,xq,tmin,ind,t)
    
    % this subroutine finds out the intersection between a quadratic 
    % and straight line segment 
    ind = 0;
    ax = x0(1) - 2*xc(1) + x1(1);
    ay = x0(2) - 2*xc(2) + x1(2);
    bx = 2.*(xc(1) - x0(1));
    by = 2.*(xc(2) - x0(2));
    cx = x0(1);
    cy = x0(2);
    dx = xp(1);
    dy = xp(2);
    ex = xq(1) - xp(1);
    ey = xq(2) - xp(2);
    a = ax*ey - ay*ex;
    b = bx*ey - by*ex;
    c = (cx - dx)*ey - (cy-dy)*ex;

    if (abs(a) > 1*10^-3)
        % the discriminant
        dis = b^2 - 4.*a*c;
        dis = sqrt(dis);
        i1 = 1;
        i2 = 1;
        root1 = (-b+dis)/(2.*a);
        if (root1 > 1 || root1 <= tmin)
            i1 = 0;
        end
        root2 = (-b-dis)/(2.*a);
        if (roo2 > 1 || root2 <= tmin)
            i2 = 0;
        end
        if (abs(ex) > 1*10^-4)
            s1 = (ax*root1*root1 + bx*root1+cx-dx)/ex;
            s2 = (ax*root2*root2 + bx*root2+cx-dx)/ex;
        else
            s1 = (ay*root1*root1 + by*root1+cy-dy)/ey;
            s2 = (ay*root2*root2 + by*root2+cy-dy)/ey;
        end
        
        if (s1 <0 || s1 > 1)
            i1 = 0;
        end
        if (s2 < 0 || s2 > 1)
            i2 = 0;
        end

        if i1 == 0 && i2 == 0 
            return
        end
        ind = 1;
        if i1 == 0
            root1 = 1*10^6;
        end
        if i2 == 0
            root2 = 1*10^6;
        end
        t = amin1([root1,root2]);
    else
        % only ont root at most
        if (abs(b) < 1*10^-5) 
            if (abs(c) > 1*10^-8)
                return
            else
                fprintf('Trouble in cross!!');
            end
        else
            i1 = 1;
            root1 = -c/b;
            if root1 > 1 || root1 < tmin
                i1 = 0;
            end
            if abs(ex) > 1*10^-3 
                s1 = (bx*root1 + cx -dx)/ex;
            else
                s1 = (by*root1 + cy -dy)/ey;
            end
            if s1 < 0 || s1 > 1
                i1 = 0;
            end
            if i1 == 0
                return
            end
            ind = 1;
            t = root1;
        end
    end
    
    return
end
