function [u_wind,v_wind, rr] = fun_TC_wind(grid_lon, grid_lat, lon, lat, vmax)
    
    grid_lon=mod(grid_lon,360);
    lon=mod(lon,360);
    
    rmw = 46.4*exp(-0.0155*vmax+0.0169*abs(lat));
    X1 = 317.1-2.026*vmax+1.915*abs(lat);
    X2 = 25;
    n = 0.4067 + 0.0144 * vmax - 0.0038*abs(lat);
    A = 0.0696 + 0.0049 * vmax - 0.0064*abs(lat);
    
    if A<0
        A=0;
    end
    x=sym("x");
    eq=126/25^5*((rmw-x)^5)-420/25^6*((rmw-x)^6)+540/25^7*((rmw-x)^7)-...
            315/25^8*((rmw-x)^8)+70/25^9*((rmw-x)^9) - n*((1-A)*X1+A*X2)/(n*((1-A)*X1+A*X2)+rmw);
    R1=double(vpasolve(eq,x));
    R1=R1(imag(R1)==0 & R1<rmw);
    if R1<0
        R1=0;
    end
    R2=R1+25;
    
    u_wind = zeros(size(grid_lon));
    v_wind = zeros(size(grid_lon));
    rr     = zeros(size(grid_lon));
    %=====================================================================%
    %  build the wind field
    %=====================================================================%
    for j = 1:size(grid_lon,1)
        for k = 1:size(grid_lon,2)
            xp=grid_lon(j,k)-lon;
            yp=grid_lat(j,k)-lat;
            if xp>180
                xp=xp-360;
            elseif xp<-180
                xp=xp+360;
            end
            tem_cos = sind(lat).*sind(grid_lat(j,k))+...
                cosd(lat).*cosd(grid_lat(j,k)).*cosd(grid_lon(j,k) - lon);
            if tem_cos > 1
                tem_cos = 1;
            elseif tem_cos < -1
                tem_cos = -1;
            end
            rr(j,k) = 6371.*acos(tem_cos);  %%%%% distance to the center of TC [km]
            
            if (rr(j,k) <= R1)
                vr=vmax*((rr(j,k)/rmw).^n);
            elseif (rr(j,k) >= R2)
                vr=vmax*((1-A)*exp(-(rr(j,k)-rmw)/X1)+A*exp(-(rr(j,k)-rmw)/X2));
            else
                tem_coef=(rr(j,k)-R1)/25;
                tem_w=126*(tem_coef.^5)-420*(tem_coef.^6)+540*(tem_coef.^7)-315*(tem_coef.^8)+70*(tem_coef.^9);
                vr=vmax*((rr(j,k)/rmw).^n).*(1-tem_w)+...
                    vmax*((1-A)*exp(-(rr(j,k)-rmw)/X1)+A*exp(-(rr(j,k)-rmw)/X2)).*tem_w;
                
            end
            
            if lat < 0
                ang=atan2(yp,xp);
                uu=vr*sin(ang);
                vv=-vr*cos(ang);
            else
                ang=atan2(yp,xp);
                uu=-vr*sin(ang);
                vv=vr*cos(ang);
            end
           
            u_wind(j, k) = uu;
            v_wind(j, k) = vv;
            
        end
    end
end
