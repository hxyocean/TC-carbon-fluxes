function [dist] = fun_TC_dist(tem_lat, tem_lon, grid_lat, grid_lon)

    tem_cos = sind(tem_lat).*sind(grid_lat)+...
                    cosd(tem_lat).*cosd(grid_lat).*cosd(grid_lon - tem_lon);
    dist = 6371 .* acos(tem_cos);

end
