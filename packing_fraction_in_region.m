function pf = packing_fraction_in_region(x, y, radii, container_width, region_bottom, region_top)
% packing_fraction_in_region computes the local packing fraction in a horizontal strip
% defined from region_bottom to region_top, accounting for partial overlaps.
%
% Inputs:
%   x, y            - particle positions
%   radii           - vector of particle radii
%   container_width - width of container (assumes full horizontal span)
%   region_bottom   - lower boundary of the region
%   region_top      - upper boundary of the region
%
% Output:
%   pf              - packing fraction in the defined region

    num_particles = length(x);
    area_in_region = 0;
    region_height = region_top - region_bottom;

    for i = 1:num_particles
        r = radii(i);
        yc = y(i);  % center y-position
        y_top = yc + r;
        y_bot = yc - r;

        % Case 1: completely inside the region
        if y_bot >= region_bottom && y_top <= region_top
            area_in_region = area_in_region + pi * r^2;

        % Case 2: fully above or below
        elseif y_top <= region_bottom || y_bot >= region_top
            continue;

        % Case 3: partial overlap
        else
            % Compute overlap with bottom boundary
            if y_bot < region_bottom && y_top > region_bottom
                h = region_bottom - yc;
                h = max(min(h, r), -r);  % clamp
                theta = acos(-h / r);
                area_bottom = r^2 * theta + h * sqrt(r^2 - h^2);
            else
                area_bottom = 0;
            end

            % Compute overlap with top boundary
            if y_bot < region_top && y_top > region_top
                h = region_top - yc;
                h = max(min(h, r), -r);
                theta = acos(h / r);
                area_top = r^2 * theta - h * sqrt(r^2 - h^2);
            else
                area_top = 0;
            end

            % Total area inside region = full - (top + bottom outside)
            area_inside = pi * r^2 - area_top - area_bottom;
            area_in_region = area_in_region + max(area_inside, 0);
        end
    end

    region_area = container_width * region_height;
    pf = area_in_region / region_area;
end
