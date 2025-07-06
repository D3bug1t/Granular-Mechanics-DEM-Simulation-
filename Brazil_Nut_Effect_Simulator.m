% Brazil Nut Effect in a container - MATLAB simulation

% Parameters
num_particles = 60;
container_width = 1;
container_height = 1;
dt = 0.015;
steps = 2000;
g = 9.81;
vibration_amplitude = 0.03;
vibration_freq = 7;
mode = 0; % 0 for Sinusoidal 1 for Dirac-Delta
% Radii
radii = 0.04 * ones(num_particles,1);
large_percentage = 1/num_particles;
brazil_nuts = large_percentage*num_particles;
small_nuts = (1-large_percentage)*num_particles;
radii(1:1:brazil_nuts) = 0.15;  % Brazil nut (big)

% Initialize particle positions in increasing vertical order
tap_interval = 20;      % every 20 steps → tap every (dt * 50) seconds
tap_accel = 200;         % spike acceleration (m/s^2), adjust this

x = container_width * rand(num_particles, 1);
y_sorted = linspace(0, container_height - max(radii), num_particles)';
[~, shuffle_idx] = sort(rand(num_particles,1));
x = x(shuffle_idx);
y = y_sorted(shuffle_idx);

for i = 1:brazil_nuts
    y(i) = radii(i); 
end
x(1:1:brazil_nuts) = radii(brazil_nuts)*ones(brazil_nuts,1) + (container_width-2*radii(brazil_nuts)) * 0.5*ones(brazil_nuts);
vx = zeros(num_particles,1);
vy = zeros(num_particles,1);
mass = radii.^2; % mass ∝ area
Gamma = vibration_amplitude*((2*pi*vibration_freq)^2)/g;
prev_base_y = 0;
base_y = 0;
% Visualization
figure;
for t = 1:steps
    y(1)
    clf;
    hold on;
    grid on;
    axis([0 container_width -vibration_amplitude-0.01 container_height]);
    axis equal;
    title(sprintf('Step %d', t));
    % Used BFS to get the approriate contact layers
    layer_index = compute_layered_contacts(x, y, radii,base_y);
    max_layer = max(layer_index);
    brazil_y(t) = y(1);
      colors = lines(10);
    for i = 1:num_particles
        layer = layer_index(i);
        if layer == -1
            color = [0.6 0.6 0.6];
        else
            color = colors(mod(layer, size(colors,1)) + 1, :);
        end
        viscircles([x(i), y(i)], radii(i), 'Color', color);
    end
    % Drawing container box
    rectangle('Position', [0, 0, container_width, container_height], ...
              'EdgeColor', 'k', 'LineWidth', 2);
    if(mode==0)
        vibration_offset = vibration_amplitude * cos(2 * pi * vibration_freq * t * dt);  % dynamic vertical offset
        line_y = -vibration_amplitude+ vibration_offset;  % center around -0.2
        base_y = vibration_amplitude * cos(2 * pi * vibration_freq * t * dt);  % moving floor position
        plot([0 container_width], [line_y line_y], 'r-', 'LineWidth', 3);
    end
    
    % Processing layer by layer upwards
    for layer = 0:max(max_layer,0)
        layer_particles = find(layer_index == layer);
        for k = 1:length(layer_particles)
            i = layer_particles(k);
            % Apply tap as Dirac delta/Sinusoidal to base layer only
            if(mode==0)
                a_base = - (2 * pi * vibration_freq)^2 * vibration_amplitude * cos(2 * pi * vibration_freq * t*dt);
               if y(i) - radii(i) <= base_y
                        a_total = a_base - g;  % only apply shaking if particle on floor
                    else
                        a_total = -g;
               end
            else
                if mod(t, tap_interval) == 0 && layer == 0
                    a_total = tap_accel - g;
                else
                    a_total = -g;
                end
            end
    
            % Update velocity and position
            vy(i) = vy(i) + a_total * dt;
            y(i) = y(i) + vy(i) * dt;
    
            % Floor collision
            if y(i) - radii(i) <= base_y
                % y(i) = radii(i);
                y(i) = base_y + radii(i);
                % vy(i) = -0.3 * vy(i);
                vy(i) = -0.3 * vy(i) + (base_y - prev_base_y) / dt; %prev_base set to 0 for Dirac Delta Tapping
            end
    
            % Ceiling collision
            if y(i) + radii(i) > container_height
                y(i) = container_height - radii(i);
                vy(i) = -0.3 * vy(i);
            end
    
            % Wall collision
            if x(i) - radii(i) < 0
                x(i) = radii(i);
                vx(i) = -0.3 * vx(i);
            elseif x(i) + radii(i) > container_width
                x(i) = container_width - radii(i);
                vx(i) = -0.3 * vx(i);
            end
        end
    end
    % Process disconnected floating particles (not part of contact structure)
    % (in BFS marked as -1)
    disconnected_particles = find(layer_index == -1);
    for k = 1:length(disconnected_particles)
        i = disconnected_particles(k);
    
        % Apply gravity only
        a_total = -g;
    
        % Update velocity and position
        vy(i) = vy(i) + a_total * dt;
        y(i) = y(i) + vy(i) * dt;
    
        % Collisions
        if y(i) - radii(i) < 0
            % y(i) = radii(i);
            y(i) = base_y + radii(i);
            % vy(i) = -0.3 * vy(i);
            vy(i) = -0.3 * vy(i) + (base_y - prev_base_y) / dt;
        end
        if y(i) + radii(i) > container_height
            y(i) = container_height - radii(i);
            vy(i) = -0.3 * vy(i);
        end
        if x(i) - radii(i) < 0
            x(i) = radii(i);
            vx(i) = -0.3 * vx(i);
        elseif x(i) + radii(i) > container_width
            x(i) = container_width - radii(i);
            vx(i) = -0.3 * vx(i);
        end
    end
    
    
    
        % Basic inter-particle repulsion 
    for i = 1:num_particles
        for j = i+1:num_particles
            dx = x(i) - x(j);
            dy = y(i) - y(j);
            dist = sqrt(dx^2 + dy^2);
            min_dist = radii(i) + radii(j);
    
            if dist < min_dist && dist > 0
                % Normalize collision normal
                nx = dx / dist;
                ny = dy / dist;
    
                % Relative velocity
                dvx = vx(i) - vx(j);
                dvy = vy(i) - vy(j);
    
               % Velocity along normal
                vn = dvx * nx + dvy * ny;
    
                if vn < 0  % they are moving toward each other
                    % Compute impulse scalar
                    m1 = mass(i);
                    m2 = mass(j);
                    impulse = (2 * vn) / (m1 + m2);
    
                    % Apply impulse to velocities
                    vx(i) = vx(i) - impulse * m2 * nx;
                    vy(i) = vy(i) - impulse * m2 * ny;
                    vx(j) = vx(j) + impulse * m1 * nx;
                    vy(j) = vy(j) + impulse * m1 * ny;
                end
    
                % Position correction (to remove overlap)
                overlap = min_dist - dist + 1e-6;
                x(i) = x(i) + 0.5 * overlap * nx;
                y(i) = y(i) + 0.5 * overlap * ny;
                x(j) = x(j) - 0.5 * overlap * nx;
                y(j) = y(j) - 0.5 * overlap * ny;
            end
        end       
    end
    % vibration_offset = vibration_amplitude * sin(2*pi*vibration_freq*t*dt); % scale as needed
    % plot([0 container_width], [-0.1 -0.1] + vibration_offset, 'r-', 'LineWidth', 2);
    % text(container_width/2, -0.3, sprintf('Vibration: %.2f', vibration_offset), ...
    %      'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 10);
    
      
    packing_fraction(t) = packing_fraction_in_region(x,y,radii,container_width,base_y,0.25*container_height+base_y);
    prev_base_y = base_y;
    pause(0.01);
end

%%
figure(2)
time_mesh = 1:1:steps;
time_mesh = time_mesh*dt;
plot(time_mesh,brazil_y)
xlabel("Time(s)");
ylabel("Height of Brazil Nut");
grid on;
title("Variation of Height of Brazil-Nut effect");
%%
figure(3)
time_mesh = 1:1:steps;
time_mesh = time_mesh*dt;
plot(time_mesh,packing_fraction);
xlabel("Time(s)");
ylabel("Packing Fraction {\phi}");
grid on;
title("Variation of Packing Fraction");
function layer_index = compute_layered_contacts(x, y, radii,base_y)
    num_particles = length(x);
    epsilon = 1e-4;

    % Initialize all particles as unassigned
    layer_index = -1 * ones(num_particles, 1);

    %Assigning Layer 0 (touching the floor)
    current_layer = 0;
    layer_particles = find(y - radii <= base_y+epsilon);
    layer_index(layer_particles) = current_layer;

    % Iteratively building layers
    while ~isempty(layer_particles)
        next_layer_particles = [];

        for i = layer_particles'
            for j = 1:num_particles
                if layer_index(j) == -1
                    dx = x(i) - x(j);
                    dy = y(i) - y(j);
                    dist = sqrt(dx.^2 + dy.^2);
                    if dist <= radii(i) + radii(j) + epsilon
                        layer_index(j) = current_layer + 1;
                        next_layer_particles(end+1) = j; 
                    end
                end
            end
        end

        current_layer = current_layer + 1;
        layer_particles = next_layer_particles;
    end
end
    
