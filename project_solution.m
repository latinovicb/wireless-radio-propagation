% Constants
frequencies = [900e6]; % UMTS, LTE, 5G frequencies
ht = 40; % Transmitter antenna height in meters
hr = 1,8;  % Receiver antenna height in meters
epsilon_r = 17.2; % Relative permittivity of the ground
c = 3e8

% Distance range


% Plotting for each frequency
figure;
for f = frequencies
    lambda = c / f; % Wavelength
    received_power = zeros(size(distances));
    
    Ev_vals = {};
    Eh_vals = {};

    distances = {}; % points from 100m to 10km

    k = (2*pi)/(c/f);
    for d = 100:10:1000
        distances  = [distances, d];

        beta_angle = atan((ht + hr) / d);
        R1 = (ht/sin(beta_angle));
        R2 = (hr/sin(beta_angle));
        Rd = sqrt((ht - hr)^2 + d^2);
        Ri = R1 + R2;
        delta = Ri - Rd;

        Rv = (-(epsilon_r*sin(beta_angle)) + sqrt(epsilon_r - cos(beta_angle).^2)) / ((epsilon_r*sin(beta_angle)) + sqrt(epsilon_r - cos(beta_angle).^2));
        Rh = (sin(beta_angle)) - sqrt(epsilon_r - cos(beta_angle).^2) / ((sin(beta_angle)) + sqrt(epsilon_r - cos(beta_angle).^2));

        E0 = (exp(-j*k*Rd)) / Rd;
        Ev = E0 * (1 + Rv*exp(-j*k*delta));
        Eh = E0 * (1 + Rh*exp(-j*k*delta));

        Ev_vals = [Ev_vals, Ev];
        Eh_vals = [Eh_vals, Eh];
        
        % theta = atan((ht - hr) / d);
        % Rh = ((cos(theta) - sqrt(epsilon_r - sin(theta)^2)) / (cos(theta) + sqrt(epsilon_r - sin(theta)^2)))^2;
        % received_power(i) = 1 / (d^2) * (1 + Rh * exp(-j * 2 * pi * d / lambda));
    end

    % Plotting
    % loglog(distances, abs(received_power), 'DisplayName', sprintf('%d MHz', f / 1e6));
    % hold on;
    plot(Ev_vals,distances);
    plot(Eh_vals, distances);

end

% title('Two-Ray Reflection Model');
% xlabel('Distance (m)');
% ylabel('Received Power');
% legend('show');
% grid on;
% hold off;

