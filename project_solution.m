% Constants
frequencies = [900e6, 1500e6, 3600e6]; % UMTS, LTE, 5G frequencies
ht = 40; % Transmitter antenna height in meters
hr = 1,8;  % Receiver antenna height in meters
epsilon_r = 17.2; % Relative permittivity of the ground

% Distance range
distances = logspace(100, 1, 100); % 100 points from 100m to 10km

k = (2*pi)/(c/f)
beta_angle = atan((ht + hr) / d);

R1 = (ht/sin(beta_angle))
R2 = (hr/sin(beta_angle))
Rd = sqrt((ht - hr)^2 + d^2)
Ri = R1 + R2

detla = Ri - Rd

% Plotting for each frequency
figure;
for f = frequencies
    lambda = 3e8 / f; % Wavelength
    received_power = zeros(size(distances));

    for i = 1:length(distances)
        d = distances(i);
        theta = atan((ht - hr) / d);
        Rh = ((cos(theta) - sqrt(epsilon_r - sin(theta)^2)) / (cos(theta) + sqrt(epsilon_r - sin(theta)^2)))^2;
        received_power(i) = 1 / (d^2) * (1 + Rh * exp(-j * 2 * pi * d / lambda));
    end

    % Plotting
    loglog(distances, abs(received_power), 'DisplayName', sprintf('%d MHz', f / 1e6));
    hold on;
end

title('Two-Ray Reflection Model');
xlabel('Distance (m)');
ylabel('Received Power');
legend('show');
grid on;
hold off;

