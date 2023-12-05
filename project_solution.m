% Constants
frequencies = [900e6, 1800e6, 3600e6]; % UMTS, LTE, 5G frequencies
model_names = {'UMTS', 'LTE', '5G'}
loss_model_names = {'plane earth loss', 'free space path loss', 'COST-231 hata'}
ht = 40; % Transmitter antenna height in meters
hr = 1.8;  % Receiver antenna height in meters
epsilon_r = 17.2; % Relative permittivity of the ground
c = 3e8;

Cm = 15; %dB

% For plotting
counter = 0;

% Plotting for each frequency
for f = frequencies
    lambda = c / f; % Wavelength
    
    counter = counter + 1;
    Ev_vals = [];
    Eh_vals = [];

    distances = []; % points from 100m to 10km

    k = (2*pi)/(c/f);
    for d = 100:1:10000
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

        Ev_vals = [Ev_vals, 10*log(abs(Ev))];
        Eh_vals = [Eh_vals, 10*log(abs(Eh))];
        
    end

    figure(1);
    subplot(3, 1, counter);
    semilogx(distances, Eh_vals)
    hold on;
    semilogx(distances, Ev_vals);
    title(model_names(counter));
    xlabel('Distance (m)');
    ylabel('Received Power');
    legend('show');
    grid on;

    figure(2);
    subplot(3, 1, counter);
    semilogx(distances, Eh_vals)
    hold on;
    semilogx(distances, Ev_vals);
    title(loss_model_names(counter));
    xlabel('Distance (m)');
    ylabel('Path loss');
    legend('show');
    grid on;
    hold off;
end

