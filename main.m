clc;
close all;
clear all;

% Constants
frequencies = [900 * 1e6, 1800 * 1e6, 3600 * 1e6]; % UMTS, LTE, 5G frequencies
model_names = {'900 Mhz', '1500 Mhz', '3.6 Ghz'};
loss_model_names = {'plane earth loss', 'free space path loss', 'COST-231 hata'};
ht = 40; % Transmitter antenna height in meters
hr = 1.8;  % Receiver antenna height in meters
epsilon_r = 17.2; % Relative permittivity of the ground
c = 3e8;

C_tunable = 45; % tunable paramter

Cm = 15; %dB

% For plotting
counter = 0;

% Plotting for each frequency
for f = frequencies
    lambda = c / f; % Wavelength
    
    distance_vector = 100:1:10000; % points from 100m to 10km
    total_dist = length(distance_vector); % total distance

    k = (2*pi)/(lambda);

    counter = counter + 1;
    Ev_vals = [];
    Eh_vals = [];

    Lray_vals = [];
    Lf_vals = [];
    L_vals = [];

    for o = 1:length(distance_vector)

        d = distance_vector(o); % distance point

        beta_angle = atan((ht + hr) / (d)); % 2 way reflection model - no epsilon_r ? 

        % -------------------------------- propagation models


        R1 = (ht/sin(beta_angle));
        R2 = (hr/sin(beta_angle));
        Rd = sqrt((ht - hr)^2 + d^2);
        Ri = R1 + R2;
        delta = Ri - Rd;

        Rv = ((-epsilon_r * sin(beta_angle)) + sqrt(epsilon_r - cos(beta_angle).^2)) / ...
        ((epsilon_r*sin(beta_angle)) + sqrt(epsilon_r - cos(beta_angle).^2));
        
        Rh = ((sin(beta_angle)) - sqrt(epsilon_r - cos(beta_angle).^2)) / ...
        ((sin(beta_angle)) + sqrt(epsilon_r - cos(beta_angle).^2));

        E0 = (exp(-j*k*Rd)) / Rd;

        Ev = E0 * (1 + Rv*exp(-j*k*delta));
        Eh = E0 * (1 + Rh*exp(-j*k*delta));

        Ev_vals = [Ev_vals, 20*log10(abs(Ev))]; % in db
        Eh_vals = [Eh_vals, 20*log10(abs(Eh))];
        
        % -------------------------------- loss models

        Lray = ((((4*pi*d))/lambda)^2) / ...
          (4 * sin((2*pi*ht*hr)/(lambda*d)).^2);
        Lray = 10*log10(Lray); % convert to db

        % free space loss
        Lf = 32.4 + 20*log10(d/1000) + 20*log10(f/1000000); % convert dist to km and frequency to Mhz

        % okomura-hata
        if (f >= 150e6) && (f <= 1500e6)
          A = 69.65;
          B = 26.16;
        end

        if (f >= 1500e6) && (f <= 3601e6) % adjusted for the sake of hata model
          A = 46.3;
          B = 33.9;
        end
        hr_small_city = (1.1*log10(f) - 0.7)*hr - (1.56*log10(f) - 0.8); % small city, f in Mhz ??
        L = A + B*log10(f/1000000) - 13.82*log10(ht) - hr_small_city + (C_tunable - 6.55*log10(ht))*log10(d/1000) + Cm;

        Lray_vals = [Lray_vals, -Lray];
        Lf_vals = [Lf_vals, -Lf]; % already in dB
        L_vals = [L_vals, -L]; % already in dB
    end
    
    break_point_d = (4*hr*ht)/(lambda);

    figure(1);
    subplot(3, 1, counter);
    semilogx(distance_vector, Eh_vals)
    hold on;
    semilogx(distance_vector, Ev_vals);
    hold off;
    title(model_names(counter));
    xlabel('Distance (m)');
    ylabel('Received Power [dBm]');
    xline(break_point_d,"r--")
    legend('Eh','Ev','Break_point');
    grid on;
  
    figure(2);
    subplot(3, 1, counter);
    semilogx(distance_vector, Lray_vals);
    hold on;
    semilogx(distance_vector, Lf_vals);
    semilogx(distance_vector, L_vals);
    hold off;
    title(model_names(counter));
    xlabel('Distance (m)');
    ylabel('Path loss [-dB]');
    legend('plane earth loss', 'free space loss','cost 231 loss');
    grid on;
end

