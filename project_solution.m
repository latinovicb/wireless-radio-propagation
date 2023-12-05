% Constants
frequencies = [900e6, 1800e6, 3600e6]; % UMTS, LTE, 5G frequencies
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
    
    counter = counter + 1;
    Ev_vals = [];
    Eh_vals = [];

    Lray_vals = [];
    Lf_vals = [];
    L_vals = [];

    distances = []; % points from 100m to 10km

    k = (2*pi)/(lambda);
    for d = 100:1:10000
        % -------------------------------- propagation models

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
        
        % -------------------------------- loss models

        % plane earth loss
        Lray = 40*log10(d) - 20*log10(ht) - 20*log10(hr); % what is r ? 

        % free space loss
        Lf = 32.4 + 20*log10(d*10e-3) + 20*log10(f); % check if Rkm is correct

        % okomura-hata
        if (f >= 150e6) && (f <= 1500e6)
          A = 69.65;
          B = 26.16;
        end

        if (f >= 1500e6) && (f <= 3601e6) % adjusted for the sake of hata model
          A = 46.3;
          B = 33.9;
        end
        hr_small_city = (1.1*log10(f) - 0.7)*hr - (1.56*log10(f) - 0.8); % small city
        L = A + B*log10(f) - 13.82*log10(ht) - hr_small_city + (C_tunable - 6.55*log10(ht))*log10(d*10-3) + Cm;

        Lray_vals = [Lray_vals, -20*log10(Lray)];
        Lf_vals = [Lf_vals, -20*log10(Lf)];
        L_vals = [L_vals, -20*log10(L)];
    end
    
    break_point_d = (4*pi*hr*hr)/(lambda);

    figure(1);
    subplot(3, 1, counter);
    semilogx(distances, Eh_vals)
    hold on;
    semilogx(distances, Ev_vals);
    hold off;
    title(model_names(counter));
    xlabel('Distance (m)');
    ylabel('Received Power [dB]');
    xline(break_point_d,"r--")
    legend('Eh','Ev','Break_point');
    grid on;
  
    figure(2);
    subplot(3, 1, counter);
    semilogx(distances, Lray_vals);
    hold on;
    semilogx(distances, Lf_vals);
    semilogx(distances, L_vals);
    hold off;
    title(model_names(counter));
    xlabel('Distance (m)');
    ylabel('Path loss [-dB]');
    legend('plane earth loss', 'free space loss','cost 231 loss');
    grid on;
end

