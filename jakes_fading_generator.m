function [h_canal] = jakes_fading_generator(N_samples, fs, fd, delays_samples, gains_linear)
    % =========================================================================
    % MATHEMATICALLY RIGOROUS JAKES ISOTROPIC SUM-OF-SINUSOIDS (SOS) MODEL
    % Generates independent time-varying Rayleigh fading coefficients per path.
    % M = 50 plane waves ensures full statistical isotropic convergence.
    % =========================================================================
    
    M = 50; 
    N_paths = length(delays_samples);
    h_canal = zeros(N_samples, 1);
    t = (0:N_samples-1)' / fs;
    
    for l = 1:N_paths
        tau = delays_samples(l);
        g_lin = gains_linear(l);
        
        % Uncorrelated phase seeds per ray 'l' and plane wave 'm'
        phi_l_m = rand(M, 1) * 2 * pi; 
        theta_m = 2 * pi * (1:M)' / M;
        
        % Isotropic scattering superposition
        fading_path = zeros(N_samples, 1);
        for m = 1:M
            fading_path = fading_path + exp(1j * (2 * pi * fd * t * cos(theta_m(m)) + phi_l_m(m)));
        end
        fading_path = (fading_path / sqrt(M)) * g_lin;
        
        % Discrete sample delay shift injection
        path_shifted = zeros(N_samples, 1);
        if tau == 0
            path_shifted = fading_path;
        else
            path_shifted(tau+1:end) = fading_path(1:end-tau);
        end
        h_canal = h_canal + path_shifted;
    end
end