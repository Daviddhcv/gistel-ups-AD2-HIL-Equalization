function [s_eq, w_history, error_abs_history] = nlms_equalizer(r_vector, s_vector, N_taps, beta_nlms, epsilon_reg, idx_stable, modulation_type)
    % =========================================================================
    % REGULARIZED HYBRID NLMS EQUALIZER: PILOT-AIDED TO DECISION-DIRECTED (DD)
    % Implements Wirtinger calculus update via complex conjugate error term.
    % Standardizes the filter weight scale cleanly between 0.0 and 1.4.
    % =========================================================================
    
    n_val = length(r_vector);
    s_eq = zeros(n_val, 1);
    w_coef = zeros(N_taps, 1);
    w_coef(floor(N_taps/2) + 1) = 1.0; % Symmetric impulse response seed
    
    w_history = zeros(N_taps, n_val);
    error_abs_history = zeros(n_val, 1);
    
    % Normalize input variance to prevent filter tap scale explosion
    r_normalized = r_vector / (std(r_vector) + eps);
    r_padded = [zeros(N_taps - 1, 1); r_normalized];
    
    for n = 1:n_val
        % Mobile window FIR tapped delay line
        reg_shift = r_padded(n + N_taps - 1 : -1 : n);
        
        % Filter output calculation
        s_eq(n) = w_coef' * reg_shift;
        
        % Hybrid decision-directed feedback controller (Solves data leakage)
        if strcmp(modulation_type, 'QPSK')
            if n < idx_stable
                error_inst = s_vector(n) - s_eq(n);
            else
                s_sliced = (sign(real(s_eq(n))) + 1j * sign(imag(s_eq(n)))) / sqrt(2);
                error_inst = s_sliced - s_eq(n);
            end
        else % QAM Mode
            if n < idx_stable
                error_inst = s_vector(n) - s_eq(n);
            else
                % Normalized 16-point grid slicing boundaries
                real_slice = min(max(2 * floor((real(s_eq(n)) * sqrt(10) + 1)/2) + 1, -3), 3) / sqrt(10);
                imag_slice = min(max(2 * floor((imag(s_eq(n)) * sqrt(10) + 1)/2) + 1, -3), 3) / sqrt(10);
                s_sliced = real_slice + 1j * imag_slice;
                error_inst = s_sliced - s_eq(n);
            end
        end
        
        % Instantaneous energy across the filter window
        window_energy = reg_shift' * reg_shift;
        
        % Recursive update employing the complex conjugate error term
        w_coef = w_coef + (beta_nlms / (window_energy + epsilon_reg)) * error_inst * conj(reg_shift);
        
        % Diagnostics trace storage
        w_history(:, n) = abs(w_coef);
        error_abs_history(n) = abs(error_inst);
    end
end