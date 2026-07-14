% =========================================================================
% MASTER RUNNER & HARDWARE CHARACTERIZATION EMULATOR FOR SENSORS
% Performs the complete SNR sweeps (0 to 20 dB, 0.2 dB step) with 1.5M bits.
% Pure alphabetic naming: 'QPSK' and 'QAM' to solve all field restrictions.
% =========================================================================
clearvars; clc; close all;

% 1. System Parameters (Table 1 Compliance)
modulation_modes = {'QPSK', 'QAM'}; 
fs = 30000;              % Sample frequency (30 kHz)
sps = 6;                 % Samples per symbol
BaudRate = fs / sps;     % 5000 Baud
K_val = 50000;           % 50,000 symbols per frame (100k bits QPSK)
idx_stable = round(0.35 * K_val); % 35% Pilot window (17.5k symbols)
fd = 30;                 % Doppler Shift (30 Hz)

% Multipath Profile (Table 1)
delays = [0, 2, 5];
gains_db = [0, -3, -9];
gains_lin = 10.^(gains_db / 20);

% SNR Sweep Parameters (0 to 20 dB in 0.2 dB steps)
SNR_sweep = 0:0.2:20;
n_snr = length(SNR_sweep);
iterations = 15;         % Total bits evaluated per SNR = 100,000 * 15 = 1.5M bits

% RRC Pulseshaping
rolloff = 0.35;
span = 6;
rrcFilter = rcosdesign(rolloff, span, sps, 'sqrt');

% Initialize metrics arrays (Fixed Structure Names)
metrics = struct();
for m_idx = 1:2
    mod_type = modulation_modes{m_idx};
    metrics.(mod_type).BER_raw = zeros(n_snr, 1);
    metrics.(mod_type).BER_reg = zeros(n_snr, 1);
    metrics.(mod_type).BER_unreg = zeros(n_snr, 1); 
    metrics.(mod_type).EVM_stable = zeros(n_snr, 1);
end

% Diagnostic captures at 18 dB SNR
diag_QPSK = struct();
diag_QAM  = struct();

%% 2. Execution of the Monte Carlo sweeps
for m_idx = 1:2
    mod_type = modulation_modes{m_idx};
    fprintf('\n---> Running Characterization for %s modulation...\n', mod_type);
    
    if strcmp(mod_type, 'QPSK')
        bits_per_sym = 2;
    else % QAM Mode (16-QAM mapping structure)
        bits_per_sym = 4;
    end
    
    for s = 1:n_snr
        snr_db = SNR_sweep(s);
        snr_lin = 10^(snr_db / 10);
        
        err_raw = 0; err_reg = 0; err_unreg = 0; evm_val = 0;
        
        for it = 1:iterations
            % TX Baseband Generation
            bits_tx = randi([0 1], K_val * bits_per_sym, 1);
            if strcmp(mod_type, 'QPSK')
                bits_mat = reshape(bits_tx, 2, [])';
                syms_tx = ((2*bits_mat(:,1)-1) + 1j*(2*bits_mat(:,2)-1)) / sqrt(2);
            else % QAM
                bits_mat = reshape(bits_tx, 4, [])';
                syms_tx = zeros(K_val, 1);
                for sym_i = 1:K_val
                    val_real = (2 * bits_mat(sym_i, 1) + bits_mat(sym_i, 2) - 1.5) * 2;
                    val_imag = (2 * bits_mat(sym_i, 3) + bits_mat(sym_i, 4) - 1.5) * 2;
                    syms_tx(sym_i) = (val_real + 1j*val_imag) / sqrt(10);
                end
            end
            
            upsampled = zeros(K_val * sps, 1);
            upsampled(1:sps:end) = syms_tx;
            tx_baseband = conv(upsampled, rrcFilter, 'same');
            
            % Rayleigh Multipath Time-Varying Jakes Channel
            tx_channel = jakes_fading_generator(length(tx_baseband), fs, fd, delays, gains_lin);
            tx_faded = tx_baseband .* tx_channel;
            
            % Controlled SNR Noise Addition
            Ps = mean(abs(tx_faded).^2);
            noise = sqrt(Ps / (2 * snr_lin)) * (randn(size(tx_faded)) + 1j*randn(size(tx_faded)));
            rx_raw = tx_faded + noise;
            
            % Symmetrical RRC Matched Filter at Receiver
            rx_filt = conv(rx_raw, rrcFilter, 'same');
            rx_symbols = rx_filt(1:sps:end); 
            rx_symbols = rx_symbols(1:K_val);
            
            % --- Baseband Direct Linear Demodulation (No Eq) ---
            if strcmp(mod_type, 'QPSK')
                bits_rx_raw = [real(rx_symbols) > 0, imag(rx_symbols) > 0];
                bits_rx_raw = reshape(bits_rx_raw', [], 1);
            else % QAM
                bits_rx_raw = zeros(K_val, 4);
                for sym_i = 1:K_val
                    bits_rx_raw(sym_i, 1) = real(rx_symbols(sym_i)) > 0;
                    bits_rx_raw(sym_i, 2) = abs(real(rx_symbols(sym_i)) * sqrt(10)) > 2;
                    bits_rx_raw(sym_i, 3) = imag(rx_symbols(sym_i)) > 0;
                    bits_rx_raw(sym_i, 4) = abs(imag(rx_symbols(sym_i)) * sqrt(10)) > 2;
                end
                bits_rx_raw = reshape(bits_rx_raw', [], 1);
            end
            err_raw = err_raw + sum(bits_tx ~= bits_rx_raw);
            
            % --- Regularized NLMS Adaptive Filtering ---
            [s_eq_reg, w_hist_reg, err_hist_reg] = nlms_equalizer(rx_symbols, syms_tx, 21, 0.08, 1e-5, idx_stable, mod_type);
            
            % --- Unregularized NLMS (Ablation Case: epsilon = 0) ---
            [s_eq_unreg, ~, ~] = nlms_equalizer(rx_symbols, syms_tx, 21, 0.08, 0, idx_stable, mod_type);
            
            % Evaluate steady-state performance (excl. pilot window)
            s_reg_ss = s_eq_reg(idx_stable:end);
            s_unreg_ss = s_eq_unreg(idx_stable:end);
            s_ideal_ss = syms_tx(idx_stable:end);
            
            % Global gain/phase alignment to avoid constellation leakage
            s_reg_ss = s_reg_ss / (std(s_reg_ss) + eps);
            s_unreg_ss = s_unreg_ss / (std(s_unreg_ss) + eps);
            
            % Bit Error Counts
            if strcmp(mod_type, 'QPSK')
                bits_rx_reg = [real(s_reg_ss) > 0, imag(s_reg_ss) > 0];
                bits_rx_reg = reshape(bits_rx_reg', [], 1);
                
                bits_rx_unreg = [real(s_unreg_ss) > 0, imag(s_unreg_ss) > 0];
                bits_rx_unreg = reshape(bits_rx_unreg', [], 1);
            else % QAM
                bits_rx_reg = zeros(length(s_reg_ss), 4);
                bits_rx_unreg = zeros(length(s_unreg_ss), 4);
                for sym_i = 1:length(s_reg_ss)
                    bits_rx_reg(sym_i, 1) = real(s_reg_ss(sym_i)) > 0;
                    bits_rx_reg(sym_i, 2) = abs(real(s_reg_ss(sym_i)) * sqrt(10)) > 2;
                    bits_rx_reg(sym_i, 3) = imag(s_reg_ss(sym_i)) > 0;
                    bits_rx_reg(sym_i, 4) = abs(imag(s_reg_ss(sym_i)) * sqrt(10)) > 2;
                    
                    bits_rx_unreg(sym_i, 1) = real(s_unreg_ss(sym_i)) > 0;
                    bits_rx_unreg(sym_i, 2) = abs(real(s_unreg_ss(sym_i)) * sqrt(10)) > 2;
                    bits_rx_unreg(sym_i, 3) = imag(s_unreg_ss(sym_i)) > 0;
                    bits_rx_unreg(sym_i, 4) = abs(imag(s_unreg_ss(sym_i)) * sqrt(10)) > 2;
                end
                bits_rx_reg = reshape(bits_rx_reg', [], 1);
                bits_rx_unreg = reshape(bits_rx_unreg', [], 1);
            end
            
            bits_ref_ss = bits_tx((idx_stable-1)*bits_per_sym + 1 : end);
            err_reg = err_reg + sum(bits_ref_ss ~= bits_rx_reg);
            err_unreg = err_unreg + sum(bits_ref_ss ~= bits_rx_unreg);
            
            % Constellation Dispersal metrics
            evm_val = evm_val + sqrt(mean(abs(s_reg_ss - s_ideal_ss).^2)) * 100;
        end
        
        metrics.(mod_type).BER_raw(s) = err_raw / (K_val * bits_per_sym * iterations);
        metrics.(mod_type).BER_reg(s) = err_reg / (length(bits_ref_ss) * iterations);
        metrics.(mod_type).BER_unreg(s) = err_unreg / (length(bits_ref_ss) * iterations);
        metrics.(mod_type).EVM_stable(s) = evm_val / iterations;
    end
    
    % Store temporal diagnostic captures at 18 dB
    if strcmp(mod_type, 'QPSK')
        diag_QPSK.s_eq_ss = s_reg_ss(1:2000); 
        diag_QPSK.w_hist = w_hist_reg;
        diag_QPSK.err_hist = err_hist_reg;
    else % QAM
        diag_QAM.s_eq_ss = s_reg_ss(1:2000);
    end
end

save('synchronized_HIL_dataset.mat', 'metrics', 'diag_QPSK', 'diag_QAM', 'SNR_sweep', 'K_val', 'idx_stable');
fprintf('\n>>> Synchronized Dataset successfully generated and locked!\n');