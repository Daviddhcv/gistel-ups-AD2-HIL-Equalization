% =========================================================================
% AUTOMATED HIGH-RESOLUTION FIGURE GENERATOR FOR SENSORS MANUSCRIPT
% Compiles Figures 1 to 7 in complete professional English language.
% =========================================================================
clearvars; clc; close all;

if ~exist('synchronized_HIL_dataset.mat', 'file')
    error('Please run Run_HIL_Simulation.m first to generate the validated datasets.');
end
load('synchronized_HIL_dataset.mat');

set(0, 'DefaultAxesFontName', 'Helvetica');
set(0, 'DefaultTextFontName', 'Helvetica');

% -------------------------------------------------------------------------
% FIGURE 1: MSE Convergence Profile and Pilot Marking
% -------------------------------------------------------------------------
fig1 = figure('Color', 'w', 'Position', [100, 100, 750, 480]);
mse_db = 10 * log10(diag_QPSK.err_hist.^2 + eps);
mse_smoothed = movmean(mse_db, 100);
plot(1:K_val, mse_smoothed, 'LineWidth', 2, 'Color', [0.0, 0.27, 0.53]); hold on;

% Shaded 35% Pilot Window Box
patch([0, idx_stable, idx_stable, 0], [-35, -35, 5, 5], [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
grid on; box on; xlim([0, K_val]); ylim([-35, 5]);
xlabel('Symbol Index'); ylabel('Mean Squared Error (dB)');
title('Filter Convergence History (SNR = 18 dB)', 'FontWeight', 'bold');

text(200, -30, 'Pilot Training Phase', 'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.3, 0.3, 0.3]);
text(idx_stable + 100, -30, 'Decision-Directed (Autonomous)', 'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.0, 0.27, 0.53]);
print(fig1, '1_MSE_Learning', '-dpng', '-r300');

% -------------------------------------------------------------------------
% FIGURE 2: Corrected Bounded Tap Trajectories
% -------------------------------------------------------------------------
fig2 = figure('Color', 'w', 'Position', [150, 150, 750, 480]);
plot(1:K_val, diag_QPSK.w_hist', 'LineWidth', 1.5); grid on; box on;
xlim([0, K_val]); ylim([0, 1.4]); 
xlabel('Symbol Index'); ylabel('Tap Magnitude |\omega_n|');
title('FIR Adaptive Tap Trajectories under Isotropic Doppler', 'FontWeight', 'bold');
legend(arrayfun(@(x) sprintf('\\omega_{%d}', x), 0:20, 'UniformOutput', false), 'NumColumns', 7, 'Location', 'southoutside');
print(fig2, '2_Tap_Trajectories', '-dpng', '-r300');

% -------------------------------------------------------------------------
% FIGURE 3: QPSK Performance with explicit legends
% -------------------------------------------------------------------------
fig3 = figure('Color', 'w', 'Position', [200, 200, 750, 480]);
yyaxis left;
semilogy(SNR_sweep, metrics.QPSK.BER_raw, 'r--', 'LineWidth', 1.8); hold on;
semilogy(SNR_sweep, metrics.QPSK.BER_reg, 'b-', 'LineWidth', 2.2);
ylabel('Bit Error Rate (BER)'); ylim([1e-5, 1.0]);
yyaxis right;
plot(SNR_sweep, metrics.QPSK.EVM_stable, 'Color', [0.93, 0.49, 0.19], 'LineStyle', '-.', 'LineWidth', 1.8);
ylabel('Steady-State EVM (%)'); ylim([0, 50]);
grid on; box on; xlabel('Signal-to-Noise Ratio (SNR dB)');
title('QPSK Macro Performance Curves', 'FontWeight', 'bold');
legend({'Unequalized Link (BER)', 'Regularized NLMS (BER)', 'Stable EVM (%)'}, 'Location', 'southwest');
print(fig3, '3_BER_EVM_Curves', '-dpng', '-r300');

% -------------------------------------------------------------------------
% FIGURE 4: Recovered QPSK Constellation
% -------------------------------------------------------------------------
fig4 = figure('Color', 'w', 'Position', [250, 250, 500, 500]);
plot(real(diag_QPSK.s_eq_ss), imag(diag_QPSK.s_eq_ss), 'r.', 'MarkerSize', 3); hold on;
plot([-1, 1]/sqrt(2), [-1, 1]/sqrt(2), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
grid on; box on; axis equal; xlim([-1.5, 1.5]); ylim([-1.5, 1.5]);
xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');
title('QPSK Equalized Constellation (18 dB)', 'FontWeight', 'bold');
legend({'Recovered Grid', 'Target Bounds'}, 'Location', 'northeast');
print(fig4, '4_QPSK_Constellation', '-dpng', '-r300');

% -------------------------------------------------------------------------
% FIGURE 5: 16-QAM Universality (using pure alphabetic field QAM)
% -------------------------------------------------------------------------
fig5 = figure('Color', 'w', 'Position', [300, 300, 750, 480]);
yyaxis left;
semilogy(SNR_sweep, metrics.QAM.BER_raw, 'm--', 'LineWidth', 1.8); hold on;
semilogy(SNR_sweep, metrics.QAM.BER_reg, 'g-', 'LineWidth', 2.2);
ylabel('Bit Error Rate (BER)'); ylim([1e-5, 1.0]);
yyaxis right;
plot(SNR_sweep, metrics.QAM.EVM_stable, 'Color', [0.85, 0.65, 0.12], 'LineStyle', '-.', 'LineWidth', 1.8);
ylabel('Steady-State EVM (%)'); ylim([0, 60]);
grid on; box on; xlabel('Signal-to-Noise Ratio (SNR dB)');
title('16-QAM Universality Validation Performance', 'FontWeight', 'bold');
legend({'Raw Multi-Level (BER)', 'Regularized NLMS (BER)', 'Stable EVM (%)'}, 'Location', 'southwest');
print(fig5, '5_16QAM_BER_EVM', '-dpng', '-r300');

% -------------------------------------------------------------------------
% FIGURE 6: Recovered 16-QAM Grid
% -------------------------------------------------------------------------
fig6 = figure('Color', 'w', 'Position', [350, 350, 500, 500]);
plot(real(diag_QAM.s_eq_ss), imag(diag_QAM.s_eq_ss), 'g.', 'MarkerSize', 3); hold on;
[X, Y] = meshgrid(-3:2:3, -3:2:3);
grid_ideal = (X + 1j*Y) / sqrt(10);
plot(real(grid_ideal), imag(grid_ideal), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
grid on; box on; axis equal; xlim([-1.5, 1.5]); ylim([-1.5, 1.5]);
xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');
title('Recovered 16-QAM Constellation (18 dB)', 'FontWeight', 'bold');
legend({'Equalized Grid', 'Target Coordinates'}, 'Location', 'northeast');
print(fig6, '6_16QAM_Constellation', '-dpng', '-r300');

% -------------------------------------------------------------------------
% FIGURE 7: Final Consolidated Summary Panel Mosaic
% -------------------------------------------------------------------------
fig7 = figure('Color', 'w', 'Position', [50, 50, 1100, 850]);

subplot(3,2,1);
plot(1:K_val, mse_smoothed, 'LineWidth', 1.5, 'Color', [0.0, 0.27, 0.53]); hold on;
patch([0, idx_stable, idx_stable, 0], [-35, -35, 5, 5], [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
grid on; box on; xlim([0, K_val]); ylim([-35, 5]); title('(a) MSE Learning Curve');

subplot(3,2,2);
plot(1:K_val, diag_QPSK.w_hist', 'LineWidth', 1.1); grid on; box on;
xlim([0, K_val]); ylim([0, 1.4]); title('(b) FIR Filter Tap History');

subplot(3,2,3);
yyaxis left; semilogy(SNR_sweep, metrics.QPSK.BER_raw, 'r--'); hold on;
semilogy(SNR_sweep, metrics.QPSK.BER_reg, 'b-', 'LineWidth', 1.5); ylim([1e-5, 1]);
yyaxis right; plot(SNR_sweep, metrics.QPSK.EVM_stable, 'Color', [0.93, 0.49, 0.19], 'LineStyle', '-.'); ylim([0, 50]);
grid on; box on; xlim([0, 20]); title('(c) QPSK Macro Curves');

subplot(3,2,4);
plot(real(diag_QPSK.s_eq_ss), imag(diag_QPSK.s_eq_ss), 'r.', 'MarkerSize', 2); hold on;
plot([-1, 1]/sqrt(2), [-1, 1]/sqrt(2), 'ks', 'MarkerFaceColor', 'k');
grid on; box on; axis equal; xlim([-1.5, 1.5]); ylim([-1.5, 1.5]); title('(d) Recovered QPSK Grid');

subplot(3,2,5);
yyaxis left; semilogy(SNR_sweep, metrics.QAM.BER_raw, 'm--'); hold on;
semilogy(SNR_sweep, metrics.QAM.BER_reg, 'g-', 'LineWidth', 1.5); ylim([1e-5, 1]);
yyaxis right; plot(SNR_sweep, metrics.QAM.EVM_stable, 'Color', [0.85, 0.65, 0.12], 'LineStyle', '-.'); ylim([0, 60]);
grid on; box on; xlim([0, 20]); title('(e) 16-QAM Macro Curves');

subplot(3,2,6);
plot(real(diag_QAM.s_eq_ss), imag(diag_QAM.s_eq_ss), 'g.', 'MarkerSize', 2); hold on;
plot(real(grid_ideal), imag(grid_ideal), 'ks', 'MarkerFaceColor', 'k');
grid on; box on; axis equal; xlim([-1.5, 1.5]); ylim([-1.5, 1.5]); title('(f) Recovered 16-QAM Grid');

sgtitle('GISTEL-UPS Hardware-in-the-Loop Validation Characterization Panel', 'FontSize', 12, 'FontWeight', 'bold');
print(fig7, '7_Panel_Resumen_Sensors', '-dpng', '-r300');
close all;