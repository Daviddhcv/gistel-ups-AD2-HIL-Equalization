%% ADVANCED METROLOGICAL SCRIPT: HIL TRANSCEIVER WITH COMPLETE ENGLISH PLOTS


%% 0. LIMPIEZA DEL ESPACIO DE TRABAJO
clearvars; % Remueve variables del espacio base para impedir contaminación cruzada de datos.
clc; % Despeja la interfaz de la consola de comandos de MATLAB.
close all; % Cierra los lienzos y ventanas gráficas secundarias activas.

try
    %% 1. PARÁMETROS ESTRUCTURALES DEL ENTORNO EXPERIMENTAL
    tipo_modulacion = 'QPSK'; % Selección del formato de modulación digital vectorial bajo análisis.
    deviceID        = "AD1";   % Identificador físico de la tarjeta Analog Discovery 2 en el sistema operativo.
    fs              = 30000;   % Frecuencia de muestreo analógica establecida rígidamente en 30 kHz.
    snr_vector      = 0:2:16;  % Vector lineal de barrido para la relación señal a ruido en decibelios.
    pruebas_por_snr = 15;      % Iteraciones para el promedio estadístico y reducción de varianza.
    
    % Parámetros estocásticos del canal variante en el tiempo (Modelo Doppler de Jakes)
    fd_Doppler   = 30;        % Frecuencia Doppler máxima a 30 Hz para emulación de condiciones severas.
    retardos_sam = [0, 2, 5]; % Vector de retardos multitrayecto discretizados en muestras temporales.
    ganancias_db = [0, -3, -9];% Atenuación de potencia en decibelios para cada rayo independiente.
    ganancias_lin= 10.^(ganancias_db / 20); % Homologación de ganancias a magnitudes de escala lineal.
    
    % Gestión automatizada del directorio de volcado metrológico Q1
    carpeta_salida = fullfile(pwd, ['Sensors_Q1_NLMS_' datestr(now, 'yyyymmdd_HHMMSS')]);
    if ~exist(carpeta_salida, 'dir')
        mkdir(carpeta_salida); % Construcción del directorio físico para preservación de trazas.
    end
    
    % Dimensionamiento extendido de la carga útil binaria para análisis de estado estable
    bits_per_trama = 2500 * 2; % Trama de 2500 símbolos para garantizar el análisis asintótico del filtro.

    %% 2. CONFIGURACIÓN DE LA INTERFAZ DE INSTRUMENTACIÓN VIRTUAL (DAQ)
    dq      = daq("digilent"); % Instanciación del controlador de adquisición de datos de Digilent.
    dq.Rate = fs;               % Asignación de la tasa de transferencia de muestras en el hardware.
    flush(dq);                  % Vaciado físico de los registros de memoria interna de la AD2.
    
    addoutput(dq, deviceID, "1", "Voltage"); % Asignación de la salida analógica W1 (Componente I).
    addoutput(dq, deviceID, "2", "Voltage"); % Asignación de la salida analógica W2 (Componente Q).
    addinput(dq, deviceID, "1", "Voltage");  % Canal de lectura CH1 para digitalización de fase I.
    addinput(dq, deviceID, "2", "Voltage");  % Canal de lectura CH2 para digitalización de cuadratura Q.

    %% 3. SÍNTESIS DEL FILTRO DE CONFORMACIÓN VECTORIAL RRC
    sps       = 6; % Muestras por símbolo (Velocidad de modulación física = 5000 Baudios).
    rolloff   = 0.5; % Factor de exceso de ancho de banda para el confinamiento espectral.
    span      = 6; % Duración del soporte temporal medido en número de símbolos.
    rrcFilter = rcosdesign(rolloff, span, sps, 'sqrt'); % Generación de la respuesta impulsional RRC.

    %% 4. PREASIGNACIÓN DE MATRICES PARA CARACTERIZACIÓN ESTADÍSTICA
    n_snr       = length(snr_vector); % Dimensión del bucle exterior de evaluación de ruido.
    ber_sin_eq  = zeros(n_snr, 1);     % Vector para la tasa de error bruta sujeta a selectividad.
    ber_con_eq  = zeros(n_snr, 1);     % Vector para la tasa de error neta corregida por el ecualizador.
    evm_final   = zeros(n_snr, 1);     % Contenedor para el error vectorial cuadrático medio estable.
    
    % Matrices de captura metrológica profunda para el último punto de SNR
    N_simbolos_estimados = (bits_per_trama / 2); % Cantidad nominal de símbolos por trama de datos.
    w_history = zeros(21, N_simbolos_estimados); % Historial de evolución temporal de pesos adaptativos.
    error_abs_history = zeros(N_simbolos_estimados, 1); % Historial del residuo escalar de error complejo.
    rx_signal_power_density = []; % Contenedor espectral para estimación de densidad de potencia.

    %% 5. BUCLE EJECUTIVO DE BARRIDO METROLÓGICO
    fprintf('\n=== EVALUACIÓN METROLÓGICA AVANZADA — ENFOQUE SENSORS Q1 ===\n');
    
    for i = 1:n_snr
        snr_db  = snr_vector(i); % Extracción de la relación señal a ruido bajo análisis.
        snr_lin = 10^(snr_db / 10); % Conversión logarítmica a parámetro escalar lineal.
        
        err_acum_sin = 0; % Inicialización del contador de colisiones binarias sin ecualizar.
        err_acum_con = 0; % Inicialización del contador de errores bajo corrección adaptativa.
        evm_acum     = 0; % Inicialización del acumulador euclidiano de dispersión de fase.
        total_bits   = 0; % Sumatoria para la normalización de la tasa de error global.
        total_bits_evm = 0; % Sumatoria exclusiva para la región de régimen permanente.

        for p = 1:pruebas_por_snr
            
            % --- TRANSMISOR DIGITAL DE SEÑALES (DSP) ---
            bits_tx = randi([0 1], bits_per_trama, 1); % Generación del vector binario pseudoaleatorio.
            bits_mat = reshape(bits_tx, 2, [])'; % Segmentación bidimensional para entrelazado ortogonal.
            syms_tx = ((2*bits_mat(:,1)-1) + 1j*(2*bits_mat(:,2)-1)) / sqrt(2); % Mapeo de la constelación QPSK.
            
            upsampled = zeros(length(syms_tx) * sps, 1); % Creación de la estructura del vector interpolado.
            upsampled(1:sps:end) = syms_tx; % Inserción de nulos muestrales entre símbolos de datos.
            tx_baseband = conv(upsampled, rrcFilter, 'same'); % Convolución discreta FIR para la conformación de pulso.

            % --- MODELADO AVANZADO DEL CANAL ESTOCÁSTICO (RAYLEIGH MULTIPATH + DOPPLER) ---
            N_muestras = length(tx_baseband); % Medición longitudinal de la trama de transmisión.
            tx_canal = zeros(N_muestras, 1); % Inicialización del contenedor para la mezcla del canal.
            t = (0:N_muestras-1)' / fs; % Eje de tiempo discreto síncrono con la tasa analógica.
            
            for path = 1:length(retardos_sam)
                tau = retardos_sam(path); % Extracción del retardo intrínseco del trayecto.
                g_lin = ganancias_lin(path); % Extracción de la magnitud lineal del trayecto.
                
                % Desvanecimiento Rayleigh variant en el tiempo basado en el espectro selectivo de Jakes
                rayleigh_gain = (randn(N_muestras, 1) + 1j*randn(N_muestras, 1)) / sqrt(2); % Proceso complejo gaussiano ortogonal.
                filtro_doppler = cos(2 * pi * fd_Doppler * t); % Componente espectral oscilatoria de la portadora.
                coef_dinamico = rayleigh_gain .* filtro_doppler * g_lin; % Coeficiente variable acoplado.
                
                % Desplazamiento lineal de muestras en el dominio del tiempo discreto
                tx_desplazada = zeros(size(tx_baseband)); % Inicialización del buffer desplazado.
                if tau == 0
                    tx_desplazada = tx_baseband; % El primer rayo directo no posee retraso en el lazo.
                else
                    tx_desplazada(tau+1:end) = tx_baseband(1:end-tau); % Desplazamiento de muestras de sub-símbolo.
                end
                tx_canal = tx_canal + tx_desplazada .* coef_dinamico; % Sumatoria lineal de trayectos de dispersión.
            end
            
            % Adición de ruido blanco gaussiano aditivo (AWGN)
            Ps_I = mean(real(tx_canal).^2); % Evaluación de la densidad de potencia rama en fase.
            noise_I = sqrt(Ps_I / snr_lin) * randn(N_muestras, 1); % Modelado del ruido de cuantización analógico I.
            tx_I = real(tx_canal) + noise_I; % Composición analógica final para el canal DAC W1.
            
            Ps_Q = mean(imag(tx_canal).^2); % Evaluación de la densidad de potencia rama en cuadratura.
            noise_Q = sqrt(Ps_Q / snr_lin) * randn(N_muestras, 1); % Modelado del ruido de cuantización analógico Q.
            tx_Q = imag(tx_canal) + noise_Q; % Composición analógica final para el canal DAC W2.

            % --- INTERFAZ DE HARDWARE EN LAZO CERRADO REAL (readwrite) ---
            tx_matriz = [tx_I, tx_Q]; % Compactación matricial síncrona de señales de instrumentación.
            data_cap = readwrite(dq, tx_matriz); % Transferencia bidireccional física a través del bus de la AD2.
            
            rx_I_raw = data_cap{:, 1}; % Digitalización de la rama I mediante el ADC CH1.
            rx_Q_raw = data_cap{:, 2}; % Digitalización de la rama Q mediante el ADC CH2.
            
            % --- RECEPTOR DIGITAL AVANZADO ---
            rx_I_filt = conv(rx_I_raw, rrcFilter, 'same'); % Filtrado acoplado RRC para contención de ruido rama I.
            rx_Q_filt = conv(rx_Q_raw, rrcFilter, 'same'); % Filtrado acoplado RRC para contención de ruido rama Q.
            rx_compleja = rx_I_filt + 1j*rx_Q_filt; % Reconstrucción analítica del espacio de fases recibido.
            
            % Algoritmo de sincronización de muestras por máxima verosimilitud temporal
            [c, lags] = xcorr(rx_I_filt, real(tx_baseband)); % Correlación de ráfagas temporales con la referencia limpia.
            [~, idx_max] = max(abs(c)); % Localización del pico energético de coincidencia de fase.
            offset = lags(idx_max); % Extracción del desfase absoluto del lazo físico.
            
            if offset > 0 && offset < length(rx_compleja)
                rx_sinc = rx_compleja(offset+1:end); % Alineación rígida del vector elminando retardo de propagación.
            else
                rx_sinc = rx_compleja; % Retención de la estructura base ante anomalías del lazo.
            end
            
            % Diezmado adaptado al centro geométrico del diagrama de ojo
            centro = floor(sps / 2) + 1; % Índice óptimo de máxima apertura espectral.
            muestras_rx = rx_sinc(centro:sps:end); % Remuestreo síncrono a tasa de símbolo original.
            
            n_val = min(length(syms_tx), length(muestras_rx)); % Normalización dimensional de vectores acoplados.
            r_vector = muestras_rx(1:n_val); % Vector de variables de canal complejas recibidas.
            s_vector = syms_tx(1:n_val); % Vector de variables ideales preestablecidas.
            
            % Decisión dura directa (Límite lineal de canal degradado)
            bits_rx_sin = [real(r_vector) > 0, imag(r_vector) > 0]; % Demapeo directo por evaluación de signo cuadrantal.
            bits_rx_sin = reshape(bits_rx_sin', [], 1); % Reestructuración vectorial para cómputo bruto de error.
            
            % --- ALGORITMO ADAPTATIVO OPTIMIZADO: FILTRADO REGULARIZADO NLMS ---
            N_taps  = 21; % Estructura de 21 coeficientes para selectividad espectral severa.
            beta_nlms = 0.08; % Factor de velocidad adaptativa normalizada.
            epsilon_reg = 1e-5; % Factor de regularización contra desvanecimientos profundos.
            w_coef  = zeros(N_taps, 1); % Inicialización del vector de pesos del ecualizador.
            w_coef(floor(N_taps/2)+1) = 1; % Forzado del impulso central simétrico.
            
            r_padded = [zeros(N_taps-1, 1); r_vector]; % Inserción de márgenes de guarda nulos.
            s_ecualizada = zeros(n_val, 1); % Inicialización del vector de salida ecualizado.
            
            for n_idx = 1:n_val
                reg_desplazamiento = r_padded(n_idx + N_taps - 1 : -1 : n_idx); % Ventana móvil FIR.
                s_ecualizada(n_idx) = w_coef' * reg_desplazamiento; % Combinación lineal adaptativa.
                
                error_instantneo = s_vector(n_idx) - s_ecualizada(n_idx); % Residuo complejo instantáneo.
                potencia_ventana = (reg_desplazamiento' * reg_desplazamiento); % Estimación energética de la ventana.
                
                % Mecanismo de actualización NLMS regularizado
                w_coef = w_coef + (beta_nlms / (potencia_ventana + epsilon_reg)) * error_instantneo * conj(reg_desplazamiento);
                
                % Conservación de trazas temporales en el último nodo de SNR
                if i == n_snr && p == pruebas_por_snr
                    w_history(:, n_idx) = abs(w_coef); % Volcado de la magnitud de los coeficientes.
                    error_abs_history(n_idx) = abs(error_instantneo); % Volcado del error absoluto.
                end
            end
            
            % --- CONTROL METROLÓGICO AVANZADO: SEGREGACIÓN DEL TRANSITORIO ---
            idx_estable = floor(n_val * 0.35) + 1; % Ventaneo estricto (Descarte del 35% inicial de la trama).
            s_estacionaria = s_ecualizada(idx_estable:end); % Segmento estable en régimen permanente.
            s_ideal_estable = s_vector(idx_estable:end); % Símbolos de referencia homólogos.
            
            % Normalización estricta de segundo orden post-filtrado
            s_estacionaria = s_estacionaria - mean(s_estacionaria); % Supresión de derivas continuas de hardware.
            s_estacionaria = s_estacionaria / (std(s_estacionaria) + eps); % Escalamiento a potencia promedio unitaria.
            
            % Decisión de símbolos bajo el efecto de la ecualización NLMS
            bits_rx_con = [real(s_estacionaria) > 0, imag(s_estacionaria) > 0]; % Demapeo cuadrantal limpio.
            bits_rx_con = reshape(bits_rx_con', [], 1); % Vectorización binaria neta.
            
            % Alineación de los vectores binarios de referencia original
            bits_ref_total = bits_tx(1:2*n_val); % Referencia para canal degradado directo.
            bits_ref_estable = bits_tx((idx_estable-1)*2 + 1 : 2*n_val); % Referencia para zona ecualizada.
            
            % Acumulación estadística agregada por subtrama
            err_acum_sin = err_acum_sin + sum(bits_ref_total ~= bits_rx_sin(1:length(bits_ref_total)));
            err_acum_con = err_acum_con + sum(bits_ref_estable ~= bits_rx_con);
            evm_acum     = evm_acum + sqrt(mean(abs(s_estacionaria - s_ideal_estable).^2));
            total_bits   = total_bits + length(bits_ref_total);
            total_bits_evm = total_bits_evm + length(bits_ref_estable);
            
            % Almacenamiento de la densidad de potencia espectral de la señal recibida
            if i == n_snr && p == pruebas_por_snr
                rx_signal_power_density = rx_compleja; % Muestra cruda para análisis por FFT.
            end
        end % Fin del ciclo de promediación estadística de tramas
        
        % Cómputo definitivo de las variables metrológicas macro del sistema
        ber_sin_eq(i) = err_acum_sin / total_bits; % Tasa de error de bit del canal degradado.
        ber_con_eq(i) = err_acum_con / total_bits_evm; % Tasa de error de bit optimizada tras ecualización NLMS.
        evm_final(i)  = (evm_acum / pruebas_por_snr) * 100; % Porcentaje de EVM neto en estado estable.
        
        fprintf('SNR: %2d dB | BER Bruta: %.4e | BER NLMS: %.4e | EVM Estable: %5.2f%%\n', ...
                snr_db, ber_sin_eq(i), ber_con_eq(i), evm_final(i));
    end

    %% 6. COMPILACIÓN GRÁFICA EDITORIAL EN INGLÉS (FONDO CLARO EXCLUSIVO)
    % Función auxiliar para forzar el guardado unificado a 300 DPI y fondo claro.
    exportar_grafica = @(h_fig, f_name) print(h_fig, fullfile(carpeta_salida, f_name), '-dpng', '-r300');

    % -----------------------------------------------------------------------
    % FIGURA 1: Tasa de Error de Bit (BER) vs SNR (En Inglés)
    % -----------------------------------------------------------------------
    fig1 = figure('Color', 'w', 'Name', 'Metrological BER Curves');
    semilogy(snr_vector, ber_sin_eq, 'r--s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r'); hold on;
    semilogy(snr_vector, ber_con_eq, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    grid on; ax_fig1 = gca; set(ax_fig1, 'Color', 'w'); xlabel('SNR (dB)'); ylabel('BER');
    title('Bit Error Rate Performance: Degraded Channel vs. NLMS Equalizer');
    legend('Degraded Channel (No Eq)', 'Corrected Signal (NLMS Filter)', 'Location', 'southwest');
    exportar_grafica(fig1, '1_BER_Metrologica');

    % -----------------------------------------------------------------------
    % FIGURA 2: Magnitud del Vector de Error (EVM %) vs SNR (En Inglés)
    % -----------------------------------------------------------------------
    fig2 = figure('Color', 'w', 'Name', 'EVM Evolution Curves');
    plot(snr_vector, evm_final, 'm-^', 'LineWidth', 1.5, 'MarkerFaceColor', 'm');
    grid on; ax_fig2 = gca; set(ax_fig2, 'Color', 'w'); xlabel('SNR (dB)'); ylabel('EVM (%)');
    title('Steady-State Error Vector Magnitude (EVM) Evolution');
    exportar_grafica(fig2, '2_EVM_Estable');

    % -----------------------------------------------------------------------
    % FIGURA 3: Diagrama de Constelación del Transceptor (En Inglés)
    % -----------------------------------------------------------------------
    fig3 = figure('Color', 'w', 'Name', 'QPSK Equalized Constellation');
    plot(real(s_estacionaria), imag(s_estacionaria), 'r.', 'MarkerSize', 4); hold on;
    plot(real(s_ideal_estable), imag(s_ideal_estable), 'k+', 'LineWidth', 2, 'MarkerSize', 12);
    grid on; ax_fig3 = gca; set(ax_fig3, 'Color', 'w'); axis([-2 2 -2 2]);
    xlabel('In-Phase Component (I)'); ylabel('Quadrature Component (Q)');
    title(['Normalized Received Constellation @ SNR = ' num2str(snr_vector(end)) ' dB']);
    legend('Hardware Samples', 'Ideal Centers', 'Location', 'northeast');
    exportar_grafica(fig3, '3_Constelacion_Q1');

    % -----------------------------------------------------------------------
    % FIGURA 4: Seguimiento Temporal de Coeficientes NLMS (En Inglés)
    % -----------------------------------------------------------------------
    fig4 = figure('Color', 'w', 'Name', 'Filter Coefficients Tracking');
    plot(w_history', 'LineWidth', 1.1); grid on;
    ax_fig4 = gca; set(ax_fig4, 'Color', 'w'); xlabel('Processed Symbol Index'); ylabel('Filter Weight Magnitude |\omega_n|');
    title('Temporal Evolution of Equalizer Coefficients under Doppler Stress');
    exportar_grafica(fig4, '4_Trayectoria_Pesos');

    % -----------------------------------------------------------------------
    % FIGURA 5: Histograma de Distribución de Potencia del Error (En Inglés)
    % -----------------------------------------------------------------------
    fig5 = figure('Color', 'w', 'Name', 'Error Residue Distribution');
    histogram(error_abs_history, 40, 'FaceColor', 'g', 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    grid on; ax_fig5 = gca; set(ax_fig5, 'Color', 'w'); xlabel('Complex Error Magnitude |e_n|'); ylabel('Absolute Frequency');
    title('Statistical Distribution of Residual Error post-NLMS Convergence');
    exportar_grafica(fig5, '5_Histograma_Error');

    % -----------------------------------------------------------------------
    % FIGURA 6: Espectro de Potencia Estimado del Canal Selectivo (En Inglés)
    % -----------------------------------------------------------------------
    fig6 = figure('Color', 'w', 'Name', 'Channel Power Spectrum');
    [psd_val, f_vec] = periodogram(rx_signal_power_density, [], 512, fs, 'centered');
    plot(f_vec/1e3, 10*log10(psd_val), 'k-', 'LineWidth', 1.2);
    grid on; ax_fig6 = gca; set(ax_fig6, 'Color', 'w'); xlabel('Frequency (kHz)'); ylabel('Power Spectral Density (dB/Hz)');
    title('Estimated Channel Spectral Response under Multipath Fading');
    exportar_grafica(fig6, '6_Espectro_Canal');

    % -----------------------------------------------------------------------
    % FIGURA 7: Panel Resumen 7-en-1 con Etiquetas Consolidadas en Inglés
    % -----------------------------------------------------------------------
    fig7 = figure('Color', 'w', 'Name', 'Sensors Global Synthesis Panel', 'Position', [80 80 1300 850]);
    
    % Panel 1: BER (En Inglés)
    subplot(3, 2, 1); semilogy(snr_vector, ber_sin_eq, 'r--s', 'LineWidth', 1.2, 'MarkerFaceColor', 'r'); hold on;
    semilogy(snr_vector, ber_con_eq, 'b-o', 'LineWidth', 1.2, 'MarkerFaceColor', 'b'); grid on;
    ax_p1 = gca; set(ax_p1, 'Color', 'w'); xlabel('SNR (dB)'); ylabel('BER'); title('BER Performance');
    
    % Panel 2: EVM (En Inglés)
    subplot(3, 2, 2); plot(snr_vector, evm_final, 'm-^', 'LineWidth', 1.2, 'MarkerFaceColor', 'm'); grid on;
    ax_p2 = gca; set(ax_p2, 'Color', 'w'); xlabel('SNR (dB)'); ylabel('EVM (%)'); title('EVM % Evolution');
    
    % Panel 3: Constelación (En Inglés)
    subplot(3, 2, 3); plot(real(s_estacionaria), imag(s_estacionaria), 'r.', 'MarkerSize', 2); hold on;
    plot(real(s_ideal_estable), imag(s_ideal_estable), 'k+', 'LineWidth', 1.5, 'MarkerSize', 8); grid on;
    ax_p3 = gca; set(ax_p3, 'Color', 'w'); axis([-2 2 -2 2]); xlabel('I'); ylabel('Q'); title('Received Constellation');
    
    % Panel 4: Coeficientes (En Inglés)
    subplot(3, 2, 4); plot(w_history', 'LineWidth', 0.9); grid on;
    ax_p4 = gca; set(ax_p4, 'Color', 'w'); xlabel('Symbols'); ylabel('|\omega_n|'); title('Coefficients Tracking');
    
    % Panel 5: Histograma (En Inglés)
    subplot(3, 2, 5); histogram(error_abs_history, 30, 'FaceColor', 'g', 'EdgeColor', 'k', 'FaceAlpha', 0.6); grid on;
    ax_p5 = gca; set(ax_p5, 'Color', 'w'); xlabel('|e_n|'); ylabel('Frequency'); title('Residual Error Distribution');
    
    % Panel 6: Espectro (En Inglés)
    subplot(3, 2, 6); plot(f_vec/1e3, 10*log10(psd_val), 'k-', 'LineWidth', 1); grid on;
    ax_p6 = gca; set(ax_p6, 'Color', 'w'); xlabel('kHz'); ylabel('dB/Hz'); title('Estimated Channel PSD');
    
    % Título institucional superior en inglés
    sgtitle(['HIL Hardware Instrumentation Characterization AD2 @ ' num2str(fs/1e3) ' kHz'], ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Guardado definitivo del mosaico
    exportar_grafica(fig7, '7_Panel_Resumen_Sensors');

    %% 7. EXPORTACIÓN DE HOJAS DE COMPILACIÓN EXCEL (.XLSX)
    archivo_excel = fullfile(carpeta_salida, 'Datos_Metrologicos_Sensors.xlsx');
    T_resultados = table(snr_vector', ber_sin_eq, ber_con_eq, evm_final, ...
        'VariableNames', {'SNR_dB', 'BER_Canal_Degradado', 'BER_Ecualizador_NLMS', 'EVM_Porcentaje'});
    writetable(T_resultados, archivo_excel, 'Sheet', 'Métricas Consolidadas');
    
    fprintf('\n✓ Proceso finalizado de forma exitosa. Archivos de datos compilados en:\n  %s\n', carpeta_salida);
    
catch ME
    % Control estructurado ante eventuales fallas en los registros del bus USB de la AD2
    fprintf('\n[FALLA DE HARDWARE] Línea %d: %s\n', ME.stack(1).line, ME.message);
end