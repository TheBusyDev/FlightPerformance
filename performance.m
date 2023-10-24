clc; clear; close all

%% DATI
% dati aria e combustibile
gamma = 1.4; % rapporto tra c_p e c_v
c_p = 1004; % calore specifico a pressione costante
R = 287; % costante specifica dell'aria
Q_f = 43e+6; % potere calorifico inferiore combustibile

% dati motore
e_d_vect = [1, 0.95, 0.95, 0.95]; % efficienza del diffusore (vettore)
beta_tot = 30; % rapporto di compressione totale
beta_f = 1.5; % rapporto di compressione del fan
beta_c = beta_tot / beta_f; % rapporto di compressione del compressore
n_f = 0.88; % rendimento isoentropico del fan
n_c = 0.85; % rendimento isoentropico del compressore
n_b = 0.99; % rendimento del combustore
e_b = 0.95; % efficienza del combustore
T_04 = 1400; % temperatura massima del ciclo
n_t = 0.90; % rendimento della turbina
n_n = 0.98; % rendimento dell'ugello

BPR_vect = 0:0.001:11; % vettore dei BPR

% dati di volo
z_vect = [0, 8, 10, 12]*1e+3; % vettore delle quote di volo
M_0_vect = [0, 0.8, 0.8, 0.8]; % vettore del Mach di volo
p_a_vect = [101.325, 36, 26, 19]*1e+3; % vettore delle pressioni esterne
T_a_vect = [288, 236, 223, 217]; % vettore delle temperature esterne

%% CICLO TERMODINAMICO
% grafico I-BPR
figure(1); hold on; grid on; legend on
title ("SPINTA SPECIFICA (I)")
ylabel("I [m/s]"); xlabel("BPR")
xlim([0, 11]); xticks(0:1:11)

% grafico TSFC-BPR
figure(2); hold on; grid on; legend ('on', 'Location', 'southeast')
title ("CONSUMO SPECIFICO (TSFC)")
ylabel("TSFC [kg/(h*kgf)]"); xlabel("BPR")
xlim([0, 11]); xticks(0:1:11)

% grafico n_th-BPR
figure(3); hold on; grid on; legend on
title ("RENDIMENTO TERMODINAMICO (\eta_{th})")
ylabel("\eta_{th}"); xlabel("BPR")
xlim([0, 11]); xticks(0:1:11)

% grafico n_p-BPR
figure(4); hold on; grid on; legend ('on', 'Location', 'southeast')
title ("RENDIMENTO PROPULSIVO (\eta_{p})")
ylabel("\eta_{p}"); xlabel("BPR")
xlim([0, 11]); xticks(0:1:11)

% grafico n_o-BPR
figure(5); hold on; grid on;  legend ('on', 'Location', 'southeast')
title ("RENDIMENTO GLOBALE (\eta_{o})")
ylabel("\eta_{o}"); xlabel("BPR")
xlim([0, 11]); xticks(0:1:11)

% grafico p_cr_1-BPR
figure(6); hold on; grid on; legend on
title ("PRESSIONE CRITICA (p^*), PUNTO FISSO")
ylabel("Pressione [bar]"); xlabel("BPR")

% grafico v_9-BPR
figure(7); hold on; grid on; legend on
title ("VELOCITÀ DI EFFLUSSO UGELLO PRIMARIO (v_9)")
ylabel("v_9 [m/s]"); xlabel("BPR")
xlim([0, 11]); xticks(0:1:11); clc

% temperatura in una trasformazione adiabatica isoentropica
is_temp = @(T, beta) T * beta^((gamma-1) / gamma);

% pressione in una trasformazione adiabatica isoentropica
is_pres = @(p, temp_ratio) p * temp_ratio^(gamma / (gamma-1));

for i = 1:length(z_vect) % vario la quota
    z = z_vect(i);
    M_0 = M_0_vect(i);
    p_a = p_a_vect(i);
    T_a = T_a_vect(i);
    e_d = e_d_vect(i);

    for j = 1:length(BPR_vect) % faccio variare il BPR
        BPR = BPR_vect(j);

        % velocità di volo
        v_0 = M_0 * sqrt (gamma * R * T_a);

        % diffusore 
        T_0a = T_a * (1 + (gamma-1)/2 * M_0^2);
        p_0a = p_a * (1 + (gamma-1)/2 * M_0^2)^(gamma / (gamma-1));

        T_02 = T_0a;
        p_02 = e_d * p_0a;

        % fan
        p_021 = beta_f * p_02;
        T_021_is = is_temp (T_02, beta_f);
        T_021 = T_02 + (T_021_is - T_02) / n_f;

        % compressore
        p_03 = beta_c * p_021;
        T_03_is = is_temp (T_021, beta_c);
        T_03 = T_021 + (T_03_is - T_021) / n_c;

        % combustore
        f = (c_p * (T_04 - T_03)) / (n_b*Q_f - c_p*T_04);
        p_04 = e_b * p_03;

        % turbina
        T_05 = T_04 - 1 / (1+f) * ((1 + BPR)*(T_021 - T_02) + (T_03 - T_021));
        T_05_is = T_04 - (T_04 - T_05) / n_t;
        p_05 = is_pres (p_04, T_05_is / T_04);

        % ugello primario
        [p_cr_1(i).val(j), p_9, T_9, v_9(i).val(j)] = nozzle (gamma, c_p, n_n, p_a, p_05, T_05);

        % ugello secondario
        [p_cr_2(i).val(j), p_19, T_19, v_19] = nozzle (gamma, c_p, n_n, p_a, p_021, T_021);

        if ~isreal(v_9(i).val(j)) || ~isreal(v_19)
            v_9(i).val(j) =    []; % cancella l'ultimo valore di 'v_9'
            p_cr_1(i).val(j) = []; % cancella l'ultimo valore di 'p_cr_1'
            p_cr_2(i).val(j) = []; % cancella l'ultimo valore di 'p_cr_2'
            j = j-1;

            fprintf ('\nMASSIMO BPR: %.3f\n', BPR)
            break
        end

        %% PRESTAZIONI
        % spinta diviso la portata del flusso primario (T / m_a1)
        t = ((1+f) * v_9(i).val(j) - v_0) + BPR * (v_19 - v_0);

        % spinta specifica
        I(i).val(j) = t / (1 + BPR);

        % consumo specifico
        TSFC(i).val(j) = f / t * 3600 * 9.8;

        % variazione di energia cinetica specifica
        e_kin = 0.5 * ( ((1+f) * v_9(i).val(j)^2 - v_0^2) + BPR * (v_19^2 - v_0^2) );

        % rendimento termodinamico
        n_th(i).val(j) = e_kin / (f * Q_f);

        % rendimento propulsivo
        n_p(i).val(j) = t * v_0 / e_kin;

        % rendimento globale
        n_o(i).val(j) = n_th(i).val(j) * n_p(i).val(j);

        if BPR == 5 % salva i risultati per un BPR scelto
            temp(:,i) = [ % salva le temperature
                T_a
                T_02
                T_021
                T_03
                T_04
                T_05
                T_9
                T_19
            ];

            pres(:,i) = [ % salva le pressioni
                p_a
                p_02
                p_021
                p_03
                p_04
                p_05
                p_9
                p_19
            ];

            speed(:,i) = [ % salva le velocità
                v_0
                v_9(i).val(j)
                v_19
            ];

            perf(:,i) = [ % salva le prestazioni
                I(i).val(j)
                TSFC(i).val(j)
                n_th(i).val(j)
                n_p(i).val(j)
                n_o(i).val(j)
            ];
        end
    end

    %% PLOTTO GRAFICI
    % prestazioni
    BPR_plot = BPR_vect (1:j);
    label = sprintf ('Quota: %d m', z);

    figure(1)
    plot (BPR_plot, I(i).val, "LineWidth", 1.5, "DisplayName", label)

    figure(2)
    plot (BPR_plot, TSFC(i).val, "LineWidth", 1.5, "DisplayName", label)

    figure(3)
    plot (BPR_plot, n_th(i).val, "LineWidth", 1.5, "DisplayName", label)

    figure(4)
    plot (BPR_plot, n_p(i).val, "LineWidth", 1.5, "DisplayName", label)

    figure(5)
    plot (BPR_plot, n_o(i).val, "LineWidth", 1.5, "DisplayName", label)

    % velocità di efflusso
    figure(7)
    plot (BPR_plot, v_9(i).val, "LineWidth", 1.5, "DisplayName", label)

    if z == 0
        % plotto grafico della della pressione critica
        p_a_plot = p_a * ones (1,j) / 1e+5;

        figure(6)
        % pressioni critiche ugello primario e secondario
        plot (BPR_plot, p_cr_1(i).val / 1e+5, "LineWidth", 1.5, "DisplayName", "Pressione critica ugello primario p^*_1")
        plot (BPR_plot, p_cr_2(i).val / 1e+5, "LineWidth", 1.5, "DisplayName", "Pressione critica ugello secondario p^*_2")

        % pressione ambiente
        plot (BPR_plot, p_a_plot, "LineWidth", 1.5, "DisplayName", "Pressione ambiente p_a")
    end
end



function [p_cr, p_out, T_out, v_out] = nozzle (gamma, c_p, n_n, p_a, p_in, T_in)

    p_cr = p_in * (2 / (gamma + 1))^(gamma / (gamma-1)); % pressione critica
    p_out = p_a; % suppongo ugello sempre adattato

    T_out_is = T_in * (p_out / p_in)^((gamma-1) / gamma);
    T_out = T_in - n_n * (T_in - T_out_is);
    
    v_out = sqrt (2 * c_p * (T_in - T_out));

end
