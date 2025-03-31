clc
clear all
close all

% MATLAB script to launch a fiber structure generation and the corresponding LBM simulation
%
%INPUT VARIABLES:
%
% SEED: integer representing the seed for initializing the random
% generator. If seed=0, automatic seed generation. If you want to reproduce
% the same fiber structure, use the same seed (fibers will be located at the same place). 
%
% MEAN_D: contains the mean fiber to be used
%
% STD_D: contains the standard deviation of the fiber diameters
%
% PORO: estimated porosity of the fiber structure to be generated
% 
% NX: domain lateral size in grid cell

seed=102;
% Paramètres fixés
NX = 200;
dx = 1e-6;
seed_base = 100;
N = 100;  % Nombre de tirages Monte Carlo
mean_poro = 0.9;
std_poro = 0.0075;
mean_fiber_d = 12.5;
std_d = 2.85;
deltaP = 0.1;
k_results = zeros(1, N);
poro_samples = normrnd(mean_poro, std_poro, [1, N]);
filename= 'fiber_mat.tiff' ;
for i = 1:N
    
    poro_i = poro_samples(i);
    % filename = sprintf('fiber_mat_poro_%d.tiff', i);
 
    % Génération de la structure et simulation
    d_equivalent = Generate_sample(i, filename, mean_fiber_d, std_d, poro_i, NX, dx);
    k_results(i) = LBM(filename, NX, deltaP, dx, d_equivalent);
    i
end

%Distribution log-normale
log_k = log(k_results);
mean_log = mean(log_k);
std_log = std(log_k);
median_k = exp(mean_log)


u_input_min = median_k - exp(mean_log - std_log);
u_input_plus = exp(mean_log + std_log) - median_k;

E = median_k - 80.6;

fprintf('u_input- = %.2f µm², u_input+ = %.2f µm²\n', u_input_min, u_input_plus);

fprintf('E (erreur de simulation) = S - D = %.2f µm²\n', E);

% Estimation de la densité (PDF) de la porosité
x_vals = linspace(0.87, 0.93, 200);  % domaine réaliste autour de la moyenne ± 3σ
[f_pdf, x_pdf] = ksdensity(poro_samples, x_vals);

% Estimation de la CDF
[f_cdf, x_cdf] = ksdensity(poro_samples, x_vals, 'Function', 'cdf');

% Tracé de la PDF
figure;
plot(x_pdf, f_pdf,'b', 'LineWidth', 2);
xlabel('Porosité \epsilon', 'FontSize', 12);
ylabel('Densité de probabilité', 'FontSize', 12);
title('PDF estimée de la porosité d’entrée', 'FontSize', 14);
xlim([0.87 0.93]);
grid on;

% Tracé de la CDF
figure;
plot(x_cdf, f_cdf, 'b', 'LineWidth', 2);
xlabel('Porosité \epsilon', 'FontSize', 12);
ylabel('Fonction de répartition cumulative (CDF)', 'FontSize', 12);
title('CDF estimée de la porosité d’entrée', 'FontSize', 14);
xlim([0.87 0.93]);
ylim([0 1]);
grid on;

% Affichage de la PDF ou histogramme
figure;
histogram(log(k_results), 'Normalization', 'pdf');
xlabel('log(k)');
ylabel('Densité de probabilité');
title('Distribution log-normale de la perméabilité');

% Tracé de la CDF de k_results
figure;
sorted_k = sort(k_results);
cdf_values = (1:length(sorted_k)) / length(sorted_k);
plot(sorted_k, cdf_values, 'LineWidth', 2);
xlabel('Perméabilité k [µm²]');
ylabel('Fonction de répartition cumulative (CDF)');
title('CDF de la perméabilité simulée');
grid on;

%Calcul de Ud

br = 10.0;     % biais systématique (µm²)
sr = 14.7;     % incertitude aléatoire (µm²)

u_D = sqrt(br^2 + sr^2);
fprintf('u_D = %.2f µm²\n', u_D);


%Calcul de Unum

%GCI =
u_num = 0.4642379591;

%Calcul de \delta_{model}
k_val = 2;  % pour un intervalle à 95.4 %

uval_minus = sqrt(u_input_min^2 + u_num^2 + u_D^2);
uval_plus  = sqrt(u_input_plus^2  + u_num^2 + u_D^2);

delta_model_min = E - k_val * uval_plus;
delta_model_max = E + k_val * uval_minus;

fprintf('Intervalle de delta_model : [%.2f , %.2f] µm²\n', delta_model_min, delta_model_max);