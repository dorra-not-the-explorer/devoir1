clc;
clear all;
close all;

% Paramètres d'entrée
deltaP = 0.1;             % pression en Pa
NX = 300;                 % taille de domaine de base
poro = 0.9;
mean_fiber_d = 12.5;      % µm
std_d = 2.85;             % µm
filename = 'fiber_mat.tiff';
NXdx = 2e-4;

% Grille raffinée :
Nx_values = linspace(50, 170, 25);  % ← 25 valeurs de raffinements
Nx_values = [Nx_values,200];
N = length(Nx_values);
dx_values = NXdx ./ Nx_values;

% Nombre de seeds (10 seeds = 127 à 136)
Nb_seed = 20;
k_computed = zeros(N, Nb_seed);  % N dx x Nb_seed
r = 2;

% Boucle sur les seeds
for j = 127:(127 + Nb_seed - 1)
    seed = j;
    for i = 1:N
        fprintf('Seed = %d | NX = %d | dx = %.2e\n', seed, Nx_values(i), dx_values(i));
        % Génération de la structure
        d_equivalent = Generate_sample(seed, filename, mean_fiber_d, std_d, poro, Nx_values(i), dx_values(i));

        % Simulation LBM
        k_computed(i, j - 126) = LBM(filename, Nx_values(i), deltaP, dx_values(i), d_equivalent);
    end
end

% Calcul de l'erreur  par rapport au maillage le plus fin (dernier dx)

k_ref = k_computed(end, :);  % dernière ligne : plus fin dx
erreur = zeros(N - 1, Nb_seed);

for l = 1:(N - 1)
    erreur(l, :) = abs(k_computed(l, :) - k_ref) ;%./ abs(k_ref);
end

% Affichage en log-log
figure(1); clf; hold on;
colors = lines(Nb_seed);

for j = 1:Nb_seed
    loglog(dx_values(1:end-1), erreur(:, j), '-o', 'Color', colors(j, :), 'DisplayName', ['Seed ', num2str(j + 126)]);
end

xlabel('dx');
ylabel('Erreur relative');
title('Erreur relative en fonction de dx (log-log)');
legend show;
grid on;







