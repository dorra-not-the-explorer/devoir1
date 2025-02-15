import numpy as np
import matplotlib.pyplot as plt

try:
    from pilier_fct import *
except:
    pass

#------------------------------------------------------------------------------
# Code principal pour l'analyse des résultats
# Il faudra faire appel aux fonctions programmées dans pilier_fct.py .
#------------------------------------------------------------------------------

# Assignation des paramètres
class parametres():
    D_eff = 10**-10       # [m^2/s] diffusion
    S = 2*(10**-8)       # [mol/m^3/s] consommation
    D = 1      # [m] Diamètre
    C_e = 20    # [mol/m^3] Concentration de surface
    N = 5       # [-] Nombre de points en r

# Nombre de points pour l'analyse de convergence
N_points = [5, 10, 20, 40, 80, 160, 320, 640]
L1_modele1, L2_modele1, Linf_modele1 = [], [], []
L1_modele2, L2_modele2, Linf_modele2 = [], [], []
dr_values = [] 
#Calcul des erreurs pour chaque modèle et chaque nombre de points

for N in N_points:
    
    prm = parametres()
    prm.N = N
    prm.dr = prm.D / (2 * (N - 1))  # Calcul de dr
    dr_values.append(prm.dr)  

    #Calcul de la solution analytique
    r = np.linspace(0, prm.D/2, prm.N)
    C_analytique = (1/4)*(prm.S/prm.D_eff)*((prm.D/2)**2) *( (r**2 / ((prm.D/2)**2)) -1 ) + prm.C_e
    
    # Modèle avec erreur de troncature d'ordre 1
    C_mdf1, _ = mdf1(prm)
    L1_modele1.append(Erreur_L1(C_analytique, C_mdf1, prm))
    L2_modele1.append(Erreur_L2(C_analytique, C_mdf1, prm))
    Linf_modele1.append(Erreur_Linf(C_analytique, C_mdf1, prm))

    # Modèle avec erreur de troncature d'ordre 2
    C_mdf2, _ = mdf2(prm)
    L1_modele2.append(Erreur_L1(C_analytique, C_mdf2, prm))
    L2_modele2.append(Erreur_L2(C_analytique, C_mdf2, prm))
    Linf_modele2.append(Erreur_Linf(C_analytique, C_mdf2, prm))


#Tracé des erreurs



# Modèle 1

plt.figure(1)

plt.loglog(dr_values, L1_modele1, 'o-', label='L1 (Modèle 1)')
plt.loglog(dr_values, L2_modele1, 's-', label='L2 (Modèle 1)')
plt.loglog(dr_values, Linf_modele1, '^-', label='Linf (Modèle 1)')


plt.xlabel('Taille des éléments (m)')
plt.ylabel('Erreur')
plt.title('Erreurs (L1, L2, Linf) en fonction de la taille des éléments pour le modèle 1')
plt.legend()
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.show()


# Modèle 2

plt.figure(2)

plt.loglog(dr_values, L1_modele2, 'o--', label='L1 (Modèle 2)')
plt.loglog(dr_values, L2_modele2, 's--', label='L2 (Modèle 2)')
plt.loglog(dr_values, Linf_modele2, '^--', label='Linf (Modèle 2)')


plt.xlabel('Taille des éléments (m)')
plt.ylabel('Erreur')
plt.title('Erreurs (L1, L2, Linf) en fonction de la taille des éléments pour le modèle 2')
plt.legend()
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.show()


# Ordres pour le modèle 1
ordre_L1_modele1 = calcul_ordre(L1_modele1, N_points)
ordre_L2_modele1 = calcul_ordre(L2_modele1, N_points)
ordre_Linf_modele1 = calcul_ordre(Linf_modele1, N_points)

# Ordres pour le modèle 2
ordre_L1_modele2 = calcul_ordre(L1_modele2, N_points)
ordre_L2_modele2 = calcul_ordre(L2_modele2, N_points)
ordre_Linf_modele2 = calcul_ordre(Linf_modele2, N_points)

print("Ordres de convergence (Modèle 1) :")
print(f"L1 : {ordre_L1_modele1}")
print(f"L2 : {ordre_L2_modele1}")
print(f"Linf : {ordre_Linf_modele1}")

print("\nOrdres de convergence (Modèle 2) :")
print(f"L1 : {ordre_L1_modele2}")
print(f"L2 : {ordre_L2_modele2}")
print(f"Linf : {ordre_Linf_modele2}")
