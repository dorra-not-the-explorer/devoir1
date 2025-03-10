import numpy as np
import matplotlib.pyplot as plt

# Assignation des paramètres
class parametres():
    D_eff = 1e-10  # Coefficient de diffusion
    R = 0.5       # Rayon maximal du domaine
    N = 10       # Pas spatial
    dt = 4e6     # Pas temporel
    tf = 4e9      # Temps final
    C_e=20
    K = 4e-9
prm = parametres()


def Cm(prm,r,t) :
    
    return (prm.C_e - prm.K * t) / prm.R**2 * r**2 + prm.K * t

# Fonction source
def S_func(prm,r,t):
    return prm.K * ((prm.C_e - prm.K * t) / prm.R**2 * r**2 + prm.K * t) - prm.D_eff * (4 * (prm.C_e - prm.K * t) / prm.R**2) + (prm.K / prm.R**2) * r**2 - prm.K#prm.D_eff*4*(prm.C_e-prm.K*t)/prm.R**2 -prm.K*((prm.C_e-prm.K*t)/prm.R**2 * r**2 +prm.K*t)-prm.K+prm.K/prm.R**2* r**2


def solve_diffusion_euler_implicite(prm, S_func):
    """
    Résout l'équation de diffusion en 1D en coordonnées radiales avec Euler implicite

    """
    r = np.linspace(0, prm.R, prm.N)
    dr = prm.R/(prm.N-1)
    prm.N = len(r)
    M = int(prm.tf / prm.dt) +1
    
    # Construction explicite de la matrice tridiagonale avec une boucle for
    A = np.zeros((prm.N, prm.N))
    
    for i in range(1, prm.N - 1):
        A[i, i - 1] = -prm.D_eff * prm.dt / dr**2 + prm.D_eff * prm.dt / (2 * dr * r[i])
        A[i, i] = 1 + 2 * prm.D_eff * prm.dt / dr**2 + prm.K*prm.dt
        A[i, i + 1] = -prm.D_eff * prm.dt / dr**2 - prm.D_eff * prm.dt / (2 * dr * r[i])
    
    # Conditions aux limites
    A[0, 0] = -3
    A[0, 1] = 4 
    A[0, 2] = -1   # Flux nul à r=0 avec Gear avant
    A[-1, -1] = 1  # Dirichlet à r=R
    
    # Conditions initiales
    #C = prm.C_e/prm.R**2 *r**2
    # Résolution temporelle
    sol = np.zeros((M, prm.N))
    #sol[0, :] = Cm(prm,r,0)
    #sol[:, prm.N-1 ] = prm.C_e
    
    for m in range(1, M):
        b = sol[m-1, :].copy()         # + prm.dt *  S_func(prm,r,(m-1)*prm.dt)
        b[0] = 0 # Flux nul à r=0
        b[-1] = prm.C_e  # Condition Dirichlet
        sol[m, :] = np.linalg.solve(A, b)  
    return sol, np.linspace(0, prm.tf, M), r

# Résolution
sol, temps, r = solve_diffusion_euler_implicite(prm, S_func)

# Graphique des concentrations à des temps spécifiques
plt.figure(figsize=(8,6))

time_indices = np.linspace(0, len(temps)-1, 15, dtype=int)
for t_idx in time_indices:
    if t_idx < len(temps):
        plt.plot(r, sol[t_idx], label=f"t = {temps[t_idx]:.1f} s")
        #plt.plot(r, (prm.C_e-prm.K*temps[t_idx])/prm.R**2 * r**2 + prm.K*temps[t_idx],'o', label=f"sol. analytique t = {temps[t_idx]:.1f} s")

plt.xlabel('Rayon r')
plt.ylabel('Concentration')
plt.title('Évolution de la concentration à différents instants de temps')
plt.legend()
plt.grid()
plt.show()

