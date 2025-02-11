# Importation des modules
import numpy as np


def mdf1(prm):
    """Fonction qui calcule le profil de concentration de sel le long du rayon d'un pilier 
    // Ordre 1 pour la dérivée simple 

    Entrées:
        - prm : Objet class parametres()
            - D : Diamètre
            - D_eff :  Coefficient de diffusion effectif
            - S : Constante de de consommation
            - C_e : Concentration du sel à la surface
            - N : Nombre de points utilisés pour la méthode

    Sorties (dans l'ordre énuméré ci-bas):
        - Vecteur (np.array) donnant la concentration de sel le long du rayon du pilier
        - Vecteur (np.array) donnant la position tout au long du pilier (axe r) en mètre
    """

    # Fonction à écrire
    dr=prm.D/(2*(prm.N-1))      #Taille des intervalles
    position=np.zeros(prm.N)    #Creation du tableau qui contient les positions
    for i in range(len(position)):
        position[i]=i*dr

    # Matrice des coefficient
    M_C=np.zeros((prm.N, prm.N))

    #Condition frontière à r=0
    M_C[0,0]=-1
    M_C[0,1]=1 

    #Condition frontière à r=R
    M_C[-1,-1]=1

    #Dans le domaine
    for i in range(1,prm.N-1):
        M_C[i,i+1]=prm.D_eff/(position[i]*dr) + prm.D_eff/(dr**2)
        M_C[i,i]=-prm.D_eff/(position[i]*dr) - 2*prm.D_eff/dr**2
        M_C[i,i-1]=prm.D_eff/dr**2

    # Vecteur du terme de droite
    V_C=np.zeros(prm.N)
    V_C[prm.N-1]=prm.C_e
    for i in range(1,len(V_C)-1):
        V_C[i]=prm.S
    V_C=np.transpose(V_C)
    C=np.linalg.solve(M_C,V_C)
    return C,position


def mdf2(prm):
    """Fonction qui calcule le profil de concentration de sel le long du rayon d'un pilier
    // Ordre 2 pour la dérivée simple 
    Entrées:
        - prm : Objet class parametres()
            - D : Diamètre
            - D_eff :  Coefficient de diffusion effectif
            - S : Constante de de consommation
            - C_e : Concentration du sel à la surface
            - N : Nombre de points utilisés pour la méthode

    Sorties (dans l'ordre énuméré ci-bas):
        - Vecteur (np.array) donnant la concentration de sel le long du rayon du pilier
        - Vecteur (np.array) donnant la position tout au long du pilier (axe r) en mètre
    """

    # Fonction à écrire
    dr=prm.D/(2*(prm.N-1))
    position=np.zeros(prm.N)
    for i in range(len(position)):
        position[i]=i*dr

    # Matrice des coefficient
    M_C=np.zeros((prm.N, prm.N))

    #Condition frontière à r=0
    M_C[0,0]=-3
    M_C[0,1]=4 
    M_C[0,2]=-1 
    #Condition frontière à r=R
    M_C[-1,-1]=1

    #Dans le domaine
    for i in range(1,prm.N-1):
        M_C[i,i+1]=prm.D_eff/(2*position[i]*dr) + prm.D_eff/(dr**2)
        M_C[i,i]= -2*prm.D_eff/dr**2
        M_C[i,i-1]=-prm.D_eff/(2*position[i]*dr) + prm.D_eff/dr**2

    # Vecteur du terme de droite
    V_C=np.zeros(prm.N)
    V_C[prm.N-1]=prm.C_e
    for i in range(1,len(V_C)-1):
        V_C[i]=prm.S
    V_C=np.transpose(V_C)
    C=np.linalg.solve(M_C,V_C)
    return C,position

def Erreur_L1(solution_analytique,mdf,prm):
    
    return 1 / prm.N * np.sum(np.abs(solution_analytique - mdf))

def Erreur_L2(solution_analytique,mdf,prm):
    
    return np.sqrt(1/prm.N * np.sum((solution_analytique - mdf)**2)) 

def Erreur_Linf(solution_analytique,mdf,prm):
    
    return np.max(np.abs(solution_analytique - mdf))


def calcul_ordre(erreurs, N_points):
    
    erreur_1, erreur_2 = erreurs[-2], erreurs[-1]
    N_1, N_2 = N_points[-2], N_points[-1]

    return np.log(erreur_2 / erreur_1) / np.log(N_1 / N_2)

