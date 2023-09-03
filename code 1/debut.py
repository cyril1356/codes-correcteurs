import numpy as np
matrice_test=np.array([[1, 0, 0, 0] ,[ 0, 1, 0, 0],[ 0 ,0 ,1, 0 ],[ 0, 0, 0, 1 ], [1 ,0 ,1 ,1 ],[ 1, 1, 0, 1 ], [1, 1, 1, 0 ]])
matrice_test=matrice_test.transpose()
print(matrice_test)
tab_valeurs=[]
def modulo_matrice(matrice):
    for i in range (len(matrice)):
        for j in range(len(matrice[0])):
            matrice[i][j]=matrice[i][j] %2 

    return (matrice)

def modulo_vecteur(matrice):
    for i in range (len(matrice)):
            matrice[i]=matrice[i]%2 

    return (matrice)

n=2
for i in range (n):
    for j in range (n):
        for k in range (n):
            for l in range (n):
                valeur=np.array([i,j,k,l])
                tab_valeurs.append([modulo_vecteur(np.matmul(valeur,matrice_test))[-3:],valeur])



               
for element1 in tab_valeurs: 
    for element2 in tab_valeurs:
        if not np.array_equiv(element1,element2):
            dist=element1-element2
            res=0
            for element in dist : 
                element= element%2
                res+=element
            if res==2: 
                print("aa")
  
