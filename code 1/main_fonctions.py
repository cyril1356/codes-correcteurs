import numpy as np
import random
from copy import deepcopy
from tqdm import tqdm

### cette partie contient toute les fonction utilisées durant le stage

def encode(self): 
    #transfore une chaine de char en bits
    L=[]
    for element in self:
        L.append('0'+format(ord(element), 'b'))
    return L

def incetion_vecteur_espace(tab):
    #remplace les valeur dans un vecteur par les représent des classes d'equivalences
    for i in range (len(tab)):
        tab[i]=tab[i]%2
    return(tab)
    
def crea_tab_sindorme(H):
    ##construit le tableau syndrome pour les codes de hamming 4-7
    v=2
    sindrome=[]
    for i in range (v):
        for j in range (v):
            for k in range (v):
                for l in range (v):
                    for m in range (v):
                        for n in range (v):
                            for o in range (v):
                                ok =True
                                ok_2=False
                                conpteur=111
                                tab=np.array([i,j,k,l,m,n,o])
                                tab=tab.transpose()
                                val=incetion_vecteur_espace( np.dot(H,tab))
                                for p in range (len(sindrome)) :
                                    if np.array_equal(sindrome[p][0],val):
                                        
                                        
                                        if poids(sindrome[p][1])<= poids(tab.T):
                                            ok=False
                                        else:
                                            ok_2=True
                                            conpteur=p

                                if ok: 

                                    sindrome.append([val, tab.T])
                                    if ok_2:
                                        del sindrome[conpteur]
    return(sindrome)

def crea_tab_sindorme_new(H,n):
    ## construi un tableau syndrome pour tout codes (généralisation de la fonction précedante )
    v=2
    sindrome=[]

    v=2
    sindrome=[]
    for i in range (2**n):##on vas faire tout les élément comme ca on est sur 
        tab=np.zeros(n)
        val_bin=change_valeur_bit(i,n)
        for j in range (len(tab)):
            tab[j]=val_bin[j]
        ok =True
        ok_2=False
        conpteur=111
        tab=tab.transpose()
        val=incetion_vecteur_espace( np.dot(H,tab))
        for p in range (len(sindrome)) :
            if np.array_equal(sindrome[p][0],val):
                
                
                if poids(sindrome[p][1])<= poids(tab.T):
                    ok=False
                else:
                    ok_2=True
                    conpteur=p

        if ok: 

            sindrome.append([val, tab.T])
            if ok_2:
                del sindrome[conpteur]
    return(sindrome)

def encode_decoupe(decoupage,G):
    ##permet d'appliqueée la fonction phi à une chaine de char déja découpée en vecteur de tailles k
    new=[]
    for i in range (len(decoupage)):
        new.append(modulo_vecteur(np.dot(decoupage[i],G)))
    return(new)



def change_bit_valeur(chaine):
    ##focntion pour transformer un nombre en bits 
    valeur=0
    compteur=0
    rev = list(reversed(chaine))
    for element in rev:
        if element=='1':
            valeur+=2**compteur
        compteur+=1
    return valeur

def change_valeur_bit(valeur,longeur):
    ##focntion pour changer une chaine de char de bit en un nombre
    L=np.array([])
    
    mot=format(valeur,'b')
    assert len(mot)<=longeur

    while len(mot)<longeur:
        mot='0'+mot

    for element in mot:
        L=np.hstack((L,np.array([int(element)])))
    return L

def calcul_valeur_n_et_k(valeurfin):
    #calcul le tableau des valeur possible entre n et k pour le code d'Hamming
    L=[]
    compteur=1
    for i in range(1,valeurfin):
        x=i-np.log2(i+1)
        if int(x)==x:
            L.append((int(x),i))

    return L 

def combinaisons_v1(p,n): 
    ## rend une liste de tuples représenants les différents sous ensemble de taille p de {0...n}
    liste_combinaisons=[]    
    indices=list(range(p)) 
  
    liste_combinaisons.append(tuple(indices))
  
    if p==n: 
        return liste_combinaisons 
  
    i = p-1 
    while (i != -1): 
  
        indices[i] += 1 
  
        for j in range(i + 1, p): 
            indices[j] = indices[j - 1] + 1            
  
        if indices[i] == (n  - p + i): 
            i -= 1 
        else: 
            i = p - 1 
  
        liste_combinaisons.append(tuple(indices)) 
  
    return liste_combinaisons


def construction_martice_generalise_forme_systematique(k,n):
    #construi les matrices pour le codes de Hamming pour tout n et k 
    assert 2**n==2**k*(n+1)
    G=np.identity(k)
    A=np.zeros((k,n-k))
    combinaisons=[]
    for i in range (2,n-k+1):
        combinaisons+=combinaisons_v1(i,n-k)
    for i in range (len(A)):
        for indices in combinaisons[i]:
            A[i][indices]=1
    G=np.concatenate((G,A),axis=1)##ok 
    return(G)

def construction_H_via_G(G):
    ##construction de H si G est sous la forme systématique 
    forme=np.shape(G)
    k=forme[0]
    n=forme[1]
    A=np.array([G[0][k:]])

    for i in range(1,k):
        A=np.concatenate((A,[G[i][k:]]),axis=0)        
    H=A.T
    H=np.concatenate((H,np.identity(n-k)),axis=1)
    return(H)



def val_addition(L,v):
    x=np.zeros(len(L[0]))
    for i in range (len(L)):
        x+=int(v[i])*L[i]
    return(x)


def calucul_distance(H):
    #calcul la distance d'un code corrécteur sachant H
    ##attention cela ne marche qu'avec la matrice H
    nb_comnones=np.shape(H)[1]
    compteur=2
    ok=True
    while ok: 
        compteur+=1
        combi=combinaisons_v1(compteur,nb_comnones)
        for i in tqdm(range(len(combi))):
            L=[]

            for element2 in combi[i]:
                L.append(H[:,element2])
            for i in range (1,2**len(L)):
                v=change_valeur_bit(i,len(L[0]))
                v=np.flipud(v)
                res=val_addition(L,v)
                res=modulo_vecteur(res)
                
                if np.array_equal(res,np.zeros(len(L[0]))):

                    return (compteur,L)


def calucul_distance_avec_poids_minimal(H_reduit,poid_minimal):
    ##Fonction utile pour calculer la distance en utilisant l'ensemble E_x introduit dans le théoreme 3.4
    nb_comnones=np.shape(H_reduit)[1]
    compteur=0
    ok=True
    while ok: 
        compteur+=1
        combi=combinaisons_v1(compteur,nb_comnones)
        for i in tqdm(range(len(combi))):
            L=[]

            for element2 in combi[i]:
                L.append(H_reduit[:,element2])
            for i in range (1,2**len(L)):
                v=change_valeur_bit(i,len(L[0]))
                v=np.flipud(v)
                res=val_addition(L,v)
                res=modulo_vecteur(res)
                
                if poids(res)>=poid_minimal:

                    return (compteur,L)



def decode(self):
    ##transfore une chaine de bit en ine chaine de char
    mot = ''
    for element in self: 

        mot+=chr(change_bit_valeur(element))
    return(mot)


def petit_decoupage(mot,n):
    ##utile pour grand découpage
    assert len(mot)%n==0
    res=[]
    val=''
    for element in mot:
        val+=element
        if len(val)%n ==0:
            res.append(val)
            val=''
    return res

def poids(T1):
    #rend le poid d'un vecteur
    new=T1
    p=0
    for element in new :
        if element!=0:
            p+=1

    return (p)



def grand_decoupage_numpy(list_mot,n):
    ##découpe un string en sous string de taille n 
    res=[]
    for element in list_mot:
        petite_matrice=petit_decoupage(element,n)
        for petit_mot in petite_matrice:
            mat=np.array([])
            for lettre in petit_mot :
                mat=np.append(mat,int(lettre))
            res.append(mat)
    return (res)


def brouillage(tab,proba):
    #change un bit de notre mot avec une chance 1/proba
    new = deepcopy(tab)
    
    for j in range (len(new)):
        
        for i in range (len(new[j])):
            v=random.randint(0,proba)
            if v==1:
                new[j][i]=(new[j][i]+1)%2
    return new

def recherche_closet_leader(synd,val):
    #effectue la recherche du mot dans tableau pour aplliquer de décodage par syndrome 
    for element in synd:
        if np.array_equal(element[0],val):
            return (element[1])

def reparation(lst_mot_casse,sindrome,H):
    #effectue tout le décodeage par syndrome 
    nouveau_mot=[]
    for mot in lst_mot_casse:

        val = np.dot(H,mot.T)
        val=modulo_vecteur(val)

        e=recherche_closet_leader(sindrome,val)
        nouveau_mot.append(modulo_vecteur(mot-e))
    return(nouveau_mot)


def modulo_matrice(matrice):
    #remplace les valeur dans une matrice par les représent des classes d'equivalences
    for i in range (len(matrice)):
        for j in range(len(matrice[0])):
            matrice[i][j]=matrice[i][j] %2 

    return (matrice)

def modulo_vecteur(matrice):
    #remplace les valeur dans un vecteur par les représent des classes d'equivalences
    for i in range (len(matrice)):
            matrice[i]=matrice[i]%2 

    return (matrice)


def decriptage(liste_mot):
    #changes les bit en mots 
    new=[]
    assert (len(liste_mot))%2==0 

    for j in range(len(liste_mot)):
        new+=[""]
        lettre1=liste_mot[j][:4]
        lettre=lettre1
        for j in lettre1:
            new[-1]=new[-1]+str(int(j))

    return(new)
        
    
def retour_mot(mot_decripter):

    mot=[]
    for j in range (len(mot_decripter)//2):
        mot.append(mot_decripter[2*j]+mot_decripter[2*j+1])
    return(decode(mot))

def sigma(vecteur):
    new=np.zeros(len(vecteur))
    for i in range (len(vecteur)-1):
        new[i+1]=vecteur[i]
    new[0]=vecteur[-1]
    return new



def construction_matrice_code_cycl(g,n):
    #genére la matrice G d'un code cyclique
    assert (len(g))<=n
    v=len(g)
    k=n-len(g)+1##attention traitre 
    while len(g)!=n:
        g=np.concatenate((g,np.array([0])))
    matrice =np.array([g])
    for i in range (2,k+1):
        new_g=sigma(g)
        g=new_g

        matrice=np.concatenate((matrice,[new_g]),axis=0)
    return matrice


def test_matrice_lign(matrice1,matrice2):#on cherche a voir si deux matrices sont equivalente a permutation des lignes prés 
    if np.shape(matrice1)!=np.shape(matrice2):#si une des deux matrice a 2 même lignes ca marche pas 
        print("les matrices n'ont pas la même taille")
        return False
    compteur=0
 
    for i in range (np.shape(matrice1)[0]):
        ok=True
        ligne1=matrice1[i,:]
        for j in range(np.shape(matrice2)[0]):
            ligne2=matrice2[j,:]
            if np.array_equal(ligne1,ligne2) and ok:
                
                compteur +=1
                ok =False
    if compteur==np.shape(matrice1)[0]:
        return True
    return False

def test_matrice_colone(matrice1,matrice2):#on cherche a voir si deux matrices sont equivalente a permutation des lignes prés 
    if np.shape(matrice1)!=np.shape(matrice2):#si une des deux matrice a 2 même lignes ca marche pas 
        print("les matrices n'ont pas la meme taille")
        return False
    compteur=0
 
    for i in range (np.shape(matrice1)[1]):
        ok=True
        ligne1=matrice1[:,i]
        for j in range(np.shape(matrice2)[1]):
            ligne2=matrice2[:,j]
            if np.array_equal(ligne1,ligne2) and ok:
                
                compteur +=1
                ok =False
    if compteur==np.shape(matrice1)[1]:
        return True
    return False

def test_matrice_equiv_ligne(mat1,mat2):
    if test_matrice_lign(mat1,mat2) and test_matrice_lign(mat2,mat1):
        return True
    else: 
        return False
    
def test_matrice_equiv_colone(mat1,mat2):
    if test_matrice_colone(mat1,mat2) and test_matrice_colone(mat2,mat1):
        return True
    else: 
        return False
    
def construction_h_code_cycl(liste,n):
    ##construction de la matrice H d'une coge cyclique via la polynome g
    liste=np.flip(liste)
    return construction_matrice_code_cycl(liste,n)

    