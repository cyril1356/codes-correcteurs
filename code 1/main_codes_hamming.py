from main_fonctions import *

mot=input("quel mot voulez vous envoyer(sans espaces)  ")
proba=int(input("chances d'erreur qu'un bit change(valeur entiére)  "))
n=7
k=4

mot=encode(mot)

print(mot)

mot_decoupe=grand_decoupage_numpy(mot,k)

G=construction_martice_generalise_forme_systematique(k,n)

H_reduit=construction_H_via_G(G)

syndorme = crea_tab_sindorme_new(H_reduit,n)
print(syndorme)

mot_encode_decoupe=encode_decoupe(mot_decoupe,G)

mot_brouiller= brouillage(mot_encode_decoupe,proba)

V=mot_brouiller

X=reparation(mot_brouiller,syndorme,H_reduit)

print(retour_mot(decriptage(V)),"-----mot reçus")
print(retour_mot(decriptage(X)),'-----mot reparer')