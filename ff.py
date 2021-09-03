def nb_detections(Testimated_with_imputed: list, T_measured:list):
    #
    nb = 0
    for i in range(len(T_measured)):
        th = abs(Testimated_with_imputed[i]-T_measured[i])
        # print(th)
        if th>=0.5:
            nb= nb + 1
            #print('nb=', nb)
    return nb

# somme de difference de signal pour chaque variable avec données imputés
# nb_detections entre deux signaux: données mesurée de T et T estimée au niveau des zones imputés