import numpy as np
import matplotlib.pyplot as plt
import math
import two_weeks_data
import scipy
from scipy.optimize import rosen, differential_evolution
import datetime
import scipy.stats
def simulate(structure_temperature:list, office_power_gains:list, Tout:list, Tcor:list): #avec x est Toffice
    #Les matrices
    matrixA_thermics = -2.1578228994330876e-06
    matrixB_thermics_phi_in = 3.1753190677875445e-09
    matrixB_thermics_Tout = 5.770162506326268e-08
    matrixB_thermics_Tn = 2.1001212743698246e-06
    matrixC_thermics = 0.12182074134801915
    matrixD_thermics_phi_in = 0.001292274424219787
    matrixD_thermics_Tout = 0.023483099717953423
    matrixD_thermics_Tn = 0.8546961589340274
    #
    R= 0.001292274424219787
    Ri= 0.010608
    sample_time=60*60
    alphak = math.exp(matrixA_thermics * sample_time)
    Kobserver = Ri * alphak / R
    betak = (alphak - 1) / matrixA_thermics
    office_temperatures_estimated=[]
    for k in range(len(two_weeks_data.current_datetimes)):
        #-----------créer des DM sur les entrées---------------
        structure_temerature=23*np.ones(336)
        structure_temperature= alphak * structure_temerature + betak * (matrixB_thermics_phi_in * office_power_gains[k] + matrixB_thermics_Tout * Tout[k] + matrixB_thermics_Tn * Tcor[k])
        #print('size_structure_temerature=', np.size(structure_temerature))
        office_temperature_estimated=matrixC_thermics * structure_temperature[k] + matrixD_thermics_phi_in * two_weeks_data.office_power_gains[k] + matrixD_thermics_Tout * two_weeks_data.outdoor_temperature[k] + matrixD_thermics_Tn * two_weeks_data.corridor_temperature[k]
        office_temperatures_estimated.append(office_temperature_estimated)
        return office_temperature_estimated #déduit Testimated

# Cette fonction consiste à calculer combien de fois on sort de seuil.
#def nb_detections(y1:list, y2:list):#y1
def nb_detections(y1:list, y2:list):
    y1, y2=[], []
    s1 = simulate(two_weeks_data.office_power_gains, two_weeks_data.outdoor_temperature,two_weeks_data.corridor_temperature, two_weeks_data.Toffice)
    y1.append(s1)
    for k in range(len(y1)):
        y2.append((max(y1)))
    nb = 0
    for i in range(len(y1)):
        th = abs(y1[i] - y2[i])
        # print(th)
        if th != 0:
            nb = nb + 1
            #print('nb=', nb)
    return nb
#Pour mieux comparer la somme des deux histogrammes, la méthode consiste à calculer la somme des écarts
def sum_of_difference(x:list, y:list):
    x = list()
    y= list()
    sum_of_difference = 0
    for i in range(len(x)):
        sum_of_difference += abs(x[i] - y[i])
    return sum_of_difference
    # sum_of_difference = 0
    # for i in range(len(x)):
    #     sum_of_difference += abs(x[i] - y[i] )
    #     print(sum_of_difference)
    # y = simulate(office_power_gains, outdoor_temperature, corridor_temperature, structure_temperature)
    # x = simulate(two_weeks_data.office_power_gains_copie, two_weeks_data.outdoor_temperature_copie,two_weeks_data.corridor_temperature_copie, two_weeks_data.Toffice)

# PDF de la dérivée
def moyenne(tableau):
    return sum(tableau, 0.0) / len(tableau)
#La variance est définie comme la moyenne des carrés des écarts à la moyenne
def variance(tableau):
    m=moyenne(tableau)
    return moyenne([(x-m)**2 for x in tableau])

#L'écart-type est défini comme la racine carrée de la variance
def ecartype(tableau):
    return variance(tableau)**0.5
def pdf_derivative(y1, y2):#mu1, mu2, mu3, mu4, mu5, mu6, sigma_1, sigma_2, sigma_3, sigma_4, sigma_5, sigma_6
    y1=two_weeks_data.corridor_temperature_copie
    y2=two_weeks_data.corridor_temperature
    mu1=variance(y1)
    sigma1=ecartype(y1)
    mu2 =variance(y2)
    sigma2=ecartype(y2)

    a1= scipy.stats.norm.pdf(y1, mu1, sigma1) #la variable Tcor suit la loi normale (j'ai fait une étude d'adéquation de PDF pour chaque variable)
    a2= scipy.stats.norm.pdf(y2, mu2, sigma2)
    pdf_d=sum_of_difference(a1, a2)
    #
    y3 = two_weeks_data.outdoor_temperature_copie
    y4 = two_weeks_data.outdoor_temperature
    mu3 = variance(y3)
    sigma3 = ecartype(y3)
    mu4 = variance(y3)
    sigma4 = ecartype(y4)
    a3 = scipy.stats.norm.pdf(y3, mu3, sigma3)
    a4 = scipy.stats.norm.pdf(y4, mu4, sigma4)
    pdf_d34 = sum_of_difference(a3, a4)
    #
    y5 = two_weeks_data.office_power_gains_copie
    y6 = two_weeks_data.office_power_gains
    mu5 = variance(y5)
    sigma5 = ecartype(y5)
    mu6 = variance(y6)
    sigma6 = ecartype(y6)
    a5 = scipy.stats.norm.pdf(y5, mu5, sigma5)
    a6 = scipy.stats.norm.pdf(y6, mu6, sigma6)
    pdf_d56 = sum_of_difference(a5, a6)
    return pdf_d#+pdf_d34+pdf_d56

def cost_function(missing_values:list=None):
    inputs=[office_power_gains,outdoor_temperature, corridor_temperature, structure_temperature]
    #=list()
    y = simulate(office_power_gains, outdoor_temperature, corridor_temperature, structure_temperature)
    x = simulate(two_weeks_data.office_power_gains_copie, two_weeks_data.outdoor_temperature_copie, two_weeks_data.corridor_temperature_copie, two_weeks_data.Toffice)
    #x, y=list(), list()
    for i in range(len(inputs)):
        if inputs[i]==np.nan:
            missing_values.append(inputs[i])
    for j in range(len(missing_values)):
            return 15*sum_of_difference(x, y)+70*nb_detections(x, y)#+15*pdf_derivative(x, y)
# somme de difference de signal pour chaque variable avec données imputés
# nb_detections entre deux signaux: données mesurée de T et T estimée au niveau des zones imputés

office_power_gains = two_weeks_data.office_power_gains
outdoor_temperature = two_weeks_data.outdoor_temperature
corridor_temperature = two_weeks_data.corridor_temperature
structure_temperature = two_weeks_data.Toffice

#bounds=[(22.5, 22.7), (12, 19.5), (12, 19.5), (12, 19.5), (446, 500), (446, 500), (446, 500), (446, 500)]
bounds=[]
# for i in range(93):
#     bounds.append((24.8, 31.29))#avec 24.8, 31.29 representent les "Bounds" pour la variable "Tout"
#     #print(bounds)
bounds_Tcor=[]
bounds_Tout=[]
bounds_office_power_gains=[]
bb=[]
for i in range(93):
    bounds_Tcor.append((25, 31))
    bounds_Tcor.append((12, 37))
    bounds_office_power_gains.append((451, 1100))#805 et 460
    bound_s=bounds_Tcor+ bounds_Tout+bounds_office_power_gains
    bb.append((451, 1000))
result = differential_evolution(cost_function, bb, constraints=, , strategy='best1bin', maxiter=100, popsize=15, tol=0.01, mutation=(0.5, 1), recombination=0.7, seed=None, callback=None, disp=False, polish=True, init='latinhypercube')
print('result_x=', list(result.x))

