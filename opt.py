from xgboost.sklearn import XGBRegressor
import numpy as np
import os
from bayes_opt import BayesianOptimization
from pymatgen import Lattice, Structure, Composition
from matminer.featurizers.composition import BandCenter
from StructureFeature import structure_symmetry, average_atomic_volume
from ComponentFeature import meredig, band_center

#enreg = XGBRegressor()

def calculate_features(id):
    features = []
    structure = Structure.from_file("../ML/cif/"+str(id) + ".cif")
    ###### Structure Feature ######
    crystal_system, centrosymmetric = structure_symmetry(structure)
    average_atom_volume, average_volume_atom = average_atomic_volume(structure) # dimension:2

    ###### Component Feature ######
    a = structure.composition
    elemental_property = meredig(Composition(a))
    bandcenter = band_center(Composition(a))
    feature = crystal_system
    feature.append(centrosymmetric)
    feature.append(average_atom_volume)
    feature.append(average_volume_atom)
    feature = feature + elemental_property + bandcenter

    features = features + feature
    features = np.array(features)
    features = features.reshape(-1,28)
    np.savetxt("features.txt", features)


def main():
    print("****************************************")
    print("*** Welcome to BayesOptimization !!! ***")
    print("* Please input your cif id:")
    cif_id = input()
    calculate_features(id=cif_id)
    print("* Features have been calculated successfully !!!")
    print("* Please input your elements:")
    global elements
    elements = input()
    elements = elements.split(" ")
    global element_number
    element_number = len(elements)
    if element_number == 2:
        print("* Two elements are obtained successfully !!!")
        print("* Please input the fraction range of each element:")
        fraction = input()
        fraction = fraction.split(" ")
        band_opt(fraction)
    elif element_number == 3:
        print("Three elements are obtained successfully !!!")
        print("* Please input the fraction range of each element:")
        fraction = input()
        fraction = fraction.split(" ")
        band_opt(fraction)
    elif element_number == 4:
        print("Four elements are obtained successfully !!!")
    elif element_number == 5:
        print("Five elements are obtained successfully !!!")
    #print("Please input ")



def model(feature):
    enreg = XGBRegressor()
    enreg.load_model('xgb-1.model')
    p_1 = enreg.predict(feature)
    enreg.load_model('xgb-2.model')
    p_2 = enreg.predict(feature)
    enreg.load_model('xgb-3.model')
    p_3 = enreg.predict(feature)
    enreg.load_model('xgb-4.model')
    p_4 = enreg.predict(feature)
    enreg.load_model('xgb-5.model')
    p_5 = enreg.predict(feature)
    enreg.load_model('xgb-6.model')
    p_6 = enreg.predict(feature)
    enreg.load_model('xgb-7.model')
    p_7 = enreg.predict(feature)
    enreg.load_model('xgb-8.model')
    p_8 = enreg.predict(feature)
    enreg.load_model('xgb-9.model')
    p_9 = enreg.predict(feature)
    enreg.load_model('xgb-10.model')
    p_10 = enreg.predict(feature)
    return (p_1+p_2+p_3+p_4+p_5+p_6+p_7+p_8+p_9+p_10)/10


def calculate_band(a,b,c):
    if element_number == 2:
        formula = elements[0]+str(a)+elements[1]+str(b)
    if element_number == 3:
        formula = elements[0]+str(a)+elements[1]+str(b)+elements[2]+str(c)
    composition = Composition(formula)
    features = np.loadtxt("features.txt").reshape(1,28)
    element_property = np.array(meredig(Composition(composition)))
    BC = BandCenter()
    bandcenter = BC.featurize(composition)
    features[0][10:27] = element_property[0:]
    features[0][27] = bandcenter[0]
    result = model(features)
    return result[0]*(-1)

def band_opt(fraction):
    #fraction = [0.9, 1.0, 0.9, 1.0]
    #element_number = 2
    if element_number == 2:
        opt = BayesianOptimization(calculate_band,\
        {'a':(fraction[0],fraction[1]),\
        'b':(fraction[2],fraction[3])})
        opt.maximize(init_points=200, n_iter=20)
    if element_number == 3:
        opt = BayesianOptimization(calculate_band,\
        {'a':(fraction[0],fraction[1]),\
        'b':(fraction[2],fraction[3]),\
        'c':(fraction[4],fraction[5])})
        opt.maximize(init_points=200, n_iter=20)



if __name__ == '__main__':
    main()
    #band_opt()

