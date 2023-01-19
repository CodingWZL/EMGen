import os
#from matminer.featurizers.structure.sites import SiteStatsFingerprint
from pymatgen import Lattice, Structure, Composition
from StructureFeature import structure_symmetry, radial_distribution, xray_diffraction, average_atomic_volume
from SiteFeature import average_neighbor_num, average_bond_length, average_bond_angle, angular_fourier_series, effective_location
from ComponentFeature import Li_fraction, meredig, band_center
import numpy as np

# Citrination API Key: K7680bck0TbSG5uyRtnzDwtt
# Materials project API Key: uTvDJ6mH2k9aXKQXJ6
# os.environ['CITRINE_KEY'] = ":Q4lUr9HO38ZAJNr21lbljgtt"

cif_number = np.loadtxt("log")

def main():
    #f = open("component.txt", "r")
    index = []
    features = []
    for i in range(1,116753):
        #if i in [69,79,81,362,380,459,573,813,879,1073,\
        #         1099,1225,1532,1536,1537,1544,1550,1565,\
        #         1571,1706,1710,2413,2761,2762,2763,2883,\
        #         2909,2945,2946,3010,3124,3129,3137,3141,\
        #         3144,3145,3147,3155,3159,3164,3167,3175]:
        #    continue
        structure = Structure.from_file("../../MP-database/cif/" + str(int(i)) + ".cif")
        #structure = Structure.from_file("../cif/" + str(i) + ".cif")
        #a = f.readline()
        #print(a)
        ###### Structure Feature ######
        #symmetry_information = structure_symmetry(structure) # dimension:4
        crystal_system, centrosymmetric = structure_symmetry(structure)
        #radial_distribution_information = list(radial_distribution(structure)) # dimension:40
        #xray_diffraction_information = list(xray_diffraction(structure)) # dimension:30
        average_atom_volume, average_volume_atom = average_atomic_volume(structure) # dimension:2

        ###### Site Feature ######
        #neighbor_num = average_neighbor_num(structure) # dimension:1
        #bond_length = average_bond_length(structure) # dimension:1
        #bond_angle = average_bond_angle(structure) # dimension:1
        #average_angular_fourier_series = list(angular_fourier_series(structure)) # dimension: 100
        #effective_ratio = effective_location(structure) # dimension:1

        ###### Component Feature ######
        a = structure.composition
        #print(int(cif_number[i]))
        elemental_property = meredig(Composition(a))
        #try:
        #    site_feature = list(Li_fraction(Composition(a))) # dimension:12
        #except:
        #    index.append(i)
        #    continue
        try:
            bandcenter = band_center(Composition(a)) # dimension:1
        except:
            index.append(i)
            continue
        if np.isnan(elemental_property[5]):
            index.append(i)
            continue
        # 7-element
        feature = crystal_system# + site_feature + bandcenter
        # 1-element
        feature.append(centrosymmetric)
        # 1-element
        feature.append(average_atom_volume)
        # 1-element
        feature.append(average_volume_atom)
        #feature.append(neighbor_num)
        #feature.append(bond_length)
        #feature.append(bond_angle)
        #feature.append(effective_ratio)
        # 17-element 1-element
        feature = feature + elemental_property + bandcenter #site_feature# + bandcenter

        #feature.append(bandcenter)
        features = features + feature
        #print(i)
        #print(i)
        #print(site_feature)
    #f.close()
        print(i)
    features = np.array(features)
    features = features.reshape(-1,28)
    np.savetxt("features-mp.txt", features)
    
    #HSE = np.loadtxt("HSE.txt")
    #b = np.delete(HSE, index)
    #np.savetxt("HSE-new.txt",b)

if __name__ == '__main__':
    main()
