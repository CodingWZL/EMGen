from matminer.featurizers.structure import GlobalSymmetryFeatures
#from matminer.featurizers.structure import RadialDistributionFunction
#from matminer.featurizers.structure import XRDPowderPattern
from matminer.featurizers.structure import SineCoulombMatrix
from pymatgen import Lattice, Structure
#from matminer.featurizers.site import SiteStatsFingerprint
import matplotlib.pyplot as plt
import numpy as np


# Extract the symmetry information for a crystal structure
def structure_symmetry(structure):
    """
    :param structure: Pymatgen Structure
    :return: space_group, crystal_system_int, is_centrosymmetric, n_symmetry_ops
    crystal_system_int:
        "triclinic": 7,
        "monoclinic": 6,
        "orthorhombic": 5,
        "tetragonal": 4,
        "trigonal": 3,
        "hexagonal": 2,
        "cubic": 1,
    """
    GSF = GlobalSymmetryFeatures()
    space_group, crystal_system, crystal_system_int, is_centrosymmetric = GSF.featurize(structure)
    if is_centrosymmetric is True:  # convert the boolean to int
        is_centrosymmetric = 1
    else:
        is_centrosymmetric = 0
    a = [0,0,0,0,0,0,0]
    if crystal_system_int == 1:
        a[0] = 1
    elif crystal_system_int == 2:
        a[1] = 1
    elif crystal_system_int == 3:
        a[2] = 1
    elif crystal_system_int == 4:
        a[3] = 1
    elif crystal_system_int == 5:
        a[4] = 1
    elif crystal_system_int == 6:
        a[5] = 1
    elif crystal_system_int == 7:
        a[6] = 1
    return a,is_centrosymmetric

# Calculate the radial distribution function for a crystal structure
def radial_distribution(structure, show=0):
    """
    :param structure: Pymatgen Structure
    :param show: whether to show the rdf figure
    :return: list of rdf information, default dimension is 40
    """
    RDF = RadialDistributionFunction(cutoff=20.0, bin_size=0.5)  # dimension = cutoff/bin_size
    radial_distribution_information = RDF.featurize(structure)
    if show == 1:
        plt.figure(1)
        x = np.arange(40)
        plt.plot(x, radial_distribution_information)
        plt.show()
    return radial_distribution_information



# Calculate the X-ray diffraction for a crystal structure
def xray_diffraction(structure, show=0):
    """
    :param structure: Pymatgen Structure
    :param show: whether to show the XRD
    :return: list of XRD information, default dimension is patter_length or two_theta_range+1
    """
    XRD = XRDPowderPattern(two_theta_range=(0, 90), bw_method=0.05, pattern_length=30)
    xray_diffraction_information = XRD.featurize(structure)
    if show == 1:
        plt.figure(1)
        x = np.arange(30)
        plt.bar(x, xray_diffraction_information)
        plt.show()
    return xray_diffraction_information



# Calculate the modified coulomb matrix for a crystal structure
def sinecoulombmatrix(X, show=0):
    """
    :param X: a list of Pymatgen Structures
    :param show: whether to show the flatten Coulomb Matrix
    :return: flatten Coulomb Matrix, dimension  = max(N_atom) of X
    """
    SCM = SineCoulombMatrix()
    SCM.fit(X)
    sinecoulombmatrix_information = SCM.featurize(X)
    if show == 1:
        plt.figure(1)
        x = np.arange(32)
        plt.bar(x, sinecoulombmatrix_information[0])
        plt.show()
    return sinecoulombmatrix_information



# Calculate the average atomic volume and average volume atoms for a crystal structure
def average_atomic_volume(structure):
    volume = structure.volume
    number_atoms = structure.num_sites
    average_atom_volume = number_atoms/volume
    average_volume_atom = volume/number_atoms
    return average_atom_volume, average_volume_atom




