import numpy as np
from pymatgen.core import Composition, Structure, Element
from matminer.featurizers.composition import ElementFraction, BandCenter
from matminer.featurizers.composition import Meredig


# calculate the fraction of Li of a composition

def meredig(composition):
    composition = Composition(composition)
    meredig = Meredig()
    elemental_property = meredig.featurize(composition)
    return elemental_property[-17:]


# calculate the band center by first ionization energy and electron affinity
def band_center(composition):
    composition = Composition(composition)
    BC = BandCenter()
    bandcenter = BC.featurize(composition)
    return bandcenter


def main():
    #structure = Structure.from_file('Li7P3S11_mp-641703_primitive.cif')

    print(meredig("Li3O"))


if __name__ == '__main__':
    main()
