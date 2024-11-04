# EMGen

The source code provides the feature calculations, model training, and directed design.
![image](https://github.com/CodingWZL/ElecML/assets/104205506/59a4c1ab-904e-4eec-9feb-7b30f4134edd)

# environmental.py
which provides the dependencies for installing and running EMGen.
One can use the "conda env create -f environmental.py -n you_name" to create the running environment.

# calculate.py
which is used for calculating the descriptors for crystal structures in current directory. Users should provide at least one .cif file in current directory, and set the filename as parameters in this program, and then run this program by "python calculate.py"

# ComponentFeature.py
which is used for calculating the component features for crystal structures in current directory, and is called by calculate.py.

# StructureFeature.py
which is used for calculating the structural features for crystal structures in current directory, and is called by calculate.py.

# opt.py
which is used to directed design semiconductors by applying active learning methods, where the well-trained models are served as the surrogate models.

# Reference:
Zhilong Wang, Sixian Liu, Jing Gao, Kehao Tao, An Chen, Hongxiao Duan, Yanqiang Han, Fengqi You*, Gang Liu*, Jinjin Li*. Interpretable Surrogate Learning for Electronic Materials Generation. ACS Nano Accepted (2024) https://doi.org/10.1021/acsnano.4c12166.
