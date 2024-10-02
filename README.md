## Functional Approximation Encoder
This repo provides examples of using MFA, Multivariate Functional Approximation, to encode 3D volume into continuous model, called micro-model, proposed in the paper of adaptive function approximation multi-resolutio (Adaptive-FAM).

###  Dependencies
- C++11 or higher compiler.
- [mfa](https://github.com/adaptive-fam/mfa), Multivariate Functional Approximations (MFA) library.
- [MPI](http://www.mpich.org), Message Passing Interface

### Build
1. Create project folder
```bash
mkdir project
cd project
```
2. Get mfa library, no need to build or install it.
```bash
git clone https://github.com/sunjianxin/Adaptive-FAM_mfa.git
```
3. Get and build code
```bash
git clone https://github.com/adaptive-fam/mfa_utility
cd mfa_utility
mkdir build
cd build
cmake ..  \
-DCMAKE_CXX_COMPILER=mpicxx \
-DMFA_INCLUDE_DIR=path_to_mfa_include_folder
make
```
*path_to_mfa_include_folder* is the folder location in the project folder in step 2.
### MFA Modeling
Encoding raw volumetric data into MFA model. *fixed* is the program encoding raw volumetric data into MFA model.
```bash
./build/src/fixed/fixed -m 3 -d 4 -i sinc -n 8 -v 8 -q 2
```
Above operation encodes 3D sinc volume data into micro-model named test.mfab. using number of control points of 8 and polynomial degree of 2. Detail options of encoding can be checked using *-h* flag.
