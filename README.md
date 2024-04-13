## example_MFA-DVR
This repo provides examples of using MFA-DVR pipeline to visualize volumetric data, including encoding raw data to MFA model and visualizing them using MFA-DVR.

###  Dependencies
- C++11 or higher compiler.
- [mfa](https://github.com/adaptive-fam/mfa), Multivariate Functional Approximations (MFA) library.
- [MPI](http://www.mpich.org)

### Build
1. Create project folder
```bash
mkdir myProject
cd myProject
```
2. Get mfa library, no need to build or install it.
```bash
git clone https://github.com/adaptive-fam/mfa
```
3. Get and build example code
```bash
git clone https://github.com/sunjianxin/example_MFA-DVR
cd example_MFA-DVR
mkdir build
cd build
cmake ..  \
-DCMAKE_CXX_COMPILER=mpicxx \
-DMFA_INCLUDE_DIR=path_to_mfa_include_folder \
-DVTK_DIR:PATH=path_to_VTK_DVR-MFA_installation_folder
```
*path_to_mfa_include_folder* is the folder location in the project folder in step 2. *path_to_VTK_DVR-MFA_installation_folder* is the installation location when you configure VTK_DVR-MFA before building, and it is normally at */usr/local/include/vtk-9.0* by default.
### MFA Modeling
Encoding raw volumetric data into MFA model. *fixed* is the program encoding raw volumetric data into MFA model.
* Synthetic data
```bash
cd data
../build/src/fixed/fixed -m 3 -d 4 -f gaussianBeam_8x8x8.xyz -i general -n 8 -v 8 -q 4
mv approx.out gaussianBeam.mfa
```
Above operation encodes raw volumetric data, *.xyz*, into mfa model named *approx.out* using number of control points of 8 and polynomial degree of 4. Detail options of encoding can be checked using *-h* flag.
* Real data
```bash
cd data
../build/src/fixed/fixed -m 3 -d 4 -f vertebra_128x128x128.xyz -i general -n 128 -v 128 -q 4
mv approx.out vertebra.mfa
```
For real data with a large size, a compressed MFA model can be encoded by using a smaller number of control points that is less or equal to the data sample size of each dimension, in another word, using an argument value smaller than 128 for *-v* option for the vertebra dataset. The current MFA encoding library generates the approx.out file including MFA model and metadata, so its size might not be smaller than the raw data. The compression ratio together with other information can be found in the output of the program. Only MFA model will be generated in the upcoming version of mfa library in the future.

### MFA Rendering
Visualize data using MFA-DVR. *simple_ray_cast_switch* is the program generating volume rendering results using MFA-DVR or default NNS or trilinear interpolation methods provided by VTK. Shading can be turn on or off through *SHADEON* macro.
* Synthetic data
```bash
cd build/src/simple_ray_cast_switch
./simple_ray_cast_switch -v ../../../data/gaussianBeam_8x8x8.vtk -m ../../../data/gaussianBeam.mfa -n gaussianbeam -d 8 -e mfa
```
![results](https://github.com/sunjianxin/example_MFA-DVR/blob/main/doc/gaussianBeam.png)
* Real data
```bash
cd build/src/simple_ray_cast_switch
./simple_ray_cast_switch -v ../../../data/vertebra_128x128x128.vtk -m ../../../data/vertebra.mfa -n vertebra -d 128 -e mfa
```
Above operation generate the volume rendering image using MFA-DVR. The *.vtk* file is the VTK data file in structured points format used for comparison. Detail options of volume rendering can be checked using *-h* flag.
![results](https://github.com/sunjianxin/example_MFA-DVR/blob/main/doc/vertebra.png)
