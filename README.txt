Environnement :
export CC=/chemin/vers/gcc/8.3/bin/gcc 
export CXX=/chemin/vers/gcc/8.3/bin/g++

Configuration cmake :
module load gnu/8.3.0
module load cmake/3.18.1
mkdir build
cd build
En Release :
   cmake -S .. -DArcane_ROOT=/chemin/vers/arcane/blabla/release -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_INIT="-g"
En Debug :
   cmake -S .. -DArcane_ROOT=/chemin/vers/arcane/blabla/dbg -DCMAKE_BUILD_TYPE=Debug 

Compilation :
module load gnu/8.3.0
cd build
cmake --build . -- -j 4
ou bien
make VERBOSE=1

Execution :
module load gnu/8.3.0
cd tests/
cd p4gpu_compute_cqs_vector_tensor/
ccc_mprun -n1 ../../build/src/Pattern4GPU -arcane_opt max_iteration 100 ComputeCqsVectorTensor10x10x10.arc

Organistaion des sources :
src/Pattern4GPUModule.cc         : les points d'entrée sur les tenseurs multi-env, vecteur, cqs ...
src/geomenv                      : le module d'initialisation des géométries multi-env
src/Pattern4GPUCartesian.cc      : les points d'entrée pour calculs de volumes par face en cartésien
src/Pattern4GPUBenchCartesian.cc : plein de boucles élémentaires sur grille cartésienne
src/cartesian/                   : les classes "à la" Arcane pour proposer des services cartésiens

