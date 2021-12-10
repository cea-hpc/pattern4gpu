# Pattern4GPU

## Environnement
```
export CC=/chemin/vers/gcc/8.3/bin/gcc 
export CXX=/chemin/vers/gcc/8.3/bin/g++
```

## Configuration cmake
```
module load gnu/8.3.0
module load cmake/3.18.1
mkdir build
cd build
```
En Release [avec profiling nvtx] :
```
cmake -S .. -DArcane_ROOT=/chemin/vers/arcane/blabla/release -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_INIT="-g" -DWANT_CUDA=TRUE -DCMAKE_CUDA_COMPILER=/chemin/vers/nvcc [-DWANT_PROF_ACC=TRUE]
```
En Debug :
```
cmake -S .. -DArcane_ROOT=/chemin/vers/arcane/blabla/dbg -DCMAKE_BUILD_TYPE=Debug -DWANT_CUDA=TRUE -DCMAKE_CUDA_COMPILER=/chemin/vers/nvcc
```

## Compilation
```
module load gnu/8.3.0
cd build
cmake --build . -- -j 4
```
ou bien
```
make VERBOSE=1
```

## Exécution
### Exécution sur CPU
```
module load gnu/8.3.0
cd tests/
cd p4gpu_compute_cqs_vector_tensor/
ccc_mprun -n1 ../../build/src/Pattern4GPU -arcane_opt max_iteration 100 ComputeCqsVectorTensor10x10x10.arc
```

### Exécution sur accéléréateur (après connexion sur noeud GPU)
```
/chemin/vers/build/src/Pattern4GPU -A,AcceleratorRuntime=cuda Test.arc
```

#### Execution avec instrumentation nsys 
[si projet configuré avec `-DWANT_PROF_ACC=TRUE`, prise en compte des points d'entrée] :
```
nsys profile --stats=true --force-overwrite true -o p4gpu /chemin/vers/build/src/Pattern4GPU -A,AcceleratorRuntime=cuda Test.arc
```

#### Execution avec instrumentation ncu 
[si projet configuré avec `-DWANT_PROF_ACC=TRUE`, prise en compte des points d'entrée NVTX] :
```
ncu --nvtx -f --set detailed -o ncu_p4gpu /chemin/vers/build/src/Pattern4GPU -A,AcceleratorRuntime=cuda Test.arc
```

#### Exécution avec instrumentation nvprof (sortie ASCII dans `p4gpu.lognvprof`) 
[si projet configuré avec `-DWANT_PROF_ACC=TRUE`, prise en compte des points d'entrée] :
```
nvprof --print-api-trace --print-gpu-trace --normalized-time-unit col --log-file p4gpu.lognvprof /chemin/vers/build/src/Pattern4GPU -A,AcceleratorRuntime=cuda Test.arc
```

