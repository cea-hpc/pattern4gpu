Organisation des sources :
Pattern4GPUModule.cc               : les points d'entrée sur les tenseurs multi-env, vecteur, cqs ...
Pattern4GPUComputeAndPrintError.cc : gestion d'une erreur dans une boucle
Pattern4GPUEnvOrder.cc             : le point d'entrée pour déterminer un ordre de traitement des environnements
Pattern4GPUMultiEnv.cc             : patterns caractéristiques sur traitements multi-environnement
geomenv/                           : le module d'initialisation des géométries multi-env
Pattern4GPUCartesian.cc            : les points d'entrée pour calculs de volumes par face en cartésien
Pattern4GPUTestCartesian.cc        : tests sur des acces par stencil par direction
Pattern4GPUBenchCartesian.cc       : plein de boucles élémentaires sur grille cartésienne
cartesian/                         : les classes "à la" Arcane pour proposer des services cartésiens

