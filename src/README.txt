Organisation des sources :
Pattern4GPUModule.cc         : les points d'entrée sur les tenseurs multi-env, vecteur, cqs ...
Pattern4GPUEnvOrder.cc       : le point d'entrée pour déterminer un ordre de traitement des environnements
geomenv                      : le module d'initialisation des géométries multi-env
Pattern4GPUCartesian.cc      : les points d'entrée pour calculs de volumes par face en cartésien
Pattern4GPUBenchCartesian.cc : plein de boucles élémentaires sur grille cartésienne
cartesian/                   : les classes "à la" Arcane pour proposer des services cartésiens

