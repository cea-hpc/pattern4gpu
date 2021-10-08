p4gpu_bench_cartesian            :  pleins de petites boucles sur des accès élémentaires cartésiens
p4gpu_cart_vol                   :  calculs par direction de volumes par face en cartésien (1 seul point d'entrée)
p4gpu_compute_cqs                :  calcul de la cqs et du vecteur à partir de la cqs (1 seul pt d'entrée)
p4gpu_compute_cqs_vector_tensor  :  calcul de la cqs et vecteur + utilisation du tenseur pour maj vecteur (2 pts d'entrée)
p4gpu_cqs_vector_tensor          :  calcul cqs et vecteur + maj vecteur par tenseur + maj tenseur en multi-env (3 pts d'entrée)
p4gpu_env_order                  :  détermination de l'ordre des environnements sur les mailles mixtes (1 seul pt d'entrée)
p4gpu_partial_impure_only        :  maj de valeurs partielles sur les mailles mixtes (impures) (1 seul pt d'entrée)
p4gpu_test_cartesian             :  tests sur des acces par stencil directionnel
p4gpu_update_tensor              :  maj tenseur en multi-env (1 seul pt d'entrée)
p4gpu_update_vector_from_tensor  :  utilisation du tenseur pour maj vecteur (1 seul pt d'entrée)
