#ifndef PATTERN_4_GPU_OPTIONS_H
#define PATTERN_4_GPU_OPTIONS_H

/*! \brief Définit les implémentations de InitTensor
 */
enum eInitTensorVersion {
  ITV_ori = 0, //! Version CPU d'origine
  ITV_arcgpu_v1 //! Implémentation API GPU Arcane version 1
};

/*! \brief Définit les implémentations de InitNodeVector
 */
enum eInitNodeVectorVersion {
  INVV_ori = 0, //! Version CPU d'origine
  INVV_mt, //! Version CPU multi-thread
  INVV_arcgpu_v1, //! Implémentation API GPU Arcane version 1
  INVV_kokkos  //! Implémentation Kokkos
};

/*! \brief Définit les implémentations de InitNodeCoordBis
 */
enum eInitNodeCoordBisVersion {
  INCBV_ori = 0, //! Version CPU d'origine
  INCBV_arcgpu_v1, //! Implémentation API GPU Arcane version 1
  INCBV_kokkos //! Implémentation Kokkos
};

/*! \brief Définit les implémentations de InitCqs
 */
enum eInitCqsVersion {
  ICQV_ori = 0, //! Version CPU d'origine
  ICQV_arcgpu_v1, //! Implémentation API GPU Arcane version 1
  ICQV_arcgpu_v5, //! Implémentation API GPU Arcane version 5 (GG)
  ICQV_kokkos //! Implémentation Kokkos
};

/*! \brief Définit les implémentations de InitCqs1
 */
enum eInitCqs1Version {
  ICQ1V_ori = 0, //! Version CPU d'origine
  ICQ1V_arcgpu_v1 //! Implémentation API GPU Arcane version 1
};

/*! \brief Définit les implémentations de InitCellArr12
 */
enum eInitCellArr12Version {
  IA12V_ori = 0, //! Version CPU d'origine
  IA12V_arcgpu_v1, //! Implémentation API GPU Arcane version 1
  IA12V_kokkos //! Implémentation Kokkos
};

/*! \brief Définit les implémentations de UpdateVectorFromTensor
 */
enum eUpdateVectorFromTensorVersion {
  UVTV_ori = 0, //! Version CPU d'origine
  UVTV_mt //! Implémentation CPU multi-thread
};

/*! \brief Définit les implémentations de UpdateTensor
 */
enum eUpdateTensorVersion {
  UVV_ori = 0, //! Version CPU d'origine
  UVV_ori_v2, //! Version CPU sans calculs inutiles des valeurs moyennes
  UVV_ori_v3, //! Version CPU avec maj des valeurs moyennes
  UVV_arcgpu_v1, //! Implémentation API GPU Arcane version 1
  UVV_arcgpu_v2a, //! Implémentation API GPU Arcane correspondant à ori_v2 avec tableau intermédiaire
  UVV_arcgpu_v2b, //! Implémentation API GPU Arcane correspondant à ori_v2 sans tableau intermédiaire
  UVV_arcgpu_v3b  //! Implémentation API GPU Arcane correspondant à ori_v3 sans tableau intermédiaire
};

/*! \brief Définit les implémentations de ComputeCqsAndVector
 */
enum eComputeCqsVectorVersion {
  CCVV_ori = 0, //! Version CPU d'origine
  CCVV_mt, //! Implémentation CPU Arcane multi-thread
  CCVV_mt_v2, //! Implémentation CPU Arcane multi-thread version 2
  CCVV_arcgpu_v1, //! Implémentation API GPU Arcane version 1
  CCVV_arcgpu_v2, //! Implémentation API GPU Arcane version 2 (avec NumArray pour CQS)
  CCVV_arcgpu_v5, //! Implémentation API GPU Arcane version 5 (GG)
  CCVV_kokkos //! Implémentation Kokkos
};

/*! \brief Définit les implémentations de ComputeAndPrintError
 */
enum eComputeAndPrintError {
  CPEV_ori = 0, //! Version CPU d'origine
  CPEV_arcgpu_v0, //! Implémentation API GPU Arcane version 0 (sans gestion d'erreur)
  CPEV_arcgpu_v1  //! Implémentation API GPU Arcane version 1
};

/*! \brief Définit les implémentations de InitMEnvVar
 */
enum eInitMEnvVar {
  IMVV_ori = 0, //! Version CPU d'origine
  IMVV_arcgpu_v1 //! Implémentation API GPU Arcane version 1
};

/*! \brief Définit les implémentations de PartialImpureOnly
 */
enum ePartialImpureOnlyVersion {
  PIOV_ori = 0, //! Version CPU d'origine
  PIOV_arcgpu_v1 //! Implémentation API GPU Arcane version 1
};

/*! \brief Définit les implémentations de PartialOnly
 */
enum ePartialOnlyVersion {
  POV_ori = 0, //! Version CPU d'origine
  POV_arcgpu_v1, //! Implémentation API GPU Arcane version 1
  POV_arcgpu_v2, //! Implémentation API GPU Arcane version 2
  POV_arcgpu_v3, //! Implémentation API GPU Arcane version 3
  POV_arcgpu_v4  //! Implémentation API GPU Arcane version 4, draft support multimat
};

/*! \brief Définit les implémentations de PartialAndMean
 */
enum ePartialAndMeanVersion {
  PMV_ori = 0, //! Version CPU d'origine
  PMV_ori_v2, //! Implem Arcane CPU version 2
  PMV_arcgpu_v1, //! Implémentation API GPU Arcane version 1
  PMV_arcgpu_v2  //! Implémentation API GPU Arcane version 2
};

/*! \brief Définit les implémentations de PartialAndMean4
 */
enum ePartialAndMean4Version {
  PM4V_ori = 0, //! Version CPU d'origine
  PM4V_ori_v2, //! Implem Arcane CPU version 2
  PM4V_arcgpu_v1, //! Implémentation API GPU Arcane version 1
  PM4V_arcgpu_v2  //! Implémentation API GPU Arcane version 2
};

/*! \brief Définit les implémentations de PartialAndGlobal5
 */
enum ePartialAndGlobal5Version {
  PG5V_ori = 0,   //! Version CPU d'origine
  PG5V_alter,     //! Version CPU alternative
  PG5V_arcgpu_v1, //! Implémentation API GPU Arcane version 1
  PG5V_arcgpu_v2, //! Implémentation API GPU Arcane version 2, draft support multimat
  PG5V_arcgpu_v3  //! Implémentation API GPU Arcane version 2, draft support multimat, queues async % nb env
};

#endif
