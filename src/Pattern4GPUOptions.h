#ifndef PATTERN_4_GPU_OPTIONS_H
#define PATTERN_4_GPU_OPTIONS_H

/*! \brief Définit les implémentations de InitNodeVector
 */
enum eInitNodeVectorVersion {
  INVV_ori = 0, //! Version CPU d'origine
  INVV_mt, //! Version CPU multi-thread
  INVV_arcgpu_v1 //! Implémentation API GPU Arcane version 1
};

/*! \brief Définit les implémentations de InitNodeCoordBis
 */
enum eInitNodeCoordBisVersion {
  INCBV_ori = 0, //! Version CPU d'origine
  INCBV_arcgpu_v1 //! Implémentation API GPU Arcane version 1
};

/*! \brief Définit les implémentations de InitCqs
 */
enum eInitCqsVersion {
  ICQV_ori = 0, //! Version CPU d'origine
  ICQV_arcgpu_v1 //! Implémentation API GPU Arcane version 1
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
  IA12V_arcgpu_v1 //! Implémentation API GPU Arcane version 1
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
  CCVV_arcgpu_v1 //! Implémentation API GPU Arcane version 1
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
  POV_arcgpu_v2  //! Implémentation API GPU Arcane version 2
};

/*! \brief Définit les implémentations de PartialAndMean
 */
enum ePartialAndMeanVersion {
  PMV_ori = 0, //! Version CPU d'origine
  PMV_ori_v2, //! Implem Arcane CPU version 2
  PMV_arcgpu_v1 //! Implémentation API GPU Arcane version 1
};

/*! \brief Définit les implémentations de PartialAndMean4
 */
enum ePartialAndMean4Version {
  PM4V_ori = 0, //! Version CPU d'origine
  PM4V_ori_v2, //! Implem Arcane CPU version 2
  PM4V_arcgpu_v1, //! Implémentation API GPU Arcane version 1
  PM4V_arcgpu_v2  //! Implémentation API GPU Arcane version 2
};

#endif
