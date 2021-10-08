#ifndef PATTERN_4_GPU_OPTIONS_H
#define PATTERN_4_GPU_OPTIONS_H

/*! \brief Définit les implémentations de InitNodeVector
 */
enum eInitNodeVectorVersion {
  INVV_ori = 0, //! Version CPU d'origine
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

/*! \brief Définit les implémentations de InitCellArr12
 */
enum eInitCellArr12Version {
  IA12V_ori = 0, //! Version CPU d'origine
  IA12V_arcgpu_v1 //! Implémentation API GPU Arcane version 1
};

/*! \brief Définit les implémentations de ComputeCqsAndVector
 */
enum eComputeCqsVectorVersion {
  CCVV_ori = 0, //! Version CPU d'origine
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

#endif
