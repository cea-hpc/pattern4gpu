#ifndef PATTERN_4_GPU_OPTIONS_H
#define PATTERN_4_GPU_OPTIONS_H

/*! \brief Définit les implémentations de InitNodeVector
 */
enum eInitNodeVectorVersion {
  INVV_ori = 0, //! Version CPU d'origine
  INVV_arcgpu_v1 //! Implémentation API GPU Arcane version 1
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

#endif
