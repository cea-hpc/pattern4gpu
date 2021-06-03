#ifndef PATTERN_4_GPU_MODULE_H
#define PATTERN_4_GPU_MODULE_H

#include "Pattern4GPU_axl.h"

#include <arcane/geometry/IGeometryMng.h>
#include <arcane/materials/IMeshMaterialMng.h>

using namespace Arcane;


class Pattern4GPUModule
: public ArcanePattern4GPUObject
{
 public:
  /** Constructeur de la classe */
  Pattern4GPUModule(const ModuleBuildInfo& mbi)
    : ArcanePattern4GPUObject(mbi) {}
  /** Destructeur de la classe */
  ~Pattern4GPUModule() {}
  
 public:

  //! points d'entrée
  void geomEnvInit() override; // GeomEnvInit
  void updateTensor() override; // UpdateTensor

 private:

  Materials::IMeshMaterialMng* m_material_mng = nullptr;
};

#endif

