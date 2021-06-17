#ifndef GEOM_ENV_MODULE_H
#define GEOM_ENV_MODULE_H

#include <arcane/geometry/IGeometryMng.h>
#include <arcane/materials/IMeshMaterialMng.h>
#include <arcane/materials/MeshMaterialVariableRef.h>

/*! \brief Définit les configurations géométriques pré-définies
 */
enum eGeomScene {
  env5m3 = 0 //! 5 environnements dont des mailles avec 3 environnements
};

#include "geomenv/GeomEnv_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;


class GeomEnvModule
: public ArcaneGeomEnvObject
{
 public:
  /** Constructeur de la classe */
  GeomEnvModule(const ModuleBuildInfo& mbi);

  /** Destructeur de la classe */
  ~GeomEnvModule();
  
 public:

  //! points d'entrée
  void initGeomEnv() override; // InitGeomEnv

 private:

  IMeshMaterialMng* m_mesh_material_mng;
  CellGroup m_active_cells;
};

#endif

