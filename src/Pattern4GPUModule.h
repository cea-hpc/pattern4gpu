#ifndef PATTERN_4_GPU_MODULE_H
#define PATTERN_4_GPU_MODULE_H

#include <arcane/materials/IMeshMaterialMng.h>
#include <arcane/materials/MeshMaterialVariableRef.h>
#include <arcane/materials/CellToAllEnvCellConverter.h>
#include <arcane/cea/ICartesianMesh.h>
#include "cartesian/ICartesianMesh.h"

#include "Pattern4GPU_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;


class Pattern4GPUModule
: public ArcanePattern4GPUObject
{
 public:
  /** Constructeur de la classe */
  Pattern4GPUModule(const ModuleBuildInfo& mbi);

  /** Destructeur de la classe */
  ~Pattern4GPUModule();
  
 public:

  //! points d'entrée
  void initP4GPU() override; // InitP4GPU
  void initTensor() override; // InitTensor
  void initNodeVector() override; // InitNodeVector
  void initNodeCoordBis() override; // InitNodeCoordBis
  void initCqs() override; // InitCqs
  void initCellArr12() override; // InitCellArr12
  void initCartesian() override; // InitCartesian

  void updateTensor() override; // UpdateTensor
  void updateVectorFromTensor() override; // UpdateVectorFromTensor
  void computeCqsAndVector() override; // ComputeCqsAndVector
  void benchCartesian() override; // benchCartesian

 private:

  void _updateVariable(const MaterialVariableCellReal& volume, MaterialVariableCellReal& f);

  template<Integer DIM>
  void _benchCartesianDim();

 private:

  IMeshMaterialMng* m_mesh_material_mng;
  CellToAllEnvCellConverter* m_allenvcell_converter=nullptr;
  CellGroup m_active_cells;
  MaterialVariableCellReal m_compxx;
  MaterialVariableCellReal m_compxy;
  MaterialVariableCellReal m_compyy;

  // Pour comparer 2 implémentations cartésiennes
  Cartesian::ICartesianMesh* m_cart_cartesian_mesh = nullptr;
  Arcane::ICartesianMesh* m_arc_cartesian_mesh = nullptr;
};

#endif

