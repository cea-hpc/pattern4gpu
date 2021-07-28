#ifndef PATTERN_4_GPU_MODULE_H
#define PATTERN_4_GPU_MODULE_H

#include <arcane/materials/IMeshMaterialMng.h>
#include <arcane/materials/MeshMaterialVariableRef.h>
#include <arcane/materials/CellToAllEnvCellConverter.h>
#include <arcane/cea/ICartesianMesh.h>
#include "cartesian/ICartesianMesh.h"
#include "cartesian/interface/ICartesianMesh.h"

// Ajout pour accélérateur
#include "arcane/UnstructuredMeshConnectivity.h"
#include "AcceleratorUtils.h"
//

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

  //! points d'entrée "build"
  void accBuild() override; // AccBuild

  //! points d'entrée "init"
  void initP4GPU() override; // InitP4GPU
  void initTensor() override; // InitTensor
  void initNodeVector() override; // InitNodeVector
  void initNodeCoordBis() override; // InitNodeCoordBis
  void initCqs() override; // InitCqs
  void initCellArr12() override; // InitCellArr12

  void initBenchCartesian() override; // InitBenchCartesian

  void initCartMesh() override; // InitCartMesh
  void initForVol() override; // InitForVol

  void initEnvOrder() override; // InitEnvOrder

  //! points d'entrée "compute-loop"
  void updateTensor() override; // UpdateTensor
  void updateVectorFromTensor() override; // UpdateVectorFromTensor
  void computeCqsAndVector() override; // ComputeCqsAndVector

  void benchCartesian() override; // BenchCartesian

  void computeVol() override; // ComputeVol

  void detEnvOrder() override; // DetEnvOrder

 private:

  void _updateVariable(const MaterialVariableCellReal& volume, MaterialVariableCellReal& f);

  template<Integer DIM>
  void _benchCartesianDim();

  template<typename CartesianMeshT, template<class> class ViewInDirReal>
  void _computeVolDir(const Integer dir, const Real dt);

  template<typename CartesianMeshT>
  void _detEnvOrder();

 private:

  IMeshMaterialMng* m_mesh_material_mng;
  CellToAllEnvCellConverter* m_allenvcell_converter=nullptr;
  CellGroup m_active_cells;
  MaterialVariableCellReal m_compxx;
  MaterialVariableCellReal m_compxy;
  MaterialVariableCellReal m_compyy;

  // CartesianInterface:: = Arcane:: ou Cartesian::
  CartesianInterface::ICartesianMesh* m_cartesian_mesh = nullptr;

  // Pour comparer 2 implémentations cartésiennes
  Cartesian::ICartesianMesh* m_cart_cartesian_mesh = nullptr;
  Arcane::ICartesianMesh* m_arc_cartesian_mesh = nullptr;

  // Pour l'utilisation des accélérateurs
  ax::Runner m_runner;

  UnstructuredMeshConnectivityView m_connectivity_view;
};

#endif

