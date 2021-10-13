#ifndef PATTERN_4_GPU_MODULE_H
#define PATTERN_4_GPU_MODULE_H

#include <arcane/materials/IMeshMaterialMng.h>
#include <arcane/materials/MeshMaterialVariableRef.h>
#include <arcane/materials/CellToAllEnvCellConverter.h>
#include <arcane/cea/ICartesianMesh.h>
#include "cartesian/ICartesianMesh.h"
#include "cartesian/interface/ICartesianMesh.h"

// Ajout pour accélérateur
#include "accenv/IAccEnv.h"
#include "accenv/AcceleratorUtils.h"
//

// Ce fichier doit être inclu avant Pattern4GPU_axl.h
#include "Pattern4GPUOptions.h"

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
  void initCqs1() override; // InitCqs1
  void initCellArr12() override; // InitCellArr12

  void initBenchCartesian() override; // InitBenchCartesian

  void initCartMesh() override; // InitCartMesh
  void initForVol() override; // InitForVol

  void initEnvOrder() override; // InitEnvOrder

  void initMEnvVar() override; // InitMEnvVar


  //! points d'entrée "compute-loop"
  void updateTensor() override; // UpdateTensor
  void updateVectorFromTensor() override; // UpdateVectorFromTensor
  void computeCqsAndVector() override; // ComputeCqsAndVector

  void testCartesian() override; // TestCartesian
  void benchCartesian() override; // BenchCartesian

  void computeVol() override; // ComputeVol

  void detEnvOrder() override; // DetEnvOrder

  void partialImpureOnly() override; // PartialImpureOnly
  void partialOnly() override; // PartialOnly
  void partialAndMean() override; // PartialAndMean
  void partialAndMean4() override; // PartialAndMean4

 public:
  // Implémentations des points d'entrées, devrait être private mais 
  // impossible car toute méthode déportée sur GPU doit être publique !
  void _computeCqsAndVector_Vori();
  void _computeCqsAndVector_Vmt();
  void _computeCqsAndVector_Vmt_v2();
  void _computeCqsAndVector_Varcgpu_v1();

  void _testCell2Cell();
  void _testNode2Node();
  void _testFace2Cell();
  void _testCell2Face();
  void _stencilCartesian();

 private:

  void _updateVariable(const MaterialVariableCellReal& volume, MaterialVariableCellReal& f);

  template<Integer DIM>
  void _benchCartesianDim();

  template<typename CartesianMeshT, template<class> class ViewInDirReal>
  void _computeVolDir(const Integer dir, const Real dt);

  template<typename CartesianMeshT>
  void _detEnvOrder();

  // Ecriture m_menv_var1 dans m_menv_var1_visu pour visualisation
  void _dumpVisuMEnvVar();

 private:

  IMeshMaterialMng* m_mesh_material_mng;
  CellToAllEnvCellConverter* m_allenvcell_converter=nullptr;
  CellGroup m_active_cells;
  MaterialVariableCellReal m_compxx;
  MaterialVariableCellReal m_compxy;
  MaterialVariableCellReal m_compyy;
  VariableCellReal m_tmp1; // un tableau temporaire de travail

  // CartesianInterface:: = Arcane:: ou Cartesian::
  CartesianInterface::ICartesianMesh* m_cartesian_mesh = nullptr;

  // Pour comparer 2 implémentations cartésiennes
  Cartesian::ICartesianMesh* m_cart_cartesian_mesh = nullptr;
  Arcane::ICartesianMesh* m_arc_cartesian_mesh = nullptr;

  // Pour l'utilisation des accélérateurs
  IAccEnv* m_acc_env=nullptr;
};

#endif

