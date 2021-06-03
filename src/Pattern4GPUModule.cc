#include "Pattern4GPUModule.h"

#include <arcane/materials/MeshBlockBuildInfo.h>
#include <arcane/materials/MeshEnvironmentBuildInfo.h>
#include <arcane/materials/IMeshBlock.h>
#include <arcane/materials/MeshMaterialModifier.h>
#include <arcane/materials/MatItemEnumerator.h>
#include <arcane/materials/ComponentPartItemVectorView.h>
#include <arcane/materials/CellToAllEnvCellConverter.h>
#include <arcane/IMesh.h>
#include <arcane/utils/ArcaneGlobal.h>
#include <arcane/utils/StringBuilder.h>

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class IShape {
 public:
  IShape(String name, IMesh* mesh) :
  m_name (name), 
  m_node_coord(mesh->nodesCoordinates()) {
  }
  virtual ~IShape() {}

  const String &name() const { return m_name; }
  virtual void isInside(VariableNodeBool &node_inside, NodeGroup node_group) = 0;

 protected:
  String m_name;
  const VariableNodeReal3& m_node_coord;
};

class ShapeLayer3D : public IShape {
 protected:
  inline bool _isInsidePt(Real x, Real y, Real z) {
    const Real z_min = m_cmin.x * x + m_cmin.y * y + m_cmin.z;
    const Real z_max = m_cmax.x * x + m_cmax.y * y + m_cmax.z;
    return (z_min<=z && z<z_max);
  }
 public:
  ShapeLayer3D(String name, IMesh* mesh, Real3 pmin, Real3 pmax) :
  IShape (name, mesh) {
    m_cmin.z = pmin.z;
    ARCANE_ASSERT(pmin.x != 0., ("pmin.x ne peut pas être nul !"));
    ARCANE_ASSERT(pmin.y != 0., ("pmin.y ne peut pas être nul !"));
    m_cmin.x = -pmin.z/pmin.x;
    m_cmin.y = -pmin.z/pmin.y;

    m_cmax.z = pmax.z;
    ARCANE_ASSERT(pmax.x != 0., ("pmax.x ne peut pas être nul !"));
    ARCANE_ASSERT(pmax.y != 0., ("pmax.y ne peut pas être nul !"));
    m_cmax.x = -pmax.z/pmax.x;
    m_cmax.y = -pmax.z/pmax.y;
  }
  virtual ~ShapeLayer3D() {}

  void isInside(VariableNodeBool &node_inside, NodeGroup node_group) override {
    ENUMERATE_NODE(inode, node_group) {
      const Real3& pt = m_node_coord[inode];
      node_inside[inode] = _isInsidePt(pt.x, pt.y, pt.z);
    }
  }
 protected:
  Real3 m_cmin;
  Real3 m_cmax;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
geomEnvInit()
{
  info() << "Dans geomEnvInit";
  // On va d'abord créer les environnements en les lisant dans le JDD

  m_material_mng = IMeshMaterialMng::getReference(defaultMesh());
  MeshBlockBuildInfo mbbi("BLOCK1",allCells());

  // Définition des différents objets qui vont composer la scene géométrique
  UniqueArray<Real3> l_pts(5); // points qui délimitent les différentes couches
  l_pts[0] = Real3(-0.5,-0.5,-0.5);
  l_pts[1] = Real3(0.5,0.5,0.5);
  l_pts[2] = Real3(0.7,0.7,0.7);
  l_pts[3] = Real3(0.9,0.9,0.9);
  l_pts[4] = Real3(3,3,3);
  Integer nb_sh=l_pts.size()-1;
  UniqueArray<IShape*> l_shape(nb_sh);
  for(Integer ish(0) ; ish<nb_sh ; ++ish) {
    StringBuilder str_build("MIL");
    str_build+=ish;
    l_shape[ish] = new ShapeLayer3D(str_build.toString(), mesh(), l_pts[ish], l_pts[ish+1]);
  }
  // Fin défintion scene géométrique

  // On a un seul matériau par environnement
  for( Integer i=0,n=l_shape.size(); i<n; ++i ){
    String mat_name = l_shape[i]->name()+String("_mat");
    debug() << "Add material name=" << mat_name;
    m_material_mng->registerMaterialInfo(mat_name);
    const String &env_name = l_shape[i]->name();
    MeshEnvironmentBuildInfo env_build(env_name);
    debug() << "Add material=" << mat_name << " in environment=" << env_name;
    env_build.addMaterial(mat_name);
    debug() << "Materiau cree";
    IMeshEnvironment* env = m_material_mng->createEnvironment(env_build);
    debug() << "Environment cree";
    debug() << "Add environment " << env_name << " to BLOCK1";
    mbbi.addEnvironment(env);
  }

  IMeshBlock* block1 = m_material_mng->createBlock(mbbi);

  m_material_mng->endCreate(subDomain()->isContinue());

  MeshMaterialModifier modifier(m_material_mng);
  VariableNodeBool node_inside(VariableBuildInfo(mesh(),"TemporaryNodeInside"));

  // tableau de travail, liste des mailles qui appartiendront à l'environnement
  Int32UniqueArray mat_indexes; 
  mat_indexes.reserve(allCells().size());

  Integer nb_env = block1->nbEnvironment();
  ENUMERATE_ENV(ienv, m_material_mng) {
    IMeshEnvironment* env = *ienv;
    Integer env_id = env->id();

    // On demande à la forme géométrique si tels ou tels noeuds lui appartient
    l_shape[env_id]->isInside(node_inside, allNodes());

    mat_indexes.resize(0); // le tableau est pré-alloué mais on remet sa taille à 0
    ENUMERATE_CELL(icell, allCells()) {
      // Pour chaque maille, si au moins un noeud est dans la forme
      // géométrique, alors la maille appartiendra à l'environnement
      Cell cell = *icell;
      bool is_cell_inside = false;
      ENUMERATE_NODE(inode, cell.nodes()) {
        if (node_inside[inode]) {
          is_cell_inside = true;
        }
      }
      if (is_cell_inside) {
        mat_indexes.add(icell.localId());
      }
    }

    if (!mat_indexes.empty()) {
      // Hypothèse : un SEUL matériau par environnement
      ARCANE_ASSERT(env->nbMaterial() == 1, ("Un environnement ne doit contenir qu'un seul matériau"));
      IMeshMaterial* mat = env->materials()[0];
      modifier.addCells(mat, mat_indexes);
    }
  }
  modifier.endUpdate(); // Pour etre sur que c'est pris en compte pour les statistiques

  for(Integer ish(0) ; ish<l_shape.size() ; ish++) {
    delete l_shape[ish];
  }

  // Statistiques
  Integer nb_tot_cells=allCells().size();
  info() << "Nb total de mailles : " << nb_tot_cells;
  ENUMERATE_ENV(ienv, m_material_mng) {
    IMeshEnvironment* env = *ienv;
    Integer nb_pure_in_env=env->pureEnvItems().nbItem();
    Integer nb_impure_in_env=env->impureEnvItems().nbItem();
    info() << "Dans environnement " << env->name() 
      << ", nb de mailles pures et mixtes : " << nb_pure_in_env << ", " << nb_impure_in_env;
  }

  // nb_cell_env[0] = nb de mailles dont le nb d'env == 0
  // nb_cell_env[1] = nb de mailles dont le nb d'env == 1
  // nb_cell_env[2] = nb de mailles dont le nb d'env >= 2
  UniqueArray<Integer> nb_cell_env(3, /*value=*/0);
  CellToAllEnvCellConverter all_env_cell_converter(m_material_mng);
  ENUMERATE_CELL(icell, allCells()) {
    AllEnvCell all_env_cell = all_env_cell_converter[(*icell)];
    Integer nb_env=all_env_cell.nbEnvironment();
    Integer nb_env_bounded=std::min(nb_env,2);
    nb_cell_env[nb_env_bounded]++;
    m_nbenv[icell]=Real(nb_env);
  }
  info() << "Nb de mailles vides  : " << nb_cell_env[0] << ", ratio=" << nb_cell_env[0]/Real(nb_tot_cells);
  info() << "Nb de mailles pures  : " << nb_cell_env[1] << ", ratio=" << nb_cell_env[1]/Real(nb_tot_cells);
  info() << "Nb de mailles mixtes : " << nb_cell_env[2] << ", ratio=" << nb_cell_env[2]/Real(nb_tot_cells);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
updateTensor()
{
  info() << "Dans updateTensor";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_PATTERN4GPU(Pattern4GPUModule);

