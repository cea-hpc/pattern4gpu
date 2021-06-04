#include "Pattern4GPUModule.h"

#include <arcane/geometry/IGeometry.h>
#include <arcane/materials/MeshBlockBuildInfo.h>
#include <arcane/materials/MeshEnvironmentBuildInfo.h>
#include <arcane/materials/IMeshBlock.h>
#include <arcane/materials/MeshMaterialModifier.h>
#include <arcane/materials/MatItemEnumerator.h>
#include <arcane/materials/ComponentPartItemVectorView.h>
#include <arcane/materials/MaterialVariableBuildInfo.h>
#include <arcane/IMesh.h>
#include <arcane/utils/ArcaneGlobal.h>
#include <arcane/utils/StringBuilder.h>

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Pattern4GPUModule::
Pattern4GPUModule(const ModuleBuildInfo& mbi)
  : ArcanePattern4GPUObject(mbi), 
  m_mesh_material_mng(IMeshMaterialMng::getReference(defaultMesh())),
  m_compxx(MaterialVariableBuildInfo(
        m_mesh_material_mng, "Compxx", IVariable::PTemporary | IVariable::PExecutionDepend)),
  m_compxy(MaterialVariableBuildInfo(
        m_mesh_material_mng, "Compxy", IVariable::PTemporary | IVariable::PExecutionDepend)),
  m_compyy(MaterialVariableBuildInfo(
        m_mesh_material_mng, "Compyy", IVariable::PTemporary | IVariable::PExecutionDepend))
{
}

Pattern4GPUModule::
~Pattern4GPUModule() {
  delete m_allenvcell_converter;
}

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
  debug() << "Dans geomEnvInit";
  // On va d'abord créer les environnements en les lisant dans le JDD

  MeshBlockBuildInfo mbbi("BLOCK1",allCells());

  // Définition des différents objets qui vont composer la scene géométrique
  UniqueArray<Real3> l_pts(5); // points qui délimitent les différentes couches
  l_pts[0] = Real3(-0.5,-0.5,-0.5);
  l_pts[1] = Real3(0.5,0.5,0.5);
  l_pts[2] = Real3(0.7,0.7,0.7);
  l_pts[3] = Real3(0.9,0.9,0.9);
  l_pts[4] = Real3(4,4,4);
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
    m_mesh_material_mng->registerMaterialInfo(mat_name);
    const String &env_name = l_shape[i]->name();
    MeshEnvironmentBuildInfo env_build(env_name);
    debug() << "Add material=" << mat_name << " in environment=" << env_name;
    env_build.addMaterial(mat_name);
    debug() << "Materiau cree";
    IMeshEnvironment* env = m_mesh_material_mng->createEnvironment(env_build);
    debug() << "Environment cree";
    debug() << "Add environment " << env_name << " to BLOCK1";
    mbbi.addEnvironment(env);
  }

  IMeshBlock* block1 = m_mesh_material_mng->createBlock(mbbi);

  m_mesh_material_mng->endCreate(subDomain()->isContinue());

  // On précalcule les volumes globaux sur les mailles
  Numerics::IGeometryMng* geom_service=options()->geometry();
  geom_service->init();
  Numerics::IGeometry* geom=geom_service->geometry();
  VariableCellReal cell_volume(VariableBuildInfo(mesh(),"TemporaryCellVolume"));

  ENUMERATE_CELL(icell, allCells()) {
    cell_volume[icell]=geom->computeMeasure(*icell);
  }

  MeshMaterialModifier modifier(m_mesh_material_mng);
  VariableNodeBool node_inside(VariableBuildInfo(mesh(),"TemporaryNodeInside"));

  Integer max_nb_env = block1->nbEnvironment();
  // tableau de travail, liste des mailles qui appartiendront aux environnements
  UniqueArray<Int32UniqueArray> mat_indexes(max_nb_env); 

  // Liste des volumes partiels en cohérence avec mat_indexes
  UniqueArray<RealUniqueArray> partial_volume(max_nb_env);

  ENUMERATE_ENV(ienv, m_mesh_material_mng) {
    IMeshEnvironment* env = *ienv;
    Integer env_id = env->id();

    // On demande à la forme géométrique si tels ou tels noeuds lui appartient
    l_shape[env_id]->isInside(node_inside, allNodes());

    auto& mat_indexes_env=mat_indexes[env_id];
    auto& partial_volume_env=partial_volume[env_id];

    ENUMERATE_CELL(icell, allCells()) {
      // Pour chaque maille, si au moins un noeud est dans la forme
      // géométrique, alors la maille appartiendra à l'environnement
      Cell cell = *icell;
      Integer nb_node_inside=0;
      ENUMERATE_NODE(inode, cell.nodes()) {
        if (node_inside[inode]) {
          nb_node_inside++;
        }
      }
      if (nb_node_inside>0) {
        mat_indexes_env.add(icell.localId());

        // On calcule le volume partiel au prorata du nb de noeuds présents
        // dans la forme géométrique
        Real part_vol_cell = nb_node_inside*cell_volume[icell]/Real(cell.nbNode());
        partial_volume_env.add(part_vol_cell);
      }
    }

    if (!mat_indexes_env.empty()) {
      // Hypothèse : un SEUL matériau par environnement
      ARCANE_ASSERT(env->nbMaterial() == 1, ("Un environnement ne doit contenir qu'un seul matériau"));
      IMeshMaterial* mat = env->materials()[0];
      modifier.addCells(mat, mat_indexes_env);
    }
  }
  modifier.endUpdate(); // Pour etre sur que c'est pris en compte pour les statistiques
  // On peut créer maintenant l'objet car la composition des environnements
  // est connue
  m_allenvcell_converter=new CellToAllEnvCellConverter(m_mesh_material_mng);

  for(Integer ish(0) ; ish<l_shape.size() ; ish++) {
    delete l_shape[ish];
  }

  CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
  // On remplit le tableau des volumes partiels m_volume
  auto internals=defaultMesh()->itemsInternal(IK_Cell).data();
  ENUMERATE_ENV(ienv, m_mesh_material_mng) {
    IMeshEnvironment* env = *ienv;
    Integer env_id = env->id();

    auto& mat_indexes_env=mat_indexes[env_id];
    auto& partial_volume_env=partial_volume[env_id];

    // Un peu compliqué car mat_indexes[env_id] n'est pas forcément trié de la
    // même façon que la liste des EnvCell pour env_id
    for(Integer cptr=0 ; cptr<mat_indexes_env.size() ; ++cptr) {
      Integer lid=mat_indexes_env[cptr];
      Cell cell(internals, lid);
      AllEnvCell all_env_cell = allenvcell_converter[cell];
      EnvCell envcell = env->findEnvCell(all_env_cell);
      ARCANE_ASSERT(mat_indexes_env[cptr]==envcell.globalCell().localId(), ("Incohérence entre les lids des mailles"));
      m_volume[envcell] = partial_volume_env[cptr];
    }
  }
  // On calcule le volume global (et on vérifie qu'il est cohérent avec celui
  // calculé par Arcane)
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell=(*icell);
    AllEnvCell all_env_cell = allenvcell_converter[cell];
    Real vol_sum = 0.;
    ENUMERATE_CELL_ENVCELL (envcell_i, all_env_cell) {
      vol_sum += m_volume[envcell_i];
    }
    Real vol_ref=cell_volume[icell];
    Real ecart=math::abs(vol_sum-vol_ref)/vol_ref;
    ARCANE_ASSERT(ecart<1.e-10, ("Ecart trop important sur volume calculé"));
    m_volume[icell]=vol_sum;
  }

  // Sortie du volume pour la visu
  if (options()->visuVolume()) {
    m_volume_visu.resize(max_nb_env);
    m_volume_visu.fill(0.);
    ENUMERATE_ENV(ienv, m_mesh_material_mng) {
      IMeshEnvironment* env = *ienv;
      Integer env_id = env->id();

      ENUMERATE_ENVCELL (envcell_i, env) {
        Cell cell = (*envcell_i).globalCell();
        m_volume_visu[cell][env_id] = m_volume[envcell_i];
      }
    }
  }

  // Statistiques
  Integer nb_tot_cells=allCells().size();
  info() << "Nb total de mailles : " << nb_tot_cells;
  ENUMERATE_ENV(ienv, m_mesh_material_mng) {
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
  ENUMERATE_CELL(icell, allCells()) {
    AllEnvCell all_env_cell = allenvcell_converter[(*icell)];
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
initTensor()
{
  debug() << "Dans initTensor";

  // Initialise m_tensor de telle sorte qu'il soit symétrique
  ENUMERATE_ENV(ienv, m_mesh_material_mng) {
    IMeshEnvironment* env = *ienv;
    Integer env_id = env->id();
    const Real d=3*env_id;
    ENUMERATE_ENVCELL (envcell_i, env) {
      Real3x3& tens3x3 = m_tensor[envcell_i];
      tens3x3.x.x=d+1.;        tens3x3.x.y=d+1.5;       tens3x3.x.z=d+1.75;
      tens3x3.y.x=tens3x3.x.y; tens3x3.y.y=d+2.;        tens3x3.y.z=d+2.5;
      tens3x3.z.x=tens3x3.x.z; tens3x3.z.y=tens3x3.y.z; tens3x3.z.z=-tens3x3.x.x-tens3x3.y.y;
    }
  }
  CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
  const Real3x3 zero_r3x3(Real3x3::zero());

  ENUMERATE_CELL (cell_i, allCells()) {
    AllEnvCell allenvcell = allenvcell_converter[*cell_i];

    Real3x3& tens3x3 = m_tensor[cell_i];
    Real3x3 tens_sum(zero_r3x3);
    ENUMERATE_CELL_ENVCELL (envcell_i, allenvcell) {
      tens_sum += m_volume[envcell_i]*m_tensor[envcell_i]; // Real3x3 += Real * Real3x3
    }
    tens_sum /= m_volume[cell_i];
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
_updateVariable(const MaterialVariableCellReal& volume, MaterialVariableCellReal& f)
{
  CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
  ENUMERATE_CELL (cell_i, allCells()) {
    AllEnvCell allenvcell = allenvcell_converter[*cell_i];

    if (allenvcell.nbEnvironment() <= 1) continue;

    Real sum = 0.;
    ENUMERATE_CELL_ENVCELL (envcell_i, allenvcell) {
      sum += f[envcell_i];
    }
    f[cell_i] = sum / volume[cell_i];
  }

  ENUMERATE_ENV (env_i, m_mesh_material_mng) {
    ENUMERATE_ENVCELL (envcell_i, *env_i) {
      f[envcell_i] /= volume[envcell_i];
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
updateTensor()
{
  debug() << "Dans updateTensor";

  // Remplir les variables composantes
  ENUMERATE_ENV (env_i, m_mesh_material_mng) {
    ENUMERATE_ENVCELL (envcell_i, *env_i) {
      const Real3x3& real3x3 = m_tensor[envcell_i];
      m_compxx[envcell_i] = real3x3.x.x;
      m_compxy[envcell_i] = real3x3.x.y;
      m_compyy[envcell_i] = real3x3.y.y;
    }
  }

  _updateVariable(m_volume, m_compxx);
  _updateVariable(m_volume, m_compxy);
  _updateVariable(m_volume, m_compyy);

  // Ranger les variables composantes
  ENUMERATE_ENV (env_i, m_mesh_material_mng) {
    ENUMERATE_ENVCELL (envcell_i, *env_i) {
      Real3x3& real3x3 = m_tensor[envcell_i];

      real3x3.x.y = m_compxy[envcell_i];
      real3x3.y.x = m_compxy[envcell_i];

      real3x3.x.x = m_compxx[envcell_i];
      real3x3.y.y = m_compyy[envcell_i];
      real3x3.z.z = - m_compxx[envcell_i] - m_compyy[envcell_i];
    }
  }

  if (defaultMesh()->dimension() == 3) {
    // Remplir les variables composantes (on utilise les memes variables tampon que pour les champs en 2D)
    ENUMERATE_ENV (env_i, m_mesh_material_mng) {
      ENUMERATE_ENVCELL (envcell_i, *env_i) {
        m_compxx[envcell_i] = m_tensor[envcell_i].x.z;
        m_compyy[envcell_i] = m_tensor[envcell_i].y.z;
      }
    }

    // Projeter les variables
    _updateVariable(m_volume, m_compxx);
    _updateVariable(m_volume, m_compyy);

    // Ranger les variables composantes
    ENUMERATE_ENV (env_i, m_mesh_material_mng) {
      ENUMERATE_ENVCELL (envcell_i, *env_i) {
        Real3x3& real3x3 = m_tensor[envcell_i];

        real3x3.x.z = m_compxx[envcell_i];
        real3x3.y.z = m_compyy[envcell_i];
        real3x3.z.x = m_compxx[envcell_i];
        real3x3.z.y = m_compyy[envcell_i];
      }
    }
  }  // end if (dim == 3)
  
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_PATTERN4GPU(Pattern4GPUModule);

