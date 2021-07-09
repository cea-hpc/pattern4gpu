#include "Pattern4GPUModule.h"
#include <arcane/materials/MeshBlockBuildInfo.h>
#include <arcane/materials/MeshEnvironmentBuildInfo.h>
#include <arcane/materials/IMeshBlock.h>
#include <arcane/materials/MeshMaterialModifier.h>
#include <arcane/materials/MatItemEnumerator.h>
#include <arcane/materials/ComponentPartItemVectorView.h>
#include <arcane/materials/MaterialVariableBuildInfo.h>
#include <arcane/IMesh.h>
#include <arcane/IParallelMng.h>
#include <arcane/IItemFamily.h>
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

void Pattern4GPUModule::
initP4GPU()
{
  debug() << "Dans initP4GPU";

  // On peut créer maintenant l'objet car la composition des environnements
  // est connue car le le pt d'entree GeomEnv.InitGeomEnv a été appelé
  m_allenvcell_converter=new CellToAllEnvCellConverter(m_mesh_material_mng);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initTensor()
{
  debug() << "Dans initTensor";

  // Initialise m_tensor de telle sorte qu'il soit symétrique
  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();
  ENUMERATE_ENV(ienv, m_mesh_material_mng) {
    IMeshEnvironment* env = *ienv;
    Integer env_id = env->id();
    const Real d=3*env_id;
    ENUMERATE_ENVCELL (envcell_i, env) {
      EnvCell envcell=(*envcell_i);
      Cell cell=envcell.globalCell();
      const Node& first_node=cell.node(0);
      const Real3& c=node_coord[first_node];
      const Real dd=d+0.5*sin(1+c.x+c.y+c.z);

      Real3x3& tens3x3 = m_tensor[envcell_i];
      tens3x3.x.x=dd+1.;       tens3x3.x.y=dd+1.5;      tens3x3.x.z=dd+1.75;
      tens3x3.y.x=tens3x3.x.y; tens3x3.y.y=dd+2.;       tens3x3.y.z=dd+2.5;
      tens3x3.z.x=tens3x3.x.z; tens3x3.z.y=tens3x3.y.z; tens3x3.z.z=-tens3x3.x.x-tens3x3.y.y;
    }
  }
  CellToAllEnvCellConverter& allenvcell_converter=*m_allenvcell_converter;
  const Real3x3 zero_r3x3(Real3x3::zero());

  ENUMERATE_CELL (cell_i, allCells()) {
    AllEnvCell allenvcell = allenvcell_converter[*cell_i];

    Real3x3& tens3x3 = m_tensor[cell_i];
    tens3x3 = zero_r3x3;
    if (m_volume[cell_i]>0.) {
      ENUMERATE_CELL_ENVCELL (envcell_i, allenvcell) {
        tens3x3 += m_volume[envcell_i]*m_tensor[envcell_i]; // Real3x3 += Real * Real3x3
      }
      tens3x3 /= m_volume[cell_i];
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initNodeVector()
{
  debug() << "Dans initNodeVector";

  // chaque noeud aura un vecteur de norme 1 mais dans des directions
  // différentes
  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();
  ENUMERATE_NODE(node_i, allNodes()) {
    const Real3& c=node_coord[node_i];
    Real cos_th=cos(c.x+c.y+c.z); // garantit une valeur dans [-1,+1]
    Real sin_th=math::sqrt(1-cos_th*cos_th); // garantit une valeur dans [0,+1]
    Real phi=(c.x+1)*(c.y+1)*(c.z+1);
    m_node_vector[node_i]=Real3(sin_th*cos(phi), sin_th*sin(phi),cos_th);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initNodeCoordBis()
{
  debug() << "Dans initNodeCoordBis";

  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();

  ENUMERATE_NODE(node_i, allNodes()) {
    m_node_coord_bis[node_i]=node_coord[node_i];
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initCqs()
{
  debug() << "Dans initCqs";

  // Valable en 3D, 8 noeuds par maille
  m_cell_cqs.resize(8);

  ENUMERATE_CELL (cell_i, allCells()) {
    for(Integer inode(0) ; inode<8 ; ++inode) {
      m_cell_cqs[cell_i][inode] = Real3::zero();
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initCellArr12()
{
  debug() << "Dans initCellArr12";

  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();
  ENUMERATE_CELL(cell_i, allCells()) {
    const Node& first_node=(*cell_i).node(0);
    const Real3& c=node_coord[first_node];
    m_cell_arr1[cell_i]=1.+math::abs(sin(c.x+1)*cos(c.y+1)*sin(c.z+2));
    m_cell_arr2[cell_i]=2.+math::abs(cos(c.x+2)*sin(c.y+1)*cos(c.z+1));
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

void Pattern4GPUModule::
updateVectorFromTensor() {

  debug() << "Dans updateVectorFromTensor";

  Integer nb_blocks = m_mesh_material_mng->blocks().size();

  for (Integer i = 0; i < nb_blocks; i++) {
    // on récupère le groupe de mailles associé au bloc
    IMeshBlock* b = (m_mesh_material_mng->blocks())[i];
    CellGroup cell_group = b->cells();

    Real3x3 cell_tensor;
    ENUMERATE_CELL (cell_i, cell_group) {
    // boucle sur les noeuds de la maille
      ENUMERATE_NODE (node_i, cell_i->nodes()) {
        cell_tensor = m_tensor[cell_i];
        m_node_vector[node_i] -= math::prodTensVec(cell_tensor,
            m_cell_cqs[cell_i][node_i.index()]);
      }
    }
  }  // end iblock loop
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
computeCqsAndVector() {

  debug() << "Dans computeCqsAndVector";

  constexpr Real k025 = 0.25;

  CellGroup active_cells = defaultMesh()->cellFamily()->findGroup("active_cells");

  UniqueArray<Real3> pos(8);
  ENUMERATE_CELL (cell_i,  allCells()) {
    for (Integer ii = 0; ii < 8; ++ii) {
      pos[ii] = m_node_coord_bis[cell_i->node(ii)];
    }

    m_cell_cqs[cell_i][0] = -k025*math::vecMul(pos[4]-pos[3], pos[1]-pos[3]);
    m_cell_cqs[cell_i][1] = -k025*math::vecMul(pos[0]-pos[2], pos[5]-pos[2]);
    m_cell_cqs[cell_i][2] = -k025*math::vecMul(pos[1]-pos[3], pos[6]-pos[3]);
    m_cell_cqs[cell_i][3] = -k025*math::vecMul(pos[7]-pos[2], pos[0]-pos[2]);
    m_cell_cqs[cell_i][4] = -k025*math::vecMul(pos[5]-pos[7], pos[0]-pos[7]);
    m_cell_cqs[cell_i][5] = -k025*math::vecMul(pos[1]-pos[6], pos[4]-pos[6]);
    m_cell_cqs[cell_i][6] = -k025*math::vecMul(pos[5]-pos[2], pos[7]-pos[2]);
    m_cell_cqs[cell_i][7] = -k025*math::vecMul(pos[6]-pos[3], pos[4]-pos[3]);
  }

  ENUMERATE_NODE (node_i, allNodes()) {
    m_node_vector[node_i].assign(0., 0., 0.);
  }

  // Calcul du gradient de pression
  ENUMERATE_CELL (cell_i, active_cells) {
    ENUMERATE_NODE (node_i, cell_i->nodes()) {
      m_node_vector[node_i] += (m_cell_arr1[cell_i] + m_cell_arr2[cell_i]) *
        m_cell_cqs[cell_i][node_i.index()];
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_PATTERN4GPU(Pattern4GPUModule);

