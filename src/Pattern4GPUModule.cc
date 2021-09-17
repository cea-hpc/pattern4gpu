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

#include <arcane/AcceleratorRuntimeInitialisationInfo.h>

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
        m_mesh_material_mng, "Compyy", IVariable::PTemporary | IVariable::PExecutionDepend)),
  m_node_index_in_cells(platform::getAcceleratorHostMemoryAllocator())
{
}

Pattern4GPUModule::
~Pattern4GPUModule() {
  delete m_allenvcell_converter;
  delete m_acc_mem_adv;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
accBuild()
{
  PROF_ACC_BEGIN(__FUNCTION__);

  info() << "Using Pattern4GPU with accelerator";
  IApplication* app = subDomain()->application();
  initializeRunner(m_runner,traceMng(),app->acceleratorRuntimeInitialisationInfo());

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* m_node_index_in_cells[nid*N + ième-maille-du-noeud] =                     */
/*                                     indice du noeud nid dans la maille    */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_computeNodeIndexInCells() {
  debug() << "_computeNodeIndexInCells";
  // Un noeud est connecté au maximum à MAX_NODE_CELL mailles
  // Calcul pour chaque noeud son index dans chacune des
  // mailles à laquelle il est connecté.
  NodeGroup nodes = allNodes();
  Integer nb_node = nodes.size();
  m_node_index_in_cells.resize(MAX_NODE_CELL*nb_node);
  m_node_index_in_cells.fill(-1);
  auto node_cell_cty = m_connectivity_view.nodeCell();
  auto cell_node_cty = m_connectivity_view.cellNode();
  ENUMERATE_NODE(inode,nodes){
    NodeLocalId node = *inode;
    Int32 index = 0; 
    Int32 first_pos = node.localId() * MAX_NODE_CELL;
    for( CellLocalId cell : node_cell_cty.cells(node) ){
      Int16 node_index_in_cell = 0; 
      for( NodeLocalId cell_node : cell_node_cty.nodes(cell) ){
        if (cell_node==node)
          break;
        ++node_index_in_cell;
      }    
      m_node_index_in_cells[first_pos + index] = node_index_in_cell;
      ++index;
    }    
  }

  // Tableau préferentiellement lu sur le device
  m_acc_mem_adv->setReadMostly(m_node_index_in_cells.view());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initP4GPU()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initP4GPU";

  // On peut créer maintenant l'objet car la composition des environnements
  // est connue car le le pt d'entree GeomEnv.InitGeomEnv a été appelé
  m_allenvcell_converter=new CellToAllEnvCellConverter(m_mesh_material_mng);

  // On impose un pas de temps (pour l'instant, non paramétrable)
  m_global_deltat = 1.e-3;

  // Pour accélérateur
  m_acc_mem_adv = new AccMemAdviser(options()->getAccMemAdvise());
  m_connectivity_view.setMesh(this->mesh());
  _computeNodeIndexInCells();

  // "Conseils" mémoire
  // CellLocalId
  m_acc_mem_adv->setReadMostly(allCells().view().localIds());
  m_acc_mem_adv->setReadMostly(ownCells().view().localIds());

  // NodeLocalId
  m_acc_mem_adv->setReadMostly(allNodes().view().localIds());
  m_acc_mem_adv->setReadMostly(ownNodes().view().localIds());

  // FaceLocalId
  m_acc_mem_adv->setReadMostly(allFaces().view().localIds());
  m_acc_mem_adv->setReadMostly(ownFaces().view().localIds());

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initTensor()
{
  PROF_ACC_BEGIN(__FUNCTION__);
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
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initNodeVector()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initNodeVector";

  // chaque noeud aura un vecteur de norme 1 mais dans des directions
  // différentes
  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();
  if (options()->getInitNodeVectorVersion() == INVV_ori) 
  {
    ENUMERATE_NODE(node_i, allNodes()) {
      const Real3& c=node_coord[node_i];
      Real cos_th=cos(c.x+c.y+c.z); // garantit une valeur dans [-1,+1]
      Real sin_th=math::sqrt(1-cos_th*cos_th); // garantit une valeur dans [0,+1]
      Real phi=(c.x+1)*(c.y+1)*(c.z+1);
      m_node_vector[node_i]=Real3(sin_th*cos(phi), sin_th*sin(phi),cos_th);
    }
  }
  else if (options()->getInitNodeVectorVersion() == INVV_arcgpu_v1)
  {
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);

    auto in_node_coord = ax::viewIn(command, node_coord);
    auto out_node_vector = ax::viewOut(command, m_node_vector);

    command << RUNCOMMAND_ENUMERATE(Node, nid, allNodes()) {
      const Real3 c=in_node_coord[nid];
      Real cos_th=cos(c.x+c.y+c.z); // garantit une valeur dans [-1,+1]
      Real sin_th=math::sqrt(1-cos_th*cos_th); // garantit une valeur dans [0,+1]
      Real phi=(c.x+1)*(c.y+1)*(c.z+1);
      out_node_vector[nid]=Real3(sin_th*cos(phi), sin_th*sin(phi),cos_th);
    };
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initNodeCoordBis()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initNodeCoordBis";

  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();

  if (options()->getInitNodeCoordBisVersion() == INCBV_ori)
  {
    ENUMERATE_NODE(node_i, allNodes()) {
      m_node_coord_bis[node_i]=node_coord[node_i];
    }
  }
  else if (options()->getInitNodeCoordBisVersion() == INCBV_arcgpu_v1)
  {
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);

    auto in_node_coord = ax::viewIn(command, node_coord);
    auto out_node_coord_bis = ax::viewOut(command, m_node_coord_bis);

    command << RUNCOMMAND_ENUMERATE(Node, nid, allNodes()) {
      out_node_coord_bis[nid]=in_node_coord[nid];
    };
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initCqs()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initCqs";

  // Valable en 3D, 8 noeuds par maille
  m_cell_cqs.resize(8);

  if (options()->getInitCqsVersion() == ICQV_ori)
  { 
    ENUMERATE_CELL (cell_i, allCells()) {
      for(Integer inode(0) ; inode<8 ; ++inode) {
        m_cell_cqs[cell_i][inode] = Real3::zero();
      }
    }
  }
  else if (options()->getInitCqsVersion() == ICQV_arcgpu_v1)
  {
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);

    auto out_cell_cqs = ax::viewOut(command, m_cell_cqs);

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
      for(Integer inode(0) ; inode<8 ; ++inode) {
        out_cell_cqs[cid][inode] = Real3::zero();
      }
    };
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
initCellArr12()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans initCellArr12";

  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();
  if (options()->getInitCellArr12Version() == IA12V_ori) 
  {
    ENUMERATE_CELL(cell_i, allCells()) {
      const Node& first_node=(*cell_i).node(0);
      const Real3& c=node_coord[first_node];
      m_cell_arr1[cell_i]=1.+math::abs(sin(c.x+1)*cos(c.y+1)*sin(c.z+2));
      m_cell_arr2[cell_i]=2.+math::abs(cos(c.x+2)*sin(c.y+1)*cos(c.z+1));
    }
  }
  else if (options()->getInitCellArr12Version() == IA12V_arcgpu_v1)
  {
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);

    auto in_node_coord = ax::viewIn(command, node_coord);

    auto out_cell_arr1 = ax::viewOut(command, m_cell_arr1);
    auto out_cell_arr2 = ax::viewOut(command, m_cell_arr2);

    auto cnc = m_connectivity_view.cellNode();

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
      NodeLocalId first_nid(cnc.nodes(cid)[0]);
      Real3 c=in_node_coord[first_nid];
      out_cell_arr1[cid]=1.+math::abs(sin(c.x+1)*cos(c.y+1)*sin(c.z+2));
      out_cell_arr2[cid]=2.+math::abs(cos(c.x+2)*sin(c.y+1)*cos(c.z+1));
    };
  }
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
_updateVariable(const MaterialVariableCellReal& volume, MaterialVariableCellReal& f)
{
  PROF_ACC_BEGIN(__FUNCTION__);
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
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
updateTensor()
{
  PROF_ACC_BEGIN(__FUNCTION__);
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
  
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
updateVectorFromTensor() {

  PROF_ACC_BEGIN(__FUNCTION__);
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
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Implémentation d'origine de computeCqsAndVector()                         */
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
_computeCqsAndVector_Vori() {

  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans _computeCqsAndVector_Vori";

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

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Implémentation API GPU Arcane version 1                                   */
/*---------------------------------------------------------------------------*/

ARCCORE_HOST_DEVICE inline void computeCQs(Real3 pos[8], Span<Real3> out_cqs) {
  constexpr Real k025 = 0.25;
  Real3 p0 = pos[0];
  Real3 p1 = pos[1];
  Real3 p2 = pos[2];
  Real3 p3 = pos[3];
  Real3 p4 = pos[4];
  Real3 p5 = pos[5];
  Real3 p6 = pos[6];
  Real3 p7 = pos[7];

  out_cqs[0] = -k025*math::vecMul(p4-p3, p1-p3);
  out_cqs[1] = -k025*math::vecMul(p0-p2, p5-p2);
  out_cqs[2] = -k025*math::vecMul(p1-p3, p6-p3);
  out_cqs[3] = -k025*math::vecMul(p7-p2, p0-p2);
  out_cqs[4] = -k025*math::vecMul(p5-p7, p0-p7);
  out_cqs[5] = -k025*math::vecMul(p1-p6, p4-p6);
  out_cqs[6] = -k025*math::vecMul(p5-p2, p7-p2);
  out_cqs[7] = -k025*math::vecMul(p6-p3, p4-p3);
}

void Pattern4GPUModule::
_computeCqsAndVector_Varcgpu_v1() {

  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans _computeCqsAndVector_Varcgpu_v1";

  {
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);

    auto in_node_coord_bis = ax::viewIn(command,m_node_coord_bis);
    auto out_cell_cqs = ax::viewInOut(command,m_cell_cqs);

    auto cnc = m_connectivity_view.cellNode();

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()){
      // Recopie les coordonnées locales (pour le cache)
      Real3 pos[8];
      Int32 index=0;
      for( NodeLocalId nid : cnc.nodes(cid) ){
        pos[index]=in_node_coord_bis[nid];
        ++index;
      }
      // Calcule les résultantes aux sommets
      computeCQs(pos, out_cell_cqs[cid]);
    };
  }
  {
    // On inverse boucle Cell <-> Node car la boucle originelle sur les mailles n'est parallélisable
    // Du coup, on boucle sur les Node
    // On pourrait construire un groupe de noeuds des mailles active_cells et boucler sur ce groupe
    // Mais ici, on a pré-construit un tableau global au mailles qui indique si une maille fait partie
    // du groupe active_cells ou pas (m_is_active_cell)
    // Rem : en décomp. de dom., pour la plupart des sous-dom. on aura : active_cells = allCells
    // Ainsi, en bouclant sur tous les noeuds allNodes(), pour un noeud donné :
    //   1- on initialise le vecteur à 0
    //   2- on calcule les constributions des mailles connectées au noeud uniquement si elles sont actives
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);

    auto in_cell_arr1 = ax::viewIn(command, m_cell_arr1);
    auto in_cell_arr2 = ax::viewIn(command, m_cell_arr2);
    auto in_cell_cqs  = ax::viewIn(command, m_cell_cqs);
    auto in_is_active_cell = ax::viewIn(command, m_is_active_cell);

    auto out_node_vector = ax::viewOut(command, m_node_vector);

    auto node_index_in_cells = m_node_index_in_cells.constSpan();
    auto nc_cty = m_connectivity_view.nodeCell();

    command << RUNCOMMAND_ENUMERATE(Node,nid,allNodes()) {
      Int32 first_pos = nid.localId() * MAX_NODE_CELL;
      Integer index = 0;
      Real3 node_vec = Real3::zero();
      for( CellLocalId cid : nc_cty.cells(nid) ){
        if (in_is_active_cell[cid]) { // la maille ne contribue que si elle est active
          Int16 node_index = node_index_in_cells[first_pos + index];
          node_vec += (in_cell_arr1[cid]+in_cell_arr2[cid]) 
            * in_cell_cqs[cid][node_index];
        }
        ++index;
      }
      out_node_vector[nid] = node_vec;
    };
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
computeCqsAndVector() {

  PROF_ACC_BEGIN(__FUNCTION__);

  switch (options()->getComputeCqsVectorVersion()) {
    case CCVV_ori: _computeCqsAndVector_Vori(); break;
    case CCVV_arcgpu_v1: _computeCqsAndVector_Varcgpu_v1(); break;
  };

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_PATTERN4GPU(Pattern4GPUModule);

