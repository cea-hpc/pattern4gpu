#include "Pattern4GPUModule.h"
#include "ViewInDir.h"
#include "cartesian/CartesianMeshProperties.h"
#include "cartesian/CartesianItemSorter.h"

#include "cartesian/FactCartDirectionMng.h"

#include "arcane/VariableView.h"

#include "cartesian/CartTypes.h"
#include "cartesian/CartesianMeshT.h"
#include "arcane/IParallelMng.h"

#include "arcane/cartesianmesh/CellDirectionMng.h"

#define P4GPU_PROFILING // Pour activer le profiling
#include "P4GPUTimer.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*!
 * \brief Initialisation du maillage cartésien
 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
initCartMesh() {
  PROF_ACC_BEGIN(__FUNCTION__);
  Cartesian::CartesianMeshProperties cart_mesh_prop(mesh());
  if (cart_mesh_prop.isPureCartesianMesh()) {
    info() << "Maillage cartésien détecté, tri cartésien des faces";
    Cartesian::CartesianItemSorter cart_sorter(mesh());
    cart_sorter.sortFaces();

    m_cartesian_mesh = CartesianInterface::ICartesianMesh::getReference(mesh(), true);
    m_cartesian_mesh->computeDirections();

    // On force l'implémentation Arcane::
    m_arc_cartesian_mesh = Arcane::arcaneCreateCartesianMesh(mesh());
    m_arc_cartesian_mesh->computeDirections();

  } else {
    info() << "Maillage non cartésien";
  }
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Initialisation du maillage cartésien
 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
initForVol() {
  PROF_ACC_BEGIN(__FUNCTION__);
  const VariableNodeReal3& node_coord = defaultMesh()->nodesCoordinates();

  Cartesian::CartesianMeshProperties cart_mesh_prop(mesh());
  ARCANE_ASSERT(cart_mesh_prop.isPureCartesianMesh(), ("Maillage cartésien obligatoire, ce n'est pas le cas"));
  // Calcul du pas d'espace dans chaque direction
  // D'abord les coordonnées min et max
  Cartesian::FactCartDirectionMng fact_cart_dm(subDomain()->defaultMesh());
  auto* cart_grid=fact_cart_dm.cartesianGrid();
  const auto& cart_num_node=cart_grid->cartNumNode(); // Numérotation cartésienne aux noeuds
  auto loc_first_id=cart_num_node.firstId(); // Id du noeud "en bas à gauche"
  auto loc_last_id=loc_first_id+cart_num_node.nbItem()-1; // Id du noeud "en haut à droite"

  auto node_dm=fact_cart_dm.nodeDirection(0);
  Node first_node=node_dm.toNode(loc_first_id);
  Node last_node=node_dm.toNode(loc_last_id);
  Real3 loc_min_c3=node_coord[first_node]; // plus petites coordonnées du sous-domaine
  Real3 loc_max_c3=node_coord[last_node]; // plus grandes coordonnées du sous-domaine

  Real3 glob_min_c3=defaultMesh()->parallelMng()->reduce(Parallel::ReduceMin, loc_min_c3);
  Real3 glob_max_c3=defaultMesh()->parallelMng()->reduce(Parallel::ReduceMax, loc_max_c3);

  // Puis il faut le nb total de mailles dans chaque direction
  const UniqueIdType3& glob_ncell3=cart_mesh_prop.globalNbCell3();
  // Hypothèse : maillage régulier
  Real3 space_step3=glob_max_c3-glob_min_c3;
  space_step3.x/=Real(glob_ncell3[0]);
  space_step3.y/=Real(glob_ncell3[1]);
  space_step3.z/=Real(glob_ncell3[2]);

  debug() << "Nb de mailles globales dans domaine calcul : (" << glob_ncell3[0] << ", " << glob_ncell3[1] << ", " << glob_ncell3[2] << ")";
  debug() << "Coordonnées min/max : " 
    << " (" << glob_min_c3.x << ", " << glob_min_c3.y << ", " << glob_min_c3.z << "),"
    << " (" << glob_max_c3.x << ", " << glob_max_c3.y << ", " << glob_max_c3.z << ")";
  debug() << "Pas d'espace par direction : " 
    << " (" << space_step3.x << ", " << space_step3.y << ", " << space_step3.z << ")";

  // m_cart_space_step est le pas d'espace
  // Hypothèse : maillage régulier
  m_cart_space_step.fill(space_step3);

  // Récupération des coordonnées des noeuds et construction des coordonnées "déformées"
  // on en profite pour calculer les vitesses aux noeuds
  m_global_deltat=1.e-6; // On impose le pas de temps qu'on relira avec globalDeltaT()
  const Real inv_dt=1./globalDeltaT();

  ENUMERATE_NODE(node_i, allNodes()) {
    const Real3& c=node_coord[node_i];

    // "Déformation" ou déplacement du noeud
    Real3 def3=0.1*space_step3; // maximum 10% d'une maille dans chaque direction
    // calcul d'un vecteur directeur unitaire en coordonnées sphériques
    Real cos_th=cos(10*c.x+15*c.y+13*c.z); // garantit une valeur dans [-1,+1]
    Real sin_th=math::sqrt(1-cos_th*cos_th); // garantit une valeur dans [0,+1]
    Real phi=(c.x+1)*(c.y+1)*(c.z+1);
    def3.x *= sin_th*cos(phi);
    def3.y *= sin_th*sin(phi);
    def3.z *= cos_th;

    // Les variables d'intérêt
    m_car_node_coord[node_i]=c;
    m_def_node_coord[node_i]=c+def3;
    m_node_velocity[node_i]=inv_dt*def3;
  }

  // Ce sont les variables à calculer
  m_dir_trans_area_left.fill(-1.);
  m_face_velocity_left.fill(-1.);
  m_dir_def_coord_left.fill(-1.);
  m_dir_car_coord_left.fill(-1.);
  m_dir_vol1_left.fill(-1.);
  m_dir_vol2_left.fill(-1.);
  m_dir_trans_area_right.fill(-2.);
  m_face_velocity_right.fill(-2.);
  m_dir_def_coord_right.fill(-2.);
  m_dir_car_coord_right.fill(-2.);
  m_dir_vol1_right.fill(-2.);
  m_dir_vol2_right.fill(-2.);
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<typename CartesianMeshT, template<class> class ViewInDirReal >
void Pattern4GPUModule::
_computeVolDir(const Integer dir, const Real dt) {
  PROF_ACC_BEGIN(__FUNCTION__);

  Integer dir_perp_0=(dir+1)%3;
  Integer dir_perp_1=(dir+2)%3;

  using CellDirectionMngType = typename CartesianMeshT::CellDirectionMngType;
  using ConnectivityCellFaceNode = typename CartesianMeshT::ConnectivityCellFaceNode;
  using CellGroupType = typename CartesianMeshT::CellGroupType;

  CartesianMeshT cart_mesh_t(m_cartesian_mesh);

  // Recuperation de toutes les mailles cartesiennes dans la direction dir
  CellDirectionMngType&& cart_cell_dm = cart_mesh_t.cellDirection(dir);
  const CellGroupType&& dirCartCellGroup = cart_cell_dm.allCells();

  ConnectivityCellFaceNode&& cart_conn_cfn = cart_mesh_t.connectivityCellFaceNode(dir);

  // Vues en lecture dans une direction pour des variables vectorielles
  // Attention, des copies peuvent exister en fonctions de l'implem. de la vue
  ViewInDirReal<Cell> in_cart_space_step_dir_perp0(m_cart_space_step, dir_perp_0);
  ViewInDirReal<Cell> in_cart_space_step_dir_perp1(m_cart_space_step, dir_perp_1);
  ViewInDirReal<Node> in_node_velocity_dir(m_node_velocity, dir);
  ViewInDirReal<Node> in_def_node_coord_dir(m_def_node_coord, dir);
  ViewInDirReal<Node> in_car_node_coord_dir(m_car_node_coord, dir);

  // Vues "a la C" de tableaux bien ordonnes
  auto v_dir_trans_area_left      = viewInOut(m_dir_trans_area_left);
  auto v_dir_trans_area_right     = viewInOut(m_dir_trans_area_right);
  auto v_face_velocity_left       = viewInOut(m_face_velocity_left);
  auto v_face_velocity_right      = viewInOut(m_face_velocity_right);
  auto v_dir_def_coord_left       = viewInOut(m_dir_def_coord_left);
  auto v_dir_car_coord_left       = viewInOut(m_dir_car_coord_left);
  auto v_dir_def_coord_right      = viewInOut(m_dir_def_coord_right);
  auto v_dir_car_coord_right      = viewInOut(m_dir_car_coord_right);
  auto v_dir_vol1_left            = viewInOut(m_dir_vol1_left);
  auto v_dir_vol1_right           = viewInOut(m_dir_vol1_right);
  auto v_dir_vol2_left            = viewInOut(m_dir_vol2_left);
  auto v_dir_vol2_right           = viewInOut(m_dir_vol2_right);

  const Integer nb_node_on_face = cart_conn_cfn.nbNode();
  const Real nb_node_inverse = 1.0 / nb_node_on_face;

  // Un calcul elementaire pour une maille pour la direction dir pour un cote donne
  auto lbd_dirvol = [&cart_conn_cfn, dt, nb_node_on_face, nb_node_inverse,
                     &in_cart_space_step_dir_perp0, &in_cart_space_step_dir_perp1,
                     &in_node_velocity_dir, &in_def_node_coord_dir, &in_car_node_coord_dir](
      const CellLocalId &cell_id, 
      const auto &cart_cell_i,
      auto &v_dir_trans_area,
      auto &v_face_velocity,
      auto &v_dir_def_coord,
      auto &v_dir_car_coord,
      auto &v_dir_vol1,
      auto &v_dir_vol2
      ) {

    // calcul de l'aire transversale
    const Real trans_area = in_cart_space_step_dir_perp0[cell_id]
      * in_cart_space_step_dir_perp1[cell_id];

    v_dir_trans_area[cell_id] = trans_area;

    // Passage maille => noeuds sur la face dans la direction
    cart_conn_cfn.initCartCell(cart_cell_i);

    // Modification de l'aire transversale
    Real face_velocity = 0.;
    Real def_coord = 0.;
    Real car_coord = 0.;
    for(Integer inode = 0 ; inode < nb_node_on_face ; inode++) {
      const auto node_id{cart_conn_cfn.node(inode)};

      face_velocity += in_node_velocity_dir[node_id];
      def_coord += in_def_node_coord_dir[node_id];
      car_coord += in_car_node_coord_dir[node_id];
    }
    face_velocity *= nb_node_inverse;
    def_coord *= nb_node_inverse;
    car_coord *= nb_node_inverse;

    v_face_velocity[cell_id] = face_velocity;
    v_dir_def_coord[cell_id] = def_coord;
    v_dir_car_coord[cell_id] = car_coord;

    v_dir_vol1[cell_id] = dt * face_velocity * trans_area;
    v_dir_vol2[cell_id] = v_dir_vol1[cell_id];
  };

  P4GPU_DECLARE_TIMER(subDomain(), Loop_Cell1); P4GPU_START_TIMER(Loop_Cell1);
  // On va d'abord effectuer tous les calculs a gauche, puis on recopiera a droite sauf pour la deniere rangee
  cart_conn_cfn.initSide(MS_previous);

  // Hypothese forte : les mailles doivent etre parcourues de facon cartesienne
  // Pour un plan donne, pour une ligne donnee, parcours des mailles d'une ligne
  ENUMERATE_AUTO_CELL(cell_i, dirCartCellGroup) {

    CellLocalId cell_id(cell_i.localId());

    // Acces "previous" et grandeurs a gauche
    lbd_dirvol(cell_id, cell_i, 
        v_dir_trans_area_left,
        v_face_velocity_left,
        v_dir_def_coord_left,
        v_dir_car_coord_left,
        v_dir_vol1_left,
        v_dir_vol2_left
        );

  }  // Fin ENUMERATE_CELL

  // Maintenant, on recupere les valeurs a droite, sauf pour la derniere rangee ou l'on effectue le calcul
  cart_conn_cfn.initSide(MS_next);

  ENUMERATE_AUTO_CELL(cell_i, dirCartCellGroup) {

    CellLocalId cell_id(cell_i.localId());

    const auto &&dir_cell = cart_cell_dm[cell_i]; // <=> cart_cell_dm.cell(cell_i...)
    CellLocalId next_cell_id(dir_cell.next()); // La maille apres la cellule courante

    if (next_cell_id >= 0) {
      // J'ai une maille a ma droite
      // J'affecte dans MA valeur de droite la valeur de gauche de ma maille de droite
      v_dir_trans_area_right[cell_id] = v_dir_trans_area_left[next_cell_id];
      v_face_velocity_right[cell_id] = v_face_velocity_left[next_cell_id];
      v_dir_def_coord_right[cell_id] = v_dir_def_coord_left[next_cell_id];
      v_dir_car_coord_right[cell_id] = v_dir_car_coord_left[next_cell_id];
      v_dir_vol1_right[cell_id] = v_dir_vol1_left[next_cell_id];
      v_dir_vol2_right[cell_id] = v_dir_vol2_left[next_cell_id];

    } else {
      // Je suis sur la derniere rangee, je calcule
      // Acces "next" et grandeurs a droite
      lbd_dirvol(cell_id, cell_i, 
          v_dir_trans_area_right,
          v_face_velocity_right,
          v_dir_def_coord_right,
          v_dir_car_coord_right,
          v_dir_vol1_right,
          v_dir_vol2_right
          );
    }

  }  // Fin ENUMERATE_CELL

  P4GPU_STOP_TIMER(Loop_Cell1);
  PROF_ACC_END;
}

void Pattern4GPUModule::
_computeVol_Varcgpu_v2() 
{
  PROF_ACC_BEGIN(__FUNCTION__);
  P4GPU_DECLARE_TIMER(subDomain(), Loop_Cell1); P4GPU_START_TIMER(Loop_Cell1);

  Real dt=globalDeltaT();
  auto rqueue = m_acc_env->refQueueAsync();

  for(Integer dir=0 ; dir<mesh()->dimension() ; ++dir) 
  {
    auto cell_dm = m_arc_cartesian_mesh->cellDirection(dir);
    auto fnc = m_acc_env->connectivityView().faceNode();

    Integer dir_perp_0(1), dir_perp_1(2);
    switch (dir) {
      case 0: dir_perp_0 = 1; dir_perp_1 = 2; break;
      case 1: dir_perp_0 = 0; dir_perp_1 = 2; break;
      case 2: dir_perp_0 = 0; dir_perp_1 = 1; break;
    }

    const Integer nb_node_on_face = 1<<(mesh()->dimension()-1);
    const Real nb_node_inverse = 1.0 / nb_node_on_face;

    // Previous (left)
    {
      auto command = makeCommand(rqueue.get());

      // Pour l'instant, on ne se préoccupe pas du rangement des données
      auto in_cart_space_step = ax::viewIn(command, m_cart_space_step);
      auto in_node_velocity   = ax::viewIn(command, m_node_velocity );
      auto in_def_node_coord  = ax::viewIn(command, m_def_node_coord);
      auto in_car_node_coord  = ax::viewIn(command, m_car_node_coord);

      auto out_dir_trans_area_left  = ax::viewOut(command, m_dir_trans_area_left);
      auto out_face_velocity_left   = ax::viewOut(command, m_face_velocity_left);
      auto out_dir_def_coord_left   = ax::viewOut(command, m_dir_def_coord_left);
      auto out_dir_car_coord_left   = ax::viewOut(command, m_dir_car_coord_left);
      auto out_dir_vol1_left        = ax::viewOut(command, m_dir_vol1_left);
      auto out_dir_vol2_left        = ax::viewOut(command, m_dir_vol2_left);

      command.addKernelName("prev") << RUNCOMMAND_ENUMERATE(Cell, cid, cell_dm.allCells()) {

	// calcul de l'aire transversale
	// TODO : utiliser la vue par direction
	const Real trans_area = in_cart_space_step[cid][dir_perp_0]
	  * in_cart_space_step[cid][dir_perp_1];

	out_dir_trans_area_left[cid] = trans_area;

	// Modification de l'aire transversale
	Real face_velocity = 0.;
	Real def_coord = 0.;
	Real car_coord = 0.;

	// On boucle sur les noeuds de la face "previous" transverse
	DirCellFaceLocalId cf(cell_dm.dirCellFaceId(cid));
	FaceLocalId pfid(cf.previousId());
	for(NodeLocalId nid : fnc.nodes(pfid)) {

	  // TODO : utiliser la vue par direction
	  face_velocity += in_node_velocity [nid][dir];
	  def_coord     += in_def_node_coord[nid][dir];
	  car_coord     += in_car_node_coord[nid][dir];
	}
	face_velocity *= nb_node_inverse;
	def_coord *= nb_node_inverse;
	car_coord *= nb_node_inverse;

	out_face_velocity_left[cid] = face_velocity;
	out_dir_def_coord_left[cid] = def_coord;
	out_dir_car_coord_left[cid] = car_coord;

	Real dir_vol1 = dt * face_velocity * trans_area;
	out_dir_vol1_left[cid] = dir_vol1;
	out_dir_vol2_left[cid] = dir_vol1;
      };
    }

    // Next (right)
    {
      auto command = makeCommand(rqueue.get());

      // Pour l'instant, on ne se préoccupe pas du rangement des données
      auto in_cart_space_step = ax::viewIn(command, m_cart_space_step);
      auto in_node_velocity   = ax::viewIn(command, m_node_velocity );
      auto in_def_node_coord  = ax::viewIn(command, m_def_node_coord);
      auto in_car_node_coord  = ax::viewIn(command, m_car_node_coord);

      auto inout_dir_trans_area_left      = ax::viewInOut(command, m_dir_trans_area_left);
      auto inout_dir_trans_area_right     = ax::viewInOut(command, m_dir_trans_area_right);
      auto inout_face_velocity_left       = ax::viewInOut(command, m_face_velocity_left);
      auto inout_face_velocity_right      = ax::viewInOut(command, m_face_velocity_right);
      auto inout_dir_def_coord_left       = ax::viewInOut(command, m_dir_def_coord_left);
      auto inout_dir_car_coord_left       = ax::viewInOut(command, m_dir_car_coord_left);
      auto inout_dir_def_coord_right      = ax::viewInOut(command, m_dir_def_coord_right);
      auto inout_dir_car_coord_right      = ax::viewInOut(command, m_dir_car_coord_right);
      auto inout_dir_vol1_left            = ax::viewInOut(command, m_dir_vol1_left);
      auto inout_dir_vol1_right           = ax::viewInOut(command, m_dir_vol1_right);
      auto inout_dir_vol2_left            = ax::viewInOut(command, m_dir_vol2_left);
      auto inout_dir_vol2_right           = ax::viewInOut(command, m_dir_vol2_right);

      command.addKernelName("next") << RUNCOMMAND_ENUMERATE(Cell, cid, cell_dm.allCells()) {

	DirCellLocalId cc(cell_dm.dirCellId(cid));
	CellLocalId ncid{cc.next()}; // La maille apres la cellule courante

	if (!ncid.isNull()) {
	  // J'ai une maille a ma droite
	  // J'affecte dans MA valeur de droite la valeur de gauche de ma maille de droite
	  inout_dir_trans_area_right[cid] = inout_dir_trans_area_left[ncid];
	  inout_face_velocity_right[cid] = inout_face_velocity_left[ncid];
	  inout_dir_def_coord_right[cid] = inout_dir_def_coord_left[ncid];
	  inout_dir_car_coord_right[cid] = inout_dir_car_coord_left[ncid];
	  inout_dir_vol1_right[cid] = inout_dir_vol1_left[ncid];
	  inout_dir_vol2_right[cid] = inout_dir_vol2_left[ncid];

	} else {
	  // Je suis sur la derniere rangee, je calcule
	  // Acces "next" et grandeurs a droite
	  // GROS COPIER-COLLER EN REMPLACANT LEFT PAR RIGHT ET EN METTANT +1 POUR NEXT
	  //
	  // calcul de l'aire transversale
	  // TODO : utiliser la vue par direction
	  const Real trans_area = in_cart_space_step[cid][dir_perp_0]
	    * in_cart_space_step[cid][dir_perp_1];

	  inout_dir_trans_area_right[cid] = trans_area;

	  // Modification de l'aire transversale
	  Real face_velocity = 0.;
	  Real def_coord = 0.;
	  Real car_coord = 0.;

	  // On boucle sur les noeuds de la face "next" transverse
	  DirCellFaceLocalId cf(cell_dm.dirCellFaceId(cid));
	  FaceLocalId pfid(cf.nextId());
	  for(NodeLocalId nid : fnc.nodes(pfid)) {

	    // TODO : utiliser la vue par direction
	    face_velocity += in_node_velocity [nid][dir];
	    def_coord     += in_def_node_coord[nid][dir];
	    car_coord     += in_car_node_coord[nid][dir];
	  }
	  face_velocity *= nb_node_inverse;
	  def_coord *= nb_node_inverse;
	  car_coord *= nb_node_inverse;

	  inout_face_velocity_right[cid] = face_velocity;
	  inout_dir_def_coord_right[cid] = def_coord;
	  inout_dir_car_coord_right[cid] = car_coord;

	  Real dir_vol1 = dt * face_velocity * trans_area;
	  inout_dir_vol1_right[cid] = dir_vol1;
	  inout_dir_vol2_right[cid] = dir_vol1;
	}
      };
    }
    rqueue->barrier();
  }

  //rqueue->barrier();

  P4GPU_STOP_TIMER(Loop_Cell1);
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Parcours toutes les directions et appele _computeVolDir                   */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
computeVol() {
  PROF_ACC_BEGIN(__FUNCTION__);
  Real dt=globalDeltaT();

  if (options()->getComputeVolVersion() == CVOV_ori)
  {
//#define SOA
#ifdef SOA
    // On recupere les valeurs par direction
#define VIEW_IN_DIR_REAL ViewInDirReal_SoA
#else
    // Vue sur les tableaux en Real3 (AoS)
#define VIEW_IN_DIR_REAL ViewInDirReal_AoS
#endif
    Cartesian::FactCartDirectionMng fact_cart_dm(subDomain()->defaultMesh());
    bool is_cartesian_mesh = fact_cart_dm.isPureCartesianMesh(); 
    //bool is_cartesian_mesh = false; 

    for(Integer dir=0 ; dir<mesh()->dimension() ; ++dir) {
      if (is_cartesian_mesh) {
        _computeVolDir<CartCartesianMeshT, VIEW_IN_DIR_REAL>(dir, dt);
      } else {
        _computeVolDir<UnstructCartesianMeshT, VIEW_IN_DIR_REAL>(dir, dt);
      }
    }
  }
  else if (options()->getComputeVolVersion() == CVOV_arcgpu_v1)
  {
    auto queue = m_acc_env->newQueue();

    using CartesianMeshT = CartCartesianMeshT;
    using ConnectivityCellFaceNode = typename CartesianMeshT::ConnectivityCellFaceNode;

    Cartesian::FactCartDirectionMng cartesian_mesh(mesh());
    CartesianMeshT cart_mesh_t(m_cartesian_mesh);
   
    for(Integer dir=0 ; dir<mesh()->dimension() ; ++dir) 
    {
      auto cell_dm = cartesian_mesh.cellDirection(dir);
      auto c2cid_stm = cell_dm.cell2CellIdStencil();
      auto cell_group = cell_dm.allCells();

      Integer dir_perp_0(1), dir_perp_1(2);
      switch (dir) {
        case 0: dir_perp_0 = 1; dir_perp_1 = 2; break;
        case 1: dir_perp_0 = 0; dir_perp_1 = 2; break;
        case 2: dir_perp_0 = 0; dir_perp_1 = 1; break;
      }

      ConnectivityCellFaceNode&& cart_conn_cfn = cart_mesh_t.connectivityCellFaceNode(dir);

      const Integer nb_node_on_face = cart_conn_cfn.nbNode();
      const Real nb_node_inverse = 1.0 / nb_node_on_face;

      // Previous (left)
      {
        auto command = makeCommand(queue);

        // Pour l'instant, on ne se préoccupe pas du rangement des données
        auto in_cart_space_step = ax::viewIn(command, m_cart_space_step);
        auto in_node_velocity   = ax::viewIn(command, m_node_velocity );
        auto in_def_node_coord  = ax::viewIn(command, m_def_node_coord);
        auto in_car_node_coord  = ax::viewIn(command, m_car_node_coord);

        auto out_dir_trans_area_left  = ax::viewOut(command, m_dir_trans_area_left);
        auto out_face_velocity_left   = ax::viewOut(command, m_face_velocity_left);
        auto out_dir_def_coord_left   = ax::viewOut(command, m_dir_def_coord_left);
        auto out_dir_car_coord_left   = ax::viewOut(command, m_dir_car_coord_left);
        auto out_dir_vol1_left        = ax::viewOut(command, m_dir_vol1_left);
        auto out_dir_vol2_left        = ax::viewOut(command, m_dir_vol2_left);

        // Pour récupérer les noeuds sur la face "previous" orthogonale à dir sur GPU
        auto cfn = cart_conn_cfn.cellFace2Node(MS_previous);

        command.addKernelName("prev") << RUNCOMMAND_LOOP(iter, cell_group.loopRanges()) {
          auto [cid, idx] = c2cid_stm.idIdx(iter);

          // calcul de l'aire transversale
          // TODO : utiliser la vue par direction
          const Real trans_area = in_cart_space_step[cid][dir_perp_0]
            * in_cart_space_step[cid][dir_perp_1];

          out_dir_trans_area_left[cid] = trans_area;

          // Modification de l'aire transversale
          Real face_velocity = 0.;
          Real def_coord = 0.;
          Real car_coord = 0.;

          // On boucle sur les noeuds de la face "previous" transverse
          NodeLocalId bnid{cfn.baseNode(idx)};
          for(Integer inode = 0 ; inode < nb_node_on_face ; inode++) {
            NodeLocalId nid{cfn.node(bnid, inode)};

            // TODO : utiliser la vue par direction
            face_velocity += in_node_velocity [nid][dir];
            def_coord     += in_def_node_coord[nid][dir];
            car_coord     += in_car_node_coord[nid][dir];
          }
          face_velocity *= nb_node_inverse;
          def_coord *= nb_node_inverse;
          car_coord *= nb_node_inverse;

          out_face_velocity_left[cid] = face_velocity;
          out_dir_def_coord_left[cid] = def_coord;
          out_dir_car_coord_left[cid] = car_coord;

          Real dir_vol1 = dt * face_velocity * trans_area;
          out_dir_vol1_left[cid] = dir_vol1;
          out_dir_vol2_left[cid] = dir_vol1;
        };
      }

      // Next (right)
      {
        auto command = makeCommand(queue);

        // Pour l'instant, on ne se préoccupe pas du rangement des données
        auto in_cart_space_step = ax::viewIn(command, m_cart_space_step);
        auto in_node_velocity   = ax::viewIn(command, m_node_velocity );
        auto in_def_node_coord  = ax::viewIn(command, m_def_node_coord);
        auto in_car_node_coord  = ax::viewIn(command, m_car_node_coord);

        auto inout_dir_trans_area_left      = ax::viewInOut(command, m_dir_trans_area_left);
        auto inout_dir_trans_area_right     = ax::viewInOut(command, m_dir_trans_area_right);
        auto inout_face_velocity_left       = ax::viewInOut(command, m_face_velocity_left);
        auto inout_face_velocity_right      = ax::viewInOut(command, m_face_velocity_right);
        auto inout_dir_def_coord_left       = ax::viewInOut(command, m_dir_def_coord_left);
        auto inout_dir_car_coord_left       = ax::viewInOut(command, m_dir_car_coord_left);
        auto inout_dir_def_coord_right      = ax::viewInOut(command, m_dir_def_coord_right);
        auto inout_dir_car_coord_right      = ax::viewInOut(command, m_dir_car_coord_right);
        auto inout_dir_vol1_left            = ax::viewInOut(command, m_dir_vol1_left);
        auto inout_dir_vol1_right           = ax::viewInOut(command, m_dir_vol1_right);
        auto inout_dir_vol2_left            = ax::viewInOut(command, m_dir_vol2_left);
        auto inout_dir_vol2_right           = ax::viewInOut(command, m_dir_vol2_right);

        // Pour récupérer les noeuds sur la face "next" orthogonale à dir sur GPU
        auto cfn = cart_conn_cfn.cellFace2Node(MS_next);

        command.addKernelName("next") << RUNCOMMAND_LOOP(iter, cell_group.loopRanges()) {
          auto [cid, idx] = c2cid_stm.idIdx(iter);

          CellLocalId ncid{c2cid_stm.cell(cid,idx).next()}; // La maille apres la cellule courante

          if (!ItemId::null(ncid)) {
            // J'ai une maille a ma droite
            // J'affecte dans MA valeur de droite la valeur de gauche de ma maille de droite
            inout_dir_trans_area_right[cid] = inout_dir_trans_area_left[ncid];
            inout_face_velocity_right[cid] = inout_face_velocity_left[ncid];
            inout_dir_def_coord_right[cid] = inout_dir_def_coord_left[ncid];
            inout_dir_car_coord_right[cid] = inout_dir_car_coord_left[ncid];
            inout_dir_vol1_right[cid] = inout_dir_vol1_left[ncid];
            inout_dir_vol2_right[cid] = inout_dir_vol2_left[ncid];

          } else {
            // Je suis sur la derniere rangee, je calcule
            // Acces "next" et grandeurs a droite
            // GROS COPIER-COLLER EN REMPLACANT LEFT PAR RIGHT ET EN METTANT +1 POUR NEXT
            //
            // calcul de l'aire transversale
            // TODO : utiliser la vue par direction
            const Real trans_area = in_cart_space_step[cid][dir_perp_0]
              * in_cart_space_step[cid][dir_perp_1];

            inout_dir_trans_area_right[cid] = trans_area;

            // Modification de l'aire transversale
            Real face_velocity = 0.;
            Real def_coord = 0.;
            Real car_coord = 0.;

            // On boucle sur les noeuds de la face "previous" transverse
            // TODO : à encapsuler
            NodeLocalId bnid{cfn.baseNode(idx)};
            for(Integer inode = 0 ; inode < nb_node_on_face ; inode++) {
              NodeLocalId nid{cfn.node(bnid, inode)};

              // TODO : utiliser la vue par direction
              face_velocity += in_node_velocity [nid][dir];
              def_coord     += in_def_node_coord[nid][dir];
              car_coord     += in_car_node_coord[nid][dir];
            }
            face_velocity *= nb_node_inverse;
            def_coord *= nb_node_inverse;
            car_coord *= nb_node_inverse;

            inout_face_velocity_right[cid] = face_velocity;
            inout_dir_def_coord_right[cid] = def_coord;
            inout_dir_car_coord_right[cid] = car_coord;

            Real dir_vol1 = dt * face_velocity * trans_area;
            inout_dir_vol1_right[cid] = dir_vol1;
            inout_dir_vol2_right[cid] = dir_vol1;
          }
        };
      }
    }
  }
  else if (options()->getComputeVolVersion() == CVOV_arcgpu_v2)
  {
    _computeVol_Varcgpu_v2();
  }

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

