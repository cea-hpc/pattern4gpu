#include "Pattern4GPUModule.h"
#include "cartesian/CartesianMeshProperties.h"
#include "cartesian/CartesianItemSorter.h"

#include "cartesian/CartesianFaceId.h"
#include "cartesian/FactCartDirectionMng.h"
#include "cartesian/CartConnectivityCellNode.h"
#include "cartesian/CellDirectionMng.h"
#include "cartesian/FaceDirectionMng.h"
#include "cartesian/CartesianConnectivity.h"

#include "arcane/VariableView.h"
#include "arcane/cea/CartesianConnectivity.h"
#include "arcane/cea/CellDirectionMng.h"
#include "arcane/cea/FaceDirectionMng.h"
#include <random>

#define P4GPU_PROFILING // Pour activer le profiling
#include "P4GPUTimer.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*!
 * \brief Initialisation des maillages cartésiens et des vraiables pour les
 * autres points d'entrées cartésiens
 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
initCartesian() {

  // Construction des maillages Cartésiens Cartesian:: et Arcane::
  Cartesian::CartesianMeshProperties cart_mesh_prop(mesh());
  if (cart_mesh_prop.isPureCartesianMesh()) {
    info() << "Maillage cartésien détecté, tri cartésien des faces";
    Cartesian::CartesianItemSorter cart_sorter(mesh());
    cart_sorter.sortFaces();

    m_cart_cartesian_mesh = Cartesian::ICartesianMesh::getReference(mesh(), true);
    m_cart_cartesian_mesh->computeDirections();

    // On force l'implémentation Arcane::
    m_arc_cartesian_mesh = Arcane::arcaneCreateCartesianMesh(mesh());
    m_arc_cartesian_mesh->computeDirections();

  } else {
    info() << "Maillage non cartésien";
  }

  // Initialisation des variables aux items
  ENUMERATE_CELL (cell_i, allCells()) {
    m_cell_arr1[cell_i] = 0;
    m_cell_arr2[cell_i] = 0;
  }

  std::minstd_rand0 gen(0); 
  std::uniform_real_distribution<> dis(-1., +1.);  

  ENUMERATE_NODE (node_i, allNodes()) {
    m_node_arr1[node_i] = (dis(gen)>0 ? +1 : -1);
  }
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Comparaison des implementations Cartesian:: et Arcane::
 */
/*---------------------------------------------------------------------------*/
template<Integer DIM>
void Pattern4GPUModule::
_benchCartesianDim() {
  Cartesian::FactCartDirectionMng fact_cart_dm(mesh());
  auto* cart_grid = fact_cart_dm.cartesianGrid();
  auto&& all_cells = fact_cart_dm.cellDirection(0).allCells();
  auto&& all_nodes = fact_cart_dm.nodeDirection(0).allNodes();

  auto in_cell_arr2 = viewIn(m_cell_arr2);
  auto in_cell_arr1 = viewIn(m_cell_arr1);
  auto out_cell_arr1 = viewOut(m_cell_arr1);

  auto in_node_arr1 = viewIn(m_node_arr1);
  auto out_node_arr1 = viewOut(m_node_arr1);

  auto out_face_arr1 = viewOut(m_face_arr1);

  constexpr Integer nrun=100;

  /* Cell => Node
   */
  // Void
  auto ldb_Void_CN = [nrun, &out_cell_arr1](const Integer nb_node, const auto& all_cells) {
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      ENUMERATE_AUTO_CELL (cell_i, all_cells) {
        CellLocalId cell_id(cell_i.localId());
        Real sum=0;
        for(Integer inode(0) ; inode<nb_node ; ++inode) {
          sum += 1;
        }
        out_cell_arr1[cell_id]=sum;
      }
    }
  };
  Integer nb_node = 1 << DIM;
  P4GPU_DECLARE_TIMER(subDomain(), C2N_Arcane____VOID); P4GPU_START_TIMER(C2N_Arcane____VOID);
  ldb_Void_CN(nb_node, allCells());
  P4GPU_STOP_TIMER(C2N_Arcane____VOID);

  P4GPU_DECLARE_TIMER(subDomain(), C2N_Cartesian_VOID); P4GPU_START_TIMER(C2N_Cartesian_VOID);
  ldb_Void_CN(nb_node, all_cells);
  P4GPU_STOP_TIMER(C2N_Cartesian_VOID);


  auto lbd_CC_CellNode = [this, nrun, &out_cell_arr1](const auto& cc) {
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      ENUMERATE_CELL (cell_i, allCells()) {
        Cell&& cell = *cell_i;
        CellLocalId cell_id(cell_i.localId());

        const Node&& node_0 = cc.upperLeft(cell);
        const Node&& node_1 = cc.upperRight(cell);
        const Node&& node_2 = cc.lowerRight(cell);
        const Node&& node_3 = cc.lowerLeft(cell);

        Real sum=0;
        if (!node_0.null()) {
          sum += m_node_arr1[node_0];
        }
        if (!node_1.null()) {
          sum += m_node_arr1[node_1];
        }
        if (!node_2.null()) {
          sum += m_node_arr1[node_2];
        }
        if (!node_3.null()) {
          sum += m_node_arr1[node_3];
        }

        if constexpr (DIM == 3) {
          const Node&& node_4 = cc.topZUpperLeft(cell);
          const Node&& node_5 = cc.topZUpperRight(cell);
          const Node&& node_6 = cc.topZLowerRight(cell);
          const Node&& node_7 = cc.topZLowerLeft(cell);

          if (!node_4.null()) {
            sum += m_node_arr1[node_4];
          }
          if (!node_5.null()) {
            sum += m_node_arr1[node_5];
          }
          if (!node_6.null()) {
            sum += m_node_arr1[node_6];
          }
          if (!node_7.null()) {
            sum += m_node_arr1[node_7];
          }
        }

        out_cell_arr1[cell_id]=sum;
      }
    }
  };

  // Arcane::CartesianConnectivity
  {
    auto arc_cc = m_arc_cartesian_mesh->connectivity();

    P4GPU_DECLARE_TIMER(subDomain(), C2N_Arcane____CartesianConnectivity); P4GPU_START_TIMER(C2N_Arcane____CartesianConnectivity);
    lbd_CC_CellNode(arc_cc);
    P4GPU_STOP_TIMER(C2N_Arcane____CartesianConnectivity);
  }

  // Cartesian::CartesianConnectivity
  {
    Cartesian::CartesianConnectivity cc = m_cart_cartesian_mesh->connectivity();

    P4GPU_DECLARE_TIMER(subDomain(), C2N_Cartesian_CartesianConnectivity); P4GPU_START_TIMER(C2N_Cartesian_CartesianConnectivity);
    lbd_CC_CellNode(cc);
    P4GPU_STOP_TIMER(C2N_Cartesian_CartesianConnectivity);
  }
  
  // ENUMERATE_NODE
  {
    P4GPU_DECLARE_TIMER(subDomain(), C2N_Arcane____ENUMERATE_NODE); P4GPU_START_TIMER(C2N_Arcane____ENUMERATE_NODE);
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      ENUMERATE_CELL (cell_i, allCells()) {
        CellLocalId cell_id(cell_i.localId());
        Real sum=0;
        ENUMERATE_NODE(node_i, cell_i->nodes()) {
          sum+=m_node_arr1[node_i];
        }
        out_cell_arr1[cell_id]=sum;
      }
    }
    P4GPU_STOP_TIMER(C2N_Arcane____ENUMERATE_NODE);
  }
 
  auto lbd_CartConnCN = [nrun, &all_cells, &in_node_arr1, &out_cell_arr1](const auto& conn_cn) {
    Integer nb_node(conn_cn.nbNode());
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      ENUMERATE_AUTO_CELL (cell_i, all_cells) {
        const auto&& cell2node = conn_cn.cellConnectivity(cell_i);
        CellLocalId cell_id(cell_i.localId());
        Real sum=0;
        for(Integer inode(0) ; inode<nb_node ; ++inode) {
          sum += in_node_arr1[cell2node.node(inode)];
        }
        out_cell_arr1[cell_id]=sum;
      }
    }
  };

  // Cartesian::CartConnectivityCellNode tri SN_cart
  {
    Cartesian::CartConnectivityCellNode conn_cn(*cart_grid, Cartesian::CartConnectivityCellNode::SN_cart);

    P4GPU_DECLARE_TIMER(subDomain(), C2N_Cartesian_CartConctvtyCellNode_SCart); P4GPU_START_TIMER(C2N_Cartesian_CartConctvtyCellNode_SCart);
    lbd_CartConnCN(conn_cn);
    P4GPU_STOP_TIMER(C2N_Cartesian_CartConctvtyCellNode_SCart);
  }
  
  // Cartesian::CartConnectivityCellNode tri SN_arc
  {
    Cartesian::CartConnectivityCellNode conn_cn(*cart_grid, Cartesian::CartConnectivityCellNode::SN_arc);

    P4GPU_DECLARE_TIMER(subDomain(), C2N_Cartesian_CartConctvtyCellNode_SArc); P4GPU_START_TIMER(C2N_Cartesian_CartConctvtyCellNode_SArc);
    lbd_CartConnCN(conn_cn);
    P4GPU_STOP_TIMER(C2N_Cartesian_CartConctvtyCellNode_SArc);
  }
  
  // Cartesian::CartConnectivityCellNode tri SN_trigo
  {
    Cartesian::CartConnectivityCellNode conn_cn(*cart_grid, Cartesian::CartConnectivityCellNode::SN_trigo);

    P4GPU_DECLARE_TIMER(subDomain(), C2N_Cartesian_CartConctvtyCellNode_STrigo); P4GPU_START_TIMER(C2N_Cartesian_CartConctvtyCellNode_STrigo);
    lbd_CartConnCN(conn_cn);
    P4GPU_STOP_TIMER(C2N_Cartesian_CartConctvtyCellNode_STrigo);
  }


  /* 
   * Node => Cell
   */
  auto lbd_CC_NodeCell = [this, nrun, &out_node_arr1](const auto& cc) {
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      ENUMERATE_NODE (node_i, allNodes()) {
        Node&& node = *node_i;

        const Cell&& cell_0 = cc.upperLeft(node);
        const Cell&& cell_1 = cc.upperRight(node);
        const Cell&& cell_2 = cc.lowerRight(node);
        const Cell&& cell_3 = cc.lowerLeft(node);

        Real sum=0;
        if (!cell_0.null()) {
          sum += m_cell_arr1[cell_0];
        }
        if (!cell_1.null()) {
          sum += m_cell_arr1[cell_1];
        }
        if (!cell_2.null()) {
          sum += m_cell_arr1[cell_2];
        }
        if (!cell_3.null()) {
          sum += m_cell_arr1[cell_3];
        }

        if constexpr (DIM == 3) {
          const Cell&& cell_4 = cc.topZUpperLeft(node);
          const Cell&& cell_5 = cc.topZUpperRight(node);
          const Cell&& cell_6 = cc.topZLowerRight(node);
          const Cell&& cell_7 = cc.topZLowerLeft(node);

          if (!cell_4.null()) {
            sum += m_cell_arr1[cell_4];
          }
          if (!cell_5.null()) {
            sum += m_cell_arr1[cell_5];
          }
          if (!cell_6.null()) {
            sum += m_cell_arr1[cell_6];
          }
          if (!cell_7.null()) {
            sum += m_cell_arr1[cell_7];
          }
        }

        NodeLocalId node_id(node_i.localId());
        out_node_arr1[node_id]=sum;
      }
    }
  };

  // Arcane::CartesianConnectivity
  {
    auto arc_cc = m_arc_cartesian_mesh->connectivity();

    P4GPU_DECLARE_TIMER(subDomain(), N2C_Arcane____CartesianConnectivity); P4GPU_START_TIMER(N2C_Arcane____CartesianConnectivity);
    lbd_CC_NodeCell(arc_cc);
    P4GPU_STOP_TIMER(N2C_Arcane____CartesianConnectivity);
  }

  // Cartesian::CartesianConnectivity
  {
    Cartesian::CartesianConnectivity cc = m_cart_cartesian_mesh->connectivity();

    P4GPU_DECLARE_TIMER(subDomain(), N2C_Cartesian_CartesianConnectivity); P4GPU_START_TIMER(N2C_Cartesian_CartesianConnectivity);
    lbd_CC_NodeCell(cc);
    P4GPU_STOP_TIMER(N2C_Cartesian_CartesianConnectivity);
  }

  auto lbd_CC_EnumNodeCell = [this, nrun, &out_node_arr1, &in_cell_arr1, &all_nodes](const auto& cc) {
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      ENUMERATE_AUTO_NODE (node_i, all_nodes) {

        const CellLocalId&& cell_0 = cc.upperLeft(node_i);
        const CellLocalId&& cell_1 = cc.upperRight(node_i);
        const CellLocalId&& cell_2 = cc.lowerRight(node_i);
        const CellLocalId&& cell_3 = cc.lowerLeft(node_i);

        Real sum=0;
        if (!ItemId::null(cell_0.localId())) {
          sum += in_cell_arr1[cell_0];
        }
        if (!ItemId::null(cell_1.localId())) {
          sum += in_cell_arr1[cell_1];
        }
        if (!ItemId::null(cell_2.localId())) {
          sum += in_cell_arr1[cell_2];
        }
        if (!ItemId::null(cell_3.localId())) {
          sum += in_cell_arr1[cell_3];
        }

        if constexpr (DIM == 3) {
          const CellLocalId&& cell_4 = cc.topZUpperLeft(node_i);
          const CellLocalId&& cell_5 = cc.topZUpperRight(node_i);
          const CellLocalId&& cell_6 = cc.topZLowerRight(node_i);
          const CellLocalId&& cell_7 = cc.topZLowerLeft(node_i);

          if (!ItemId::null(cell_4.localId())) {
            sum += in_cell_arr1[cell_4];
          }
          if (!ItemId::null(cell_5.localId())) {
            sum += in_cell_arr1[cell_5];
          }
          if (!ItemId::null(cell_6.localId())) {
            sum += in_cell_arr1[cell_6];
          }
          if (!ItemId::null(cell_7.localId())) {
            sum += in_cell_arr1[cell_7];
          }
        }

        NodeLocalId node_id(node_i.localId());
        out_node_arr1[node_id]=sum;
      }
    }
  };

  // Cartesian::CartesianConnectivity avec CartNodeEnumerator
  {
    Cartesian::CartesianConnectivity cc = m_cart_cartesian_mesh->connectivity();

    P4GPU_DECLARE_TIMER(subDomain(), N2C_Cartesian_CartesianConnectivity_CaEmNd); P4GPU_START_TIMER(N2C_Cartesian_CartesianConnectivity_CaEmNd);
    lbd_CC_EnumNodeCell(cc);
    P4GPU_STOP_TIMER(N2C_Cartesian_CartesianConnectivity_CaEmNd);
  }


  auto lbd_CartConnNC = [nrun, &in_cell_arr1, &out_node_arr1, &all_nodes](const auto& conn_nc) {
    Integer max_nb_cell(conn_nc.maxNbCell());
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      ENUMERATE_AUTO_NODE (node_i, all_nodes) {
        const auto&& node2cell = conn_nc.nodeConnectivity(node_i);
        NodeLocalId node_id(node_i.localId());
        Real sum=0;
        for(Integer icell(0) ; icell<max_nb_cell ; ++icell) {
          CellLocalId cell_id(node2cell.cell(icell));
          if (!ItemId::null(cell_id.localId())) {
            sum += in_cell_arr1[cell_id];
          }
        }
        out_node_arr1[node_id]=sum;
      }
    }
  };

  // Cartesian::CartConnectivityNodeCell tri SC_cart
  {
    Cartesian::CartConnectivityNodeCell conn_nc(*cart_grid, Cartesian::CartConnectivityNodeCell::SC_cart);

    P4GPU_DECLARE_TIMER(subDomain(), N2C_Cartesian_CartConctvtyNodeCell_SCart); P4GPU_START_TIMER(N2C_Cartesian_CartConctvtyNodeCell_SCart);
    lbd_CartConnNC(conn_nc);
    P4GPU_STOP_TIMER(N2C_Cartesian_CartConctvtyNodeCell_SCart);
  }

  // Cartesian::CartConnectivityNodeCell tri SC_arc
  {
    Cartesian::CartConnectivityNodeCell conn_nc(*cart_grid, Cartesian::CartConnectivityNodeCell::SC_arc);

    P4GPU_DECLARE_TIMER(subDomain(), N2C_Cartesian_CartConctvtyNodeCell_SArc); P4GPU_START_TIMER(N2C_Cartesian_CartConctvtyNodeCell_SArc);
    lbd_CartConnNC(conn_nc);
    P4GPU_STOP_TIMER(N2C_Cartesian_CartConctvtyNodeCell_SArc);
  }

  // Cartesian::CartConnectivityNodeCell tri SC_trigo
  {
    Cartesian::CartConnectivityNodeCell conn_nc(*cart_grid, Cartesian::CartConnectivityNodeCell::SC_trigo);

    P4GPU_DECLARE_TIMER(subDomain(), N2C_Cartesian_CartConctvtyNodeCell_STrigo); P4GPU_START_TIMER(N2C_Cartesian_CartConctvtyNodeCell_STrigo);
    lbd_CartConnNC(conn_nc);
    P4GPU_STOP_TIMER(N2C_Cartesian_CartConctvtyNodeCell_STrigo);
  }

  /*
   *  INNER NODES
   */

  // Inner nodes dans toutes les directions
  const auto& cart_numb_node = cart_grid->cartNumNode();
  Cartesian::LocalIdType3 inner_beg={0,0,0}, inner_end={1,1,1};
  for(Integer d(0) ; d<DIM ; ++d) {
    inner_beg[d] = 1;
    inner_end[d] = cart_numb_node.nbItemDir(d)-1;
  }
  Cartesian::CartNodeGroup inner_nodes(mesh()->itemsInternal(IK_Node).data(), 0, *cart_grid, cart_numb_node, inner_beg, inner_end);

  auto lbd_InnerCartConnNC = [nrun, &in_cell_arr1, &out_node_arr1, &inner_nodes](const auto& conn_nc) {
    Integer max_nb_cell(conn_nc.maxNbCell());
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      ENUMERATE_AUTO_NODE (node_i, inner_nodes) {
        const auto&& node2cell = conn_nc.innerNodeConnectivity(node_i);
        NodeLocalId node_id(node_i.localId());
        Real sum=0;
        for(Integer icell(0) ; icell<max_nb_cell ; ++icell) {
          sum += in_cell_arr1[node2cell.cell(icell)];
        }
        out_node_arr1[node_id]=sum;
      }
    }
  };

  // Cartesian::CartConnectivityNodeCell tri SC_cart noeuds INTERNES
  {
    Cartesian::CartConnectivityNodeCell conn_nc(*cart_grid, Cartesian::CartConnectivityNodeCell::SC_cart);

    P4GPU_DECLARE_TIMER(subDomain(), N2C_Cartesian_CartConctvtyNodeCell_InNods); P4GPU_START_TIMER(N2C_Cartesian_CartConctvtyNodeCell_InNods);
    lbd_InnerCartConnNC(conn_nc);
    P4GPU_STOP_TIMER(N2C_Cartesian_CartConctvtyNodeCell_InNods);
  }

  /*
   * CellDirectionMng avec des Cells, per Cell previous/next
   */
  auto lbd_CellDirectionMng = [this, nrun](auto cartesian_mesh)
  {
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      for(Integer dir(0) ; dir < DIM ; ++dir) {
        auto&& cell_dm = cartesian_mesh->cellDirection(dir);
        ENUMERATE_CELL(cell_i, cell_dm.allCells()) {
          auto dir_cell = cell_dm[cell_i];

          Cell prev_cell = dir_cell.previous();
          Cell next_cell = dir_cell.next();

          Real sum=0;
          if (!prev_cell.null())
            sum += m_cell_arr2[prev_cell];
          if (!next_cell.null())
            sum += m_cell_arr2[next_cell];

          m_cell_arr1[cell_i]=sum;
        }
      }
    }
  };

  // Arcane::CellDirectionMng 
  {
    P4GPU_DECLARE_TIMER(subDomain(), C2C_Arcane____CellDirectionMng_StdCPN); P4GPU_START_TIMER(C2C_Arcane____CellDirectionMng_StdCPN);
    lbd_CellDirectionMng(m_arc_cartesian_mesh);
    P4GPU_STOP_TIMER(C2C_Arcane____CellDirectionMng_StdCPN);
  }

  /*
   * CellDirectionMng  per Cell previous/next
   */
  auto lbd_CellDirectionMng_CPN = [nrun, &in_cell_arr2, &out_cell_arr1](auto cartesian_mesh)
  {
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      for(Integer dir(0) ; dir < DIM ; ++dir) {
        auto&& cell_dm = cartesian_mesh->cellDirection(dir);
        ENUMERATE_AUTO_CELL(cell_i, cell_dm.allCells()) {
          const auto&& dir_cell = cell_dm[cell_i];

          auto&& prev_cell = dir_cell.previous();
          auto&& next_cell = dir_cell.next();

          CellLocalId prev_cell_id(prev_cell.localId());
          CellLocalId next_cell_id(next_cell.localId());

          Real sum=0;
          if (!ItemId::null(prev_cell_id))
            sum += in_cell_arr2[prev_cell_id];
          if (!ItemId::null(next_cell_id))
            sum += in_cell_arr2[next_cell_id];

          CellLocalId cell_id(cell_i.localId());
          out_cell_arr1[cell_id]=sum;
        }
      }
    }
  };

  // Arcane::CellDirectionMng  per Cell previous/next
  {
    P4GPU_DECLARE_TIMER(subDomain(), C2C_Arcane____CellDirectionMng_CPN); P4GPU_START_TIMER(C2C_Arcane____CellDirectionMng_CPN);
    lbd_CellDirectionMng_CPN(m_arc_cartesian_mesh);
    P4GPU_STOP_TIMER(C2C_Arcane____CellDirectionMng_CPN);
  }

  // Cartesian::CellDirectionMng  per Cell previous/next
  {
    P4GPU_DECLARE_TIMER(subDomain(), C2C_Cartesian_CellDirectionMng_CPN); P4GPU_START_TIMER(C2C_Cartesian_CellDirectionMng_CPN);
    lbd_CellDirectionMng_CPN(m_cart_cartesian_mesh);
    P4GPU_STOP_TIMER(C2C_Cartesian_CellDirectionMng_CPN);
  }

  // Cartesian::CartCellDirectionMng  per Cell previous/next
  {
    P4GPU_DECLARE_TIMER(subDomain(), C2C_Cartesian_CartCellDirectionMng_CPN); P4GPU_START_TIMER(C2C_Cartesian_CartCellDirectionMng_CPN);
    lbd_CellDirectionMng_CPN(&fact_cart_dm);
    P4GPU_STOP_TIMER(C2C_Cartesian_CartCellDirectionMng_CPN);
  }

  /*
   * CellDirectionMng  Cell previous puis Cell next
   */
  auto lbd_CellDirectionMng_CPCN = [nrun, &in_cell_arr2, &in_cell_arr1, &out_cell_arr1](auto cartesian_mesh)
  {
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      for(Integer dir(0) ; dir < DIM ; ++dir) {
        auto&& cell_dm = cartesian_mesh->cellDirection(dir);
        // Tout previous
        ENUMERATE_AUTO_CELL(cell_i, cell_dm.allCells()) {
          const auto&& dir_cell = cell_dm[cell_i];

          auto&& prev_cell = dir_cell.previous();

          CellLocalId prev_cell_id(prev_cell.localId());

          Real sum=0;
          if (!ItemId::null(prev_cell_id))
            sum += in_cell_arr2[prev_cell_id];

          CellLocalId cell_id(cell_i.localId());
          out_cell_arr1[cell_id]=sum;
        }
        // Tout next
        ENUMERATE_AUTO_CELL(cell_i, cell_dm.allCells()) {
          const auto&& dir_cell = cell_dm[cell_i];

          auto&& next_cell = dir_cell.next();

          CellLocalId next_cell_id(next_cell.localId());

          Real sum=0;
          if (!ItemId::null(next_cell_id))
            sum += in_cell_arr2[next_cell_id];

          CellLocalId cell_id(cell_i.localId());
          out_cell_arr1[cell_id]=in_cell_arr1[cell_id]+sum;
        }
      }
    }
  };

  // Arcane::CellDirectionMng  Cell previous puis Cell next
  {
    P4GPU_DECLARE_TIMER(subDomain(), C2C_Arcane____CellDirectionMng_CPCN); P4GPU_START_TIMER(C2C_Arcane____CellDirectionMng_CPCN);
    lbd_CellDirectionMng_CPCN(m_arc_cartesian_mesh);
    P4GPU_STOP_TIMER(C2C_Arcane____CellDirectionMng_CPCN);
  }

  // Cartesian::CellDirectionMng  Cell previous puis Cell next
  {
    P4GPU_DECLARE_TIMER(subDomain(), C2C_Cartesian_CellDirectionMng_CPCN); P4GPU_START_TIMER(C2C_Cartesian_CellDirectionMng_CPCN);
    lbd_CellDirectionMng_CPCN(m_cart_cartesian_mesh);
    P4GPU_STOP_TIMER(C2C_Cartesian_CellDirectionMng_CPCN);
  }

  // Cartesian::CartCellDirectionMng  Cell previous puis Cell next
  {
    P4GPU_DECLARE_TIMER(subDomain(), C2C_Cartesian_CartCellDirectionMng_CPCN); P4GPU_START_TIMER(C2C_Cartesian_CartCellDirectionMng_CPCN);
    lbd_CellDirectionMng_CPCN(&fact_cart_dm);
    P4GPU_STOP_TIMER(C2C_Cartesian_CartCellDirectionMng_CPCN);
  }

  /*
   * Arcane::CellDirectionMng  per Inner Cell previous/next
   */
  auto lbd_CellDirectionMng_AICPN = [nrun, &in_cell_arr2, &out_cell_arr1](auto cartesian_mesh)
  {
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      for(Integer dir(0) ; dir < DIM ; ++dir) {
        auto&& cell_dm = cartesian_mesh->cellDirection(dir);
        // Inner
        ENUMERATE_CELL(cell_i, cell_dm.innerCells()) {
          const auto&& dir_cell = cell_dm[cell_i];

          auto&& prev_cell = dir_cell.previous();
          auto&& next_cell = dir_cell.next();

          CellLocalId prev_cell_id(prev_cell.localId());
          CellLocalId next_cell_id(next_cell.localId());

          Real sum=0;
          sum += in_cell_arr2[prev_cell_id];
          sum += in_cell_arr2[next_cell_id];

          CellLocalId cell_id(cell_i.localId());
          out_cell_arr1[cell_id]=sum;
        }
        // Outer 
        ENUMERATE_AUTO_CELL(cell_i, cell_dm.outerCells()) {
          const auto&& dir_cell = cell_dm[cell_i];

          auto&& prev_cell = dir_cell.previous();
          auto&& next_cell = dir_cell.next();

          CellLocalId prev_cell_id(prev_cell.localId());
          CellLocalId next_cell_id(next_cell.localId());

          Real sum=0;
          if (!prev_cell.null())
            sum += in_cell_arr2[prev_cell_id];
          if (!next_cell.null())
            sum += in_cell_arr2[next_cell_id];

          CellLocalId cell_id(cell_i.localId());
          out_cell_arr1[cell_id]=sum;
        }
      }
    }
  };

  // Arcane::CellDirectionMng  per Inner Cell previous/next
  {
    P4GPU_DECLARE_TIMER(subDomain(), C2C_Arcane____CellDirectionMng_InnerCPN); P4GPU_START_TIMER(C2C_Arcane____CellDirectionMng_InnerCPN);
    lbd_CellDirectionMng_AICPN(m_arc_cartesian_mesh);
    P4GPU_STOP_TIMER(C2C_Arcane____CellDirectionMng_InnerCPN);
  }

  /*
   * CartCellDirectionMng  per Inner Cell previous/next
   */
  auto lbd_CellDirectionMng_ICPN = [nrun, &in_cell_arr2, &out_cell_arr1](auto cartesian_mesh)
  {
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      for(Integer dir(0) ; dir < DIM ; ++dir) {
        auto&& cell_dm = cartesian_mesh->cellDirection(dir);
        // Inner
        ENUMERATE_AUTO_CELL(cell_i, cell_dm.innerCells()) {
//          const auto&& dir_cell = cell_dm[cell_i];
          const auto&& dir_cell = cell_dm.innerCell(cell_i);

          auto&& prev_cell = dir_cell.previous();
          auto&& next_cell = dir_cell.next();

          CellLocalId prev_cell_id(prev_cell.localId());
          CellLocalId next_cell_id(next_cell.localId());

          Real sum=0;
          sum += in_cell_arr2[prev_cell_id];
          sum += in_cell_arr2[next_cell_id];

          CellLocalId cell_id(cell_i.localId());
          out_cell_arr1[cell_id]=sum;
        }
        // Outer gauche
        ENUMERATE_AUTO_CELL(cell_i, cell_dm.previousOuterCells()) {
          const auto&& dir_cell = cell_dm.innerCell(cell_i);

          auto&& next_cell = dir_cell.next();

          CellLocalId next_cell_id(next_cell.localId());
          CellLocalId cell_id(cell_i.localId());
          out_cell_arr1[cell_id]=in_cell_arr2[next_cell_id];
        }
        // Outer droite
        ENUMERATE_AUTO_CELL(cell_i, cell_dm.nextOuterCells()) {
          const auto&& dir_cell = cell_dm.innerCell(cell_i);

          auto&& prev_cell = dir_cell.previous();

          CellLocalId prev_cell_id(prev_cell.localId());
          CellLocalId cell_id(cell_i.localId());
          out_cell_arr1[cell_id]=in_cell_arr2[prev_cell_id];
        }
      }
    }
  };

  // Cartesian::CartCellDirectionMng  per Inner Cell previous/next
  {
    P4GPU_DECLARE_TIMER(subDomain(), C2C_Cartesian_CartCellDirectionMng_InnerCPN); P4GPU_START_TIMER(C2C_Cartesian_CartCellDirectionMng_InnerCPN);
    lbd_CellDirectionMng_ICPN(&fact_cart_dm);
    P4GPU_STOP_TIMER(C2C_Cartesian_CartCellDirectionMng_InnerCPN);
  }

  /*
   * Face => Cell
   */
  auto lbd_FaceDirectionMng_FPN = [this, nrun](auto cartesian_mesh)
  {
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      for(Integer dir(0) ; dir < DIM ; ++dir) {
        auto&& face_dm = cartesian_mesh->faceDirection(dir);
        ENUMERATE_FACE(face_i, face_dm.innerFaces()) {

          auto&& cart_face = face_dm[face_i];
          Cell&& prev_cell = cart_face.previousCell();
          Cell&& next_cell = cart_face.nextCell();

          m_face_arr1[face_i] = m_cell_arr2[prev_cell] + m_cell_arr2[next_cell];
        }  // Fin ENUMERATE_FACE
      }
    }
  };

  // Arcane::FaceDirectionMng
  {
    P4GPU_DECLARE_TIMER(subDomain(), F2C_Arcane____FaceDirectionMng_InnerFPN); P4GPU_START_TIMER(F2C_Arcane____FaceDirectionMng_InnerFPN);
    lbd_FaceDirectionMng_FPN(m_arc_cartesian_mesh);
    P4GPU_STOP_TIMER(F2C_Arcane____FaceDirectionMng_InnerFPN);
  }

  // Cartesian::FaceDirectionMng
  {
    P4GPU_DECLARE_TIMER(subDomain(), F2C_Cartesian_FaceDirectionMng_InnerFPN); P4GPU_START_TIMER(F2C_Cartesian_FaceDirectionMng_InnerFPN);
    lbd_FaceDirectionMng_FPN(m_cart_cartesian_mesh);
    P4GPU_STOP_TIMER(F2C_Cartesian_FaceDirectionMng_InnerFPN);
  }

  auto lbd_CartFaceDirectionMng_FPN = [nrun, &in_cell_arr2, &out_face_arr1](auto cartesian_mesh)
  {
    for(Integer irun(0) ; irun<nrun ; ++irun) {
      for(Integer dir(0) ; dir < DIM ; ++dir) {
        auto&& face_dm = cartesian_mesh->faceDirection(dir);
        ENUMERATE_AUTO_FACE(face_i, face_dm.innerFaces()) {

          auto&& cart_face = face_dm[face_i];
          CellLocalId prev_cell_id(cart_face.previousCell());
          CellLocalId next_cell_id(cart_face.nextCell());

          FaceLocalId face_id(face_i.localId());
          out_face_arr1[face_id] = in_cell_arr2[prev_cell_id] + in_cell_arr2[next_cell_id];
        }  // Fin ENUMERATE_FACE
      }
    }
  };

  // Cartesian::CartFaceDirectionMng
  {
    P4GPU_DECLARE_TIMER(subDomain(), F2C_Cartesian_CartFaceDirectionMng_InnerFPN); P4GPU_START_TIMER(F2C_Cartesian_CartFaceDirectionMng_InnerFPN);
    lbd_CartFaceDirectionMng_FPN(&fact_cart_dm);
    P4GPU_STOP_TIMER(F2C_Cartesian_CartFaceDirectionMng_InnerFPN);
  }

}

void Pattern4GPUModule::
benchCartesian() {
  info() << "------------------------------------------------";
  info() << "Comparaison des implementations Cartesian:: et";
  info() << "                                   Arcane:: ";
  info() << "------------------------------------------------";
  if (mesh()->dimension() == 2) {
    _benchCartesianDim<2>();
  } else if (mesh()->dimension() == 3) {
    _benchCartesianDim<3>();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

