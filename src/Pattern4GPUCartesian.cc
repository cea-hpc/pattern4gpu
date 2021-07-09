#include "Pattern4GPUModule.h"
#include "ViewInDir.h"
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

#include "cartesian/CartTypes.h"
#include "cartesian/CartesianMeshT.h"
#include "arcane/IParallelMng.h"

#define P4GPU_PROFILING // Pour activer le profiling
#include "P4GPUTimer.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*!
 * \brief Initialisation des maillages cartésiens et des variables pour le
 * point d'entrée BenchCartesian
 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
initBenchCartesian() {

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

/*---------------------------------------------------------------------------*/
/*!
 * \brief Initialisation du maillage cartésien
 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
initCartMesh() {
  Cartesian::CartesianMeshProperties cart_mesh_prop(mesh());
  if (cart_mesh_prop.isPureCartesianMesh()) {
    info() << "Maillage cartésien détecté, tri cartésien des faces";
    Cartesian::CartesianItemSorter cart_sorter(mesh());
    cart_sorter.sortFaces();

    m_cartesian_mesh = CartesianInterface::ICartesianMesh::getReference(mesh(), true);
    m_cartesian_mesh->computeDirections();

  } else {
    info() << "Maillage non cartésien";
  }
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Initialisation du maillage cartésien
 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
initForVol() {
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
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<typename CartesianMeshT, template<class> class ViewInDirReal >
void Pattern4GPUModule::
_computeVolDir(const Integer dir, const Real dt) {

  Integer dir_perp_0=(dir+1)%3;
  Integer dir_perp_1=(dir+2)%3;

  using CellDirectionMngType = typename CartesianMeshT::CellDirectionMngType;
  using FaceDirectionMngType = typename CartesianMeshT::FaceDirectionMngType;
  using ConnectivityCellNode = typename CartesianMeshT::ConnectivityCellNode;
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
  auto v_dir_trans_area_left      = viewOut(m_dir_trans_area_left);
  auto v_dir_trans_area_right     = viewOut(m_dir_trans_area_right);
  auto v_face_velocity_left       = viewOut(m_face_velocity_left);
  auto v_face_velocity_right      = viewOut(m_face_velocity_right);
  auto v_dir_def_coord_left       = viewOut(m_dir_def_coord_left);
  auto v_dir_car_coord_left       = viewOut(m_dir_car_coord_left);
  auto v_dir_def_coord_right      = viewOut(m_dir_def_coord_right);
  auto v_dir_car_coord_right      = viewOut(m_dir_car_coord_right);
  auto out_dir_vol1_left          = viewOut(m_dir_vol1_left);
  auto out_dir_vol1_right         = viewOut(m_dir_vol1_right);
  auto out_dir_vol2_left          = viewOut(m_dir_vol2_left);
  auto out_dir_vol2_right         = viewOut(m_dir_vol2_right);

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
      const auto node_id(cart_conn_cfn.node(inode));

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
        out_dir_vol1_left,
        out_dir_vol2_left
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
      out_dir_vol1_right[cell_id] = out_dir_vol1_left[next_cell_id];
      out_dir_vol2_right[cell_id] = out_dir_vol2_left[next_cell_id];

    } else {
      // Je suis sur la derniere rangee, je calcule
      // Acces "next" et grandeurs a droite
      lbd_dirvol(cell_id, cell_i, 
          v_dir_trans_area_right,
          v_face_velocity_right,
          v_dir_def_coord_right,
          v_dir_car_coord_right,
          out_dir_vol1_right,
          out_dir_vol2_right
          );
    }

  }  // Fin ENUMERATE_CELL

  P4GPU_STOP_TIMER(Loop_Cell1);
}

/*---------------------------------------------------------------------------*/
/* Parcours toutes les directions et appele _computeVolDir                   */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
computeVol() {
  Real dt=globalDeltaT();

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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

