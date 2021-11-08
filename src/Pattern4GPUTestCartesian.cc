#include "Pattern4GPUModule.h"
#include "arcane/Assertion.h"

#include "cartesian/CartesianMeshProperties.h"
#include "cartesian/CartesianItemSorter.h"

#include "cartesian/CartesianFaceId.h"
#include "cartesian/FactCartDirectionMng.h"

#include "arcane/VariableView.h"
#include "arcane/cea/CellDirectionMng.h"
#include "arcane/cea/NodeDirectionMng.h"
#include "arcane/cea/FaceDirectionMng.h"

/*---------------------------------------------------------------------------*/
/*!
 * \brief Test des passages de Cell => Cell (avec stencil) par direction
 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_testCell2Cell() {
  VariableCellInteger a_pcid(VariableBuildInfo(mesh(), "TemporaryCellPcid"));
  VariableCellInteger a_ccid(VariableBuildInfo(mesh(), "TemporaryCellCcid"));
  VariableCellInteger a_ncid(VariableBuildInfo(mesh(), "TemporaryCellNcid"));

  VariableCellInteger a_cid_m3(VariableBuildInfo(mesh(), "TemporaryCellCidM3"));
  VariableCellInteger a_cid_m2(VariableBuildInfo(mesh(), "TemporaryCellCidM2"));
  VariableCellInteger a_cid_m1(VariableBuildInfo(mesh(), "TemporaryCellCidM1"));
  VariableCellInteger a_cid_0 (VariableBuildInfo(mesh(), "TemporaryCellCid0"));
  VariableCellInteger a_cid_p1(VariableBuildInfo(mesh(), "TemporaryCellCidP1"));
  VariableCellInteger a_cid_p2(VariableBuildInfo(mesh(), "TemporaryCellCidP2"));
  VariableCellInteger a_cid_p3(VariableBuildInfo(mesh(), "TemporaryCellCidP3"));

  Cartesian::FactCartDirectionMng cartesian_mesh(mesh());

  auto queue = m_acc_env->newQueue();

  for(Integer dir(0) ; dir < mesh()->dimension() ; ++dir) {
    
    // C2C
    auto command = makeCommand(queue);


    auto out_pcid = ax::viewOut(command, a_pcid);
    auto out_ccid = ax::viewOut(command, a_ccid);
    auto out_ncid = ax::viewOut(command, a_ncid);

    auto out_cid_m3 = ax::viewOut(command, a_cid_m3);
    auto out_cid_m2 = ax::viewOut(command, a_cid_m2);
    auto out_cid_m1 = ax::viewOut(command, a_cid_m1);
    auto out_cid_0  = ax::viewOut(command, a_cid_0);
    auto out_cid_p1 = ax::viewOut(command, a_cid_p1);
    auto out_cid_p2 = ax::viewOut(command, a_cid_p2);
    auto out_cid_p3 = ax::viewOut(command, a_cid_p3);

    auto cell_dm = cartesian_mesh.cellDirection(dir);
    auto c2cid_stm = cell_dm.cell2CellIdStencil();
    auto cell_group = cell_dm.allCells();

    command << RUNCOMMAND_LOOP(iter, cell_group.loopRanges()) {
      auto [cid, idx] = c2cid_stm.idIdx(iter);

      auto c2cid = c2cid_stm.stencilCell<3>(cid,idx);
      out_cid_m3[cid] = c2cid(-3).localId();
      out_cid_m2[cid] = c2cid(-2).localId();
      out_cid_m1[cid] = c2cid(-1).localId();
      out_cid_0 [cid] = c2cid( 0).localId();
      out_cid_p1[cid] = c2cid(+1).localId();
      out_cid_p2[cid] = c2cid(+2).localId();
      out_cid_p3[cid] = c2cid(+3).localId();

      out_pcid[cid] = c2cid.previousId().localId();
      out_ccid[cid] = c2cid.centralId().localId();
      out_ncid[cid] = c2cid.nextId().localId();
    };

    auto arc_cell_dm = m_arc_cartesian_mesh->cellDirection(dir);

    auto in_pcid = viewIn(a_pcid);
    auto in_ccid = viewIn(a_ccid);
    auto in_ncid = viewIn(a_ncid);

    auto in_cid_m3 = viewIn(a_cid_m3);
    auto in_cid_m2 = viewIn(a_cid_m2);
    auto in_cid_m1 = viewIn(a_cid_m1);
    auto in_cid_0  = viewIn(a_cid_0);
    auto in_cid_p1 = viewIn(a_cid_p1);
    auto in_cid_p2 = viewIn(a_cid_p2);
    auto in_cid_p3 = viewIn(a_cid_p3);

    Assertion ass;

    ENUMERATE_CELL(cell_i, arc_cell_dm.allCells()) {
      CellLocalId cid(cell_i.localId());
      auto dir_cell = arc_cell_dm.cell(cid);
      CellLocalId pcid(dir_cell.previousId());
      CellLocalId ncid(dir_cell.nextId());

      ass.ASSERT_EQUAL(pcid.localId(), in_pcid[cid]);
      ass.ASSERT_EQUAL(cid.localId() , in_ccid[cid]);
      ass.ASSERT_EQUAL(ncid.localId(), in_ncid[cid]);

      ass.ASSERT_EQUAL(pcid.localId(), in_cid_m1[cid]);
      ass.ASSERT_EQUAL(cid.localId() , in_cid_0[cid]);
      ass.ASSERT_EQUAL(ncid.localId(), in_cid_p1[cid]);

      // previous ...
      if (!ItemId::null(pcid)) {
        auto dir_cell_m1 = arc_cell_dm.cell(pcid);
        CellLocalId pcid_m2(dir_cell_m1.previousId());
        ass.ASSERT_EQUAL(pcid_m2.localId(), in_cid_m2[cid]);
        if (!ItemId::null(pcid_m2)) {
          auto dir_cell_m2 = arc_cell_dm.cell(pcid_m2);
          CellLocalId pcid_m3(dir_cell_m2.previousId());
          ass.ASSERT_EQUAL(pcid_m3.localId(), in_cid_m3[cid]);
        }
      }

      // next ...
      if (!ItemId::null(ncid)) {
        auto dir_cell_m1 = arc_cell_dm.cell(ncid);
        CellLocalId ncid_p2(dir_cell_m1.nextId());
        ass.ASSERT_EQUAL(ncid_p2.localId(), in_cid_p2[cid]);
        if (!ItemId::null(ncid_p2)) {
          auto dir_cell_m2 = arc_cell_dm.cell(ncid_p2);
          CellLocalId ncid_p3(dir_cell_m2.nextId());
          ass.ASSERT_EQUAL(ncid_p3.localId(), in_cid_p3[cid]);
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Test des passages de Node => Node (avec stencil) par direction
 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_testNode2Node() {
  VariableNodeInteger a_pnid(VariableBuildInfo(mesh(), "TemporaryNodePnid"));
  VariableNodeInteger a_cnid(VariableBuildInfo(mesh(), "TemporaryNodeCnid"));
  VariableNodeInteger a_nnid(VariableBuildInfo(mesh(), "TemporaryNodeNnid"));

  VariableNodeInteger a_nid_m3(VariableBuildInfo(mesh(), "TemporaryNodeNidM3"));
  VariableNodeInteger a_nid_m2(VariableBuildInfo(mesh(), "TemporaryNodeNidM2"));
  VariableNodeInteger a_nid_m1(VariableBuildInfo(mesh(), "TemporaryNodeNidM1"));
  VariableNodeInteger a_nid_0 (VariableBuildInfo(mesh(), "TemporaryNodeNid0"));
  VariableNodeInteger a_nid_p1(VariableBuildInfo(mesh(), "TemporaryNodeNidP1"));
  VariableNodeInteger a_nid_p2(VariableBuildInfo(mesh(), "TemporaryNodeNidP2"));
  VariableNodeInteger a_nid_p3(VariableBuildInfo(mesh(), "TemporaryNodeNidP3"));

  Cartesian::FactCartDirectionMng cartesian_mesh(mesh());

  auto queue = m_acc_env->newQueue();

  for(Integer dir(0) ; dir < mesh()->dimension() ; ++dir) {
    
    // N2N
    auto command = makeCommand(queue);


    auto out_pnid = ax::viewOut(command, a_pnid);
    auto out_cnid = ax::viewOut(command, a_cnid);
    auto out_nnid = ax::viewOut(command, a_nnid);

    auto out_nid_m3 = ax::viewOut(command, a_nid_m3);
    auto out_nid_m2 = ax::viewOut(command, a_nid_m2);
    auto out_nid_m1 = ax::viewOut(command, a_nid_m1);
    auto out_nid_0  = ax::viewOut(command, a_nid_0);
    auto out_nid_p1 = ax::viewOut(command, a_nid_p1);
    auto out_nid_p2 = ax::viewOut(command, a_nid_p2);
    auto out_nid_p3 = ax::viewOut(command, a_nid_p3);

    auto node_dm = cartesian_mesh.nodeDirection(dir);
    auto n2nid_stm = node_dm.node2NodeIdStencil();
    auto node_group = node_dm.allNodes();

    command << RUNCOMMAND_LOOP(iter, node_group.loopRanges()) {
      auto [nid, idx] = n2nid_stm.idIdx(iter);

      auto n2nid = n2nid_stm.stencilNode<3>(nid,idx);
      out_nid_m3[nid] = n2nid(-3).localId();
      out_nid_m2[nid] = n2nid(-2).localId();
      out_nid_m1[nid] = n2nid(-1).localId();
      out_nid_0 [nid] = n2nid( 0).localId();
      out_nid_p1[nid] = n2nid(+1).localId();
      out_nid_p2[nid] = n2nid(+2).localId();
      out_nid_p3[nid] = n2nid(+3).localId();

      out_pnid[nid] = n2nid.previousId().localId();
      out_cnid[nid] = n2nid.centralId().localId();
      out_nnid[nid] = n2nid.nextId().localId();
    };

    auto arc_node_dm = m_arc_cartesian_mesh->nodeDirection(dir);

    auto in_pnid = viewIn(a_pnid);
    auto in_cnid = viewIn(a_cnid);
    auto in_nnid = viewIn(a_nnid);

    auto in_nid_m3 = viewIn(a_nid_m3);
    auto in_nid_m2 = viewIn(a_nid_m2);
    auto in_nid_m1 = viewIn(a_nid_m1);
    auto in_nid_0  = viewIn(a_nid_0);
    auto in_nid_p1 = viewIn(a_nid_p1);
    auto in_nid_p2 = viewIn(a_nid_p2);
    auto in_nid_p3 = viewIn(a_nid_p3);

    Assertion ass;

    ENUMERATE_NODE(node_i, arc_node_dm.allNodes()) {
      NodeLocalId nid(node_i.localId());
      auto dir_node = arc_node_dm.node(nid);
      NodeLocalId pnid(dir_node.previousId());
      NodeLocalId nnid(dir_node.nextId());

      ass.ASSERT_EQUAL(pnid.localId(), in_pnid[nid]);
      ass.ASSERT_EQUAL(nid.localId() , in_cnid[nid]);
      ass.ASSERT_EQUAL(nnid.localId(), in_nnid[nid]);

      ass.ASSERT_EQUAL(pnid.localId(), in_nid_m1[nid]);
      ass.ASSERT_EQUAL(nid.localId() , in_nid_0[nid]);
      ass.ASSERT_EQUAL(nnid.localId(), in_nid_p1[nid]);

      // previous ...
      if (!ItemId::null(pnid)) {
        auto dir_node_m1 = arc_node_dm.node(pnid);
        NodeLocalId pnid_m2(dir_node_m1.previousId());
        ass.ASSERT_EQUAL(pnid_m2.localId(), in_nid_m2[nid]);
        if (!ItemId::null(pnid_m2)) {
          auto dir_node_m2 = arc_node_dm.node(pnid_m2);
          NodeLocalId pnid_m3(dir_node_m2.previousId());
          ass.ASSERT_EQUAL(pnid_m3.localId(), in_nid_m3[nid]);
        }
      }

      // next ...
      if (!ItemId::null(nnid)) {
        auto dir_node_m1 = arc_node_dm.node(nnid);
        NodeLocalId nnid_p2(dir_node_m1.nextId());
        ass.ASSERT_EQUAL(nnid_p2.localId(), in_nid_p2[nid]);
        if (!ItemId::null(nnid_p2)) {
          auto dir_node_m2 = arc_node_dm.node(nnid_p2);
          NodeLocalId nnid_p3(dir_node_m2.nextId());
          ass.ASSERT_EQUAL(nnid_p3.localId(), in_nid_p3[nid]);
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Test des passages de Face => Cell (avec stencil) par direction
 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_testFace2Cell() {
  VariableFaceInteger a_cid_m3(VariableBuildInfo(mesh(), "TemporaryFaceCidM3"));
  VariableFaceInteger a_cid_m2(VariableBuildInfo(mesh(), "TemporaryFaceCidM2"));
  VariableFaceInteger a_cid_m1(VariableBuildInfo(mesh(), "TemporaryFaceCidM1"));
  VariableFaceInteger a_cid_p1(VariableBuildInfo(mesh(), "TemporaryFaceCidP1"));
  VariableFaceInteger a_cid_p2(VariableBuildInfo(mesh(), "TemporaryFaceCidP2"));
  VariableFaceInteger a_cid_p3(VariableBuildInfo(mesh(), "TemporaryFaceCidP3"));

  Cartesian::FactCartDirectionMng cartesian_mesh(mesh());

  auto queue = m_acc_env->newQueue();

  for(Integer dir(0) ; dir < mesh()->dimension() ; ++dir) {
    
    // F2C
    auto command = makeCommand(queue);


    auto out_cid_m3 = ax::viewOut(command, a_cid_m3);
    auto out_cid_m2 = ax::viewOut(command, a_cid_m2);
    auto out_cid_m1 = ax::viewOut(command, a_cid_m1);
    auto out_cid_p1 = ax::viewOut(command, a_cid_p1);
    auto out_cid_p2 = ax::viewOut(command, a_cid_p2);
    auto out_cid_p3 = ax::viewOut(command, a_cid_p3);

    auto face_dm = cartesian_mesh.faceDirection(dir);
    auto f2cid_stm = face_dm.face2CellIdStencil();
    auto face_group = face_dm.allFaces();

    command << RUNCOMMAND_LOOP(iter, face_group.loopRanges()) {
      auto [fid, idx] = f2cid_stm.idIdx(iter);

      auto f2cid = f2cid_stm.stencilFace2Cell<3>(fid,idx);
      out_cid_m3[fid] = f2cid.previousId(-3).localId();
      out_cid_m2[fid] = f2cid.previousId(-2).localId();
      out_cid_m1[fid] = f2cid.previousId(-1).localId();
      out_cid_p1[fid] = f2cid.nextId(+1).localId();
      out_cid_p2[fid] = f2cid.nextId(+2).localId();
      out_cid_p3[fid] = f2cid.nextId(+3).localId();
    };

    auto arc_face_dm = m_arc_cartesian_mesh->faceDirection(dir);
    auto arc_cell_dm = m_arc_cartesian_mesh->cellDirection(dir);

    auto in_cid_m3 = viewIn(a_cid_m3);
    auto in_cid_m2 = viewIn(a_cid_m2);
    auto in_cid_m1 = viewIn(a_cid_m1);
    auto in_cid_p1 = viewIn(a_cid_p1);
    auto in_cid_p2 = viewIn(a_cid_p2);
    auto in_cid_p3 = viewIn(a_cid_p3);

    Assertion ass;

    ENUMERATE_FACE(face_i, arc_face_dm.allFaces()) {
      FaceLocalId fid(face_i.localId());
      auto dir_face = arc_face_dm.face(fid);
      CellLocalId pcid(dir_face.previousCellId());
      CellLocalId ncid(dir_face.nextCellId());

      ass.ASSERT_EQUAL(pcid.localId(), in_cid_m1[fid]);
      ass.ASSERT_EQUAL(ncid.localId(), in_cid_p1[fid]);

      // previous ...
      if (!ItemId::null(pcid)) {
        auto dir_cell_m1 = arc_cell_dm.cell(pcid);
        CellLocalId pcid_m2(dir_cell_m1.previousId());
        ass.ASSERT_EQUAL(pcid_m2.localId(), in_cid_m2[fid]);
        if (!ItemId::null(pcid_m2)) {
          auto dir_cell_m2 = arc_cell_dm.cell(pcid_m2);
          CellLocalId pcid_m3(dir_cell_m2.previousId());
          ass.ASSERT_EQUAL(pcid_m3.localId(), in_cid_m3[fid]);
        }
      }

      // next ...
      if (!ItemId::null(ncid)) {
        auto dir_cell_m1 = arc_cell_dm.cell(ncid);
        CellLocalId ncid_p2(dir_cell_m1.nextId());
        ass.ASSERT_EQUAL(ncid_p2.localId(), in_cid_p2[fid]);
        if (!ItemId::null(ncid_p2)) {
          auto dir_cell_m2 = arc_cell_dm.cell(ncid_p2);
          CellLocalId ncid_p3(dir_cell_m2.nextId());
          ass.ASSERT_EQUAL(ncid_p3.localId(), in_cid_p3[fid]);
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Test des passages de Cell => Face (avec stencil) par direction
 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_testCell2Face() {
  VariableCellInteger a_fid_m3(VariableBuildInfo(mesh(), "TemporaryCellFidM3"));
  VariableCellInteger a_fid_m2(VariableBuildInfo(mesh(), "TemporaryCellFidM2"));
  VariableCellInteger a_fid_m1(VariableBuildInfo(mesh(), "TemporaryCellFidM1"));
  VariableCellInteger a_fid_p1(VariableBuildInfo(mesh(), "TemporaryCellFidP1"));
  VariableCellInteger a_fid_p2(VariableBuildInfo(mesh(), "TemporaryCellFidP2"));
  VariableCellInteger a_fid_p3(VariableBuildInfo(mesh(), "TemporaryCellFidP3"));

  Cartesian::FactCartDirectionMng cartesian_mesh(mesh());

  auto queue = m_acc_env->newQueue();

  for(Integer dir(0) ; dir < mesh()->dimension() ; ++dir) {
    
    // C2F
    auto command = makeCommand(queue);


    auto out_fid_m3 = ax::viewOut(command, a_fid_m3);
    auto out_fid_m2 = ax::viewOut(command, a_fid_m2);
    auto out_fid_m1 = ax::viewOut(command, a_fid_m1);
    auto out_fid_p1 = ax::viewOut(command, a_fid_p1);
    auto out_fid_p2 = ax::viewOut(command, a_fid_p2);
    auto out_fid_p3 = ax::viewOut(command, a_fid_p3);

    auto cell_dm = cartesian_mesh.cellDirection(dir);
    auto c2fid_stm = cell_dm.cell2FaceIdStencil();
    auto cell_group = cell_dm.allCells();

    command << RUNCOMMAND_LOOP(iter, cell_group.loopRanges()) {
      auto [cid, idx] = c2fid_stm.idIdx(iter);

      auto c2fid = c2fid_stm.stencilCell2Face<3>(cid,idx);
      out_fid_m3[cid] = c2fid.previousId(-3).localId();
      out_fid_m2[cid] = c2fid.previousId(-2).localId();
      out_fid_m1[cid] = c2fid.previousId(-1).localId();
      out_fid_p1[cid] = c2fid.nextId(+1).localId();
      out_fid_p2[cid] = c2fid.nextId(+2).localId();
      out_fid_p3[cid] = c2fid.nextId(+3).localId();
    };

    auto arc_face_dm = m_arc_cartesian_mesh->faceDirection(dir);
    auto arc_cell_dm = m_arc_cartesian_mesh->cellDirection(dir);

    auto in_fid_m3 = viewIn(a_fid_m3);
    auto in_fid_m2 = viewIn(a_fid_m2);
    auto in_fid_m1 = viewIn(a_fid_m1);
    auto in_fid_p1 = viewIn(a_fid_p1);
    auto in_fid_p2 = viewIn(a_fid_p2);
    auto in_fid_p3 = viewIn(a_fid_p3);

    Assertion ass;

    ENUMERATE_CELL(cell_i, arc_cell_dm.allCells()) {
      CellLocalId cid(cell_i.localId());
      auto dir_cell_face = arc_cell_dm.cellFace(cid);
      FaceLocalId pfid(dir_cell_face.previousId());
      FaceLocalId nfid(dir_cell_face.nextId());

      ass.ASSERT_EQUAL(pfid.localId(), in_fid_m1[cid]);
      ass.ASSERT_EQUAL(nfid.localId(), in_fid_p1[cid]);

      // previous ...
      if (!ItemId::null(pfid)) {
        auto dir_face_m1 = arc_face_dm.face(pfid);
        CellLocalId pcid_m1(dir_face_m1.previousCellId());
        if (!ItemId::null(pcid_m1)) {
          auto dir_cell_face_m1 = arc_cell_dm.cellFace(pcid_m1);
          FaceLocalId pfid_m2(dir_cell_face_m1.previousId());
          ass.ASSERT_EQUAL(pfid_m2.localId(), in_fid_m2[cid]);

          auto dir_face_m2 = arc_face_dm.face(pfid_m2);
          CellLocalId pcid_m2(dir_face_m2.previousCellId());
          if (!ItemId::null(pcid_m2)) {
            auto dir_cell_face_m2 = arc_cell_dm.cellFace(pcid_m2);
            FaceLocalId pfid_m3(dir_cell_face_m2.previousId());
            ass.ASSERT_EQUAL(pfid_m3.localId(), in_fid_m3[cid]);
          }
        }
      }

      // next ...
      if (!ItemId::null(nfid)) {
        auto dir_face_p1 = arc_face_dm.face(nfid);
        CellLocalId ncid_p1(dir_face_p1.nextCellId());
        if (!ItemId::null(ncid_p1)) {
          auto dir_cell_face_p1 = arc_cell_dm.cellFace(ncid_p1);
          FaceLocalId nfid_p2(dir_cell_face_p1.nextId());
          ass.ASSERT_EQUAL(nfid_p2.localId(), in_fid_p2[cid]);

          auto dir_face_p2 = arc_face_dm.face(nfid_p2);
          CellLocalId ncid_p2(dir_face_p2.nextCellId());
          if (!ItemId::null(ncid_p2)) {
            auto dir_cell_face_p2 = arc_cell_dm.cellFace(ncid_p2);
            FaceLocalId nfid_p3(dir_cell_face_p2.nextId());
            ass.ASSERT_EQUAL(nfid_p3.localId(), in_fid_p3[cid]);
          }
        }
      }
    }
  }
}

void Pattern4GPUModule::
_stencilCartesian() {
  PROF_ACC_BEGIN(__FUNCTION__);

#define DO_ASSERT 

  Cartesian::FactCartDirectionMng cartesian_mesh(mesh());

  auto queue = m_acc_env->newQueue();

  for(Integer dir(0) ; dir < mesh()->dimension() ; ++dir) {
    
    // C2C
    auto command = makeCommand(queue);

    auto in_cell_arr2 = ax::viewIn(command, m_cell_arr2);
    auto inout_cell_arr1 = ax::viewInOut(command, m_cell_arr1);

    auto cart_cell_dm = cartesian_mesh.cellDirection(dir);
    auto c2cid_stm = cart_cell_dm.cell2CellIdStencil();
    //auto cell_group = cart_cell_dm.innerCells();
    auto cell_group = cart_cell_dm.allCells();

    command << RUNCOMMAND_LOOP(iter, cell_group.loopRanges()) {
      auto [cid, idx] = c2cid_stm.idIdx(iter);

      // Acces mailles gauche/droite
      auto c2cid = c2cid_stm.cell(cid, idx);
      CellLocalId pcid(c2cid.previous());
      CellLocalId ncid(c2cid.next());

      Real sum=0.;
      if (!ItemId::null(pcid))
        sum+=in_cell_arr2[pcid];
      if (!ItemId::null(ncid))
        sum+=in_cell_arr2[ncid];

      auto c2cid_st3 = c2cid_stm.stencilCell<3>(cid, idx);
      // Acces mailles stencil - façon 1
      Real sum_st1=0.;
      for(Integer ilayer=-3/*-c2cid_st3.nLayer()*/ ; ilayer<=3/*c2cid_st3.nLayer()*/ ; ilayer++) {
        // acces à la maille de la couche ilayer
        // Rem1 : ilayer=0 => cid ;  Rem2 : la maille peut ne pas être valide si ext. au domaine
        CellLocalId cid_st(c2cid_st3(ilayer)); // acces à la maille de la couche ilayer
        if (!ItemId::null(cid_st))
          sum_st1 += in_cell_arr2[cid_st];
      }

      // Acces mailles stencil - façon 2
      Real sum_st2=0.;
      for(Integer ilayer=c2cid_st3.validMin() ; ilayer<=c2cid_st3.validMax() ; ilayer++) {
        // acces à la maille de la couche ilayer, maille valide car dans [validMin(),validMax()]
        sum_st2 += in_cell_arr2[ c2cid_st3(ilayer) ];
      }

#if 0 // TODO
      // Acces mailles stencil - façon 3
      Real sum_st3=0.;
      for(Integer ilayer : c2cid_st3.validLayers()) {
        // acces à la maille de la couche ilayer, maille valide car dans [validMin(),validMax()]
        sum_st2 += in_cell_arr2[ c2cid_st3(ilayer) ];
      }
#endif

#if !defined(ARCCORE_DEVICE_CODE) && defined(DO_ASSERT)
      Assertion ass;
      ass.ASSERT_NEARLY_EQUAL_EPSILON(sum_st1,sum_st2,Real(1.e-12));
#endif

      inout_cell_arr1[cid]+=sum+sum_st2;
    };

    // N2N
    auto command2 = makeCommand(queue);

    auto in_node_arr1 = ax::viewIn(command2, m_node_arr1);
    auto inout_node_arr2 = ax::viewInOut(command2, m_node_arr2);

    auto cart_node_dm = cartesian_mesh.nodeDirection(dir);
    auto n2nid_stm = cart_node_dm.node2NodeIdStencil();
    //auto node_group = cart_node_dm.innerNodes();
    auto node_group = cart_node_dm.allNodes();

    command2 << RUNCOMMAND_LOOP(iter, node_group.loopRanges()) {
      auto [nid, idx] = n2nid_stm.idIdx(iter);

      auto n2nid_st3 = n2nid_stm.stencilNode<3>(nid, idx);
      // Acces noeuds stencil - façon 1
      Real sum_st1=0.;
      for(Integer ilayer=-3/*-n2nid_st3.nLayer()*/ ; ilayer<=3/*n2nid_st3.nLayer()*/ ; ilayer++) {
        // acces au noeud de la couche ilayer
        // Rem1 : ilayer=0 => nid ;  Rem2 : le noeud peut ne pas être valide si ext. au domaine
        NodeLocalId nid_st(n2nid_st3(ilayer)); // acces au noeud de la couche ilayer
        if (!ItemId::null(nid_st))
          sum_st1 += in_node_arr1[nid_st];
      }

      // Acces noeuds stencil - façon 2
      Real sum_st2=0.;
      for(Integer ilayer=n2nid_st3.validMin() ; ilayer<=n2nid_st3.validMax() ; ilayer++) {
        // acces au noeud de la couche ilayer, noeud valide car dans [validMin(),validMax()]
        sum_st2 += in_node_arr1[ n2nid_st3(ilayer) ];
      }

#if !defined(ARCCORE_DEVICE_CODE) && defined(DO_ASSERT)
      Assertion ass;
      ass.ASSERT_NEARLY_EQUAL_EPSILON(sum_st1,sum_st2,Real(1.e-12));
#endif

      inout_node_arr2[nid]+=sum_st1+sum_st2;
    };

    // C2F
    auto command3 = makeCommand(queue);

    auto in_face_arr1 = ax::viewIn(command3, m_face_arr1);
    auto inout_cell_arr2 = ax::viewInOut(command3, m_cell_arr2);

    //auto cart_cell_dm = cartesian_mesh.cellDirection(dir);
    auto c2fid_stm = cart_cell_dm.cell2FaceIdStencil();
    //auto cell_group = cart_cell_dm.innerCells();
    //auto cell_group = cart_cell_dm.allCells();

    command3 << RUNCOMMAND_LOOP(iter, cell_group.loopRanges()) {
      auto [cid, idx] = c2fid_stm.idIdx(iter);

      // Acces faces gauche/droite qui existent forcement
      auto c2fid = c2fid_stm.cellFace(cid, idx);
      FaceLocalId pfid(c2fid.previousId());
      FaceLocalId nfid(c2fid.nextId());

      Real sum=in_face_arr1[pfid]+in_face_arr1[nfid];

      auto c2fid_st3 = c2fid_stm.stencilCell2Face<3>(cid, idx);
      // Acces faces "previous" dans stencil - façon 1
      Real sum_st1=0.;
      for(Integer ilayer=-3/*-c2fid_st3.nLayer()*/ ; ilayer<=-1 ; ilayer++) {
        // acces à la face de la couche ilayer
        // Rem1 : ilayer=0 => fid ;  Rem2 : la face peut ne pas être valide si ext. au domaine
        FaceLocalId fid_st(c2fid_st3.previousId(ilayer)); // acces à la face de la couche ilayer
        if (!ItemId::null(fid_st))
          sum_st1 += in_face_arr1[fid_st];
      }
      // Acces faces "next" dans stencil - façon 1
      for(Integer ilayer=+1 ; ilayer<=+3/*c2fid_st3.nLayer()*/ ; ilayer++) {
        // acces à la face de la couche ilayer
        // Rem1 : ilayer=0 => fid ;  Rem2 : la face peut ne pas être valide si ext. au domaine
        FaceLocalId fid_st(c2fid_st3.nextId(ilayer)); // acces à la face de la couche ilayer
        if (!ItemId::null(fid_st))
          sum_st1 += in_face_arr1[fid_st];
      }

      // Acces faces stencil - façon 2
      Real sum_st2=0.;
      for(Integer ilayer=c2fid_st3.validMin() ; ilayer<=-1 ; ilayer++) {
        // acces à la face de la couche ilayer, face valide car dans [validMin(),-1]
        sum_st2 += in_face_arr1[ c2fid_st3.previousId(ilayer) ];
      }
      for(Integer ilayer=+1 ; ilayer<=c2fid_st3.validMax() ; ilayer++) {
        // acces à la face de la couche ilayer, face valide car dans [+1,validMax()]
        sum_st2 += in_face_arr1[ c2fid_st3.nextId(ilayer) ];
      }

#if !defined(ARCCORE_DEVICE_CODE) && defined(DO_ASSERT)
      Assertion ass;
      ass.ASSERT_NEARLY_EQUAL_EPSILON(sum_st1,sum_st2,Real(1.e-12));
#endif

      inout_cell_arr2[cid]+=sum+sum_st2;
    };

    // F2C
    auto command4 = makeCommand(queue);

    auto in_cell_arr1 = ax::viewIn(command4, m_cell_arr1);
    auto inout_face_arr1 = ax::viewInOut(command4, m_face_arr1);

    auto cart_face_dm = cartesian_mesh.faceDirection(dir);
    auto f2cid_stm = cart_face_dm.face2CellIdStencil();
    //auto face_group = cart_face_dm.innerFaces();
    auto face_group = cart_face_dm.allFaces();

    command4 << RUNCOMMAND_LOOP(iter, face_group.loopRanges()) {
      auto [fid, idx] = f2cid_stm.idIdx(iter);

      // Acces mailles gauche/droite 
      auto f2cid = f2cid_stm.face(fid, idx);
      CellLocalId pcid(f2cid.previousCell());
      CellLocalId ncid(f2cid.nextCell());

      Real sum=0.;
      if (!ItemId::null(pcid))
        sum+=in_cell_arr1[pcid];
      if (!ItemId::null(ncid))
        sum+=in_cell_arr1[ncid];


      auto f2cid_st3 = f2cid_stm.stencilFace2Cell<3>(fid, idx);
      // Acces mailles "previous" dans stencil - façon 1
      Real sum_st1=0.;
      for(Integer ilayer=-3/*-f2cid_st3.nLayer()*/ ; ilayer<=-1 ; ilayer++) {
        // acces à la maille de la couche ilayer
        // Rem1 : ilayer=0 => fid ;  Rem2 : la maille peut ne pas être valide si ext. au domaine
        CellLocalId cid_st(f2cid_st3.previousId(ilayer)); // acces à la maille de la couche ilayer
        if (!ItemId::null(cid_st))
          sum_st1 += in_cell_arr1[cid_st];
      }
      // Acces mailles "next" dans stencil - façon 1
      for(Integer ilayer=+1 ; ilayer<=+3/*f2cid_st3.nLayer()*/ ; ilayer++) {
        // acces à la maille de la couche ilayer
        // Rem1 : ilayer=0 => fid ;  Rem2 : la maille peut ne pas être valide si ext. au domaine
        CellLocalId cid_st(f2cid_st3.nextId(ilayer)); // acces à la maille de la couche ilayer
        if (!ItemId::null(cid_st))
          sum_st1 += in_cell_arr1[cid_st];
      }

      // Acces faces stencil - façon 2
      Real sum_st2=0.;
      for(Integer ilayer=f2cid_st3.validMin() ; ilayer<=-1 ; ilayer++) {
        // acces à la maille de la couche ilayer, maille valide car dans [validMin(),-1]
        sum_st2 += in_cell_arr1[ f2cid_st3.previousId(ilayer) ];
      }
      for(Integer ilayer=+1 ; ilayer<=f2cid_st3.validMax() ; ilayer++) {
        // acces à la maille de la couche ilayer, maille valide car dans [+1,validMax()]
        sum_st2 += in_cell_arr1[ f2cid_st3.nextId(ilayer) ];
      }

#if !defined(ARCCORE_DEVICE_CODE) && defined(DO_ASSERT)
      Assertion ass;
      ass.ASSERT_NEARLY_EQUAL_EPSILON(sum_st1,sum_st2,Real(1.e-12));
#endif

      inout_face_arr1[fid]+=sum+sum_st2;
    };

  }
  PROF_ACC_END;
}

void Pattern4GPUModule::
testCartesian() {
  PROF_ACC_BEGIN(__FUNCTION__);

  _testCell2Cell();
  _testNode2Node();
  _testFace2Cell();
  _testCell2Face();
  _stencilCartesian();

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

