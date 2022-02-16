#include "Pattern4GPUModule.h"

#define P4GPU_PROFILING // Pour activer le profiling
#include "P4GPUTimer.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
updateVectorFromTensor() {

  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "Dans updateVectorFromTensor";

  if (options()->getUpdateVectorFromTensorVersion() == UVTV_ori)
  {
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
  else if (options()->getUpdateVectorFromTensorVersion() == UVTV_mt)
  {
    Integer nb_blocks = m_mesh_material_mng->blocks().size();

    for (Integer i = 0; i < nb_blocks; i++) {
      // on récupère le groupe de mailles associé au bloc
      IMeshBlock* b = (m_mesh_material_mng->blocks())[i];
      CellGroup cell_group = b->cells();

      auto node_index_in_cells = m_acc_env->nodeIndexInCells();
      const Integer max_node_cell = m_acc_env->maxNodeCell();

      ParallelLoopOptions options;
      options.setPartitioner(ParallelLoopOptions::Partitioner::Auto);

      NodeGroup node_group = allNodes(); // TODO : passer aux noeuds du blocks
      arcaneParallelForeach(node_group, options, [&](NodeVectorView nodes) {
      ENUMERATE_NODE (node_i, nodes) {
        Int32 first_pos = node_i.localId() * max_node_cell;
        ENUMERATE_CELL(cell_i, node_i->cells()) {
          if (true) { // TODO : vrai ssi cell_i est dans cell_group
            Int16 node_index = node_index_in_cells[first_pos + cell_i.index()];
            const Real3x3& cell_tensor = m_tensor[cell_i];
            m_node_vector[node_i] -= math::prodTensVec(cell_tensor,
              m_cell_cqs[cell_i][node_index]);
          }
        }
      }
      });
    }  // end iblock loop
  }
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

