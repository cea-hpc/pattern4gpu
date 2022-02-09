#ifndef PATTERN_4_GPU_4_KOKKOS_H
#define PATTERN_4_GPU_4_KOKKOS_H

#include <algorithm>

#include "Kokkos_Core.hpp"

#include "arcane/VariableTypedef.h"
#include "arcane/MeshVariableScalarRef.h"
#include "arcane/MeshVariableArrayRefT.H"
#include "arcane/Item.h"
#include "arcane/ItemGroup.h"
#include "arcane/utils/UtilsTypes.h"
#include "arcane/utils/Real3.h"
#include "arcane/utils/ArrayView.h"
#include "arcane/MathUtils.h"

#define FOR_KOKKOS
#include "SentinelArray.h"

// Pour tester de plomber la coalescence des Kokkos::View
// #define TargetMem Kokkos::LayoutRight,Kokkos::CudaUVMSpace

/*---------------------------------------------------------------------------*/
/* Classe encapsulante pour les méthodes qui vont utilisées Kokkos           */
/*---------------------------------------------------------------------------*/

struct KokkosWrapper {
  using TargetMem = Kokkos::CudaSpace;
  // using TargetMem = Kokkos::CudaHostPinnedSpace
  // using TargetMem = Kokkos::CudaUVMSpace;
  // using TargetMem = Kokkos::HostSpace;

  using TargetExec = Kokkos::Cuda;
  // using TargetExec = Kokkos::HPX;  // not possible ATM
  // using TargetExec = Kokkos::OpenMP;
  // using TargetExec = Kokkos::Serial;
  // using TargetExec = Kokkos::Threads;
  
  void init(const Arcane::CellGroup& all_cells, const Arcane::NodeGroup& all_nodes, const Arcane::VariableCellByte& is_active_cell, const Arcane::Span<const Arcane::Int16>& node_index_in_cells);
  
  // void initNodeVector(const Arcane::ConstArrayView<Arcane::Real3>& node_coord);
  void initNodeVector(const Arcane::VariableNodeReal3& node_coord, const Arcane::NodeGroup& all_nodes);
  
  // void initNodeCoordBis(const Arcane::ConstArrayView<Arcane::Real3>& node_coord);
  void initNodeCoordBis(const Arcane::VariableNodeReal3& node_coord, const Arcane::NodeGroup& all_nodes);
    
  void initCqs();
  
  void initCellArr12(const Arcane::CellGroup& all_cells, const Arcane::VariableNodeReal3& node_coord);
  
  void computeCqsAndVector();
  void computeCqsAndVectorV2();
  
  void syncHostData(const Arcane::CellGroup& all_cells, const Arcane::NodeGroup& all_nodes,
                    Arcane::VariableNodeReal3 node_vector, Arcane::VariableNodeReal3 node_coord_bis,
                    Arcane::VariableCellArrayReal3 cell_cqs, Arcane::VariableCellReal cell_arr1,
                    Arcane::VariableCellReal cell_arr2);
  
  static void end() {Kokkos::finalize();}
  
  
  // Pour les active_cells:  // Arcane::VariableCellByte comme flag sur les allCells dans le module
  Kokkos::View<Arcane::Int16*, TargetMem> m_is_active_cell;
  // Pour un minimum de connectivité active cell Id -> Node Ids  // Maximum 8 noeuds par maille
  Kokkos::View<Arcane::Int32*[8], TargetMem> m_cell_node_id;
  // Pour de la connectivité Node -> Cells  // Maximum 8 mailles par noeuds
  Kokkos::View<SentinelArray<Arcane::Int32, 8>*, TargetMem> m_node_cell_id;
  
  // Pb de recalage numerique, clairement un pb d'equivalence entre les idx des noeuds dans les cellules... j'abdique 
  Kokkos::View<Arcane::Int16*, TargetMem> m_node_index_in_cells;
  
  
  // Pour les coordonnées des noeuds (bis)  //   Arcane::VariableNodeReal3 dans le module
  Kokkos::View<Arcane::Real3*, TargetMem> m_node_coord_bis;
  
  // pour les cqs  // Arcane::VariableCellArrayReal3 dans le module
  Kokkos::View<Arcane::Real3*[8], TargetMem> m_cell_cqs;
  
  // pour les node vector  //   Arcane::VariableNodeReal3 dans le module
  Kokkos::View<Arcane::Real3*, TargetMem> m_node_vector;
  
  // pour deux quantités aux cells  // Arcane::VariableCellReal dans le module
  Kokkos::View<Arcane::Real*, TargetMem> m_cell_arr1;
  Kokkos::View<Arcane::Real*, TargetMem> m_cell_arr2;
  
  // Nb cells
  size_t m_nb_cells;
  size_t m_nb_nodes;
};

#endif  // PATTERN_4_GPU_4_KOKKOS_H

