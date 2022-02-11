#ifdef WITH_KOKKOS

#include "Pattern4GPU4Kokkos.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void KokkosWrapper::init(const Arcane::CellGroup& all_cells,
                         const Arcane::NodeGroup& all_nodes,
                         const Arcane::VariableCellByte& is_active_cell,
                         const Arcane::Span<const Arcane::Int16>& node_index_in_cells)
{
  Kokkos::initialize();
  
  m_nb_cells = all_cells.size();
  m_nb_nodes = all_nodes.size();
  
  Kokkos::resize(m_cell_node_id, m_nb_cells);
  Kokkos::resize(m_node_cell_id, m_nb_nodes);
  Kokkos::resize(m_is_active_cell, m_nb_cells);
  Kokkos::resize(m_node_index_in_cells, m_nb_nodes*8);

  Kokkos::View<Arcane::Int32*[8], TargetMem>::HostMirror cell_node_id_host =
    Kokkos::create_mirror_view(m_cell_node_id);
  Kokkos::View<SentinelArray<Arcane::Int32, 8>*, TargetMem>::HostMirror node_cell_id_host =
    Kokkos::create_mirror_view(m_node_cell_id);
  Kokkos::View<Arcane::Int16*, TargetMem>::HostMirror is_active_cell_host =
    Kokkos::create_mirror_view(m_is_active_cell);
  Kokkos::View<Arcane::Int16*, TargetMem>::HostMirror node_index_in_cells_host =
    Kokkos::create_mirror_view(m_node_index_in_cells);
  
  ENUMERATE_CELL(icell, all_cells) {
    cell_node_id_host(icell->localId(), 0) = icell->node(0).localId();
    cell_node_id_host(icell->localId(), 1) = icell->node(1).localId();
    cell_node_id_host(icell->localId(), 2) = icell->node(2).localId();
    cell_node_id_host(icell->localId(), 3) = icell->node(3).localId();
    cell_node_id_host(icell->localId(), 4) = icell->node(4).localId();
    cell_node_id_host(icell->localId(), 5) = icell->node(5).localId();
    cell_node_id_host(icell->localId(), 6) = icell->node(6).localId();
    cell_node_id_host(icell->localId(), 7) = icell->node(7).localId();
    
    is_active_cell_host(icell->localId()) = is_active_cell[icell];
  }
  
  ENUMERATE_NODE(inode, all_nodes) {
    node_cell_id_host(inode->localId()) =
      SentinelArray<Arcane::Int32, 8>(inode->cells().localIds().begin(),
                                      inode->cells().localIds().end());
  }
  
  for (int i(0); i < node_index_in_cells.size(); ++i)
    node_index_in_cells_host(i) = node_index_in_cells[i];
  
  Kokkos::deep_copy(m_cell_node_id, cell_node_id_host);  // H2D
  Kokkos::deep_copy(m_node_cell_id, node_cell_id_host);  // H2D
  Kokkos::deep_copy(m_is_active_cell, is_active_cell_host);  // H2D
  Kokkos::deep_copy(m_node_index_in_cells, node_index_in_cells_host);  // H2D
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void KokkosWrapper::initNodeVector(const Arcane::VariableNodeReal3& node_coord,
                                   const Arcane::NodeGroup& all_nodes)
{
  Kokkos::resize(m_node_vector, m_nb_nodes);
  // Create host mirror of node_vector
  Kokkos::View<Arcane::Real3*, TargetMem>::HostMirror node_vector_host =
    Kokkos::create_mirror_view(m_node_vector);
  
  ENUMERATE_NODE(node_i, all_nodes) {
    const Arcane::Real3& c=node_coord[node_i];
    Arcane::Real cos_th=cos(c.x+c.y+c.z);  // garantit une valeur dans [-1,+1]
    Arcane::Real sin_th=Arcane::math::sqrt(1-cos_th*cos_th);  // garantit une valeur dans [0,+1]
    Arcane::Real phi=(c.x+1)*(c.y+1)*(c.z+1);
    node_vector_host(node_i->localId())=Arcane::Real3(sin_th*cos(phi), sin_th*sin(phi),cos_th);
  }

  Kokkos::deep_copy(m_node_vector, node_vector_host);  // H2D
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

//TODO: essayez de copier la var arcane qui reside deja en device directement, sans hostmirror
void KokkosWrapper::initNodeCoordBis(const Arcane::VariableNodeReal3& node_coord, const Arcane::NodeGroup& all_nodes)
{
  Kokkos::resize(m_node_coord_bis, m_nb_nodes);
  Kokkos::View<Arcane::Real3*, TargetMem>::HostMirror node_coord_bis_host =
    Kokkos::create_mirror_view(m_node_coord_bis);
  
  ENUMERATE_NODE(node_i, all_nodes) {
    node_coord_bis_host(node_i->localId()) = node_coord[node_i];
  }

  Kokkos::deep_copy(m_node_coord_bis, node_coord_bis_host);  // H2D
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void KokkosWrapper::initCqs()
{
  Kokkos::resize(m_cell_cqs, m_nb_cells);
  Kokkos::View<Arcane::Real3*[8], TargetMem>::HostMirror cell_cqs_host =
    Kokkos::create_mirror_view(m_cell_cqs);

  for (Arcane::Integer icell(0); icell < m_nb_cells; ++icell) {
    for(Arcane::Integer inode(0) ; inode < 8 ; ++inode) {
      cell_cqs_host(icell, inode) = Arcane::Real3::zero();
    }
  }

  Kokkos::deep_copy(m_cell_cqs, cell_cqs_host);  // H2D
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void KokkosWrapper::initCellArr12(const Arcane::CellGroup& all_cells,
                                  const Arcane::VariableNodeReal3& node_coord)
{
  Kokkos::resize(m_cell_arr1, m_nb_cells);
  Kokkos::resize(m_cell_arr2, m_nb_cells);
  
  Kokkos::View<Arcane::Real*, TargetMem>::HostMirror cell_arr1_host =
    Kokkos::create_mirror_view(m_cell_arr1);
  Kokkos::View<Arcane::Real*, TargetMem>::HostMirror cell_arr2_host =
    Kokkos::create_mirror_view(m_cell_arr2);

  ENUMERATE_CELL(icell, all_cells) {
    const Arcane::Node& first_node = (*icell).node(0);
    const Arcane::Real3& c = node_coord[first_node];
    cell_arr1_host[icell->localId()] =
      1. + Arcane::math::abs(sin(c.x + 1) * cos(c.y + 1) * sin(c.z + 2));
    cell_arr2_host[icell->localId()] =
      2. + Arcane::math::abs(cos(c.x + 2) * sin(c.y + 1) * cos(c.z + 1));
  }
  
  Kokkos::deep_copy(m_cell_arr1, cell_arr1_host);  // H2D
  Kokkos::deep_copy(m_cell_arr2, cell_arr2_host);  // H2D
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void KokkosWrapper::computeCqsAndVector()
{
  constexpr Arcane::Real k025 = 0.25;
  
  Kokkos::parallel_for(m_nb_cells, KOKKOS_CLASS_LAMBDA(const size_t& cell_i)  // allcells
	{
    //Kokkos::View<Arcane::Real3[8], Kokkos::Cuda> pos;  // pb à la desallocation
    std::array<Arcane::Real3, 8> pos;
    
    for (Arcane::Integer ii = 0; ii < 8; ++ii) {
      pos[ii] = m_node_coord_bis(m_cell_node_id(cell_i, ii));
    }
    m_cell_cqs(cell_i, 0) = -k025 * Arcane::math::cross(pos[4] - pos[3], pos[1] - pos[3]);
    m_cell_cqs(cell_i, 1) = -k025 * Arcane::math::cross(pos[0] - pos[2], pos[5] - pos[2]);
    m_cell_cqs(cell_i, 2) = -k025 * Arcane::math::cross(pos[1] - pos[3], pos[6] - pos[3]);
    m_cell_cqs(cell_i, 3) = -k025 * Arcane::math::cross(pos[7] - pos[2], pos[0] - pos[2]);
    m_cell_cqs(cell_i, 4) = -k025 * Arcane::math::cross(pos[5] - pos[7], pos[0] - pos[7]);
    m_cell_cqs(cell_i, 5) = -k025 * Arcane::math::cross(pos[1] - pos[6], pos[4] - pos[6]);
    m_cell_cqs(cell_i, 6) = -k025 * Arcane::math::cross(pos[5] - pos[2], pos[7] - pos[2]);
    m_cell_cqs(cell_i, 7) = -k025 * Arcane::math::cross(pos[6] - pos[3], pos[4] - pos[3]);
  });

  Kokkos::parallel_for(m_nb_nodes, KOKKOS_CLASS_LAMBDA(const size_t& node_i)  // allnodes
  {
    m_node_vector(node_i) = Arcane::Real3(0., 0., 0.);
  });
  
  // Useless ?
  Kokkos::fence();

  // /!\ Pas correct en parallèle à cause des accès en écriture qui peuvent être concurrents => data race condition

  // Calcul du gradient de pression
  //  Kokkos::RangePolicy<Kokkos::Cuda, int> range(0, m_nb_cells);  // equivalent
  Kokkos::parallel_for(m_nb_cells, KOKKOS_CLASS_LAMBDA(const int cell_i)  // allcells
	{
    if (m_is_active_cell(cell_i)) {
      for (int ii = 0; ii < 8; ++ii) {
        m_node_vector(m_cell_node_id(cell_i, ii)) += (m_cell_arr1(cell_i) + m_cell_arr2(cell_i)) * m_cell_cqs(cell_i, ii);
      }
    }
  });
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void KokkosWrapper::syncHostData(const Arcane::CellGroup& all_cells, const Arcane::NodeGroup& all_nodes,
                                 Arcane::VariableNodeReal3 node_vector, Arcane::VariableNodeReal3 node_coord_bis,
                                 Arcane::VariableCellArrayReal3 cell_cqs, Arcane::VariableCellReal cell_arr1,
                                 Arcane::VariableCellReal cell_arr2)
{
  // Ancienne version, empeche le basculement simple d'execution policy
/*
  Kokkos::View<Arcane::Real3*, TargetMem>::HostMirror node_vector_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), m_node_vector);
  Kokkos::View<Arcane::Real3*, TargetMem>::HostMirror node_coord_bis_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), m_node_coord_bis);
  Kokkos::View<Arcane::Real3*[8], TargetMem>::HostMirror cell_cqs_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), m_cell_cqs);
  Kokkos::View<Arcane::Real*, TargetMem>::HostMirror cell_arr1_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), m_cell_arr1);
  Kokkos::View<Arcane::Real*, TargetMem>::HostMirror cell_arr2_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), m_cell_arr2);
*/

  Kokkos::View<Arcane::Real3*, TargetMem>::HostMirror node_vector_host =
    Kokkos::create_mirror_view(m_node_vector);
  Kokkos::View<Arcane::Real3*, TargetMem>::HostMirror node_coord_bis_host =
    Kokkos::create_mirror_view(m_node_coord_bis);
  Kokkos::View<Arcane::Real3*[8], TargetMem>::HostMirror cell_cqs_host =
    Kokkos::create_mirror_view(m_cell_cqs);
  Kokkos::View<Arcane::Real*, TargetMem>::HostMirror cell_arr1_host =
    Kokkos::create_mirror_view(m_cell_arr1);
  Kokkos::View<Arcane::Real*, TargetMem>::HostMirror cell_arr2_host =
    Kokkos::create_mirror_view(m_cell_arr2);
  // copy D2H
  Kokkos::deep_copy(node_vector_host, m_node_vector);
  Kokkos::deep_copy(node_coord_bis_host, m_node_coord_bis);
  Kokkos::deep_copy(cell_cqs_host, m_cell_cqs);
  Kokkos::deep_copy(cell_arr1_host, m_cell_arr1);
  Kokkos::deep_copy(cell_arr2_host, m_cell_arr2);

  ENUMERATE_NODE(inode, all_nodes) {
    node_vector[inode] = node_vector_host(inode->localId());
    node_coord_bis[inode] = node_coord_bis_host(inode->localId());
  }

  ENUMERATE_CELL(icell, all_cells) {
    for (Arcane::Integer i = 0; i < 8; ++i)
      cell_cqs[icell][i] = cell_cqs_host(icell->localId(), i);
    cell_arr1[icell] = cell_arr1_host(icell->localId());
    cell_arr2[icell] = cell_arr2_host(icell->localId());
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void KokkosWrapper::computeCqsAndVectorV2()
{
  constexpr Arcane::Real k025 = 0.25;

  // Kokkos::parallel_for(m_nb_cells, KOKKOS_CLASS_LAMBDA(const size_t& cell_i)  // Si on lance le parallel comme ça, il force le run en cuda
  Kokkos::RangePolicy<TargetExec, int> all_cells_range(0, m_nb_cells);  // equivalent mais permet de specifier l'execution space
  Kokkos::parallel_for(all_cells_range, KOKKOS_CLASS_LAMBDA(const size_t& cell_i)  // allcells
	{
    //Kokkos::View<Arcane::Real3[8], Kokkos::Cuda> pos;  // pb à la desallocation
    std::array<Arcane::Real3, 8> pos;
    
    for (Arcane::Integer ii = 0; ii < 8; ++ii) {
      pos[ii] = m_node_coord_bis(m_cell_node_id(cell_i, ii));
    }
    m_cell_cqs(cell_i, 0) = -k025 * Arcane::math::cross(pos[4] - pos[3], pos[1] - pos[3]);
    m_cell_cqs(cell_i, 1) = -k025 * Arcane::math::cross(pos[0] - pos[2], pos[5] - pos[2]);
    m_cell_cqs(cell_i, 2) = -k025 * Arcane::math::cross(pos[1] - pos[3], pos[6] - pos[3]);
    m_cell_cqs(cell_i, 3) = -k025 * Arcane::math::cross(pos[7] - pos[2], pos[0] - pos[2]);
    m_cell_cqs(cell_i, 4) = -k025 * Arcane::math::cross(pos[5] - pos[7], pos[0] - pos[7]);
    m_cell_cqs(cell_i, 5) = -k025 * Arcane::math::cross(pos[1] - pos[6], pos[4] - pos[6]);
    m_cell_cqs(cell_i, 6) = -k025 * Arcane::math::cross(pos[5] - pos[2], pos[7] - pos[2]);
    m_cell_cqs(cell_i, 7) = -k025 * Arcane::math::cross(pos[6] - pos[3], pos[4] - pos[3]);
	});

  // Useless ? Pas sûr, si je l'enlève, les stats d'execution diffèrent
  // (mais le temps global reste quasiment identique)
  Kokkos::fence();

  // Rebelote avec les horreurs de index in cells
  // Kokkos::parallel_for(m_nb_nodes, KOKKOS_CLASS_LAMBDA(const size_t& nid)  // allnodes
  Kokkos::RangePolicy<TargetExec, int> all_nodes_range(0, m_nb_nodes);  // equivalent
  Kokkos::parallel_for(all_nodes_range, KOKKOS_CLASS_LAMBDA(const size_t& nid)  // allnodes
  {
    Arcane::Int32 first_pos = nid * 8;
    Arcane::Integer index = 0;
    Arcane::Real3 node_vec = Arcane::Real3::zero();
    for (auto cid : m_node_cell_id(nid)) {
      if (m_is_active_cell(cid)) {
        Arcane::Int16 node_index = m_node_index_in_cells(first_pos + index);
        node_vec += (m_cell_arr1(cid) + m_cell_arr2(cid)) * m_cell_cqs(cid, node_index);
      }
      ++index;
    }
    m_node_vector(nid) = node_vec;
  });

/*
// Version moins optimale car il y a une boucle de recherche d'idx des noeuds
// c'est presque 30% plus lent que de stocker la table d'idx des noeuds
  Kokkos::RangePolicy<TargetExec, int> all_nodes_range(0, m_nb_nodes);  // equivalent
  Kokkos::parallel_for(all_nodes_range, KOKKOS_CLASS_LAMBDA(const size_t& node_i)  // allnodes
  {
    m_node_vector(node_i) = Arcane::Real3(0., 0., 0.);
    
    for (auto cell_i : m_node_cell_id(node_i)) {
      if (m_is_active_cell(cell_i)) {
        int node_idx_in_cell(-1);
        for (int ii = 0; ii < 8; ++ii) {
          if (m_cell_node_id(cell_i, ii) == node_i) {
            node_idx_in_cell = ii;
            break;
          }
        }
        // if (node_idx_in_cell >7) {Kokkos::abort("DAMN");}
        // if (node_idx_in_cell ==-1) {Kokkos::abort("DAMN");}
        m_node_vector(node_i) += (m_cell_arr1(cell_i) + m_cell_arr2(cell_i)) * m_cell_cqs(cell_i, node_idx_in_cell);
      }
    }
  });
*/
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif  // WITH_KOKKOS