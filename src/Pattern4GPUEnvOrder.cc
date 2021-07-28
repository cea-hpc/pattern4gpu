#include "Pattern4GPUModule.h"
#include "cartesian/CartesianMeshProperties.h"
#include "cartesian/CartesianMeshT.h"

#include "arcane/VariableView.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*!
 * \brief Prépare l'utilisation de la méthode DetEnvOrder 
 *  pour la boucle de calcul
 */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
initEnvOrder() {
  // Peut-être allouer m_cell_env_order ici pour éviter de le faire à chaque
  // appel dans la boucle de calcul
}


typedef std::pair<Real, Integer> mypair;
bool operator<(const mypair& p1, const mypair& p2)
{
  return p1.first < p2.first;
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Détermine un ordre de traitement des environnements par maille mixte
 *  à partir des fractions volumiques des environnements des mailles voisines
 */
/*---------------------------------------------------------------------------*/
template<typename CartesianMeshT>
void Pattern4GPUModule::
_detEnvOrder()
{
  ARCANE_ASSERT(m_cartesian_mesh!=nullptr, ("Maillage cartésien inexistant ou non initialisé"));
  CartesianMeshT cart_mesh_t(m_cartesian_mesh);

  using NeighCellsType = typename CartesianMeshT::template NeighCells<3>;
  using ArrayCellType = typename NeighCellsType::ArrayCellType;
  using ConnectivityCellNode = typename CartesianMeshT::ConnectivityCellNode;

  NeighCellsType&& neigh_cells = cart_mesh_t.template neighCells<3>();
  ConnectivityCellNode&& conn_cn = cart_mesh_t.connectivityCellNode();
  auto&& cell_dm = cart_mesh_t.cellDirection(0);
  auto&& all_cells = cell_dm.allCells();

  // Pour le parcours des milieux
  CellToAllEnvCellConverter allenvcell_converter(m_mesh_material_mng);

  m_cell_env_order.resize(m_mesh_material_mng->environments().size());
  m_cell_env_order.fill(-1);

  // Pre-allocation de tableaux temporaires dimensionnes au nb max d'environnements possibles par maille
  const Integer max_nb_env = m_mesh_material_mng->environments().size();

  UniqueArray<bool> has_null_volume(max_nb_env);
  UniqueArray<Integer> env_order(max_nb_env);
  UniqueArray<Real> x_alpha(max_nb_env), y_alpha(max_nb_env);
  UniqueArray<Real> proj_x_alpha(max_nb_env), proj_y_alpha(max_nb_env);
  UniqueArray<Real> alpha_total(max_nb_env);
  UniqueArray<Integer> nb_pure(max_nb_env);
  UniqueArray<Real> alpha(max_nb_env);
  UniqueArray<EnvCell> env_cell_array;
  UniqueArray<mypair> dist(max_nb_env);

  env_cell_array.reserve(max_nb_env);

  const Integer stencil_sz = NeighCellsType::stencil_sz;
  ArrayCellType _cells(stencil_sz);

  UniqueArray<Real> z_alpha(max_nb_env), proj_z_alpha(max_nb_env);
  UniqueArray<mypair> dd(3);

  const Integer nb_node=conn_cn.nbNode();

  auto in_node_coord = viewIn(m_node_coord_bis);

  ENUMERATE_AUTO_CELL (cell_i, all_cells) {
    const Cell& cell(*cell_i);
    AllEnvCell allenvcell = allenvcell_converter[cell];
    const Integer nb_env(allenvcell.nbEnvironment());

    if (nb_env > 1) {
      // Pour identifier les environnements de volume nul
      has_null_volume.resize(nb_env);
      for (Integer env_i(0) ; env_i < nb_env ; env_i++) {
        has_null_volume[env_i] = false;
      }

      // Mailles voisines avec conditions aux limites
      neigh_cells.neighCellsBC(cell_i, _cells);

      // On récupère les milieux présents dans la maille

      env_cell_array.clear();
      ENUMERATE_CELL_ENVCELL (envcell_i, allenvcell) {
        env_cell_array.add(*envcell_i);
      }

      // Pour chacun de ces milieux, on calcule un barycentre
      for (Integer mat_i = 0 ; mat_i < nb_env ; ++mat_i) {
        EnvCell envcell = env_cell_array[mat_i];
        // On récupère l'environnement en question
        IMeshEnvironment* env(envcell.environment());

        if (m_volume[envcell] == 0) {
          has_null_volume[mat_i] = true;
        }

        Real cumul_alpha = 0.;
        alpha_total[mat_i] = 0.;
        x_alpha[mat_i] = 0.;
        y_alpha[mat_i] = 0.;
        z_alpha[mat_i] = 0.;
        nb_pure[mat_i] = 0;

        if (!has_null_volume[mat_i]) {
          // On parcourt le stencil de mailles
          for (Integer jj(0); jj < _cells.size(); jj ++) {
            auto cell_jj = _cells[jj];
            // On calcule le centre de la maille
            Real x_jj = 0., y_jj = 0., z_jj = 0.;
            // Centre de la maille consideree
            const auto&& cell_conn = conn_cn.cellConnectivity(cell_jj);
            for(Integer inode(0) ; inode < nb_node ; inode++) {
              const auto&& node_jj(cell_conn.node(inode));
              x_jj += in_node_coord[node_jj].x / ((Real) nb_node);
              y_jj += in_node_coord[node_jj].y / ((Real) nb_node);
              z_jj += in_node_coord[node_jj].z / ((Real) nb_node);
            }

            // On recherche l'environnement envcell dans cette maille
            AllEnvCell allenvcell_jj = allenvcell_converter[CellLocalId(cell_jj)];
            EnvCell env_cell_jj = env->findEnvCell(allenvcell_jj);
            // S'il est present
            if (!env_cell_jj.null()) {
              x_alpha[mat_i] += m_frac_vol[env_cell_jj] * x_jj;
              y_alpha[mat_i] += m_frac_vol[env_cell_jj] * y_jj;
              z_alpha[mat_i] += m_frac_vol[env_cell_jj] * z_jj;
              cumul_alpha += m_frac_vol[env_cell_jj];

              if (allenvcell_jj.nbEnvironment() == 1) {
                nb_pure[mat_i]++;
              }
            }
          }
          // On divise les centroides par le cumul des fractions de volume
          x_alpha[mat_i] /= cumul_alpha;
          y_alpha[mat_i] /= cumul_alpha;
          z_alpha[mat_i] /= cumul_alpha;
          alpha_total[mat_i] = cumul_alpha;
        }
      }

      Real x_bar(0.0), y_bar(0.0), z_bar(0.0);
      for (Integer mat_i(0) ; mat_i < nb_env ; ++mat_i) {
        alpha[mat_i] = 1. / ((Real) nb_env);
        x_bar += x_alpha[mat_i] * alpha[mat_i];
        y_bar += y_alpha[mat_i] * alpha[mat_i];
        z_bar += z_alpha[mat_i] * alpha[mat_i];
      }

      // On forme les composantes de la matrice re covariance
      Real axx(0.0), ayy(0.0), azz(0.0), axy(0.0), axz(0.0), ayz(0.0);
      for (Integer mat_i(0) ; mat_i < nb_env ; ++mat_i) {
        axx += (x_alpha[mat_i] - x_bar) * (x_alpha[mat_i] - x_bar) * alpha[mat_i];
        ayy += (y_alpha[mat_i] - y_bar) * (y_alpha[mat_i] - y_bar) * alpha[mat_i];
        azz += (z_alpha[mat_i] - z_bar) * (z_alpha[mat_i] - z_bar) * alpha[mat_i];
        axy += (x_alpha[mat_i] - x_bar) * (y_alpha[mat_i] - y_bar) * alpha[mat_i];
        axz += (x_alpha[mat_i] - x_bar) * (z_alpha[mat_i] - z_bar) * alpha[mat_i];
        ayz += (y_alpha[mat_i] - y_bar) * (z_alpha[mat_i] - z_bar) * alpha[mat_i];
      }
      Real nx(0.), ny(0.), nz(0.);

      if (nb_env == 2) {
        // Les points sont alignes
        nx = x_alpha[1] - x_alpha[0];
        ny = y_alpha[1] - y_alpha[0];
        nz = z_alpha[1] - z_alpha[0];
        Real norm = math::sqrt(nx * nx + ny * ny + nz * nz);
        if (norm != 0.) {
          nx /= norm, ny /= norm, nz /= norm;
        } else {
          for (Integer mat_i(0) ; mat_i < nb_env ; ++mat_i) {
            env_order[nb_env -1 - mat_i] = env_cell_array[mat_i].environmentId();
          }
        }
      } else {
        Real3x3 a, rotation_matrix;
        Real3 eigen_values;

        a[0][0] = axx, a[0][1] = axy, a[0][2] = axz;
        a[1][0] = axy, a[1][1] = ayy, a[1][2] = ayz;
        a[2][0] = axz, a[2][1] = ayz, a[2][2] = azz;

        // TODO : diogonaliser a => rotation_matrix, eigen_values
        {
          rotation_matrix=a;
          eigen_values=Real3(axx,axy,axz);
        }

        // On classe les valeurs propres
        for (Integer kk = 0 ; kk < 3 ; ++kk) {
          dd[kk] = mypair(eigen_values[kk], kk);
        }

        std::sort(dd.begin(), dd.end());
        Integer i_lambda = dd[0].second;
        nx = rotation_matrix[0][i_lambda];
        ny = rotation_matrix[1][i_lambda];
        nz = rotation_matrix[2][i_lambda];
        Real norm = math::sqrt(nx * nx + ny * ny + nz * nz);
        if (norm != 0.) {
          nx /= norm, ny /= norm, nz /= norm;
        } else {
          for (Integer mat_i(0) ; mat_i < nb_env ; ++mat_i) {
            env_order[nb_env -1 - mat_i] = env_cell_array[mat_i].environmentId();
          }
        }
      }

      for (Integer mat_i(0) ; mat_i < nb_env ; ++mat_i) {
        Real xn = nx * (x_alpha[mat_i] - x_bar) + ny * (y_alpha[mat_i] - y_bar) +
          nz * (z_alpha[mat_i] - z_bar);
        proj_x_alpha[mat_i] = x_bar +nx * xn;
        proj_y_alpha[mat_i] = y_bar +ny * xn;
        proj_z_alpha[mat_i] = z_bar +nz * xn;
      }

      for (Integer mat_i(0) ; mat_i < nb_env ; ++mat_i) {
        dist[mat_i] = mypair(nx * (proj_x_alpha[mat_i] - x_bar) + ny * (proj_y_alpha[mat_i] - y_bar) +
            nz * (proj_z_alpha[mat_i] - z_bar), mat_i);
      }

      for (Integer mat_i(0) ; mat_i < nb_env ; ++mat_i) {
        env_order[mat_i] = env_cell_array[dist[mat_i].second].environmentId();
      }

      // -----------------------------------------------------------
      // Pour définir l'ordre :
      // On prend come premier milieu celui qui a le plus de mailles pures dans le stencil
      // En cas d'égalité, celui qui a le plus grand volume dans le stencil
      if ((nb_pure[nb_env-1] == nb_pure[0] &&
            alpha_total[nb_env-1] > alpha_total[0]) ||
          (nb_pure[nb_env-1] > nb_pure[0]) ) {
        std::reverse(env_order.begin(), env_order.begin()+nb_env);
      }

      // Les env. de volume vide (i.e. à éliminer) sont mis à la fin
      // -----------------------------------------------------------

      // il faut recalculer ici has_null_volume ici car la correspondance
      // envcell.environmentId() <-> env_order[mat_i]) a changé depuis le calcul initial plus haut
      for (Integer mat_i(0) ; mat_i < nb_env; ++mat_i) {
        has_null_volume[mat_i] = false;
        ENUMERATE_CELL_ENVCELL (envcell_i, allenvcell) {
          const EnvCell& envcell = *envcell_i;
          if (envcell.environmentId() == env_order[mat_i]) {
            if (m_volume[envcell] == 0) {
              has_null_volume[mat_i] = true;
            }
          }
        }
      }
      for (Integer mat_i(0) ; mat_i < nb_env - 1 ; ++mat_i) {
        if (has_null_volume[mat_i]) {
          pinfo() << " cell = " << (*cell_i).uniqueId() << "mat_i = " << mat_i
            << "env_order[mat_i]=" << env_order[mat_i];
          Integer sav_id(env_order[mat_i]);
          for (Integer mat_j(mat_i) ; mat_j < nb_env - 1; ++mat_j) {
            env_order[mat_j] = env_order[mat_j + 1];
            has_null_volume[mat_j] = has_null_volume[mat_j + 1];
          }
          env_order[nb_env -1] = sav_id;
          has_null_volume[nb_env -1] = true;
        }
      }

      // -----------------------------------------------------------
      // On recopie env_order[0:nb_env) dans m_cell_env_order[*cell_i]
      m_cell_env_order[*cell_i].copy(env_order.subConstView(0, nb_env));
    }  // nb_env > 1
  }  // ENUMERATE_CELL
}

void Pattern4GPUModule::
detEnvOrder() {
  PROF_ACC_BEGIN(__FUNCTION__);
  ARCANE_ASSERT(subDomain()->defaultMesh()->dimension()==3, ("Seul le 3D est supporté"));

  Cartesian::CartesianMeshProperties cart_mesh_prop(mesh());
  //bool is_cartesian_mesh = cart_mesh_prop.isPureCartesianMesh(); 
  bool is_cartesian_mesh = false; 

  if (is_cartesian_mesh) {
    _detEnvOrder<CartCartesianMeshT>();
  } else {
    _detEnvOrder<UnstructCartesianMeshT>();
  }
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
