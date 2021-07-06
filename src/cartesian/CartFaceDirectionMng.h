#ifndef CARTESIAN_CART_FACE_DIRECTION_MNG_H
#define CARTESIAN_CART_FACE_DIRECTION_MNG_H

#include "cartesian/CartTypes.h"
#include "cartesian/CartItemGroup.h"
#include "cartesian/CartesianGridT.h"

namespace Cartesian {
  
/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Meme interface que DirFace mais en cartésien, 
 * permet de passer de l'identifiant de la face aux mailles adjacentes
 * selon la direction orthogonale à la face
 */
/*---------------------------------------------------------------------------*/
class CartDirFace {
 public:
  using CellType = CellLocalId;

  //! Type pour la numérotation cartésienne locale aux mailles
  using CartesianNumbering = CartesianNumberingT<LocalIdType>;

 public:
  CartDirFace(LocalIdType face_idx_dir, LocalIdType nfacesm1_dir, LocalIdType next_cell_id, LocalIdType delta_prev_cell)
  : m_face_idx_dir (face_idx_dir),
  m_nfacesm1_dir (nfacesm1_dir),
  m_next_cell_id (next_cell_id),
  m_delta_prev_cell (delta_prev_cell) {
  }

  //! Identifiant de la maille adjacente à la face et juste avant cette face selon la direction
  CellType previousCell() const {
    return CellType(m_face_idx_dir == 0 ? -1 : /*m_face_id + m_delta_conv*/m_next_cell_id - m_delta_prev_cell);
  }

  //! Identifiant de la maille adjacente à la face et juste après cette face selon la direction
  CellType nextCell() const {
    // La maille suivante se repère par les mêmes indices que la face
    return CellType(m_face_idx_dir == m_nfacesm1_dir ? -1 : m_next_cell_id /*m_face_id + m_delta_conv*/);
  }

 private:

  LocalIdType m_face_idx_dir; //! m_face_ijk[m_dir] pré-calculé
  LocalIdType m_nfacesm1_dir; //! Nb de faces -1 dans la direction m_dir
  LocalIdType m_next_cell_id; //! L'identifiant de la maille qui suit la face
  LocalIdType m_delta_prev_cell; //! Delta à retirer pour passer de la maille suivante à précédente
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Meme interface que FaceDirectionMng
 */
/*---------------------------------------------------------------------------*/
class CartFaceDirectionMng {
 public:
  using FaceType = FaceLocalId;
  using NodeType = NodeLocalId;
  using FaceEnumeratorType = CartFaceEnumerator;

  //! Type pour une grille cartésienne avec identifiants locaux
  using CartesianGrid = CartesianGridT<LocalIdType>;

  //! Type pour la numérotation cartésienne sur des identifiants locaux
  using CartesianNumbering = typename CartesianGrid::CartesianNumbering;

  //! Type tableau sur pointeurs d'ItemInternal (implémentation d'un Item)
  using ItemInternalPtr = Item::ItemInternalPtr;

 public:
  CartFaceDirectionMng(const ItemInternalPtr* internals,
    Integer dir, const CartesianGrid &cart_grid)
  : m_internals (internals), 
  m_dir (dir),
  m_cart_grid (cart_grid),
  m_cart_cell_numbering (m_cart_grid.cartNumCell()),
  m_cart_face_numbering (m_cart_grid.cartNumFace(m_dir)),
  m_nfaces_dir (m_cart_face_numbering.nbItem3()),
  m_nfacesm1_dir (m_nfaces_dir[m_dir]-1),
  m_delta_prev_cell (m_cart_cell_numbering.deltaDir(m_dir)) {

    for(Integer d(0) ; d < 3 ; ++d) {
      if (d == m_dir) {
        // on retire la première et la dernière couche de faces selon m_dir
        m_inner_faces_beg[m_dir] = 1; 
        m_inner_faces_end[m_dir] = m_nfacesm1_dir;
        // on ne retient que la première couche de faces selon m_dir
        m_prev_outer_faces_end[m_dir] = 1; 
        // on ne retient que la dernière couche de faces selon m_dir
        m_next_outer_faces_beg[m_dir] = m_nfacesm1_dir;
      } else {
        // On prend l'intégralité du domaine selon m_dir
        m_inner_faces_beg[d] = 0;
        m_inner_faces_end[d] = m_nfaces_dir[d];
        // On prend l'intégralité du domaine selon m_dir
        m_prev_outer_faces_end[d] = m_nfaces_dir[d];
        // On prend l'intégralité du domaine selon m_dir
        m_next_outer_faces_beg[d] = 0;
      }
    }
  }

  CartFaceDirectionMng(const CartFaceDirectionMng& rhs)
  : m_internals (rhs.m_internals),
  m_dir (rhs.m_dir),
  m_cart_grid (rhs.m_cart_grid),
  m_cart_cell_numbering (rhs.m_cart_cell_numbering),
  m_cart_face_numbering (rhs.m_cart_face_numbering),
  m_nfaces_dir (rhs.m_nfaces_dir),
  m_nfacesm1_dir (rhs.m_nfacesm1_dir),
  m_delta_prev_cell (rhs.m_delta_prev_cell)
  {
    for(Integer d(0) ; d < 3 ; ++d) {
      m_inner_faces_beg[d] = rhs.m_inner_faces_beg[d];
      m_inner_faces_end[d] = rhs.m_inner_faces_end[d];
      m_prev_outer_faces_end[d] = rhs.m_prev_outer_faces_end[d];
      m_next_outer_faces_beg[d] = rhs.m_next_outer_faces_beg[d];
    }
  }

  //! Items adjacents à la face f
  CartDirFace face(const FaceEnumeratorType &f) const {
    return CartDirFace(f.itemIdxDir(m_dir), m_nfacesm1_dir, f.localIdConv(m_type_cell), m_delta_prev_cell);
  }

  //! Items adjacents à la face f
  CartDirFace operator[](const FaceEnumeratorType &f) const {
    return CartDirFace(f.itemIdxDir(m_dir), m_nfacesm1_dir, f.localIdConv(m_type_cell), m_delta_prev_cell);
  }

  //! Retourne le groupe de toutes les faces cartesiennes
  CartFaceGroup allFaces() const {
    return CartFaceGroup(m_internals, m_dir, m_cart_grid, m_cart_face_numbering, {0, 0, 0}, m_nfaces_dir);
  }

  //! Groupe de toutes les faces cartesiennes internes à la direction
  CartFaceGroup innerFaces() const {
    return CartFaceGroup(m_internals, m_dir, m_cart_grid, m_cart_face_numbering, m_inner_faces_beg, m_inner_faces_end);
  }

  //! Groupe de toutes les faces cartesiennes externes à la direction à gauche (la face "avant" est nulle)
  CartFaceGroup previousOuterFaces() const {
    return CartFaceGroup(m_internals, m_dir, m_cart_grid, m_cart_face_numbering, {0, 0, 0}, m_prev_outer_faces_end);
  }

  //! Groupe de toutes les faces cartesienness externes à la direction à droite (la face "après" est nulle)
  CartFaceGroup nextOuterFaces() const {
    return CartFaceGroup(m_internals, m_dir, m_cart_grid, m_cart_face_numbering, m_next_outer_faces_beg, m_nfaces_dir);
  }

  eMeshDirection direction() const {
    return eMeshDirection(m_dir);
  }

 private:
  
 private:
  const ItemInternalPtr* m_internals;  //! Tableau dimensionne au nb total de Face, chaque case pointe vers un ItemInternal
  Integer m_dir;  //! Direction des faces

  const CartesianGrid &m_cart_grid;  //! Grille cartésienne en numérotation locale (cohérence entre mailles, faces, noeuds)
  const CartesianNumbering &m_cart_cell_numbering;  //! Permet de numeroter les mailles a partir de (i,j,k) et reciproquement
  const CartesianNumbering &m_cart_face_numbering;  //! Permet de numeroter les faces de direction m_dir a partir de (i,j,k) et reciproquement

  const LocalIdType3 &m_nfaces_dir;  //! Nb de faces (de direction m_dir) par direction
  LocalIdType m_nfacesm1_dir;  //! Nb de faces-1 dans la direction m_dir (ie m_nfaces_dir[m_dir]-1)
  LocalIdType m_delta_prev_cell;  //! Delta pour passer de next_cell_id à prev_cell_id dans la direction m_dir

  LocalIdType3 m_inner_faces_beg;  //! Triplet inclu "en bas à gauche" délimitant les faces intérieures à la direction m_dir
  LocalIdType3 m_inner_faces_end;  //! Triplet exclu "en haut à droite" délimitant les faces intérieures à la direction m_dir
  LocalIdType3 m_prev_outer_faces_end;  //! Triplet exclu "en haut à droite" délimitant les faces de bord au début de la direction m_dir
  LocalIdType3 m_next_outer_faces_beg;  //! Triplet inclu "en bas à gauche" délimitant les faces de bord à la fin de la direction m_dir

  constexpr static Cell* m_type_cell = nullptr;  //! Permet d'utiliser la méthode localIdConv sur les Cell
};

}

#endif

