#ifndef CARTESIAN_CART_NODE_DIRECTION_MNG_H
#define CARTESIAN_CART_NODE_DIRECTION_MNG_H

#include "cartesian/CartTypes.h"
#include "cartesian/CartItemGroup.h"
#include "cartesian/CartesianGridT.h"

namespace Cartesian {
  
/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Meme interface que NodeDirectionMng
 */
/*---------------------------------------------------------------------------*/
class CartNodeDirectionMng {
 public:

  //! Type pour grille cartésienne sur les identifiants locaux
  using CartesianGrid = CartesianGridT<LocalIdType>;

  //! Type pour la numérotation cartésienne sur des identifiants locaux
  using CartesianNumbering = CartesianGrid::CartesianNumbering;

  //! Type tableau sur pointeurs d'ItemInternal (implémentation d'un Item)
  using ItemInternalPtr = Item::ItemInternalPtr;

 public:
  CartNodeDirectionMng(const ItemInternalPtr* internals, 
    Integer dir, const CartesianGrid &cart_grid) 
  : m_internals (internals), 
  m_dir (dir),
  m_cart_grid (cart_grid),
  m_cart_numb_node (m_cart_grid.cartNumNode()),
  m_nnodes_dir (m_cart_numb_node.nbItem3()) {

    // Pour les groupes
    for(Integer d(0) ; d < 3 ; ++d) {
      if (d == m_dir) {
        // on retire la première et la dernière couche de nodes selon m_dir
        m_inner_nodes_beg[m_dir] = 1; 
        m_inner_nodes_end[m_dir] = m_nnodes_dir[m_dir]-1;
        // on ne retient que la première couche de nodes selon m_dir
        m_prev_outer_nodes_end[m_dir] = 1; 
        // on ne retient que la dernière couche de nodes selon m_dir
        m_next_outer_nodes_beg[m_dir] = m_nnodes_dir[m_dir]-1;
      } else {
        // On prend l'intégralité du domaine selon m_dir
        m_inner_nodes_beg[d] = 0;
        m_inner_nodes_end[d] = m_nnodes_dir[d];
        // On prend l'intégralité du domaine selon m_dir
        m_prev_outer_nodes_end[d] = m_nnodes_dir[d];
        // On prend l'intégralité du domaine selon m_dir
        m_next_outer_nodes_beg[d] = 0;
      }
    }
  }

  CartNodeDirectionMng(const CartNodeDirectionMng& rhs)
  : m_internals (rhs.m_internals),
  m_dir (rhs.m_dir),
  m_cart_grid (rhs.m_cart_grid),
  m_cart_numb_node (rhs.m_cart_numb_node),
  m_nnodes_dir (rhs.m_nnodes_dir)
  {
    for(Integer d(0) ; d < 3 ; ++d) {
      m_inner_nodes_beg[d] = rhs.m_inner_nodes_beg[d];
      m_inner_nodes_end[d] = rhs.m_inner_nodes_end[d];
      m_prev_outer_nodes_end[d] = rhs.m_prev_outer_nodes_end[d];
      m_next_outer_nodes_beg[d] = rhs.m_next_outer_nodes_beg[d];
    }
  }

  //! Création d'une instance de Node à partir de son local_id
  Node toNode(LocalIdType node_id) const {
    if (ItemId::null(node_id)) {
      return Node();
    } else {
      return Node(m_internals, node_id);
    }
  }

  //! Retourne le groupe de toutes les noeuds cartésiens
  CartNodeGroup allNodes() const {
    return CartNodeGroup(m_internals, m_dir, m_cart_grid, m_cart_numb_node, {0, 0, 0}, m_nnodes_dir);
  }

  //! Groupe de tous les noeuds cartesiens internes à la direction
  CartNodeGroup innerNodes() const {
    return CartNodeGroup(m_internals, m_dir, m_cart_grid, m_cart_numb_node, m_inner_nodes_beg, m_inner_nodes_end);
  }

  //! Groupe de tous les noeuds cartesiens externes à la direction à gauche (le noeud avant est nul)
  CartNodeGroup previousOuterNodes() const {
    return CartNodeGroup(m_internals, m_dir, m_cart_grid, m_cart_numb_node, {0, 0, 0}, m_prev_outer_nodes_end);
  }

  //! Groupe de tous les noeuds cartesiens externes à la direction à droite (le noeud après est nul)
  CartNodeGroup nextOuterNodes() const {
    return CartNodeGroup(m_internals, m_dir, m_cart_grid, m_cart_numb_node, m_next_outer_nodes_beg, m_nnodes_dir);
  }

  eMeshDirection direction() const {
    return eMeshDirection(m_dir);
  }

 private:
  const ItemInternalPtr* m_internals;  //! Tableau dimensionne au nb total de Node, chaque case pointe vers un ItemInternal
  Integer m_dir;  //! Direction privilegiee

  const CartesianGrid &m_cart_grid;  //! Grille cartésienne contenant toutes les numérotations des différents items
  const CartesianNumbering &m_cart_numb_node; //! Numérotation cartésienne des noeuds (ids locaux)

  const LocalIdType3 &m_nnodes_dir; //! Nb de noeuds par direction

  LocalIdType3 m_inner_nodes_beg;  //! Triplet inclu "en bas à gauche" délimitant les noeuds intérieurs à la direction m_dir
  LocalIdType3 m_inner_nodes_end;  //! Triplet exclu "en haut à droite" délimitant les noeuds intérieurs à la direction m_dir
  LocalIdType3 m_prev_outer_nodes_end;  //! Triplet exclu "en haut à droite" délimitant les noeuds de bord au début de la direction m_dir
  LocalIdType3 m_next_outer_nodes_beg;  //! Triplet inclu "en bas à gauche" délimitant les noeuds de bord à la fin de la direction m_dir
};

}

#endif

