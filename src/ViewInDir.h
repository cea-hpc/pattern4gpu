#ifndef P4GPU_VIEW_IN_DIR_H
#define P4GPU_VIEW_IN_DIR_H

#include "arcane/VariableView.h"
#include "arcane/IMesh.h"
#include "arcane/Item.h"
#include "arcane/MeshVariableScalarRef.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*!
 * \brief Expose une vue d'une composante directionnelle d'un tableau de coordonnées 3D aux <items>
 * L'implémentation peut être AoS (Array of Structs), c'est l'implémentation native Arcane
 * Ou bien SoA (Struct of Arrays)
 *
 * Pour l'instant, sur des tableaux de Real sont implémentés
 */

//
// AOS AOS AOS
//

//! Vue d'une direction pour un array of Real3, implémentation à partir d'un Real3
template<typename ItemType>
class ViewInDirReal_AoS {
 public:
  using LocalIdType = typename ItemType::LocalIdType;

  ViewInDirReal_AoS(const MeshVariableScalarRefT<ItemType, Real3> &var_in, Integer dir) 
  : m_dir(dir), m_var_in(viewIn(var_in)) {
  }

  // En lecture
  Real operator[](LocalIdType item_id) const {
    return m_var_in[item_id][m_dir];
  }

 private:
  Integer m_dir;
  ItemVariableScalarInViewT<ItemType, Real3> m_var_in;
};

//
// SOA SOA SOA
//

// Copie des valeurs dir d'un tableau dans un tableau aux Items
template<typename ItemType>
void copy_dir(
    MeshVariableScalarRefT<ItemType, Real> &var_dir_out, 
    const MeshVariableScalarRefT<ItemType, Real3> &var_in, 
    Integer dir, IMesh *mesh) {

  // TODO : a specialiser
  abort();
}

// Specialisation Cell
template<>
void copy_dir(
    MeshVariableScalarRefT<Cell, Real> &var_dir_out, 
    const MeshVariableScalarRefT<Cell, Real3> &var_in, 
    Integer dir, IMesh *mesh) {

  ENUMERATE_CELL (cell_i, mesh->allCells()) {

    // Passage AoS => SoA
    var_dir_out[cell_i] = var_in[cell_i][dir];

  }  // Fin ENUMERATE_CELL
}

// Specialisation Node
template<>
void copy_dir(
    MeshVariableScalarRefT<Node, Real> &var_dir_out, 
    const MeshVariableScalarRefT<Node, Real3> &var_in, 
    Integer dir, IMesh *mesh) {

  ENUMERATE_NODE (node_i, mesh->allNodes()) {

    // Passage AoS => SoA
    var_dir_out[node_i] = var_in[node_i][dir];

  }  // Fin ENUMERATE_NODE
}


//! Vue d'une direction pour un array of Real3, implémentation à partir d'un tableau direct
template<typename ItemType>
class ViewInDirReal_SoA {
 public:
  using LocalIdType = typename ItemType::LocalIdType;

  ViewInDirReal_SoA(const MeshVariableScalarRefT<ItemType, Real3> &var_in, Integer dir) 
  : m_dir(dir),
  m_var_dir ( VariableBuildInfo(var_in.subDomain()->defaultMesh(),String("TemporaryCell")+var_in.name()+String::fromNumber(dir)) ),
  m_var_in (viewIn(m_var_dir)) {

    // Recopie des valeurs de var_in de la direction dir dans le tableau m_var_dir
    // Cette fonction est specialisee selon ItemType = Cell, Node, ...
    copy_dir(m_var_dir, var_in, dir, var_in.subDomain()->defaultMesh());
  }

  MeshVariableScalarRefT<ItemType, Real> &varDir() {
    return m_var_dir;
  }

  // En lecture
  Real operator[](LocalIdType item_id) const {
    return m_var_in[item_id];
  }

 private:
  Integer m_dir;
  MeshVariableScalarRefT<ItemType, Real> m_var_dir; //! Tableau aux items stockant les valeurs de la direction m_dir
  ItemVariableScalarInViewT<ItemType, Real> m_var_in;  //! Vue en lecture sur m_var_dir
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


#endif

