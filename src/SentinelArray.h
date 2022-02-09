/*******************************************************************************
 * Copyright (c) 2022 CEA
 * This program and the accompanying materials are made available under the 
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 * Contributors: Francois Letierce
 *******************************************************************************/

#include <array>
#include <initializer_list>
#include <memory>
#include <type_traits>
#include <iterator>
#include <cassert>
#include <string>
#include <stdexcept>
#include <vector>
 
#ifdef FOR_KOKKOS
#define __COMPILE_FLAG__ KOKKOS_INLINE_FUNCTION
#else
#define __COMPILE_FLAG__
#endif

/**
 * @class SentinelArray
 * @brief Simple static 1D array which adds variable end support to std::array.
 *
 *                   -----------------
 * std::array<int,4> | 0 | 1 | 2 | 3 |
 *                   -----------------
 * 
 *                                   -------------------------
 * SentinelArray<int, 4>({42, 1337}) | 42 | 1337 | xxx | xxx |
 *                                   -------------------------
 *                                                ^^^^^^^^^^^
 *                                  uninitialized and wasted memory space
 * Sentinel value is 2 is this case.
 * Using iterators, you will have access to 42 and 1337.
 * 
 * 
 * @note GPU kernel friendly. 
 *       Handy for variable data size use when vector cannot be used
 * 
 * 
 * @todo - Add tests !!!!!
 *       - Add __cplusplus macro check and constexpr methods for >c++17
 *         - Reminder :
 *           - C++ pre-C++98: __cplusplus is 1.
 *           - C++98: __cplusplus is 199711L.
 *           - C++11: __cplusplus is 201103L.
 *           - C++14: __cplusplus is 201402L.
 *           - C++17: __cplusplus is 201703L.
 *           - C++20: __cplusplus is 202002L.
 *       - Implement a better GPU macros support (__global__, KOKKOS_*, etc.)
 *       - Add override method or support for :
 *         - fill
 *         - swap
 *         - operator==
 *         - operator!=
 *         - operator<
 *         - operator<=
 *         - operator>
 *         - operator>=
 *         - operator<=>
 *         - std::get
 *         - std::swap
 *         - to_array
 *         - std::tuple_size
 * 
 * 
 * @tparam T Data type
 * @tparam N Maximum container size
 */
template <typename T, size_t N>
class SentinelArray : public std::array<T,N>
{
 public:
  /***************************** Ctors ****************************************/
  // Default Ctor
  __COMPILE_FLAG__
  SentinelArray() = default;
  
  /// Ctor from initializer list of convertible elements
  template <typename U, typename std::enable_if_t<std::is_convertible_v<T, U>>* = nullptr>
  __COMPILE_FLAG__
  SentinelArray(std::initializer_list<U> init) : m_sentinel(init.size())
  {
    assert(m_sentinel <= N);
    // std::uninitialized_copy(init.begin(), init.end(), this->data());  // host function only
    size_t i(0);
    for (auto it(init.begin()); it != init.end(); ++it, ++i)
      this->data()[i] = *it;
  }
  
  /// Ctor from begin and end iterators
  template <typename Iter, typename std::enable_if_t<std::is_same_v<typename std::iterator_traits<Iter>::iterator_category, std::random_access_iterator_tag>>* = nullptr>
  __COMPILE_FLAG__
  SentinelArray(Iter begin, Iter end) : m_sentinel(std::distance(begin, end))
  {
    assert(m_sentinel <= N);
    // std::uninitialized_copy(begin, end, this->data());  // host function only :(
    size_t i(0);
    for (Iter it(begin); it != end; ++it, ++i)
      this->data()[i] = *it;
  }

  /// Ctor from std::array of convertible elements
  template <typename U, size_t M, typename std::enable_if_t<std::is_convertible_v<T, U>>* = nullptr>
  __COMPILE_FLAG__
  SentinelArray(std::array<U,M> a) : m_sentinel(M)
  {
    static_assert(M <= N);
    // std::uninitialized_copy(a.begin(), a.end(), this->data());  // host function only
    size_t i(0);
    for (auto it(a.begin()); it != a.end(); ++it, ++i)
      this->data()[i] = *it;
  }
  
  /// Ctor from std::vector of convertible elements
  template <typename U, typename std::enable_if_t<std::is_convertible_v<T, U>>* = nullptr>
  __COMPILE_FLAG__
  SentinelArray(std::vector<U> v) : m_sentinel(v.size())
  {
    assert(m_sentinel <= N);
    // std::uninitialized_copy(a.begin(), a.end(), this->data());  // host function only
    size_t i(0);
    for (auto it(v.begin()); it != v.end(); ++it, ++i)
      this->data()[i] = *it;
  }
  
  /***************************** Copy operator ********************************/
  /// Copy operator from an initializer list of convertible elements
  template <typename U, typename std::enable_if_t<std::is_convertible_v<T, U>>* = nullptr>
  __COMPILE_FLAG__
  SentinelArray& operator=(std::initializer_list<U> init)
  {
    assert(init.size() <= m_sentinel);
    // std::uninitialized_copy(init.begin(), init.end(), this->data());  // host function only
    size_t i(0);
    for (auto it(init.begin()); it != init.end(); ++it, ++i)
      this->data()[i] = *it;
    return *this;
  }
  
  /// Copy operator from std::array of convertible elements
  template <typename U, size_t M, typename std::enable_if_t<std::is_convertible_v<T, U>>* = nullptr>
  __COMPILE_FLAG__
  SentinelArray& operator=(std::array<U,M> a)
  {
    assert(M <= m_sentinel);
    // std::uninitialized_copy(a.begin(), a.end(), this->data());  // host function only
    size_t i(0);
    for (auto it(a.begin()); it != a.end(); ++it, ++i)
      this->data()[i] = *it;
    return *this;
  }
  
  /// Copy operator from std::vector of convertible elements
  template <typename U, typename std::enable_if_t<std::is_convertible_v<T, U>>* = nullptr>
  __COMPILE_FLAG__
  SentinelArray& operator=(std::vector<U> v)
  {
    assert(v.size() <= m_sentinel);
    // std::uninitialized_copy(a.begin(), a.end(), this->data());  // host function only
    size_t i(0);
    for (auto it(v.begin()); it != v.end(); ++it, ++i)
      this->data()[i] = *it;
    return *this;
  }  
  
  /***************************** End iterators ********************************/
  /// End iterator (override to support sentinel value)
  __COMPILE_FLAG__
  typename std::array<T,N>::iterator end() noexcept
  {
      return typename std::array<T,N>::iterator(this->data()+m_sentinel);
  }
  /// Const end iterator (override to support sentinel value)
  __COMPILE_FLAG__
  typename std::array<T,N>::const_iterator end() const noexcept
  {
    return typename std::array<T,N>::const_iterator(this->data()+m_sentinel);
  }
  /// Cend iterator (override to support sentinel value)
  __COMPILE_FLAG__
  typename std::array<T,N>::const_iterator cend() const noexcept
  {
    return typename std::array<T,N>::const_iterator(this->data()+m_sentinel);
  }

  /***************************** Reverse begin iterators **********************/
  /// Rbegin iterator (override to support sentinel value)
  __COMPILE_FLAG__
  typename std::array<T,N>::reverse_iterator rbegin() noexcept
  {
    return typename std::array<T,N>::reverse_iterator(this->data()+m_sentinel);
  }
  /// Const rbegin iterator (override to support sentinel value)
  __COMPILE_FLAG__
  typename std::array<T,N>::const_reverse_iterator rbegin() const noexcept
  {
    return typename std::array<T,N>::const_reverse_iterator(this->data()+m_sentinel);
  }
  /// Crbegin iterator (override to support sentinel value)
  __COMPILE_FLAG__
  typename std::array<T,N>::const_reverse_iterator crbegin() const noexcept
  {
    return typename std::array<T,N>::const_reverse_iterator(this->data()+m_sentinel);
  }

  /***************************** Size operator ********************************/
  /// Size operator (override to support sentinel value)
  __COMPILE_FLAG__
  typename std::array<T,N>::size_type size() const noexcept
  {
    return m_sentinel;
  }
  
  /***************************** Bound checking accessors *********************/
  /// Boundaries checking accessor (override to support sentinel value)
  typename std::array<T,N>::reference at(typename std::array<T,N>::size_type pos)
  {
      if (!(pos < m_sentinel))
          throw std::out_of_range("array::at: pos (which is " + std::to_string(pos) + ") >= m_sentinel (which is " + std::to_string(m_sentinel) + ")");
      return (*this)[pos];
  }
  /// Boundaries checking const accessor (override to support sentinel value)
  typename std::array<T,N>::const_reference at(typename std::array<T,N>::size_type pos) const
  {
      if (!(pos < m_sentinel))
          throw std::out_of_range("array::at: pos (which is " + std::to_string(pos) + ") >= m_sentinel (which is " + std::to_string(m_sentinel) + ")");
      return (*this)[pos];
  }
  
  /***************************** Back accessor ********************************/
  /// Last value accessor (override to support sentinel value)
  typename std::array<T,N>::reference back() noexcept
  {
    return m_sentinel ? *(end() - 1) : *end();
  }
  /// Last value const accessor (override to support sentinel value)
  /* FIXME: cannot override
  typename std::array<T,N>::const_reference back() noexcept
  {
    return m_sentinel ? *(end() - 1) : *end();
  }
  */
  
  /***************************** Sentinel setter ******************************/
  /// Sentinel setter. Defines artificial new end of array.
  __COMPILE_FLAG__
  void setSentinel(typename std::array<T,N>::size_type s)
  {
    m_sentinel = s;
  }

  /***************************** Private attribute ****************************/
 private:
  typename std::array<T,N>::size_type m_sentinel;
};
