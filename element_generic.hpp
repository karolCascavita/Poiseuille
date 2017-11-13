/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#pragma once

#include <vector>

#include "mesh/mesh.hpp"
#include "mesh/mesh_storage.hpp"
#include "mesh/ident.hpp"
#include "mesh/point.hpp"

namespace disk {

template<size_t DIM, size_t CODIM>
class generic_element;

/* type traits for the elements */
template<typename T>
struct generic_element_traits;

template<size_t DIM, size_t CODIM>
struct generic_element_traits<generic_element<DIM, CODIM>>
{
    static_assert(CODIM < DIM, "generic_element_traits: CODIM must be less than DIM");

    typedef generic_element<DIM, CODIM+1>                           subelement_type;
    typedef identifier<generic_element<DIM,CODIM>, ident_impl_t, 0> id_type;
    static const size_t dimension = DIM;
    static const size_t codimension = CODIM;
};

template<size_t DIM>
struct generic_element_traits<generic_element<DIM, DIM>>
{
    typedef identifier<generic_element<DIM,DIM>, ident_impl_t, 0>   id_type;
    static const size_t dimension = DIM;
    static const size_t codimension = DIM;
};


/* element base class (for subelements) */
template<typename T>
class generic_element_base
{
    typedef typename generic_element_traits<T>::id_type                     id_type;
    typedef typename generic_element_traits<T>::subelement_type::id_type    sub_id_type;
    typedef point_identifier<generic_element_traits<T>::dimension>          point_id_type;
protected:
    id_type                             m_element_id;
    std::vector<sub_id_type>            m_sids_ptrs;
    std::vector<point_id_type>          m_pts_ptrs;
public:
    typedef typename std::vector<sub_id_type>::iterator         sid_iterator;
    typedef typename std::vector<sub_id_type>::const_iterator   const_sid_iterator;

    generic_element_base()
    {
        //std::cout << "/* generic_element_base constructor void */" << std::endl;
    }

    generic_element_base(std::initializer_list<sub_id_type> l)
        : m_sids_ptrs(l)
    {
        //std::cout << "/* generic_element_base constructor not void */" << std::endl;
    }

    explicit generic_element_base(const std::vector<sub_id_type>& sids_ptrs)
        : m_sids_ptrs(sids_ptrs)
    {
        //std::cout << "/* generic_element_base explicit constructor */" << std::endl;
    }

    template<typename Itor>
    void set_subelement_ids(const Itor e_begin, const Itor e_end)
    {
        m_sids_ptrs = std::vector<sub_id_type>(e_begin, e_end);
        //std::cout << "/* set_subelement_ids */" << std::endl;
    }



    /* All this point-related stuff is simply horrible, because points should
     * be automatically derived from lower-dimensional elements. However,
     * this is not possible because points must be stored in a specific order,
     * otherwise geometry-related computations fail.
     *
     * No, in 2D they won't fail! But what about 3D?
     */
    template<typename Itor>
    void set_point_ids(const Itor p_begin, const Itor p_end)
    {
        m_pts_ptrs = std::vector<point_id_type>(p_begin, p_end);
        //std::cout << "/* set_point_ids function */" << std::endl;
    }

    std::vector<point_id_type>
    point_ids(void) const
    {
        assert( m_pts_ptrs.size() != 0 );
        return m_pts_ptrs;
        //std::cout << "/* point_ids function*/" << std::endl;
    }

    std::vector<sub_id_type>
    faces_ids(void) const
    {
        assert( m_sids_ptrs.size() != 0 );
        return m_sids_ptrs;
        //std::cout << "/* point_ids function*/" << std::endl;
    }

    std::vector<sub_id_type>
    get_m_sids_ptrs(void) const
    {
        assert( m_sids_ptrs.size() != 0 );
        return m_sids_ptrs;
        //std::cout << "/* point_ids function*/" << std::endl;
    }

    void set_element_id(const id_type id)
    {
        m_element_id = id;
    }

    id_type
    get_id(void) const
    {
        //assert( m_element_id < 0 );
        return m_element_id;
    }


    sid_iterator        subelement_id_begin() { return m_sids_ptrs.begin(); }
    sid_iterator        subelement_id_end()   { return m_sids_ptrs.end(); }
    const_sid_iterator  subelement_id_begin() const { return m_sids_ptrs.begin(); }
    const_sid_iterator  subelement_id_end()   const { return m_sids_ptrs.end(); }

    size_t              subelement_size() const { return m_sids_ptrs.size(); }
};

/* element class, CODIM < DIM case */
template<size_t DIM, size_t CODIM>
class generic_element : public generic_element_base<generic_element<DIM, CODIM>>
{
    static_assert(CODIM <= DIM, "generic_element: CODIM > DIM does not make sense");

public:
    typedef typename generic_element_traits<generic_element>::subelement_type  subelement_type;
    typedef typename generic_element_traits<generic_element>::id_type          id_type;

    generic_element()
    {
        //std::cout << "/* generic_element      constructor void CODIM*/" << std::endl;

    }
    generic_element(std::initializer_list<typename subelement_type::id_type> l)
        : generic_element_base<generic_element<DIM, CODIM>>(l)
    {
        //std::cout << "/* generic_element      constructor not void CODIM*/" << std::endl;
    }

    explicit generic_element(const std::vector<typename subelement_type::id_type>& sids_ptrs)
        : generic_element_base<generic_element<DIM, CODIM>>(sids_ptrs)
    {
        //std::cout << "/* generic_element      constructor explicit CODIM*/" << std::endl;
    }
    bool operator<(const generic_element& other) const
    {
        return this->m_sids_ptrs < other.m_sids_ptrs;
    }

    bool operator==(const generic_element& other) const
    {
        return std::equal(this->m_sids_ptrs.begin(), this->m_sids_ptrs.end(),
                          other.m_sids_ptrs.begin());
    }
};

/* element class, CODIM == DIM case */
template<size_t DIM>
class generic_element<DIM, DIM>
{
    point_identifier<DIM>    m_assoc_point;

public:
    typedef typename generic_element_traits<generic_element>::id_type          id_type;

    generic_element()
    {
        //std::cout << "/* generic_element      constructor void DIM*/" << std::endl;

    }

    explicit generic_element(point_identifier<DIM> assoc_point)
    {
        m_assoc_point = assoc_point;
        //std::cout << "/* generic_element      constructor explicit DIM*/" << std::endl;

    }
    std::vector<point_identifier<DIM>>
    point_ids(void) const
    {
        return std::vector<point_identifier<DIM>>{{m_assoc_point}};
    }

    bool operator<(const generic_element& other) const
    {
        return this->m_assoc_point < other.m_assoc_point;
    }

    bool operator==(const generic_element& other) const
    {
        return this->m_assoc_point == other.m_assoc_point;
    }
};

/* Output streaming operator for elements */
template<size_t DIM, size_t CODIM>
std::ostream&
operator<<(std::ostream& os, const generic_element<DIM, CODIM>& e)
{
    os << "generic_element<" << DIM << "," << CODIM << ">: ";
    for (auto itor = e.subelement_id_begin();
         itor != e.subelement_id_end();
         itor++)
        os << *itor << " ";

    return os;
}

template<size_t DIM>
std::ostream&
operator<<(std::ostream& os, const generic_element<DIM, DIM>& e)
{
    os << "element<" << DIM << "," << DIM << ">: ";
    return os;
}

} // namespace disk
