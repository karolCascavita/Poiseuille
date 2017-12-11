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

 /*
  * This source file is part of EMT, the ElectroMagneticTool.
  *
  * Copyright (C) 2013-2015, Matteo Cicuttin - matteo.cicuttin@uniud.it
  * Department of Electrical Engineering, University of Udine
  * All rights reserved.
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions are met:
  *     * Redistributions of source code must retain the above copyright
  *       notice, this list of conditions and the following disclaimer.
  *     * Redistributions in binary form must reproduce the above copyright
  *       notice, this list of conditions and the following disclaimer in the
  *       documentation and/or other materials provided with the distribution.
  *     * Neither the name of the University of Udine nor the
  *       names of its contributors may be used to endorse or promote products
  *       derived from this software without specific prior written permission.
  *
  * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
  * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
  * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  */


#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include <numeric>
#include <cassert>
#include <iterator>

#include "ident.hpp"
#include "point.hpp"

namespace disk {

namespace priv {
    template<typename Iterator, typename Predicate>
    class filter_iterator
    {
        Predicate       predicate_;
        Iterator        itor_{}, end_{};

        void advance(void) { itor_++; }

        void find_next(void)
        {
            while ( itor_ != end_ )
            {
                if ( predicate_(*itor_) )
                    return;
                advance();
            }
        }

    public:
        using value_type = typename std::iterator_traits<Iterator>::value_type;
        using reference = typename std::iterator_traits<Iterator>::reference;

        filter_iterator() = default;

        filter_iterator(Predicate pred, Iterator itor, Iterator end)
            : predicate_(pred), itor_(itor), end_(end)
        {
            if (itor != end)
                find_next();
        }

        reference operator*() { return *itor_; }
        value_type *operator->() const { return &(*itor_); }

        filter_iterator& operator++() {
            if ( itor_ != end_ )
                advance();
            find_next();
            return *this;
        }

        filter_iterator operator++(int) {
            auto it = *this;
            ++(*this);
            return it;
        }

        bool operator==(const filter_iterator& other) { return (itor_ == other.itor_); }
        bool operator!=(const filter_iterator& other) { return (itor_ != other.itor_); }
    };

    template<typename mesh_type>
    class is_boundary_pred
    {
        const mesh_type&    msh_;
    public:
        is_boundary_pred(const mesh_type& msh) : msh_(msh) {}

        template<typename T>
        bool operator()(const T& elem) { return msh_.is_boundary(elem); }
    };

    template<typename mesh_type>
    class is_internal_pred
    {
        const mesh_type&    msh_;
    public:
        is_internal_pred(const mesh_type& msh) : msh_(msh) {}

        template<typename T>
        bool operator()(const T& elem) { return !msh_.is_boundary(elem); }
    };

} //namespace priv

template<typename T>
[[deprecated]]
std::pair<bool, typename T::id_type>
find_element_id(const std::vector<T>& elements, const T& element)
{
    auto itor = std::lower_bound(elements.begin(), elements.end(), element);

    if (itor != elements.end() && !(element < *itor))
    {
        typename T::id_type id(std::distance(elements.begin(), itor));
        return std::make_pair(true, id);
    }

    return std::make_pair(false, typename T::id_type());
}

template<typename T, typename Iterator>
std::pair<bool, typename T::id_type>
find_element_id(const Iterator& begin, const Iterator& end, const T& element)
{
    auto itor = std::lower_bound(begin, end, element);

    if (itor != end && !(element < *itor))
    {
        typename T::id_type id(std::distance(begin, itor));
        return std::make_pair(true, id);
    }

    return std::make_pair(false, typename T::id_type());
}



/****************************************************************************/
namespace priv {

/* Mesh bones.
 *
 * @DIM     Space dimension
 * @T       Type that the mesh uses to represent points.
 */
template<typename T, size_t DIM, typename Storage>
class mesh_bones
{
    static_assert(DIM > 0 && DIM <= 3, "mesh: Allowed dimensions are 1, 2 and 3");

    std::shared_ptr<Storage>    m_storage;

public:
    mesh_bones() { m_storage = std::make_shared<Storage>(); }

    /* Return a shared_ptr to the backend storage. */
    std::shared_ptr<Storage>
    backend_storage(void)
    {
        return m_storage;
    }

    /* Return a shared_ptr to the backend storage. */
    const std::shared_ptr<Storage>
    backend_storage(void) const
    {
        return m_storage;
    }
};

/* Generic template for a mesh.
 *
 * This template has to be specialized for the 1D, 2D and 3D cases and it
 * represents the actual interface between the user and the mesh. It is in
 * `priv` and `mesh` inherits from `mesh_base` to provide an additional
 * decoupling layer.
 * The user should interact with the mesh in terms of cells and faces only.
 *
 * @DIM     Space dimension
 * @T       Type that the mesh uses to represent points.
 */
template<typename T, size_t DIM, typename Storage>
class mesh_base
{
    static_assert(DIM > 0 && DIM <= 3, "mesh: Allowed dimensions are 1, 2 and 3");
};



/* Template specialization for 3D meshes.
 *
 * @T       Type that the mesh uses to represent points.
 */
template<typename T, typename Storage>
class mesh_base<T,3,Storage> : public mesh_bones<T,3,Storage>
{
public:
    typedef typename Storage::volume_type                       volume_type;
    typedef typename Storage::surface_type                      surface_type;
    typedef typename Storage::edge_type                         edge_type;
    typedef typename Storage::node_type                         node_type;
    typedef T                                                   value_type;
    typedef typename Storage::point_type                        point_type;

    typedef volume_type                                         cell;
    typedef surface_type                                        face;
    const static size_t dimension = 3;

    /* cell iterators */
    typedef typename std::vector<volume_type>::iterator         cell_iterator;
    typedef typename std::vector<volume_type>::const_iterator   const_cell_iterator;

    cell_iterator           cells_begin() { return this->backend_storage()->volumes.begin(); }
    cell_iterator           cells_end()   { return this->backend_storage()->volumes.end(); }
    const_cell_iterator     cells_begin() const { return this->backend_storage()->volumes.begin(); }
    const_cell_iterator     cells_end()   const { return this->backend_storage()->volumes.end(); }

    /* face iterators */
    typedef typename std::vector<surface_type>::iterator        face_iterator;
    typedef typename std::vector<surface_type>::const_iterator  const_face_iterator;

    face_iterator           faces_begin() { return this->backend_storage()->surfaces.begin(); }
    face_iterator           faces_end()   { return this->backend_storage()->surfaces.end(); }
    const_face_iterator     faces_begin() const { return this->backend_storage()->surfaces.begin(); }
    const_face_iterator     faces_end()   const { return this->backend_storage()->surfaces.end(); }

    size_t  cells_size() const { return this->backend_storage()->volumes.size(); }
    size_t  faces_size() const { return this->backend_storage()->surfaces.size(); }

    bool is_boundary(typename face::id_type id) const
    {
        return this->backend_storage()->boundary_surfaces.at(id);
    }

    bool is_boundary(const face& f) const
    {
        auto e = find_element_id(faces_begin(), faces_end(), f);
        if (e.first == false)
            throw std::invalid_argument("Cell not found");

        return this->backend_storage()->boundary_surfaces.at(e.second);
    }

    bool is_boundary(const face_iterator& itor) const
    {
        auto ofs = std::distance(faces_begin(), itor);
        return this->backend_storage()->boundary_surfaces.at(ofs);
    }

    size_t  boundary_faces_size() const
    {
        return std::count(this->backend_storage()->boundary_surfaces.begin(),
                          this->backend_storage()->boundary_surfaces.end(),
                          true);
    }

    size_t  internal_faces_size() const
    {
        return std::count(this->backend_storage()->boundary_surfaces.begin(),
                          this->backend_storage()->boundary_surfaces.end(),
                          false);
    }
};



/* Template specialization for 2D meshes.
 *
 * @T       Type that the mesh uses to represent points.
 */
template<typename T, typename Storage>
class mesh_base<T,2,Storage> : public mesh_bones<T,2,Storage>
{
public:
    typedef typename Storage::surface_type                      surface_type;
    typedef typename Storage::edge_type                         edge_type;
    typedef typename Storage::node_type                         node_type;
    typedef T                                                   value_type;
    typedef typename Storage::point_type                        point_type;

    typedef surface_type                                        cell;
    typedef edge_type                                           face;


    const static size_t dimension = 2;

    /* cell iterators */
    typedef typename std::vector<surface_type>::iterator        cell_iterator;
    typedef typename std::vector<surface_type>::const_iterator  const_cell_iterator;

    cell_iterator           cells_begin() { return this->backend_storage()->surfaces.begin(); }
    cell_iterator           cells_end()   { return this->backend_storage()->surfaces.end(); }
    const_cell_iterator     cells_begin() const { return this->backend_storage()->surfaces.begin(); }
    const_cell_iterator     cells_end()   const { return this->backend_storage()->surfaces.end(); }

    /* face iterators */
    typedef typename std::vector<edge_type>::iterator           face_iterator;
    typedef typename std::vector<edge_type>::const_iterator     const_face_iterator;

    face_iterator           faces_begin() { return this->backend_storage()->edges.begin(); }
    face_iterator           faces_end()   { return this->backend_storage()->edges.end(); }
    const_face_iterator     faces_begin() const { return this->backend_storage()->edges.begin(); }
    const_face_iterator     faces_end()   const { return this->backend_storage()->edges.end(); }

    size_t  cells_size() const { return this->backend_storage()->surfaces.size(); }
    size_t  faces_size() const { return this->backend_storage()->edges.size(); }

    bool
    is_special_cell(typename cell::id_type& id) const
    {
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //This is only for triangles adaptation this must to be leave out for quadrilaterals
        #if 0
        bool check = false;
        if( this->backend_storage()->surfaces.at(id).size() > 3 )
            check = true;
        if(check != this->backend_storage()->special_surfaces.at(id).first)
        throw std::logic_error(Check is_special_cell, since cell is not considered has triangle.
        This is only for testing the right performance oh this functon for triangles);
        #endif
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        return this->backend_storage()->special_surfaces.at(id).first;
    }
    bool
    is_special_cell(const cell& cl) const
    {
        typename cell::id_type id = cl.get_id();

        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //This is only for triangles adaptation this must to be leave out for quadrilaterals
        #if 0
        bool check = false;
        if(cl.subelement_size() > 3 ) // Only when working with triangles !!!
            check = true;
        if(check != this->backend_storage()->special_surfaces.at(id).first)
            throw std::logic_error(Check is_special_cell, since cell is not considered has triangle.
            This is only for testing the right performance oh this functon for triangles);
        #endif
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        return this->backend_storage()->special_surfaces.at(id).first;
    }
    std::vector<typename point_type::id_type>
    get_vertices_pos(typename cell::id_type& id) const
    {
        auto sp =  this->backend_storage()->special_surfaces.at(id);

        #if 0
        if(sp.second.size() > 3 )
            throw std::logic_error(vertices ids are bigger than 3. Check get_vertices
            _pos or take this out when starting with quadrilaterals);
        #endif

        return sp.second;
    }

    std::vector<typename point_type::id_type>
    get_vertices_pos(const cell& cl) const
    {
        auto id = cl.get_id();
        auto sp =  this->backend_storage()->special_surfaces.at(id);

        #if 0
        if(sp.second.size() > 3 )
            throw std::logic_error("vertices ids are bigger than 3. Check get_vertices_pos or take this out when starting with quadrilaterals");
        #endif

        return sp.second;
    }

    std::vector<point_type>
    get_vertices(typename cell::id_type& id , const std::vector<point_type>& pts) const
    {
        typedef std::pair<bool, std::vector<typename point_type::id_type>> pair;
        pair sp =  this->backend_storage()->special_surfaces.at(id);
        std::vector<typename point_type::id_type> vertices_ids = sp.second;
        std::vector<point_type> vertices(vertices_ids.size());

        size_t i = 0;
        for(auto& id: vertices_ids)
            vertices.at(i++) = pts.at(id);
        #if 0
        if(vertices.size() > 3 )
                throw std::logic_error("vertices are bigger than 3. Check get_vertices or take this out when starting with quadrilaterals");
        #endif

        return vertices;
    }
    std::vector<point_type>
    get_vertices(const cell& cl, const std::vector<point_type>& pts) const
    {
        typedef std::pair<bool, std::vector<typename point_type::id_type>> pair;

        auto id = cl.get_id();
        auto sp =  this->backend_storage()->special_surfaces.at(id);
        auto vertices_ids = sp.second;
        std::vector<point_type> vertices(vertices_ids.size());

        size_t i = 0;
        assert(vertices_ids.size()<= pts.size());

        for(auto id: vertices_ids)
            vertices.at(i++) = pts.at(id);
        #if 0
        if(vertices.size() > 3 )
            throw std::logic_error("vertices are bigger than 3. Check get_vertices or take this out when starting with quadrilaterals");
        #endif

        return vertices;
    }
    cell
    neighbor(const  cell & cl,
                const  face & fc)
    {
        auto f = find_element_id(faces_begin(), faces_end(), fc);
        if (f.first == false)
            throw std::invalid_argument("Face not found");
        auto c = find_element_id(cells_begin(), cells_end(), cl);
        if (c.first == false)
            throw std::invalid_argument("Cell not found");

        auto cells = this->backend_storage()->edges_owners.at(f.second);
        //second part just to confirm cell is also owner
        if(cells.first != c.second && cells.second == c.second)
            return *std::next(this->backend_storage()->surfaces.begin(), cells.first);
        else if(cells.second != c.second && cells.first == c.second)
        return *std::next(this->backend_storage()->surfaces.begin(), cells.second);
        else
            throw std::logic_error("Neighbors not found or cell is not owner");
    }
    cell
    neighbor(const  cell & cl,
                const  face & fc) const
    {
        auto f = find_element_id(faces_begin(), faces_end(), fc);
        if (f.first == false)
            throw std::invalid_argument("Face not found");
        auto c = find_element_id(cells_begin(), cells_end(), cl);
        if (c.first == false)
            throw std::invalid_argument("Cell not found");

        auto cells = this->backend_storage()->edges_owners.at(f.second);
        //second part just to confirm cell is also owner
        if(cells.first != c.second && cells.second == c.second)
            return *std::next(this->backend_storage()->surfaces.begin(), cells.first);
        else if(cells.second != c.second && cells.first == c.second)
        return *std::next(this->backend_storage()->surfaces.begin(), cells.second);
        else
            throw std::logic_error("Neighbors not found or cell is not owner");
    }
    auto
    neighbor_id(const  cell & cl,
                const  face & fc)
    {
        auto f = find_element_id(faces_begin(), faces_end(), fc);
        if (f.first == false)
            throw std::invalid_argument("Face not found");
        auto c = find_element_id(cells_begin(), cells_end(), cl);
        if (c.first == false)
            throw std::invalid_argument("Cell not found");

        auto cells = this->backend_storage()->edges_owners.at(f.second);
        //second part just to confirm cell is also owner
        if(cells.first != c.second && cells.second == c.second)
            return cells.first;
        else if(cells.second != c.second && cells.first == c.second)
            return cells.second;
        else
            throw std::logic_error("Neighbors not found or cell is not owner");
    }

    auto
    neighbor_id(const  cell & cl,
                const  face & fc) const
    {
        auto f = find_element_id(faces_begin(), faces_end(), fc);
        if (f.first == false)
            throw std::invalid_argument("Face not found");
        auto c = find_element_id(cells_begin(), cells_end(), cl);
        if (c.first == false)
            throw std::invalid_argument("Cell not found");

        std::cout << "this->backend_storage()->edges_owners"<<this->backend_storage()->edges_owners.size() << std::endl;
        auto cells = this->backend_storage()->edges_owners.at(f.second);
        //second part just to confirm cell is also owner
        if(cells.first != c.second && cells.second == c.second)
            return cells.first;
        else if(cells.second != c.second && cells.first == c.second)
            return cells.second;
        else
            throw std::logic_error("Neighbors not found or cell is not owner");
    }

    std::pair<int,int>
    owners(const face & fc )
    {
        auto f = find_element_id(faces_begin(), faces_end(), fc);
        if (f.first == false)
            throw std::invalid_argument("Face not found");

        return this->backend_storage()->edges_owners.at(f.second);
    }
    std::pair<int,int>
    owners(const face_iterator & itor )
    {
        auto f = find_element_id(faces_begin(), faces_end(), *itor);
        if (f.first == false)
            throw std::invalid_argument("Face not found");

        return this->backend_storage()->edges_owners.at(f.second);
    }
    std::pair<int,int>
    owners(const const_face_iterator & itor )
    {
        auto f = find_element_id(faces_begin(), faces_end(), *itor);
        if (f.first == false)
            throw std::invalid_argument("Face not found");

        return this->backend_storage()->edges_owners.at(f.second);
    }

    bool is_boundary(typename face::id_type id) const
    {
        return this->backend_storage()->boundary_edges.at(id);
    }

    bool is_boundary(const face& f) const
    {
        auto e = find_element_id(faces_begin(), faces_end(), f);
        if (e.first == false)
            throw std::invalid_argument("Cell not found");

        return this->backend_storage()->boundary_edges.at(e.second);
    }
    auto boundary_information(const face& f) const
    {
        auto e = find_element_id(faces_begin(), faces_end(), f);
        if (e.first == false)
            throw std::invalid_argument("Face not found");

        return this->backend_storage()->boundary_info.at(e.second);
    }
    bool is_boundary(const face_iterator& itor) const
    {
        auto ofs = std::distance(faces_begin(), itor);
        return this->backend_storage()->boundary_edges.at(ofs);
    }

    size_t  boundary_faces_size() const
    {
        return std::count(this->backend_storage()->boundary_edges.begin(),
                          this->backend_storage()->boundary_edges.end(),
                          true);
    }
    size_t  boundary_type_faces_size() const
    {
        return std::count(this->backend_storage()->boundary_info.begin(),
                          this->backend_storage()->boundary_info.end(),
                          true);
    }


    size_t  internal_faces_size() const
    {
        return std::count(this->backend_storage()->boundary_edges.begin(),
                          this->backend_storage()->boundary_edges.end(),
                          false);
    }
};

/* mesh base class defining the data arrays for the 1D case */
template<typename T, typename Storage>
class mesh_base<T,1,Storage> : public mesh_bones<T,1,Storage>
{
public:
    typedef typename Storage::edge_type                     edge_type;
    typedef typename Storage::node_type                     node_type;
    typedef T                                               value_type;
    typedef typename Storage::point_type                    point_type;

    typedef edge_type                                       cell;
    typedef node_type                                       face;
    const static size_t dimension = 1;

    /* cell iterators */
    typedef typename std::vector<edge_type>::iterator       cell_iterator;
    typedef typename std::vector<edge_type>::const_iterator const_cell_iterator;

    cell_iterator           cells_begin() { return this->backend_storage()->edges.begin(); }
    cell_iterator           cells_end()   { return this->backend_storage()->edges.end(); }
    const_cell_iterator     cells_begin() const { return this->backend_storage()->edges.begin(); }
    const_cell_iterator     cells_end()   const { return this->backend_storage()->edges.end(); }

    /* face iterators */
    typedef typename std::vector<node_type>::iterator       face_iterator;
    typedef typename std::vector<node_type>::const_iterator const_face_iterator;

    face_iterator           faces_begin() { return this->backend_storage()->nodes.begin(); }
    face_iterator           faces_end()   { return this->backend_storage()->nodes.end(); }
    const_face_iterator     faces_begin() const { return this->backend_storage()->nodes.begin(); }
    const_face_iterator     faces_end()   const { return this->backend_storage()->nodes.end(); }

    size_t  cells_size() const { return this->backend_storage()->edges.size(); }
    size_t  faces_size() const { return this->backend_storage()->nodes.size(); }

    bool is_boundary(typename face::id_type id) const
    {
        return this->backend_storage()->boundary_nodes.at(id);
    }

    bool is_boundary(const face& f) const
    {
        auto e = find_element_id(faces_begin(), faces_end(), f);
        if (e.first == false)
            throw std::invalid_argument("Cell not found");

        return this->backend_storage()->boundary_nodes.at(e.second);
    }

    bool is_boundary(const face_iterator& itor) const
    {
        auto ofs = std::distance(faces_begin(), itor);
        return this->backend_storage()->boundary_nodes.at(ofs);
    }

    size_t  boundary_faces_size() const
    {
        return std::count(this->backend_storage()->boundary_nodes.begin(),
                          this->backend_storage()->boundary_nodes.end(),
                          true);
    }

    size_t  internal_faces_size() const
    {
        return std::count(this->backend_storage()->boundary_nodes.begin(),
                          this->backend_storage()->boundary_nodes.end(),
                          false);
    }

    #if 0
    auto boundary_type(const face& f) const
    {
        auto e = find_element_id(faces_begin(), faces_end(), f);
        if (e.first == false)
            throw std::invalid_argument("Cell not found");
        auto  ret = int(0);
        return ret;
    }
    #endif
};

} // namespace priv

template<typename T, size_t DIM, typename Storage>
class mesh : public priv::mesh_base<T,DIM,Storage>
{
public:
    static const size_t dimension = DIM;
    typedef T scalar_type;

    typedef typename priv::mesh_base<T, DIM, Storage>::point_type   point_type;
    typedef typename priv::mesh_base<T, DIM, Storage>::cell         cell;
    typedef typename priv::mesh_base<T, DIM, Storage>::face         face;

    /* point iterators */
    typedef typename std::vector<point_type>::iterator              point_iterator;
    typedef typename std::vector<point_type>::const_iterator        const_point_iterator;

    point_iterator          points_begin() { return this->backend_storage()->points.begin(); }
    point_iterator          points_end()   { return this->backend_storage()->points.end(); }
    const_point_iterator    points_begin() const { return this->backend_storage()->points.begin(); }
    const_point_iterator    points_end()   const { return this->backend_storage()->points.end(); }

    size_t  points_size() const { return this->backend_storage()->points.size(); }

    typedef priv::filter_iterator<typename mesh::face_iterator,
                                  priv::is_boundary_pred<mesh>>
                                  boundary_face_iterator;

    typedef priv::filter_iterator<typename mesh::face_iterator,
                                  priv::is_internal_pred<mesh>>
                                  internal_face_iterator;

    boundary_face_iterator  boundary_faces_begin()
    {
        typedef priv::is_boundary_pred<mesh> ibp;
        return boundary_face_iterator(ibp(*this), this->faces_begin(), this->faces_end());
    }

    boundary_face_iterator  boundary_faces_end()
    {
        typedef priv::is_boundary_pred<mesh> ibp;
        return boundary_face_iterator(ibp(*this), this->faces_end(), this->faces_end());
    }

    internal_face_iterator  internal_faces_begin()
    {
        typedef priv::is_internal_pred<mesh> iip;
        return internal_face_iterator(iip(*this), this->faces_begin(), this->faces_end());
    }

    internal_face_iterator  internal_faces_end()
    {
        typedef priv::is_internal_pred<mesh> iip;
        return internal_face_iterator(iip(*this), this->faces_end(), this->faces_end());
    }

    /* Apply a transformation to the mesh. Transform should be a functor or
     * a lambda function of type
     *      mesh_type::point_type -> mesh_type::point_type
     */
    template<typename Transform>
    void transform(const Transform& tr)
    {
        std::transform(points_begin(), points_end(), points_begin(), tr);
    }

    /* Returns the numerial ID of a cell. */
    typename cell::id_type lookup(const cell& cl) const
    {
        auto ci = find_element_id(this->cells_begin(), this->cells_end(), cl);
        if (!ci.first)
            throw std::invalid_argument("Cell not present in mesh");

        return ci.second;
    }

    /* Returns the numerial ID of a face. */
    typename face::id_type lookup(const face& fc) const
    {
        auto fi = find_element_id(this->faces_begin(), this->faces_end(), fc);
        if (!fi.first)
            throw std::invalid_argument("Face not present in mesh");

        return fi.second;
    }

    /* Th->maximumNumberOfFaces() */
    [[deprecated("subelement_size works only on generic_element")]]
    size_t  max_faces_per_element(void) const
    {
        size_t mfpe = 0;
        for (auto itor = this->cells_begin(); itor != this->cells_end(); itor++)
        {
            auto cell = *itor;
            mfpe = std::max(mfpe, cell.subelement_size());
        }

        return mfpe;
    }
};

template<typename T, size_t DIM, typename Storage>
typename mesh<T, DIM, Storage>::cell_iterator
begin(mesh<T, DIM, Storage>& msh)
{
    return msh.cells_begin();
}

template<typename T, size_t DIM, typename Storage>
typename mesh<T, DIM, Storage>::cell_iterator
end(mesh<T, DIM, Storage>& msh)
{
    return msh.cells_end();
}

template<typename T, size_t DIM, typename Storage>
typename mesh<T, DIM, Storage>::const_cell_iterator
begin(const mesh<T, DIM, Storage>& msh)
{
    return msh.cells_begin();
}

template<typename T, size_t DIM, typename Storage>
typename mesh<T, DIM, Storage>::const_cell_iterator
end(const mesh<T, DIM, Storage>& msh)
{
    return msh.cells_end();
}



template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
dump_to_matlab(const Mesh<T, 2, Storage>& msh, const std::string& filename)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file_dump_to_matlab"<<std::endl;

    ofs << " figure;"<<std::endl;
    ofs << " hold on"<<std::endl;

    for (auto cl : msh)
    {
        auto pts = points(msh, cl);
        auto b  = barycenter(msh,cl);
        auto id = msh.lookup(cl);
        auto fcs = faces(msh, cl);

        ofs << "strName = strtrim(cellstr(num2str("<< id <<",'(%d)')));"<<std::endl;
        ofs << "text("<<b.x()<< ","<< b.y() <<",strName,'VerticalAlignment','bottom');"<<std::endl;

        for (auto fc : fcs)
        {
            auto pts = points(msh, fc);
            if(pts.size() < 1)
               ofs << " \%display(\' no tiene pts\') ;"<<std::endl;


            auto fid = msh.lookup(fc);
            auto fb  = barycenter(msh, fc);
            ofs << "strName = strtrim(cellstr(num2str("<< fid <<",'(%d)')));"<<std::endl;
            ofs << "text("<<fb.x()<< ","<< fb.y() <<",strName,'VerticalAlignment','bottom', 'Color', 'r');"<<std::endl;

            if ( msh.is_boundary(fc) )
            {
                ofs << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
                ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'r');";
                ofs << std::endl;
            }
            else
            {
                ofs << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
                ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'k');";
                ofs << std::endl;
            }
        }
    }
    ofs.close();
}
template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
dump_to_matlab(const Mesh<T, 2, Storage>& msh, const std::string& filename,
                const std::vector<size_t> &vec, const std::string& some) // solve this, since this is declare in the same way the function below, thus I had to put a new variable some, taht do nothing
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file dump_to_matlab"<<std::endl;

    ofs << " figure;"<<std::endl;
    ofs << " hold on"<<std::endl;

    for (auto cl : msh)
    {
        auto pts = points(msh, cl);
        auto b  = barycenter(msh,cl);
        auto id = msh.lookup(cl);
        auto fcs = faces(msh, cl);

        ofs << "strName = strtrim(cellstr(num2str("<< id <<",'(%d)')));"<<std::endl;
        ofs << "text("<<b.x()<< ","<< b.y() <<",strName,'VerticalAlignment','bottom');"<<std::endl;

        for (auto fc : fcs)
        {
            auto pts = points(msh, fc);
            if(pts.size() < 1)
               ofs << " \%display(\' no tiene pts\') ;"<<std::endl;


            auto fid = msh.lookup(fc);
            auto fb  = barycenter(msh, fc);
            ofs << "strName = strtrim(cellstr(num2str("<< fid <<",'(%d)')));"<<std::endl;
            ofs << "text("<<fb.x()<< ","<< fb.y() <<",strName,'VerticalAlignment','bottom');"<<std::endl;

            if ( msh.is_boundary(fc) )
            {
                ofs << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
                ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'r');";
                ofs << std::endl;
            }
            else if(vec.at(fid) == 1)
            {
                ofs << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
                ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'm');";
                ofs << std::endl;
            }
            else
            {
                ofs << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
                ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'k');";
                ofs << std::endl;
            }

        }
    }
    ofs.close();
}

template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
dump_to_matlab(const Mesh<T, 2, Storage>& msh, const std::string& filename, const std::vector<size_t>& vec)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file dump_to_matlab"<<std::endl;

    ofs << " figure;"<<std::endl;
    ofs << " hold on"<<std::endl;

    for (auto cl : msh)
    {

        auto b  = barycenter(msh,cl);
        auto id = msh.lookup(cl);
        ofs<< "plot( "<< b.x() << ", " << b.y() <<",'r');"<<std::endl;
        ofs<< "str = strtrim(cellstr(num2str("<< vec.at(id) <<",'(%d)')));"<<std::endl;
        ofs<< "text("<<b.x()<< ","<< b.y() <<",str,'VerticalAlignment','bottom');"<<std::endl;

        auto fcs = faces(msh, cl);
        for (auto fc : fcs)
        {
            auto pts = points(msh, fc);
            auto fid = msh.lookup(fc);
            auto fb  = barycenter(msh, fc);
            ofs << "strName = strtrim(cellstr(num2str("<< fid <<",'(%d)')));"<<std::endl;
            ofs << "text("<<fb.x()<< ","<< fb.y() <<",strName,'VerticalAlignment','bottom');"<<std::endl;

            if ( msh.is_boundary(fc) )
            {
                ofs << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
                ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'b');";
                ofs << std::endl;
            }
            else
            {
                auto cells_ids = face_owner_cells_ids(msh,fc);
                auto cl1_mark = vec.at(cells_ids.first);
                auto cl2_mark = vec.at(cells_ids.second);
                auto mark = std::max(cl1_mark, cl2_mark);

                ofs << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
                if(mark == 3)
                    ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'r');";
                if(mark == 2)
                    ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'g');";
                if(mark == 1)
                    ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'c');";
                if(mark == 0)
                    ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'k');";
                ofs << std::endl;

            }
        }
    }
    ofs.close();
}
template<template<typename, size_t, typename> class Mesh,
         typename T, typename Storage>
void
dump_to_matlab(const Mesh<T, 2, Storage>& msh, const std::string& filename,
                const std::vector<size_t>& levels, const size_t imsh)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file dump_to_matlab"<<std::endl;

    ofs << " figure;"<<std::endl;
    ofs << " hold on"<<std::endl;
    for (auto cl : msh)
    {
        ofs << " hold on"<<std::endl;

        auto pts = points(msh, cl);
        ofs << " \% display(\'  celda "<< msh.lookup(cl)<<" \') ;"<<std::endl;

        if(pts.size() < 1)
           ofs << " \%display(\' no tiene pts\') ;"<<std::endl;

        for(auto& p : pts)
            ofs<< " \%plot( "<< p.x() << ", " << p.y() <<",'o');"<<std::endl;

        auto b  = barycenter(msh,cl);
        auto id = msh.lookup(cl);
        ofs<< "\%plot( "<< b.x() << ", " << b.y() <<",'r');"<<std::endl;
        //ofs<< "\%strLevel = strtrim(cellstr(num2str("<< levels.at(id) <<",'(%d)')));"<<std::endl;
        //ofs<< "\%text("<<b.x()<< ","<< b.y() <<",strLevel,'VerticalAlignment','bottom');"<<std::endl;

        ofs<< "\%strName = strtrim(cellstr(num2str("<< id <<",'(%d)')));"<<std::endl;
        ofs<< "\%text("<<b.x()<< ","<< b.y() <<",strName,'VerticalAlignment','bottom');"<<std::endl;


        auto fcs = faces(msh, cl);
        if(fcs.size() < 1)
            ofs << " \% display(\'cell "<< msh.lookup(cl)<<" no tiene caras\'); "<<std::endl;

        for (auto fc : fcs)
        {
            auto pts = points(msh, fc);
            if(pts.size() < 1)
               ofs << " \%display(\' face "<< msh.lookup(fc)<<" en la celda "<< msh.lookup(cl)<<" no tiene pts\'); "<<std::endl;

            if ( msh.is_boundary(fc) )
            {
                ofs << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
                ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'r');";
                ofs << std::endl;
            }
            else
            {
                ofs << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
                ofs << pts[0].y() << " " << pts[1].y() << "], 'Color', 'k');";
                ofs << std::endl;
            }
        }
    }
    ofs.close();
}


} // namespace disk
