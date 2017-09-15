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
#include <algorithm>    // std::find

#include "geometry/geometry_all.hpp"
#include "geometry/element_generic.hpp"

namespace disk {

namespace generic_priv {

template<size_t DIM>
struct element_types
{
    static_assert(DIM > 0 && DIM <= 3, "element_types: CODIM must be less than DIM");
};

template<>
struct element_types<3> {
        typedef generic_element<3,0>    volume_type;
        typedef generic_element<3,1>    surface_type;
        typedef generic_element<3,2>    edge_type;
        typedef generic_element<3,3>    node_type;
};

template<>
struct element_types<2> {
        typedef generic_element<2,0>    surface_type;
        typedef generic_element<2,1>    edge_type;
        typedef generic_element<2,2>    node_type;
};

template<>
struct element_types<1> {
        typedef generic_element<1,0>    edge_type;
        typedef generic_element<1,1>    node_type;
};

} // namespace priv

template<typename T, size_t DIM>
using generic_mesh_storage = mesh_storage<T, DIM, generic_priv::element_types<DIM>>;

template<typename T, size_t DIM>
using generic_mesh = mesh<T, DIM, generic_mesh_storage<T, DIM>>;

/* Return the number of elements of the specified cell */
template<typename T, size_t DIM>
size_t
number_of_faces(const generic_mesh<T,DIM>& msh, const typename generic_mesh<T,DIM>::cell& cl)
{
    return cl.subelement_size();
}

/* Return the actual faces of the specified element */
template<typename T, size_t DIM>
std::vector<typename generic_mesh<T, DIM>::face>
faces(const generic_mesh<T, DIM>& msh, const typename generic_mesh<T, DIM>::cell& cl)
{
    auto id_to_face = [&](const typename generic_mesh<T, DIM>::face::id_type& id) -> auto {
        return *std::next(msh.faces_begin(), size_t(id));
    };

    std::vector<typename generic_mesh<T, DIM>::face> ret;
    ret.resize( cl.subelement_size() );

    std::transform(cl.subelement_id_begin(), cl.subelement_id_end(),
                   ret.begin(), id_to_face);

    return ret;
}

template<typename T, size_t DIM>
std::pair<int,int>
face_owner_cells_ids(const generic_mesh<T, DIM>& msh,
                    const typename generic_mesh<T, DIM>::face& fc)
{
    std::pair<int,int> ret;

    std::vector<typename generic_mesh<T, DIM>::cell::id_type> vec;
    for(auto& cl:msh)
    {
        auto fcs  = faces(msh,cl);
        std::sort(fcs.begin(), fcs.end());
        auto itor = std::lower_bound(fcs.begin(), fcs.end(), fc);

        if (itor != fcs.end() && !(fc < *itor))
        {
            auto cl_id = cl.get_id();
            vec.push_back(cl_id);
        }
    }
    if(vec.size() == 0)
        throw std::logic_error(" No owners found for this face");
    else if(vec.size() == 1)
    {
        if(msh.is_boundary(fc))
        {
            ret.first = vec[0];
            ret.second = -1;
            return ret;
        }
        else
            throw std::logic_error("Just one owner found for a inner face.");
    }
    else if(vec.size() == 2)
    {
        if(msh.is_boundary(fc))
            throw std::logic_error("Two owners found for a Boundary edge.");
        else
        {
            ret.first  = vec[0];
            ret.second = vec[1];;
            return ret;
        }
    }
    else
        throw std::logic_error(" Finding more than 2 cell for the face");

}


template<typename T, size_t DIM>
void
check_older_msh(const generic_mesh<T, DIM>& msh)
{
    auto storage = msh.backend_storage();
    auto points  = storage ->points;
    std::cout << "/ **************** POINTS        **************** /" << std::endl;
    for(auto& p: points)
        std::cout << "      "<< p << std::endl;
    std::cout << "/ **************** FACES:  cells **************** /" << std::endl;
    //typename generic_mesh<T, DIM>::cell::id_type
    for(auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
    {
        auto fc = *itor;
        std::cout << "fc_id  : "<< msh.lookup(fc)<< "      ; owner_cells:";

        for(auto& cl:msh)
        {
            auto fcs  = faces(msh,cl);
            std::sort(fcs.begin(), fcs.end());
            auto itor = std::lower_bound(fcs.begin(), fcs.end(), fc);

            if (itor != fcs.end() && !(fc < *itor))
            {
                std::cout << "     "<< msh.lookup(cl);
            }
        }
        std::cout<< std::endl;
    }
    std::cout << "/ **************** FACES:  nodes **************** /" << std::endl;

    //typename generic_mesh<T, DIM>::cell::id_type
    for(auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
    {
        auto fc = *itor;
        std::cout << "fc_id  : "<< msh.lookup(fc)<< "      ; nodes :";
        auto nodes = fc.point_ids();
        for(auto& n: nodes)
            std::cout << "     "<< n;
        std::cout<< std::endl;
    }
    std::cout << "/ **************** CELSS:  faces **************** /" << std::endl;

    for(auto& cl : msh)
    {
        std::cout << "cell : "<< msh.lookup(cl)<< "      ; faces :";
        auto fcs  = faces(msh,cl);

        for(auto& fc : fcs)
            std::cout <<"     "<< msh.lookup(fc);
        std::cout<< std::endl;
    }
    std::cout<< std::endl;

    std::cout << "/ **************** CELSS:  nodes **************** /" << std::endl;

    for(auto& cl : msh)
    {
        auto id = msh.lookup(cl);
        std::cout << "cell : "<< id << "      ; nodes :  ("<< storage->special_surfaces.at(id).first;
        std::cout << ")" ;
        auto nodes  = cl.point_ids();

        for(auto& p : nodes)
            std::cout <<"     "<< p;
        std::cout<< std::endl;
    }
    std::cout<< std::endl;


}



/* Compute the measure of a 2-cell (= area) */
template<typename T>
T
measure(const generic_mesh<T,2>& msh, const typename generic_mesh<T,2>::cell& cl)
{
    auto pts = points(msh, cl);



    T acc{};
    for (size_t i = 1; i < pts.size() - 1; i++)
    {
        auto u = (pts.at(i) - pts.at(0)).to_vector();
        auto v = (pts.at(i+1) - pts.at(0)).to_vector();
        auto n = cross(u, v);
        acc += n.norm() / T(2);
    }

    auto vts_ids = msh.get_vertices_pos(cl);
    std::vector<point<T,2>> vts(vts_ids.size());
    size_t i = 0;
    for(auto id: vts_ids)
        vts.at(i++) = pts.at(id);

    T acc_v{};
    for (size_t i = 1; i < vts.size() - 1; i++)
    {
        auto u = (vts.at(i) - vts.at(0)).to_vector();
        auto v = (vts.at(i+1) - vts.at(0)).to_vector();
        auto n = cross(u, v);
        acc_v += n.norm() / T(2);
    }

    if(std::abs(acc_v -acc)/acc > 1.e-4)
    {
        std::cout << "cell :"<< cl.get_id() << std::endl;
        std::cout << "diff :"<< std::abs(acc_v -acc)<< std::endl;
        std::cout << "area :"<< acc << std::endl;
        std::cout << "area without hanging_nodes :"<< acc_v << std::endl;
        throw std::logic_error("Areas are not the same, review get_vertices!!");
    }
    return acc;
}

/* Compute the measure of a 2-face (= length) */
template<typename T>
T
measure(const generic_mesh<T,2>& msh, const typename generic_mesh<T,2>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);
    return (pts[1] - pts[0]).to_vector().norm();
}

/* Compute the measure of a 1-cell (= area) */
template<typename T>
T
measure(const generic_mesh<T,1>& msh, const typename generic_mesh<T,1>::cell& cl)
{
    auto pts = points(msh, cl);
    assert(pts.size() == 2);
    return (pts[1] - pts[0]).to_vector().norm();
}

/* Compute the measure of a 1-face (= length) */
template<typename T>
T
measure(const generic_mesh<T,1>& msh, const typename generic_mesh<T,1>::face& fc)
{
    return T(1);
}

template<typename Mesh, typename Element>
typename Mesh::scalar_type
diameter(const Mesh& msh, const Element& elem)
{
    auto pts = points(msh, elem);

    typename Mesh::scalar_type diam = 0.;

    for (size_t i = 0; i < pts.size(); i++)
        for (size_t j = i+1; j < pts.size(); j++)
            diam = std::max((pts[i] - pts[j]).to_vector().norm(), diam);

    return diam;
}
template<typename Mesh>
typename Mesh::scalar_type
diameter(const Mesh& msh, const std::vector<typename Mesh::point_type>& pts)
{
    typename Mesh::scalar_type diam = 0.;

    for (size_t i = 0; i < pts.size(); i++)
        for (size_t j = i+1; j < pts.size(); j++)
            diam = std::max((pts[i] - pts[j]).to_vector().norm(), diam);

    return diam;
}

/* Compute the barycenter of a 2-face */
template<typename T>
point<T,2>
barycenter(const generic_mesh<T,2>& msh, const typename generic_mesh<T,2>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);
    auto bar = (pts[0] + pts[1]) / T(2);
    return bar;
}

template<typename T>
static_vector<T, 2>
normal(const generic_mesh<T,2>& msh,
       const typename generic_mesh<T,2>::cell& cl,
       const typename generic_mesh<T,2>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);

    auto v = pts[1] - pts[0];
    auto n = (point<T,2>({-v.y(), v.x()})).to_vector();

    auto cell_bar = barycenter(msh, cl);
    auto face_bar = barycenter(msh, fc);
    auto outward_vector = (face_bar - cell_bar).to_vector();

    if ( n.dot(outward_vector) < T(0) )
        return -n/n.norm();

    return n/n.norm();
}

template<typename T>
T
normal(const generic_mesh<T,1>& msh,
       const typename generic_mesh<T,1>::cell& cl,
       const typename generic_mesh<T,1>::face& fc)
{
    auto fcs = faces(msh, cl);
    assert(fcs.size() == 2);

    if (fc == fcs[0])
        return -1.;

    if (fc == fcs[1])
        return 1.;

    throw std::logic_error("shouldn't have arrived here");
}

template<typename T>
T
normal_factor(const generic_mesh<T,2>& msh,
       const typename generic_mesh<T,2>::cell& cl,
       const typename generic_mesh<T,2>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);


    auto v = pts[1] - pts[0];
    auto n = (point<T,2>({-v.y(), v.x()})).to_vector();
    auto cell_bar = barycenter(msh, cl);
    auto face_bar = barycenter(msh, fc);
    auto outward_vector = (face_bar - cell_bar).to_vector();

    //std::cout << "value dot product ("<< cl.get_id()<<"): "<<  n.dot(outward_vector) << std::endl;
    //std::cout << " normal : "<< n/n.norm() <<std::endl;
    if ( n.dot(outward_vector) < T(0) )
        return -1.;

    return 1.;
}

template<typename T>
static_vector<T, 2>
normal_face(const generic_mesh<T,2>& msh,
       const typename generic_mesh<T,2>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);

    auto v = pts[1] - pts[0];
    auto n = (point<T,2>({-v.y(), v.x()})).to_vector();

    return n/n.norm();
}



} // namespace disk
