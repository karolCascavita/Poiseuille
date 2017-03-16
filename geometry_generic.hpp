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
        return *(msh.faces_begin() + size_t(id));
    };

    std::vector<typename generic_mesh<T, DIM>::face> ret;
    ret.resize( cl.subelement_size() );

    std::transform(cl.subelement_id_begin(), cl.subelement_id_end(),
                   ret.begin(), id_to_face);

    return ret;
}

template<typename T, size_t DIM>
std::vector<typename generic_mesh<T, DIM>::cell::id_type>
face_owner_cells_ids(const generic_mesh<T, DIM>& msh,
                        const typename generic_mesh<T, DIM>::face& fc)
{
    #if 0
    std::vector<typename generic_mesh<T, DIM>::cell::id_type> ret;
    /* Search the face in the other cells*/
    for(auto& cl:msh)
    {
        auto fcs = faces(msh,cl);
        std::sort(fcs.begin(), fcs.end());
        if(std::binary_search(fcs.begin(),fcs.end(),fc))
        {
            auto cl_id = msh.lookup(cl);
            //std::cout << "/* cell_id = "<< cl_id<<" */" << std::endl;
            ret.push_back(cl_id);
        }
    }
    if(ret.size() > 2)
        throw std::logic_error("face have more than 2 owners");
    // WK:This should be redefine as a pair, instead of a vector, so it will always have just 2 owners
    // and I dont need to check the size
    return ret;
    #endif

    std::vector<typename generic_mesh<T, DIM>::cell::id_type> ret;
    for(auto& cl:msh)
    {
        auto fcs  = faces(msh,cl);
        std::sort(fcs.begin(), fcs.end());
        auto itor = std::lower_bound(fcs.begin(), fcs.end(), fc);

        if (itor != fcs.end() && !(fc < *itor))
        {
            auto cl_id = cl.get_id();
            ret.push_back(cl_id);
        }
    }
    if(ret.size() > 2)
        throw std::logic_error(" Finding more than 2 cell for the face");
    else if(ret.size() == 0)
        throw std::logic_error(" No owners found for the face");
    else
        return ret;
}


template<typename T, size_t DIM>
typename generic_mesh<T, DIM>::cell::id_type
face_owner_cells_ids(const generic_mesh<T, DIM>& msh,
                        const typename generic_mesh<T, DIM>::face& face,
                        const typename generic_mesh<T, DIM>::cell& cell)
{
    #if 0
    typedef typename generic_mesh<T, DIM>::cell::id_type cell_id_type;
    cell_id_type ret;
    /* Search the face in the other cells*/
    auto cell_id = msh.lookup(cell);

    for(auto& cl:msh)
    {
        auto fcs = faces(msh,cl);
        std::sort(fcs.begin(), fcs.end());
        if(std::binary_search(fcs.begin(),fcs.end(),fc))
        {
            cell_id_type cl_id = msh.lookup(cl);

            std::cout << "  inside  cl_id "<< cl_id << std::endl;

            if(!(cl_id == cell_id))
            {
                ret = cl_id;
                std::cout << "  inside  ngh_id "<< cl_id << std::endl;
                break;
            }
        }
    }
    #endif
    if(msh.is_boundary(face))
    {
        throw std::logic_error("Do not call face_owner_cells_ids(msh, fc, cl) for boundaray edges; since it gives the neighbor of cl for this face.");
    }
    else
    {
        auto cell_id =cell.get_id();

        for(auto& cl : msh)
        {
            auto fcs   = faces(msh,cl);
            std::sort(fcs.begin(), fcs.end());
            auto itor  = std::lower_bound(fcs.begin(), fcs.end(), face);
            auto neighbor_id = cl.get_id();

            if( itor != fcs.end() && !(face < *itor) && (neighbor_id != cell_id))
                return neighbor_id;
        }
        throw std::logic_error(" Not neighbor found");
    }
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
        std::cout << "cell : "<< msh.lookup(cl)<< "      ; nodes :";
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

template<typename T>
T
diameter(const generic_mesh<T,2>& msh, const typename generic_mesh<T,2>::cell& cl)
{

    T c_meas = measure(msh, cl);
    T af_meas = T(0);
    auto fcs = faces(msh, cl);
    for (auto& f : fcs)
        af_meas += measure(msh, f);

    return c_meas/af_meas;

    /*
    auto bar = barycenter(msh, cl);
    auto pts = points(msh, cl);

    T diam = 0.;

    for (auto& pt : pts)
    {
        auto d = (pt-bar).to_vector().norm();
        diam = std::max(d, diam);
    }

    return diam;
    */
}

template<typename T>
T
diameter(const generic_mesh<T,2>& msh, const typename generic_mesh<T,2>::face& fc)
{
    return measure(msh, fc);
}

template<typename T>
T
diameter(const generic_mesh<T,1>& msh, const typename generic_mesh<T,1>::cell& cl)
{
    return measure(msh, cl);
}

template<typename T>
T
diameter(const generic_mesh<T,1>&, const typename generic_mesh<T,1>::face&)
{
    return 1.;
}
template<typename T,  size_t DIM>
T
diameter(const generic_mesh<T,DIM>& msh)
{
    T diam = T(0);
    for(auto& cl : msh)
    {
        auto cl_diam = diameter(msh, cl);
        diam = std::max(diam, cl_diam);
    }
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




} // namespace disk
