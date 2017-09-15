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

#include <algorithm>
#include <vector>

#include "mesh/mesh.hpp"
#include "common/util.h"

namespace disk {

/* Compute an estimate of the mesh discretization step 'h' */
template<typename Mesh>
typename Mesh::scalar_type
mesh_h_max(const Mesh& msh)
{
    typename Mesh::scalar_type h{};
    for (auto itor = msh.cells_begin(); itor != msh.cells_end(); itor++)
    {
        auto cell = *itor;
        auto pts = points(msh,cell);
        auto vts = msh.get_vertices(cell, pts);
        auto cell_diameter   = diameter(msh, vts);
        //auto cell_diameter = diameter(msh , cell);

        h = std::max(h, cell_diameter);
    }
    return h;
}
template<typename Mesh>
typename Mesh::scalar_type
mesh_h_min(const Mesh& msh)
{
    typename Mesh::scalar_type h{};
    h = 1.;
    for (auto itor = msh.cells_begin(); itor != msh.cells_end(); itor++)
    {
        auto cell = *itor;
        //auto cell_diameter = diameter(msh , cell);
        auto pts = points(msh,cell);
        auto vts = msh.get_vertices(cell, pts);
        auto cell_diameter   = diameter(msh, vts);
        h = std::min(h, cell_diameter);
    }

    return h;
}

template<typename Mesh>
typename Mesh::scalar_type
average_diameter(const Mesh& msh)
{
    typename Mesh::scalar_type h{};
    for (auto& cl : msh)
    {
        auto pts = points(msh,cl);
        auto vts = msh.get_vertices(cl, pts);
        auto cell_diameter   = diameter(msh, vts);
        h += cell_diameter;
        //h += diameter(msh, cl);

    }

    return h/msh.cells_size();
}
template<typename Mesh, typename Element>
std::vector<typename Mesh::point_type>
points(const Mesh& msh, const Element& elem)
{
    auto ptids = elem.point_ids();

    auto ptid_to_point = [&](const point_identifier<Mesh::dimension>& pi) -> auto {
        auto itor = msh.points_begin();
        std::advance(itor, pi);
        return *itor;
    };

    std::vector<typename Mesh::point_type> pts(ptids.size());
    std::transform(ptids.begin(), ptids.end(), pts.begin(), ptid_to_point);

    return pts;
}

template<typename T, typename Storage>
point<T,2>
//barycenter(const generic_mesh<T,2>& msh, const typename generic_mesh<T,2>::cell& cl)
barycenter(const mesh<T, 2, Storage>& msh, const typename mesh<T,2, Storage>::cell& cl)
{
   auto pts = points(msh, cl);
   auto ptsnum = pts.size();

  // typedef typename generic_mesh<T,2>::point_type point_type;
  typedef typename mesh<T,2,Storage>::point_type point_type;

   point_type  bar = point_type({0.0, 0.0});
   T           area = 0.0;


   for (size_t i = 0; i < ptsnum; i++)
   {
       auto p0 = pts[i];
       auto p1 = pts[(i+1)%ptsnum];

       auto a = p0.x()*p1.y() - p1.x()*p0.y();

       bar = bar + (p0 + p1) * a;
       area += a;
   }

   area *= 0.5;

   return bar/(6.0 * area);
}

} // namespace disk
