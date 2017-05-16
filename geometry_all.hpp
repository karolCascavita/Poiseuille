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
        auto cell_diameter = diameter(msh , cell);

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
        auto cell_diameter = diameter(msh , cell);

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
        h += diameter(msh, cl);
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

/* Compute the barycenter of a cell */
template<typename Mesh, typename Element>
point<typename Mesh::value_type, Mesh::dimension>
barycenter(const Mesh& msh, const Element& elm)
{
    auto pts = points(msh, elm);
    auto bar = std::accumulate(std::next(pts.begin()), pts.end(), pts.front());
    return bar / typename Mesh::value_type( pts.size() );
}


} // namespace disk
