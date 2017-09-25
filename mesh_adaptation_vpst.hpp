/*
 *       /\
 *      /__\
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
#include <iomanip>

#include "common/eigen.hpp"
#include "bases/bases_utils.hpp"
#include "bases/bases_ranges.hpp"

namespace disk {

template<typename T>
struct mesh_parameters
{
    mesh_parameters()
    {
        start_exact = 0;
        mesh_name   = 0;
        num_remesh  = 0;
        marker_name = 2;
        initial_imsh = 0;
        percent     = 0.1;
        recycle     = true;
        diff        = false;
        mark_all    = false; // There is a problem here read comment
        hanging_nodes = true;
        call_mesher   = false;
    }
    // If marl_all ==true, it cannot be used  in circular for iter == 0, that used quadrilaterals
    // Since, we have triangles in the center, that will make hanging_nodes to appear
    // However, since the quadrilaterals are split in triangles and doesn't support
    // hanging nodes, the triangles made by quadrilaterals won-t see the hanging_nodes

    T       percent;
    int     mesh_name;
    int     num_remesh;
    bool    call_mesher;
    int     initial_imsh;
    int     marker_name;
    int     start_exact;
    bool    recycle;
    bool    hanging_nodes;
    bool    diff;
    bool    mark_all;
    std::string     short_mesh_name;
    std::string     directory;
    std::string     summary;
    std::string     summary_old;

    friend std::ostream& operator<<(std::ostream& os, const mesh_parameters<T>& mp) {
        os << "Mesh Parameters: "<<std::endl;
        os << "* percent      : "<< mp.percent<< std::endl;
        os << "* mesh_name    : "<< mp.mesh_name<< std::endl;
        os << "* num_remesh   : "<< mp.num_remesh<< std::endl;
        os << "* call_mesher  : "<< mp.call_mesher<< std::endl;
        os << "* initial_imsh : "<< mp.initial_imsh<< std::endl;
        os << "* marker_name  : "<< mp.marker_name<< std::endl;
        os << "* recycle      : "<< mp.recycle<< std::endl;
        os << "* hanging_nodes: "<< mp.hanging_nodes<< std::endl;
        os << "* diffusion    : "<< mp.diff<< std::endl;
        os << "* mark_all     : "<< mp.mark_all<< std::endl;
        os << "* short_mesh_name: "<< mp.short_mesh_name<< std::endl;
        os << "* directory    : "<< mp.directory<< std::endl;
        os << "* summary      : "<< mp.summary<< std::endl;
        os << "* summary_old  : "<< mp.summary_old<< std::endl;
        return os;
    }
};

template<typename MeshType, typename CellType, typename FaceType>
size_t
face_position(const MeshType& msh, const CellType& cell, const FaceType& face)
{
    auto c_faces = faces(msh, cell);
    size_t   j = 0;
    for(auto& fc : c_faces)
    {
        if(fc == face)
            return j;
        j++;
    }
    throw std::invalid_argument("This is a bug: face not found");
};

template<typename T, size_t DIM, typename Storage>
std::vector<dynamic_vector<T>>
solution_zero_vector(const disk::mesh<T,DIM,Storage>& msh, const size_t degree)
{
    typedef disk::mesh<T,DIM,Storage>              mesh_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, typename mesh_type::cell>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, typename mesh_type::face>    face_basis_type;

    std::vector<dynamic_vector<T>> U(msh.cells_size());
    size_t i = 0;
    for (auto& cl : msh)
    {
        auto fcs = faces(msh, cl);
        auto num_faces       = fcs.size();
        auto num_cell_dofs   = cell_basis_type(degree).size();
        auto num_face_dofs   = face_basis_type(degree).size();
        disk::dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);
        U[i] =  dynamic_vector<T>::Zero(dsr.total_size());
        ++i;
    }
    return U;
};
}//Disk
