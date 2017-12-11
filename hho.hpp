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

#include "common/eigen.hpp"
#include "bases/bases_utils.hpp"
#include "bases/bases_ranges.hpp"
#include "timecounter.h"

namespace disk {

    template<typename T>
    dynamic_matrix<T>
    make_vectorial_matrix(const Eigen::Matrix<T, Eigen::Dynamic, 1>& phi, const int DIM)
    {
        dynamic_matrix<T> ret = dynamic_matrix<T>::Zero(DIM, DIM * phi.size());
        auto cbs = phi.size();

        for (size_t d = 0; d < DIM; d++)
            ret.block(d, cbs* d, 1, cbs) = phi.transpose();

        return ret;
    }


template<typename CellBasisType, typename CellQuadType, typename Mesh,
         typename Function>
dynamic_vector<typename Mesh::scalar_type>
compute_rhs(const Mesh& msh, const typename Mesh::cell& cl,
            const Function& f, size_t degree)
{
    typedef dynamic_vector<typename Mesh::scalar_type> vector_type;

    auto cell_basis     = CellBasisType(degree);
    auto cell_quad      = CellQuadType(2 * degree);
    vector_type ret = vector_type::Zero(cell_basis.size());

    auto cell_quadpoints = cell_quad.integrate(msh, cl);
    for (auto& qp : cell_quadpoints)
    {
        auto phi = cell_basis.eval_functions(msh, cl, qp.point());
        auto fval = f(qp.point());
        ret += qp.weight() * fval * phi;
    }
    return ret;
}


template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class projector_nopre
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef CellBasisType                       cell_basis_type;
    typedef CellQuadType                        cell_quadrature_type;
    typedef FaceBasisType                       face_basis_type;
    typedef FaceQuadType                        face_quadrature_type;

    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    cell_basis_type                             cell_basis;
    cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;

    size_t                                      m_degree_cell;
    size_t                                      m_degree;

public:

    projector_nopre()
    {
        auto degree = 1;

        m_degree_cell = degree;
        m_degree = degree;

        cell_basis          = cell_basis_type(m_degree_cell);
        cell_quadrature     = cell_quadrature_type(2*m_degree_cell + 2);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree + 2);
    }
    #if 0
    projector_nopre(size_t degree)
    {
        m_degree_cell = degree;
        m_degree= degree;

        cell_basis          = cell_basis_type(m_degree_cell);
        cell_quadrature     = cell_quadrature_type(2*m_degree_cell + 2);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree + 2);
    }
    #endif
    projector_nopre(size_t degree, const bool& change_degree)
    {
        m_degree =  degree;

        m_degree_cell = degree + change_degree;

        cell_basis          = cell_basis_type(m_degree_cell);
        cell_quadrature     = cell_quadrature_type(2*m_degree_cell + 2);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree + 2);
    }

    matrix_type cell_mm;

    template<typename Function>
    vector_type
    compute_vec_cell(const mesh_type& msh, const cell_type& cl, const Function& f)
    {
        //if(mp.diff)
        //auto number = set_cell_number(msh, cl);
        auto number = 0;
        auto DIM  = mesh_type::dimension;
        auto zero_range = cell_basis.range(0, m_degree);
        auto cbs  = zero_range.size() * DIM;
        matrix_type mm  = matrix_type::Zero(cbs, cbs);
        vector_type rhs = vector_type::Zero(cbs);

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto phi = cell_basis.eval_functions(msh, cl, qp.point(), 0, m_degree);

            matrix_type vec_phi = make_vectorial_matrix(phi, DIM);
            auto fval = f(qp.point(), number);

#ifdef USE_BLAS
    #if 0
            int sz = phi.size();
            double w = qp.weight();
            int one = 1;
            dger(&sz, &sz, &w, phi.data(), &one, phi.data(), &one, mm.data(), &sz);
            w *= fval;
            daxpy(&sz, &w, phi.data(), &one, rhs.data(), &one);
    #endif
#else

            mm  += qp.weight() * vec_phi.transpose() * vec_phi;
            rhs += qp.weight() * vec_phi.transpose() * fval;
#endif
        }
        return mm.ldlt().solve(rhs);
    }

    template<typename Function>
    vector_type
    compute_cell(const mesh_type& msh, const cell_type& cl, const Function& f)
    {
        //if(mp.diff)
        //auto number = set_cell_number(msh, cl);
        auto number = 0;
        matrix_type mm = matrix_type::Zero(cell_basis.size(), cell_basis.size());
        vector_type rhs = vector_type::Zero(cell_basis.size());

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto phi = cell_basis.eval_functions(msh, cl, qp.point());
            auto fval = f(qp.point(), number);

#ifdef USE_BLAS
            int sz = phi.size();
            double w = qp.weight();
            int one = 1;
            dger(&sz, &sz, &w, phi.data(), &one, phi.data(), &one, mm.data(), &sz);
            w *= fval;
            daxpy(&sz, &w, phi.data(), &one, rhs.data(), &one);
#else

            mm  += qp.weight() * phi * phi.transpose();
            rhs += qp.weight() * phi * fval;
#endif
        }

        cell_mm = mm;
        return mm.llt().solve(rhs);
    }

    template<typename Function>
    vector_type
    compute_whole(const mesh_type& msh, const cell_type& cl, const Function& f)
    {
        //if(mp.diff)
        //auto  number = set_cell_number(msh, cl);

        auto fcs = faces(msh, cl);
        vector_type ret(cell_basis.size() + fcs.size()*face_basis.size());

        ret.block(0, 0, cell_basis.size(), 1) = compute_cell(msh, cl, f);

        size_t face_offset = cell_basis.size();
        for (auto& fc : fcs)
        {
            matrix_type mm = matrix_type::Zero(face_basis.size(), face_basis.size());
            vector_type rhs = vector_type::Zero(face_basis.size());

            auto face_quadpoints = face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
                auto phi = face_basis.eval_functions(msh, fc, qp.point());
                auto fval = f(qp.point(), 0);

                mm  += qp.weight() * phi * phi.transpose();
                rhs += qp.weight() * phi * fval;
            }

            ret.block(face_offset, 0, face_basis.size(), 1) = mm.llt().solve(rhs);
            face_offset += face_basis.size();
        }

        return ret;
    }
};


template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class gradient_reconstruction_nopre
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef CellBasisType                       cell_basis_type;
    typedef CellQuadType                        cell_quadrature_type;
    typedef FaceBasisType                       face_basis_type;
    typedef FaceQuadType                        face_quadrature_type;

    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    cell_basis_type                             cell_basis;
    cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;

    size_t                                      m_degree, m_degree_cell;

public:
    matrix_type     oper;
    matrix_type     gradrecoper;

    matrix_type     data;

    gradient_reconstruction_nopre()
        : m_degree(1)
    {
        m_degree_cell = m_degree;
        //if()        m_degree_cell = m_degree + 1;

        cell_basis          = cell_basis_type(m_degree + 1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }
    #if 0
    gradient_reconstruction_nopre(const size_t degree)
        : m_degree(degree)
    {
        m_degree_cell = m_degree;
        //if()        m_degree_cell = m_degree + 1;

        cell_basis          = cell_basis_type(m_degree + 1);
        cell_quadrature     = cell_quadrature_type(2*(degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }
    #endif
    gradient_reconstruction_nopre(const size_t degree, const bool& change_degree)
        : m_degree(degree), m_degree_cell(degree)
    {
        if(change_degree)
            m_degree_cell = degree + 1;

        cell_basis          = cell_basis_type(m_degree + 1);
        cell_quadrature     = cell_quadrature_type(2*(degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }



    void compute(const mesh_type& msh, const cell_type& cl)
    {
        auto zero_mindeg = 0;
        auto zero_maxdeg = m_degree;

        auto cell_mindeg = 0;
        auto cell_maxdeg = m_degree_cell;

        auto one_mindeg = 1;
        auto one_maxdeg = m_degree + 1;

        auto zero_range = cell_basis.range( zero_mindeg,  zero_maxdeg);
        auto cell_range = cell_basis.range( cell_mindeg,  cell_maxdeg);
        auto one_range  = cell_basis.range( one_mindeg ,  one_maxdeg);

        matrix_type stiff_mat = matrix_type::Zero(cell_basis.size(), cell_basis.size());

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point());
            stiff_mat  += qp.weight() * c_dphi * c_dphi.transpose();
        }

        /* RHS, volumetric part. */
        auto BG_row_range = one_range;
        auto BG_col_range = cell_range;

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_range.size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        /* LHS: take basis functions derivatives from degree 1 to K+1 */
        auto MG_rowcol_range = one_range;

        matrix_type MG = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

        matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

        BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
                                take(stiff_mat, BG_row_range, BG_col_range);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = face_quadrature.integrate(msh, fc);

	        matrix_type CG = matrix_type::Zero(BG_row_range.size(),num_cell_dofs);
	        matrix_type FG = matrix_type::Zero(BG_row_range.size(),num_face_dofs);

            //matrix_type m1 = matrix_type::Zeros(BG_row_range.size(), BG_col_range.size());
            for (auto& qp : face_quadpoints)
            {
                auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point(),
                                            one_mindeg, one_maxdeg);
                auto c_phi  = cell_basis.eval_functions(msh, cl, qp.point(),
                                            cell_mindeg, cell_maxdeg);
                // ( \nabla w * n , v_T )
                CG += qp.weight() * c_dphi * n * c_phi.transpose();

                auto f_phi = face_basis.eval_functions(msh, fc, qp.point());

                FG += qp.weight() * c_dphi * n * f_phi.transpose();
            }

            BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) -= CG;
            BG.block(0, current_face_range.min(), BG_row_range.size(),
                    current_face_range.size()) += FG;

        }

        oper  = MG.llt().solve(BG);    // GT
        data  = BG.transpose() * oper;  // A
        return;
    }

    void
    compute_low_order(const mesh_type& msh, const cell_type& cl)
    {
        auto DIM = mesh_type::dimension;
        /* LHS: take vectorial basis functions degree K */
        auto zero_mindeg = 0;
        auto zero_maxdeg = m_degree;
        auto zero_range  = cell_basis.range( zero_mindeg,  zero_maxdeg);

        auto cell_mindeg = 0;
        auto cell_maxdeg = m_degree_cell;
        auto cell_range  = cell_basis.range( cell_mindeg,  cell_maxdeg);

        size_t cbs           = zero_range.size(); // cell_basis are until degree k!!!!
        size_t vec_cbs       = DIM * cbs;
        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);

        matrix_type MG  = matrix_type::Zero(vec_cbs, vec_cbs);

        for (auto& qp : cell_quadpoints)
        {
            auto  c_phi = cell_basis.eval_functions(msh, cl, qp.point(),
                                            zero_mindeg, zero_maxdeg);

            matrix_type   vec_phi = make_vectorial_matrix(c_phi, DIM);
            MG   += qp.weight() * vec_phi.transpose() * vec_phi;
        }

        /* RHS, volumetric part. */
        auto BG_row_size = vec_cbs;
        auto BG_col_size = cell_range.size();

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_range.size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        matrix_type CM = matrix_type::Zero(BG_row_size, BG_col_size);

        for (auto& qp : cell_quadpoints)
        {
            auto  c_dphi = cell_basis.eval_gradients(msh, cl, qp.point(),
                                            cell_mindeg, cell_maxdeg);
            auto  c_phi  = cell_basis.eval_functions(msh, cl, qp.point(),
                                            zero_mindeg, zero_maxdeg);

            matrix_type   vec_phi  = make_vectorial_matrix(c_phi, DIM);

            CM  += qp.weight() * (c_dphi * vec_phi).transpose();
        }

        // RHS
        matrix_type BG = matrix_type::Zero(BG_row_size, dsr.total_size());
        //(\nabla v_T , w)
        BG.block(0, 0, BG_row_size, BG_col_size) = CM;

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = face_quadrature.integrate(msh, fc);

            //matrix_type m1 = matrix_type::Zeros(BG_row_range.size(), BG_col_range.size());
            for (auto& qp : face_quadpoints)
            {
                // ( w * n , v_T )
                auto c_phi = cell_basis.eval_functions(msh, cl, qp.point(),
                                                cell_mindeg, cell_maxdeg);
                // w basis
                vector_type  w_phi   = c_phi.head(zero_range.size());
                matrix_type  vec_phi = make_vectorial_matrix(w_phi, DIM);
                matrix_type  CG = qp.weight() * vec_phi.transpose() * n * c_phi.transpose();
                BG.block(0, 0, BG_row_size, cell_range.size()) -= CG;

                // (  w * n , v_F )
                auto f_phi = face_basis.eval_functions(msh, fc, qp.point());
                matrix_type FG = qp.weight() * vec_phi.transpose() * n * f_phi.transpose() ;
                BG.block(0, current_face_range.min(), BG_row_size,
                                    current_face_range.size()) += FG;
            }
        }

        gradrecoper  = MG.llt().solve(BG);    // GT
        data  = BG.transpose() * gradrecoper;  // A
        return;
    }
};

#if 0
template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType,
                        typename DivCellBasisType, typename DivCellQuadType>
class divergence_reconstruction_nopre
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef CellBasisType                       cell_basis_type;
    typedef CellQuadType                        cell_quadrature_type;
    typedef FaceBasisType                       face_basis_type;
    typedef FaceQuadType                        face_quadrature_type;
    typedef DivCellBasisType                    div_cell_basis_type;
    typedef DivCellQuadType                     div_cell_quadrature_type;

    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    cell_basis_type                             cell_basis;
    cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;

    div_cell_basis_type                         div_cell_basis;
    div_cell_quadrature_type                    div_cell_quadrature;

    size_t                                      m_degree;

public:
    matrix_type     oper;
    matrix_type     data;

    divergence_reconstruction_nopre()
        : m_degree(1)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
        div_cell_basis      = div_cell_basis_type(m_degree);
        div_cell_quadrature = div_cell_quadrature_type(2*m_degree);
    }

    divergence_reconstruction_nopre(size_t degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
        div_cell_basis      = div_cell_basis_type(m_degree);
        div_cell_quadrature = div_cell_quadrature_type(2*m_degree);
    }

    void compute(const mesh_type& msh, const cell_type& cl)
    {
        auto dcbs = div_cell_basis.size();
        matrix_type MD = matrix_type::Zero(dcbs, dcbs);

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto phi_d = div_cell_basis.eval_functions(msh, cl, qp.point());
            for (size_t i = 0; i < dcbs; i++)
                for (size_t j = 0; j < dcbs; j++)
                    MD(i,j) += qp.weight() * mm_prod(phi_d[i], phi_d[j]);
        }

        /* RHS, volumetric part. */
        auto fcs = faces(msh, cl);
        auto num_cell_dofs = cell_basis.range(0, m_degree).size();
        auto num_face_dofs = face_basis.size();
        auto num_faces = fcs.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        matrix_type BD = matrix_type::Zero(dcbs, dsr.total_size());
        for (auto& qp : cell_quadpoints)
        {
            auto phi = cell_basis.eval_functions(msh, cl, qp.point());
            auto dphi_d = div_cell_basis.eval_gradients(msh, cl, qp.point());

            for (size_t i = 0; i < dphi_d.size(); i++)
                for (size_t j = 0; j < num_cell_dofs; j++)
                    BD(i,j) -= qp.weight() * mm_prod(dphi_d[i], phi[j]);
        }

        size_t face_offset = num_cell_dofs;

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = face_quadrature.integrate(msh, fc);

            for (auto& qp : face_quadpoints)
            {
                auto phi_d = div_cell_basis.eval_functions(msh, cl, qp.point());
                auto phi = face_basis.eval_functions(msh, fc, qp.point());
                for (size_t i = 0; i < phi_d.size(); i++)
                {
                    for (size_t j = 0; j < face_basis.size(); j++)
                    {
                        auto p1 = mm_prod(phi[j], n);
                        scalar_type p2 = mm_prod(p1, phi_d[i]);
                        BD(i,face_offset+j) += qp.weight() * p2;
                    }
                }
            }

            face_offset += face_basis.size();
        }

        oper = MD.partialPivLu().solve(BD);
        data = BD.transpose() * oper;
    }
};
#endif
template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class diffusion_like_stabilization_nopre
{
    typedef Mesh                            mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;

    typedef CellBasisType                   cell_basis_type;
    typedef CellQuadType                    cell_quadrature_type;
    typedef FaceBasisType                   face_basis_type;
    typedef FaceQuadType                    face_quadrature_type;

    typedef dynamic_matrix<scalar_type>     matrix_type;
    typedef dynamic_vector<scalar_type>     vector_type;

    cell_basis_type                         cell_basis;
    cell_quadrature_type                    cell_quadrature;

    face_basis_type                         face_basis;
    face_quadrature_type                    face_quadrature;

    size_t                                  m_degree, m_degree_cell;

public:
    matrix_type     low_order_data;
    matrix_type     data;
    matrix_type     oper;
    diffusion_like_stabilization_nopre()
        : m_degree(1)
    {
        m_degree_cell = m_degree;
        //if()        m_degree_cell = m_degree + 1;

        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }
    #if 0
    diffusion_like_stabilization_nopre(size_t degree)
        : m_degree(degree), m_degree_cell(degree)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }
    #endif
    diffusion_like_stabilization_nopre(size_t degree, const bool& change_degree)
        : m_degree(degree),  m_degree_cell(degree)
    {
        if(change_degree)
            m_degree_cell = m_degree + 1;

        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    void compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
    {
        auto zero_mindeg = 0;
        auto zero_maxdeg = m_degree;
        auto zero_range  = cell_basis.range( zero_mindeg,  zero_maxdeg);

        auto cell_mindeg = 0;
        auto cell_maxdeg = m_degree_cell;
        auto cell_range  = cell_basis.range( cell_mindeg,  cell_maxdeg);

        auto one_mindeg = 1;
        auto one_maxdeg = m_degree + 1;
        auto one_range  = cell_basis.range( one_mindeg, one_maxdeg);

        auto pts = points(msh,cl);

        matrix_type mass_mat = matrix_type::Zero(cell_basis.size(), cell_basis.size());

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
            mass_mat  += qp.weight() * c_phi * c_phi.transpose();
        }

        // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

        //Step 1: compute \pi_T^k p_T^k v (third term).
        matrix_type M1 = take(mass_mat, cell_range, cell_range);
        matrix_type M2 = take(mass_mat, cell_range, one_range);
        matrix_type proj1 = -M1.llt().solve(M2 * gradrec_oper);

        //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
        matrix_type I_T = matrix_type::Identity(cell_range.size(), cell_range.size());
        proj1.block(0, 0, cell_range.size(), cell_range.size()) += I_T;

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_range.size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        data = matrix_type::Zero(dsr.total_size(), dsr.total_size());
        oper = matrix_type::Zero(num_face_dofs, num_faces * dsr.total_size());

        //std::cout << "cell: "<< msh.lookup(cl)  << std::endl;
        // Step 3: project on faces (eqn. 21)
        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fbs = current_face_range.size();

            auto vts = msh.get_vertices(cl, pts);
            auto h   = diameter(msh, vts);
            //auto h = diameter(msh, /*fcs[face_i]*/cl);
            auto fc = fcs[face_i];

            matrix_type face_mass_matrix    = matrix_type::Zero(fbs, fbs);
            matrix_type face_trace_matrix   = matrix_type::Zero(fbs, cell_basis.size());

            auto face_quadpoints = face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
                auto f_phi = face_basis.eval_functions(msh, fc, qp.point());
                auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());

                face_mass_matrix  += qp.weight() * f_phi * f_phi.transpose();
                face_trace_matrix += qp.weight() * f_phi * c_phi.transpose();
            }

            Eigen::LLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix);

            // Step 3a:  \pi_F^k( v_F - p_T^k v )
            auto face_range = current_face_range.remove_offset();
            matrix_type MR1 = take(face_trace_matrix, face_range, one_range);
            matrix_type proj2 = piKF.solve(MR1 * gradrec_oper);
            matrix_type I_F = matrix_type::Identity(fbs, fbs);
            auto block_offset = current_face_range.min();
            proj2.block(0, block_offset, fbs, fbs) -= I_F;

            // Step 3b:  \pi_F^k( v_T - \pi_T^k p_T^k v )
            matrix_type MR2   = take(face_trace_matrix, face_range, cell_range);
            matrix_type proj3 = piKF.solve(MR2 * proj1);

            matrix_type BRF = proj2 + proj3;

            oper.block(0, face_i * dsr.total_size(), fbs, dsr.total_size()) = BRF;
            data += BRF.transpose() * face_mass_matrix * BRF / h;
        }
    }

    void
    compute_low_order(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
    {
        auto pts = points(msh,cl);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto zero_mindeg = 0;
        auto zero_maxdeg = m_degree;
        auto zero_range  = cell_basis.range( zero_mindeg, zero_maxdeg);

        auto cell_mindeg = 0;
        auto cell_maxdeg = m_degree_cell;
        auto cell_range  = cell_basis.range( cell_mindeg, cell_maxdeg);

        auto one_mindeg  = 1;
        auto one_maxdeg  = m_degree + 1;
        auto one_range   = cell_basis.range( one_mindeg, one_maxdeg);

        auto num_cell_dofs = cell_range.size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        // Build \pi_F^k (v_F -  v_T) equations (21) and (22)

        //Step 1: store v_T (third term).
        matrix_type proj1 = matrix_type::Zero(cell_range.size() , dsr.total_size());
        matrix_type I_T   = matrix_type::Identity(cell_range.size(), cell_range.size());
        proj1.block(0, 0, cell_range.size(), cell_range.size()) += I_T;

        low_order_data  = matrix_type::Zero(dsr.total_size(), dsr.total_size());
        oper = matrix_type::Zero(num_face_dofs, num_faces * dsr.total_size());

        // Step 3: project on faces (eqn. 21)
        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fbs = current_face_range.size();
            auto vts = msh.get_vertices(cl, pts);
            //auto h = diameter(msh, /*fcs[face_i]*/cl);
            auto h = diameter(msh, /*fcs[face_i]*/vts);
            auto fc = fcs[face_i];

            matrix_type face_mass_matrix    = matrix_type::Zero(fbs, fbs);
            matrix_type face_trace_matrix   = matrix_type::Zero(fbs, cell_basis.size());

            auto face_quadpoints = face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
                auto f_phi = face_basis.eval_functions(msh, fc, qp.point());
                auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());

                face_mass_matrix  += qp.weight() * f_phi * f_phi.transpose();
                face_trace_matrix += qp.weight() * f_phi * c_phi.transpose();
            }

            Eigen::LLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix);

            // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
            auto face_range   = current_face_range.remove_offset();
            matrix_type MR2   = take(face_trace_matrix, face_range, cell_range);
            matrix_type proj3 = piKF.solve(MR2 * proj1);

            // Step 3a: \pi_F^k( v_F - p_T^k v )
            matrix_type I_F = matrix_type::Identity(fbs, fbs);
            auto block_offset = current_face_range.min();
            proj3.block(0, block_offset, fbs, fbs) -= I_F;

            matrix_type BRF = proj3;
            oper.block(0, face_i * dsr.total_size(), fbs, dsr.total_size()) = BRF;
            low_order_data += BRF.transpose() * face_mass_matrix * BRF / h;


        }
        return;
    }

};

template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class diffusion_like_static_condensation_nopre
{
    typedef Mesh                            mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;

    typedef CellBasisType                   cell_basis_type;
    typedef CellQuadType                    cell_quadrature_type;
    typedef FaceBasisType                   face_basis_type;
    typedef FaceQuadType                    face_quadrature_type;

    typedef dynamic_matrix<scalar_type>     matrix_type;
    typedef dynamic_vector<scalar_type>     vector_type;

    cell_basis_type                         cell_basis;
    cell_quadrature_type                    cell_quadrature;

    face_basis_type                         face_basis;
    face_quadrature_type                    face_quadrature;

    size_t                                  m_degree, m_degree_cell;

public:
    diffusion_like_static_condensation_nopre()
        : m_degree(1)
    {
        m_degree_cell = m_degree;

        cell_basis          = cell_basis_type(m_degree_cell);
        cell_quadrature     = cell_quadrature_type(2*m_degree_cell);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }
    #if 0
    diffusion_like_static_condensation_nopre(size_t degree)
        : m_degree(degree), m_degree_cell(degree)
    {
        cell_basis          = cell_basis_type(m_degree_cell);
        cell_quadrature     = cell_quadrature_type(2*m_degree_cell);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }
    #endif

    diffusion_like_static_condensation_nopre(size_t degree, const bool& change_degree)
        : m_degree(degree),  m_degree_cell(degree)
    {
        if(change_degree)
            m_degree_cell = m_degree + 1;

        cell_basis          = cell_basis_type(m_degree_cell);
        cell_quadrature     = cell_quadrature_type(2*m_degree_cell);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }


    auto
    compute(const mesh_type& msh,
            const cell_type& cl,
            const matrix_type& local_mat,
            const vector_type& cell_rhs)
    {
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_basis.size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        size_t cell_size = dsr.cell_range().size();
        size_t face_size = dsr.all_faces_range().size();

        matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
        matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
        matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
        matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);

        assert(K_TT.cols() == cell_size);
        assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
        assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
        assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
        assert(K_FT.cols() + K_FF.cols() == local_mat.cols());

        auto K_TT_ldlt = K_TT.llt();
        matrix_type AL = K_TT_ldlt.solve(K_TF);
        vector_type bL = K_TT_ldlt.solve(cell_rhs);

        matrix_type AC = K_FF - K_FT * AL;
        vector_type bC = /* no projection on faces, eqn. 26*/ - K_FT * bL;

        return std::make_pair(AC, bC);
    }

    vector_type
    recover(const mesh_type& msh,
            const cell_type& cl,
            const matrix_type& local_mat,
            const vector_type& cell_rhs,
            const vector_type& solF)
    {
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_basis.size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        size_t cell_size        = dsr.cell_range().size();
        size_t all_faces_size   = dsr.all_faces_range().size();

        vector_type ret( dsr.total_size() );

        matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
        matrix_type K_TF = local_mat.topRightCorner(cell_size, all_faces_size);

        vector_type solT = K_TT.llt().solve(cell_rhs - K_TF * solF);

        ret.head(cell_size)         = solT;
        ret.tail(all_faces_size)    = solF;

        return ret;
    }
};

template<typename Mesh, //typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class assembler_nopre
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    //typedef CellBasisType                       cell_basis_type;
    //typedef CellQuadType                        cell_quadrature_type;
    typedef FaceBasisType                       face_basis_type;
    typedef FaceQuadType                        face_quadrature_type;

    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    //cell_basis_type                             cell_basis;
    //cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;

    size_t                                      m_degree;

    typedef Eigen::SparseMatrix<scalar_type>    sparse_matrix_type;
    typedef Eigen::Triplet<scalar_type>         triplet_type;

    std::vector<triplet_type>                   m_triplets;
    size_t                                      m_num_unknowns;

public:

    sparse_matrix_type      matrix;
    vector_type             rhs;

    assembler_nopre()               = delete;
    assembler_nopre(const assembler_nopre&) = delete;
    assembler_nopre(assembler_nopre&&)      = delete;

    assembler_nopre(const mesh_type& msh, const size_t degree)
        : m_degree(degree)
    {
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);

        m_num_unknowns = face_basis.size() * (msh.faces_size() + msh.boundary_faces_size());
        matrix = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
        rhs = vector_type::Zero(m_num_unknowns);
    }

    template<typename LocalContrib>
    void
    assemble(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc)
    {
        auto fcs = faces(msh, cl);
        std::vector<size_t> l2g(fcs.size() * face_basis.size());
        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            auto face_offset = face_id * face_basis.size();

            auto pos = face_i * face_basis.size();

            for (size_t i = 0; i < face_basis.size(); i++)
                l2g[pos + i] = face_offset + i;
        }

        assert(lc.first.rows() == lc.first.cols());
        assert(lc.first.rows() == lc.second.size());
        assert(lc.second.size() == l2g.size());

        //std::cout << lc.second.size() << " " << l2g.size() << std::endl;

        for (size_t i = 0; i < lc.first.rows(); i++)
        {
            for (size_t j = 0; j < lc.first.cols(); j++)
            {
                //std::cout <<  l2g.at(i) <<","<< l2g.at(j)<<","<< lc.first(i,j)  << std::endl;
                m_triplets.push_back( triplet_type( l2g.at(i), l2g.at(j), lc.first(i,j) ) );
            }

            rhs(l2g.at(i)) += lc.second(i);
        }
    }

    template<typename Function>
    void
    impose_boundary_conditions(mesh_type& msh, const Function& bc)
    {
        size_t fbs = face_basis.size();
        size_t face_i = 0;
        for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
        {
            auto bfc = *itor;

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            auto face_offset = face_id * fbs;
            auto face_offset_lagrange = (msh.faces_size() + face_i) * fbs;

            auto fqd = face_quadrature.integrate(msh, bfc);
            //auto binfo = bfc.boundary_info();

            matrix_type MFF     = matrix_type::Zero(fbs, fbs);
            vector_type rhs_f   = vector_type::Zero(fbs);

            for (auto& qp : fqd)
            {
                auto f_phi = face_basis.eval_functions(msh, bfc, qp.point());

                MFF   += qp.weight() * f_phi * f_phi.transpose();
                rhs_f += qp.weight() * f_phi * bc(qp.point(), 0);
            }

            for (size_t i = 0; i < MFF.rows(); i++)
            {
                for (size_t j = 0; j < MFF.cols(); j++)
                {
                    m_triplets.push_back( triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i,j)) );
                    m_triplets.push_back( triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i,j)) );
                }
                rhs(face_offset_lagrange+i) = rhs_f(i);
            }

            face_i++;
        }
    }

    void
    finalize()
    {
        matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
        m_triplets.clear();
    }
};

template< typename T, size_t DIM, typename Storage>
auto
high_order_reconstruction_nopre(const mesh<T, DIM, Storage>&  msh,
                                const typename mesh<T, DIM, Storage>::cell& cl,
                                const dynamic_matrix<T>& rec_oper,
                                const dynamic_vector<T>& v,
                                const size_t& degree)
{
    typedef mesh<T, DIM, Storage>                           mesh_type;
    typedef typename mesh_type::cell                        cell_type;
    typedef typename mesh_type::scalar_type                 scalar_type;
    typedef disk::quadrature<mesh_type, cell_type>          cell_quadrature_type;
    typedef dynamic_matrix<scalar_type>                     matrix_type;
    typedef dynamic_vector<scalar_type>                     vector_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;

    cell_basis_type          cell_basis(degree+1);
    auto cell_quadrature     = cell_quadrature_type(2*(degree+1));
    auto cbs = cell_basis.size();

    matrix_type mm = matrix_type::Zero(cbs, cbs);
    auto cell_quadpoints = cell_quadrature.integrate(msh, cl);

    for (auto& qp : cell_quadpoints)
    {
        auto phi = cell_basis.eval_functions(msh, cl, qp.point());

        mm  += qp.weight() * phi * phi.transpose();
    }

    // Use eqn. (22) to do high order reconstruction.
    auto zero_range = cell_basis.range(0, degree);
    auto one_range = cell_basis.range(1, degree+1);

    vector_type P = vector_type::Zero(cbs);
    vector_type vT = v.head(zero_range.size());

    P.tail(one_range.size()) = rec_oper * v;

    matrix_type M1 = take(mm, zero_range, zero_range);//cell_mass_matrix.block(0, 0, basis_k_size, basis_k_size);
    matrix_type M2 = take(mm, zero_range, one_range);//cell_mass_matrix.block(0, 0, basis_k_size, cbs);
    matrix_type R = vT - M1.ldlt().solve(M2*P);

    P.head(zero_range.size()) += R;

    return P;
};

} // namespace disk
