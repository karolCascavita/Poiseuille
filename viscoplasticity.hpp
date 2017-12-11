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
#include "timecounter.h"

namespace disk {

template<typename T>
struct plasticity_data
{
    plasticity_data(): Lref(1.), Vref(1.), Bn(0.), mu(1.), alpha(1.), beta(1.), f(1.),
                        method(false)
    {}
    T f;                    //WK: Cuidado porque f deberia ser el valor externo de la fuente.
    T Lref;                 /* Charactetistic length */
    T Vref;                 /* Reference velocity */
    T Bn;                   /* Bingham number */
    T mu;                   /* viscosity */
    T alpha;
    T beta;
    T yield;
    bool   method;

    friend std::ostream& operator<<(std::ostream& os, const plasticity_data<T>& p){
        os << "Plastic data: "<<std::endl;
        os << "* method : "<< p.method<< std::endl;
        os << "* f      : "<< p.f<< std::endl;
        os << "* Lref   : "<< p.Lref<< std::endl;
        os << "* Vref   : "<< p.Vref<< std::endl;
        os << "* yield  : "<< p.yield<< std::endl;
        os << "* mu     : "<< p.mu<< std::endl;
        os << "* Bingham: "<< p.Bn<< std::endl;
        os << "* alpha  : "<< p.alpha<< std::endl;
        os << "* beta   : "<< p.beta<< std::endl;

        return os;
    }
};

template<typename MeshType, typename ElementType>
class tensor_bones
{
    typedef ElementType                         element_type;
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;

    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    vector_type vec;
    size_t m_degree;
};


template<typename T, typename Storage>
class tensor_bones<disk::mesh<T,2,Storage>, typename disk::mesh<T,2,Storage>::face>
{
public:
    typedef disk::mesh<T,2,Storage>                     mesh_type;
    typedef typename mesh_type::face                    face_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quad_type;
    typedef dynamic_matrix<T>         matrix_type;
    typedef dynamic_vector<T>         vector_type;
    typedef T scalar_type;

    vector_type vec;
    size_t m_degree;

    //WK: try to do this in the constructor.(I had problems in tensor_zero_vector,
    // since it s asking for the input of the constructor)
    tensor_bones(){};
    tensor_bones(const mesh_type& msh,
                 const face_type& face,
                 const size_t& degree,
                 const size_t& num_dofs): m_degree(degree)
    {
        vec = vector_type::Zero(num_dofs);
    };

    size_t      degree()        { return m_degree;}

    friend tensor_bones<mesh_type, face_type> operator-
                            (const tensor_bones<mesh_type, face_type>& tp1,
                            const tensor_bones< mesh_type, face_type>& tp2)
    {
        tensor_bones<mesh_type, face_type> ret;
        if(tp1.m_degree != tp2.m_degree)
            throw std::logic_error("(degree) is different for these tensors, so it is not possible to do (-)operation.");
        ret.m_degree  = tp1.m_degree;
        ret.vec  = tp1.vec - tp2.vec;

        return ret;
    }
    friend tensor_bones<mesh_type, face_type> operator+
                            (const tensor_bones<mesh_type, face_type>& tp1,
                            const tensor_bones<mesh_type, face_type>& tp2)
    {
        tensor_bones<mesh_type, face_type> ret;
        if(tp1.m_degree != tp2.m_degree)
            throw std::logic_error("(degree) is different for these tensors, so it is not possible to do (-)operation.");
        ret.m_degree  = tp1.m_degree;
        ret.vec = tp1.vec + tp2.vec;

        return ret;
    }
    friend tensor_bones<mesh_type, face_type> operator*(const scalar_type& scalar,
                            const tensor_bones<mesh_type, face_type>& tp)
    {
        tensor_bones<mesh_type, face_type> ret;
        ret.m_degree  = tp.m_degree;
        ret.vec = scalar * tp.vec ;

        return ret;
    }

    #if 0
    friend vector_type operator*(const vector_type& vec,
                            const tensor_bones<mesh_type, face_type>& tp)
    {
        vector_type ret = vector_type::Zero( tp.num_eval_pts);
        ret.m_quad_degree  = tp.m_quad_degree;
        ret.num_eval_pts   = tp.num_eval_pts;

        size_t i, j;
        for( i = 0; i < tp.quad_evals_mat.cols(); i++)
            ret(i)= make_prod(tp.quad_evals_mat, vec);

        for(j= 0; j < tp.nodes_evals_mat.cols(); j++)
            ret(j + i)= make_prod(tp.nodes_evals_mat, vec);

        ret.bar_eval_vec(j + i + 1)   = make_prod( vec,tp.bar_eval_vec);

        return ret;
    }
    #endif
    friend std::ostream& operator<<(std::ostream& os,
                                const tensor_bones<mesh_type, face_type>& tp)
    {
        os << " * Tensor at face "<<std::endl;
        os << " * degree : "<< tp.m_degree<<std::endl;
        os << " * vector : size = "<< tp.vec.rows()<< std::endl;

        for(size_t i = 0 ; i < tp.vec.rows(); i++)
            os<< "  " << tp.vec(i);
        return os;
    }

    template<typename FaceBasisType>
    void
    Zero()
    {
        FaceBasisType fb(m_degree);
        auto DIM = mesh_type::dimension;
        vec = vector_type::Zero(DIM, fb.size());
        return;
    }
    #if 0
    matrix_type
    get_all()
    {
        matrix_type ret = matrix_type::Zero(2, num_eval_pts);

        ret.block(0, 0, 2, num_eval_quad_pts) = quad_evals_mat;
        ret.block(0, num_eval_quad_pts, 2, 1) = bar_eval_vec;

        return ret;
    }
    #endif
};

template<typename T, typename Storage>
class tensor_bones<disk::mesh<T,2,Storage>, typename disk::mesh<T,2,Storage>::cell>
{
public:
    typedef disk::mesh<T,2,Storage>                     mesh_type;
    typedef typename mesh_type::cell                    cell_type;
    typedef disk::quadrature<mesh_type, cell_type>      cell_quad_type;
    typedef dynamic_matrix<T>         matrix_type;
    typedef dynamic_vector<T>         vector_type;
    typedef T scalar_type;

    vector_type vec;
    size_t m_degree;

    //WK: try to do this in the constructor.(I had problems in tensor_zero_vector,
    // since it s asking for the input of the constructor)
    tensor_bones(){};
    tensor_bones(const mesh_type& msh,
                 const cell_type& cell,
                 const size_t& degree,
                 const size_t& num_dofs): m_degree(degree)
    {
        vec = vector_type::Zero(num_dofs);
    };

    size_t      degree()        { return m_degree;}

    friend tensor_bones<mesh_type, cell_type> operator-
                            (const tensor_bones<mesh_type, cell_type>& tp1,
                            const tensor_bones< mesh_type, cell_type>& tp2)
    {
        tensor_bones<mesh_type, cell_type> ret;
        if(tp1.m_degree != tp2.m_degree)
            throw std::logic_error("(degree) is different for these tensors, so it is not possible to do (-)operation.");
        ret.m_degree  = tp1.m_degree;
        ret.vec  = tp1.vec - tp2.vec;

        return ret;
    }
    friend tensor_bones<mesh_type, cell_type> operator+
                            (const tensor_bones<mesh_type, cell_type>& tp1,
                            const tensor_bones<mesh_type, cell_type>& tp2)
    {
        tensor_bones<mesh_type, cell_type> ret;
        if(tp1.m_degree != tp2.m_degree)
            throw std::logic_error("(degree) is different for these tensors, so it is not possible to do (-)operation.");
        ret.m_degree  = tp1.m_degree;
        ret.vec = tp1.vec + tp2.vec;

        return ret;
    }
    friend tensor_bones<mesh_type, cell_type> operator*(const scalar_type& scalar,
                            const tensor_bones<mesh_type, cell_type>& tp)
    {
        tensor_bones<mesh_type, cell_type> ret;
        ret.m_degree  = tp.m_degree;
        ret.vec = scalar * tp.vec ;

        return ret;
    }
    friend std::ostream& operator<<(std::ostream& os,
                                const tensor_bones<mesh_type, cell_type>& tp)
    {
        os << " * Tensor at cell "<<std::endl;
        os << " * degree : "<< tp.m_degree<<std::endl;
        os << " * vector : size = "<< tp.vec.rows()<< std::endl;

        for(size_t i = 0 ; i < tp.vec.rows(); i++)
            os<< "  " << tp.vec(i);
        return os;
    }
};


template<typename MeshType>
class tensors
{
public:

    typedef MeshType   mesh_type;
    typedef typename mesh_type::cell      cell_type;
    typedef typename mesh_type::face      face_type;
    typedef typename mesh_type::scalar_type    scalar_type;
    typedef dynamic_matrix<scalar_type>        matrix_type;
    typedef dynamic_vector<scalar_type>        vector_type;
    typedef std::vector<vector_type>           vec_vec_type;
    typedef std::vector<matrix_type>           vec_mat_type;
    vec_vec_type     cell_tensor_vector;
    vec_mat_type     c_face_tensor_vector;

    tensors()
    {}

    tensors(const MeshType& msh)
    {
        cell_tensor_vector = vec_vec_type(msh.cells_size());
        c_face_tensor_vector = vec_mat_type(msh.cells_size());
    }

    tensors(const size_t cells_size)
    {
        cell_tensor_vector = vec_vec_type(cells_size);
        c_face_tensor_vector = vec_mat_type(cells_size);
    }
    template<typename CellBasisType, typename FaceBasisType>
    void
    zero_vector(const mesh_type msh, const size_t& degree)
    {
        CellBasisType cell_basis(degree);
        FaceBasisType face_basis(degree);

        auto DIM = mesh_type::dimension;
        size_t num_cell_dofs = DIM * cell_basis.size();
        size_t num_face_dofs = face_basis.size();

        cell_tensor_vector = vec_vec_type(msh.cells_size());
        c_face_tensor_vector = vec_mat_type(msh.cells_size());

        for (auto& cl : msh)
        {
            auto fcs = faces(msh, cl);
            auto num_faces = number_of_faces(msh, cl);
            auto id = cl.get_id();
            cell_tensor_vector.at(id)   = vector_type::Zero(num_cell_dofs);
            c_face_tensor_vector.at(id) = matrix_type::Zero(num_face_dofs, num_faces);
        }

        return;
    };
    vec_vec_type  at_all_cells(void) {  return cell_tensor_vector;}
    vec_mat_type  at_all_faces(void) {  return c_face_tensor_vector;}

    vec_vec_type  at_all_cells(void) const {  return cell_tensor_vector;}
    vec_mat_type  at_all_faces(void) const {  return c_face_tensor_vector;}

    typedef typename vec_vec_type::iterator         cell_iterator;
    typedef typename vec_vec_type::const_iterator   const_cell_iterator;

    cell_iterator           at_cells_begin() { return cell_tensor_vector.begin(); }
    cell_iterator           at_cells_end()   { return cell_tensor_vector.end(); }
    const_cell_iterator     at_cells_begin() const { return cell_tensor_vector.begin(); }
    const_cell_iterator     at_cells_end()   const { return cell_tensor_vector.end(); }

    /* face iterators */
    typedef typename vec_mat_type::iterator        face_iterator;
    typedef typename vec_mat_type::const_iterator  const_face_iterator;

    face_iterator           at_faces_begin() { return c_face_tensor_vector.begin(); }
    face_iterator           at_faces_end()   { return c_face_tensor_vector.end(); }
    const_face_iterator     at_faces_begin() const { return c_face_tensor_vector.begin(); }
    const_face_iterator     at_faces_end()   const { return c_face_tensor_vector.end(); }

    size_t  at_all_cells_size() const { return cell_tensor_vector.size(); }
    size_t  at_all_faces_size() const { return c_face_tensor_vector.size(); }

    vector_type
    at_cell(const mesh_type & msh, const cell_type& cl)
    {
        auto id = msh.lookup(cl);
        return cell_tensor_vector.at(id);
    };

    vector_type
    at_cell(const mesh_type & msh, const cell_iterator & itor)
    {
        auto id = msh.lookup(*itor);
        return cell_tensor_vector.at(id);
    };

    vector_type
    at_cell(const mesh_type & msh, const cell_type& cl) const
    {
        auto id = msh.lookup(cl);
        return cell_tensor_vector.at(id);
    };

    vector_type
    at_cell(const mesh_type & msh, const cell_iterator & itor)  const
    {
        auto id = msh.lookup(*itor);
        return cell_tensor_vector.at(id);
    };

    vector_type
    at_face(const mesh_type & msh, const cell_type& cl, const face_type& fc)
    {
        auto id  = msh.lookup(cl);
        auto pos = face_position(msh, cl, fc);
        return c_face_tensor_vector.at(id).col(pos);
    };

    vector_type
    at_face(const mesh_type & msh, const cell_type& cl, const face_type& fc)
    const
    {
        auto id  = msh.lookup(cl);
        auto pos = face_position(msh, cl, fc);
        return c_face_tensor_vector.at(id).col(pos);
    };

    matrix_type
    at_element_faces(const mesh_type & msh, const cell_type& cl)
    {
        auto id  = msh.lookup(cl);
        return c_face_tensor_vector.at(id);
    };

    matrix_type
    at_element_faces(const mesh_type & msh, const cell_type& cl)
    const
    {
        auto id  = msh.lookup(cl);
        return c_face_tensor_vector.at(id);
    };

    vector_type
    at_face(const mesh_type & msh, const cell_type& cl, const face_iterator & itor)
    {
        auto fc  = *itor;
        auto id  = msh.lookup(id);
        auto pos = face_position(msh, cl, fc);
        return c_face_tensor_vector.at(id).col(pos);
    };

    void
    save(const mesh_type & msh, const cell_type & cl, const vector_type& vec)
    {
        auto id = msh.lookup(cl);
        cell_tensor_vector.at(id) = vec;
        return;
    }
    void
    save(const mesh_type & msh, const cell_type& cl, const face_type & fc,
         const vector_type& vec)
    {
        auto id  = msh.lookup(cl);
        auto pos = face_position(msh, cl, fc);
        c_face_tensor_vector.at(id).col(pos) = vec;
        return;
    }
    void
    instantiate(const size_t i)
    {
        vec_vec_type().swap(cell_tensor_vector);
        vec_mat_type().swap(c_face_tensor_vector);
        cell_tensor_vector = vec_vec_type(i);
        c_face_tensor_vector = vec_mat_type(i);
        return;
    }
};

template<typename T,typename Mesh,typename CellBasisType, typename CellQuadType,
                       typename FaceBasisType, typename FaceQuadType, typename Tensors>
class plasticity2
{
    typedef Mesh                            mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef dynamic_matrix<scalar_type>     matrix_type;
    typedef dynamic_vector<scalar_type>     vector_type;
    typedef Tensors                         tensors_type;

    typedef CellBasisType                   cell_basis_type;
    typedef FaceBasisType                   face_basis_type;
    typedef CellQuadType                    cell_quadrature_type;
    typedef FaceQuadType                    face_quadrature_type;

    typedef std::vector<vector_type>           vec_vec_type;

    cell_basis_type                         cell_basis;
    face_basis_type                         face_basis;
    cell_quadrature_type                    cell_quadrature;
    face_quadrature_type                    face_quadrature;


    plasticity_data<T>          pst;
    size_t                      m_degree, m_degree_cell;
    T                           m_method_coef;
    std::vector<matrix_type>    m_gradrec_opers;
    std::vector<matrix_type>    m_stabilization_opers;

public:

    tensors_type                m_multiplier;
    tensors_type                m_decoupled_var;

    plasticity2()
    {}
    plasticity2(const size_t& degree,
                const plasticity_data<T>& pst_data,
                const std::vector<matrix_type>& gradrec_opers,
                const std::vector<matrix_type>& stabilization_opers,
                const bool& change_degree):
                pst(pst_data),
                m_degree(degree),
                m_degree_cell(degree),
                m_gradrec_opers(gradrec_opers),
                m_stabilization_opers(stabilization_opers)
                //WK: for mu as tensor one should check this part
    {
        if(change_degree)
            m_degree_cell = m_degree + 1;

        cell_basis  =  cell_basis_type(m_degree + 1);
        face_basis  =  face_basis_type(m_degree);

        cell_quadrature  = cell_quadrature_type(2 * degree);//+ 2);
        face_quadrature  = face_quadrature_type(2 * degree);

        if(pst.method == true)
            m_method_coef = 1. / pst.alpha ;
        else
            m_method_coef = 1. / (pst.alpha + pst.mu);
    }

    template<typename Errors>
    void
    update_multiplier(const mesh_type    & msh,
                         tensors_type       & mtp,
                         const std::vector<vector_type>& velocity,
                         Errors       & error)
    {
        auto zero_mindeg = 0;
        auto zero_maxdeg = m_degree;
        auto zero_range  = cell_basis.range( zero_mindeg,  zero_maxdeg);

        for(auto& cl:msh)
        {
            auto id     = msh.lookup(cl);
            vector_type vel_TF     = velocity.at(id);
            vector_type Gt_uT      = grad_rec_uTF(msh, cl, vel_TF );
            vector_type gamma      = m_decoupled_var.at_cell(msh, cl);
            vector_type old_sigma  = mtp.at_cell(msh, cl);
            vector_type diff_sigma = pst.alpha * ( Gt_uT - gamma );
            vector_type sigma      = old_sigma + diff_sigma;

            //Updating mtp in the cell
            mtp.save(msh, cl, sigma);

            //Error computation for sigma (cell)
            auto DIM = mesh_type::dimension;
            auto cbs = zero_range.size();
            auto vec_cbs = DIM * zero_range.size();


            matrix_type mass_mat  = matrix_type::Zero(vec_cbs, vec_cbs);

            auto row_range = disk::dof_range( 0, DIM);
            auto cell_quadpoints   = cell_quadrature.integrate(msh, cl);
            for (auto& qp : cell_quadpoints)
            {
                auto c_phi = cell_basis.eval_functions(msh, cl, qp.point(),
                                                    zero_mindeg, zero_maxdeg);
                matrix_type vec_phi = make_vectorial_matrix(c_phi, DIM);

                mass_mat += qp.weight() * vec_phi.transpose() * vec_phi;
            }

            error.sigma += diff_sigma.transpose() * mass_mat * diff_sigma;

            // Updating mtp on the faces
            auto fbs = face_basis.size();
            auto fcs = faces(msh , cl);

            matrix_type   face_mass_mat = matrix_type::Zero(fbs, fbs);

            for (auto& fc :fcs)
            {
                vector_type STF_u = stab_oper(msh, cl, fc) * vel_TF ;
                vector_type Psi   = m_decoupled_var.at_face(msh, cl, fc);
                vector_type old_varsigma  = mtp.at_face(msh, cl, fc);
                vector_type diff_varsigma = pst.alpha * ( STF_u - Psi );
                vector_type varsigma =  old_varsigma + diff_varsigma;

                //Error computation for varsigma (cells)
                auto face_quadpoints = face_quadrature.integrate(msh, fc);
                for (auto& qp : face_quadpoints)
                {
                    auto f_phi = face_basis.eval_functions(msh, fc, qp.point());
                    face_mass_mat  += qp.weight() * f_phi * f_phi.transpose();
                }

                error.varsigma += diff_varsigma.transpose() * face_mass_mat * diff_varsigma;

                mtp.save( msh, cl, fc, varsigma);
            }
        }
        return;
    }

    void
    eval_plasticity(const mesh_type  & msh,
                  const cell_type    & cl,
                  const face_type    & fc,
                  const vector_type  & varsigma,
                  const vector_type  & velocity)
    {
        vector_type STF_u = stab_oper(msh, cl, fc) * velocity;

        auto fbs = face_basis.size();
        matrix_type   face_mass_matrix = matrix_type::Zero(fbs, fbs);
        vector_type   lhs_vec = vector_type::Zero(fbs);

        auto pts = points(msh, cl);
        auto vts = msh.get_vertices(cl, pts);
        auto h   = diameter(msh, /*fcs[face_i]*/vts);

        auto face_quadpoints = face_quadrature.integrate(msh, fc);
        for (auto& qp : face_quadpoints)
        {
            auto f_phi = face_basis.eval_functions(msh, fc, qp.point());

            face_mass_matrix  += qp.weight() * f_phi * f_phi.transpose();

            // Eval point-wise gamma function at test points
            vector_type theta = varsigma  +  pst.alpha * STF_u / h ;
            T     theta_eval  =  f_phi.dot(theta);
            auto  theta_eval_norm  =  theta_eval;

            T Psi_eval = T(0);
            if(theta_eval_norm > pst.yield)
                Psi_eval = m_method_coef * (theta_eval_norm - pst.yield) *
                                                (theta_eval / theta_eval_norm);
            lhs_vec += qp.weight() * f_phi * Psi_eval;
        }

        //Recover Psi function: DOFs in  P^(k)_(d-1)
        vector_type Psi = face_mass_matrix.llt().solve(lhs_vec);
        m_decoupled_var.save(msh, cl, fc, Psi);
        return;
    }
    void
    eval_plasticity(const mesh_type    & msh,
                  const cell_type    & cl,
                  const vector_type  & sigma,
                  const vector_type& velocity)
    {
        auto id = msh.lookup(cl);
        vector_type Gt_uT  = grad_rec_uTF(msh, cl, velocity);

        auto zero_mindeg = 0;
        auto zero_maxdeg = m_degree;
        auto zero_range = cell_basis.range( zero_mindeg,  zero_maxdeg);

        auto DIM = mesh_type::dimension;
        auto row_range = disk::dof_range( 0, DIM);
        auto cbs = zero_range.size();
        auto vec_cbs = DIM * zero_range.size();
        size_t i = 0;

        vector_type lhs_vec    = vector_type::Zero(vec_cbs);
        matrix_type mass_mat   = matrix_type::Zero(vec_cbs, vec_cbs);

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto c_phi = cell_basis.eval_functions(msh, cl, qp.point(),
                                            zero_mindeg, zero_maxdeg);
            matrix_type vec_phi = make_vectorial_matrix(c_phi, DIM);

            mass_mat += qp.weight() * vec_phi.transpose() * vec_phi;

            // Eval point-wise gamma function at test points
            vector_type theta      = sigma  +  pst.alpha * Gt_uT;
            vector_type theta_eval = vec_phi * theta;
            auto  theta_eval_norm  = theta_eval.norm();

            vector_type gamma_eval = vector_type::Zero(mesh_type::dimension);
            if(theta_eval_norm > pst.yield)
            {
                gamma_eval = m_method_coef * (theta_eval_norm - pst.yield) *
                                                    (theta_eval / theta_eval_norm);
            }
            lhs_vec += qp.weight() * vec_phi.transpose() * gamma_eval;

        }

        //Recover gamma function: DOFs in  [P^k]^d
        vector_type gamma = mass_mat.llt().solve(lhs_vec);
        m_decoupled_var.save( msh, cl, gamma);

        return;
    }
    template<typename Timing>
    void
    compute_decoupling(const mesh_type     & msh,
                       const tensors_type  & multiplier,
                       const vec_vec_type  & velocity,
                        Timing & timings)
    {
        m_decoupled_var.instantiate(msh.cells_size());
        m_decoupled_var.template zero_vector<CellBasisType, FaceBasisType>(msh, m_degree);

        for(auto& cl : msh)
        {
            //timecounter tc_detail;
            //tc_detail.tic();
            auto id  = msh.lookup(cl);
            vector_type sigma   = multiplier.at_cell( msh, cl);
            vector_type vel_TF  = velocity.at(id);
            //tc_detail.toc();
            //timings["Plasticity cells - grbg"] += tc_detail.to_double();

            //tc_detail.tic();
            eval_plasticity( msh, cl, sigma, vel_TF);
            //tc_detail.toc();
            //timings["Plasticity cells"] += tc_detail.to_double();

            auto fcs = faces(msh, cl);

            for (auto & fc :fcs)
            {
                //tc_detail.tic();
                vector_type varsigma = multiplier.at_face( msh, cl, fc );
                //tc_detail.toc();
                //timings["Plasticity faces - grbg"] += tc_detail.to_double();

                //tc_detail.tic();
                eval_plasticity(msh, cl, fc, varsigma, vel_TF);
                //tc_detail.toc();
                //timings["Plasticity faces"] += tc_detail.to_double();
            }
        }
        return;
    }

    vector_type
    compute_integral(const mesh_type     & msh,
                     const cell_type     & cl,
                     const tensors_type  & multiplier)
    {
        auto num_faces  = number_of_faces(msh, cl);
        auto cell_range = cell_basis.range(0, m_degree_cell);
        auto face_range = face_basis.range(0, m_degree);

        dofspace_ranges dsr(cell_range.size(), face_range.size(), num_faces);
        auto a_faces_range = dsr.all_faces_range();
        vector_type ret = vector_type::Zero(dsr.total_size());

        ret  = compute_cell_integral( msh, cl, multiplier);
        ret += compute_faces_integral(msh, cl, dsr, multiplier);
        return ret;
    }

    vector_type
    compute_cell_integral(const mesh_type    & msh,
                          const cell_type    & cl,
                          const tensors_type & multiplier)
    {

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);

        auto zero_mindeg = 0;
        auto zero_maxdeg = m_degree;
        auto DIM = mesh_type::dimension;
        auto cbs = cell_basis.range(zero_mindeg, zero_maxdeg).size();
        auto vec_cbs = DIM * cbs;

        matrix_type mass_mat = matrix_type::Zero(vec_cbs, vec_cbs);
        vector_type sigma   = multiplier.at_cell( msh, cl );
        vector_type gamma   = m_decoupled_var.at_cell( msh, cl );

        for (auto& qp : cell_quadpoints)
        {
            auto c_phi = cell_basis.eval_functions(msh, cl, qp.point(),
                                                zero_mindeg, zero_maxdeg );
            matrix_type vec_phi = make_vectorial_matrix(c_phi, DIM);
            mass_mat += qp.weight() * vec_phi.transpose
            () * vec_phi;
        }

        auto id = msh.lookup(cl);

        matrix_type rec_oper = m_gradrec_opers.at(id);
        matrix_type RMG = rec_oper.transpose() * mass_mat;

        return RMG * (sigma - pst.alpha * gamma);
    }

    vector_type
    compute_faces_integral( const mesh_type   & msh,
                            const cell_type   & cl,
                            const dofspace_ranges & dsr,
                            const tensors_type    & multiplier)
    {
        auto pts = points(msh, cl);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        vector_type ret = vector_type::Zero(dsr.total_size());


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

            auto face_quadpoints = face_quadrature.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
                auto f_phi = face_basis.eval_functions(msh, fc, qp.point());

                face_mass_matrix  += qp.weight() * f_phi * f_phi.transpose();
            }

            vector_type varsigma = multiplier.at_face(msh, cl, fc);
            vector_type Psi  = m_decoupled_var.at_face(msh, cl, fc);
            matrix_type STF  = stab_oper(msh, cl, fc);

            ret += STF.transpose() * face_mass_matrix *
                                    (varsigma - (pst.alpha / h) * Psi );
        }
        return ret;
    }

    #if 0
    vector_type
    grad_rec_uTF_eval( const point<T,2>& pt, //Check for this to be generalized to other dimensions.
                   const mesh_type&  msh,
                   const cell_type&  cl,
                   const vector_type& u_TF)
    {
        auto id = cl.get_id();
        matrix_type  rec_op = m_gradrec_opers.at(id);

        auto col_range      = cell_basis.range(1,m_degree+1);
        auto row_range      = disk::dof_range( 0,mesh_type::dimension);

        auto dphi   = cell_basis.eval_gradients(msh, cl, pt);
        matrix_type dphi_matrix = make_gradient_matrix(dphi, row_range, col_range);
        vector_type dphi_ruh    = dphi_matrix * rec_op * u_TF;

        return dphi_ruh;
    }

    std::pair<vector_type, vector_type>
    grad_uh_eval( const point<T,2>& pt, //Check for this to be generalized to other dimensions.
                   const mesh_type&  msh,
                   const cell_type&  cl,
                   const vector_type& u_TF)
    {
        auto id = msh.lookup(cl);
        matrix_type  rec_op = m_gradrec_opers.at(id);

        auto col_range = cell_basis.range(1,m_degree+1);
        auto row_range = disk::dof_range( 0,mesh_type::dimension);

        auto cell = *std::next(msh.cells_begin(), id);
        auto dphi = cell_basis.eval_gradients(msh, cell, pt);
        matrix_type dphi_matrix = make_gradient_matrix(dphi);
        matrix_type dphi_taken  = take(dphi_matrix, row_range, col_range);
        vector_type dphi_ruh    = dphi_taken * rec_op * u_TF;

        auto cell_range     = cell_basis.range(0,m_degree);
        vector_type  uh_T   = u_TF.head(cell_range.size());
        matrix_type  dphi_T = take(dphi_matrix, row_range, cell_range);
        vector_type  dphi_uh_T = dphi_T * uh_T;

        return std::make_pair(dphi_ruh, dphi_uh_T);
    }
    #endif

    vector_type
    grad_rec_uTF_eval( const point<T,2>& pt, //Check for this to be generalized to other dimensions.
                   const mesh_type&  msh,
                   const cell_type&  cl,
                   const vector_type& u_TF)
    {
        auto id = cl.get_id();
        matrix_type  grec_op = m_gradrec_opers.at(id);
        auto zero_mindeg = 0;
        auto zero_maxdeg = m_degree;

        auto col_range      = cell_basis.range(zero_mindeg, zero_maxdeg);
        auto row_range      = disk::dof_range( 0,mesh_type::dimension);

        auto c_phi   = cell_basis.eval_functions(msh, cl, pt);
        matrix_type vec_phi = make_vectorial_matrix(c_phi, mesh_type::dimension);
        vector_type G_uT    = vec_phi * grec_op * u_TF;

        return G_uT;
    }

    vector_type
    grad_rec_uTF(  const mesh_type  &  msh,
                   const cell_type  &  cl,
                   const vector_type&  u_TF)
    {
        // We keep the name as grad_rec_oper and not only rec_oper since we
        // are still working with only D P^(k+1), i.e, we havent specified
        // the mean value
        auto id = msh.lookup(cl);
        matrix_type  grec_oper = m_gradrec_opers.at(id);
        //std::cout << "grec_oper size"<< grec_oper.rows() <<" x "<< grec_oper.cols() << std::endl;
        //std::cout << "u_TF size"<< u_TF.size() << std::endl;

        return grec_oper * u_TF;
    }

    matrix_type
    stab_oper( const mesh_type&  msh,
                   const cell_type &  cl,
                   const face_type&  fc)
    {
        auto cl_id = msh.lookup(cl);
        auto pos = face_position(msh, cl, fc);
        auto fbs = face_basis.size();
        auto cbs = cell_basis.range(0, m_degree_cell).size();
        auto num_faces = number_of_faces(msh, cl);
        auto num_total_dofs = cbs + fbs * num_faces;

        matrix_type ret = m_stabilization_opers.at(cl_id);
        return ret.block( 0, pos * num_total_dofs, fbs, num_total_dofs);
    }


};// plasticity

}//Disk
