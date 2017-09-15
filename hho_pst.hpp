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
#include "../plasticity/qps_ranges.hpp"
#include "plasticity.hpp"
namespace disk {

template<typename MeshType>
class punctual_tensor;

template<typename T, typename Storage>
class punctual_tensor<disk::mesh<T,2,Storage>>
{
public:
    typedef dynamic_matrix<T>                           matrix_type;
    typedef disk::mesh<T,2,Storage>                     mesh_type;
    typedef typename mesh_type::cell                    cell_type;
    typedef typename mesh_type::face                    face_type;
    typedef typename mesh_type::point_type              point_type;
    typedef disk::quadrature<mesh_type, cell_type>      cell_quad_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quad_type;

    matrix_type cell_quad_pts;
    matrix_type face_quad_pts;
    matrix_type grid_points;
    matrix_type face_bars;
    matrix_type cell_bar;
    qpspace_ranges qsr;
    size_t m_quad_degree;
    size_t cell_qps_size;
    size_t face_qps_size;
    size_t grid_pts_size;

    //WK: try to do this in the constructor.(I had problems in tensor_zero_vector,
    // since it s asking for the input of the constructor)

    punctual_tensor(){m_quad_degree = 10000; };
    //Warning if m_quad-degree has not a value,  all others functions wont be ok
    punctual_tensor(const mesh_type& msh, const cell_type& cl, const size_t& quad_degree):m_quad_degree(quad_degree)
    {

        auto cq  = cell_quad_type(m_quad_degree);
        auto fq  = face_quad_type(m_quad_degree);
        auto fcs = faces(msh, cl);
        auto num_fcs =  number_of_faces(msh, cl);

        auto num_face_qps = fq.integrate(msh, fcs[0]).size();
        auto num_cell_qps = cq.integrate(msh, cl).size();

        qsr = qpspace_ranges(num_cell_qps, num_face_qps, num_fcs);

        cell_qps_size = qsr.cell_range().size();
        face_qps_size = qsr.all_faces_range().size();
        grid_pts_size = num_fcs;
    };

    size_t cell_qps_cols(){   return cell_qps_size;}
    size_t face_qps_cols(){   return face_qps_size;}
    size_t grid_pts_cols(){   return grid_pts_size;}
    size_t facebars_cols(){   return grid_pts_size;}
    size_t cols(){ return   cell_qps_size + face_qps_size + 2 * grid_pts_size + 1;}

    void
    Zero(const size_t rows)
    {
        if(m_quad_degree >= 100 )
            throw std::logic_error("m_quad_degree >100, it has not been defined or it is taking impossible values");

        //std::cout << "ZERO QUAD DEGREE :"<< quad_degree << std::endl;
        cell_quad_pts = matrix_type::Zero(rows, cell_qps_size);
        face_quad_pts = matrix_type::Zero(rows, face_qps_size);
        grid_points   = matrix_type::Zero(rows, grid_pts_size);
        face_bars = matrix_type::Zero(rows, grid_pts_size);
        cell_bar  = matrix_type::Zero(rows, 1);
        return;
    }

    matrix_type cell_qp_matrix() { return cell_quad_pts; }
    matrix_type face_qp_matrix() { return face_quad_pts; }
    matrix_type grid_pt_matrix() { return grid_points; }
    matrix_type facebar_matrix() { return face_bars; }
    matrix_type cellbar_matrix() { return face_bars; }

    matrix_type
    join_all()
    {

        if(m_quad_degree >= 100 )
            throw std::logic_error("m_quad_degree >100, it has not been defined or it is taking impossible values");

        size_t cols = cell_quad_pts.cols() + face_quad_pts.cols() + grid_points.cols()
                        + face_bars.cols() + 1;

        auto fqs_size = face_quad_pts.cols();
        auto offset_fqs      = cell_quad_pts.cols();
        auto offset_vertex   = cell_quad_pts.cols() + fqs_size;
        auto offset_face_bar = offset_vertex   + grid_points.cols();
        auto offset_cell_bar = offset_face_bar + grid_points.cols();

        assert(cell_quad_pts.rows() == face_quad_pts.rows());
        assert(cell_quad_pts.rows() == grid_points.rows());
        assert(face_bars.rows() == grid_points.rows());
        assert(face_bars.rows() == cell_bar.rows());

        auto rows = cell_quad_pts.rows();
        matrix_type  ret = matrix_type::Zero(rows, cols);
        ret.block(0, 0, rows, offset_fqs)        = cell_quad_pts;
        ret.block(0, offset_fqs, rows, fqs_size) = face_quad_pts;
        ret.block(0, offset_vertex,   rows, grid_points.cols()) = grid_points;
        ret.block(0, offset_face_bar, rows, grid_points.cols()) = face_bars;
        ret.block(0, offset_cell_bar, rows, 1)   = cell_bar;
        return ret;
    }
    void
    split_all(const matrix_type &mat, const size_t& rows)
    {
        if(m_quad_degree >= 100 )
            throw std::logic_error("m_quad_degree >10E5, it has not been defined or it is taking impossible values");

        Zero(rows);

        auto offset_fqs      = cell_qps_size;
        auto offset_vertex   = offset_fqs + face_qps_size;
        auto offset_face_bar = offset_vertex   + grid_pts_size;
        auto offset_cell_bar = offset_face_bar + grid_pts_size;

        cell_quad_pts = mat.block(0, 0, rows, offset_fqs);
        face_quad_pts = mat.block(0, offset_fqs, rows, face_qps_size);
        grid_points   = mat.block(0, offset_vertex,   rows, grid_pts_size);
        face_bars     = mat.block(0, offset_face_bar, rows, grid_pts_size);
        cell_bar      = mat.block(0, offset_cell_bar, rows, 1);

        assert(cell_quad_pts.rows() == face_quad_pts.rows());
        assert(cell_quad_pts.rows() == grid_points.rows());
        assert(face_bars.rows() == grid_points.rows());
        assert(face_bars.rows() == cell_bar.rows());

        //m_quad_degree = tsr.sigma.m_quad_degree;
        return;
    }
    friend punctual_tensor<mesh_type> operator-(const punctual_tensor< mesh_type>& tp1,
                                                const punctual_tensor< mesh_type>& tp2)
    {
        punctual_tensor<mesh_type> ret;
        if(tp1.m_quad_degree != tp2.m_quad_degree)
            throw std::logic_error("(quad_degree) is different for these tensors, so it is not possible to do (-)operation.");
        ret.m_quad_degree = tp1.m_quad_degree;
        ret.cell_quad_pts = tp1.cell_quad_pts - tp2.cell_quad_pts;
        ret.face_quad_pts = tp1.face_quad_pts - tp2.face_quad_pts;
        ret.face_bars = tp1.face_bars - tp2.face_bars;
        ret.cell_bar  = tp1.cell_bar  - tp2.cell_bar;
        ret.grid_points  = tp1.grid_points - tp2.grid_points;

        ret.cell_qps_size = ret.cell_quad_pts.cols();
        ret.face_qps_size = ret.face_quad_pts.cols();
        ret.grid_pts_size = ret.grid_points.cols();

        return ret;
    }
    friend punctual_tensor<mesh_type> operator+(const punctual_tensor<mesh_type>& tp1,
                                               const punctual_tensor<mesh_type>& tp2)
    {
        if(tp1.m_quad_degree != tp2.m_quad_degree)
            throw std::logic_error("(quad_degree) is different for these tensors, so it is not possible to do (+)operation.");

        punctual_tensor<mesh_type> ret;
        ret.m_quad_degree = tp1.m_quad_degree;

        ret.cell_quad_pts = tp1.cell_quad_pts + tp2.cell_quad_pts;
        ret.face_quad_pts = tp1.face_quad_pts + tp2.face_quad_pts;
        ret.face_bars = tp1.face_bars + tp2.face_bars;
        ret.cell_bar  = tp1.cell_bar  + tp2.cell_bar;
        ret.grid_points  = tp1.grid_points  + tp2.grid_points;

        ret.cell_qps_size = ret.cell_quad_pts.cols();
        ret.face_qps_size = ret.face_quad_pts.cols();
        ret.grid_pts_size = ret.grid_points.cols();
        return ret;
    }
    friend punctual_tensor<mesh_type> operator*(const T& scalar, const punctual_tensor<mesh_type>& tp)
    {
        punctual_tensor<mesh_type> ret;
        ret.m_quad_degree = tp.m_quad_degree;
        ret.cell_quad_pts = scalar * tp.cell_quad_pts;
        ret.face_quad_pts = scalar * tp.face_quad_pts;
        ret.face_bars     = scalar * tp.face_bars;
        ret.cell_bar      = scalar * tp.cell_bar;
        ret.grid_points      = scalar * tp.grid_points;

        ret.cell_qps_size = ret.cell_quad_pts.cols();
        ret.face_qps_size = ret.face_quad_pts.cols();
        ret.grid_pts_size = ret.grid_points.cols();

        return ret;
    }

    friend std::ostream& operator<<(std::ostream& os, const punctual_tensor<mesh_type>& tp){
        os << "* m_quad_degree : "<< tp.m_quad_degree<<std::endl;
        os << "* cell_quad_pts : size = "<< tp.cell_quad_pts.rows()<< "x";
        os <<  tp.cell_qps_cols()<< std::endl;
        for(size_t j = 0 ; j < tp.cell_quad_pts.rows(); j++)
        {
            for(size_t i = 0 ; i < tp.cell_qps_cols(); i++)
                    os << tp.cell_quad_pts(i,j) << "  ";
            os <<std::endl;
        }
        os << "* face_quad_pts : size = "<< tp.face_quad_pts.rows();
        os << "x"<<  tp.face_qps_cols()<< std::endl;
        for(size_t j = 0 ; j < tp.face_quad_pts.rows(); j++)
        {
            for(size_t i = 0 ; i < tp.face_qps_cols(); i++)
                    os << tp.face_quad_pts(i,j) << "  ";
            os <<std::endl;
        }

        os << "* face_bar : size = "<< tp.face_bars.rows();
        os << "x"<<  tp.facebars_cols()<< std::endl;
        for(size_t j = 0 ; j < tp.face_quad_pts.rows(); j++)
        {
            for(size_t i = 0 ; i < tp.facebars_cols(); i++)
                    os << tp.face_quad_pts(i,j) << "  ";
            os <<std::endl;
        }

        os << "* vertices : size = "<< tp.grid_points.rows();
        os << "x"<<  tp.grid_pts_cols()<< std::endl;
        for(size_t j = 0 ; j < tp.grid_points.rows(); j++)
        {
            for(size_t i = 0 ; i < tp.grid_pts_cols(); i++)
                    os << tp.grid_points(i,j) << "  ";
            os <<std::endl;
        }
        os << "* cell_bar : size = "<< tp.cell_bar.rows();
        os << "x"<<  tp.cell_bar.cols()<< std::endl;
        for(size_t j = 0 ; j < tp.cell_bar.rows(); j++)
        {
            for(size_t i = 0 ; i < tp.cell_bar.cols(); i++)
                    os << tp.cell_bar(i,j) << "  ";
            os <<std::endl;
        }

        return os;
    }

    std::vector<point_type>
    tensor_points(const mesh_type& msh,
                  const cell_type& cell)
    {
        auto cq  = cell_quad_type(m_quad_degree);
        auto fq  = face_quad_type(m_quad_degree);
        auto cqs = cq.integrate(msh, cell);
        auto fcs = faces(msh, cell);
        auto pts = points(msh, cell);
        auto total_size = cell_qps_size + face_qps_size + 2 * grid_pts_size + 1;
        auto ret = std::vector<point_type>(total_size);

        size_t cont = 0;
        for(auto& cq : cqs)
            ret.at(cont++) = cq.point();

        for(auto& fc :fcs)
        {
            auto fqs = fq.integrate(msh, fc);
            for(auto& fq : fqs)
                ret.at(cont++) = fq.point();
        }

        for(auto& p : pts)
            ret.at(cont++) = p;

        for(auto& fc :fcs)
            ret.at(cont++) = barycenter(msh, fc);

        ret.at(cont) = barycenter(msh, cell);

        return ret;
    }
};

template<typename MeshType>
struct tensors
{
public:
    typedef typename MeshType::scalar_type      scalar_type;
    typedef typename MeshType::cell             cell_type;
    typedef Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic> tensor_matrix;

    punctual_tensor<MeshType>   sigma;
    punctual_tensor<MeshType>   gamma;
    punctual_tensor<MeshType>   xi_norm;

    tensors() //quad_degree(2*degree+2)
    {}

    void
    set_quad_degree(const size_t& degree)
    {}

    void
    zero_tensor(const MeshType& msh, const cell_type& cl, const size_t& degree)
    {
        sigma   = punctual_tensor<MeshType>(msh, cl, 2 * degree + 2);
        gamma   = punctual_tensor<MeshType>(msh, cl, 2 * degree + 2);
        xi_norm = punctual_tensor<MeshType>(msh, cl, 2 * degree + 2);

        sigma.Zero( 2);
        gamma.Zero( 2);
        xi_norm.Zero(1);

        #if 0
        std::cout << " SIZE SIGMA POINTS ("<< cl.get_id()<<") : "<< cols << std::endl;

        size_t cont = 0;
        auto cqs = cq.integrate(msh, cl);
        for(auto & cq: cqs)
        {
            auto b =  cq.point();
            std::cout<< "strName = strtrim(cellstr(num2str("<< cont <<",'(%d)')));"<<std::endl;
            std::cout<< "text("<< b.x()<< ","<< b.y() <<",strName,'VerticalAlignment','bottom');"<<std::endl;
            std::cout<< "plot("<< b.x()<< ","<< b.y() << ",'*')" << std::endl;
            cont++;
        }
        for(size_t ifc = 0;  ifc < fcs.size(); ifc++)
        {
            auto fqs  = fq.integrate(msh, fcs[ifc]); //WK: This should change if k different for each face
            auto face = *std::next(msh.faces_begin(), ifc);
            auto pts  = points(msh, face);

            std::cout << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
            std::cout << pts[0].y() << " " << pts[1].y() << "], 'Color', 'r');";
            std::cout << std::endl;


            for(auto & fpq: fqs)
            {
                auto b =  fpq.point();
                std::cout << "strName = strtrim(cellstr(num2str("<< cont <<",'(%d)')));"<<std::endl;
                std::cout << "text("<< b.x()<< ","<< b.y() <<",strName,'VerticalAlignment','bottom');"<<std::endl;
                std::cout << "plot("<< b.x()<< ","<< b.y() << ",'*')" << std::endl;
                cont++;
            }
        }
        #endif

    }
};

template<typename T, size_t DIM, typename Storage>
std::vector<tensors< mesh<T,DIM,Storage> >>
tensor_zero_vector(const mesh<T,DIM,Storage>& msh, const size_t& degree)
{
    typedef   mesh<T,DIM,Storage> mesh_type;
    // WK: degree of quadrature assumed as for (Dru,Dru)
    auto    vec = std::vector<tensors<mesh_type>>(msh.cells_size());

    for (auto& cl : msh)
    {
        auto  id = cl.get_id();
        tensors<mesh_type>  tsr;
        tsr.zero_tensor(msh, cl, degree);
        vec.at(id)   =  tsr;
    }
    return vec;
};
template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class gradient_reconstruction_pst
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

    size_t                                      m_degree;

public:

    gradient_reconstruction_pst()
    {}

    gradient_reconstruction_pst(const size_t& degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }
    gradient_reconstruction_pst(const size_t& degree, const size_t& quad_degree)
        : m_degree(degree)
    {
        std::cout << "degree : "<< degree<< "  ; m_degree :"<< m_degree << std::endl;
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(quad_degree);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }
    std::pair<matrix_type, matrix_type>
    compute(const mesh_type& msh, const cell_type& cl)
    {
        matrix_type stiff_mat = matrix_type::Zero(cell_basis.size(), cell_basis.size());

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point());
            for (size_t i = 0; i < cell_basis.size(); i++)
                for (size_t j = 0; j < cell_basis.size(); j++)
                    stiff_mat(i,j) += qp.weight() * mm_prod(c_dphi[i], c_dphi[j]);
        }

        /* LHS: take basis functions derivatives from degree 1 to K+1 */
        auto MG_rowcol_range = cell_basis.range(1, m_degree+1);
        matrix_type MG = take(stiff_mat, MG_rowcol_range, MG_rowcol_range);

        /* RHS, volumetric part. */
        auto BG_row_range = cell_basis.range(1, m_degree+1);
        auto BG_col_range = cell_basis.range(0, m_degree);

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_basis.range(0, m_degree).size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        matrix_type BG = matrix_type::Zero(BG_row_range.size(), dsr.total_size());

        BG.block(0, 0, BG_row_range.size(), BG_col_range.size()) =
                                take(stiff_mat, BG_row_range, BG_col_range);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs[face_i];
            auto n = normal(msh, cl, fc);
            auto face_quadpoints = face_quadrature.integrate(msh, fc);

            //matrix_type m1 = matrix_type::Zeros(BG_row_range.size(), BG_col_range.size());
            for (auto& qp : face_quadpoints)
            {
                auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
                auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point());
                for (size_t i = BG_row_range.min(), ii = 0; i < BG_row_range.max(); i++, ii++)
                {
                    for (size_t j = BG_col_range.min(), jj = 0; j < BG_col_range.max(); j++, jj++)
                    {
                        auto p1 = mm_prod(c_dphi[i], n);
                        auto p2 = mm_prod(c_phi[j], p1);
                        BG(ii,jj) -= qp.weight() * p2;
                    }
                }

                auto f_phi = face_basis.eval_functions(msh, fc, qp.point());

                for (size_t i = BG_row_range.min(), ii = 0; i < BG_row_range.max(); i++, ii++)
                {
                    for (size_t j = 0, jj = current_face_range.min();
                                j < current_face_range.size();
                                j++, jj++)
                    {
                        auto p1 = mm_prod(c_dphi[i], n);
                        auto p2 = mm_prod(f_phi[j], p1);
                        BG(ii,jj) += qp.weight() * p2;
                    }
                }
            }
        }

        matrix_type oper  = MG.llt().solve(BG);    // GT
        matrix_type data  = BG.transpose() * oper;  // A
        return std::make_pair(oper,data);
    }
};

template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class diffusion_like_stabilization_pst
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

    size_t                                  m_degree;

public:


    diffusion_like_stabilization_pst()
        : m_degree(1)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    diffusion_like_stabilization_pst(size_t degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    matrix_type
    compute(const mesh_type& msh, const cell_type& cl, const matrix_type& gradrec_oper)
    {
        auto pts = points(msh,cl);

        matrix_type mass_mat = matrix_type::Zero(cell_basis.size(), cell_basis.size());

        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
        for (auto& qp : cell_quadpoints)
        {
            auto c_phi = cell_basis.eval_functions(msh, cl, qp.point());
            for (size_t i = 0; i < cell_basis.size(); i++)
                for (size_t j = 0; j < cell_basis.size(); j++)
                    mass_mat(i,j) += qp.weight() * mm_prod(c_phi[i], c_phi[j]);
        }

        auto zero_range         = cell_basis.range(0, m_degree);
        auto one_range          = cell_basis.range(1, m_degree+1);

        // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

        //Step 1: compute \pi_T^k p_T^k v (third term).
        matrix_type M1 = take(mass_mat, zero_range, zero_range);
        matrix_type M2 = take(mass_mat, zero_range, one_range);
        matrix_type proj1 = -M1.llt().solve(M2*gradrec_oper);

        //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
        matrix_type I_T = matrix_type::Identity(zero_range.size(), zero_range.size());
        proj1.block(0, 0, zero_range.size(), zero_range.size()) += I_T;

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_basis.range( 0, m_degree).size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        matrix_type data = matrix_type::Zero(dsr.total_size(), dsr.total_size());

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

                for (size_t i = 0; i < face_basis.size(); i++)
                    for (size_t j = 0; j < face_basis.size(); j++)
                        face_mass_matrix(i,j)  += qp.weight() * mm_prod(f_phi[i], f_phi[j]);

                for (size_t i = 0; i < face_basis.size(); i++)
                    for (size_t j = 0; j < cell_basis.size(); j++)
                        face_trace_matrix(i,j) += qp.weight() * mm_prod(f_phi[i], c_phi[j]);
            }

            Eigen::LLT<matrix_type> piKF;
            piKF.compute(face_mass_matrix);

            // Step 3a: \pi_F^k( v_F - p_T^k v )
            auto face_range = current_face_range.remove_offset();
            matrix_type MR1 = take(face_trace_matrix, face_range, one_range);
            matrix_type proj2 = piKF.solve(MR1*gradrec_oper);
            matrix_type I_F = matrix_type::Identity(fbs, fbs);
            auto block_offset = current_face_range.min();
            proj2.block(0, block_offset, fbs, fbs) -= I_F;

            // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
            matrix_type MR2 = take(face_trace_matrix, face_range, zero_range);
            matrix_type proj3 = piKF.solve(MR2*proj1);

            matrix_type BRF = proj2 + proj3;

            data += BRF.transpose() * face_mass_matrix * BRF / h;
        }
        return data;
    }
};

template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class diffusion_like_static_condensation_pst
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

    size_t                                  m_degree;

public:
    diffusion_like_static_condensation_pst()
        : m_degree(1)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(2*m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }

    diffusion_like_static_condensation_pst(size_t degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(2*m_degree);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
        face_basis          = face_basis_type(m_degree);
    }


    matrix_type
    compute_lhs(const mesh_type& msh, const cell_type& cl,
                         const matrix_type& local_mat)
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
        matrix_type AC = K_FF - K_FT * AL;

        return AC;
    }

    #if 0
    vector_type
    plasticity_rhs(const mesh_type& msh,
                    );
    {

    }
    #endif
    vector_type
    compute_rhs(const mesh_type& msh, const cell_type& cl,
            const matrix_type& local_mat,
            const vector_type& cell_rhs,
            const vector_type& plastic_rhs)
    {
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        auto num_cell_dofs = cell_basis.size();
        auto num_face_dofs = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        size_t cell_size = dsr.cell_range().size();
        size_t face_size = dsr.all_faces_range().size();

        matrix_type K_TT = local_mat.topLeftCorner(cell_size, cell_size);
        matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
        vector_type bP_F = plastic_rhs.tail(face_size);
        vector_type bP_T = plastic_rhs.head(cell_size);


        assert(K_TT.cols() == cell_size);
        assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
        assert(bP_T.rows() + bP_F.rows() == plastic_rhs.rows());

        auto K_TT_ldlt = K_TT.llt();
        vector_type bL = K_TT_ldlt.solve(cell_rhs - bP_T);
        vector_type bC = /* no projection on faces, eqn. 26*/ - bP_F  - K_FT * bL;

        return bC;
    }
    vector_type
    recover(const mesh_type& msh, const cell_type& cl,
            const matrix_type& local_mat,
            const vector_type& cell_rhs,
            const vector_type& solF,
            const vector_type& plastic_rhs)
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

        vector_type bP_T = plastic_rhs.head(cell_size);

        vector_type solT = K_TT.llt().solve(cell_rhs - bP_T - K_TF*solF);

        assert(cell_size + all_faces_size == dsr.total_size() );
        ret.head(cell_size)         = solT;
        ret.tail(all_faces_size)    = solF;

        return ret;
    }
};

template<typename Mesh, //typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class assembler_pst
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

    //assembler_pst(){}               = delete;
    //assembler_pst(const assembler_pst&) = delete;
    //assembler_pst(assembler_pst&&)      = delete;
    assembler_pst()
    {}
    assembler_pst(const mesh_type& msh, const size_t& degree)
        : m_degree(degree)
    {

        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);

        m_num_unknowns = face_basis.size() * (msh.faces_size() + msh.boundary_faces_size());
    }
    void
    initialize_rhs(void)
    {
        rhs = vector_type::Zero(m_num_unknowns);
    }
    void
    initialize_lhs(void)
    {
        matrix = sparse_matrix_type(m_num_unknowns, m_num_unknowns);
    }

    template<typename LocalContrib>
    void
    assemble_lhs(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc)
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
                l2g.at(pos + i) = face_offset + i;
        }

        assert(lc.rows() == lc.cols());
        assert(lc.rows() == l2g.size());

        for (size_t i = 0; i < lc.rows(); i++)
        {
            for (size_t j = 0; j < lc.cols(); j++)
                m_triplets.push_back( triplet_type( l2g.at(i), l2g.at(j), lc(i,j) ) );
        }
    }
    template<typename LocalContrib>
    void
    assemble_rhs(const mesh_type& msh, const cell_type& cl, const LocalContrib& lc)
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
            auto pos   = face_i * face_basis.size();

            for (size_t i = 0; i < face_basis.size(); i++)
                l2g.at(pos + i) = face_offset + i;
        }

        assert(lc.size() == l2g.size());

        for (size_t i = 0; i < lc.size(); i++)
            rhs(l2g.at(i)) += lc(i);
    }

    void
    impose_boundary_conditions_lhs( mesh_type& msh)
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

            matrix_type MFF     = matrix_type::Zero(fbs, fbs);
            vector_type rhs_f   = vector_type::Zero(fbs);

            for (auto& qp : fqd)
            {
                auto f_phi = face_basis.eval_functions(msh, bfc, qp.point());

                for (size_t i = 0; i < f_phi.size(); i++)
                {
                    for (size_t j = 0; j < f_phi.size(); j++)
                        MFF(i,j) += qp.weight() * mm_prod(f_phi[i], f_phi[j]);
                }
            }

            for (size_t i = 0; i < MFF.rows(); i++)
            {
                for (size_t j = 0; j < MFF.cols(); j++)
                {
                    m_triplets.push_back( triplet_type(face_offset + i, face_offset_lagrange + j, MFF(i,j)) );
                    m_triplets.push_back( triplet_type(face_offset_lagrange + j, face_offset + i, MFF(i,j)) );
                }
            }

            face_i++;
        }
    }
    template<typename Function>
    void
    impose_boundary_conditions_rhs(mesh_type& msh, const Function& bc)
    {
        //std::cout << "Imposing_boundary_conditions_rhs" << std::endl;
        size_t fbs = face_basis.size();
        size_t face_i = 0;
        for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
        {
            auto bfc = *itor;

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id     = eid.second;
            auto face_offset = face_id * fbs;
            auto face_offset_lagrange = (msh.faces_size() + face_i) * fbs;
            auto fqd   = face_quadrature.integrate(msh, bfc);
            auto binfo = msh.boundary_info(bfc);

            vector_type rhs_f   = vector_type::Zero(fbs);
            #if 0
            auto pts = points(msh, bfc);
            auto b  = barycenter(msh,bfc);
            std::cout << "line([" << pts[0].x() << " " << pts[1].x() << "], [";
            std::cout << pts[0].y() << " " << pts[1].y() << "], 'Color', 'r');";
            std::cout << std::endl;
            std::cout << "strName = strtrim(cellstr(num2str("<< binfo.boundary_id <<",'(%d)')));"<<std::endl;
            std::cout << "text("<< b.x()<< ","<<b.y() <<",'strName','VerticalAlignment','bottom');"<<std::endl;
            #endif

            for (auto& qp : fqd)
            {
                auto f_phi = face_basis.eval_functions(msh, bfc, qp.point());

                for (size_t i = 0; i < f_phi.size(); i++)
                    rhs_f(i) += qp.weight() * mm_prod(f_phi[i], bc(qp.point(),binfo.boundary_id));
                    //rhs_f(i) += qp.weight() * mm_prod(f_phi[i], bc(qp.point()));

            }

            for (size_t i = 0; i < fbs; i++)
                rhs(face_offset_lagrange + i) = rhs_f(i);

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

template<typename Mesh, typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class projector_pst
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

    size_t                                      m_degree;

public:

    projector_pst()
        : m_degree(1)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(2*m_degree+ 2);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree+ 2);
    }

    projector_pst(size_t degree)
        : m_degree(degree)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(2*m_degree + 2);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree + 2);
    }

    matrix_type cell_mm;

    template<typename Function>
    vector_type
    compute_cell(const mesh_type& msh, const cell_type& cl, const Function& f)
    {

        auto x0 = 0.5;
        auto pts = points(msh, cl);
        size_t number = 0;
        for(auto& p : pts)
        {
            if(p.x() < x0)
            number = 1;
            if(p.x() > x0)
            number = 2;
        }
        if(number == 0)
            throw std::invalid_argument("Invalid number domain.");



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

#ifdef FILL_COLMAJOR
            for (size_t j = 0; j < phi.size(); j++)
            {
                for (size_t i = 0; i < phi.size(); i++)
                    mm(i,j) += qp.weight() * mm_prod(phi[i], phi[j]);

                rhs(j) += qp.weight() * mm_prod(fval, phi[j]);
            }
#else
            for (size_t i = 0; i < phi.size(); i++)
            {
                for (size_t j = 0; j < phi.size(); j++)
                    mm(i,j) += qp.weight() * mm_prod(phi[i], phi[j]);

                rhs(i) += qp.weight() * mm_prod(fval, phi[i]);
            }
#endif
#endif
        }

        cell_mm = mm;
        return mm.llt().solve(rhs);
    }

    template<typename Function>
    vector_type
    compute_whole(const mesh_type& msh, const cell_type& cl, const Function& f)
    {
        auto x0 = 0.5;
        auto pts = points(msh, cl);
        size_t number = 0;
        for(auto& p : pts)
        {
            if(p.x() < x0)
            number = 1;
            if(p.x() > x0)
            number = 2;
        }
        if(number == 0)
            throw std::invalid_argument("Invalid number domain.");

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
                auto fval = f(qp.point(), number);

#ifdef FILL_COLMAJOR
                for (size_t j = 0; j < phi.size(); j++)
                {
                    for (size_t i = 0; i < phi.size(); i++)
                        mm(i,j) += qp.weight() * mm_prod(phi[i], phi[j]);

                    rhs(j) += qp.weight() * mm_prod(fval, phi[j]);
                }
#else
                for (size_t i = 0; i < phi.size(); i++)
                {
                    for (size_t j = 0; j < phi.size(); j++)
                        mm(i,j) += qp.weight() * mm_prod(phi[i], phi[j]);

                    rhs(i) += qp.weight() * mm_prod(fval, phi[i]);
                }
#endif
            }

            ret.block(face_offset, 0, face_basis.size(), 1) = mm.llt().solve(rhs);
            face_offset += face_basis.size();
        }

        return ret;
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

    cell_basis_type                             cell_basis(degree+1);
    auto cell_quadrature     = cell_quadrature_type(2*(degree+1));
    auto cbs = cell_basis.size();

    matrix_type mm = matrix_type::Zero(cbs, cbs);
    auto cell_quadpoints = cell_quadrature.integrate(msh, cl);

    for (auto& qp : cell_quadpoints)
    {
        auto phi = cell_basis.eval_functions(msh, cl, qp.point());
        for (size_t i = 0; i < cbs; i++)
            for (size_t j = 0; j < cbs; j++)
                mm(i,j) += qp.weight() * mm_prod(phi[i], phi[j]);
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
template<typename T, int DIM>
dynamic_matrix<T>
make_gradient_matrix(const std::vector<static_vector<T,DIM>>& dphi)
{
    dynamic_matrix<T> ret(DIM, dphi.size());

    for (size_t i = 0; i < dphi.size(); i++)
        ret.col(i) = dphi[i];

    return ret;
}

template<typename T>
dynamic_matrix<T>
make_gradient_matrix(const std::vector<T>& dphi)
{
    dynamic_matrix<T> ret(1, dphi.size());

    for (size_t i = 0; i < dphi.size(); i++)
        ret(i) = dphi[i];

    return ret;
}
template<typename T, int DIM>
T
make_prod(const  dynamic_vector<T>& st, const static_vector<T,DIM>& n )
{
    T p1(0.);
    //for(size_t i = 0; i < mesh_type::dimension; i++)
    //    p1  += st(i) * n(i);
    p1 = st.dot(n);
    return p1;
}
template<typename T>
T
make_prod(const  dynamic_vector<T>& st, const T& n )
{
    std::cout << "******************** THIS IS ONLY FOR 1D ********************" << std::endl;
    std::cout << "******************** THIS IS ONLY FOR 1D ********************" << std::endl;
    std::cout << "******************** THIS IS ONLY FOR 1D ********************" << std::endl;

    return  st(0) * n;
}

template<typename T, typename MeshType, typename CellType, typename CellBasisType>
dynamic_vector<T>
nabla_ruh_at_point( const point<T,2>& pt, //Check for this to be generalized to other dimensions.
               const MeshType&  msh,
               const CellType&  cl,
               const CellBasisType& cell_basis,
               const dynamic_vector<T>& ruh,
               const size_t m_degree)
{
    auto col_range      = cell_basis.range(1,m_degree+1);
    auto row_range      = disk::dof_range(0,MeshType::dimension);

    auto dphi   = cell_basis.eval_gradients(msh, cl, pt);
    dynamic_matrix<T> dphi_matrix = make_gradient_matrix(dphi);
    dynamic_matrix<T> dphi_taken  = take(dphi_matrix, row_range, col_range);
    dynamic_vector<T> dphi_ruh    = dphi_taken * ruh;
    return dphi_ruh;
}

template<typename T,typename Mesh,typename CellBasisType, typename CellQuadType,
                       typename FaceBasisType, typename FaceQuadType>
class plasticity
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
    typedef static_vector<T, mesh_type::dimension>   m_static_vector;
    typedef typename CellBasisType::gradient_value_type gradient_value_type_pr;

    cell_basis_type                         cell_basis;
    cell_quadrature_type                    cell_quadrature;

    face_basis_type                         face_basis;
    face_quadrature_type                    face_quadrature;

    size_t                                  m_quad_degree;
    size_t                                  m_degree;
    T                                       m_method_coef;
public:
    //vector_type     rhs;

    //plasticity()                = delete;
    plasticity()
    {}
    plasticity(const size_t& degree, const size_t& quad_degree)

    //WK: for mu as tensor one should check this part
        : m_degree(degree), m_quad_degree(quad_degree)
    {

        cell_basis  =  cell_basis_type(m_degree + 1); // This is k + 1  just if I want to use \nabla r^{k + 1}
        face_basis  =  face_basis_type(m_degree);

        cell_quadrature     = cell_quadrature_type(m_quad_degree);
        face_quadrature     = face_quadrature_type(m_quad_degree);

    }

    std::pair< vector_type,T>
    eval_plasticity(const point<T,2> & pt,
                  const mesh_type    & msh,
                  const cell_type    & cl,
                  const vector_type  & ruh,
                  const vector_type  & uh_T,
                  const plasticity_data<T> & pst,
                  vector_type     & sigma)
    {
        auto dphi   = cell_basis.eval_gradients(msh, cl, pt);

        auto cell_range = cell_basis.range(0,m_degree);
        auto col_range  = cell_basis.range(1,m_degree + 1);
        auto row_range  = disk::dof_range(0,mesh_type::dimension);

        matrix_type dphi_matrix = make_gradient_matrix(dphi);
        matrix_type dphi_taken  = take(dphi_matrix, row_range, col_range);
        vector_type dphi_ruh    = dphi_taken * ruh;

        vector_type xi =  vector_type::Zero(mesh_type::dimension);
        if(pst.hho == 1)
             xi = sigma  +  pst.alpha * dphi_ruh;
        if(pst.hho == 2)
        {
            matrix_type dphi_  = take(dphi_matrix, row_range, cell_range);
            xi = sigma  +  pst.alpha * dphi_* uh_T;
        }

        auto xi_norm = xi.norm();

        vector_type gamma = vector_type::Zero(mesh_type::dimension);
        if(xi_norm > pst.yield)
            gamma = m_method_coef * (xi_norm - pst.yield) * (xi / xi_norm);

        sigma += pst.alpha * ( dphi_ruh - gamma );

        return std::make_pair(gamma, xi_norm);
    }

    void
    compute_other(const mesh_type& msh,
                 const cell_type& cl,
                 const vector_type& ruh,
                 const vector_type  & uh_T,
                 const plasticity_data<T>& pst,
                 tensors<mesh_type> & tsr)
    {
        auto fcs = faces(msh, cl);
        auto pts = points(msh, cl);
        // Values at the vertices
        for(size_t i = 0;  i < fcs.size(); i++)
        {
            vector_type v_sigma = tsr.sigma.grid_points.col(i);
            auto pair  =  eval_plasticity(pts.at(i), msh, cl, ruh, uh_T,pst, v_sigma);

            tsr.sigma.grid_points.col(i) = v_sigma;
            tsr.gamma.grid_points.col(i) = pair.first;
            tsr.xi_norm.grid_points(i)   = pair.second;
        }

        // Values at the barycenter of the faces
        for(size_t i = 0; i < pts.size(); i++)
        {
            auto face  = fcs.at(i);
            auto fbar  = barycenter(msh, face);

            vector_type fb_sigma = tsr.sigma.face_bars.col(i);
            auto pair  = eval_plasticity(fbar, msh, cl, ruh, uh_T,pst, fb_sigma);

            tsr.sigma.face_bars.col(i)  = fb_sigma;
            tsr.gamma.face_bars.col(i)  = pair.first;
            tsr.xi_norm.face_bars(i)    = pair.second;
        }

        // Values at the barycenter of the cell
        vector_type cb_sigma = tsr.sigma.cell_bar;
        auto cbar  =  barycenter(msh, cl);
        auto pair  =  eval_plasticity(cbar, msh, cl, ruh, uh_T, pst, cb_sigma);

        tsr.sigma.cell_bar   = cb_sigma;
        tsr.gamma.cell_bar   = pair.first;
        tsr.xi_norm.cell_bar(0) = pair.second;

        return;
    }

    vector_type
    compute(const mesh_type& msh,
                 const cell_type& cl,
                 const matrix_type& rec_oper,
                 const vector_type& uh_TF,
                 const plasticity_data<T>& pst,
                 tensors<mesh_type> & tsr)
    {
        if(pst.method == true)
            m_method_coef = 1. / pst.alpha ;
        else
            m_method_coef = 1. / (pst.alpha + pst.mu);

        //WK: for p-adaptation quad_degree is not the same for all cells
        auto fcs = faces(msh, cl);

        auto num_faces      = fcs.size();
        auto cell_range     = cell_basis.range(0,m_degree);
        auto col_range      = cell_basis.range(1,m_degree+1);
        auto row_range      = disk::dof_range(0,mesh_type::dimension);
        auto num_cell_dofs  = cell_range.size();
        auto num_face_dofs  = face_basis.size();
        dofspace_ranges    dsr(num_cell_dofs, num_face_dofs, num_faces);

        vector_type    uh_T = uh_TF.head(num_cell_dofs);
        vector_type    ruh  = rec_oper * uh_TF;
        vector_type    rhs  = vector_type::Zero(dsr.total_size());

        auto   cqs   =  cell_quadrature.integrate(msh, cl);
        size_t cont  =  0;

        for (auto& qp : cqs)
        {
            if(cont > cqs.size())
                throw std::logic_error("Wrong access to columns");
            vector_type sigma =  tsr.sigma.cell_quad_pts.col(cont);
            auto pair  =  eval_plasticity(qp.point(), msh, cl, ruh, uh_T, pst, sigma);
            vector_type gamma = pair.first;
            tsr.sigma.cell_quad_pts.col(cont)   =  sigma;
            tsr.gamma.cell_quad_pts.col(cont)   =  gamma;
            tsr.xi_norm.cell_quad_pts(cont) =  pair.second;

            //( Dv_T, sigma - alpha* gamma)_T
            auto dphi   = cell_basis.eval_gradients(msh, cl, qp.point());
            matrix_type   dphi_matrix = make_gradient_matrix(dphi);

            if(pst.hho == 2)
            {
                matrix_type dphi_  = take(dphi_matrix, row_range, cell_range);
                rhs.head(num_cell_dofs)  += qp.weight() * dphi_.transpose() *
                                                    (sigma - pst.alpha * gamma);
            }
            else if(pst.hho ==  1)
            {
                matrix_type dphi_taken  = take(dphi_matrix, row_range, col_range);
                matrix_type dphi_rec = dphi_taken * rec_oper;
                rhs  += qp.weight() * dphi_rec.transpose() * (sigma - pst.alpha * gamma);
            }
            else
                throw std::invalid_argument("Invalid discretization name for the plasticity term. Review the name in parametes.txt or pst.hho");

            ++cont;
        }

        // Boundary term
        size_t fcont = 0;
        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fc   = fcs.at(face_i);
            auto fqs  = face_quadrature.integrate(msh, fc);

            static_vector<T,2> n = normal(msh, cl, fc);

            for (auto& qp : fqs)
            {
                vector_type f_sigma =  tsr.sigma.face_quad_pts.col(fcont);
                auto pair  =  eval_plasticity(qp.point(), msh, cl, ruh, uh_T, pst, f_sigma);
                vector_type f_gamma  = pair.first;

                tsr.sigma.face_quad_pts.col(fcont)  = f_sigma;
                tsr.gamma.face_quad_pts.col(fcont)  = f_gamma;
                tsr.xi_norm.face_quad_pts(fcont)    = pair.second;

                assert( fcont < tsr.sigma.face_qps_cols() );
                assert( fcont < tsr.gamma.face_qps_cols() );
                assert( fcont < tsr.xi_norm.face_qps_cols() );

                auto f_phi  = face_basis.eval_functions(msh, fc, qp.point());
                auto c_phi  = cell_basis.eval_functions(msh, cl, qp.point());

                // ((sigma - alpha*gamma) * n
                vector_type sigma_agamma =  f_sigma - pst.alpha * f_gamma;
                scalar_type p1  = make_prod( sigma_agamma, n);
                if(pst.hho == 2)
                {
                    // ( (sigma - alpha*gamma) * n, vT )_F
                    for (size_t i = 0; i < num_cell_dofs; i++)
                        rhs(i)   -= qp.weight() * p1 * c_phi.at(i);

                    // ( (sigma - alpha*gamma) * n, vF )_F
                    for (size_t i = 0, ii = current_face_range.min();
                                        i < current_face_range.size(); i++, ii++)
                        rhs(ii)  += qp.weight() * p1 * f_phi.at(i);
                }
                ++fcont;
            }
        }

        compute_other(msh, cl, ruh, uh_T, pst, tsr);

        return rhs;
    }

    template<typename Solution>
    dynamic_vector<T>
    compute(const mesh_type& msh,
                 const cell_type& cl,
                 const matrix_type& rec_oper,
                 const vector_type& uh_TF,
                 const plasticity_data<T>& pst,
                 const Solution& solution )
    {
        //WK: for p-adaptation quad_degree is not the same for all cells
        auto fcs = faces(msh, cl);

        auto num_faces      = fcs.size();
        auto cell_range     = cell_basis.range(0,m_degree);
        auto col_range      = cell_basis.range(1,m_degree+1);
        auto row_range      = disk::dof_range(0,mesh_type::dimension);
        auto num_cell_dofs  = cell_range.size();
        auto num_face_dofs  = face_basis.size();
        dofspace_ranges    dsr(num_cell_dofs, num_face_dofs, num_faces);
        vector_type    uh_T = uh_TF.head(num_cell_dofs);
        vector_type    ruh  = rec_oper * uh_TF;
        vector_type    rhs  = vector_type::Zero(dsr.total_size());

        auto cqs  =  cell_quadrature.integrate(msh, cl);
        auto number = set_cell_number(msh, cl);
        size_t cont  =  0;
        for (auto& qp : cqs)
        {
            if(cont > cqs.size())
                throw std::logic_error("Wrong access to columns");

            //( Dv_T, sigma - alpha* gamma)_T
            auto dphi   = cell_basis.eval_gradients(msh, cl, qp.point());
            matrix_type   dphi_matrix = make_gradient_matrix(dphi);

            if(pst.hho == 2)
            {
                matrix_type dphi_  = take(dphi_matrix, row_range, cell_range);
                rhs.head(num_cell_dofs)  += qp.weight() * dphi_.transpose()
                                                * solution.st( qp.point(), number);
            }
            else if(pst.hho ==  1)
            {
                matrix_type dphi_taken  = take(dphi_matrix, row_range, col_range);
                matrix_type dphi_rec = dphi_taken * rec_oper;
                rhs  += qp.weight() * dphi_rec.transpose() * solution.st(qp.point(), number);
            }
            else
                throw std::invalid_argument("Invalid discretization name for the plasticity term. Review the name in parametes.txt or pst.hho");

        }

        // Boundary term
        if(pst.hho == 2)
        {
            for (size_t face_i = 0; face_i < num_faces; face_i++)
            {
                auto current_face_range = dsr.face_range(face_i);
                auto fc   = fcs.at(face_i);
                auto fqs  = face_quadrature.integrate(msh, fc);

                static_vector<T,2> n = normal(msh, cl, fc);


                for (auto& qp : fqs)
                {

                    auto f_phi  = face_basis.eval_functions(msh, fc, qp.point());
                    auto c_phi  = cell_basis.eval_functions(msh, cl, qp.point());
                    size_t  number_f;
                    if(!msh.is_boundary(fc))
                    {
                        auto cell_id =  cl.get_id();
                        auto cells   =  face_owner_cells_ids(msh, fc);
                        auto neighbor_id =  cells.first;
                        if(cells.first != cell_id && cells.second == cell_id) //second part just to confirm cell is also owner
                            neighbor_id =  cells.first;
                        else if(cells.second != cell_id && cells.first == cell_id) //second part just to confirm cell is also owner
                            neighbor_id = cells.second;
                        else
                            throw std::logic_error("Neighbors not found");

                        auto ngh =  *std::next(msh.cells_begin(), neighbor_id);
                        auto number_nc = set_cell_number(msh, ngh);

                        //Esto escoge el valor de la derecha para tau_F
                        number_f = number;
                        if(number != number_nc )
                            number_f == 1;
                    }
                    else
                        number_f = number;


                    // ((sigma - alpha*gamma) * n
                    vector_type str   =  solution.st( qp.point(), number);
                    vector_type str_f =  solution.st( qp.point(), number_f);
                    scalar_type p1_c  =  make_prod( str , n);
                    scalar_type p1_f  =  make_prod( str_f , n);

                    // ( (sigma - alpha*gamma) * n, vT )_F
                    for (size_t i = 0; i < num_cell_dofs; i++)
                        rhs(i)   -= qp.weight() * p1_c * c_phi.at(i);

                    // ( (sigma - alpha*gamma) * n, vF )_F
                    for (size_t i = 0, ii = current_face_range.min();
                                        i < current_face_range.size(); i++, ii++)
                        rhs(ii)  += qp.weight() * p1_f * f_phi.at(i);
                }
            }
        }

        //compute_other(msh, cl, ruh, pst, tsr);
        return rhs;
    }
};

template< typename T, typename Storage, typename CellType>
size_t
set_cell_number(const  disk::mesh<T,2,Storage>& msh, const CellType& cl)
{
    //std::cout << "CELL :"<<cl.get_id() << std::endl;
    auto x0 = 0.5;
    auto pts = points(msh, cl);
    size_t number = 0;
    for(auto& p : pts)
    {
        if(p.x() < x0)
            number = 1;
        if(p.x() > x0)
            number = 2;
        //std::cout << " * point : "<< p.x() << "  "<< p.y()  << std::endl;
        //std::cout << " * number: "<< number << std::endl;
    }
    if(number == 0)
        throw std::invalid_argument("Invalid number domain.");
    //std::cout << " * number_final: "<< number << std::endl;
    return number;
};


template<typename T, size_t DIM>
struct solution
{
    typedef point<T, DIM>               point_type;
    typedef Eigen::Matrix< T, DIM, 1>   gradient_vector;

    typedef std::function<T (const point_type & p)>    function;
    typedef std::function<gradient_vector(const point_type & p)>   dfunction;
};

template<typename T, size_t DIM>
struct poiseuille;

template<typename T>
struct poiseuille<T,1>: public solution<T,1>
{
    typedef typename solution<T,1>::function     function;
    typedef typename solution<T,1>::dfunction    dfunction;
    typedef typename solution<T,1>::point_type   point_type;
    typedef typename solution<T,1>::gradient_vector gradient_vector;

    bool            is_exact;
    function f,sf;
    dfunction df;
    std::string     name;

    poiseuille(plasticity_data<T>& pst)
    {
        is_exact = true;
        name = "pouseille";
        T fvalue(1.);
        pst.f = fvalue;

        f  =  [fvalue](const point_type& p) -> T {
            T ret = fvalue;
            return ret;
         };
        sf =  [pst](const point_type& p)-> T
        {
            T  ret, xO(0), xL(pst.Lref);
            T  xc = pst.Lref/2.0 - pst.Bn*pst.Lref*pst.f;

            if( p.x() >= xO & p.x() <=  xc)
                ret = (xc*xc - (p.x() - xc)*(p.x() - xc));
            if( p.x() > xc & p.x() <  (xL- xc))
                ret = xc*xc;
            if( p.x() <= xL & p.x() >= (xL - xc))
                ret = (xc*xc - (p.x() - (pst.Lref - xc))*(p.x() - (pst.Lref -xc)));
            return ret*(0.5 * pst.f / pst.mu);
        };

        df =  [pst](const point_type& p)-> gradient_vector
        {
            gradient_vector ret   = gradient_vector::Zero();
            T  xO(0);
            T  xL(pst.Lref);
            T  xc = pst.Lref / 2.0  -  pst.Bn * pst.Lref * pst.f;

            if( p.x() >= xO & p.x() <=  xc)
                ret(0) = - 2.*(p.x() - xc);
            if( p.x() > xc & p.x() <  (xL- xc))
                ret(0) = 0.;
            if( p.x() <= xL & p.x() >= (xL - xc))
                ret(0) = - 2.*(p.x() - (pst.Lref - xc));
            ret = (0.5*pst.f/pst.mu)*ret;
            return ret;
        };
    };
};
template<typename T, size_t DIM>
struct tuyau;

template<typename T>
struct tuyau<T,2>: public solution<T,2>
{
    typedef typename solution<T,2>::function        function;
    typedef typename solution<T,2>::dfunction       dfunction;
    typedef typename solution<T,2>::point_type      point_type;
    typedef typename solution<T,2>::gradient_vector gradient_vector;


    typedef std::function<T (const point_type & p, const size_t n)> nfunction;
    typedef std::function<gradient_vector(const point_type & p,const size_t n)>
                                                                         ndfunction;
#if 0
    // not using the number
    bool            is_exact;
    function        f,sf;
    dfunction       df;
    std::string     name;

    tuyau(plasticity_data<T>& pst)
    {
        is_exact = false;

        //name = "square/31_EST_STR/PER01";
        //name = "Estimator/AEM/square";
        //name = "square/21_GAUSS_PTS/PER01";
        //name = "square/11_MARKER_4/PER01";
        //name = "square/12_MARKER_4/PER01/5734";
        name = "square/test/M3";
        pst.f = 1.;

        f  = [](const point_type& p) -> T {
            T ret = 1.;
            return ret;
        };

        sf = [](const point_type& p) -> T {
            T ret = 0.;
            return ret;
        };
        df = [](const point_type& p) -> gradient_vector {
            return gradient_vector::Zero();
        };
    };
    #endif
    bool            is_exact;
    function        f;
    nfunction       sf;
    ndfunction      df,st;
    std::string     name;
    std::string     identity;

    tuyau(plasticity_data<T>& pst, const std::string& name)
    {
        is_exact = false;
        pst.f = 1.;
        identity = "square";

        f = [&](const point_type& p) -> T {
            T ret = 1.;
            return ret;
        };

        sf = [&](const point_type& p, const size_t number) -> T {
            T ret = 0.;
            return ret;
        };

        df = [](const point_type& p,  const size_t number) -> gradient_vector {
            return gradient_vector::Zero();
        };

        st = [&](const point_type& p, const size_t number) -> gradient_vector {

            throw std::invalid_argument("st function is not intended for square problem.Review solution definition and parameters.txt");
            return gradient_vector::Zero();
        };
    };

};

template<typename T, size_t DIM>
struct circular_tuyau;

template<typename T>
struct circular_tuyau<T,2>: public solution<T,2>
{
    typedef typename solution<T,2>::function        function;
    typedef typename solution<T,2>::dfunction       dfunction;
    typedef typename solution<T,2>::point_type      point_type;
    typedef typename solution<T,2>::gradient_vector gradient_vector;

    typedef std::function<T (const point_type & p, const size_t n)> nfunction;
    typedef std::function<gradient_vector(const point_type & p,const size_t n)>
                                                                         ndfunction;

    #if 0
    bool            is_exact;
    function        f,sf;
    dfunction       df;
    std::string     name;
    std::string     identity;

    circular_tuyau(plasticity_data<T>& pst)
    {
        is_exact = true;
        identity = "circular";
        //name = "circular/31_EST_TOT";
        //name = "circular/22_GAUSS_PTS/PER01";
        //name = "circular/11_MARKER_4/PER01";
        //name = "circular/12_MARKER_4/PER01/5734";
        name  = "circular/ERROR/M1";
        pst.f = 1.;
        T   R = 1.; // this should be read from the mesh

        f  = [pst](const point_type& p) -> T {
            T ret = 1.;
            return ret;
        };
        sf = [R,pst](const point_type& p) -> T {
            T  ret = 0;
            T    r = std::sqrt(p.x()*p.x() + p.y()*p.y());
            T   Bi = 2. * pst.yield/(pst.f * R);

            if(r/R > Bi)
                ret = 0.5*(1. - (r/R)*(r/R)) - Bi*(1. - r/R);
            else
                ret = 0.5*(1. - Bi)*(1. - Bi);

            ret = (0.5*pst.f*R*R/pst.mu)*ret;
            return ret;
        };
        df = [R,pst](const point_type& p) -> gradient_vector {

            gradient_vector ret   = gradient_vector::Zero();

            T    r = std::sqrt(p.x()*p.x() + p.y()*p.y());
            T   Bi = 2.*pst.yield/(pst.f*R);

            if(r/R >= Bi)
            {
                T dfdr =  - r/(R*R) + Bi/R;
                T drdx =  p.x()/r;
                T drdy =  p.y()/r;
                ret(0) = dfdr*drdx;
                ret(1) = dfdr*drdy;
            }

            ret = (0.5*pst.f*R*R/pst.mu)*ret;
            return ret;
        };
    };


    circular_tuyau(plasticity_data<T>& pst, const std::string& name)
    {
        is_exact = true;
        identity = "circular";

        pst.f = 1.;
        T   R = 1.; // this should be read from the mesh

        f  = [pst](const point_type& p) -> T {
            T ret = 1.;
            return ret;
        };
        sf = [R,pst](const point_type& p) -> T {
            T  ret = 0;
            T    r = std::sqrt(p.x()*p.x() + p.y()*p.y());
            T   Bi = 2. * pst.yield/(pst.f * R);

            if(r/R > Bi)
                ret = 0.5*(1. - (r/R)*(r/R)) - Bi*(1. - r/R);
            else
                ret = 0.5*(1. - Bi)*(1. - Bi);

            ret = (0.5*pst.f*R*R/pst.mu)*ret;
            return ret;
        };
        df = [R,pst](const point_type& p) -> gradient_vector {

            gradient_vector ret   = gradient_vector::Zero();

            T    r = std::sqrt(p.x()*p.x() + p.y()*p.y());
            T   Bi = 2.*pst.yield/(pst.f*R);

            if(r/R >= Bi)
            {
                T dfdr =  - r/(R*R) + Bi/R;
                T drdx =  p.x()/r;
                T drdy =  p.y()/r;
                ret(0) = dfdr*drdx;
                ret(1) = dfdr*drdy;
            }

            ret = (0.5*pst.f*R*R/pst.mu)*ret;
            return ret;
        };
    };
    #endif

    bool            is_exact;
    function        f;
    nfunction       sf;
    ndfunction      df,st;
    std::string     name;
    std::string     identity;


    circular_tuyau(plasticity_data<T>& pst, const std::string& name)
    {
        is_exact = true;
        identity = "circular";
        pst.f = 1.;
        T   R = 1.; // this should be read from the mesh
        f  = [pst](const point_type& p) -> T {
            T ret = 1.;
            return ret;
        };
        sf = [R,pst](const point_type& p, const size_t number) -> T {
            T  ret = 0;
            T    r = std::sqrt(p.x()*p.x() + p.y()*p.y());
            T   Bi = 2. * pst.yield/(pst.f * R);

            if(r/R > Bi)
                ret = 0.5*(1. - (r/R)*(r/R)) - Bi*(1. - r/R);
            else
                ret = 0.5*(1. - Bi)*(1. - Bi);

            ret = (0.5*pst.f*R*R/pst.mu)*ret;
            return ret;
        };
        df = [R,pst](const point_type& p, const size_t number) -> gradient_vector {

            gradient_vector ret   = gradient_vector::Zero();

            T    r = std::sqrt(p.x()*p.x() + p.y()*p.y());
            T   Bi = 2.*pst.yield/(pst.f*R);

            if(r/R >= Bi)
            {
                T dfdr =  - r/(R*R) + Bi/R;
                T drdx =  p.x()/r;
                T drdy =  p.y()/r;
                ret(0) = dfdr*drdx;
                ret(1) = dfdr*drdy;
            }

            ret = (0.5*pst.f*R*R/pst.mu)*ret;
            return ret;
        };
        st = [&](const point_type& p, const size_t number) -> gradient_vector {

            throw std::invalid_argument("st function is not intended for circular_tuyau problem. Review solution definition and parameters.txt");
            return gradient_vector::Zero();
        };
    };

};

template<typename T, size_t DIM>
struct diffusion;

template<typename T>
struct diffusion<T,2>: public solution<T,2>
{
    typedef typename solution<T,2>::function        function;
    typedef typename solution<T,2>::dfunction       dfunction;
    typedef typename solution<T,2>::point_type      point_type;
    typedef typename solution<T,2>::gradient_vector gradient_vector;

    typedef std::function<T (const point_type & p, const size_t n)> nfunction;
    typedef std::function<gradient_vector(const point_type & p,const size_t n)>
                                                                     ndfunction;

    bool            is_exact;
    function        f,dst;
    nfunction       sf;
    ndfunction      df, st;
    std::string     name, identity;
    #if 0
    diffusion(plasticity_data<T>& pst)
    {
        is_exact = false;
        pst.f  = 1.;
        name = "Estimator/AEM/square";

        f  = [](const point_type& p) -> T {
            return  2.* M_PI * M_PI * sin(p.x() * M_PI) * sin(p.y() * M_PI);
        };
        sf = [](const point_type& p) -> T {
            return sin(p.x() * M_PI) * sin(p.y() * M_PI);
        };
        df = [](const point_type& p) -> gradient_vector {
            gradient_vector ret =  gradient_vector::Zero();
            ret(0) = M_PI * cos(p.x() * M_PI) * sin(p.y() * M_PI);
            ret(1) = M_PI * sin(p.x() * M_PI) * cos(p.y() * M_PI);
            return ret;
        };
        st = [](const point_type& p) -> gradient_vector {
            gradient_vector ret =  gradient_vector::Zero();
            ret(0) = -0.5 * M_PI * cos(p.x() * M_PI) * sin(p.y() * M_PI);
            ret(1) = -0.5 * M_PI * sin(p.x() * M_PI) * cos(p.y() * M_PI);
            return ret;
        };
        dst = [](const point_type& p) -> T {
            return  M_PI * M_PI * sin(p.x() * M_PI) * sin(p.y() * M_PI);
        };


    };
    #endif
    diffusion(plasticity_data<T>& pst , const std::string& name)
    {
        is_exact = true;
        pst.f    = 1.;
        identity = "diffusion";
        f  = [&](const point_type& p) -> T {
            T x0    = 0.5;
            T c1    = 1.;
            #if 0
            if( p.x() <= x0)
            {
                return 1.;
            }
            else
            {
                return 0.5;
            }
            #endif
            return 0.;
        };

        sf = [&](const point_type& p, const size_t number) -> T {
            T ret;
            T x0    = 0.5;

            switch (number)
            {
                case 1:
                    ret = -0.25 * p.x() * p.x() -  0.25 * p.y() * p.y();
                    //std::cout << " ("<< p.x() << " " <<p.y() << ") = " << ret << std::endl;
                    //std::cout << " bdn = 1" << std::endl;
                    break;
                case 2:
                    ret = -0.25 * x0 * x0 -  0.25 * p.y() * p.y();
                    //std::cout << " ("<< p.x() << " " <<p.y() << ") = " << ret << std::endl;
                    //std::cout << " bdn  = 2" << std::endl;
                    break;
                default:
                std::cout << "number: "<<number << std::endl;
                    throw std::invalid_argument("Not given values for this number, sf");
                    break;
            }
            return ret;
        };
        df = [&](const point_type& p, const size_t number) -> gradient_vector {
            gradient_vector ret =  gradient_vector::Zero();
            switch (number)
            {
                case 1:
                    ret(0) = -0.5 * p.x() ;
                    ret(1) = -0.5 * p.y() ;
                    break;
                case 2:
                    ret(0) =  0.;
                    ret(1) = -0.5 * p.y() ;
                    break;
                default:
                    throw std::invalid_argument("Not given values for this number, df");
                    break;
            }
            return ret;
        };
        st = [&](const point_type& p, const size_t number) -> gradient_vector {

            gradient_vector ret =  gradient_vector::Zero();
            // Esto funciona para el metodo 1, pero hace falta saber que hacer con el metodo dos.
            T c1    = 1.;
            switch (number)
            {
                case 1:
                    ret(0) = p.x();
                    ret(1) = p.x();//0.;
                    break;
                case 2:
                    ret(0) = 0.5 * p.x();
                    ret(1) = 0.5 * p.x();//c1;
                    break;
                default:
                    throw std::invalid_argument("Not given values for this number, st");
                    break;
            }
            return ret;// gradient_vector::Zero();
        };
        dst = [&](const point_type& p) -> T {
            T x0    = 0.5;
            T c1    = 1.;

            if( p.x() <= x0)
            {
                return 1.;
            }
            else
            {
                return 0.5;
            }
        };

    };

};


}//disk
