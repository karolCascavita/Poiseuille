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
struct plasticity_data
{
    plasticity_data(): Lref(1.), Vref(1.), Bn(0.), mu(1.), alpha(1.), beta(1.), f(1.),
                        method(false), hho(2)
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
    size_t hho;

    friend std::ostream& operator<<(std::ostream& os, const plasticity_data<T>& p){
        os << "Plastic data: "<<std::endl;
        os << "* method : "<< p.method<< std::endl;
        os << "* hho    : "<< p.hho<< std::endl;
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


template<typename T>
struct tensors
{
public:
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>     tensor_matrix;
    tensor_matrix      siglam, gamma, xi_norm;
    size_t quad_degree;

    tensors() //quad_degree(2*degree+2)
    {}

    void
    set_quad_degree(const size_t& degree)
    {
    }
    void
    Zero(const size_t cols)
    {
        auto siglam   = tensor_matrix::Zero(2, cols);
        auto gamma    = tensor_matrix::Zero(2, cols);
        auto xi_norm  = tensor_matrix::Zero(1, cols);
    }
    template< typename MeshType>
    void
    zero_tensor(const MeshType& msh, const typename MeshType::cell& cl, const size_t& degree)
    {
        typedef MeshType                                    mesh_type;
        typedef typename MeshType::cell                     cell_type;
        typedef typename MeshType::face                     face_type;
        typedef disk::quadrature<mesh_type, cell_type>      cell_quad_type;
        typedef disk::quadrature<mesh_type, face_type>      face_quad_type;

        //WK: try to do this in the constructor.(I had problems in tensor_zero_vector, since it s asking for the input of the constructor)
        quad_degree = 2 * degree + 2;
        //std::cout << "ZERO QUAD DEGREE :"<< quad_degree << std::endl;

        auto DIM = MeshType::dimension;
        auto cq  = cell_quad_type(quad_degree);
        auto fq  = face_quad_type(quad_degree);

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        auto fqsz = num_faces * fq.integrate(msh, fcs[0]).size(); //WK: This should change if k different for each face
        auto cqsz = cq.integrate(msh, cl).size();


        //#if 0
        auto cols = cqsz + fqsz + 2 * num_faces + 1;
        siglam   = tensor_matrix::Zero(DIM, cols);
        gamma    = tensor_matrix::Zero(DIM, cols);
        xi_norm  = tensor_matrix::Zero(1  , cols);


        //#endif
        #if 0
        siglam   = tensor_matrix::Zero(DIM, cqs + ps );
        gamma    = tensor_matrix::Zero(DIM, cqs + ps );
        xi_norm  = tensor_matrix::Zero(1  , cqs + ps );
        #endif
        T a;
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
std::vector<tensors<T>>
tensor_zero_vector(const mesh<T,DIM,Storage>& msh, const size_t& degree)
{
    // WK: degree of quadrature assumed as for (Dru,Dru)
    auto    vec = std::vector<tensors<T>>(msh.cells_size());

    for (auto& cl : msh)
    {
        auto  id = cl.get_id();
        tensors<T>  tsr;
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
            auto h = diameter(msh, /*fcs[face_i]*/cl);
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
make_prod_stress_n(const  dynamic_vector<T>& st, const static_vector<T,DIM>& n )
{
    T p1(0.);
    //for(size_t i = 0; i < mesh_type::dimension; i++)
    //    p1  += st(i) * n(i);
    p1 = st.dot(n);
    return p1;
}
template<typename T>
T
make_prod_stress_n(const  dynamic_vector<T>& st, const T& n )
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

    cell_basis_type                             cell_basis;
    cell_quadrature_type                        cell_quadrature;

    face_basis_type                             face_basis;
    face_quadrature_type                        face_quadrature;

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



    vector_type
    compute(const mesh_type& msh,
                 const cell_type& cl,
                 const matrix_type& rec_oper,
                 const vector_type& uh_TF,
                 const plasticity_data<T>& pst,
                 tensors<T> & tsr)
    {

        if(pst.method == true)
            m_method_coef = 1. / pst.alpha ;
        else
            m_method_coef = 1. / (pst.alpha + pst.mu);


        //WK: for p-adaptation quad_degree is not the same for all cells

        timecounter tc;

        std::map<std::string, double> timings, timings2, timings1;

        tc.tic();

        auto fcs = faces(msh, cl);
        auto num_faces      = fcs.size();
        auto cell_range     = cell_basis.range(0,m_degree);
        auto col_range      = cell_basis.range(1,m_degree+1);
        auto row_range      = disk::dof_range(0,mesh_type::dimension);
        auto num_cell_dofs  = cell_range.size();
        auto num_face_dofs  = face_basis.size();
        dofspace_ranges    dsr(num_cell_dofs, num_face_dofs, num_faces);
        vector_type    ruh  = rec_oper * uh_TF;
        vector_type    rhs  = vector_type::Zero(dsr.total_size());
        vector_type    rhs_c  = vector_type::Zero(num_cell_dofs);
        size_t cont  =  0;
        auto   cqs   =  cell_quadrature.integrate(msh, cl);

        for (auto& qp : cqs)
        {
            timecounter tc_detail;
            tc_detail.tic();

            auto dphi       =   cell_basis.eval_gradients(msh, cl, qp.point());
            matrix_type dphi_matrix =   make_gradient_matrix(dphi);
            matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
            vector_type dphi_ruh    =   dphi_taken * ruh;
            matrix_type dphi_rec    =   dphi_taken * rec_oper;

            vector_type  xi    =  tsr.siglam.col(cont)  +  pst.alpha * dphi_ruh;
            auto xi_norm = xi.norm();

            vector_type  gamma =  vector_type::Zero(mesh_type::dimension);

            if(xi_norm > pst.yield)
                gamma   = m_method_coef * (xi_norm - pst.yield) * (xi / xi_norm);

            tsr.siglam.col(cont) +=  pst.alpha * ( dphi_ruh - gamma );
            tsr.gamma.col(cont)   =  gamma;
            tsr.xi_norm(cont)     =  xi_norm;

            //( Dv_T, sigma - alpha* gamma)_T
            if(pst.hho == 2)
            {
                matrix_type dphi_    = take(dphi_matrix, row_range, cell_range);
                //rhs_c += qp.weight() * dphi_.transpose() * (tsr.siglam.col(cont) - pst.alpha * gamma) ;
                //rhs += qp.weight() * dphi_.transpose() * (tsr.siglam.col(cont) - pst.alpha * gamma) ;
                rhs.head(num_cell_dofs)  += qp.weight() * dphi_.transpose() * (tsr.siglam.col(cont) - pst.alpha * gamma) ;
            }
            else if(pst.hho ==  1)
                rhs  += qp.weight() * dphi_rec.transpose() * (tsr.siglam.col(cont) - pst.alpha * gamma) ;
            else
                throw std::invalid_argument("Invalid discretization name for the plasticity term. Review the name in parametes.txt or pst.hho");

            ++cont;
        }

        size_t face_basis_offset = num_cell_dofs;
        size_t offset    = cqs.size();
        //#if 0
        // Boundary term
        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fc   = fcs.at(face_i);
            auto fqs  = face_quadrature.integrate(msh, fc);

            //#if 0
            auto n = normal(msh, cl, fc);
            vector_type    rhs_f  = vector_type::Zero(face_basis.size());

            for (auto& qp : fqs)
            {
                auto f_phi  = face_basis.eval_functions(msh, fc, qp.point());
                auto c_phi  = cell_basis.eval_functions(msh, cl, qp.point());
                auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point());

                matrix_type c_dphi_matrix =   make_gradient_matrix(c_dphi);
                matrix_type c_dphi_taken  =   take(c_dphi_matrix, row_range, col_range);
                vector_type c_dphi_ruh    =   c_dphi_taken * ruh;

                // xi = sigma + alpha *Dr(u)
                vector_type  xi    =  tsr.siglam.col(cont)  +  pst.alpha * c_dphi_ruh;
                auto xi_norm = xi.norm();

                vector_type  gamma =  vector_type::Zero(mesh_type::dimension);
                if(xi_norm > pst.yield)
                    gamma  = m_method_coef * (xi_norm - pst.yield)* (xi / xi_norm);

                tsr.siglam.col(cont) += pst.alpha * ( c_dphi_ruh - gamma );
                tsr.gamma.col(cont)   = gamma;
                tsr.xi_norm(cont)     = xi_norm;

                // ((sigma - alpha*gamma) * n
                vector_type str =  tsr.siglam.col(cont) - pst.alpha * gamma;
                scalar_type p1  =  make_prod_stress_n( str , n);
                if(pst.hho == 2)
                {
                    // ( (sigma - alpha*gamma) * n, vT )_F
                    for (size_t i = 0; i < num_cell_dofs; i++)
                        rhs(i)   -= qp.weight() * p1 * c_phi.at(i) ;

                    // ( (sigma - alpha*gamma) * n, vF )_F
                    for (size_t i = 0, ii = current_face_range.min(); i < current_face_range.size(); i++, ii++)
                        rhs(ii)  += qp.weight() * p1 * f_phi.at(i);
                    #if 0
                    // ( (sigma - alpha*gamma) * n, vT )_F
                    for (size_t i = 0; i < num_cell_dofs; i++)
                        rhs_c(i)   -= qp.weight() * p1 * c_phi.at(i) ;

                    // ( (sigma - alpha*gamma) * n, vF )_F
                    for (size_t i = 0; i < num_face_dofs; i++)
                        rhs_f(i)  += qp.weight() * p1 * f_phi.at(i);
                    #endif
                }
                ++cont;
            }

            //rhs.block(face_basis_offset, 0, num_face_dofs,1) = rhs_f;
            face_basis_offset += num_face_dofs;
            //#endif
            offset   += fqs.size();
        }

        //rhs.block(0, 0, num_cell_dofs,1) = rhs_c;



        //#endif

        //#if 0

        // tensor values in vertices
        auto pts = points(msh,cl);


        if(offset != cont)
            throw std::invalid_argument("Invalid indices in storing values at Gauss points");

        // Value at the vertices
        for(size_t i = 0, j = offset; i < pts.size(); i++, j++)
        {
            //std::cout << "j = "<< j<< std::endl;
            auto dphi   = cell_basis.eval_gradients(msh, cl, pts.at(i));
            matrix_type dphi_matrix = make_gradient_matrix(dphi);
            matrix_type dphi_taken  = take(dphi_matrix, row_range, col_range);
            vector_type dphi_ruh    = dphi_taken * ruh;
            vector_type xi = tsr.siglam.col(j)  +  pst.alpha * dphi_ruh;
            auto  xi_norm  = xi.norm();
            vector_type gamma =   vector_type::Zero(mesh_type::dimension);

            if(xi_norm > pst.yield)
                gamma = m_method_coef * (xi_norm - pst.yield) * (xi / xi_norm);

            tsr.siglam.col(j) += pst.alpha * ( dphi_ruh - gamma );
            tsr.gamma.col(j)   = gamma;
            tsr.xi_norm(j)     = xi_norm;
        }

        // Value at the barycenter of the faces
        offset += pts.size();
        for(size_t i = 0, j = offset; i < num_faces; i++, j++)
        {
            //std::cout << "j = "<< j<< std::endl;
            auto face   = fcs.at(i);
            auto fbar   = barycenter(msh, face);
            auto dphi   = cell_basis.eval_gradients(msh, cl, fbar);
            matrix_type dphi_matrix = make_gradient_matrix(dphi);
            matrix_type dphi_taken  = take(dphi_matrix, row_range, col_range);
            vector_type dphi_ruh    = dphi_taken * ruh;
            vector_type xi = tsr.siglam.col(j)  +  pst.alpha * dphi_ruh;
            auto  xi_norm  = xi.norm();
            vector_type gamma =   vector_type::Zero(mesh_type::dimension);

            if(xi_norm > pst.yield)
                gamma = m_method_coef * (xi_norm - pst.yield) * (xi / xi_norm);

            tsr.siglam.col(j) += pst.alpha * ( dphi_ruh - gamma );
            tsr.gamma.col(j)   = gamma;
            tsr.xi_norm(j)     = xi_norm;
        }

        // Value at the barycenter of the cell
        offset += num_faces;
        auto cbar   = barycenter(msh, cl);
        auto dphi   = cell_basis.eval_gradients(msh, cl, cbar);
        matrix_type dphi_matrix = make_gradient_matrix(dphi);
        matrix_type dphi_taken  = take(dphi_matrix, row_range, col_range);
        vector_type dphi_ruh    = dphi_taken * ruh;
        vector_type xi = tsr.siglam.col(offset)  +  pst.alpha * dphi_ruh;
        auto  xi_norm  = xi.norm();
        vector_type gamma =   vector_type::Zero(mesh_type::dimension);

        if(xi_norm > pst.yield)
            gamma = m_method_coef * (xi_norm - pst.yield) * (xi / xi_norm);

        tsr.siglam.col(offset) += pst.alpha * ( dphi_ruh - gamma );
        tsr.gamma.col(offset)   = gamma;
        tsr.xi_norm(offset)     = xi_norm;

        assert(offset + 1 == tsr.xi_norm.cols());

        //#endif
        return rhs;
    }


template< typename Solution>
vector_type
compute(const mesh_type& msh,
             const cell_type& cl,
             const matrix_type& rec_oper,
             const vector_type& uh_TF,
             const plasticity_data<T>& pst,
             const Solution& solution)
{
    // Diffusion

    //std::cout << "/**********   CELL = "<<msh.lookup(cl)<<"*********/" << std::endl;
    //WK: for p-adaptation quad_degree is not the same for all cells
    timecounter tc;

    std::map<std::string, double> timings, timings2, timings1;

    tc.tic();

    auto fcs = faces(msh, cl);
    auto num_faces      = fcs.size();
    auto cell_range     = cell_basis.range(0,m_degree);
    auto col_range      = cell_basis.range(1,m_degree+1);
    auto row_range      = disk::dof_range(0,mesh_type::dimension);
    auto num_cell_dofs  = cell_range.size();
    auto num_face_dofs  = face_basis.size();
    dofspace_ranges    dsr(num_cell_dofs, num_face_dofs, num_faces);
    matrix_type    ruh  = rec_oper * uh_TF;
    vector_type    rhs_m1  = vector_type::Zero(dsr.total_size());
    vector_type    rhs_m2  = vector_type::Zero(dsr.total_size());
    vector_type    rhs_c  = vector_type::Zero(num_cell_dofs);
    size_t cont   =  0;
    auto   cqs    =  cell_quadrature.integrate(msh, cl);
    auto   pts    =  points( msh , cl);


    T x0 = 0.5;
    size_t number = 0;
    for(auto& p : pts)
    {
        if(p.x() < x0)
        number = 1;
        if(p.x() > x0)
        number = 2;
    }
    std::cout << "number( "<< cl.get_id()<<") = "<< number    << std::endl;
    if(!(number ==2 || number ==1) )
        throw std::invalid_argument("Invalid number domain.");

    for (auto& qp : cqs)
    {
        timecounter tc_detail;
        tc_detail.tic();

        auto dphi   =   cell_basis.eval_gradients(msh, cl, qp.point());
        matrix_type dphi_matrix =   make_gradient_matrix(dphi);
        matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
        vector_type dphi_ruh    =   dphi_taken * ruh;
        matrix_type dphi_rec    =   dphi_taken * rec_oper;

        //( Dv_T, sigma - alpha* gamma)_T
        if(pst.hho == 2)
        {
            matrix_type dphi_    = take(dphi_matrix, row_range, cell_range);
            rhs_c.head(num_cell_dofs)  += qp.weight() * dphi_.transpose() * solution.st( qp.point(), number);
        }
        else if(pst.hho == 1)
            rhs_m1 += qp.weight() * dphi_rec.transpose()  * solution.st( qp.point(), number);
        else
            throw std::invalid_argument("Invalid discretization name for the plasticity term. Review the name in parametes.txt or pst.hho");
        ++cont;
    }

    //std::cout << "rhs_m1 : [ "<< rhs_m1 << " ]"<< std::endl;
    //std::cout << "rhs_m2 : [ "<< rhs_m2 << " ]"<< std::endl;

    size_t face_basis_offset = num_cell_dofs;
    size_t offset    = cqs.size();

    //#if 0
    // Boundary term
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        auto current_face_range = dsr.face_range(face_i);
        auto fc   = fcs.at(face_i);
        auto fqs  = face_quadrature.integrate(msh, fc);

        //#if 0
        auto n = normal(msh, cl, fc);
        vector_type    rhs_f  = vector_type::Zero(face_basis.size());

        for (auto& qp : fqs)
        {
            auto f_phi  = face_basis.eval_functions(msh, fc, qp.point());
            auto c_phi  = cell_basis.eval_functions(msh, cl, qp.point());
            auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point());

            matrix_type c_dphi_matrix =   make_gradient_matrix(c_dphi);
            matrix_type c_dphi_taken  =   take(c_dphi_matrix, row_range, col_range);
            vector_type c_dphi_ruh    =   c_dphi_taken * ruh;

            size_t  number_f;
            if(!msh.is_boundary(fc))
            {
                size_t number_nc;
                auto   neighbor_id =  face_owner_cells_ids(msh, fc, cl);
                auto   ngh   =  *std::next(msh.cells_begin(), neighbor_id);
                auto   pts_ngh    =  points( msh , ngh);

                for(auto& p : pts_ngh)
                {
                    if(p.x() < x0)
                        number_nc = 1;
                    if(p.x() > x0)
                        number_nc = 2;
                }
                //std::cout << "number_nc( "<< ngh.get_id()<<") = "<< number    << std::endl;
                if(!(number ==2 || number ==1) )
                throw std::invalid_argument("Invalid number domain.");

                //Esto escoge el valor de la derecha para tau_F
                number_f = number;
                if(number != number_nc )
                {
                    number_f == 1;
                    std::cout << " changing number in cell "<< cl.get_id() << std::endl;
                }
            }
            else
            {
                number_f = number;
                //std::cout << "number_f( "<< cl.get_id()<<") = "<< number_f << std::endl;
            }


                // ((sigma - alpha*gamma) * n
            vector_type str   =  solution.st( qp.point(), number);
            vector_type str_f =  solution.st( qp.point(), number_f);
            scalar_type p1_c  =  make_prod_stress_n( str , n);
            scalar_type p1_f  =  make_prod_stress_n( str_f , n);

            if(pst.hho == 2)
            {
                //#if 0
                // ( (sigma - alpha*gamma) * n, vT )_F
                //for (size_t i = 0; i < num_cell_dofs; i++)
                //    rhs_m2(i)   -= qp.weight() * p1_c * c_phi.at(i) ;

                // ( (sigma - alpha*gamma) * n, vF )_F
                //for (size_t i = 0, ii = current_face_range.min(); i < current_face_range.size(); i++, ii++)
                //    rhs_m2(ii)  += qp.weight() * p1_f * f_phi.at(i);
                //#endif
                //#if 0
                // ( (sigma - alpha*gamma) * n, vT )_F
                for (size_t i = 0; i < num_cell_dofs; i++)
                    rhs_c(i)   -= qp.weight() * p1_c * c_phi.at(i) ;

                // ( (sigma - alpha*gamma) * n, vF )_F
                for (size_t i = 0; i < num_face_dofs; i++)
                    rhs_f(i)  += qp.weight() * p1_f * f_phi.at(i);
                //#endif
            }

            //std::cout << " ** rhs_c : [ "<< rhs_c.transpose() << " ]"<< std::endl;
            //std::cout << " ** rhs_f : [ "<< rhs_f.transpose() << " ]"<< std::endl;

            ++cont;
        }
        //std::cout << "face_basis_offset : "<< face_basis_offset << std::endl;
        rhs_m2.block(face_basis_offset, 0, num_face_dofs,1) += rhs_f;
        face_basis_offset += num_face_dofs;
        //#endif
        offset   += fqs.size();

        std::cout << " * rhs_m1 : [ "<< rhs_m1.transpose() << " ]"<< std::endl;
        std::cout << " * rhs_m2 : [ "<< rhs_m2.transpose() << " ]"<< std::endl;
    }

    rhs_m2.block(0, 0, num_cell_dofs,1) = rhs_c;

    std::cout << " - rhs_m1 : [ "<< rhs_m1.transpose() << " ]"<< std::endl;
    std::cout << " - rhs_m2 : [ "<< rhs_m2.transpose() << " ]"<< std::endl;

    //#endif
    if(offset != cont)
        throw std::invalid_argument("Invalid indices in storing values at Gauss points");

    if(pst.hho == 1)
        return rhs_m1;
    else if(pst.hho == 2)
        return rhs_m2;
    else
        throw std::invalid_argument("Invalid discretization name for the plasticity term. Review the name in parametes.txt or pst.hho");
}

};
#if 0
template< typename T, typename Storage, typename CellType>
T
set_cell_number(const  disk::mesh<T,2,Storage>& msh, const CellType& cl)
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
}
#endif


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
            T  ret, xO(0.), xL(pst.Lref);
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
            T  xO(0.);
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
        is_exact = false;
        identity = "square";
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
        is_exact = false;
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
                    ret(1) = 0.;
                    break;
                case 2:
                    ret(0) = 0.5 * p.x();
                    ret(1) = c1;
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
