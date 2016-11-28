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

#include "common/eigen.hpp"
#include "bases/bases_utils.hpp"
#include "bases/bases_ranges.hpp"

namespace disk {

template<typename T>
struct plasticity_data
{
    plasticity_data(): Lref(1.), Vref(1.), Bn(0.1),mu(1.), alpha(1.), f(1.),method(true)
    {}
    T f;                    //WK: Cuidado porque f deberia ser el valor externo de la fuente.
    T Lref;                 /* Charactetistic length */
    T Vref;                 /* Reference velocity */
    T Bn;                   /* Bingham number */
    T mu;                   /* viscosity */
    T alpha;
    T yield;
    bool   method;
};

template<typename T>
struct tensors
{
public:
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>     tensor_matrix;
    tensor_matrix      siglam, gamma, xi_norm;

    tensors(){};

    template< typename MeshType>
    void
    zero_tensor(const MeshType& msh, const typename MeshType::cell& cl, const size_t degree)
    {
        typedef MeshType                                    mesh_type;
        typedef typename MeshType::cell                     cell_type;
        typedef disk::quadrature<mesh_type, cell_type>      cell_quad_type;

        auto DIM = MeshType::dimension;
        auto cq  = cell_quad_type(2*degree+2);
        auto cqs = cq.integrate(msh, cl).size();
        siglam   = tensor_matrix::Zero(DIM, cqs);
        gamma    = tensor_matrix::Zero(DIM, cqs);
        xi_norm  = tensor_matrix::Zero(1  , cqs);
        T a;
    }
};

template<typename T, size_t DIM, typename Storage>
std::vector<tensors<T>>
zero_tensor_vector(const mesh<T,DIM,Storage>& msh, const size_t degree)
{
    // WK: degree of quadrature assumed as for (Dru,Dru)
    auto    vec = std::vector<tensors<T>>(msh.cells_size());

    for (auto& cl : msh)
    {
        auto  id = msh.lookup(cl);
        tensors<T>  tsr;
        tsr.zero_tensor(msh,cl,degree);
        vec[id]   =  tsr;
    }
    return vec;
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
    auto
    compute(const mesh_type& msh, const cell_type& cl,
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
        matrix_type K_TF = local_mat.topRightCorner(cell_size, face_size);
        matrix_type K_FT = local_mat.bottomLeftCorner(face_size, cell_size);
        matrix_type K_FF = local_mat.bottomRightCorner(face_size, face_size);
        vector_type bP_F = plastic_rhs.tail(face_size);
        vector_type bP_T = plastic_rhs.head(cell_size);


        assert(K_TT.cols() == cell_size);
        assert(K_TT.cols() + K_TF.cols() == local_mat.cols());
        assert(K_TT.rows() + K_FT.rows() == local_mat.rows());
        assert(K_TF.rows() + K_FF.rows() == local_mat.rows());
        assert(K_FT.cols() + K_FF.cols() == local_mat.cols());
        assert(bP_T.rows() + bP_F.rows() == plastic_rhs.rows());

        auto K_TT_ldlt = K_TT.llt();
        matrix_type AL = K_TT_ldlt.solve(K_TF);
        vector_type bL = K_TT_ldlt.solve(cell_rhs - bP_T);


        matrix_type AC = K_FF - K_FT * AL;
        vector_type bC = /* no projection on faces, eqn. 26*/ - bP_F  - K_FT * bL;

        return std::make_pair(AC, bC);
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

template<typename MeshType, //typename CellBasisType, typename CellQuadType,
                        typename FaceBasisType, typename FaceQuadType>
class assembler_pst
{
    typedef MeshType                            mesh_type;
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

    assembler_pst()               = delete;
    assembler_pst(const assembler_pst&) = delete;
    assembler_pst(assembler_pst&&)      = delete;

    assembler_pst(const MeshType& msh, const size_t degree)
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
                l2g[pos+i] = face_offset+i;
        }

        assert(lc.first.rows() == lc.first.cols());
        assert(lc.first.rows() == lc.second.size());
        assert(lc.second.size() == l2g.size());

        //std::cout << lc.second.size() << " " << l2g.size() << std::endl;

        for (size_t i = 0; i < lc.first.rows(); i++)
        {
            for (size_t j = 0; j < lc.first.cols(); j++)
                m_triplets.push_back( triplet_type( l2g.at(i), l2g.at(j), lc.first(i,j) ) );

            rhs(l2g.at(i)) += lc.second(i);
        }
    }
    template<typename Function, typename PlasticData>
    void
    impose_boundary_conditions(mesh_type& msh, const Function& bc, const PlasticData& pst)
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

                    rhs_f(i) += qp.weight() * mm_prod(f_phi[i], bc(qp.point(), pst));
                }
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
high_order_reconstruction_nopre(const mesh<T, DIM, Storage>&  msh,  const typename mesh<T, DIM, Storage>::cell& cl,
                                const dynamic_matrix<T>& rec_oper,
                                const dynamic_vector<T>& v,
                                const size_t degree)
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

template<typename T,typename Mesh,typename CellBasisType, typename CellQuadType,
                       typename FaceBasisType>
class plasticity
{
    typedef Mesh                            mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;

    typedef CellBasisType                   cell_basis_type;
    typedef CellQuadType                    cell_quadrature_type;
    typedef FaceBasisType                   face_basis_type;

    typedef dynamic_matrix<scalar_type>     matrix_type;
    typedef dynamic_vector<scalar_type>     vector_type;
    typedef typename CellBasisType::gradient_value_type gradient_value_type_pr;

    cell_basis_type                         cell_basis;
    cell_quadrature_type                    cell_quadrature;

    face_basis_type                         face_basis;
    size_t                                  m_degree;
    T                                       m_yield;
    T                                       m_alpha;
    T                                       m_method_coef;
public:
    vector_type     rhs;


    plasticity()                = delete;
    plasticity(const size_t degree, const plasticity_data<T>& pst)
    //WK: for mu as tensor one should check this part
        : m_degree(degree), m_alpha(pst.alpha), m_yield(pst.yield)
    {
        cell_basis          = cell_basis_type(m_degree+1);
        cell_quadrature     = cell_quadrature_type(2*(m_degree+1));
        face_basis          = face_basis_type(m_degree);

        if(pst.method == true)
            m_method_coef = 1. / m_alpha ;
        else
            m_method_coef = 1. / (m_alpha + pst.mu);
    }

    void compute(const mesh_type& msh,
                 const cell_type& cl,
                 const matrix_type& rec_oper,
                 const vector_type& uh_TF,
                 tensors<T> & tsr)
    {
        auto fcs = faces(msh, cl);
        auto num_faces       = fcs.size();
        auto num_cell_dofs   = cell_basis.range(0,m_degree).size();
        auto num_face_dofs   = face_basis.size();

        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        rhs  = vector_type::Zero(dsr.total_size());

        size_t cont     = 0;
        auto cell_id    = msh.lookup(cl);
        auto cell_quadpoints = cell_quadrature.integrate(msh, cl);

        for (auto& qp : cell_quadpoints)
        {
            auto dphi       =   cell_basis.eval_gradients(msh, cl, qp.point());
            auto col_range  =   cell_basis.range(1,m_degree+1);
            auto row_range  =   disk::dof_range(0,mesh_type::dimension);

            matrix_type dphi_matrix =   make_gradient_matrix(dphi);
            matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
            matrix_type dphi_rec    =   dphi_taken * rec_oper;

            vector_type dphi_rec_uh =   dphi_rec*uh_TF;

            vector_type siglam   =   tsr.siglam.col(cont); // WK: Tener cuidado con la sacada de sigma con 2D
            vector_type    xi    =   siglam  +  m_alpha * dphi_rec_uh;
            T xi_norm      =   xi.norm();
            tsr.xi_norm(cont) = xi_norm;

            if(xi_norm > m_yield)
            {
                vector_type  gamma;
                gamma   = m_method_coef * (xi_norm - m_yield)* (xi / xi_norm);
                rhs    +=   qp.weight() * (siglam  - m_alpha * gamma).transpose() * dphi_rec ;
                tsr.siglam.col(cont) += m_alpha * ( dphi_rec_uh - gamma );
                tsr.gamma.col(cont)   = gamma;
            }
            else
            {
                rhs +=  qp.weight() * siglam * dphi_rec;
                tsr.siglam.col(cont) += m_alpha * dphi_rec_uh;
                tsr.gamma.col(cont)   = vector_type::Zero(cell_quadpoints.size());
            }
            ++cont;
        }
    }
};


template<typename MeshType>
class stress_based_mesh
{};

template<typename T, typename Storage>
class stress_based_mesh<mesh<T,1,Storage>>
{
public:
    typedef mesh<T,1,Storage> mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::cell_iterator       cell_itor;

    std::vector<bool> cl_marks_vector;


    stress_based_mesh(const mesh_type& msh,const std::vector<tensors<T>>& tsr_vec,const T yield)
    {
        std::vector<bool> vec(msh.cells_size(),false);

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto mat = tsr_vec.at(i).xi_norm;
            if(std::abs(mat.maxCoeff()) >=  yield)
            {
                if (std::abs(mat.minCoeff())< yield)
                    vec.at(i) = true;
            }
        }

         cl_marks_vector = vec;
    }

    void re_populate_mesh(mesh_type & msh,
                            const std::vector<tensors<T>>& tsr_vec,
                            const T yield,
                            const int degree)
    {
        mesh_type re_msh;

        auto storage    = msh.backend_storage();
        auto storage_rm = re_msh.backend_storage();
        auto num_cells  = msh.cells_size();
        size_t nds_counter = 0;

        for(auto& cl:msh)
        {
            size_t id  = msh.lookup(cl);
            auto xi_norms_storage =  tsr_vec[id].xi_norm;
            auto cmkr  = cl_marks_vector.at(id);
            //WK: This is only for 1D
            auto lmkr  = (id != 0)?  cl_marks_vector.at(id-1) : false;
            auto rmkr  = (id != msh.cells_size()-1)?  cl_marks_vector.at(id+1) : false;

            auto cell_pts   = points(msh, cl);
            storage_rm->points.push_back(cell_pts.at(0));

            if(cmkr || (lmkr ||rmkr))
            {
                auto new_pt = (cell_pts.at(1) + cell_pts.at(0))/2.;
                /* Points */
                storage_rm->points.push_back(new_pt);
            }
            if(id == msh.cells_size() -1 )
                 storage_rm->points.push_back(cell_pts.at(1));
        } //for cells

        /* nodes */
        auto num_rm_pts = storage_rm->points.size();
        for(size_t i = 0; i < num_rm_pts -1; i++)
        {
            auto n0 = typename node_type::id_type(i);
            auto n1 = typename node_type::id_type(i + 1);
            auto e = edge_type{{n0, n1}};

            std::vector<point_identifier<1>> pts(2);
            pts[0] = point_identifier<1>(i);
            pts[1] = point_identifier<1>(i + 1);
            e.set_point_ids(pts.begin(), pts.end());
            storage_rm->edges.push_back(e);
        }

        for (size_t i = 0; i < num_rm_pts; i++)
            storage_rm->nodes.push_back(node_type(point_identifier<1>(i)));

        storage_rm->boundary_nodes.resize(num_rm_pts);
        storage_rm->boundary_nodes.at(0) = true;
        storage_rm->boundary_nodes.at(num_rm_pts - 1) = true;

        msh = re_msh;
    }

};

template<typename T, typename Storage>
class stress_based_mesh<mesh<T,2,Storage>>
{

public:
    typedef mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::cell                cell_type;
    typedef typename mesh_type::cell_iterator       cell_itor;

    bool m_hanging_nodes;
    std::vector<size_t> cl_marks_vector;
    //std::vector<std::pair<bool,typename point<T,2>::id_type>> fc_marks_vector;
    std::vector<std::pair<bool,int>> fc_marks_vector;
    std::vector<std::array<ident_impl_t, 4>>        m_edges;
    std::vector<std::array<ident_impl_t, 2>>        m_boundary_edges;

    std::vector<std::array<ident_impl_t, 3>>        m_triangles;
    std::vector<std::array<ident_impl_t, 4>>        m_quadrangles;
    std::vector<std::array<ident_impl_t, 5>>        m_pentagons;
    std::vector<std::array<ident_impl_t, 6>>        m_hexagons;

    stress_based_mesh(const mesh_type& msh,
                        const std::vector<tensors<T>>& tsr_vec,
                        const T yield,
                        bool hanging_ndoes):m_hanging_nodes(hanging_ndoes)
    {
        if(m_hanging_nodes)
            cl_marks_vector = cells_marker_with_hanging_nodes(msh,tsr_vec,yield);
        else
            cl_marks_vector = cells_marker(msh,tsr_vec,yield);
            for(size_t cm = 0; cm < cl_marks_vector.size(); cm++)
            {
                size_t b = cl_marks_vector.at(cm);
                std::cout<<b<<"  ";
            }
            std::cout<<" ]"<<std::endl;

        fc_marks_vector.resize(msh.faces_size());
        for(size_t i = 0; i < fc_marks_vector.size(); i++)
            fc_marks_vector.at(i).first = false;
    }
    template< size_t N>
    struct n_gon
    {
       std::array<size_t, N>  p;
       std::array<bool, N>    b;
    };

    std::vector<size_t>
    cells_marker_with_hanging_nodes(const mesh_type& msh,
                        const std::vector<tensors<T>>& tsr_vec,
                        const T yield)
    {
        std::vector<size_t> vec(msh.cells_size());

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto mat = tsr_vec.at(i).xi_norm;
            if(std::abs(mat.maxCoeff()) >=  yield)
            {
                if (std::abs(mat.minCoeff()) < yield)
                    vec.at(i) = 1;
            }
        }
        return vec;
    }

    std::vector<size_t>
    cells_marker(const mesh_type& msh,
            const std::vector<tensors<T>>& tsr_vec,
            const T yield)
    {
        std::vector<size_t> vec(msh.cells_size());

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto mat = tsr_vec.at(i).xi_norm;
            if(std::abs(mat.maxCoeff()) >=  yield)
            {
                if (std::abs(mat.minCoeff())< yield)
                    vec.at(i) = 3;
            }
        }
        return vec;
    }
   template<typename U, size_t N, size_t M>
    void store(const n_gon<N> & t,  std::vector<std::array<U, M>>& m_n_gons)
    {
        if(M == N)
            m_n_gons.push_back(t.p);
        else
            throw std::invalid_argument("Polygon cannot be store");

         // Mejorar todo esto!!! para que no quede repetido
         for(size_t i = 0; i < N - 1; i++)
         {
             auto a = t.p[i];
             auto b = t.p[i+1];
             if(b < a)
                 std::swap(a,b);
             m_edges.push_back({a, b});
             if(t.b[i])
                 m_boundary_edges.push_back({a , b});
         }
         auto a = t.p[0];
         auto b = t.p[N-1];
         if(b < a)
             std::swap(a,b);
         m_edges.push_back({a , b});
         if(t.b[N -1])
             m_boundary_edges.push_back({a , b });
     }

    void
    marker_with_hanging_nodes(mesh_type & msh)
    {
        auto cmv_temp = cl_marks_vector;

        auto storage = msh.backend_storage();
        for(auto& cl : msh)
        {
            auto cl_id   = msh.lookup(cl);
            auto cl_mark = cl_marks_vector.at(cl_id);

            if(cl_mark == 1)
            {
                auto fcs = faces(msh,cl);
                for(auto& fc:fcs)
                {
                    if(!msh.is_boundary(fc))
                    {
                        auto  neighbor_id   = face_owner_cells_ids(msh,fc,cl);
                        cmv_temp.at(neighbor_id) = 1;
                        auto  neighbor      = *(msh.cells_begin() + size_t(neighbor_id));
                        auto  neighbor_fcs  = faces(msh,neighbor);
                        for(auto& nfc : neighbor_fcs)
                        {
                            auto  nfc_id    = msh.lookup(nfc);

                            if(!fc_marks_vector.at(nfc_id).first)
                            {
                                fc_marks_vector.at(nfc_id).first  = true;
                                auto bar = barycenter(msh,nfc);
                                storage ->points.push_back(bar);
                                fc_marks_vector.at(nfc_id).second = storage->points.size() -1;
                                //WK: this should be with an id_type instead of using ints

                                if(!msh.is_boundary(nfc))
                                {
                                    auto  nn_id   = face_owner_cells_ids(msh,nfc,cl);
                                    cmv_temp.at(nn_id) = 1;
                                }
                            }
                        }
                    }
                    else
                    {
                        auto  fc_id    = msh.lookup(fc);
                        //std::cout << "inside marker, marking bundary = "<<fc_id << std::endl;
                        fc_marks_vector.at(fc_id).first  = true;
                        auto bar = barycenter(msh,fc);
                        storage ->points.push_back(bar);
                        fc_marks_vector.at(fc_id).second = storage->points.size() -1;
                    }
                }
            }
        }
    }

    void
    set_marks(mesh_type & msh, const cell_type & cl)
    {
        typedef typename cell_type::id_type         cell_id_type;

        auto storage = msh.backend_storage();
        auto cl_id   = msh.lookup(cl);
        auto cl_mark = cl_marks_vector.at(cl_id);

        if(cl_mark > 1  &  cl_mark < 4)
        {
            auto fcs = faces(msh,cl);
            for(auto& fc:fcs)
            {
                if(!msh.is_boundary(fc))
                {
                    auto  neighbor_id = face_owner_cells_ids(msh,fc,cl);
                    auto  ngh_mark    = cl_marks_vector.at(neighbor_id);
                    auto  fc_id    = msh.lookup(fc);
                    if(!fc_marks_vector.at(fc_id).first)
                    {
                        fc_marks_vector.at(fc_id).first  = true;
                        auto bar = barycenter(msh,fc);
                        storage ->points.push_back(bar);
                        fc_marks_vector.at(fc_id).second = storage->points.size() -1;
                    }
                    if(ngh_mark < cl_mark)
                    {
                        if(ngh_mark == 1)
                            cl_marks_vector.at(neighbor_id) = 2;
                        else
                            cl_marks_vector.at(neighbor_id) = cl_mark - 1;
                        auto  ngh = *(msh.cells_begin() + size_t(neighbor_id));
                        set_marks(msh,ngh);
                    }

                        #if 0
                        auto  ngh      = *(msh.cells_begin() + size_t(ngh_id));
                        auto  ngh_fcs  = faces(msh,ngh);

                        for(auto& nfc : ngh_fcs)
                        {
                            auto  nfc_id  = msh.lookup(nfc);

                            if(!fc_marks_vector.at(nfc_id).first)
                            {
                                fc_marks_vector.at(nfc_id).first  = true;
                                auto bar = barycenter(msh,nfc);
                                storage ->points.push_back(bar);
                                fc_marks_vector.at(nfc_id).second = storage->points.size() ;
                            }
                        }
                        #endif

                }
                else
                {
                    auto  fc_id    = msh.lookup(fc);
                    fc_marks_vector.at(fc_id).first  = true;
                    auto bar = barycenter(msh,fc);
                    storage ->points.push_back(bar);
                    fc_marks_vector.at(fc_id).second = storage->points.size() -1;
                }
            }
        }

        //else
        //    throw std::logic_error("shouldn't have arrived here");
    }
    void
    marker(mesh_type& msh)
    {
        auto cmv_initial = cl_marks_vector;
        auto storage = msh.backend_storage();
        for(auto& cl : msh)
        {
            auto cl_id   = msh.lookup(cl);
            bool cl_mark = cmv_initial.at(cl_id);
            if(cl_mark)
                set_marks(msh, cl);
        }
    }

   template<size_t N, typename NGon, typename FaceIds>
   n_gon<N>
   n_gon_with_hanging_node(const size_t num_fcs,const FaceIds& fcs_ids,const NGon& t)
   {
        n_gon<N> nt;

       for(size_t i = 0, j = 0; i < num_fcs; i++, j++)
       {
           auto fc_mark = fc_marks_vector.at(fcs_ids[i]).first;
           nt.p[j] = t.p[i];
           nt.b[j] = t.b[i];
           if(fc_mark)
            {
                j++;
                nt.p[j] = fc_marks_vector.at(fcs_ids[i]).second;
                nt.b[j] = t.b[i];
            }
       }
       return nt;
   }


    void refine_single(const mesh_type& msh, const cell_type& cl)
    {
        n_gon<3> t;

        auto storage = msh.backend_storage();
        auto pts_ids = cl.point_ids();
        auto fcs_ids = cl.faces_ids();
        auto num_fcs = fcs_ids.size();
        auto cl_id   = msh.lookup(cl);
        auto cl_mark = cl_marks_vector.at(cl_id);

        for(size_t i = 0; i < num_fcs; i++)
        {
            auto id  = fcs_ids.at(i);
            auto fc_mark = fc_marks_vector.at(id).first;
            t.p[i]   = pts_ids.at(i);
            t.b[i]   = storage->boundary_edges.at(id);
        }

        if(cl_mark > 1 )
        {
            n_gon<3> t0,t1,t2,t3;
            auto bar0pos = fc_marks_vector.at(fcs_ids[0]).second;
            auto bar1pos = fc_marks_vector.at(fcs_ids[1]).second;
            auto bar2pos = fc_marks_vector.at(fcs_ids[2]).second;

            t0.p[0] = t.p[0];   t0.p[1] = bar0pos;  t0.p[2] = bar2pos;
            t1.p[0] = bar0pos;  t1.p[1] = t.p[1];   t1.p[2] = bar1pos;
            t2.p[0] = bar2pos;  t2.p[1] = bar1pos;  t2.p[2] = t.p[2];
            t3.p[0] = bar0pos;  t3.p[1] = bar1pos;  t3.p[2] = bar2pos;

            t0.b[0] = t.b[0];   t0.b[1] = false;    t0.b[2] = t.b[2];
            t1.b[0] = t.b[0];   t1.b[1] = t.b[1];   t1.b[2] = false;
            t2.b[0] = false;    t2.b[1] = t.b[1];   t2.b[2] = t.b[2];
            t3.b[0] = false;    t3.b[1] = false;    t3.b[2] = false;

            store(t0,m_triangles);
            store(t1,m_triangles);
            store(t2,m_triangles);
            store(t3,m_triangles);
        }

        if(cl_mark == 1)
        {
            n_gon<3> t0,t1;

            t0.p[0] = t.p[0];   t0.p[1] = t.p[1];   t0.p[2] = t.p[2];
            t1.p[0] = t.p[1];   t1.p[1] = t.p[2];   t1.p[2] = t.p[0];

            t0.b[0] = t.b[0];   t0.b[1] = t.b[1];   t0.b[2] = t.b[2];
            t1.b[0] = t.b[1];   t1.b[1] = t.b[2];   t1.b[2] = t.b[0];

            std::array<int, 3> permut = {2,0,1};

            for(size_t i = 0; i < num_fcs; i++)
            {
                auto id  = fcs_ids.at(i);
                auto fc_mark = fc_marks_vector.at(id).first;

                if(fc_mark)
                {
                    auto barpos = fc_marks_vector.at(fcs_ids[i]).second;
                    t0.p[i] = barpos;           t1.p[i] = barpos;
                    t0.b[permut[i]] = false;    t1.b[i] = false;
                }
            }
            store(t0,m_triangles);
            store(t1,m_triangles);
        }
    }


    void refine_single_with_hanging_nodes(const mesh_type& msh, const cell_type& cl)
    {
            auto storage = msh.backend_storage();

            auto count   = 0;
            auto pts_ids = cl.point_ids();
            auto fcs_ids = cl.faces_ids();
            auto num_fcs = fcs_ids.size();

            //for triangles
            n_gon<3> t;

            for(size_t i = 0; i < num_fcs; i++)
            {

                auto id  = fcs_ids.at(i);
                auto fc_mark = fc_marks_vector.at(id).first;
                t.p[i]   = pts_ids.at(i);
                t.b[i]   = storage->boundary_edges.at(id);
                if(fc_mark)
                    count++;
            }

            auto num_split_fcs =  count;

            if( num_split_fcs == 3)
            {
                n_gon<3> t0,t1,t2,t3;
                auto bar0pos = fc_marks_vector.at(fcs_ids[0]).second;
                auto bar1pos = fc_marks_vector.at(fcs_ids[1]).second;
                auto bar2pos = fc_marks_vector.at(fcs_ids[2]).second;

                t0.p[0] = t.p[0];   t0.p[1] = bar0pos;  t0.p[2] = bar2pos;
                t1.p[0] = bar0pos;  t1.p[1] = t.p[1];   t1.p[2] = bar1pos;
                t2.p[0] = bar2pos;  t2.p[1] = bar1pos;  t2.p[2] = t.p[2];
                t3.p[0] = bar0pos;  t3.p[1] = bar1pos;  t3.p[2] = bar2pos;

                t0.b[0] = t.b[0];   t0.b[1] = false;    t0.b[2] = t.b[2];
                t1.b[0] = t.b[0];   t1.b[1] = t.b[1];   t1.b[2] = false;
                t2.b[0] = false;    t2.b[1] = t.b[1];   t2.b[2] = t.b[2];
                t3.b[0] = false;    t3.b[1] = false;    t3.b[2] = false;

                store(t0,m_triangles);
                store(t1,m_triangles);
                store(t2,m_triangles);
                store(t3,m_triangles);

            }
            //Esto hay que revisarlo para una segunda adaptacion, puesto que ya hay polygonos
            if(num_split_fcs == 1)
            {
                auto nt = n_gon_with_hanging_node<4>(num_fcs,fcs_ids,t);
                store(nt,m_quadrangles);
            }
            if(num_split_fcs == 2)
            {
                auto nt = n_gon_with_hanging_node<5>(num_fcs,fcs_ids,t);
                store(nt,m_pentagons);
            }
    }

    void
    sort_uniq(std::vector<T>& v)
    {
        std::sort(v.begin(), v.end());
        auto uniq_iter = std::unique(v.begin(), v.end());
        v.erase(uniq_iter, v.end());
    }

    template<typename LoaderType>
    void re_populate_mesh( mesh_type & msh,
                            const std::vector<tensors<T>>& tsr_vec,
                            const T yield,
                            const int degree)
    {
        mesh_type re_msh, re_msh_2;


        re_msh = msh;

        auto storage    = msh.backend_storage();
        auto re_storage = re_msh.backend_storage();
        auto num_cells  = msh.cells_size();
        size_t nds_counter = 0;

        if(m_hanging_nodes)
            marker_with_hanging_nodes(re_msh);
        else
            marker(re_msh);

        dump_to_matlab(msh,"mesh.m",cl_marks_vector);

        for(auto& cl:msh)
        {
            auto cl_id   = msh.lookup(cl);
            auto cl_mark = cl_marks_vector.at(cl_id);

            if(cl_mark > 0)
            {
                if(m_hanging_nodes)
                    refine_single_with_hanging_nodes(msh,cl);
                else
                    refine_single(msh,cl);
            }
            else
            {
                auto pts_ids = cl.point_ids();
                auto fcs_ids = cl.faces_ids();
                auto num_fcs = fcs_ids.size();
                n_gon<3> t;
                for(size_t i = 0; i < num_fcs; i++)
                {
                    auto pts_ids = cl.point_ids();
                    auto id  = fcs_ids.at(i);
                    t.p[i]   = pts_ids.at(i); // Since later in populate_mesh values are subtract  by 1
                    t.b[i]   = storage->boundary_edges.at(id);
                }
                store(t,m_triangles);
            }
            //else
        }

        std::sort(m_edges.begin(), m_edges.end());
        auto uniq_iter = std::unique(m_edges.begin(), m_edges.end());
        m_edges.erase(uniq_iter, m_edges.end());

        LoaderType      loader;
        loader.m_edges        = m_edges;
        loader.m_triangles    = m_triangles;
        loader.m_quadrangles  = m_quadrangles;
        loader.m_pentagons    = m_pentagons;
        loader.m_hexagons     = m_hexagons;
        loader.m_boundary_edges = m_boundary_edges;
        //loader.m_ennagons  = m_ennagons;

        loader.m_points     = re_storage->points;

        auto storage_rm2  = re_msh_2.backend_storage();

        loader.populate_mesh(re_msh_2);
        msh = re_msh_2;

    }

};


}//disk
