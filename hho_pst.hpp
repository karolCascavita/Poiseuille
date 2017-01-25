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

template<typename T,size_t DIM>
std::string
directory(const T Bi, const std::string& name)
{
    auto Bi_str  = to_string(int(10* Bi));
    auto DIM_str = to_string(DIM);

    std::string dir_name    =   name + "/" + DIM_str + "D_Bi" + Bi_str;
    std::cout << "dir_name  = "<< dir_name << std::endl;
    return dir_name;
}

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
    size_t quad_degree;

    tensors() //quad_degree(2*degree+2)
    {}

    void
    set_quad_degree(const size_t degree)
    {
    }

    template< typename MeshType>
    void
    zero_tensor(const MeshType& msh, const typename MeshType::cell& cl, const size_t degree)
    {
        typedef MeshType                                    mesh_type;
        typedef typename MeshType::cell                     cell_type;
        typedef disk::quadrature<mesh_type, cell_type>      cell_quad_type;

        //WK: try to do this in the constructor.(I had problems in zero_tensor_vector, since it s asking for the input of the constructor)
        quad_degree = 2*degree;

        auto DIM = MeshType::dimension;
        auto cq  = cell_quad_type(quad_degree);
        auto cqs = cq.integrate(msh, cl).size();
        auto p   = points(msh,cl);
        auto ps  = p.size();
        siglam   = tensor_matrix::Zero(DIM, cqs + ps );
        gamma    = tensor_matrix::Zero(DIM, cqs + ps );
        xi_norm  = tensor_matrix::Zero(1  , cqs + ps );
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
        tsr.zero_tensor(msh, cl, degree);
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
        //WK: for p-adaptation quad_degree is not the same for all cells
        size_t quad_degree = tsr.quad_degree;

        cell_quadrature     = cell_quadrature_type(quad_degree);
        face_basis          = face_basis_type(m_degree);

        auto fcs = faces(msh, cl);
        auto num_faces      = fcs.size();
        auto num_cell_dofs  = cell_basis.range(0,m_degree).size();
        auto num_face_dofs  = face_basis.size();
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

            vector_type  xi   =   tsr.siglam.col(cont)  +  m_alpha * dphi_rec_uh;
            T xi_norm =  xi.norm();
            tsr.xi_norm(cont) = xi_norm;

            if(xi_norm > m_yield)
            {
                vector_type  gamma;
                gamma = m_method_coef * (xi_norm - m_yield)* (xi / xi_norm);
                tsr.siglam.col(cont) += m_alpha * ( dphi_rec_uh - gamma );
                tsr.gamma.col(cont)   = gamma;

                rhs  += qp.weight() * dphi_rec.transpose() * (tsr.siglam.col(cont) - m_alpha * gamma) ;
            }
            else
            {
                tsr.siglam.col(cont) += m_alpha * dphi_rec_uh;
                tsr.gamma.col(cont)   = vector_type::Zero(mesh_type::dimension);

                rhs +=  qp.weight() * dphi_rec.transpose() * tsr.siglam.col(cont);
            }
            ++cont;
        }

        auto pts = points(msh,cl);
        auto cqs = cell_quadrature.integrate(msh, cl).size();

        for(size_t i = 0, j = cqs; i < pts.size(); i++, j ++)
        {
            auto dphi       =   cell_basis.eval_gradients(msh, cl, pts.at(i));
            auto col_range  =   cell_basis.range(1,m_degree+1);
            auto row_range  =   disk::dof_range(0,mesh_type::dimension);

            matrix_type dphi_matrix =   make_gradient_matrix(dphi);
            matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
            matrix_type dphi_rec    =   dphi_taken * rec_oper;

            vector_type dphi_rec_uh =   dphi_rec*uh_TF;

            vector_type  xi   =   tsr.siglam.col(j)  +  m_alpha * dphi_rec_uh;
            T xi_norm =  xi.norm();
            tsr.xi_norm(j) = xi_norm;

            if(xi_norm > m_yield)
            {
                vector_type  gamma;
                gamma = m_method_coef * (xi_norm - m_yield)* (xi / xi_norm);
                tsr.siglam.col(j) += m_alpha * ( dphi_rec_uh - gamma );
                tsr.gamma.col(j)   = gamma;
            }
            else
            {
                tsr.siglam.col(j) += m_alpha * dphi_rec_uh;
                tsr.gamma.col(j)   = vector_type::Zero(mesh_type::dimension);
            }
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

    std::vector<size_t> cl_marks_vector;


    stress_based_mesh(const mesh_type& msh,const std::vector<tensors<T>>& tsr_vec, const T yield, bool& do_refinement)
    {
        std::vector<size_t> vec(msh.cells_size(),0);

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto mat = tsr_vec.at(i).xi_norm;
            std::cout << "mat = [" << std::endl;

            for(size_t m = 0; m < mat.size(); m++)
            {
                std::cout << "  "<< m;

                if(std::isnan(mat(0,m)))
                    throw std::logic_error("The norm of the constrain is NaN");
            }
            std::cout<<"]" << std::endl;

            if(std::abs(mat.minCoeff() - yield) < 1.e-10 || std::abs(mat.minCoeff()) < yield )
            {
                if(std::abs(mat.maxCoeff() - yield) > 1.e-10 && std::abs(mat.maxCoeff()) > yield )
                    vec.at(i) = 1;
            }
        }

        cl_marks_vector = vec;
        do_refinement = false;
        for(size_t i = 0; i < cl_marks_vector.size(); i++)
        {
            std::cout<<"cl_marks_vector["<<i<<"]"<<cl_marks_vector[i]<<std::endl;
            if(cl_marks_vector[i] == 1)
            {
                do_refinement = true;
                break;
            }
        }
        std::cout << do_refinement << std::endl;

    }
    void
    re_populate_mesh(mesh_type & msh,
                            const std::vector<tensors<T>>& tsr_vec,
                            const T yield,
                            const std::string&  directory,
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
            auto lmkr  = (id != 0)?  cl_marks_vector.at(id-1) : 0;
            auto rmkr  = (id != msh.cells_size()-1)?  cl_marks_vector.at(id+1) : 0;

            auto cell_pts   = points(msh, cl);
            storage_rm->points.push_back(cell_pts.at(0));

            if(cmkr == 1 || (lmkr == 1 || rmkr == 1))
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
    typedef typename mesh_type::face                face_type;
    typedef typename mesh_type::cell_iterator       cell_itor;

    bool m_hanging_nodes;
    std::vector<size_t> cl_marks_vector;
    //std::vector<std::pair<bool,typename point<T,2>::id_type>> fc_marks_vector;
    std::vector<std::pair<bool, int>> fc_marks_vector;
    std::vector<std::array<ident_impl_t, 4>>        m_edges;
    std::vector<std::array<ident_impl_t, 2>>        m_boundary_edges;

    std::vector<std::array<ident_impl_t, 3>>        m_triangles;
    std::vector<std::array<ident_impl_t, 4>>        m_quadrangles;
    std::vector<std::array<ident_impl_t, 5>>        m_pentagons;
    std::vector<std::array<ident_impl_t, 6>>        m_hexagons;
    std::vector<std::array<ident_impl_t, 7>>        m_heptagons;
    std::vector<std::array<ident_impl_t, 8>>        m_octagons;
    std::vector<std::array<ident_impl_t, 9>>        m_enneagons;
    std::vector<std::array<ident_impl_t,10>>        m_decagons;
    std::vector<std::array<ident_impl_t,11>>        m_hendecagons;
    std::vector<std::array<ident_impl_t,12>>        m_dodecagons;

    stress_based_mesh(const mesh_type& msh,
                        const std::vector<tensors<T>>& tsr_vec,
                        const T yield,
                        bool& do_refinement,
                        const bool hanging_nodes):m_hanging_nodes(hanging_nodes)
    {
        cl_marks_vector = cells_marker(msh,tsr_vec,yield);
        for(size_t cm = 0; cm < cl_marks_vector.size(); cm++)
        {
            size_t b = cl_marks_vector.at(cm);
            std::cout<<b<<"  ";
        }
        std::cout<<std::endl;

        fc_marks_vector.resize(msh.faces_size());
        for(size_t i = 0; i < fc_marks_vector.size(); i++)
            fc_marks_vector.at(i).first = false;

        do_refinement = false;
        for(size_t i = 0; i < cl_marks_vector.size(); i++)
        {
            if(cl_marks_vector[i] > 0)
            {
                do_refinement = true;
                break;
            }
        }
    }
    template< size_t N>
    struct n_gon
    {
       std::array<size_t, N>  p;
       std::array<bool, N>    b;
    };

    std::vector<size_t>
    cells_marker(const mesh_type& msh,
                        const std::vector<tensors<T>>& tsr_vec,
                        const T yield)
    {
        std::vector<size_t> vec(msh.cells_size());

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto mat = tsr_vec.at(i).xi_norm;

            for(size_t m = 0; m < mat.size(); m++)
            {
                if(std::isnan(mat(0,m)))
                    throw std::logic_error("The norm of the constrain is NaN");
            }

            if(std::abs(mat.minCoeff() - yield) < 1.e-10 || std::abs(mat.minCoeff()) < yield )
            {
                if(std::abs(mat.maxCoeff() - yield) > 1.e-10 && std::abs(mat.maxCoeff()) > yield )
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
             auto a = t.p.at(i);
             auto b = t.p.at(i+1);
             if(b < a)
                 std::swap(a,b);
             m_edges.push_back({a, b});
             if(t.b.at(i))
                 m_boundary_edges.push_back({a , b});
         }
         auto a = t.p.at(0);
         auto b = t.p.at(N-1);
        if(b < a)
            std::swap(a,b);
        m_edges.push_back({a , b});
        if(t.b[N -1])
            m_boundary_edges.push_back({a , b });
    }

     #if 0
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
    #endif
    void
    face_mark_hanging_nodes(mesh_type & msh, const face_type & fc)
    {
        auto storage = msh.backend_storage();
        auto  fc_id  = msh.lookup(fc);
        fc_marks_vector.at(fc_id).first  = true;
        auto bar = barycenter(msh,fc);
        storage ->points.push_back(bar);
        fc_marks_vector.at(fc_id).second = storage->points.size() -1;
        //std::cout << "inside marker, marking bundary = "<<fc_id << std::endl;
    }

    void
    set_marks_hanging_nodes(mesh_type & msh, const cell_type & cl)
    {
        auto cl_id   = msh.lookup(cl);
        auto cl_mark = cl_marks_vector.at(cl_id);
        //To adapt also the neighbors take into account 3 and 2, otherwise just 3
        if(cl_mark == 3)
        //if(cl_mark > 1) /
        {
            auto fcs     = faces(msh,cl);
            for(auto & fc: fcs)
            {
                if(msh.is_boundary(fc))
                    face_mark_hanging_nodes(msh, fc);
                else
                {
                    auto  fc_id    = msh.lookup(fc);
                    if(!fc_marks_vector.at(fc_id).first)
                    {
                        face_mark_hanging_nodes(msh, fc);
                        std::cout << "cell = "<< cl_id <<"    ; fc = "<< fc_id << std::endl;

                        size_t neighbor_id = face_owner_cells_ids(msh, fc, cl);

                        std::cout << "   ngh_id_1 = "<< neighbor_id<<std::endl;

                        typename cell_type::id_type id(neighbor_id);

                        std::cout << "   ngh_id_2 = "<< neighbor_id<<std::endl;
                        std::cout << "   id     = "<< id << std::endl;

                        auto  ngh_mark      = cl_marks_vector.at(id);

                        if(cl_mark > ngh_mark)
                        {
                            cl_marks_vector.at(neighbor_id) = cl_mark - 1;
                            auto  ngh  = *(msh.cells_begin() + size_t(neighbor_id));
                            set_marks_hanging_nodes(msh, ngh);
                        }
                    }
                }
            }
        }
        return;
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
        std::cout << "/ *******  review old mesh  ********/" << std::endl;
        borrar(msh);

        for(auto& cl : msh)
        {
            auto cl_id   = msh.lookup(cl);
            bool cl_mark = cl_marks_vector.at(cl_id);
            if(cl_mark)
            {
                if(m_hanging_nodes)
                    set_marks_hanging_nodes(msh, cl);
                else
                    set_marks(msh, cl);
            }
        }
    }

    template<size_t N, size_t M, typename FaceIds>
    n_gon<N>
    n_gon_base(const n_gon<M>& t, const FaceIds& fcs_ids)
    {
       n_gon<N> nt;
       auto num_fcs = fcs_ids.size();
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

   template<size_t N>
   n_gon<N>
   put_n_gon(const mesh_type & msh, const cell_type& cl)
   {
      n_gon<N> t;
      auto storage = msh.backend_storage();
      auto pts_ids = cl.point_ids();
      auto fcs_ids = cl.faces_ids();
      for(size_t i = 0; i < N; i++)
      {
          auto id  = fcs_ids.at(i);
          t.p[i]   = pts_ids.at(i); // Since later in populate_mesh values are subtract  by 1
          t.b[i]   = storage->boundary_edges.at(id);
      }
      return t;
   }

   template<size_t N>
   n_gon<N>
   put_n_gon(const mesh_type & msh, const std::vector<ident_impl_t>& fcs_ids,
                                    const std::vector<ident_impl_t>& pts_ids)
   {
      n_gon<N> t;
      auto storage = msh.backend_storage();
      for(size_t i = 0; i < N; i++)
      {
          auto id  = fcs_ids.at(i);
          t.p[i]   = pts_ids.at(i); // Since later in populate_mesh values are subtract  by 1
          t.b[i]   = storage->boundary_edges.at(id);
      }
      return t;
   }

   template<size_t N>
   n_gon<N>
   put_n_gon_with_hanging_node(mesh_type & msh, const cell_type& cl)
   {

       auto fcs_ids = cl.faces_ids();
       auto num_fcs = fcs_ids.size();

        if(num_fcs == 3)
        {
           auto t = put_n_gon<3>(msh, cl);
           return  n_gon_base<N>(t, fcs_ids);
        }
        if(num_fcs == 4)
        {
            auto t = put_n_gon<4>(msh, cl);
            return  n_gon_base<N>(t, fcs_ids);
        }

        if(num_fcs == 5)
        {
            auto t = put_n_gon<5>(msh, cl);
            return  n_gon_base<N, 5>(t, fcs_ids);
        }
        if(num_fcs == 6)
        {
            auto t = put_n_gon<6>(msh, cl);
            return  n_gon_base<N, 6>(t, fcs_ids);
        }
        if(num_fcs == 7)
        {
            auto t = put_n_gon<7>(msh, cl);
            return  n_gon_base<N, 7>(t, fcs_ids);
        }
        if(num_fcs == 8)
        {
            auto t = put_n_gon<8>(msh, cl);
            return  n_gon_base<N, 8>(t, fcs_ids);
        }
        if(num_fcs == 9)
        {
            auto t = put_n_gon<9>(msh, cl);
            return  n_gon_base<N, 9>(t, fcs_ids);
        }
        if(num_fcs == 10)
        {
            auto t = put_n_gon<10>(msh, cl);
            return  n_gon_base<N, 10>(t, fcs_ids);
        }
        if(num_fcs == 11)
        {
            auto t = put_n_gon<11>(msh, cl);
            return  n_gon_base<N, 11>(t, fcs_ids);
        }
        if(num_fcs == 12)
        {
            auto t = put_n_gon<12>(msh, cl);
            return  n_gon_base<N, 12>(t, fcs_ids);
        }
        throw std::logic_error("This shouldn't come to this point");
    }

    void refine_single(const mesh_type& msh, const cell_type& cl)
    {
        std::cout << "WARNING INSIDE REFINE_SINGLE!!" << std::endl;
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


    void
    refine_reference_triangle(const n_gon<3> & t , const std::array<int, 3> & bar_ids)
    {

        n_gon<3> t0,t1,t2,t3;
        auto bar0pos = bar_ids.at(0);
        auto bar1pos = bar_ids.at(1);
        auto bar2pos = bar_ids.at(2);

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
    void refine_triangle(mesh_type& msh, const cell_type& cl)
    {
        auto fcs_ids = cl.faces_ids();
        auto bar0pos = fc_marks_vector.at(fcs_ids[0]).second;
        auto bar1pos = fc_marks_vector.at(fcs_ids[1]).second;
        auto bar2pos = fc_marks_vector.at(fcs_ids[2]).second;
        std::array<int, 3> bar_ids = {bar0pos, bar1pos, bar2pos};

        auto t  = put_n_gon<3>(msh, cl);

        refine_reference_triangle( t,  bar_ids);
    }

    void
    refine_other(mesh_type& msh, const cell_type& cl)
    {

        //if( hanging nodes)
        auto pts      = points(msh,cl);
        auto storage  = msh.backend_storage();
        auto pts_ids  = cl.point_ids();
        auto fcs_ids  = cl.faces_ids();
        auto num_fcs  = fcs_ids.size();
        auto fbar_ids = std::vector<size_t>(num_fcs);
        auto cbar     = barycenter(msh,cl);
        storage ->points.push_back(cbar);
        auto cbar_id = storage ->points.size() - 1;
        for (size_t i = 0; i < num_fcs; i++)
        {
           auto fbar = ( pts.at(i) + cbar ) / 2.;
           storage ->points.push_back(fbar);
           fbar_ids.at(i)  = storage ->points.size() - 1;
        }
        for(size_t i = 0; i < num_fcs; i++)
        {
            n_gon<3>   nt;
            std::array<int, 3> tbar_ids;
            nt.p[0] = pts_ids[i];
            nt.p[1] = pts_ids[(i+1)%num_fcs];
            nt.p[2] = cbar_id;

            nt.b[0] = storage->boundary_edges.at(fcs_ids.at(i));
            nt.b[1] = false;
            nt.b[2] = false;

            tbar_ids.at(0) = fc_marks_vector.at(fcs_ids.at(i)).second;
            tbar_ids.at(1) = fbar_ids.at((i+1)%num_fcs);
            tbar_ids.at(2) = fbar_ids.at(i);

            refine_reference_triangle(nt , tbar_ids);
        }
    }

    template<typename PtA>
    std::pair<bool,std::vector<point<T,2>>>
    is_special_polygon(const mesh_type& msh, const cell_type& cl,
                        const PtA& pts)
    {
         auto fcs = faces(msh, cl);
         std::vector<std::array<T,2>> ns(fcs.size());
         size_t i = 0;
         for(auto& fc : fcs)
         {
             auto n = normal(msh,cl,fc);
             ns.at(i)[0] = n(0);
             ns.at(i)[1] = n(1);
             i++;
         }
        std::sort(ns.begin(), ns.end());
        auto uniq_iter = std::unique(ns.begin(), ns.end(),[](std::array<T,2>& l, std::array<T,2>& r)
            {return std::sqrt(std::pow((l[0] - r[0]),2.) + std::pow((l[1] - r[1]),2.)) < 1.e-10; });
        ns.erase(uniq_iter, ns.end());

         //Identify vertices
         std::vector<point<T,2>> vertices(ns.size());
         auto num_pts = pts.size();
         size_t vcount = 0;

         for(size_t i = 0; i < num_pts; i++)
         {
             size_t idx = (i == 0)? num_pts - 1 : i - 1 ;
             auto pb = pts.at(idx);
             auto p  = pts.at(i);
             auto pf = pts.at((i + 1) % num_pts);

             auto u  = (pb - p).to_vector();
             auto v  = (pf - p).to_vector();
             auto uxv_norm = cross(u, v).norm();

             if(uxv_norm > 1.e-10)
             {
                 vertices.at(vcount) = p;
                 vcount++;
             }
         }
        if(vcount != ns.size())
             std::logic_error(" Incorrect procedure to find vertices");

        bool has_hang_nodes(false);
        if(vertices.size() != pts.size())
            has_hang_nodes = true;

        auto ret = std::make_pair(has_hang_nodes, vertices);

        return ret;
    }

    void
    refine_special_triangle(mesh_type& msh, const cell_type& cl,
                                        const std::vector<point<T,2>>& vertices)
    {
        bool do_bar;
        int i(0), j(0), m(0), pos(0), npb(0);

        auto pts      = points(msh,cl);
        auto pts_ids  = cl.point_ids();
        auto fcs_ids  = cl.faces_ids();
        auto storage  = msh.backend_storage();
        std::vector<std::array<int,5>> owner_cells;
        owner_cells.reserve(pts_ids.size());
        point<T,2> bar;
        for(auto & v:vertices)
        {
            auto w   = vertices.at((i + 1)% 3);
            bar = (v + w)/2.;
            auto b_distance = (((v - w)/2).to_vector()).norm();
            std::cout << "b_distance" << b_distance<< std::endl;
            if( b_distance < 1.e-10)
                throw std::logic_error("Flat triangle");

            pos = 0;
            do_bar = true;

            for(size_t ip = m; ip < pts.size();ip++, j++)
            {
                auto p = pts.at(ip);
                auto vp_distance = ( p - v).to_vector().norm();
                auto dist = (vp_distance - b_distance) / b_distance;

                auto fc_id       = fcs_ids.at(ip);

                if( std::abs(dist) < 1.e-10)
                {
                    owner_cells.push_back({ j, i , (i + 1) % 3 , int(pts_ids.at(ip)), int(fc_id)});
                    do_bar = false;
                }
                else
                {
                    if(dist < 0.)
                    {
                        owner_cells.push_back({ j, i , -1, int(pts_ids.at(ip)), int(fc_id)});
                        if(std::abs(dist + b_distance) > 1.e-10)
                            pos++;
                    }
                    if(dist > 0.)
                    {
                        owner_cells.push_back({ j, (i + 1) % 3, -1, int(pts_ids.at(ip)), int(fc_id)});
                        if(std::abs(dist - 1.0) < 1.e-10)
                        {
                            if(do_bar)
                            {
                                storage ->points.push_back(bar);
                                auto bar_id = storage ->points.size() - 1;
                                for(size_t ii = npb + pos + 1; ii < owner_cells.size(); ii ++)
                                    ++owner_cells.at(ii)[0];

                                owner_cells.push_back({ int(npb + pos + 1), i, (i + 1) % 3, int(bar_id), -1});
                                j++;
                            }
                            m = ip + 1;
                            npb = j;
                            j++;
                            break;
                        }
                    }
                }
            }
            i++;
        }

        //for the last face
        if(do_bar)
        {
            storage ->points.push_back(bar);
            auto bar_id = storage ->points.size() - 1;
            i = i - 1;
            for(size_t ii = npb + pos + 1; ii < owner_cells.size(); ii ++)
              ++owner_cells.at(ii)[0];
            owner_cells.push_back({int(npb + pos + 1), i, (i + 1) % 3, int(bar_id),-1});
        }
        std::sort(owner_cells.begin(), owner_cells.end());

        std::cout << " cell_points = " << std::endl;
        for(auto& v : pts)
            std::cout<<"             "<<  v <<std::endl;

        std::cout << " owner_cells = " << std::endl;
        for(auto& v : owner_cells)
            std::cout<<"            "<< v[0] <<" " <<v[1]<<" " <<v[2]<<" " <<v[3]<<" "<<v[4]<<std::endl;

        std::cout << " total_points = " << std::endl;
        for(auto& v : owner_cells)
        {
            auto id = v[3];
            std::cout<<"    point( "<< id <<") = "<<  storage ->points.at(id) <<std::endl;
        }



        size_t ii = 0;
        n_gon<3> t3;

        for (auto i = owner_cells.begin(), toofar = owner_cells.end(); i != toofar; ++i)
        {
                auto ptr    = (*i).begin();
                auto is_barycenter = (*(ptr + 2) != -1)?  true : false;
                if(is_barycenter)
                {
                    auto  bar_id   =  *(ptr + 3);
                    t3.p.at(ii) = bar_id;
                    t3.b.at(ii) = false;
                    ii++;
                    if(*(ptr + 4) == -1)
                    {
                        auto it = (*(i-1)).begin() + 4;
                        *(ptr + 4) = *it;
                    }
                }
        }
        store(t3,m_triangles);
        std::cout << " owner_cells = " << std::endl;
        for(auto& v : owner_cells)
            std::cout<<"            "<< v[0] <<" " <<v[1]<<" " <<v[2]<<" " <<v[3]<<" "<<v[4]<<std::endl;

        for(size_t j = 0; j < 3; j++)
        {
            std::vector<size_t> pts_ids;
            std::vector<size_t> fcs_ids;

            pts_ids.reserve(3);
            fcs_ids.reserve(3);

            for (auto i = owner_cells.begin(), toofar = owner_cells.end(); i != toofar; ++i)
            {
                auto ptr    = (*i).begin();
                auto in_triangle = (*(ptr + 1) == j) | (*(ptr + 2) == j)?  true : false;
                if(in_triangle)
                {
                    auto  pt_id = *(ptr + 3);
                    auto  fc_id = *(ptr + 4);
                    pts_ids.push_back(pt_id);
                    fcs_ids.push_back(fc_id);

                }
            }

            auto new_num_fcs = fcs_ids.size();
            if(new_num_fcs == 3)
            {
                std::cout << "store in 3angles" << std::endl;
                auto nt = put_n_gon<3>(msh, fcs_ids, pts_ids);
                store(nt,m_triangles);
            }
            if(new_num_fcs == 4)
            {
                std::cout << "store in 4angles" << std::endl;
                auto nt = put_n_gon<4>(msh, fcs_ids, pts_ids);
                store(nt,m_quadrangles);
            }
            if(new_num_fcs == 5)
            {
                std::cout << "store in 5angles" << std::endl;
                auto nt = put_n_gon<5>(msh, fcs_ids, pts_ids);
                store(nt,m_pentagons);
            }
            if(new_num_fcs == 6)
            {
                std::cout << "store in 6angles" << std::endl;

                auto nt = put_n_gon<6>(msh, fcs_ids, pts_ids);
                store(nt,m_hexagons);
            }
            if(new_num_fcs == 7)
            {
                std::cout << "store in 7angles" << std::endl;
                auto nt = put_n_gon<7>(msh, fcs_ids, pts_ids);
                store(nt,m_heptagons);
            }
            if(new_num_fcs == 8)
            {
                std::cout << "store in 8angles" << std::endl;
                auto nt = put_n_gon<8>(msh, fcs_ids, pts_ids);
                store(nt,m_octagons);
            }
            if(new_num_fcs == 9)
            {
                std::cout << "store in 9angles" << std::endl;
                auto nt = put_n_gon<9>(msh, fcs_ids, pts_ids);
                store(nt,m_enneagons);
            }
            if(new_num_fcs == 10)
            {
                std::cout << "store in 10angles" << std::endl;
                auto nt = put_n_gon<10>(msh, fcs_ids, pts_ids);
                store(nt,m_decagons);
            }
            if(new_num_fcs == 11)
            {
                std::cout << "store in 11angles" << std::endl;

                auto nt = put_n_gon<11>(msh, fcs_ids, pts_ids);
                store(nt,m_hendecagons);
            }
            if(new_num_fcs == 12)
            {
                std::cout << "store in 12angles" << std::endl;

                auto nt = put_n_gon<12>(msh, fcs_ids, pts_ids);
                store(nt,m_dodecagons);
            }

            if(new_num_fcs > 12)
            {
                std::cout << "number of faces = "<< new_num_fcs << std::endl;
                throw std::logic_error("number of faces exceeds maximum. Add new array to store it.");
            }
        }
    }

    void
    refine_single_with_hanging_nodes(mesh_type& msh, const cell_type& cl, const size_t cl_mark)
    {
        auto fcs_ids = cl.faces_ids();
        auto num_fcs = fcs_ids.size();
        auto pts     = points(msh, cl);
        /* Neighbors Adaptation:
        To adapt marked cells(3) and/or their neighbors(2)
        cells and neihgbors => if(cl_mark > 1)
        Only cells          => if(cl_mark == 3)*/

        if(cl_mark == 3)
        {
            std::cout << " cell = "<< msh.lookup(cl) << std::endl;

            auto  e   = is_special_polygon(msh, cl, pts);

            if(e.first)
            {
                std::vector<point<T,2>>  vertices = e.second;
                std::cout << "cell with hanging_nodes ("<< pts.size() <<") v("<< vertices.size()<< ")" << std::endl;

                switch (vertices.size())
                {
                    case 1:
                    case 2:
                        throw std::logic_error("Number of faces cannot be less than 3");
                        break;
                    case 3:
                        std::cout << "refine special triangle ("<<pts.size() <<")" << std::endl;
                        std::cout << "vertices = " << std::endl;
                        for(auto& v: vertices)
                            std::cout <<"        "<< v << std::endl;

                        refine_special_triangle(msh, cl, vertices);
                        break;
                    //case 4:
                    //    refine_quadrangle();
                    //    break;
                    default:
                        std::cout << "refine other" << std::endl;
                        refine_other(msh ,cl); //WK: try to do the same for quadrangles
                        break;
                }
            }
            else
            {
                switch (num_fcs)
                {
                    std::cout << "regular cells" << std::endl;
                    case 1:
                    case 2:
                        throw std::logic_error("Number of faces cannot be less than 3");
                        break;
                    case 3:
                        std::cout << "refine triangle" << std::endl;
                        refine_triangle(msh, cl);
                        break;
                    //case 4:
                    //    refine_quadrangle();
                    //    break;
                    default:
                        std::cout << "refine other" << std::endl;
                        refine_other(msh ,cl);
                        break;
                }

            }
        }
        else
        {
            auto count   = 0;
            /* Neighbors Adaptation:
            hanging_nodes on neighbors(2)               => if( cl_mark == 2)
            hanging_nodes on neihgbors' neihgbors (1)   => if( cl_mark == 1)*/
            if( cl_mark == 2)
            {
                for(size_t i = 0; i < num_fcs; i++)
                {
                auto id  = fcs_ids.at(i);
                auto fc_mark = fc_marks_vector.at(id).first;
                if(fc_mark)
                    count++;
                }
            }
            auto new_num_fcs = count + num_fcs;
            if(new_num_fcs == 3)
            {
                auto nt = put_n_gon_with_hanging_node<3>(msh, cl);
                store(nt,m_triangles);
            }
            if(new_num_fcs == 4)
            {
                auto nt = put_n_gon_with_hanging_node<4>(msh, cl);
                store(nt,m_quadrangles);
            }
            if(new_num_fcs == 5)
            {
                auto nt = put_n_gon_with_hanging_node<5>(msh, cl);
                store(nt,m_pentagons);
            }
            if(new_num_fcs == 6)
            {

                auto nt = put_n_gon_with_hanging_node<6>(msh, cl);
                store(nt,m_hexagons);
            }
            if(new_num_fcs == 7)
            {
                auto nt = put_n_gon_with_hanging_node<7>(msh, cl);
                store(nt,m_heptagons);
            }
            if(new_num_fcs == 8)
            {
                auto nt = put_n_gon_with_hanging_node<8>(msh, cl);
                store(nt,m_octagons);
            }
            if(new_num_fcs == 9)
            {
                auto nt = put_n_gon_with_hanging_node<9>(msh, cl);
                store(nt,m_enneagons);
            }
            if(new_num_fcs == 10)
            {
                auto nt = put_n_gon_with_hanging_node<10>(msh, cl);
                store(nt,m_decagons);
            }
            if(new_num_fcs == 11)
            {
                auto nt = put_n_gon_with_hanging_node<11>(msh, cl);
                store(nt,m_hendecagons);
            }
            if(new_num_fcs == 12)
            {
                auto nt = put_n_gon_with_hanging_node<12>(msh, cl);
                store(nt,m_dodecagons);
            }

            if(new_num_fcs > 12)
            {
                std::cout << "number of faces = "<< new_num_fcs << std::endl;
                throw std::logic_error("number of faces exceeds maximum. Add new array to store it.");
            }
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
                            const std::string & directory,
                            const std::string & other_info,
                            const int degree)
    {
        mesh_type re_msh, re_msh_2;

        re_msh = msh;

        auto storage    = msh.backend_storage();
        auto re_storage = re_msh.backend_storage();
        auto num_cells  = msh.cells_size();
        size_t nds_counter = 0;

        marker(re_msh);
        dump_to_matlab(msh, directory + "/mesh" + other_info + ".m",cl_marks_vector);

        for(auto& cl : msh)
        {
            auto cl_id   = msh.lookup(cl);
            auto cl_mark = cl_marks_vector.at(cl_id);

            if(cl_mark > 0)
            {
                if(m_hanging_nodes)
                    refine_single_with_hanging_nodes(msh, cl, cl_mark);
                else
                    refine_single(msh,cl);
            }
            else
            {
                auto fcs_ids  =  faces(msh,cl); //WK: There should be direct way to know the number of faces of the cell
                int  num_fcs  =  fcs_ids.size();
                switch(num_fcs)
                {
                    case 3:
                        store(put_n_gon<3>(msh,cl), m_triangles);
                        break;
                    case 4:
                        store(put_n_gon<4>(msh,cl), m_quadrangles);
                        break;
                    case 5:
                        store(put_n_gon<5>(msh,cl), m_pentagons);
                        break;
                    case 6:
                        store(put_n_gon<6>(msh,cl), m_hexagons);
                        break;
                    case 7:
                        store(put_n_gon<7>(msh,cl), m_heptagons);
                        break;
                    case 8:
                        store(put_n_gon<8>(msh,cl), m_octagons);
                        break;
                    case 9:
                        store(put_n_gon<9>(msh,cl), m_enneagons);
                        break;
                    case 10:
                        store(put_n_gon<10>(msh,cl), m_decagons);
                        break;
                    case 11:
                        store(put_n_gon<11>(msh,cl), m_hendecagons);
                        break;
                    case 12:
                        store(put_n_gon<12>(msh,cl), m_dodecagons);
                        break;

                    default:
                        std::cout << "number of faces = "<< num_fcs << std::endl;
                        throw std::logic_error("Polygon not stored, number of faces exceeds maximum. Add an array to store it");
                        break;
                }
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
        loader.m_heptagons    = m_heptagons;
        loader.m_octagons     = m_octagons;
        loader.m_enneagons    = m_enneagons;
        loader.m_decagons     = m_decagons;
        loader.m_hendecagons  = m_hendecagons;
        loader.m_dodecagons   = m_dodecagons;

        loader.m_points       = re_storage->points;
        loader.m_boundary_edges = m_boundary_edges;

        auto storage_rm2  = re_msh_2.backend_storage();
        loader.populate_mesh(re_msh_2);
        msh = re_msh_2;

    }
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

    function f,sf;
    dfunction df;
    std::string     name;

    poiseuille(plasticity_data<T>& pst)
    {
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

    function        f,sf;
    dfunction       df;
    std::string     name;

    tuyau(plasticity_data<T>& pst)
    {
        name = "square";

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

    function        f,sf;
    dfunction       df;
    std::string     name;
    circular_tuyau(plasticity_data<T>& pst)
    {
        name = "circular";
        pst.f = 1.;
        T   R = 1.; // this should be read from the mesh

        f  = [pst](const point_type& p) -> T {
            T ret = 1.;
            return ret;
        };
        sf = [R,pst](const point_type& p) -> T {

            T  ret = 0;
            T    r = std::sqrt(p.x()*p.x() + p.y()*p.y());
            T   Bi = 2.*pst.yield/(pst.f*R);

            if(r/R >= Bi)
                0.5*(1. - (r/R)*(r/R)) - Bi*(1. - r/R);
            else
                0.5*(1. - Bi)*(1. - Bi);

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
};


}//disk
