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
        typedef typename MeshType::face                     face_type;
        typedef disk::quadrature<mesh_type, cell_type>      cell_quad_type;
        typedef disk::quadrature<mesh_type, face_type>      face_quad_type;

        //WK: try to do this in the constructor.(I had problems in zero_tensor_vector, since it s asking for the input of the constructor)
        quad_degree = 2*degree;

        auto DIM = MeshType::dimension;
        auto cq  = cell_quad_type(quad_degree);
        auto fq  = face_quad_type(quad_degree);

        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();
        auto fqs = num_faces*fq.integrate(msh, fcs[0]).size(); //WK: This should change if k differente for each face
        auto cqs = cq.integrate(msh, cl).size();

        auto p   = points(msh,cl);
        auto ps  = p.size();
        siglam   = tensor_matrix::Zero(DIM, cqs + fqs + ps );
        gamma    = tensor_matrix::Zero(DIM, cqs + fqs + ps );
        xi_norm  = tensor_matrix::Zero(1  , cqs + fqs + ps );
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
        cell_basis          = cell_basis_type(m_degree + 1); // This is k + 1  just if I want to use \nabla r^{k + 1}
        face_basis          = face_basis_type(m_degree);
        if(pst.method == true)
            m_method_coef = 1. / m_alpha ;
        else
            m_method_coef = 1. / (m_alpha + pst.mu);
    }

    template<int DIM>
    T
    make_prod_stress_n(const  dynamic_vector<T>& st, const static_vector<T,DIM>& n )
    {
        T p1;
        for(size_t i = 0; i < mesh_type::dimension; i++)
            p1  += st(i) * n(i,1);
        return p1;
    }
    T
    make_prod_stress_n(const  dynamic_vector<T>& st, const T n )
    {
        return  st(0) * n;
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
        face_quadrature     = face_quadrature_type(quad_degree);

        auto fcs = faces(msh, cl);
        auto num_faces      = fcs.size();
        auto cell_range     = cell_basis.range(0,m_degree);
        auto num_cell_dofs  = cell_range.size();
        auto num_face_dofs  = face_basis.size();
        dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);
        size_t cell_size = dsr.cell_range().size();

        rhs  = vector_type::Zero(dsr.total_size());

        size_t cont  =  0;
        auto   cqs   =  cell_quadrature.integrate(msh, cl);

        for (auto& qp : cqs)
        {
            auto dphi       =   cell_basis.eval_gradients(msh, cl, qp.point());
            auto col_range  =   cell_basis.range(1,m_degree+1);
            auto row_range  =   disk::dof_range(0,mesh_type::dimension);

            matrix_type dphi_matrix =   make_gradient_matrix(dphi);
            matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
            matrix_type dphi_rec    =   dphi_taken * rec_oper;
            vector_type dphi_rec_uh =   dphi_rec*uh_TF;
            matrix_type dphi_       =   take(dphi_matrix, cell_range, cell_range);

            vector_type  xi   =  tsr.siglam.col(cont)  +  m_alpha * dphi_rec_uh;
            vector_type  gamma;

            if(xi.norm() > m_yield)
                gamma   = m_method_coef * (xi.norm()- m_yield) * (xi / xi.norm());
            else
                gamma   = vector_type::Zero(mesh_type::dimension);

            tsr.siglam.col(cont) +=  m_alpha * ( dphi_rec_uh - gamma );
            tsr.gamma.col(cont)   =  gamma;
            tsr.xi_norm(cont)     =  xi.norm();

            rhs.head(cell_size)  += qp.weight() * dphi_.transpose() * (tsr.siglam.col(cont) - m_alpha * gamma) ;
            //rhs  += qp.weight() * dphi_rec.transpose() * (tsr.siglam.col(cont) - m_alpha * gamma) ;
            ++cont;
        }

        // Boundary term
        for (size_t face_i = 0, jj = 0; face_i < num_faces; face_i++)
        {
            auto current_face_range = dsr.face_range(face_i);
            auto fc = fcs.at(face_i);
            auto fqs= face_quadrature.integrate(msh, fc);

            auto n = normal(msh, cl, fc);

            //matrix_type m1 = matrix_type::Zeros(BG_row_range.size(), BG_col_range.size());
            for (auto& qp : fqs)
            {
                auto c_phi  = cell_basis.eval_functions(msh, cl, qp.point());
                auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point());
                auto f_phi  = face_basis.eval_functions(msh, fc, qp.point());
                auto col_range  =   cell_basis.range(1,m_degree+1);
                auto row_range  =   disk::dof_range(0,mesh_type::dimension);

                matrix_type c_dphi_matrix =   make_gradient_matrix(c_dphi);
                matrix_type c_dphi_taken  =   take(c_dphi_matrix, row_range, col_range);
                matrix_type c_dphi_rec    =   c_dphi_taken * rec_oper;
                vector_type c_dphi_rec_uh =   c_dphi_rec*uh_TF;

                vector_type  xi  =   tsr.siglam.col(cont)  +  m_alpha * c_dphi_rec_uh;
                vector_type  gamma;
                scalar_type  p1;

                if(xi.norm() > m_yield)
                    gamma  = m_method_coef * (xi.norm()- m_yield)* (xi / xi.norm());
                else
                    gamma  = vector_type::Zero(mesh_type::dimension);


                tsr.siglam.col(cont) += m_alpha * ( c_dphi_rec_uh - gamma );
                tsr.gamma.col(cont)   = gamma;
                tsr.xi_norm(cont) = xi.norm();

                // ((sigma - alpha*gamma) * n
                p1 = make_prod_stress_n(tsr.siglam.col(cont) - m_alpha * gamma , n);

                // ( (sigma - alpha*gamma) * n, vT )_F
                for (size_t i = cell_range.min(); i < cell_range.max(); i++)
                    rhs(i)   -= qp.weight() * mm_prod(c_phi.at(i), p1);

                // ( (sigma - alpha*gamma) * n, vF )_F
                for (size_t i = 0, ii = current_face_range.min();
                                      i < current_face_range.size(); i++, ii++)
                    rhs(ii)  += qp.weight() * mm_prod(f_phi.at(i), p1);

                ++cont;
            }
        }

        // tensor values in vertices
        auto pts = points(msh,cl);

        for(size_t i = 0, j = cont; i < pts.size(); i++, j ++)
        {
            auto dphi       =   cell_basis.eval_gradients(msh, cl, pts.at(i));
            auto col_range  =   cell_basis.range(1,m_degree+1);
            auto row_range  =   disk::dof_range(0,mesh_type::dimension);

            matrix_type dphi_matrix =   make_gradient_matrix(dphi);
            matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
            matrix_type dphi_rec    =   dphi_taken * rec_oper;
            vector_type dphi_rec_uh =   dphi_rec*uh_TF;

            vector_type xi   =   tsr.siglam.col(j)  +  m_alpha * dphi_rec_uh;
            vector_type gamma;

            if(xi.norm() > m_yield)
                gamma = m_method_coef * (xi.norm() - m_yield)* (xi / xi.norm());
            else
                gamma = vector_type::Zero(mesh_type::dimension);

            tsr.siglam.col(j) += m_alpha * ( dphi_rec_uh - gamma );
            tsr.gamma.col(j)   = gamma;
            tsr.xi_norm(j)     = xi.norm();

        }
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

    bool            is_exact;
    function        f,sf;
    dfunction       df;
    std::string     name;

    tuyau(plasticity_data<T>& pst)
    {
        is_exact = false;

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

    bool            is_exact;
    function        f,sf;
    dfunction       df;
    std::string     name;
    circular_tuyau(plasticity_data<T>& pst)
    {
        is_exact = true;

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
};


}//disk
