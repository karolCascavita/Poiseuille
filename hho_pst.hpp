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
#include "viscoplasticity.hpp"
namespace disk {

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
    nabla_ruh_at_point(const typename MeshType::point_type& pt,
                   const MeshType&  msh,
                   const CellType&  cl,
                   const CellBasisType& cell_basis,
                   const dynamic_vector<T>& ruh,
                   const size_t m_degree)
    {
        auto one_mindeg = 1;
        auto one_maxdeg = m_degree + 1;

        auto col_range = cell_basis.range( one_mindeg, one_maxdeg);
        auto row_range = disk::dof_range( 0, MeshType::dimension);
        auto dphi = cell_basis.eval_gradients(msh, cl, pt, one_mindeg, one_maxdeg);

        dynamic_vector<T> dphi_ruh    = dphi.transpose() * ruh;
        return dphi_ruh;
    }

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

    size_t                                  m_degree, m_degree_cell;

public:
    diffusion_like_static_condensation_pst()
        : m_degree(1), m_degree_cell(1)
    {
        cell_basis          = cell_basis_type(m_degree_cell);
        cell_quadrature     = cell_quadrature_type(2*m_degree_cell);
        face_quadrature     = face_quadrature_type(2*m_degree);
        face_basis          = face_basis_type(m_degree);
    }
    #if 0
    diffusion_like_static_condensation_pst(size_t degree)
        : m_degree(degree),  m_degree_cell(degree)
    {
        cell_basis          = cell_basis_type(m_degree_cell);
        cell_quadrature     = cell_quadrature_type(2*m_degree_cell);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
    }
    #endif
    diffusion_like_static_condensation_pst(size_t degree, const bool& change_degree)
        : m_degree(degree),  m_degree_cell(degree)
    {
        if(change_degree)
            m_degree_cell = m_degree + 1;

        cell_basis          = cell_basis_type(m_degree_cell);
        cell_quadrature     = cell_quadrature_type(2*m_degree_cell);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree);
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

        #if 0
        std::cout << "cell_rhs : "<< cell_rhs.rows()<<cell_rhs.cols() << std::endl;
        std::cout << "pst_rhs  : "<< plastic_rhs.rows()<<plastic_rhs.cols() << std::endl;
        std::cout << "bF_F  : "<< bP_F.rows()<< " x "<<bP_F.cols() << std::endl;
        std::cout << "bF_T  : "<< bP_T.rows()<< " x "<<bP_T.cols() << std::endl;
        #endif

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


        size_t num_dirichlet_edges(0);
        for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
        {
            auto bfc = *itor;
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), bfc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");
            auto face_id     = eid.second;
            auto binfo = msh.boundary_information(bfc);

            if(binfo.bndtype == boundary::DIRICHLET)
                num_dirichlet_edges++;
        }
        m_num_unknowns = face_basis.size() * (msh.faces_size() + num_dirichlet_edges);
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

            auto binfo = msh.boundary_information(bfc);

            if(binfo.bndtype == boundary::DIRICHLET)
            {

                auto face_offset = face_id * fbs;
                auto face_offset_lagrange = (msh.faces_size() + face_i) * fbs;

                auto fqd = face_quadrature.integrate(msh, bfc);

                matrix_type MFF     = matrix_type::Zero(fbs, fbs);
                vector_type rhs_f   = vector_type::Zero(fbs);

                for (auto& qp : fqd)
                {
                    auto f_phi = face_basis.eval_functions(msh, bfc, qp.point());

                    MFF += qp.weight() * f_phi * f_phi.transpose();
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
            auto binfo = msh.boundary_information(bfc);

            if(binfo.bndtype == boundary::DIRICHLET)
            {
                auto face_offset = face_id * fbs;
                auto face_offset_lagrange = (msh.faces_size() + face_i) * fbs;
                auto fqd   = face_quadrature.integrate(msh, bfc);

                vector_type rhs_f = vector_type::Zero(fbs);

                for (auto& qp : fqd)
                {
                    auto f_phi = face_basis.eval_functions(msh, bfc, qp.point());
                    rhs_f += qp.weight() * f_phi * bc(qp.point(), binfo.boundary_id);
                }

                for (size_t i = 0; i < fbs; i++)
                    rhs(face_offset_lagrange + i) = rhs_f(i);

                face_i++;
            }
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

    size_t                                      m_degree, m_degree_cell;

public:

    projector_pst()
        : m_degree(1), m_degree_cell(1)
    {
        cell_basis          = cell_basis_type(m_degree_cell);
        cell_quadrature     = cell_quadrature_type(2*m_degree_cell+ 2);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree+ 2);
    }
    #if 0
    projector_pst(size_t degree)
        : m_degree(degree), m_degree_cell(degree)
    {
        cell_basis          = cell_basis_type(m_degree);
        cell_quadrature     = cell_quadrature_type(2*m_degree + 2);
        face_basis          = face_basis_type(m_degree);
        face_quadrature     = face_quadrature_type(2*m_degree + 2);
    }
    #endif
    projector_pst(size_t degree, const bool& change_degree)
    : m_degree(degree), m_degree_cell(degree)
    {
        if(change_degree)
            m_degree_cell = degree + 1;

        cell_basis          = cell_basis_type(m_degree_cell);
        cell_quadrature     = cell_quadrature_type(2*m_degree_cell + 2);
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

                mm  += qp.weight() * phi * phi.transpose();
                rhs += qp.weight() * phi * fval;
            }

            ret.block(face_offset, 0, face_basis.size(), 1) = mm.llt().solve(rhs);
            face_offset += face_basis.size();
        }

        return ret;
    }

    template<typename Function>
    vector_type
    compute_face(const mesh_type& msh, const cell_type& cl,
                    const face_type & fc, const Function& f)
    {
        //if(mp.diff)
        auto number = set_cell_number(msh, cl);

        matrix_type mm = matrix_type::Zero(face_basis.size(), face_basis.size());
        vector_type rhs = vector_type::Zero(face_basis.size());

        auto face_quadpoints = face_quadrature.integrate(msh, fc);
        for (auto& qp : face_quadpoints)
        {
            auto phi = face_basis.eval_functions(msh, fc, qp.point());
            auto fval = f(qp.point(), number);

            mm  += qp.weight() * phi * phi.transpose();
            rhs += qp.weight() * phi * fval;
        }
        return mm.llt().solve(rhs);
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

    typedef dynamic_matrix<T>   matrix_type;
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

            ret = (0.5*pst.f * R * R/pst.mu)*ret;
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
                gradient_vector ret =  gradient_vector::Zero();

                ret(0) = - 0.5 * pst.f * p.x();
                ret(1) = - 0.5 * pst.f * p.y();
            return ret;
            //throw std::invalid_argument("st function is not intended for circular_tuyau problem. Review solution definition and parameters.txt");
            //return gradient_vector::Zero();
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
