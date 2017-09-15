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

template<typename MeshType, typename ElementType>
class tensor_bones
{
    typedef ElementType                         element_type;
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;

    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    matrix_type quad_evals_mat;
    vector_type bar_eval_vec;
    size_t m_quad_degree;
    size_t num_eval_pts;
};

template<typename T, typename Storage>
class tensor_bones<disk::mesh<T,2,Storage>, typename disk::mesh<T,2,Storage>::cell>
{
public:
    typedef disk::mesh<T,2,Storage>                 mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef disk::quadrature<mesh_type, cell_type>  cell_quad_type;
    typedef dynamic_matrix<T>         matrix_type;
    typedef dynamic_vector<T>         vector_type;
    typedef T scalar_type;

    matrix_type quad_evals_mat;
    vector_type bar_eval_vec;
    size_t m_quad_degree;
    size_t num_eval_pts;

    //WK: try to do this in the constructor.(I had problems in tensor_zero_vector,
    // since it s asking for the input of the constructor)

    tensor_bones(){};
    tensor_bones(const mesh_type& msh,
                 const cell_type& cell,
                 const size_t& quad_degree)
    {
        auto cq  = cell_quad_type(quad_degree);
        num_eval_pts = cq.integrate(msh, cell).size();
    };
    size_t      quad_degree()        { return m_quad_degree;}
    size_t      quad_evals_size()    { return num_eval_pts;}
    size_t      total_evals_size()   { return num_eval_pts + 1;}
    matrix_type quad_pts_matrix()    { return quad_evals_mat;}
    vector_type eval_at_barycenter() { return bar_eval_vec;}
    vector_type eval_at_quadpoint(size_t i) { return quad_evals_mat.col(i);}
    friend tensor_bones<mesh_type, cell_type> operator-
                            (const tensor_bones<mesh_type, cell_type>& tp1,
                            const tensor_bones< mesh_type, cell_type>& tp2)
    {
        tensor_bones<mesh_type, cell_type> ret;
        if(tp1.m_quad_degree != tp2.m_quad_degree)
            throw std::logic_error("(quad_degree) is different for these tensors, so it is not possible to do (-)operation.");
        ret.m_quad_degree  = tp1.m_quad_degree;
        ret.num_eval_pts   = ret.num_eval_pts;
        ret.quad_evals_mat = tp1.quad_evals_mat - tp2.quad_evals_mat;
        ret.bar_eval_vec   = tp1.bar_eval_vec   - tp2.bar_eval_vec;

        return ret;
    }
    friend tensor_bones<mesh_type, cell_type> operator+
                            (const tensor_bones<mesh_type, cell_type>& tp1,
                            const tensor_bones<mesh_type, cell_type>& tp2)
    {
        tensor_bones<mesh_type, cell_type> ret;
        if(tp1.m_quad_degree != tp2.m_quad_degree)
            throw std::logic_error("(quad_degree) is different for these tensors, so it is not possible to do (-)operation.");
        ret.m_quad_degree  = tp1.m_quad_degree;
        ret.num_eval_pts   = ret.num_eval_pts;
        ret.quad_evals_mat = tp1.quad_evals_mat + tp2.quad_evals_mat;
        ret.bar_eval_vec   = tp1.bar_eval_vec   + tp2.bar_eval_vec;

        return ret;
    }
    friend tensor_bones<mesh_type, cell_type> operator*(const scalar_type& scalar,
                            const tensor_bones<mesh_type, cell_type>& tp)
    {
        tensor_bones<mesh_type, cell_type> ret;
        ret.m_quad_degree  = tp.m_quad_degree;
        ret.num_eval_pts   = ret.num_eval_pts;
        ret.quad_evals_mat = scalar * tp.quad_evals_mat;
        ret.bar_eval_vec   = scalar * tp.bar_eval_vec;

        return ret;
    }

    friend vector_type operator*(const vector_type& vec,
                            const tensor_bones<mesh_type, cell_type>& tp)
    {
        vector_type ret;
        ret.m_quad_degree  = tp.m_quad_degree;
        ret.num_eval_pts   = ret.num_eval_pts;

        size_t i;
        for(i = 0; i < tp.quad_evals_mat.cols(); i++)
            ret(i)= make_prod(tp.quad_evals_mat, vec);

        ret.bar_eval_vec(i+1)   = make_prod( vec,tp.bar_eval_vec);

        return ret;
    }

    friend std::ostream& operator<<(std::ostream& os,
                                const tensor_bones<mesh_type, cell_type>& tp)
    {
        os << " * m_quad_degree : "<< tp.m_quad_degree<<std::endl;
        os << " * quad_tensor_matrix : size = "<< tp.quad_evals_mat.rows()<< "x";
        os <<  tp.cols()<< std::endl;
        for(size_t i = 0 ; i < tp.quad_evals_mat.rows(); i++)
        {
            for(size_t j = 0 ; j < tp.cols(); j++)
                os<< "  " << tp.values_matrix(i,j);
            os <<std::endl;
        }
        os << " * bar_tensor_vector : size = "<< tp.quad_evals_mat.rows()<< "x";
        os <<  1 << std::endl;
        for(size_t j = 0 ; j < tp.bar_eval_vec.rows(); j++)
            os << "  " << tp.bar_eval_vec(j,0) << std::endl;

        return os;
    }

    void
    Zero()
    {
        auto DIM = mesh_type::dimension;
        quad_evals_mat = matrix_type::Zero(DIM, num_eval_pts);
        bar_eval_vec   = vector_type::Zero(DIM);
        return;
    }

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

    matrix_type quad_evals_mat;
    matrix_type nodes_evals_mat;
    vector_type bar_eval_vec;

    size_t m_quad_degree;
    size_t num_eval_pts;
    size_t num_eval_quad_pts;

    //WK: try to do this in the constructor.(I had problems in tensor_zero_vector,
    // since it s asking for the input of the constructor)
    tensor_bones(){};
    tensor_bones(const mesh_type& msh,
                 const face_type& face,
                 const size_t& quad_degree)//: m_quad_degree(quad_degree)
    {
        auto fq  = face_quad_type(quad_degree);
        auto num_eval_quad_pts = fq.integrate(msh, face).size();
        num_eval_pts = num_eval_quad_pts + 2 + 1; //2 nodes ; 1 bar
    };
    size_t      quad_degree()        { return m_quad_degree;}
    size_t      quad_evals_size()    { return num_eval_quad_pts;}
    size_t      total_evals_size()   { return num_eval_pts;}
    matrix_type eval_quad_matrix()   { return quad_evals_mat;}
    matrix_type evals_at_nodes()     { return nodes_evals_mat;}
    vector_type eval_at_barycenter() { return bar_eval_vec;}
    vector_type eval_at_quadpoint(size_t i) { return quad_evals_mat.col(i);}
    vector_type eval_at_node(size_t i) { assert( i < 2 ); return nodes_evals_mat.col(i);}
    friend tensor_bones<mesh_type, face_type> operator-
                            (const tensor_bones<mesh_type, face_type>& tp1,
                            const tensor_bones< mesh_type, face_type>& tp2)
    {
        tensor_bones<mesh_type, face_type> ret;
        if(tp1.m_quad_degree != tp2.m_quad_degree)
            throw std::logic_error("(quad_degree) is different for these tensors, so it is not possible to do (-)operation.");
        ret.m_quad_degree  = tp1.m_quad_degree;
        ret.num_eval_pts   = ret.num_eval_pts;
        ret.quad_evals_mat = tp1.quad_evals_mat - tp2.quad_evals_mat;
        ret.bar_eval_vec   = tp1.bar_eval_vec   - tp2.bar_eval_vec;
        ret.nodes_evals_mat = tp1.nodes_evals_mat - tp2.nodes_evals_mat;

        return ret;
    }
    friend tensor_bones<mesh_type, face_type> operator+
                            (const tensor_bones<mesh_type, face_type>& tp1,
                            const tensor_bones<mesh_type, face_type>& tp2)
    {
        tensor_bones<mesh_type, face_type> ret;
        if(tp1.m_quad_degree != tp2.m_quad_degree)
            throw std::logic_error("(quad_degree) is different for these tensors, so it is not possible to do (-)operation.");
        ret.m_quad_degree  = tp1.m_quad_degree;
        ret.num_eval_pts   = ret.num_eval_pts;
        ret.quad_evals_mat = tp1.quad_evals_mat + tp2.quad_evals_mat;
        ret.bar_eval_vec   = tp1.bar_eval_vec   + tp2.bar_eval_vec;
        ret.nodes_evals_mat = tp1.nodes_evals_mat + tp2.nodes_evals_mat;

        return ret;
    }
    friend tensor_bones<mesh_type, face_type> operator*(const scalar_type& scalar,
                            const tensor_bones<mesh_type, face_type>& tp)
    {
        tensor_bones<mesh_type, face_type> ret;
        ret.m_quad_degree  = tp.m_quad_degree;
        ret.num_eval_pts   = ret.num_eval_pts;
        ret.quad_evals_mat = scalar * tp.quad_evals_mat;
        ret.bar_eval_vec   = scalar * tp.bar_eval_vec;
        ret.nodes_evals_mat = scalar * tp.nodes_evals_mat;

        return ret;
    }

    friend vector_type operator*(const vector_type& vec,
                            const tensor_bones<mesh_type, face_type>& tp)
    {
        vector_type ret = vector_type::Zero( num_eval_pts);
        ret.m_quad_degree  = tp.m_quad_degree;
        ret.num_eval_pts   = ret.num_eval_pts;

        size_t i;
        for(i = 0; i < tp.quad_evals_mat.cols(); i++)
            ret(i)= make_prod(tp.quad_evals_mat, vec);

        for(j= 0; j < tp.nodes_evals_mat.cols(); j++)
            ret(j + i)= make_prod(tp.nodes_evals_mat, vec);

        ret.bar_eval_vec(j + i + 1)   = make_prod( vec,tp.bar_eval_vec);

        return ret;
    }

    friend std::ostream& operator<<(std::ostream& os,
                                const tensor_bones<mesh_type, face_type>& tp)
    {
        os << " * m_quad_degree : "<< tp.m_quad_degree<<std::endl;
        os << " * quad_tensor_matrix : size = "<< tp.quad_evals_mat.rows()<< "x";
        os <<  tp.cols()<< std::endl;
        for(size_t i = 0 ; i < tp.quad_evals_mat.rows(); i++)
        {
            for(size_t j = 0 ; j < tp.cols(); j++)
                os<< "  " << tp.values_matrix(i,j);
            os <<std::endl;
        }
        os << " * bar_tensor_vector : size = "<< tp.quad_evals_mat.rows()<< "x";
        os <<  1 << std::endl;
        for(size_t j = 0 ; j < tp.bar_eval_vec.rows(); j++)
            os << "  " << tp.bar_eval_vec(j,0) << std::endl;

        return os;
    }

    void
    Zero()
    {
        auto DIM = mesh_type::dimension;
        quad_evals_mat = matrix_type::Zero(DIM, num_eval_pts);
        bar_eval_vec   = vector_type::Zero(DIM);
        return;
    }
    matrix_type
    get_all()
    {
        matrix_type ret = matrix_type::Zero(2, num_eval_pts);

        auto ret.block(0, 0, 2, num_eval_quad_pts) = quad_evals_mat;
        auto ret.block(0, num_eval_quad_pts, 2, 1) = bar_eval_vec;

        return ret;
    }
};

template<typename T, typename Storage>
std::vector<point<T,2>>
tensor_points(const disk::mesh<T, 2, Storage>& msh,
              const typename disk::mesh<T, 2, Storage>::cell & cell,
              const size_t m_quad_degree )
{
    typedef disk::mesh<T,2,Storage>                 mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef typename mesh_type::face                face_type;
    typedef disk::quadrature<mesh_type, cell_type>  cell_quad_type;
    typedef disk::quadrature<mesh_type, face_type>  face_quad_type;

    auto cq  = cell_quad_type(m_quad_degree);
    auto fq  = face_quad_type(m_quad_degree);
    auto cqs = cq.integrate(msh, cell);

    auto fcs = faces(msh, cell);
    auto fqs = fq.integrate(msh, fcs[0]);

    auto pts = points(msh, cell);

    auto total_size = cqs.size()+ pts.size() * fqs.size() + 2 * pts.size() + 1;
    auto ret = std::vector<typename mesh_type::point_type>(total_size);

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
};


template<typename MeshType>
class tensors2
{
public:

    typedef MeshType   mesh_type;
    typedef typename mesh_type::cell_type      cell_type;
    typedef typename mesh_type::face_type      face_type;
    typedef typename mesh_type::scalar_type    scalar_type;
    typedef std::vector<tensor_bones<MeshType, cell_type>>  cell_tensor_vector;
    typedef std::vector<tensor_bones<MeshType, face_type>>  face_tensor_vector;
    typedef dynamic_matrix<scalar_type>        matrix_type;
    typedef dynamic_vector<scalar_type>        vector_type;

    cell_tensor_vector  cells;
    face_tensor_vector  faces;


    tensors2() //quad_degree(2*degree+2)
    {}

    tensors2(const MeshType& msh) //quad_degree(2*degree+2)
    {
        cells = cell_tensor_vector(msh.cells_size());
        faces = face_tensor_vector(msh.faces_size());
    }

    tensors2(const size_t cells_size, const size_t faces_size) //quad_degree(2*degree+2)
    {
        cells = cell_tensor_vector(cells_size);
        faces = face_tensor_vector(faces_size);
    }
    void
    tensors2(const mesh_type msh, const size_t& degree)
    {
        cells = cell_tensor_vector(msh.cells_size());
        faces = face_tensor_vector(msh.faces_size());

        for (auto& cl : msh)
        {
            auto  id = cl.get_id();
            tensor_bones<mesh_type, cell_type>  ptsr(msh, cl, degree);
            ptsr.Zero();
            cells.at(id)  =  ptsr;
        }

        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            auto face = *itor;
            auto id = msh.lookup(face);
            tensor_bones<mesh_type, face_type>  ptsr(msh, face, degree);;
            ptsr.Zero();
            faces.at(id) =  ptsr;
        }
        return;
    };
    void
    zero_vector(const mesh_type msh, const size_t& degree)
    {
        cells = cell_tensor_vector(msh.cells_size());
        faces = face_tensor_vector(msh.faces_size());

        for (auto& cl : msh)
        {
            auto  id = cl.get_id();
            tensor_bones<mesh_type, cell_type>  ptsr(msh, cl, degree);
            ptsr.Zero();
            cells.at(id)  =  ptsr;
        }

        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            auto face = *itor;
            auto id = msh.lookup(face);
            tensor_bones<mesh_type, face_type>  ptsr(msh, face, degree);
            ptsr.Zero();
            faces.at(id) =  ptsr;
        }
        return;
    };
    auto  tensor_at_cells() {  return cells;}
    auto  tensor_at_faces() {  return faces;}

    typedef typename cell_tensor_vector::iterator         cell_iterator;
    typedef typename cell_tensor_vector::const_iterator   const_cell_iterator;

    cell_iterator           at_cells_begin() { return cells.begin(); }
    cell_iterator           at_cells_end()   { return cells.end(); }
    const_cell_iterator     at_cells_begin() const { return cells.begin(); }
    const_cell_iterator     at_cells_end()   const { return cells.end(); }

    /* face iterators */
    typedef typename face_tensor_vector::iterator        face_iterator;
    typedef typename face_tensor_vector::const_iterator  const_face_iterator;

    face_iterator           at_faces_begin() { return faces.begin(); }
    face_iterator           at_faces_end()   { return faces.end(); }
    const_face_iterator     at_faces_begin() const { return faces.begin(); }
    const_face_iterator     at_faces_end()   const { return faces.end(); }
    matrix_type

    #if 0

    matrix_type
    at_cell(const mesh_type & msh, const cell_type& cl)
    {
        auto id = msh.lookup(cl);
        return cells.at(id);
    };

    matrix_type
    at_cell(const mesh_type & msh, const cell_iterator & itor)
    {
        auto id = msh.lookup(*itor);
        return cells.at(id);
    };

    matrix_type
    at_cell(const mesh_type & msh, const cell_type& cl, const size_t& col)
    {
        auto id = msh.lookup(cl);
        return cells.at(id).col(col);
    };


    matrix_type
    at_face(const mesh_type & msh, const face_type& fc,  const size_t& col)
    {
        auto id = msh.lookup(fc);
        return faces.at(id).col(col);
    };

    matrix_type
    at_face(const mesh_type & msh, const face_iterator & itor)
    {
        auto id = msh.lookup(*itor);
        return faces.at(id);
    };

    void
    put_at_cell(const mesh_type & msh, const cell_type& cl, const matrix_type& mat)
    {
        auto id = msh.lookup(cl);
        assert(cells.at(id).cols() == mat.cols());
        assert(cells.at(id).rows() == mat.rows());

        cells.at(id) = mat;
        return;
    };
    void
    put_at_face(const mesh_type & msh, const face_type& fc, const matrix_type& mat)
    {
        auto id = msh.lookup(fc);
        assert(faces.at(id).cols() == mat.cols());
        assert(faces.at(id).rows() == mat.rows());

        faces.at(id) = mat;
        return;
    };
    template<typename T>
    void()
    put_tensor_at_cell(const size_t i,
               const dynamic_matrix<T>& mat,
               const dynamic_vector<T>& vec
               const size_t quad_degree)
    {
        cells.at(i).quad_evals_mat = mat;
        cells.at(i).bar_eval_vec   = vec;
        cells.at(i).num_eval_pts   = mat.cols() + 1;
        cells.at(i).m_quad_degree  = m_quad_degree;
    }

    template<typename T>
    void()
    put_tensor_at_face(const size_t i,
               const dynamic_matrix<T>& mat,
               const dynamic_vector<T>& vec
               const size_t quad_degree)
    {
        faces.at(i).quad_evals_mat = mat;
        faces.at(i).nodes_evals_mat = ;
        faces.at(i).bar_eval_vec   = vec;
        faces.at(i).num_eval_pts   = mat.cols() + 1;
        faces.at(i).m_quad_degree  = m_quad_degree;
    }
    #endif
    void()
    put_at_quad_pts(const mesh_type& msh,
                    const cell_type & cl,
                    const matrix_type& mat)
    {
        auto id = msh.lookup(cl);
        cells.at(id).quad_evals_mat = mat;
    }
    void()
    put_at_quad_pts(const mesh_type& msh,
                    const face_type & fc,
                    const matrix_type& mat)
    {
        auto id = msh.lookup(fc);
        faces.at(id).quad_evals_mat = mat;
    }
    void()
    put_at_barycenter(const mesh_type& msh,
                    const cell_type & cl,
                    const vector_type& vec)
    {
        auto id = msh.lookup(cl);
        cells.at(id).bar_eval_vec = vec;
    }
    void()
    put_at_barycenter(const mesh_type& msh,
                      const face_type & fc,
                      const vector_type& vec)
    {
        auto id = msh.lookup(fc);
        faces.at(id).bar_eval_vec = vec;
    }
    void()
    put_at_nodes(const mesh_type& msh,
                 const face_type & fc,
                 const matrix_type& mat)
    {
        auto id = msh.lookup(fc);
        faces.at(id).nodes_evals_mat = mat;
    }
    at_quad_pts(const mesh_type & msh, const cell_type& cl)
    {
        auto id = msh.lookup(cl);
        return cells.at(id).quad_evals_mat;
    }
    at_quad_pts(const mesh_type & msh, const face_type& fc)
    {
        auto id = msh.lookup(fc);
        return faces.at(id).quad_evals_mat;
    }
    at_quad_point(const mesh_type & msh, const cell_type& cl, const size_t col)
    {
        auto id = msh.lookup(cl);
        assert(col < cells.at(id).quad_evals_mat.cols());
        return cells.at(id).quad_evals_mat.col(col);
    }
    at_quad_point(const mesh_type & msh, const face_type& fc, const size_t col)
    {
        auto id = msh.lookup(fc);
        assert(col < faces.at(id).quad_evals_mat.cols());
        return faces.at(id).quad_evals_mat.col(col);
    }
    matrix_type
    at_face_nodes(const mesh_type & msh, const face_type& face)
    {
        auto id = msh.lookup(fc);
        return faces.at(id).nodes_evals_mat;
    }
    vector_type
    at_barycenter(const mesh_type & msh, const face_type& face)
    {
        auto id = msh.lookup(fc);
        return faces.at(id).bar_eval_vec;
    }
    auto
    quad_degree(const mesh_type & msh, const cell_type& cell)
    {
        auto id = msh.lookup(cell);
        return cells.at(id).m_quad_degree;
    }
    auto
    quad_degree(const mesh_type & msh, const face_type& face)
    {
        auto id = msh.lookup(fc);
        return faces.at(id).m_quad_degree;
    }
};

template<typename T,typename Mesh,typename CellBasisType, typename CellQuadType,
                       typename FaceBasisType, typename FaceQuadType>
class plasticity2
{
    typedef Mesh                            mesh_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef dynamic_matrix<scalar_type>     matrix_type;
    typedef dynamic_vector<scalar_type>     vector_type;

    typedef CellBasisType                   cell_basis_type;
    typedef FaceBasisType                   face_basis_type;
    cell_basis_type                         cell_basis;
    face_basis_type                         face_basis;

    typedef CellQuadType                    cell_quadrature_type;
    typedef FaceQuadType                    face_quadrature_type;
    cell_quadrature_type                    cell_quadrature;
    face_quadrature_type                    face_quadrature;

    plasticity_data<T>          pst;
    size_t                      m_degree;
    T                           m_method_coef;
    std::vector<matrix_type>    m_rec_oper_Th;
    std::vector<vector_type>    m_velocity_Th;
    tensors2<mesh_type>         m_stress_Th;

    plasticity2(const size_t& degree, const size_t& quad_degree,
                const plasticity_data<T> pst_data): pst(pst_data),
                                    m_degree(degree)
    //WK: for mu as tensor one should check this part
    {
        cell_basis  =  cell_basis_type(m_degree + 1); // This is k + 1  just if I want to use \nabla r^{k + 1}
        face_basis  =  face_basis_type(m_degree);

        cell_quadrature  = cell_quadrature_type(quad_degree);
        face_quadrature  = face_quadrature_type(quad_degree);
        if(pst.method == true)
            m_method_coef = 1. / pst.alpha ;
        else
            m_method_coef = 1. / (pst.alpha + pst.mu);
    }

    vector_type
    grad_ruh_at( const point<T,2>& pt, //Check for this to be generalized to other dimensions.
                   const mesh_type&  msh,
                   const cell_type&  cl)
    {
        auto id = cl.get_id();
        matrix_type  rec_op = m_rec_oper_Th.at(id);
        vector_type  uh_TF  = m_velocity_Th.at(id);

        auto col_range      = cell_basis.range(1,m_degree+1);
        auto row_range      = disk::dof_range( 0,mesh_type::dimension);

        auto dphi   = cell_basis.eval_gradients(msh, cl, pt);
        matrix_type dphi_matrix = make_gradient_matrix(dphi);
        matrix_type dphi_taken  = take(dphi_matrix, row_range, col_range);
        vector_type dphi_ruh    = dphi_taken * rec_op * uh_TF;

        return dphi_ruh;
    }
    std::pair<vector_type, vector_type>
    grad_uh_at( const point<T,2>& pt, //Check for this to be generalized to other dimensions.
                   const mesh_type&  msh,
                   const size_t&  id)
    {
        matrix_type  rec_op = m_rec_oper_Th.at(id);
        vector_type  uh_TF  = m_velocity_Th.at(id);

        auto col_range      = cell_basis.range(1,m_degree+1);
        auto row_range      = disk::dof_range( 0,mesh_type::dimension);

        auto cell  =  *std::next(msh.cells_begin(), id);
        auto dphi  =  cell_basis.eval_gradients(msh, cell, pt);
        matrix_type dphi_matrix = make_gradient_matrix(dphi);
        matrix_type dphi_taken  = take(dphi_matrix, row_range, col_range);
        vector_type dphi_ruh    = dphi_taken * rec_op * uh_TF;

        auto cell_range     = cell_basis.range(0,m_degree);
        vector_type  uh_T   = uh_TF.head(cell_range.size());
        matrix_type  dphi_T = take(dphi_matrix, row_range, cell_range.size());
        vector_type  dphi_uh_T = dphi_T * uh_T;

        return std::make_pair(dphi_ruh, dphi_uh_T);
    }

    vector_type
    eval_plasticity(const point<T,2> & pt,
                  const mesh_type    & msh,
                  const cell_type    & cl,
                  const vector_type  & sigma)
    {
        vector_type  theta =  vector_type::Zero(mesh_type::dimension);
        auto grad_uh = grad_uh_at(pt, msh, cl);

        vector_type dphi_ruh  = grad_uh.first;
        vector_type dphi_uh_T = grad_uh.second;
        if(pst.hho == 1)
             theta = sigma  +  pst.alpha * dphi_ruh;
        if(pst.hho == 2)
            theta = sigma  +  pst.alpha * dphi_uh_T;

        auto theta_norm = theta.norm();

        vector_type gamma = vector_type::Zero(mesh_type::dimension);
        if(theta_norm > pst.yield)
            gamma = m_method_coef * (theta_norm - pst.yield) * (theta / theta_norm);

        return sigma + pst.alpha * ( dphi_ruh - gamma );
    }

    vector_type
    eval_plasticity(const point<T,2> & pt,
                  const mesh_type    & msh,
                  const face_type    & face,
                  const std::pair<int,int> & owners_ids,
                  const  vector_type  & sigma)
    {
        size_t cell_id_1 = owners_ids.first;
        size_t cell_id_2 = owners_ids.first;

        if(!msh.is_boundary(face))
            cell_id_2 = owners_ids.second;

        auto grad_uh_1 = grad_uh_at(pt, msh, size_t(cell_id_1));
        auto grad_uh_2 = grad_uh_at(pt, msh, size_t(cell_id_2));

        vector_type  dphi_ruh_T1 = grad_uh_1.first;
        vector_type  dphi_ruh_T2 = grad_uh_2.first;

        vector_type  dphi_ruh    = T(0.5) * (dphi_ruh_T1 + dphi_ruh_T2);

        vector_type  theta =  vector_type::Zero(mesh_type::dimension);

        if(pst.hho == 1)
             theta = sigma  +  pst.alpha * dphi_ruh;
        if(pst.hho == 2)
        {
            vector_type dphi_uh_T1 = grad_uh_1.second;
            vector_type dphi_uh_T2 = grad_uh_2.second;

            theta = sigma  +  pst.alpha * T(0.5) *(dphi_uh_T1 + dphi_uh_T2);
        }
        auto theta_norm = theta.norm();

        vector_type gamma = vector_type::Zero(mesh_type::dimension);
        if(theta_norm > pst.yield)
            gamma = m_method_coef * (theta_norm - pst.yield) * (theta / theta_norm);

        return sigma + pst.alpha * ( dphi_ruh - gamma );
    }
    void
    compute(const mesh_type  & msh,
            const std::vector<vector_type>& velocity_Th,
            const std::vector<matrix_type>& rec_opers_Th,
            tensors2<mesh_type>& stress_Th)
    {
        m_velocity_Th = velocity_Th;
        m_stress_Th   = stress_Th;
        m_rec_oper_Th = rec_opers_Th;

        compute_at_cells(msh);
        compute_at_faces(msh);

        stress_Th = m_stress_Th;
    }
    void
    compute_at_cells(const mesh_type  & msh)
    {
        for(auto& cl : msh)
        {
            //WK: for p-adaptation quad_degree is not the same for all cells
            size_t cont  =  0;
            auto   cqs   =  cell_quadrature.integrate(msh, cl);
            vector_type new_sigma = vector_type::Zero(cqs.size());

            for (auto& qp : cqs)
            {
                if(cont > cqs.size())
                    throw std::logic_error("Wrong access to columns");

                vector_type sigma = m_stress_Th.at_quad_point(msh, face, col);
                new_sigma.col(cont) = eval_plasticity(qp.point(), msh, cl, sigma);
                ++cont;
            }
            vector_type old_sigma = m_stress_Th.at_barycenter(msh, face);
            vector_type bar_sigma = eval_plasticity(msh, cl, old_sigma);

            m_stress_Th.put_at_quad_pts(msh, cl, new_sigma);
            m_stress_Th.put_at_barycenter(msh, face, bar_sigma);
        }
        return;
    }
    void
    compute_at_faces(const mesh_type    & msh)
    {
        for(auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            //WK: for p-adaptation quad_degree is not the same for all cells
            auto face   = *itor;
            auto pts    = points(msh, face);
            auto owners = owners_cells(face);

            size_t col =  0;
            vector_type quad_sigma = vector_type::Zero(pts.size());
            for(auto& p: pts)
            {
                assert(col < pts.size())
                vector_type  sigma = m_stress_Th.at_quad_point(msh, face, col);
                nodes_sigma.col(p_cont) = eval_plasticity(p, msh, face, owners, sigma);
                col++;
            }

            col  =  0;
            auto fqs = face_quadrature.integrate(msh, face);
            vector_type qps_sigma = vector_type::Zero(fqs.size());

            for (auto& qp : fqs)
            {
                assert(col < fqs.size())
                vector_type  sigma = m_stress_Th.at_quad_point(msh, face, col);
                qps_sigma.col(col) = eval_plasticity(qp.point(), msh, face, owners, col);
                ++col;
            }
            vector_type  sigma = m_stress_Th.at_barycenter(msh, face);
            vector_type bar_sigma = eval_plasticity(msh, face, owners, sigma);

            m_stress_Th.put_at_nodes(msh, face, nodes_sigma);
            m_stress_Th.put_at_quad_pts(msh, face, qps_sigma);
            m_stress_Th.put_at_barycenter(msh, face, bar_sigma);
        }
        return;
    }

};// plasticity


#if 0
{

    auto   cqs   =  cell_quadrature.integrate(msh, cl);
    size_t cont  =  0;

    for (auto& qp : cqs)
    {
        if(cont > cqs.size())
            throw std::logic_error("Wrong access to columns");
        vector_type sigma =
        vector_type gamma =

        //( Dv_T, sigma - alpha* gamma)_T
        auto dphi   = cell_basis.eval_gradients(msh, cl, qp.point());
        matrix_type   dphi_matrix = make_gradient_matrix(dphi);

        if(pst.hho == 2)
        {
            matrix_type dphi_T  = take(dphi_matrix, row_range, cell_range);
            rhs.head(num_cell_dofs)  += qp.weight() * dphi_T.transpose() *
                                                (sigma - pst.alpha * gamma);
        }
        else if(pst.hho ==  1)
        {
            matrix_type dphi_TF = take(dphi_matrix, row_range, col_range);
            matrix_type dphi_rec = dphi_TF * rec_oper;
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

        auto sigma_agamma =  f_sigma - pst.alpha * f_gamma;
        auto p1  =  n * sigma_agamma ;

        vector_type f_sigma =  ;
        vector_type f_gamma  = ;

        for (auto& qp : fqs)
        {
            auto f_phi  = face_basis.eval_functions(msh, fc, qp.point());
            auto c_phi  = cell_basis.eval_functions(msh, cl, qp.point());

            // ((sigma - alpha*gamma) * n
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
        }
    }

    return rhs;

}
#endif
}//Disk
