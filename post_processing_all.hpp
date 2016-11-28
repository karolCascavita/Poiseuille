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

 #include "hho/hho.hpp"
 #include "hho/hho_pst.hpp"
namespace disk
{

template<typename T, size_t DIM, typename Storage>
class post_processing_base
{
    static_assert(DIM > 0 && DIM <= 3, "mesh: Allowed dimensions are 1, 2 and 3");
};

template<typename T, typename Storage>
class post_processing_base<T,1, Storage>
{
    public:

    typedef mesh<T,1,Storage>                   mesh_type;
    typedef typename mesh_type::point_type      point_type;
    typedef std::vector<point_type>             point_vector_type;
    size_t m_degree;
    size_t m_num_sub_nodes,   m_num_sub_elems;
    size_t m_total_sub_nodes, m_total_sub_elems;

    point_vector_type   m_sub_nodes;
    point_vector_type   m_vertices;

//    post_processing_base(const mesh<T,1,Storage>& msh)
    post_processing_base()
    {}

    void
    plot_nodes(size_t n, const std::vector<point_type> & vers, std::vector<point_type> & nodes)
    {
        auto p1 = vers[0];
        auto p2 = vers[1];

        if(n == 0)
        {
            nodes.push_back(p1);
            nodes.push_back(p2);
        }
        else
        {
            T h = std::abs(p1.x()-p2.x())/n;
            for(size_t i = 0; i < n + 1;  i++)
            {
                point_type p = point<T,1>({h*i});
                nodes.push_back(p);
            }
        }
    }
    std::vector<point_type>
    make_vtk_points(const mesh_type& msh,const size_t degree)
    {
        m_degree = degree;
        if(m_degree == 0)
            m_degree = 1; /* Same treatment as linear functions*/

        m_num_sub_nodes   = m_degree+1;
        m_num_sub_elems   = 1;
        m_total_sub_nodes = 2*m_num_sub_nodes*msh.cells_size();
        m_total_sub_elems = 2*msh.cells_size();
        m_sub_nodes.reserve(m_num_sub_nodes);
        m_vertices.reserve(2); /*WK: do it without reserve*/
        m_vertices[0]     = point<T,1>({0.});
        m_vertices[1]     = point<T,1>({1.});

        plot_nodes(m_degree, m_vertices, m_sub_nodes);
        std::vector<point_type> test_points(m_total_sub_nodes);
        size_t cont = 0;

        for(auto& cl : msh)
        {
            auto cl_faces = faces(msh, cl);
            auto bar = barycenter(msh, cl);
            auto pts = points(msh, cl);
            auto h   = std::abs(pts[1].x() - pts[0].x())/2.;
            size_t ifc = 0;
            for(auto& fc:cl_faces)
            {
                for(size_t  i = 0; i < m_num_sub_nodes; i++)
                {
                    auto p     = m_sub_nodes[i];
                    auto c     = point<T,1> ({h*ifc});
                    int  idx   = cont*m_num_sub_nodes + i;

                    test_points[idx] = h * p + pts[0] + c;
                }
                ++cont;
                ++ifc;
            }
        }
        return test_points;
    }
};

template<typename T, typename Storage>
class post_processing_base<T,2,Storage>
{
public:
    typedef mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef std::vector<point_type>                 point_vector_type;

    size_t m_degree;
    size_t m_num_sub_nodes,   m_num_sub_elems;
    size_t m_total_sub_nodes, m_total_sub_elems;
    point_vector_type   m_sub_nodes;
    point_vector_type   m_vertices;


    //post_processing_base(const mesh<T,2,Storage>& msh)
    post_processing_base()
    {}

    void
    plot_nodes(const size_t n, const std::vector<point_type>& vers, std::vector<point_type> & nodes)
    {
        /*WK: chechk for degree 0*/
        auto p1 = vers[0];
        auto p2 = vers[1];
        auto p3 = vers[2];

        if(n == 0)
        {
            nodes.push_back(p1);
            nodes.push_back(p2);
            nodes.push_back(p3);
        }
        else
        {
            auto p4 = (p1 + p2)/2.;
            auto p5 = (p2 + p3)/2.;
            auto p6 = (p1 + p3)/2.;

            point_vector_type vertices_T1 = point_vector_type({p1, p4, p6});
            point_vector_type vertices_T2 = point_vector_type({p4, p5, p6});
            point_vector_type vertices_T3 = point_vector_type({p4, p2, p5});
            point_vector_type vertices_T4 = point_vector_type({p6, p5, p3});

            plot_nodes(n-1, vertices_T1, nodes);
            plot_nodes(n-1, vertices_T2, nodes);
            plot_nodes(n-1, vertices_T3, nodes);
            plot_nodes(n-1, vertices_T4, nodes);
        }
    }

    std::vector<point_type>
    make_vtk_points(const mesh_type& msh, const size_t degree)
    {
        m_degree = degree;
        if(m_degree == 0)
            m_degree = 1; /* Same treatment as linear functions*/

        m_num_sub_nodes   = 3*std::pow(4,m_degree);
        m_num_sub_elems   = std::pow(4,m_degree);
        m_total_sub_nodes = m_num_sub_nodes*(msh.boundary_faces_size() +  2*msh.internal_faces_size() );
        m_total_sub_elems = m_num_sub_elems*(msh.boundary_faces_size() +  2*msh.internal_faces_size() );
        m_sub_nodes.reserve(m_num_sub_nodes);
        m_vertices.reserve(3);
        m_vertices[0]   =   point<T,2>({0.,0.});
        m_vertices[1]   =   point<T,2>({1.,0.});
        m_vertices[2]   =   point<T,2>({0.,1.});

        plot_nodes(m_degree, m_vertices, m_sub_nodes);

        std::vector<point_type> test_points(m_total_sub_nodes);
        size_t cont = 0;

        for(auto& cl : msh)
        {
            auto cl_faces = faces(msh, cl);
            auto bar = barycenter(msh, cl);
            for(auto& fc : cl_faces)
            {
                auto pts = points(msh, fc);
                for(size_t  i = 0; i < m_num_sub_nodes; i++)
                {
                    auto p    = m_sub_nodes[i];
                    int idx   = cont*m_num_sub_nodes;
                    test_points[idx+i] = ((bar + (pts[1] - bar) * p.x() ) + (pts[0] - bar) * p.y());
                }
            ++cont;
            }
        }
        return test_points;
    }
};

template<typename T, size_t DIM, typename Storage>
class post_processing: public post_processing_base<T,DIM,Storage>
{
public:

    typedef typename post_processing_base<T,DIM,Storage>::mesh_type    mesh_type;
    typedef typename mesh_type::point_type              point_type;
    typedef typename mesh_type::scalar_type             scalar_type;
    typedef typename mesh_type::cell                    cell_type;
    typedef typename mesh_type::face                    face_type;
    typedef post_processing_base<T,DIM,Storage>         pp_base;
    typedef dynamic_vector<scalar_type>                 vector_type;
    typedef dynamic_matrix<scalar_type>                 matrix_type;
    typedef std::vector<point_type>                     point_vector_type;

    size_t m_degree;

    point_vector_type sub_nodes;
    point_vector_type vertices;


//    post_processing(const mesh<T,DIM,Storage>& msh): pp_base(msh)
    post_processing():
     pp_base()
    { }

    void vtk_writer(const mesh_type& msh, const std::string& name, const size_t degree,
               const std::vector<dynamic_vector<T>> & Uh_Th, const int iter)   /* storage of velocities*/
    {
        typedef quadrature<mesh_type, cell_type>     cell_quadrature_type;
        typedef quadrature<mesh_type, face_type>     face_quadrature_type;

        typedef scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
        typedef scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

        gradient_reconstruction_nopre<mesh_type,
                                            cell_basis_type,
                                            cell_quadrature_type,
                                            face_basis_type,
                                            face_quadrature_type> gradrec_nopre(degree);


        cell_basis_type                 cb(degree+1);
        std::stringstream               vtk_points, vtk_elems, vtk_scalar_data, vtk_vector_data;

        size_t num_sub_elems, num_sub_nodes;
        size_t total_sub_nodes, total_sub_elems;
        auto   test_points  = pp_base::make_vtk_points(msh,degree);

        total_sub_nodes = pp_base::m_total_sub_nodes;
        total_sub_elems = pp_base::m_total_sub_elems;
        num_sub_nodes   = pp_base::m_num_sub_nodes;
        sub_nodes       = pp_base::m_sub_nodes;
        vertices        = pp_base::m_vertices;

        size_t cl_cont = 0;
        size_t fc_cont = 0;

        char numstr[4];
        std::string result;
        sprintf(numstr, "%i",iter);
        std::string file_name = name + numstr;

        for(auto& cel : msh)
        {
            gradrec_nopre.compute(msh, cel);
            vector_type uh_TF =  Uh_Th[cl_cont];
            dynamic_vector<scalar_type> rec(cb.size());
            rec.tail(cb.size()-1) = gradrec_nopre.oper * uh_TF;
            rec(0) = uh_TF(0);

            auto cl_faces =faces(msh,cel);
            for(auto& fc : cl_faces)
            {
                for (size_t itp = 0; itp < num_sub_nodes; itp++)
                {
                    auto pot  = 0.;
                    auto idx  = itp + num_sub_nodes*fc_cont;
                    auto tp   = test_points[idx];
                    auto phi  = cb.eval_functions(msh, cel, tp);
                    auto dphi = cb.eval_gradients(msh, cel, tp);

                    for (size_t d = 0; d < DIM; d++)
                        vtk_points << tp[d] << " ";
                    for (size_t d = 0; d < 3 - DIM; d++) /* VTK only supports until 3D */
                        vtk_points<< 0. << " ";
                    vtk_points<< '\n';

                    for (size_t i = 0; i < cb.size(); i++)
                        pot  +=  phi[i] * rec(i);
                    vtk_scalar_data << pot <<' ';

                    auto col_range  =   cb.range(1,degree+1);
                    auto row_range  =   dof_range(0,DIM);

                    matrix_type dphi_matrix =  make_gradient_matrix(dphi);
                    matrix_type dphi_high_order =   take(dphi_matrix, row_range, col_range);
                    vector_type dpot =   dphi_high_order*rec.tail(col_range.size());
                    for (size_t d = 0; d < DIM; d++)
                        vtk_vector_data << dpot(d) << " ";
                    for (size_t d = 0; d < 3 - DIM; d++) /* VTK only supports until 3D */
                        vtk_vector_data<< 0. << " ";
                        vtk_vector_data<< '\n';
                }
                ++fc_cont;
            }
            ++cl_cont;
        }

        size_t el = 0;
        std::vector<size_t> vtk_cell_type = {4,5,10}; /*WK: check this for others VTK types*/
        std::vector<size_t> vtk_nodes_inside_cell = {degree+1,3,4};/* WK: all cells are the same type */
        for (size_t i = 0; i < total_sub_elems; i++)
        {
            vtk_elems << vtk_nodes_inside_cell[DIM-1] << " ";
            for(size_t j=0; j < vtk_nodes_inside_cell[DIM-1]; j++, el++)
                vtk_elems << el<< " ";
                vtk_elems <<  '\n';
        }

        /*VTK file*/
        std::string ext = ".vtk";
        std::ofstream  ofs(file_name + ext);
        ofs << "# vtk DataFile Version 3.0"<< '\n'
            << "#This file was generated by the DISK++ library"<< '\n'
            << "ASCII"<< '\n';

        /* Data header */
        ofs << "DATASET UNSTRUCTURED_GRID\n"
            << '\n';

        /* Points */
        ofs << "POINTS " << total_sub_nodes<<" double"<<'\n'
            << vtk_points.str()
            << '\n';

        /* Cells */
        ofs << "CELLS " << total_sub_elems <<' '<<  total_sub_elems *(vtk_nodes_inside_cell[DIM-1] + 1)<< '\n'
            << vtk_elems.str()
            << '\n';

        /* Types of cells*/
        ofs << "CELL_TYPES " << total_sub_elems << '\n';
            for(size_t i = 0; i < total_sub_elems; i++)
                ofs << ' ' << vtk_cell_type[DIM-1];
        ofs << '\n';

        /* Data */
        ofs << "POINT_DATA " << total_sub_nodes
            << '\n';

        ofs << "VECTORS "
            << " Gr_U"
            << " double"
            << '\n'
            << vtk_vector_data.str()
            << '\n';

        ofs << "SCALARS "
            << " U "
            << " double 1"
            << '\n'
            << "LOOKUP_TABLE default"
            << '\n'
            << vtk_scalar_data.str()
            << '\n';
        ofs.close();
    }

    template<typename TensorMatrix>
    void tensor_norm_vtk_writer(const mesh_type& msh, const std::string& name, const size_t degree,
               const std::vector<TensorMatrix> & tensor, const int iter)   /* storage of velocities*/
    {

        typedef quadrature<mesh_type, cell_type>     cell_quadrature_type;
        typedef quadrature<mesh_type, face_type>     face_quadrature_type;

        typedef scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
        typedef scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

        cell_quadrature_type    cq(2*degree+2);
        std::stringstream       vtk_points, vtk_elems, vtk_data_x, vtk_data_y, vtk_data_norm;

        size_t cl_cont = 0;
        size_t total_num_qp  = 0;
        char numstr[4];
        std::string result;
        sprintf(numstr, "%i",iter);
        std::string file_name = name + numstr;

        for(auto& cel : msh)
        {
            size_t qp_cont = 0;
            auto cell_id   = msh.lookup(cel);
            auto cell_quadpoints     =  cq.integrate(msh, cel);
            auto local_tensor_storage =  tensor[cell_id];

            total_num_qp += cell_quadpoints.size();

            for (auto& qp : cell_quadpoints)
            {
                auto tp  = qp.point();

                for (size_t d = 0; d < DIM; d++)
                        vtk_points << tp[d] << " ";
                for (size_t d = 0; d < 3 - DIM; d++) /* VTK only supports until 3D */
                        vtk_points<< 0. << " ";
                vtk_points<< '\n';

                vector_type pot = local_tensor_storage.col(qp_cont);
                if( DIM == 1)
                    vtk_data_norm << pot << " ";
                if (DIM == 2)
                    vtk_data_norm << pot.norm() << " ";
                ++qp_cont;
            }
            ++cl_cont;
        }

        size_t el = 0;

        for (size_t i = 0; i < total_num_qp; i++)
            vtk_elems << 1 <<" "<< i <<  '\n';

        /*VTK file*/
        std::string ext = ".vtk";
        std::ofstream  ofs(file_name + ext);
        ofs << "# vtk DataFile Version 3.0"<< '\n'
            << "#This file was generated by the DISK++ library"<< '\n'
            << "ASCII"<< '\n';

        /* Data header */
        ofs << "DATASET UNSTRUCTURED_GRID\n"
            << '\n';

        /* Points */
        ofs << "POINTS " << total_num_qp <<" double"<<'\n'
            << vtk_points.str()
            << '\n';

        /* Cells */
        ofs << "CELLS " << total_num_qp <<' '<<  total_num_qp *(1 + 1)<< '\n'
            << vtk_elems.str()
            << '\n';

        /* Types of cells*/
        ofs << "CELL_TYPES " << total_num_qp << '\n';
            for(size_t i = 0; i < total_num_qp; i++)
                ofs << " " << 1;
        ofs << '\n';

        /* Data */
        ofs << "POINT_DATA " << total_num_qp
            << '\n';

        ofs << "SCALARS "
            << " Lambda_norm "
            << " double 1"
            << '\n'
            << "LOOKUP_TABLE default"
            << '\n'
            << vtk_data_norm.str()
            << '\n';

        ofs.close();
    }
};
template<typename T>
struct errors
{
    errors():Iu_uh(0.), Du_Guh(0.)
    {}
    T Iu_uh;
    T Du_Guh;
};


template<typename T, size_t DIM, typename Storage, typename CellQuadType>
class stress_based_errors
{};


template<typename T, typename Storage, typename CellQuadType>
class stress_based_errors<T,1,Storage,CellQuadType>
{
    typedef CellQuadType                                cell_quad_type;
    typedef typename cell_quad_type::quadpoint_type     quadpoint_type;
    typedef typename cell_quad_type::weight_type        weight_type;
    typedef mesh<T,1,Storage>                           mesh_type;
    typedef typename mesh_type::cell                    cell_type;
    typedef typename mesh_type::face                    face_type;
    typedef typename mesh_type::point_type              point_type;
    typedef typename mesh_type::scalar_type             scalar_type;
    typedef dynamic_vector<scalar_type>                 vector_type;
    typedef dynamic_matrix<scalar_type>                 matrix_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;

    size_t  m_degree;
    std::vector<quadpoint_type>     mmmmm;
    plasticity_data<scalar_type>    m_pst;

    public:

        stress_based_errors(const plasticity_data<scalar_type>& pst,const size_t  degree): m_degree(degree), m_pst(pst)
        {};

        template<typename Solution, typename Gradient>
        errors<scalar_type>
        integrate(  const mesh_type &   msh,
                    const cell_type &   cl,
                    const Solution  &   sf,
                    const Gradient  &   df,
                    const vector_type&  x)
        {

            errors<scalar_type> er;
            T  xc = m_pst.Lref / 2.0  -  0.5 * m_pst.Bn * m_pst.Lref * m_pst.f;
            auto pts = points(msh, cl);
            auto lp  = pts[0];
            auto rp  = pts[1];
            auto lyc = point<T,1>({xc});
            auto ryc = point<T,1>({m_pst.Lref -xc});

            if( lp.x() < lyc.x() &  rp.x() > lyc.x() )
            {
                compute(msh,cl,sf,df,lp,lyc,x,er);
                lp = lyc;
            }
            if ( lp.x() < ryc.x()  &   rp.x() > 10.e-15 + ryc.x() )
            {
                compute(msh,cl,sf,df,lp,ryc,x,er);
                lp    = ryc;
            }

            compute(msh,cl,sf,df,lp,rp,x,er);
            return er;
        }
private:
        std::vector<quadpoint_type>
        special_integrate(const point<T,1>& lp,
                          const point<T,1>& rp,
                          const std::vector<std::pair<point<T,1>, T>>& m_quadrature_data)
        {
            auto measure = [](const point<T,1>& left_p,const point<T,1>& right_p)-> auto{
                        assert(pts.size() == 2);
                        return (right_p - left_p).to_vector().norm();
            };

            auto meas       = measure(lp,rp);

            std::vector<quadpoint_type> ret(m_quadrature_data.size());

            auto tr = [&](const std::pair<point<T,1>, T>& qd) -> auto {
                    auto point = (rp - lp) * qd.first.x() + lp;
                    auto weight = qd.second * meas;
                    return make_qp(point, weight);
            };

            auto retbegin = ret.begin();

            std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                            retbegin, tr);
            return ret;

        }
        template<typename Solution, typename Gradient>
        void
        compute(const mesh_type &   msh,
                const cell_type &   cl,
                const Solution  &   solution,
                const Gradient  &   gradient,
                const point<T,1>&   lp,
                const point<T,1>&   rp,
                const vector_type&  x,
                errors<scalar_type>& er)
        {

            typedef typename mesh_type::face             face_type;
            typedef quadrature<mesh_type, face_type>     face_quad_type;

            std::vector<std::pair<point<T,1>, T>> m_quadrature_data;
            m_quadrature_data = edge_quadrature<T>(2*m_degree + 2);

            cell_basis_type         cb1(m_degree);
            cell_basis_type         cb(m_degree+1);

            gradient_reconstruction_nopre<mesh_type,
                                            cell_basis_type,
                                            cell_quad_type,
                                            face_basis_type,
                                            face_quad_type> gradrec_nopre(m_degree);

            matrix_type mass_matrix = matrix_type::Zero(cb1.size(), cb1.size());
            vector_type usource     = vector_type::Zero(cb1.size()) ;

            auto quadpoints =  special_integrate(lp, rp, m_quadrature_data);

            for (auto& qp : quadpoints)
            {
                auto phi = cb1.eval_functions(msh, cl, qp.point());
                for (size_t i = 0; i < cb1.size(); i++)
                    for (size_t j = 0; j < cb1.size(); j++)
                        mass_matrix(i,j) += qp.weight() * mm_prod(phi[i], phi[j]);

                auto ueval    =   solution(qp.point(),m_pst);

                for (size_t i = 0; i < cb1.size(); i++)
                    usource(i)  +=  qp.weight() * ueval  * phi[i];
            }

            //projection on the cell
            Eigen::LDLT<matrix_type> mass_matrix_ldlt = mass_matrix.ldlt();
            vector_type Iu_T = mass_matrix_ldlt.solve(usource);

            for (auto& qp : quadpoints)
            {
                auto phi        =   cb1.eval_functions(msh, cl, qp.point());
                auto dphi       =   cb.eval_gradients( msh, cl, qp.point());
                auto ueval      =   solution(qp.point(), m_pst);
                auto deval      =   gradient(qp.point(), m_pst);
                auto col_range  =   cb.range(1,m_degree+1);
                auto row_range  =   dof_range(0,mesh_type::dimension);

                scalar_type Iu_uh   = 0.0;
                scalar_type uh_T    = 0.0;
                auto xT         =   x.head(cb1.size());
                for (size_t i = 0; i < cb1.size(); i++)
                    Iu_uh += phi[i] * (Iu_T(i) - xT(i));
                for (size_t i = 0; i < cb1.size(); i++)
                    uh_T  += phi[i] * xT(i);

                matrix_type dphi_matrix =   make_gradient_matrix(dphi);
                matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
                matrix_type Guh  =   dphi_taken * gradrec_nopre.oper * x;

                //u_uh    += qp.weight() * ( ueval - uh_T)*(ueval  -  uh_T);
                er.Iu_uh   += qp.weight() *  Iu_uh * Iu_uh;
                er.Du_Guh  += qp.weight() * (deval - Guh).dot(deval - Guh);
            }
        }

    };

    template<typename T, typename Storage, typename CellQuadType>
    class stress_based_errors<T,2,Storage,CellQuadType>
    {
    public:
        typedef mesh<T,2,Storage>                           mesh_type;
        typedef typename mesh_type::cell                    cell_type;
        typedef typename mesh_type::scalar_type             scalar_type;
        typedef dynamic_vector<scalar_type>                 vector_type;

        stress_based_errors(){}
        stress_based_errors(const plasticity_data<scalar_type>& pst,const size_t  degree){}

        template<typename Solution, typename Gradient>
        errors<scalar_type>
        integrate(  const mesh_type &   msh,
                    const cell_type &   cl,
                    const Solution  &   sf,
                    const Gradient  &   df,
                    const vector_type&  x)
        {
            errors<scalar_type> er;
            return er;
        }
    };

    template<typename T, typename Storage, typename CellQuadType>
    class stress_based_errors<T,3,Storage,CellQuadType>
    {
    public:
        typedef mesh<T,3,Storage>                           mesh_type;
        typedef typename mesh_type::cell                    cell_type;
        typedef typename mesh_type::scalar_type             scalar_type;
        typedef dynamic_vector<scalar_type>                 vector_type;

        stress_based_errors(){}
        stress_based_errors(const plasticity_data<scalar_type>& pst,const size_t  degree){}

        template<typename Solution, typename Gradient>
        errors<scalar_type>
        integrate(  const mesh_type &   msh,
                    const cell_type &   cl,
                    const Solution  &   sf,
                    const Gradient  &   df,
                    const vector_type&  x)
        {
            errors<scalar_type> er;
            return er;
        }
    };

template<typename T, size_t DIM, typename Storage, typename TensorsType,
            typename Function, typename Gradient>
void
plasticity_post_processing( mesh<T,DIM,Storage>&  msh,
                            const Function &            sf,
                            const Gradient &            df,
                            const std::vector<TensorsType> &    tsr_vec,
                            std::vector<dynamic_vector<T>> &    Uh_Th,
                            const plasticity_data<T> &  pst,
                            const std::string &         directory,
                            const size_t                degree,
                            std::ofstream &             efs,
                            int iter)
{
    typedef mesh<T, DIM, Storage>                       mesh_type;
    typedef typename mesh_type::scalar_type             scalar_type;
    typedef typename mesh_type::cell                    cell_type;
    typedef quadrature<mesh_type, cell_type>            cell_quad_type;
    typedef dynamic_vector<scalar_type>                 vector_type;
    typedef dynamic_matrix<scalar_type>                 matrix_type;

    errors<scalar_type> er;

    for (auto& cl : msh)
    {
        auto id =   msh.lookup(cl);
        auto x  =   Uh_Th.at(id);
        stress_based_errors<T,DIM,Storage,cell_quad_type>   sber(pst,degree);
        er      =   sber.integrate(msh,cl,sf,df,x);
    }

    //u_uh_error  = std::sqrt(u_uh_error);
    er.Iu_uh    = std::sqrt(er.Iu_uh);
    er.Du_Guh   = std::sqrt(er.Du_Guh);
    std::cout  << "l2-norm error,  Iu_uh:" << er.Iu_uh  << std::endl;
    std::cout  << "l2-norm error, Du_Guh:" << er.Du_Guh << std::endl;
    efs<< iter << " " << er.Iu_uh << " " << er.Du_Guh << std::endl;

    post_processing<T,DIM,Storage> pp;
    std::string constraint_file;

    if (pst.method == true)
        constraint_file = directory + "/lambda";
    else
        constraint_file = directory + "/sigma";

    std::string solution_file  = directory + "/solution";
    std::string gamma_file     = directory + "/gamma";
    std::string xi_norm_file   = directory + "/xi_norm";
    std::string xi_funct_file  = directory + "/xi_function";


    pp.vtk_writer(msh, solution_file ,degree,Uh_Th,iter);

    typedef std::vector<Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>>  dyn_mat_vec;
    dyn_mat_vec     siglam_Th(msh.cells_size()), gamma_Th(msh.cells_size()), xi_norm_Th(msh.cells_size()), xi_function(msh.cells_size());


    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>     tensor_matrix;
    tensor_matrix vec;

    for(auto& cl:msh)
    {
        size_t id  = msh.lookup(cl);
        siglam_Th.at(id)  = tsr_vec.at(id).siglam;
        gamma_Th.at(id)   = tsr_vec.at(id).gamma;
        xi_norm_Th.at(id) = tsr_vec.at(id).xi_norm;

        auto xi_norm_vec  = tsr_vec.at(id).xi_norm;
        vec = tensor_matrix::Zero(1,xi_norm_vec.size());
        for(size_t i = 0; i < xi_norm_vec.size(); i++)
        {
            if (xi_norm_vec(i) >= pst.yield)
                vec(0,i) = 1.;
        }
        xi_function.at(id) = vec;
    }
    pp.tensor_norm_vtk_writer(msh, constraint_file, degree, siglam_Th, iter);
    pp.tensor_norm_vtk_writer(msh, gamma_file,      degree, gamma_Th, iter);
    pp.tensor_norm_vtk_writer(msh, xi_norm_file,    degree, xi_norm_Th, iter);
    pp.tensor_norm_vtk_writer(msh, xi_funct_file,   degree, xi_function, iter);
};


} // namespace disk
