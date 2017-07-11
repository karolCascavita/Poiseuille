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
// #include "visualization/gnuplot-iostream.h"
 #include "hho/hho.hpp"
 #include "hho_pst.hpp"
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
        /*WK: check for degree 0*/
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

    struct sizes
    {
        sizes(){};
        size_t num_nodes;
        size_t num_elems;
        size_t num_scalars;
        size_t num_vectors;
    };

    void paraview(const mesh_type& msh, const std::string& name, const size_t degree,
                    const std::vector<dynamic_vector<T>> & Uh_Th)
    {
    #if 0
        if(degree == 1)
            vtk_second_order(msh, name, degree, Uh_Th, iter);
        else
            vtk_higher_order(msh, name, degree, Uh_Th, iter);
    #endif
    vtk_higher_order(msh, name, degree, Uh_Th);
    }
    #if 0
    void vtk_writer(const size_t degree,
                    const std::string & file_name,
                    const std::stringstream &  vtk_points,
                    const std::stringstream &  vtk_elems,
                    const std::stringstream &  vtk_scalar_data,
                    const std::stringstream &  vtk_vector_data,const sizes& sz)
    {
        size_t el = 0;
        std::vector<size_t> vtk_cell_type = {4,5,10}; /*WK: check this for others VTK types*/
        std::vector<size_t> vtk_nodes_inside_cell = {degree+1,3,4};/* WK: all cells are the same type */
        for (size_t i = 0; i <  sz.num_elems; i++)
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
        ofs << "POINTS " << sz.num_nodes<<" double"<<'\n'
            << vtk_points.str()
            << '\n';

        /* Cells */
        ofs << "CELLS " << sz.num_elems <<' '<<  sz.num_elems *(vtk_nodes_inside_cell[DIM-1] + 1)<< '\n'
            << vtk_elems.str()
            << '\n';

        /* Types of cells*/
        ofs << "CELL_TYPES " << sz.num_elems << '\n';
            for(size_t i = 0; i < sz.num_elems; i++)
                ofs << ' ' << vtk_cell_type[DIM-1];
        ofs << '\n';

        /* Data */
        ofs << "POINT_DATA " << sz.num_nodes
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

    void vtk_second_order(const mesh_type& msh, const std::string& name, const size_t degree,
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

        size_t cl_cont = 0;

        for(auto& cl : msh)
        {
            auto fcs    = faces(msh,cl);
            auto num_fcs = fcs.size();
            auto bar    = barycenter(msh,cl);
            if(num_fcs == 3 | num_fcs == 4)
            {
                auto points = points(msh,cl);
                for(auto& fc : fcs)
                {
                    auto bar = barycenter(msh,fc);
                    points.push_back(bar);
                }
            }
            else
            {
                for( auto& fc : fcs)
                {
                    auto pts = points(msh,fc);
                    auto it  = pts.begin();
                    std::vector<point_type> new_pts(6);
                    new_pts.at(0)  = pts.at(0);
                    new_pts.at(1)  = pts.at(1);
                    new_pts.at(2)  = bar;
                    new_pts.at(3)  = (pts.at(0) +  bar )/2.;
                    new_pts.at(4)  = (pts.at(1) +  bar )/2.;
                    new_pts.at(5)  = (pts.at(0) +  pts.at(1) )/2.;
                }
            }
            #if 0
            gradrec_nopre.compute(msh, cel);
            vector_type     uh_TF =  Uh_Th[cl_cont];
            dynamic_vector<scalar_type> rec(cb.size());
            rec.tail(cb.size()-1) = gradrec_nopre.oper * uh_TF;
            rec(0) = uh_TF(0);
            auto col_range  =   cb.range(1,degree+1);
            auto row_range  =   dof_range(0,DIM);

            for(auto& tp : points)
            {
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

                matrix_type  dphi_matrix     =  make_gradient_matrix(dphi);
                matrix_type  dphi_high_order =  take(dphi_matrix, row_range, col_range);
                vector_type  dpot =   dphi_high_order * rec.tail(col_range.size());

                for (size_t d = 0; d < DIM; d++)
                    vtk_vector_data << dpot(d) << " ";
                for (size_t d = 0; d < 3 - DIM; d++) /* VTK only supports until 3D */
                    vtk_vector_data<< 0. << " ";
                vtk_vector_data<< '\n';
            }
            cl_cont++;
            #endif
        }

        //vtk_writer(vtk_points, vtk_elems, vtk_scalar_data, vtk_vector_data, sz);
    }
#endif
    void vtk_higher_order(const mesh_type& msh,
                            const std::string& name,
                            const size_t degree,
                            const std::vector<dynamic_vector<T>> & Uh_Th)
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


        cell_basis_type             cb(degree+1);
        std::stringstream           vtk_points;
        std::stringstream           vtk_elems;
        std::stringstream           vtk_scalar_data;
        std::stringstream           vtk_vector_data;

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
        //std::string file_name = name + "_i" + numstr;
        std::string file_name = name;

        std::vector<size_t> vtk_cell_type = {4,5,10}; /*WK: check this for others VTK types*/
        std::vector<size_t> vtk_nodes_inside_cell = {degree+1,3,4};/* WK: all cells are the same type */

        //vtk_writer(vtk_points, vtk_elems, vtk_scalar_data, vtk_vector_data, sz);

        /*VTK file*/
        std::string ext = ".vtk";
        std::ofstream  ofs(file_name + ext);
        ofs << "# vtk DataFile Version 3.0"<< std::endl;
        ofs << "#This file was generated by the DISK++ library"<< std::endl;
        ofs << "ASCII"<< std::endl;

        /* Data header */
        ofs << "DATASET UNSTRUCTURED_GRID\n" << std::endl;

        /* Points */
        ofs << "POINTS " << total_sub_nodes<<" double"<<std::endl;

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
                    auto tp   = test_points.at(idx);

                    for (size_t d = 0; d < DIM; d++)
                        ofs << tp.at(d) << " ";
                    for (size_t d = 0; d < 3 - DIM; d++) /* VTK only supports until 3D */
                        ofs<< 0. << " ";
                    ofs<< std::endl;
                }
                ++fc_cont;
            }
            ++cl_cont;
        }
        ofs<<std::endl;

        /* Cells */
        size_t el = 0;
        ofs << "CELLS " << total_sub_elems <<' '<<  total_sub_elems *(vtk_nodes_inside_cell[DIM-1] + 1)<< std::endl;
        for (size_t i = 0; i < total_sub_elems; i++)
        {
            ofs << vtk_nodes_inside_cell[DIM-1] << " ";
            for(size_t j=0; j < vtk_nodes_inside_cell[DIM-1]; j++, el++)
                ofs << el<< " ";
            ofs<< std::endl;
        }

        /* Types of cells*/
        ofs << "CELL_TYPES " << total_sub_elems << std::endl;
            for(size_t i = 0; i < total_sub_elems; i++)
                ofs << ' ' << vtk_cell_type[DIM-1];
        ofs << std::endl;

        /* Data */
        ofs << "POINT_DATA " << total_sub_nodes << std::endl;

        ofs << "VECTORS  Gr_U  double"<< std::endl;

        cl_cont = 0;
        fc_cont = 0;
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
                    auto tp   = test_points.at(idx);
                    auto dphi = cb.eval_gradients(msh, cel, tp);
                    auto col_range  =   cb.range(1,degree+1);
                    auto row_range  =   dof_range(0,DIM);

                    matrix_type dphi_matrix =  make_gradient_matrix(dphi);
                    matrix_type dphi_high_order =   take(dphi_matrix, row_range, col_range);
                    vector_type dpot =   dphi_high_order * rec.tail(col_range.size());
                    for (size_t d = 0; d < DIM; d++)
                        ofs <<  dpot(d) << " ";

                    for (size_t d = 0; d < 3 - DIM; d++) /* VTK only supports until 3D */
                        ofs<< 0. << " ";
                    ofs<< std::endl;
                }
                ++fc_cont;
            }
            ++cl_cont;
        }

        cl_cont = 0;
        fc_cont = 0;

        ofs << "SCALARS  U  double 1"<< std::endl;
        ofs << "LOOKUP_TABLE default"<< std::endl;
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
                    auto tp   = test_points.at(idx);
                    auto phi  = cb.eval_functions(msh, cel, tp);
                    for (size_t i = 0; i < cb.size(); i++)
                        pot  +=  phi[i] * rec(i);
                    ofs << pot <<' ';
                }
                ++fc_cont;
            }
            ++cl_cont;
        }
        ofs.close();
    }


    template<typename TensorMatrix>
    void tensor_norm_vtk_writer(const mesh_type& msh,
                const std::string& name,
                const size_t quad_degree,
                const std::vector<TensorMatrix> & tensor)   /* storage of velocities*/
    {

        typedef quadrature<mesh_type, cell_type>     cell_quadrature_type;
        typedef quadrature<mesh_type, face_type>     face_quadrature_type;

        typedef scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
        typedef scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

        cell_quadrature_type    cq(quad_degree);

        size_t total_num_pts  = 0;
        std::string result;
        //std::string file_name = name + "_i" + numstr;
        std::string file_name = name ;

        for(auto& cel : msh)
        {
            auto pts = points(msh, cel);
            auto cq_pts  =  cq.integrate(msh, cel);
            total_num_pts +=  cq_pts.size() + pts.size();
        }

        /*VTK file*/
        std::string ext = ".vtk";
        std::ofstream  ofs(file_name + ext);
        ofs << "# vtk DataFile Version 3.0"<< std::endl;
        ofs << "#This file was generated by the DISK++ library"<< std::endl;
        ofs << "ASCII"<< std::endl;

        /* Data header */
        ofs << "DATASET UNSTRUCTURED_GRID"<<std::endl;
        ofs << std::endl;

        /* Points */
        ofs << "POINTS " << total_num_pts <<" double"<<std::endl;
        for(auto& cel : msh)
        {
            size_t qp_cont = 0;
            auto cell_id   = cel.get_id();
            auto cell_quadpoints  =  cq.integrate(msh, cel);
            auto cqs  = cell_quadpoints.size();
            auto pts  = points(msh,cel);
            auto ps   = pts.size();
            for (auto& qp : cell_quadpoints)
            {
                auto tp  = qp.point();

                for (size_t d = 0; d < DIM; d++)
                    ofs << tp[d] << " ";
                for (size_t d = 0; d < 3 - DIM; d++) /* VTK only supports until 3D */
                    ofs<< 0. << " ";
                ofs<< std::endl;
            }
            for (size_t i = 0; i < ps; i++)
            {
                auto tp = pts.at(i);
                for (size_t d = 0; d < DIM; d++)
                    ofs << tp[d] << " ";
                for (size_t d = 0; d < 3 - DIM; d++) /* VTK only supports until 3D */
                    ofs<< 0. << " ";
                ofs<< std::endl;
            }
        }

        /* Cells */
        ofs << "CELLS " << total_num_pts <<' '<<  total_num_pts *(1 + 1)<<std::endl;
        for (size_t i = 0; i < total_num_pts; i++)
            ofs << 1 <<" "<< i << std::endl;
        ofs <<  std::endl;

        /* Types of cells*/
        ofs << "CELL_TYPES " << total_num_pts <<std::endl;
        for(size_t i = 0; i < total_num_pts; i++)
            ofs << " " << 1;
        ofs <<std::endl;

        /* Data */
        ofs << "POINT_DATA " << total_num_pts<<std::endl;
        ofs << std::endl;

        ofs << "SCALARS "<<std::endl;
        ofs << " tensor_norm "<<std::endl;
        ofs << " double 1"<<std::endl;
        ofs << std::endl;
        ofs << "LOOKUP_TABLE default"<<std::endl;
        ofs <<std::endl;

        for(auto& cel : msh)
        {
            size_t qp_cont = 0;
            auto cell_id   = cel.get_id();
            auto cell_quadpoints     =  cq.integrate(msh, cel);
            auto tsr =  tensor.at(cell_id);

            for (auto& qp : cell_quadpoints)
            {
                vector_type pot = tsr.cell_quad_pts.col(qp_cont);
                if( DIM == 1)
                    ofs<< pot << " ";
                if (DIM == 2)
                    ofs << pot.norm() << " ";
                ++qp_cont;
            }
            auto pts  = points(msh,cel);
            for (size_t i = 0; i < pts.size(); i++)
            {
                vector_type pot = tsr.grid_points.col(i);
                if( DIM == 1)
                    ofs << pot << " ";
                if (DIM == 2)
                    ofs << pot.norm() << " ";
            }
        }
        ofs.close();
    }
    template<typename PunctualTensor>
    void
    tensor_points_to_matlab(const mesh_type& msh, const std::string& name, const size_t quad_degree,
               const std::vector<PunctualTensor> & tensor, const std::string& other_info)   /* storage of velocities*/
    {

        typedef quadrature<mesh_type, cell_type>     cell_quadrature_type;
        typedef quadrature<mesh_type, face_type>     face_quadrature_type;
        cell_quadrature_type    cq(quad_degree);
        face_quadrature_type    fq(quad_degree);
        std::stringstream       vtk_points, vtk_elems, vtk_data_x, vtk_data_y, vtk_data_norm;


        std::ofstream cofs(name);
        cofs<< "V = cell("<< msh.cells_size()<<", 2);"<< std::endl;

        size_t i = 0;
        if (!cofs.is_open())
            std::cout << "Error opening cofs"<<std::endl;

        for(auto cell : msh)
        {
            auto cell_id = cell.get_id();
            auto tsr =  tensor[cell_id];

            cofs<<" PTS = [";
            auto pts =  tsr.tensor_points( msh, cell);
            for(auto& p: pts )
            {
                for (size_t d = 0; d < DIM; d++)
                    cofs << p[d] << " ";
                cofs<< std::endl;
            }
            cofs<<" ]"<< std::endl;


            matrix_type tensor_vec = tsr.join_all();

            cofs<<" TSR = [";
            for(size_t i = 0; i < tensor_vec.cols(); i++)
            {

                if( DIM == 1) //Temrinar esto
                    // cofs << tp[d] << " ";
                if (DIM == 2)
                {
                    vector_type tp = tensor_vec.col(i);
                    cofs << tp.norm() << std::endl;
                }

                cofs << std::endl;
            }
            cofs<<" ]"<< std::endl;
            cofs<< " V{"<<cell_id + 1 <<",1} = PTS;"<<std::endl;
            cofs<< " V{"<<cell_id + 1 <<",2} = TSR;"<<std::endl;
        }
        cofs<< "X = ["<<std::endl;
        for (auto& cl : msh)
        {
            auto bar = barycenter(msh,cl);
            cofs<< i+1 << " " << bar.x()<< " "<< bar.y()<<std::endl;
            i++;
        }
        cofs<< " ];"<<std::endl;
        cofs<< "info = 'mesh"<< other_info<<"'"<<std::endl;
        cofs<< "mesh" << other_info<<std::endl;
        cofs<< "strValues = strtrim(cellstr(num2str([X(:,1)],'(%d)')));"<<std::endl;
        cofs<< "text(X(:,2),X(:,3),strValues,'VerticalAlignment','bottom');"<<std::endl;
        cofs<< "hold on; plot(X(:,2),X(:,3),'o')"<< std::endl;
        cofs.close();
    }
};

template<typename MeshType>
void
quiver_matlab(const MeshType& msh,
              const std::string& filename,
              const size_t quad_degree,
              const std::vector<punctual_tensor<MeshType>> & global_tensor)
{
    typedef MeshType                           mesh_type;
    typedef typename mesh_type::scalar_type             scalar_type;
    typedef typename mesh_type::cell                    cell_type;
    typedef typename mesh_type::face                    face_type;
    typedef disk::quadrature<mesh_type, cell_type>      cell_quad_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quad_type;

    cell_quad_type cell_quadrature(quad_degree);
    face_quad_type face_quadrature(quad_degree);

    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file :"<< filename<<std::endl;

    ofs << " hold on;"<<std::endl;
    for (size_t i = 0; i < msh.cells_size(); i++ )
    {
        auto cl  = *std::next(msh.cells_begin(),i);
        punctual_tensor<MeshType> tensor = global_tensor.at(i);
        auto pts_to_eval = tensor.tensor_points(msh, cl);

        ofs << "pts = [";
        for (auto& p : pts_to_eval)
            ofs << p.x() << "  "<< p.y() << std::endl;
        ofs << "];"<<std::endl;

        dynamic_matrix<scalar_type> tsr_all = tensor.join_all();

        ofs << "tau = [";
        for(size_t i = 0; i < tsr_all.cols(); i++)
            ofs<< tsr_all(0, i)<< " "<< tsr_all(1, i) <<std::endl;
        ofs << "];"<<std::endl;

        ofs << "quiver(pts(:,1),pts(:,2),tau(:,1),tau(:,2));"<< std::endl;
    }
    ofs.close();
}


template<typename T>
struct solution_errors
{
    solution_errors():u_uh(0.), Iu_uh(0.), Du_Guh(0.)
    {}
    T  u_uh;
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

        template<typename Solution, typename Gradient>//,typename LoaderType>
        void
        integrate(  const mesh_type &   msh,
                    const cell_type &   cl,
                    const Solution  &   sf,
                    const Gradient  &   df,
                    const vector_type&  x,
                    solution_errors<scalar_type>& er)
        {
            T  xc = m_pst.Lref / 2.0  -  m_pst.yield * m_pst.f;
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
        }
private:
        std::vector<quadpoint_type>
        special_integrate(const point<T,1>& lp,
                          const point<T,1>& rp,
                          const std::vector<std::pair<point<T,1>, T>>& m_quadrature_data)
        {
            auto measure = [](const point<T,1>& left_p,const point<T,1>& right_p)-> auto{
                        //assert(pts.size() == 2);
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
                solution_errors<scalar_type>& er)
        {

            typedef typename mesh_type::face             face_type;
            typedef quadrature<mesh_type, face_type>     face_quad_type;

            std::vector<std::pair<point<T,1>, T>> m_quadrature_data;
            m_quadrature_data = edge_quadrature<T>(2*m_degree + 4); // why "+ 4"?

            cell_basis_type         cb1(m_degree);
            cell_basis_type         cb(m_degree+1);

            gradient_reconstruction_nopre<mesh_type,
                                            cell_basis_type,
                                            cell_quad_type,
                                            face_basis_type,
                                            face_quad_type> gradrec_nopre(m_degree);

            matrix_type mass_matrix = matrix_type::Zero(cb1.size(), cb1.size());
            vector_type u_rhs       = vector_type::Zero(cb1.size()) ;

            auto quadpoints =  special_integrate(lp, rp, m_quadrature_data);

            for (auto& qp : quadpoints)
            {
                auto phi = cb1.eval_functions(msh, cl, qp.point());
                for (size_t i = 0; i < phi.size(); i++)
                    for (size_t j = 0; j < phi.size(); j++)
                        mass_matrix(i,j) += qp.weight() * mm_prod(phi[i], phi[j]);

                auto ueval    =   solution(qp.point());

                for (size_t i = 0; i < phi.size(); i++)
                    u_rhs(i)  +=  qp.weight() * ueval  * phi[i];
            }

            //projection on the cell
            Eigen::LDLT<matrix_type> mass_matrix_ldlt = mass_matrix.ldlt();
            vector_type Iu_T = mass_matrix_ldlt.solve(u_rhs);

            auto Iu    = mass_matrix_ldlt.solve(u_rhs);
            auto xT    = x.head(cb1.size());
            auto Iu_uh = Iu_T - xT;
            //std::cout <<"/* ____cell = "<< msh.lookup(cl)<<"____*/"<<std::endl;
            //std::cout << "lp = "<< lp.x()<< ";  rp = "<<  rp.x() << " ; xT = "<< xT.transpose() << std::endl;
            //std::cout << "Iu_uh = "<< Iu_uh<< std::endl;
            //std::cout << "before er.Iu_uh"<<er.Iu_uh << std::endl;

            er.Iu_uh  += Iu_uh.transpose() * mass_matrix * Iu_uh;
            //std::cout << "cell_id: " <<msh.lookup(cl);

            for (auto& qp : quadpoints)
            {

                auto phi        =   cb1.eval_functions(msh, cl, qp.point());
                auto dphi       =   cb.eval_gradients( msh, cl, qp.point());
                auto ueval      =   solution(qp.point());
                vector_type deval = gradient(qp.point());

                auto col_range  =   cb.range(1,m_degree+1);
                auto row_range  =   dof_range(0,mesh_type::dimension);

                scalar_type uh_T    = 0.0;
                auto xT         =   x.head(cb1.size());

                for (size_t i = 0; i < cb1.size(); i++)
                    uh_T  += phi[i] * xT(i);

                gradrec_nopre.compute(msh, cl);

                matrix_type dphi_matrix =   make_gradient_matrix(dphi);
                matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
                vector_type Guh  =   dphi_taken * gradrec_nopre.oper * x;

                er.u_uh    += qp.weight() * (ueval - uh_T) * (ueval - uh_T);
                er.Du_Guh  += qp.weight() * (deval - Guh).dot(deval - Guh);
            }
            //for (size_t i = 0; i < cb1.size(); i++)
            //    Iu_uh += phi[i] * (Iu_T(i) - xT(i));
        }
    };

    template<typename T, typename Storage, typename CellQuadType>
    class stress_based_errors<T,2,Storage,CellQuadType>
    {
    public:
        typedef mesh<T,2,Storage>                           mesh_type;
        typedef typename mesh_type::cell                    cell_type;
        typedef typename mesh_type::face                    face_type;
        typedef typename mesh_type::point_type              point_type;
        typedef typename mesh_type::scalar_type             scalar_type;
        typedef dynamic_vector<scalar_type>                 vector_type;
        typedef dynamic_matrix<scalar_type>                 matrix_type;
        typedef CellQuadType                                cell_quad_type;
        typedef scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
        typedef scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

        template<size_t N> using polygon = std::array<ident_impl_t, N>;

        std::vector<point_type>                         m_points;
        std::vector<polygon<3>>                         m_triangles;
        std::vector<polygon<4>>                         m_quadrangles;
        std::vector<polygon<5>>                         m_pentagons;
        std::vector<polygon<6>>                         m_hexagons;
        std::vector<polygon<7>>                         m_ennagons;
        std::vector<std::array<ident_impl_t, 2>>        m_boundary_edges;
        std::vector<std::array<ident_impl_t, 4>>        m_edges;

        size_t  m_degree;
        plasticity_data<scalar_type>    m_pst;

        stress_based_errors(const plasticity_data<scalar_type>& pst,const size_t  degree): m_degree(degree), m_pst(pst)
        {};

        bool
        cl_mark(const mesh_type & msh,const cell_type & cl)
        {
            bool mark = false;
            auto r_yield  = m_pst.Bn*m_pst.Lref;
            auto pts      = points(msh,cl);
            std::vector<T> r(pts.size());
            for(size_t i = 0; i < pts.size(); i++)
                r.at(i)  = std::sqrt(pts[i].x() * pts[i].x() + pts[i].y() * pts[i].y());
            auto result = std::minmax_element(r.begin(), r.end());

            if(*(result.second) >  r_yield)
                if (*(result.first) < r_yield)
                    mark = true;

            std::cout << "r = "<< std::endl;
            for(auto& mm : r)
                std::cout << mm<<" ";
            std::cout << std::endl;
            std::cout << "mark = "<< mark << std::endl;
            return mark;
        }
        void
        store(const std::vector<size_t> polygon_x)
        {
            size_t N = polygon_x.size();
            if(N == 3)
            {
                std::array<ident_impl_t, 3> pol;
                for(size_t i = 0; i < N; i++)
                    pol[i] = polygon_x.at(i);
                m_triangles.push_back(pol);
            }
            if(N == 4)
            {
                std::array<ident_impl_t , 4> pol;
                for(size_t i = 0; i < N; i++)
                    pol[i] = polygon_x.at(i);
                m_quadrangles.push_back(pol);
            }
            if(N == 5)
            {
                std::array<ident_impl_t, 5> pol;
                for(size_t i = 0; i < N; i++)
                    pol[i] = polygon_x.at(i);
                m_pentagons.push_back(pol);
            }
            if(N == 6)
            {
                std::array<ident_impl_t, 6> pol;
                for(size_t i = 0; i < N; i++)
                    pol[i] = polygon_x.at(i);
                m_hexagons.push_back(pol);
            }
            if(N == 7)
            {
                std::array<ident_impl_t, 7> pol;
                for(size_t i = 0; i < N; i++)
                    pol[i] = polygon_x.at(i);
                m_ennagons.push_back(pol);
            }

            for(size_t i = 0; i < N - 1; i++)
            {
                auto a = polygon_x.at(i);
                auto b = polygon_x.at(i + 1);
                if(b < a)
                    std::swap(a,b);
                m_edges.push_back({a, b});
            }
        }
        //template<typename LoaderType>
        std::vector<cell_type>
        cl_partition(const mesh_type & msh,const cell_type & cl)
        {
            typedef generic_mesh<T,2>                       mesh_type;
            typedef typename mesh_type::scalar_type         scalar_type;
            typedef typename mesh_type::point_type          point_type;
            typedef typename mesh_type::node_type           node_type;
            typedef typename mesh_type::edge_type           edge_type;
            typedef typename mesh_type::surface_type        surface_type;
            //typedef LoaderType loader_type;
            typedef fvca5_mesh_loader<scalar_type, 2>       loader_type;

            std::vector<cell_type>      vec;
            std::vector<size_t>         polygon_1,  polygon_2, *polygon_ptr;

            polygon_ptr  = &polygon_1;
            auto fcs     = faces(msh, cl);
            auto fcs_ids = cl.faces_ids();
            auto num_fcs = fcs_ids.size();
            auto cl_pts  = points(msh,cl);
            auto i = 0, j = 0;
            auto r_yield  = m_pst.Bn*m_pst.Lref;
            std::cout << "r_yield = "<< r_yield << std::endl;
            std::cout << "cl      = "<< cl.get_id() << std::endl;
            std::cout << "cl_pts.size = "<< cl_pts.size() << std::endl;

            for(auto& fc : fcs)
            {
                auto pts = points(msh, fc);
                auto r1  = std::sqrt(pts[0].x() * pts[0].x() + pts[0].y() * pts[0].y());
                auto r2  = std::sqrt(pts[1].x() * pts[1].x() + pts[1].y() * pts[1].y());
                polygon_ptr->push_back(j + i);
                m_points.push_back(cl_pts.at(i));

                std::cout << "face = "<< msh.lookup(fc) << std::endl;
                std::cout << "fc_pts.size = "<< pts.size() << std::endl;
                std::cout << "i  = "<<i <<"   j  = "<< j<< std::endl;
                std::cout << "r1 = "<< r1<<"  r2 = "<< r2 << std::endl;
                if(std::max(r1,r2) >  r_yield & std::min(r1,r2) < r_yield)
                {
                    std::cout << "inside: i = "<<i <<"  j = "<< j<< std::endl;
                    j++;
                    auto bar    = barycenter(msh,fc);
                    polygon_ptr->push_back(j + i);
                    polygon_1.swap(polygon_2);
                    polygon_ptr->push_back(j + i);
                    m_points.push_back(bar);
                }
                i++;
            }
            std::cout << "/* message */" << std::endl;
            store(polygon_1);
            store(polygon_2);
            std::sort(m_edges.begin(), m_edges.end());
            auto uniq_iter = std::unique(m_edges.begin(), m_edges.end());
            m_edges.erase(uniq_iter, m_edges.end());

            loader_type      loader;
            loader.m_edges        = m_edges;
            loader.m_triangles    = m_triangles;
            loader.m_quadrangles  = m_quadrangles;
            loader.m_pentagons    = m_pentagons;
            loader.m_hexagons     = m_hexagons;
            loader.m_boundary_edges = m_boundary_edges;
            loader.m_ennagons   = m_ennagons;
            loader.m_points     = m_points;

            mesh_type temp_msh;
            loader.populate_mesh(temp_msh);
            size_t k = 0;
            for(auto& cel : temp_msh)
            {
                vec.at(k) = cel;
                k++;
            }
            return vec;
        }
        template<typename Solution, typename Gradient> //,typename LoaderType>
        void
        integrate(  const mesh_type &   msh,
                    const cell_type &   cl,
                    const Solution  &   sf,
                    const Gradient  &   df,
                    const vector_type&  x,
                    solution_errors<scalar_type>& er)
        {
            typedef point<T,2>          point_type;

            #if 0
            auto r_yield = m_pst.Bn*m_pst.Lref;

            std::cout << "cl_mark(msh,"<< msh.lookup(cl)<<") = "<< cl_mark(msh,cl)<< std::endl;
            if(cl_mark(msh,cl))
            {
                //auto cl_vector = cl_partition<LoaderType>(msh, cl);
                auto cl_vector = cl_partition(msh, cl);
                compute(msh, cl, cl_vector.at(0), sf, df, x, er);
                compute(msh, cl, cl_vector.at(1), sf, df, x, er);
            }
            #endif

            compute(msh, cl, cl, sf, df, x, er);
        }

        template<typename Solution, typename Gradient>
        void
        compute(const mesh_type &   msh,
                const cell_type &   cl,
                const cell_type &   sub_cl,
                const Solution  &   solution,
                const Gradient  &   gradient,
                const vector_type&  x,
                solution_errors<scalar_type>& er)
        {

            typedef typename mesh_type::face             face_type;
            typedef quadrature<mesh_type, face_type>     face_quad_type;

            cell_quad_type          cq(2*m_degree + 4);
            cell_basis_type         cb1(m_degree);
            cell_basis_type         cb(m_degree+1);

            gradient_reconstruction_nopre<mesh_type,
                                            cell_basis_type,
                                            cell_quad_type,
                                            face_basis_type,
                                            face_quad_type> gradrec_nopre(m_degree);

            matrix_type mass_matrix = matrix_type::Zero(cb1.size(), cb1.size());
            vector_type u_rhs       = vector_type::Zero(cb1.size()) ;

            gradrec_nopre.compute(msh, cl);

            auto quadpoints = cq.integrate(msh, sub_cl);


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






            for (auto& qp : quadpoints)
            {
                auto phi        =   cb1.eval_functions(msh, cl, qp.point());
                auto dphi       =   cb.eval_gradients( msh, cl, qp.point());
                auto ueval      =   solution(qp.point(), number);
                vector_type deval = gradient(qp.point(), number);

                auto col_range  =   cb.range(1,m_degree+1);
                auto row_range  =   dof_range(0,mesh_type::dimension);

                scalar_type uh_T    = 0.0;
                auto xT         =   x.head(cb1.size());

                for (size_t i = 0; i < cb1.size(); i++)
                    uh_T  += phi[i] * xT(i);

                matrix_type dphi_matrix =   make_gradient_matrix(dphi);
                matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
                vector_type Guh  =   dphi_taken * gradrec_nopre.oper * x;

                er.u_uh    += qp.weight() * (ueval - uh_T) * (ueval - uh_T);
                er.Du_Guh  += qp.weight() * (deval - Guh).dot(deval - Guh);
            }
            //for (size_t i = 0; i < cb1.size(); i++)
            //    Iu_uh += phi[i] * (Iu_T(i) - xT(i));
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
        stress_based_errors(const plasticity_data<scalar_type>& pst,const size_t  degree)
        {}

        template<typename Solution, typename Gradient>
        void
        integrate(  const mesh_type &   msh,
                    const cell_type &   cl,
                    const Solution  &   sf,
                    const Gradient  &   df,
                    const vector_type&  x,
                    solution_errors<scalar_type>& er)
        {
            int a;
        }
    };

template<typename T, size_t DIM, typename Storage, typename TensorsType,
            typename Solution>
solution_errors<T>
postprocess(const  mesh<T,DIM,Storage>&  msh,
                            const Solution &     solution,
                            const std::vector<TensorsType> &    tsr_vec,
                            const std::vector<dynamic_vector<T>> &    Uh_Th,
                            const plasticity_data<T> &  pst,
                            const mesh_parameters<T>  &  mp,
                            const size_t                degree,
                            //std::ofstream &             efs,
                            const size_t                imsh )

{


    std::cout << "/* plasticity_post_processing */" << std::endl;

    typedef mesh<T, DIM, Storage>                       mesh_type;
    typedef typename mesh_type::scalar_type             scalar_type;
    typedef typename mesh_type::cell                    cell_type;
    typedef quadrature<mesh_type, cell_type>            cell_quad_type;
    typedef dynamic_vector<scalar_type>                 vector_type;
    typedef dynamic_matrix<scalar_type>                 matrix_type;


    #if 0
    for(auto& cl:msh)
    {

        size_t id  = cl.get_id();

        std::cout << "tensor ("<<id <<")" << std::endl;
        std::cout << "SIGMA: "<<std::endl;
        std::cout << tsr_vec.at(id).sigma << std::endl;

        std::cout << "GAMMA: "<<std::endl;
        std::cout << tsr_vec.at(id).gamma << std::endl;
    }
    #endif

    solution_errors<scalar_type> er;

    std::string  msh_str = tostr(imsh);


    for (auto& cl : msh)
    {
        auto id =   cl.get_id();
        auto x  =   Uh_Th.at(id);
        stress_based_errors<T,DIM,Storage,cell_quad_type>   sber(pst,degree);
        //sber.integrate<LoaderType>(msh,cl,sf,df,x,er);
        sber.integrate(msh, cl, solution.sf, solution.df, x, er);
    }

    er.u_uh     = std::sqrt(er.u_uh);
    er.Iu_uh    = std::sqrt(er.Iu_uh);
    er.Du_Guh   = std::sqrt(er.Du_Guh);
//            auto vts       =  msh.get_vertices(cl, pts);
//            auto h_T       =  diameter(msh, vts);

//    std::cout << "Mesh diameter: " << diam << std::endl;
    #if 0
    std::cout  << "l2-norm error,   u_uh:" << er.u_uh  << std::endl;
    std::cout  << "l2-norm error,  Iu_uh:" << er.Iu_uh  << std::endl;
    std::cout  << "l2-norm error, Du_Guh:" << er.Du_Guh << std::endl;
    //efs<<  er.u_uh  <<" "<< er.Iu_uh << " " << er.Du_Guh ;
    #endif
    post_processing<T,DIM,Storage> pp;

    auto info =  mp.summary + "_R" + msh_str;

    std::string solution_file  = mp.directory + "/solution" + info;
    std::string xi_norm_file   = mp.directory + "/xi_norm"  + info;
    std::string xi_funct_file  = mp.directory + "/xi_function" + info;
    std::string gauss_pts_file = mp.directory + "/cell_gp"  + info;
    std::string quiver_gamma_file = mp.directory + "/quiver_gamma" + info;
    std::string quiver_sigma_file = mp.directory + "/quiver_sigma" + info;
    std::string quiver_sigma_agamma_file = mp.directory + "/quiver_sigma_agamma" + info;

    pp.paraview(msh, solution_file, degree, Uh_Th);

    typedef std::vector<punctual_tensor<mesh_type>>  dyn_mat_vec;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>     tensor_matrix;

    dyn_mat_vec     sigma_Th(msh.cells_size());
    dyn_mat_vec     gamma_Th(msh.cells_size());
    dyn_mat_vec     xi_norm_Th(msh.cells_size());
    dyn_mat_vec     xi_function(msh.cells_size());
    dyn_mat_vec     diff_sig_agam(msh.cells_size());
    std::vector<size_t>  quadeg(msh.cells_size());
    tensor_matrix   vec;

    get_from_tensor(sigma_Th,  tsr_vec, "sigma");
    get_from_tensor(gamma_Th,  tsr_vec, "gamma");
    get_from_tensor(xi_norm_Th,  tsr_vec, "xi_norm");
    get_from_tensor(quadeg,  tsr_vec, "quad_degree");

    //WK: for p-adaptation quad_degree is not the same for all cells
    size_t quad_degree = tsr_vec.at(0).sigma.m_quad_degree;

    for(auto& cl:msh)
    {

        size_t id  = cl.get_id();

        diff_sig_agam.at(id) = tsr_vec.at(id).sigma - pst.alpha * tsr_vec.at(id).gamma;

        auto xi_norm  = tsr_vec.at(id).xi_norm;

        vec = tensor_matrix::Zero(1,xi_norm.cols());

        matrix_type xi_norm_vec = xi_norm.join_all();
        for(size_t i = 0; i < xi_norm_vec.cols(); i++)
        {
            if (xi_norm_vec(i) >= pst.yield)
                vec(0,i) = 1.;
        }
        xi_function.at(id) = punctual_tensor<mesh_type>(msh, cl, quad_degree);
        xi_function.at(id).split_all(vec, 1);
    }
    std::cout << "QUAD DEGREE : "<< quad_degree << std::endl;
    pp.tensor_norm_vtk_writer(msh, xi_norm_file,    quad_degree, xi_norm_Th);
    pp.tensor_norm_vtk_writer(msh, xi_funct_file,   quad_degree, xi_function);
    pp.tensor_norm_vtk_writer(msh, mp.directory +"/sigma"+ info, quad_degree, sigma_Th);
    pp.tensor_norm_vtk_writer(msh, mp.directory +"/gamma"+ info, quad_degree, gamma_Th);
    pp.tensor_norm_vtk_writer(msh, mp.directory +"/sig_agam"+ info, quad_degree, diff_sig_agam);

    pp.tensor_points_to_matlab(msh, gauss_pts_file  + ".m",   quad_degree, xi_function, info);
    std::cout << "QUIVER SIGMA" << std::endl;
    quiver_matlab(msh, quiver_sigma_file + ".m", quad_degree, sigma_Th);
    std::cout << "QUIVER GAMMA" << std::endl;
    quiver_matlab(msh, quiver_gamma_file + ".m", quad_degree, gamma_Th);
    std::cout << "QUIVER SIGMA - A * GAMMA" << std::endl;
    quiver_matlab(msh, quiver_sigma_agamma_file + ".m", quad_degree, diff_sig_agam);

    //auto xi_norm_vec = get_tensor<TensorsType,T>(tsr_vec, "xi_norm");

    //pp.quiver_matlab(msh, name, quad_degree, vec);

    save_data(Uh_Th, mp.directory + "/Uh" + info + ".txt");
    save_data(quadeg,   mp.directory + "/QuadDegree" + info +".txt");
    save_data(sigma_Th, mp.directory + "/Sigma" + info + ".txt");


    auto exfunc_name = mp.directory + "/exFunction.m";

    std::ofstream effs(exfunc_name);
    if (!effs.is_open())
        std::cout << "Error opening file :"<<exfunc_name <<std::endl;

    effs << "col = 2"<<std::endl;
    effs << "Rconvergence = false"<<std::endl;
    effs << "if Rconvergence ==true"<<std::endl;
    effs << "   col = 5"<<std::endl;
    effs << "end"<<std::endl;

    effs << "draw_mesh  = true"<<std::endl;
    effs << "print_mesh = false"<<std::endl;

    effs << "draw_error = true"<<std::endl;
    effs << "print_error= false"<<std::endl;

    effs << "draw_estimator  = true"<<std::endl;
    effs << "print_estimator = false"<<std::endl;
    effs << "error_name = 'nothing'"<<std::endl;

    effs << "execute" + mp.summary<<std::endl;
    effs << "if print_error == true "<<std::endl;
    effs << "   name = strcat(error_name,'"<<mp.summary + ".png')"<<std::endl;
    effs << "   legend('-DynamicLegend')" <<std::endl;
    effs << "   set(gca,'box','on'); set(gcf,'color','w');"<<std::endl;
    effs << "   set(gca,'fontsize',12); set(gca,'xgrid', 'on', 'ygrid', 'on');"<<std::endl;
    effs << "   print -dpng name" <<std::endl;
    effs << "end"<<std::endl;
    effs.close();
    return er;
};

template<typename MeshParameters, typename PlasticData>
void
execute(std::ofstream & exfs,
        const MeshParameters& mp,
        const PlasticData& pst,
        const size_t imsh)
{
    if(imsh > 0)
    {
        exfs<< "if draw_estimator == true"<<std::endl;
        exfs<< "figure;estimator_tot"<<mp.summary<<"_RC"<<imsh<<";";
        exfs<< "if print_estimator == true;";
        exfs<< "set(gca,'box','on'); set(gcf,'color','w');set(gca,'fontsize',12);"<<std::endl;
        exfs<< "circle(0,0,"<<pst.Bn<<",'w')"<<std::endl;
        exfs<< "print('-depsc2', '-loose', 'estimator_tot"<<mp.summary<<"_RC"<<imsh<<".eps')"<<std::endl;
        exfs<< "figure;estimator_res"<<mp.summary<<"_RC"<<imsh<<";";
        exfs<< "if print_estimator == true;";
        exfs<< "set(gca,'box','on'); set(gcf,'color','w');set(gca,'fontsize',12);"<<std::endl;
        exfs<< "circle(0,0,"<<pst.Bn<<",'w')"<<std::endl;
        exfs<< "print('-depsc2', '-loose', 'estimator_res"<<mp.summary<<"_RC"<<imsh<<".eps')"<<std::endl;
        exfs<< "figure;estimator_str"<<mp.summary<<"_RC"<<imsh<<";";
        exfs<< "if print_estimator == true;";
        exfs<< "set(gca,'box','on'); set(gcf,'color','w');set(gca,'fontsize',12);"<<std::endl;
        exfs<< "circle(0,0,"<<pst.Bn<<",'w')"<<std::endl;
        exfs<< "print('-depsc2', '-loose', 'estimator_str"<<mp.summary<<"_RC"<<imsh<<".eps')"<<std::endl;
        exfs<< " end;"<<std::endl;
        exfs<< "end"<<std::endl;
    }

    exfs<< "if draw_error == true "<<std::endl;
    exfs<< "e"<< imsh <<"= load('error"<<mp.summary<<"_R"<<imsh<<".dat');"<<std::endl;
    exfs<< "figure(100); hold on; loglog(e"<<imsh<<"(:,1),e"<<imsh<<"(:,col), 'DisplayName', 'step "<<imsh<<"' , 'color', rand(1,3));"<<std::endl;
    exfs<< "set(gca,'box','on'); set(gcf,'color','w');"<<std::endl;
    exfs<< "set(gca,'fontsize',12); set(gca,'xgrid', 'on', 'ygrid', 'on');"<<std::endl;
    exfs<< "set(gca, 'xscale', 'log')"<<std::endl;
    exfs<< "set(gca, 'yscale', 'log')"<<std::endl;
    exfs<< "end" << std::endl;

    exfs<< "if draw_mesh == true "<<std::endl;
    exfs<< "figure; mesh"<<mp.summary<<"_R" <<imsh<<";";
    exfs<< "if print_mesh == true;";
    exfs<< "set(gca,'box','on'); set(gcf,'color','w');set(gca,'fontsize',12);"<<std::endl;
    exfs<< "circle(0,0,"<<pst.Bn<<",'r')"<<std::endl;
    exfs<< "print('-depsc2', '-loose', 'mesh_B"<<tostr(10*pst.Bn)<<mp.summary<<"_R" <<imsh<<".eps')"<<std::endl;
    exfs << "end"<<std::endl;
    exfs << "end"<<std::endl;
}

#if 0
template<typename T, size_t DIM, typename Storage, typename TensorsType,
            typename Function, typename Gradient>
void
gnu_plot( mesh<T,DIM,Storage>&  msh,
                            const Function &            sf,
                            const Gradient &            df,
                            std::vector<dynamic_vector<T>> &    Uh_Th,
                            const size_t                degree)
{
    typedef mesh<T, DIM, Storage>                       mesh_type;
    typedef typename mesh_type::scalar_type             scalar_type;
    typedef typename mesh_type::cell                    cell_type;
    typedef quadrature<mesh_type, cell_type>            cell_quad_type;
    typedef dynamic_vector<scalar_type>                 vector_type;
    typedef dynamic_matrix<scalar_type>                 matrix_type;

    /*****************************************************************/
    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    cell_quad_type      cq(2*degree );
    cell_basis_type     cb1(degree);

    std::vector<dynamic_vector<T>> Iu_Th(Uh_Th.size());
    size_t i = 0;

    for(auto& cl: msh)
    {
        matrix_type mass_matrix = matrix_type::Zero(cb1.size(), cb1.size());
        vector_type u_rhs       = vector_type::Zero(cb1.size()) ;

        auto quadpoints =  cq.integrate(msh,cl);

        for (auto& qp : quadpoints)
        {
            auto phi = cb1.eval_functions(msh, cl, qp.point());
            for (size_t i = 0; i < cb1.size(); i++)
                for (size_t j = 0; j < cb1.size(); j++)
                    mass_matrix(i,j) += qp.weight() * mm_prod(phi[i], phi[j]);

            auto ueval    =   sf(qp.point());
            for (size_t i = 0; i < cb1.size(); i++)
                u_rhs(i)  +=  qp.weight() * ueval  * phi[i];
            }

            //projection on the cell
            Eigen::LDLT<matrix_type> mass_matrix_ldlt = mass_matrix.ldlt();
            Iu_Th.at(i) = mass_matrix_ldlt.solve(u_rhs);
            i++;

            auto Iu = mass_matrix_ldlt.solve(u_rhs);
            std::cout<< " Iu = "<< Iu <<std::endl;

        }

        /*****************************************************************/
        post_processing_base<T,DIM,Storage>  pp_base;

        size_t num_sub_elems, num_sub_nodes;
        size_t total_sub_nodes, total_sub_elems;
        auto   test_points  = pp_base.make_vtk_points(msh,degree);
        std::vector<T> pot_Iu(test_points.size()), pot_u(test_points.size());

        total_sub_nodes = pp_base.m_total_sub_nodes;
        total_sub_elems = pp_base.m_total_sub_elems;
        num_sub_nodes   = pp_base.m_num_sub_nodes;
        auto sub_nodes       = pp_base.m_sub_nodes;
        auto vertices        = pp_base.m_vertices;


        size_t cl_cont = 0;
        size_t fc_cont = 0;
        size_t cont = 0;

        for(auto& cel : msh)
        {
            vector_type uh_TF =  Uh_Th[cl_cont];
            dynamic_vector<scalar_type> rec(cb1.size()), rec_Iu(cb1.size()) ;
            rec    = uh_TF.head(cb1.size());
            rec_Iu = Iu_Th.at(cl_cont);
            auto cl_faces =faces(msh,cel);
            for(auto& fc : cl_faces)
            {
                for (size_t itp = 0; itp < num_sub_nodes; itp++)
                {
                    auto pot  = 0.;
                    auto idx  = itp + num_sub_nodes*fc_cont;
                    auto tp   = test_points[idx];
                    auto phi  = cb1.eval_functions(msh, cel, tp);

                    for (size_t i = 0; i < cb1.size(); i++)
                    {
                        pot_u.at(cont)  +=  phi[i] * rec(i);
                        pot_Iu.at(cont) +=  phi[i] * rec_Iu(i);
                    }

                    ++cont;
                }
                ++fc_cont;
            }
            ++cl_cont;
        }

        std::vector<T> y_val = pot_u;

        std::vector<T> x_val( test_points.size());

        for(size_t i = 0; i <  test_points.size(); i++ )
        {
            auto p = test_points.at(i) ;
            x_val.at(i) = p.x();
        }
        /*****************************************************************/

        /* Plot it */
        //if (rp.draw)

        Gnuplot gp;
        gp << "set grid" << std::endl;
    //    gp << "plot '-' with lines title 'projection'" << std::endl;
        //gp.send1d(std::make_pair(x_val, y_val));

        /*****************************************************************/
}
#endif

} // namespace disk
