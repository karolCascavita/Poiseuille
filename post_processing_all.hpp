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
        //if(m_degree == 0)
        //    m_degree = 1; /* Same treatment as linear functions*/

        m_num_sub_nodes   = m_degree+1;
        m_num_sub_elems   = 1;
        m_total_sub_nodes = 2 * m_num_sub_nodes * msh.cells_size();
        m_total_sub_elems = 2 * msh.cells_size();
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
                    int  idx   = cont * m_num_sub_nodes + i;

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
    post_processing_base(const size_t degree): m_degree(degree)
    { }

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
        //if(m_degree == 0)
        //    m_degree = 1; /* Same treatment as linear functions*/

        m_num_sub_nodes   = 3 * std::pow(4,m_degree);
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

    typedef quadrature<mesh_type, cell_type>     cell_quadrature_type;
    typedef quadrature<mesh_type, face_type>     face_quadrature_type;

    typedef scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    point_vector_type sub_nodes;
    point_vector_type vertices;
    std::vector<matrix_type> m_gradrec_opers;
    cell_basis_type   cb;
    face_basis_type   fb;

    std::vector<point_type> test_points;
    size_t num_sub_elems;
    size_t num_sub_nodes;
    size_t total_sub_nodes;
    size_t total_sub_elems;
    size_t m_degree;


//    post_processing(const mesh<T,DIM,Storage>& msh): pp_base(msh)
    post_processing(const std::vector<matrix_type>& gradrec_opers,
                    const size_t degree):
     pp_base(degree)
    {
        m_gradrec_opers = gradrec_opers;
        cb = cell_basis_type(degree);
        fb = face_basis_type(degree);

        m_degree = degree;
    }

    post_processing(const size_t degree):
     pp_base(degree)
    {
        cb = cell_basis_type(degree);
        m_degree = degree;
    }

    struct sizes
    {
        sizes(){};
        size_t num_nodes;
        size_t num_elems;
        size_t num_scalars;
        size_t num_vectors;
    };

    void paraview_faces(const mesh_type   & msh,
                        const std::string & name,
                        const std::vector< std::pair<dynamic_vector<T>,
                                            dynamic_vector<T>>> & vec,
                        const std::string & type)
    {
        test_points  = pp_base::make_vtk_points(msh, m_degree);
        /*VTK file*/
        std::string ext = ".vtk";
        std::ofstream  ofs(name + ext);
        ofs << "# vtk DataFile Version 3.0"<< std::endl;
        ofs << "#This file was generated by the DISK++ library"<< std::endl;
        ofs << "ASCII"<< std::endl;

        /* Data header */
        ofs << "DATASET UNSTRUCTURED_GRID\n" << std::endl;


        auto N   = m_degree;
        if(m_degree == 0)
            N = 1;
        auto total_sub_nodes = msh.faces_size() * (N + 1);
        auto total_sub_elems = msh.faces_size();

        /* Points */
        ofs << "POINTS " << total_sub_nodes<<" double"<<std::endl;

        for(auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            auto fc = *itor;
            auto pts = points(msh, fc);
            auto h = (pts[1] - pts[0]) / N;

            for(auto i = 0; i < N + 1; i++)
            {
                auto tp = pts[0] + i * h;

                for (size_t d = 0; d < DIM; d++)
                    ofs << tp.at(d) << " ";
                for (size_t d = 0; d < 3 - DIM; d++) /* VTK only supports until 3D */
                    ofs<< 0. << " ";
                ofs<< std::endl;
            }
        }
        ofs<<std::endl;

        /* Cells */
        size_t el = 0;
        ofs << "CELLS " << total_sub_elems <<' '<<  total_sub_elems *(N+1 + 1)<< std::endl;
        for (size_t i = 0; i < total_sub_elems; i++)
        {
            ofs << N + 1 << " ";
            for(size_t j=0; j < N+1; j++, el++)
                ofs << el<< " ";
            ofs<< std::endl;
        }

        /* Types of cells*/
        ofs << "CELL_TYPES " << total_sub_elems << std::endl;
            for(size_t i = 0; i < total_sub_elems; i++)
                ofs << ' ' << 4;
        ofs << std::endl;

        /* Data */
        ofs << "POINT_DATA " << total_sub_nodes << std::endl;
        ofs << "VECTORS  Vector  double"<< std::endl;

        size_t fcount = 0;
        for(auto& data: vec)
        {
            auto fc = *std::next(msh.faces_begin(), fcount++);
            auto pts = points(msh, fc);
            auto h = (pts[1] - pts[0]) / N;

            for(auto i = 0; i < N + 1; i++)
            {
                auto tp = pts[0] + i * h;

                auto phi = fb.eval_functions(msh, fc, tp);
                ofs << phi.dot(data.first)<< " "<<phi.dot(data.second) <<" "<< 0.<<std::endl;
            }
        }
        ofs.close();
        return;
    }


    void paraview(const mesh_type   & msh,
                  const std::string & name,
                  const std::vector<dynamic_vector<T>> & vector_Th,
                  const std::string & type)
    {
        test_points  = pp_base::make_vtk_points(msh, m_degree);

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
        std::vector<size_t> vtk_nodes_inside_cell = {m_degree+1,3,4};/* WK: all cells are the same type */

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
            auto cl_faces = faces(msh, cel);

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

        if(type == "scalar")
            vtk_scalar(ofs, msh, name, vector_Th);
        if(type == "vector")
            vtk_vector(ofs, msh, name, vector_Th);

        ofs.close();

        return;
    }

    void vtk_vector(std::ofstream &ofs,
                    const mesh_type& msh,
                    const std::string& name,
                    const std::vector<dynamic_vector<T>> & vec)
    {
        ofs << "VECTORS  Vector  double"<< std::endl;

        size_t cl_cont = 0;
        size_t fc_cont = 0;

        auto mindeg = 0;
        auto maxdeg = m_degree;
        auto zero_range = cb.range( mindeg,  maxdeg);

        for(auto& cl : msh)
        {
            auto cl_id = msh.lookup( cl);
            auto fcs   = faces(msh, cl);
            vector_type v = vec.at(cl_id);

            for(auto& fc : fcs)
            {
                for (size_t itp = 0; itp < num_sub_nodes; itp++)
                {
                    auto pot  = 0.;
                    auto idx  = itp + num_sub_nodes * fc_cont;
                    auto tp   = test_points.at(idx);
                    auto c_phi = cb.eval_functions(msh, cl, tp, mindeg, maxdeg);
                    matrix_type vec_phi = make_vectorial_matrix(c_phi, DIM);
                    vector_type dpot    = vec_phi * v;
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
        return;
    }

    void
    vtk_scalar( std::ofstream &ofs,
                const mesh_type& msh,
                const std::string& name,
                const std::vector<dynamic_vector<T>> & vec)
    {
        size_t cl_cont = 0;
        size_t fc_cont = 0;

        ofs << "SCALARS  Scalar  double 1"<< std::endl;
        ofs << "LOOKUP_TABLE default"<< std::endl;
        auto mindeg = 0;
        auto maxdeg = m_degree;
        auto zero_range = cb.range( mindeg,  maxdeg);

        for(auto& cel : msh)
        {
            vector_type dofs =  vec.at(cl_cont);

            auto fcs = faces(msh, cel);

            for(auto& fc : fcs)
            {
                for (size_t itp = 0; itp < num_sub_nodes; itp++)
                {
                    auto pot  = 0.;
                    auto idx  = itp + num_sub_nodes * fc_cont;
                    auto tp   = test_points.at(idx);
                    auto phi  = cb.eval_functions(msh, cel, tp, mindeg, maxdeg);

                    assert(cb.size() == dofs.size());
                    for (size_t i = 0; i < cb.size(); i++)
                        pot  +=  phi[i] * dofs(i);
                    ofs << pot <<' ';
                }
                ++fc_cont;
            }
            ++cl_cont;
        }
        return;
    }

    // Return polar angle of p with respect to origin o
    auto
    to_angle(const point_type& p, const point_type& o)
    {
      return atan2(p.y() - o.y(), p.x() - o.x());
    }

    auto
    sort_by_polar_angle(const mesh_type & msh,
                             const cell_type& cl,
                             const std::vector<point_type>& pts)
    {
        auto h = barycenter(msh, cl);
        auto sorted_pts = pts;

        std::sort( sorted_pts.begin(), sorted_pts.end(),
         [&](const point_type& va, const point_type& vb )
            {

                auto theta_a = to_angle(va, h);
                auto theta_b = to_angle(vb, h);

                return (theta_a < theta_b);
            }
        );
        //std::cout << "/ * sorted points */" << std::endl;
        //for(auto p : sorted_pts)
        //    std::cout<< "    "<<p.x()  << " "<< p.y()<< std::endl;
        return sorted_pts;
    }

    double dot(const point_type& p1, const point_type& p2)
    {
        return p1.x()*p2.x() + p1.y()*p2.y();
    }

    double distance(const point_type& p1, const point_type& p2)
    {
        return sqrt(dot(p1, p2));
    }
    bool is_inside(const mesh_type& msh, const cell_type & cl, const point_type& pt)
    {

        auto pts = points(msh, cl);
        auto vts = msh.get_vertices(cl, pts);

        if(vts.size() > 3)
            throw std::invalid_argument("Extend is_inside for other polygons.");

        auto v0 = vts[1] - vts[0];
        auto v1 = vts[2] - vts[0];
        auto v2 = pt - vts[0];

        auto dot00 = dot(v0, v0);
        auto dot01 = dot(v0, v1);
        auto dot02 = dot(v0, v2);
        auto dot11 = dot(v1, v1);
        auto dot12 = dot(v1, v2);

        auto invDenom = 1. / (dot00 * dot11 - dot01 * dot01);
        auto u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        auto v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        /* The threshold we're discussing this afternoon may be needed here */
        //return (u >= 0) && (v >= 0) && (u + v <= 1);
        /* The threshold we're discussing this afternoon may be needed here */
        bool upos  = (std::abs(u) < 1.e-5 * invDenom) || ( u > T(0));
        bool vpos  = (std::abs(v) < 1.e-5 * invDenom) || ( v > T(0));
        bool uvpos = (std::abs(u + v - 1.) < 1.e-10) || (u + v < 1.);

        //std::cout << pt.x()<<" " << pt.y() << std::endl;
        //std::cout << v0.x()<<" " << v0.y() << std::endl;
        //std::cout << v1.x()<<" " << v1.y() << std::endl;
        //std::cout << v2.x()<<" " << v2.y() << std::endl;

        return (upos && vpos && uvpos);
    }

    // Copyright 2000 softSurfer, 2012 Dan Sunday
    // This code may be freely used and modified for any purpose
    // providing that this copyright notice is included with it.
    // SoftSurfer makes no warranty for this code, and cannot be held
    // liable for any real or imagined damage resulting from its use.
    // Users of this code must verify correctness for their application.


    // a Point is defined by its coordinates {int x, y;}
    //===================================================================


    // isLeft(): tests if a point is Left|On|Right of an infinite line.
    //    Input:  three points P0, P1, and P2
    //    Return: >0 for P2 left of the line through P0 and P1
    //            =0 for P2  on the line
    //            <0 for P2  right of the line
    //    See: Algorithm 1 "Area of Triangles and Polygons"
    auto
    isLeft( const point_type& P0, const point_type& P1, const point_type& P2 )
    {

        auto ret = ( (P1.x() - P0.x()) * (P2.y() - P0.y())
                - (P2.x() -  P0.x()) * (P1.y() - P0.y()) );
        #if 0
        std::cout << " is_Left_inside number : "<< ret << std::endl;

        std::cout << "P0: "<< P0.x() << P0.y()<< std::endl;
        std::cout << "P1: "<< P1.x() << P1.y()<< std::endl;
        std::cout << "P2: "<< P2.x() << P2.y()<< std::endl;
        #endif
        return ret;
    }
    //===================================================================


    // cn_PnPoly(): crossing number test for a point in a polygon
    //      Input:   P = a point,
    //               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
    //      Return:  0 = outside, 1 = inside
    // This code is patterned after [Franklin, 2000]
    auto
    cn_PnPoly( const point_type  & P, const std::vector<point_type>& vts, const int n )
    {
        auto V = vts;
        V.push_back(vts.at(0));

        int    cn = 0;    // the  crossing number counter

        // loop through all edges of the polygon
        for (int i = 0; i < n; i++)
        {
            // edge from V[i]  to V[i+1]
           if ( ((V[i].y() <= P.y()) && (V.at(i+1).y() > P.y() ))     // an upward crossing
                    || ((V.at(i).y() > P.y()) && (V.at(i+1).y() <=  P.y() )))
            {
                // a downward crossing
                // compute  the actual edge-ray intersect x-coordinate
                T vt = (P.y()  - V.at(i).y()) / (V.at(i+1).y() - V.at(i).y());
                if (P.x() <  V.at(i).x() + vt * (V.at(i+1).x() - V.at(i).x() )) // P.x < intersect
                     ++cn;   // a valid crossing of y=P.y right of P.x
            }
        }
        return (cn&1);    // 0 if even (out), and 1 if  odd (in)
    }
    //===================================================================


    // wn_PnPoly(): winding number test for a point in a polygon
    //      Input:   P = a point,
    //               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
    //      Return:  wn = the winding number (=0 only when P is outside)
    bool
    wn_PnPoly(const mesh_type& msh,
              const cell_type& cl,
              const point_type& P,
              const std::vector<point_type>& vts )
    {
        auto spts = sort_by_polar_angle(msh, cl, vts);
        spts.push_back(spts.at(0));

        int    wn = 0;    // the  winding number counter
        //#if 0

        // loop through all edges of the polygon
        for (int i = 0 ; i < spts.size() - 1; i++)
        {   // edge from V[i] to  V[i+1]
            if (spts.at(i).y() <= P.y() )
            {          // start y <= P.y
                if (spts.at(i+1).y()  > P.y())      // an upward crossing
                {
                    //std::cout << "***** upward crossing edge : "<< i <<"********** "<< std::endl;
                    if (isLeft( spts.at(i), spts.at(i+1), P) > T(0))  // P left of  edge
                    {
                        //std::cout << "** valid up intersect--------------------------" << std::endl;
                        ++wn;
                    }               // have  a valid up intersect
                }
            }
            else
            {
                // start y > P.y (no test needed)
                if ( spts.at(i+1).y()  <= P.y())     // a downward crossing
                {
                    //std::cout << "**** downward crossing  edge : "<< i <<"********** "<< std::endl;

                    if (isLeft( spts.at(i), spts.at(i+1), P) < T(0))  // P right of  edge
                    {
                        //std::cout << "** valid down intersect-----------------------" << std::endl;
                        --wn;
                    }               // have  a valid down intersect
                }
            }
            //std::cout << "wn = "<< wn << std::endl;
        }
        //#endif
        if(wn != 0)
            return true;
        return false;
    }

    void
    plot_over_line(const mesh_type    & msh,
                   const std::vector<vector_type>   & rec_vel_Th,
                   const std::string & filename)
    {
        std::ofstream pfs(filename);
        if(!pfs.is_open())
            std::cout << "Error opening file :"<< filename <<std::endl;

        point_type p1 = point_type({-1.0, 0.0});
        point_type p2 = point_type({ 1.0, 0.0});

        auto N = 400;
        auto h = (p2 - p1)/N;
        size_t i = 0 ;
        auto pts = std::vector<point_type>(N);

        auto mindeg = 0;
        auto maxdeg = m_degree + 1;

        cell_basis_type   cb_max(maxdeg);

        auto all_range = cb_max.range( mindeg,  maxdeg);

        for(auto& p: pts)
        {
            p = p1 + i * h;

            for(auto cl: msh)
            {
                auto cell_pts = points(msh, cl);
                auto vts = msh.get_vertices(cl, cell_pts);

                if(wn_PnPoly( msh, cl, p, vts ))
                {
                    if(!is_inside(msh, cl, p))
                    {
                        std::cout << " ** Cell_points: " << std::endl;
                        for(auto cp : cell_pts)
                            std::cout << "  P: "<< cp.x() << cp.y()<< std::endl;
                        std::cout << " ** Cell_vertices: " << std::endl;
                        for(auto cp : vts)
                            std::cout << "  V: "<< cp.x() << cp.y()<< std::endl;

                        std::cout << " ** Point : " << std::endl;
                        std::cout << "  P: "<< p.x() << p.y()<< std::endl;

                    }

                    auto cl_id  = cl.get_id();
                    vector_type u_T   = rec_vel_Th.at(cl_id);
                    auto phi  = cb_max.eval_functions(msh, cl, p, mindeg, maxdeg);
                    pfs<< p.x() << " "<< p.y() << " "<< phi.dot(u_T)<<std::endl;
                }
            }
            i++;
        }
        pfs.close();
        return;
    }


    void
    plot_over_line(const mesh_type    & msh,
                   //const point_type   & p1,
                   //const point_type   & p2,
                   const std::vector<vector_type>   & vel_cell_Th,
                   const std::vector<vector_type>   & grec_vel_Th,
                   const std::string & filename)
    {
        std::ofstream pfs(filename);
        if(!pfs.is_open())
            std::cout << "Error opening file :"<< filename <<std::endl;

        point_type p1 = point_type({-1.0, 0.0});
        point_type p2 = point_type({ 1.0, 0.0});

        auto N = 400;
        auto h = (p2 - p1)/N;
        size_t i = 0 ;
        auto pts = std::vector<point_type>(N);

        auto mindeg = 0;
        auto maxdeg = m_degree;
        auto zero_range = cb.range( mindeg,  maxdeg);

        for(auto& p: pts)
        {
            p = p1 + i * h;

            for(auto cl: msh)
            {
                auto cell_pts = points(msh, cl);
                auto vts = msh.get_vertices(cl, cell_pts);

                if(wn_PnPoly( msh, cl, p, vts ))
                {
                    if(!is_inside(msh, cl, p))
                    {
                        std::cout << " ** Cell_points: " << std::endl;
                        for(auto cp : cell_pts)
                            std::cout << "  P: "<< cp.x() << cp.y()<< std::endl;
                        std::cout << " ** Cell_vertices: " << std::endl;
                        for(auto cp : vts)
                            std::cout << "  V: "<< cp.x() << cp.y()<< std::endl;

                        std::cout << " ** Point : " << std::endl;
                        std::cout << "  P: "<< p.x() << p.y()<< std::endl;

                    }

                    auto cl_id  = cl.get_id();
                    vector_type u_T   = vel_cell_Th.at(cl_id);
                    vector_type Gu_T  = grec_vel_Th.at(cl_id);

                    auto phi  = cb.eval_functions(msh, cl, p, mindeg, maxdeg);
                    matrix_type vec_phi  =  make_vectorial_matrix(phi, mesh_type::dimension);
                    vector_type grad_fun =  vec_phi * Gu_T;

                    pfs<< p.x() << " "<< p.y() << " "<< phi.dot(u_T)<< "  ";
                    pfs<< grad_fun(0)<<" "<<grad_fun(1) <<" " << cl_id<<std::endl;
                }
            }
            i++;
        }
        pfs.close();
        return;
    }

};


template<typename T>
struct solution_errors
{
    solution_errors():u_uh(0.), Iu_uh(0.), Du_Guh(0.)
    {}
    T  u_uh;
    T Iu_uh;
    T Du_Guh;
};
template<typename T, size_t DIM, typename Storage, typename MeshParameters>
void
postdiff(const  mesh<T,DIM,Storage>&  msh,
            const std::vector<dynamic_vector<T>> & vel_cell_Th,
            const MeshParameters  &  mp,
            const size_t     degree,
            const size_t     imsh)
{
    typedef mesh<T,DIM,Storage> mesh_type;
    typedef typename mesh_type::cell cell_type;
    typedef dynamic_vector<T> vector_type;
    typedef dynamic_matrix<T> matrix_type;

    post_processing<T, DIM, Storage> pp(degree);

    std::string  msh_str = tostr(imsh);
    std::string  info    =  mp.summary + "_R" + msh_str;
    std::string  solution_file  = mp.directory + "/solution" + info;

    pp.paraview(msh, solution_file, vel_cell_Th, "scalar");
}

template<typename T, size_t DIM, typename Storage, typename MeshParameters,
        typename PlasticData>
void
postprocess(const  mesh<T,DIM,Storage>&  msh,
            const tensors<mesh<T,DIM,Storage>>   & stresses_Th,
            const std::vector<dynamic_vector<T>> & velocity_Th,
            const std::vector<dynamic_matrix<T>> & gradrec_opers,
            const std::vector<dynamic_matrix<T>> & rec_opers,
            const MeshParameters  &  mp,
            const PlasticData   & pst,
            const size_t     degree,
            const size_t     imsh)

{
    typedef mesh<T,DIM,Storage> mesh_type;
    typedef typename mesh_type::cell cell_type;
    typedef dynamic_vector<T> vector_type;
    typedef dynamic_matrix<T> matrix_type;

    post_processing<T, DIM, Storage> pp(gradrec_opers, degree);

    std::vector<vector_type> grad_rec_Th(msh.cells_size());
    std::vector<vector_type> vel_cell_Th(msh.cells_size());
    std::vector<vector_type> rec_vel_Th(msh.cells_size());

    std::vector<std::pair<vector_type, vector_type>> stresses_Fh(msh.faces_size());

    typedef scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;


    auto cell_basis  = cell_basis_type(degree + 1);
    auto cbs = cell_basis.range(0,degree).size();
    auto all_range_size = cell_basis.size();
    auto one_range_size = cell_basis.range(1, degree + 1).size();

    size_t fcount = 0;

    for(auto& cl : msh)
    {
        auto cl_id = msh.lookup( cl);

        matrix_type grec_op =  gradrec_opers.at(cl_id);
        matrix_type rec_op  =  rec_opers.at(cl_id);
        vector_type uh_TF   =  velocity_Th.at(cl_id);

        grad_rec_Th.at(cl_id) = grec_op * uh_TF;
        vel_cell_Th.at(cl_id) = uh_TF.head(cbs);

        vector_type rec = vector_type::Zero(all_range_size);
        rec.tail(one_range_size) = rec_op * uh_TF;
        rec(0) = uh_TF(0);
        rec_vel_Th.at(cl_id) = rec;
    }

    for(auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
    {
        auto fc = *itor;
        auto pts = points(msh, fc);

        auto cells_ids = face_owner_cells_ids(msh, fc);

        auto cl_1 = *std::next( msh.cells_begin(), cells_ids.first);
        auto cl_2 = cl_1;

        if(!msh.is_boundary(fc))
            cl_2 = *std::next( msh.cells_begin(), cells_ids.second);

        vector_type v1 = stresses_Th.at_face(msh, cl_1, fc);
        vector_type v2 = stresses_Th.at_face(msh, cl_2, fc);

        stresses_Fh.at(fcount++) = std::make_pair(v1, v2);
    }


    std::string  msh_str = tostr(imsh);
    std::string  info =  mp.summary + "_R" + msh_str;

    std::string solution_file  = mp.directory + "/solution" + info;
    std::string gradient_file  = mp.directory + "/gradient" + info;
    std::string stresses_file  = mp.directory + "/stresses" + info;
    std::string stresses_Fh_file  = mp.directory + "/stresses_Fh" + info;


    pp.paraview(msh, solution_file, vel_cell_Th, "scalar");
    pp.paraview(msh, gradient_file, grad_rec_Th, "vector");
    pp.paraview(msh, stresses_file, stresses_Th.at_all_cells(), "vector");
    pp.paraview_faces(msh, stresses_Fh_file, stresses_Fh, "vector");


    std::string plotline_file  = mp.directory + "/plot_over_line" + info;
    std::string rec_plotline_file  = mp.directory + "/plot_over_line_rec" + info;

    pp.plot_over_line(msh, vel_cell_Th, grad_rec_Th, plotline_file+ ".dat");
    pp.plot_over_line(msh, rec_vel_Th, rec_plotline_file +".dat");

    save_data(velocity_Th, mp.directory + "/velocities" + info + ".txt");
    save_data(stresses_Th, mp.directory + "/multiplier" + info + ".txt");


    auto exfunc_name = mp.directory + "/exFunction.m";

    std::ofstream effs(exfunc_name);
    if (!effs.is_open())
        std::cout << "Error opening file :"<<exfunc_name <<std::endl;

    effs << "col = 2"<<std::endl;
    effs << "Rconvergence = false"<<std::endl;
    effs << "if Rconvergence ==true"<<std::endl;
    effs << "   col = 5"<<std::endl;
    effs << "end"<<std::endl;

    effs << "\% Mesh "<<std::endl;
    effs << "draw_mesh  = false"<<std::endl;
    effs << "print_mesh = false"<<std::endl;

    effs << "\% Error "<<std::endl;
    effs << "draw_error = true"<<std::endl;
    effs << "print_error= false"<<std::endl;

    effs << "\% Solution "<<std::endl;
    effs << "draw_solution   = true"<<std::endl;
    effs << "print_solution  = false"<<std::endl;

    effs << "\% Estimator "<<std::endl;
    effs << "draw_estimator  = false"<<std::endl;
    effs << "print_estimator = false"<<std::endl;
    effs << "error_name = 'nothing'"<<std::endl;
    effs << "execute" + mp.summary<<std::endl;
    effs << "if print_error == true "<<std::endl;
    effs << "   legend('-DynamicLegend')" <<std::endl;
    effs << "   set(gca,'box','on'); set(gcf,'color','w');"<<std::endl;
    effs << "   set(gca,'fontsize',12); set(gca,'xgrid', 'on', 'ygrid', 'on');"<<std::endl;
    effs << "   if(col == 2)"<<std::endl;
    effs << "     title(' Total residual with k = "<< degree<<"; Bi = "<<pst.Bn<<"');"<<std::endl;
    effs << "     print('-depsc2', '-loose', 'Total_residual_"<<mp.summary<<"_R"<<imsh<<".eps')"<<std::endl;
    effs << "   end"<<std::endl;
    effs << "   if(col == 3)"<<std::endl;
    effs << "     title(' Residual with k = "<< degree<<"; Bi = "<<pst.Bn<<"');"<<std::endl;
    effs << "     print('-depsc2', '-loose', 'Residual_"<<mp.summary<<"_R"<<imsh<<".eps')"<<std::endl;
    effs << "   end"<<std::endl;
    effs << "   if(col == 4)"<<std::endl;
    effs << "     title(' |G_h^k(u_h^{n+1}- u_h^n)| with k = "<< degree<<"; Bi = "<<pst.Bn<<"');"<<std::endl;
    effs << "     print('-depsc2', '-loose', 'GUerror_"<<mp.summary<<"_R"<<imsh<<".eps')"<<std::endl;
    effs << "   end"<<std::endl;
    effs << "   if(col == 5)"<<std::endl;
    effs << "     title(' |S_h^k(u_h^{n+1}- u_h^n)| with k = "<< degree<<"; Bi = "<<pst.Bn<<"');"<<std::endl;
    effs << "     print('-depsc2', '-loose', 'SUerror_"<<mp.summary<<"_R"<<imsh<<".eps')"<<std::endl;
    effs << "   end"<<std::endl;
    effs << "   if(col == 6)"<<std::endl;
    effs << "     title(' |\\sigma_h^{n+1}- \\sigma_h^n)| with k = "<< degree<<"; Bi = "<<pst.Bn<<"');"<<std::endl;
    effs << "     print('-depsc2', '-loose', 'Sigma_"<<mp.summary<<"_R"<<imsh<<".eps')"<<std::endl;
    effs << "   end"<<std::endl;
    effs << "   if(col == 8)"<<std::endl;
    effs << "     title(' |\\varsigma_h^{n+1}- \\varsigma_h^n)| with k = "<< degree<<"; Bi = "<<pst.Bn<<"');"<<std::endl;
    effs << "     print('-depsc2', '-loose', 'Varigma_"<<mp.summary<<"_R"<<imsh<<".eps')"<<std::endl;
    effs << "   end"<<std::endl;
    effs << "   if(col == 10)"<<std::endl;
    effs << "     title(' |u_h^{n+1}- I_hu)| with k = "<< degree<<"; Bi = "<<pst.Bn<<"');"<<std::endl;
    effs << "     print('-depsc2', '-loose', 'Uerror_"<<mp.summary<<"_R"<<imsh<<".eps')"<<std::endl;
    effs << "   end"<<std::endl;
    effs << "end"<<std::endl;
    effs.close();

}
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
        exfs<< "    figure;estimator_tot"<<mp.summary<<"_RC"<<imsh<<";";
        exfs<< "    if print_estimator == true";
        exfs<< "        set(gca,'box','on'); set(gcf,'color','w');set(gca,'fontsize',12);"<<std::endl;
        exfs<< "        circle(0,0,"<<pst.Bn<<",'w')"<<std::endl;
        exfs<< "        print('-depsc2', '-loose', 'estimator_tot"<<mp.summary<<"_RC"<<imsh<<".eps')"<<std::endl;
        exfs<< "    end;"<<std::endl;
        exfs<< "    figure;estimator_res"<<mp.summary<<"_RC"<<imsh<<";";
        exfs<< "    if print_estimator == true";
        exfs<< "        set(gca,'box','on'); set(gcf,'color','w');set(gca,'fontsize',12);"<<std::endl;
        exfs<< "        circle(0,0,"<<pst.Bn<<",'w')"<<std::endl;
        exfs<< "        print('-depsc2', '-loose', 'estimator_res"<<mp.summary<<"_RC"<<imsh<<".eps')"<<std::endl;
        exfs<< "    end;"<<std::endl;
        exfs<< "    figure;estimator_str"<<mp.summary<<"_RC"<<imsh<<";";
        exfs<< "    if print_estimator == true";
        exfs<< "        set(gca,'box','on'); set(gcf,'color','w');set(gca,'fontsize',12);"<<std::endl;
        exfs<< "        circle(0,0,"<<pst.Bn<<",'w')"<<std::endl;
        exfs<< "        print('-depsc2', '-loose', 'estimator_str"<<mp.summary<<"_RC"<<imsh<<".eps')"<<std::endl;
        exfs<< "    end;"<<std::endl;
        exfs<< "end"<<std::endl;
    }

    exfs<< "if draw_error == true "<<std::endl;
    exfs<< "    e"<< imsh <<"= load('error"<<mp.summary<<"_R"<<imsh<<".dat');"<<std::endl;
    exfs<< "    figure(100); hold on; "<< std::endl;
    exfs<< "    loglog(e"<<imsh<<"(:,1),e"<<imsh<<"(:,col), 'DisplayName', 'step "<<imsh<<"' , 'color', rand(1,3));"<<std::endl;
    exfs<< "    set(gca,'box','on'); set(gcf,'color','w');"<<std::endl;
    exfs<< "    set(gca,'fontsize',12); set(gca,'xgrid', 'on', 'ygrid', 'on');"<<std::endl;
    exfs<< "    set(gca, 'xscale', 'log'); set(gca, 'yscale', 'log')"<<std::endl;
    exfs<< "end" << std::endl;
    exfs<< std::endl;

    exfs<< "if draw_mesh == true "<<std::endl;
    exfs<< "    figure; mesh"<<mp.summary<<"_R" <<imsh<<";";
    exfs<< "    if print_mesh == true";
    exfs<< "        set(gca,'box','on'); set(gcf,'color','w');set(gca,'fontsize',12);"<<std::endl;
    exfs<< "        circle(0,0,"<<pst.Bn<<",'r')"<<std::endl;
    exfs<< "        print('-depsc2', '-loose', 'mesh_B"<<tostr(10*pst.Bn)<<mp.summary<<"_R" <<imsh<<".eps')"<<std::endl;
    exfs <<"    end"<<std::endl;
    exfs <<"end"<<std::endl;
    exfs<< std::endl;

    exfs<< "if draw_solution == true "<<std::endl;
    exfs<< "    m = 101; n = 102; "<<std::endl;
    exfs<< "    e"<< imsh <<"= load('plot_over_line"<<mp.summary<<"_R"<<imsh<<".dat');"<<std::endl;
    exfs<< "    figure(m); "<<std::endl;
    exfs<< "    hold on; plot(e"<<imsh<<"(:,1),e"<<imsh<<"(:,3), 'DisplayName', 'step "<<imsh<<"' , 'color', rand(1,3));"<<std::endl;
    exfs<< "    set(gca,'box','on'); set(gcf,'color','w');"<<std::endl;
    exfs<< "    set(gca,'fontsize',12); set(gca,'xgrid', 'on', 'ygrid', 'on');"<<std::endl;
    exfs<< "    xlabel('x');  ylabel('u')" << std::endl;
    exfs<< std::endl;
    exfs<< "    er"<< imsh <<"= load('plot_over_line_rec"<<mp.summary<<"_R"<<imsh<<".dat');"<<std::endl;
    exfs<< "    hold on; plot(er"<<imsh<<"(:,1),er"<<imsh<<"(:,3), 'DisplayName', 'step "<<imsh<<"' , 'color', rand(1,3));"<<std::endl;
    exfs<< "    set(gca,'box','on'); set(gcf,'color','w');"<<std::endl;
    exfs<< "    set(gca,'fontsize',12); set(gca,'xgrid', 'on', 'ygrid', 'on');"<<std::endl;
    exfs<< "    xlabel('x');  ylabel('u')" << std::endl;
    exfs<< "    if print_solution == true";
    exfs<< "        set(gca,'box','on'); set(gcf,'color','w');set(gca,'fontsize',12);"<<std::endl;
    exfs<< "        print('-depsc2', '-loose', 'plotoline_sol_"<<mp.summary<<"_R"<<imsh<<".eps')"<<std::endl;
    exfs<< "    end;"<<std::endl;
    exfs<< std::endl;
    exfs<< "    figure(n);"<<std::endl;
    exfs<< "     hold on; plot(e"<<imsh<<"(:,1),sqrt(e"<<imsh<<"(:,5).^2 + e"<<imsh<<"(:,4).^2), 'DisplayName', 'step "<<imsh<<"' , 'color', rand(1,3));"<<std::endl;
    exfs<< "    set(gca,'box','on'); set(gcf,'color','w');"<<std::endl;
    exfs<< "    set(gca,'fontsize',12); set(gca,'xgrid', 'on', 'ygrid', 'on');"<<std::endl;
    exfs<< "    xlabel('x');  ylabel('|G_h u_h|')" << std::endl;
    exfs<< "    if print_solution == true";
    exfs<< "        set(gca,'box','on'); set(gcf,'color','w');set(gca,'fontsize',12);"<<std::endl;
    exfs<< "        print('-depsc2', '-loose', 'plotoline_grad_"<<mp.summary<<"_R"<<imsh<<".eps')"<<std::endl;
    exfs<< "    end;"<<std::endl;
    exfs<< "    exact("<<pst.Bn <<","<< 30<<", m, n )" << std::endl;
    exfs<< "end" << std::endl;
}

} // namespace disk
