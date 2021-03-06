/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
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

 /*
  * Copyright (C) 2013-2016, Matteo Cicuttin - matteo.cicuttin@uniud.it
  * Department of Electrical Engineering, University of Udine
  * All rights reserved.
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions are met:
  *     * Redistributions of source code must retain the above copyright
  *       notice, this list of conditions and the following disclaimer.
  *     * Redistributions in binary form must reproduce the above copyright
  *       notice, this list of conditions and the following disclaimer in the
  *       documentation and/or other materials provided with the distribution.
  *     * Neither the name of the University of Udine nor the
  *       names of its contributors may be used to endorse or promote products
  *       derived from this software without specific prior written permission.
  *
  * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
  * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
  * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  */

#pragma once

#include <vector>
#include <array>
#include <fstream>
#include <cassert>
#include <thread>

#include "geometry/geometry_generic.hpp"
#include "geometry/geometry_simplicial.hpp"

#include "mapped_file.h"
#include "strtot.hpp"

namespace disk {

template<typename mesh_type>
class mesh_loader
{
public:
    mesh_loader() {}

    virtual bool    read_mesh(const std::string&) { return false; }
    virtual bool    populate_mesh(mesh_type&)    = 0;

    virtual ~mesh_loader() {}
};

template<typename T, size_t N>
class uniform_mesh_loader
{
    static_assert(N == 1, "at the moment only 1D uniform meshes are available");
};

template<typename T>
class uniform_mesh_loader<T,1> : public mesh_loader<generic_mesh<T, 1>>
{
    typedef generic_mesh<T,1>                       mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;

    T       m_h, m_x_min, m_x_max;
    size_t  m_number_of_elements;

public:
    uniform_mesh_loader()
        : m_x_min(T(0)), m_x_max(T(1)), m_number_of_elements(8)
    {
        m_h = fabs(m_x_max - m_x_min) / m_number_of_elements;
    }

    uniform_mesh_loader(T x_min, T x_max, size_t N)
        : m_x_min(x_min), m_x_max(x_max), m_number_of_elements(N)
    {
        m_h = fabs(m_x_max - m_x_min) / m_number_of_elements;
    }

    bool populate_mesh(mesh_type& msh)
    {
        std::cout << " *** POPULATING UNIFORM 1D MESH ***" << std::endl;
        auto storage = msh.backend_storage();

        auto num_edges = m_number_of_elements;
        auto num_nodes = m_number_of_elements + 1;

        storage->points.resize(num_nodes);
        storage->nodes.resize(num_nodes);
        storage->edges.resize(num_edges);

        for (size_t i = 0; i < num_edges; i++)
        {
            storage->points.at(i)   = point_type({m_x_min + (m_h * i)});
            storage->points.at(i+1) = point_type({m_x_min + (m_h * (i + 1))});

            auto n0 = typename node_type::id_type(i);
            auto n1 = typename node_type::id_type(i+1);
            auto e = edge_type{{n0, n1}};

            std::vector<point_identifier<1>> pts(2);
            pts[0] = point_identifier<1>(i);
            pts[1] = point_identifier<1>(i+1);
            e.set_point_ids(pts.begin(), pts.end());
            auto id = typename edge_type::id_type(i);
            e.set_element_id(id);

            storage->edges.at(i) = e;
        }

        for (size_t i = 0; i < num_nodes; i++)
            storage->nodes.at(i) = node_type(point_identifier<1>(i));

        storage->boundary_nodes.resize(num_nodes);
        storage->boundary_nodes.at(0) = true;
        storage->boundary_nodes.at(num_nodes - 1) = true;

        return true;
    }
};

template<typename T, size_t N>
class fvca5_mesh_loader
{
    static_assert(N == 2, "fvca5 is a 2D-only mesh format");
};

template<typename T>
class fvca5_mesh_loader<T,2> : public mesh_loader<generic_mesh<T, 2>>
{
    //const char* keywords[] = {"triangles", "quadrangles", "pentagons", "hexagons",
    //                     "heptagons", "octagons", "enneagons" };

    typedef generic_mesh<T,2>                       mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;

    template<size_t N> using polygon = std::array<ident_impl_t, N>;
public:
    std::vector<point_type>                         m_points;
    std::vector<polygon<3>>                         m_triangles;
    std::vector<polygon<4>>                         m_quadrangles;
    std::vector<polygon<5>>                         m_pentagons;
    std::vector<polygon<6>>                         m_hexagons;
    std::vector<polygon<7>>                         m_heptagons;
    std::vector<polygon<8>>                         m_octagons;
    std::vector<polygon<9>>                         m_enneagons;
    std::vector<polygon<10>>                        m_decagons;
    std::vector<polygon<11>>                        m_hendecagons;
    std::vector<polygon<12>>                        m_dodecagons;
    std::vector<polygon<13>>                        m_triadecagons;
    std::vector<polygon<14>>                        m_tesseradecagons;
    std::vector<polygon<15>>                        m_pentadecagons;

    std::vector<std::array<ident_impl_t, 2>>        m_boundary_edges;
    std::vector<std::array<ident_impl_t, 4>>        m_edges;
    std::vector<int> m_index_transf;

    //std::vector<std::set<typename edge_type::id_type>>      m_edge_node_connectivity;
private:
    bool fvca5_read_points(std::ifstream& ifs)
    {
        size_t      elements_to_read;
        T           x, y;

        ifs >> elements_to_read;
        std::cout << "Attempting to read " << elements_to_read << " points" << std::endl;
        m_points.reserve(elements_to_read);

        for (size_t i = 0; i < elements_to_read; i++)
        {
            ifs >> x >> y;
            m_points.push_back(point_type{x,y});
        }

        return true;
    }

    template<typename U, size_t N>
    bool fvca5_read_tuples(std::ifstream& ifs,
                           std::vector<std::array<U, N>>& tuples)
    {
        size_t      elements_to_read;

        ifs >> elements_to_read;
        std::cout << "Attempting to read " << elements_to_read << " " << N << "-ples " << std::endl;

        tuples.reserve(elements_to_read);

        for (size_t i = 0; i < elements_to_read; i++)
        {
            std::array<U, N> tuple;
            for (size_t j = 0; j < N; j++)
            {
                U val;
                ifs >> val;
                tuple[j] = val-1;
            }

            tuples.push_back(tuple);
        }

        return true;
    }

    bool fvca5_read(const std::string& filename)
    {
        std::ifstream   ifs(filename);
        std::string     keyword;

        std::cout << " filename : "<<filename  << std::endl;
        if (!ifs.is_open())
        {
            std::cout << "Error opening " << filename << std::endl;
            return false;
        }

        ifs >> keyword;
        if ( keyword != "vertices" )
        {
            std::cout << "Expected keyword \"vertices\"" << std::endl;
            return false;
        }

        fvca5_read_points(ifs);

        ifs >> keyword;
        if ( keyword == "triangles" )
        {
            m_triangles.clear();
            fvca5_read_tuples(ifs, m_triangles);
            ifs >> keyword;
        }

        if ( keyword == "quadrangles" )
        {
            m_quadrangles.clear();
            fvca5_read_tuples(ifs, m_quadrangles);
            ifs >> keyword;
        }

        if ( keyword == "pentagons" )
        {
            m_pentagons.clear();
            fvca5_read_tuples(ifs, m_pentagons);
            ifs >> keyword;
        }

        if ( keyword == "hexagons" )
        {
            m_hexagons.clear();
            fvca5_read_tuples(ifs, m_hexagons);
            ifs >> keyword;
        }

        if ( keyword == "heptagons" )
        {
            m_heptagons.clear();
            fvca5_read_tuples(ifs, m_heptagons);
            ifs >> keyword;
        }
        if ( keyword == "octagons" )
        {
            m_heptagons.clear();
            fvca5_read_tuples(ifs, m_octagons);
            ifs >> keyword;
        }

        if ( keyword == "edges" )
        {
            std::getline(ifs, keyword); //drop the rest of the line
            m_boundary_edges.clear();
            fvca5_read_tuples(ifs, m_boundary_edges);
            ifs >> keyword;
        }
        else
        {
            std::cout << "Error parsing FVCA5 file boundary" << std::endl;
            return false;
        }

        if ( keyword == "all" )
        {
            std::getline(ifs, keyword); //drop the rest of the line
            m_edges.clear();
            fvca5_read_tuples(ifs, m_edges);
        }
        else
        {
            std::cout << "Error parsing FVCA5 file edges" << std::endl;
            return false;
        }

        ifs.close();
        return true;
    }

    template<size_t N>
    void put_polygons(mesh_type& msh,
                      const std::vector<polygon<N>>& polys,
                      std::vector<surface_type>& surfedg)
    {
        auto storage = msh.backend_storage();

        for (auto& p : polys)
        {
            //std::cout << "/* message 0 */" << std::endl;
            std::vector<typename edge_type::id_type> surface_edges(N);
            //std::cout << "/* message 1 */" << std::endl;

            for (size_t i = 0; i < N; i++)
            {
                //std::cout << "/* message 2 */" << std::endl;

                auto pt1 = typename node_type::id_type(p[i]);
                auto pt2 = typename node_type::id_type(p[(i+1) % N]);

                if (pt2 < pt1)
                    std::swap(pt1, pt2);

                edge_type edge{{pt1, pt2}};

                auto edge_id = find_element_id(storage->edges.begin(),
                                               storage->edges.end(), edge);

                if (!edge_id.first)
                    throw std::invalid_argument("Edge not found (hanging nodes?)");

                surface_edges[i] = edge_id.second;
                //std::cout << "/* message 7 */" << std::endl;

            }
            //std::cout << "/* message 8 */" << std::endl;
            auto surface = surface_type(surface_edges);
            //std::cout << "/* message 9 */" << std::endl;

            surface.set_point_ids(p.begin(), p.end()); /* XXX: crap */
            surfedg.push_back( surface );
        }

    }
    template<typename ElementType>
    void
    index_transf(const std::vector<ElementType>& elements)
    {
        std::vector<int> index(elements.size(), 0);
        for(int i = 0 ; i != index.size() ; i++)
            index.at(i) = i;

        std::sort(index.begin(), index.end(),[&](const int& a, const int& b)
            {
                surface_type s1  = elements.at(a);
                surface_type s2  = elements.at(b);
                return (s1 < s2);
            }
        );
        m_index_transf = index;
        std::cout << "INDEX_TRANSF:" << std::endl;
        for(auto& id: index)
            std::cout << id<<"  ";
        std::cout<< std::endl;
    }
    std::pair<bool, std::vector<typename point_type::id_type>>
    is_special_polygon(const mesh_type& msh, const surface_type& cl)
    {
        auto pts  = points(msh,cl);
        auto nonc_pts = pts;
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
         std::vector<typename point_type::id_type> vertices(ns.size());
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
                 vertices.at(vcount++) =  typename point_type::id_type(i);
         }
        if(vcount != ns.size())
             std::logic_error(" Incorrect procedure to find vertices");

        bool has_hang_nodes(false);
        if(vertices.size() != pts.size())
            has_hang_nodes = true;

        return std::make_pair(has_hang_nodes, vertices);
    }
    /*Set boundary number*/
    template<typename EdgeType, typename Storage>
    size_t
    set_bnd_number(const EdgeType& edge, const Storage& storage)
    {
        T x0 = 0.5;
        auto pts_ids  = edge.point_ids();
        auto p1   = storage->points.at(pts_ids.at(0));
        auto p2   = storage->points.at(pts_ids.at(1));

        if( p1.x() <= x0 && p2.x() <= x0 )
            return 1;
        else
            return 2;
    }
    #if 0
    set_neighbors(mesh_type& msh, m_edges, std::vector<surface_type>& surfaces);
    {
        auto storage = msh.backend_storage();
        auto edges   = storage->edges;

        for(auto& surf:surfaces)
        {
            auto edges = surf.subelement_size();

            for(size_t i = 0 ; i < nfaces; i++)
            {
                auto edge_id = surf.m_sids_ptrs(i);
                (storage->edge.at(edge_id).m_nghs_ptrs).push_back(srf.id);
            }
        }
    }
    #endif

public:
    fvca5_mesh_loader() = default;

    bool read_mesh(const std::string& filename)
    {
        std::cout << " *** READING FVCA5 MESH ***" << std::endl;
        return fvca5_read(filename);
    }

    bool populate_mesh(mesh_type& msh)
    {



        std::cout << " *** POPULATING FVCA5 MESH ***" << std::endl;
        auto storage = msh.backend_storage();

        //std::cout << "/* message Points */" << std::endl;
        /* Points */
        size_t nodes_size = m_points.size();
        storage->points = std::move(m_points);

        //std::cout << "/* message Nodes */" << std::endl;
        /* Nodes */
        std::vector<node_type> nodes(nodes_size);
        for (size_t i = 0; i < nodes_size; i++)
            nodes[i] = node_type(point_identifier<2>(i));

        storage->nodes = std::move(nodes);

        //std::cout << "/* message Edges */" << std::endl;
        /* Edges */
        /* Make the vector containing the edges */
        std::vector<edge_type> edges;
        edges.reserve(m_edges.size());
        for (size_t i = 0; i < m_edges.size(); i++)
        {
            auto node1 = typename node_type::id_type(m_edges[i][0]);
            auto node2 = typename node_type::id_type(m_edges[i][1]);

            if (node2 < node1)
                std::swap(node1, node2);

            auto e = edge_type{{node1, node2}};
            e.set_point_ids(m_edges[i].begin(), m_edges[i].begin() + 2); /* XXX: crap */
            edges.push_back(e);

            auto pts_ids  = e.point_ids();
        }
        /* Sort them */
        std::sort(edges.begin(), edges.end());

        /* Detect which ones are boundary edges */
        storage->boundary_edges.resize(m_edges.size());
        storage->boundary_info.resize(m_edges.size());
        for (size_t i = 0; i < m_boundary_edges.size(); i++)
        {
            auto node1 = typename node_type::id_type(m_boundary_edges[i][0]);
            auto node2 = typename node_type::id_type(m_boundary_edges[i][1]);

            if (node2 < node1)
                std::swap(node1, node2);

            auto e = edge_type{{node1, node2}};
            auto position = find_element_id(edges.begin(), edges.end(), e);
            auto edge = *next(edges.begin(), position.second);
            if (position.first == false)
            {
                std::cout << "Bad bug at " << __FILE__ << "("
                          << __LINE__ << ")" << std::endl;
                return false;
            }

            auto bnd_number = set_bnd_number(edge, storage);
            bnd_info bi{bnd_number, true};
            storage->boundary_info.at(position.second)  = bi;
            storage->boundary_edges.at(position.second) = true;
        }
        //std::cout << "/* boundary_edges */" << std::endl;
        for (size_t i = 0; i < m_edges.size(); i++)
        {
            bool be = storage->boundary_edges.at(i);
        //    std::cout<<" boundary_edges("<< i <<")="<< be<<std::endl;
        }
        storage->edges = std::move(edges);

        //std::cout << "/* message Surface */" << std::endl;
        /* Surfaces */
        std::vector<surface_type> surfaces;
        surfaces.reserve( m_triangles.size() +
                          m_quadrangles.size() +
                          m_pentagons.size() +
                          m_hexagons.size()  +
                          m_heptagons.size() +
                          m_octagons.size()  +
                          m_enneagons.size() +
                          m_decagons.size()  +
                          m_hendecagons.size() +
                          m_dodecagons.size() +
                          m_triadecagons.size()  +
                          m_pentadecagons.size() +
                          m_tesseradecagons.size())
                          ;

        //std::cout << "/* message Put Polygons*/" << std::endl;
        //std::cout << "/* triangles */" << std::endl;
        put_polygons(msh, m_triangles, surfaces);
        m_triangles.clear();
        //std::cout << "/* quadrangles */" << std::endl;

        put_polygons(msh, m_quadrangles, surfaces);
        m_quadrangles.clear();

        put_polygons(msh, m_pentagons, surfaces);
        m_pentagons.clear();

        put_polygons(msh, m_hexagons, surfaces);
        m_hexagons.clear();

        put_polygons(msh, m_heptagons, surfaces);
        m_heptagons.clear();

        put_polygons(msh, m_octagons, surfaces);
        m_octagons.clear();

        put_polygons(msh, m_enneagons, surfaces);
        m_enneagons.clear();

        put_polygons(msh, m_decagons, surfaces);
        m_decagons.clear();

        put_polygons(msh, m_hendecagons, surfaces);
        m_hendecagons.clear();

        put_polygons(msh, m_dodecagons, surfaces);
        m_dodecagons.clear();

        put_polygons(msh, m_triadecagons, surfaces);
        m_triadecagons.clear();

        put_polygons(msh, m_tesseradecagons, surfaces);
        m_tesseradecagons.clear();

        put_polygons(msh, m_pentadecagons, surfaces);
        m_pentadecagons.clear();

        index_transf(surfaces);

        std::sort(surfaces.begin(), surfaces.end());

        for(size_t i = 0; i < surfaces.size(); i++)
        {
            auto id = typename surface_type::id_type(i);
            (surfaces.at(i)).set_element_id(id);
        }

        storage->surfaces = std::move(surfaces);

        /* Detect which surfaces have hanging nodes */
        storage->special_surfaces.resize(storage->surfaces.size());
        for (size_t i = 0; i < storage->surfaces.size(); i++)
        {
            auto s = storage->surfaces.at(i);
            auto p = is_special_polygon(msh, s);

            storage->special_surfaces.at(i) = p;
        }
        //std::cout << "/* message Print stats*/" << std::endl;
        /* Print stats */
        storage->statistics();

        return true;
    }
};



template<typename T, size_t N>
class netgen_mesh_loader
{
    static_assert(N == 3, "netgen is a 3D-only mesh format");
};

namespace priv {

template<typename T>
std::tuple<T, T, T>
read_point_line(const char *str, char **endptr, T scalefactor)
{
    T t1, t2, t3;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1*scalefactor, t2*scalefactor, t3*scalefactor);
}

template<typename T>
std::tuple<T, T, T, T, T>
read_tetrahedron_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4, t5;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);
    t5 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2-1, t3-1, t4-1, t5-1);
}

template<typename T>
std::tuple<T, T, T, T>
read_triangle_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4;

    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);

    return std::make_tuple(t1, t2-1, t3-1, t4-1);
}

template<typename T>
void
sort_uniq(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    auto uniq_iter = std::unique(v.begin(), v.end());
    v.erase(uniq_iter, v.end());
}

} //namespace priv

#define THREADED
#ifdef THREADED
    #define THREAD(name, body) std::thread name([&]{body})
    #define WAIT_THREAD(name) name.join()
#else
    #define THREAD(name, body) {body}
    #define WAIT_THREAD(name)
#endif

template<typename T>
class netgen_mesh_loader<T,3> : public mesh_loader<simplicial_mesh<T,3>>
{
    typedef simplicial_mesh<T,3>                    mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::surface_type        surface_type;
    typedef typename mesh_type::volume_type         volume_type;

    std::vector<point_type>                         points;
    std::vector<node_type>                          nodes;
    std::vector<edge_type>                          edges;
    std::vector<surface_type>                       surfaces, boundary_surfaces;
    std::vector<volume_type>                        volumes;


    bool netgen_read(const std::string& filename)
    {
        /* Open file */
        if (filename.size() == 0)
        {
            std::cout << "Invalid mesh file name" << std::endl;
            return false;
        }

        size_t lines, linecount;

        mapped_file mf(filename);

        //std::cout << green << " * * * Reading NETGEN format mesh * * * ";
        //std::cout << nocolor << std::endl;

        /************************ Read points ************************/
        linecount = 0;

        const char *data = mf.mem();
        char *endptr;

        lines = strtot<size_t>(data, &endptr);

        points.reserve(lines);
        nodes.reserve(lines);

        while (linecount < lines)
        {
            if ( (linecount%100000) == 0 )
            {
                std::cout << "Reading points: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_point_line<T>(endptr, &endptr, 1.0);

            auto point = point_type( { std::get<0>(t),
                                       std::get<1>(t),
                                       std::get<2>(t) } );

            points.push_back( point );

            auto point_id = point_identifier<3>( linecount );
            auto node = node_type( { point_id } );

            nodes.push_back(node);
            /* Do something with that point */

            linecount++;
        }

        std::cout << "Reading points: " << linecount;
        std::cout << "/" << lines << std::endl;

        /************************ Read tetrahedra ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        edges.reserve(lines*6);
        surfaces.reserve(lines*4);
        volumes.reserve(lines);

        while (linecount < lines)
        {
            if ( (linecount%100000) == 0 )
            {
                std::cout << "Reading tetrahedra: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_tetrahedron_line<size_t>(endptr, &endptr);

            point_identifier<3>     p0(std::get<1>(t));
            point_identifier<3>     p1(std::get<2>(t));
            point_identifier<3>     p2(std::get<3>(t));
            point_identifier<3>     p3(std::get<4>(t));
            //domain_id_type      d(std::get<0>(t));

            edges.push_back( edge_type( { p0, p1 } ) );
            edges.push_back( edge_type( { p0, p2 } ) );
            edges.push_back( edge_type( { p0, p3 } ) );
            edges.push_back( edge_type( { p1, p2 } ) );
            edges.push_back( edge_type( { p1, p3 } ) );
            edges.push_back( edge_type( { p2, p3 } ) );

            surfaces.push_back( surface_type( { p0, p1, p2 } ) );
            surfaces.push_back( surface_type( { p0, p1, p3 } ) );
            surfaces.push_back( surface_type( { p0, p2, p3 } ) );
            surfaces.push_back( surface_type( { p1, p2, p3 } ) );

            //auto tuple = std::make_tuple(volume_type(p0, p1, p2, p3), d);
            //temp_tet.push_back( tuple );

            volumes.push_back( volume_type( { p0, p1, p2, p3 } ) );

            linecount++;
        }

        std::cout << "Reading tetrahedra: " << linecount;
        std::cout << "/" << lines << std::endl;

        /************************ Read boundary surfaces ************************/
        linecount = 0;

        lines = strtot<size_t>(endptr, &endptr);

        boundary_surfaces.reserve(lines);

        while (linecount < lines)
        {
            if ( (linecount%50000) == 0 )
            {
                std::cout << "Reading triangle: " << linecount;
                std::cout << "/" << lines << "\r";
                std::cout.flush();
            }

            auto t = priv::read_triangle_line<size_t>(endptr, &endptr);

            point_identifier<3>     p0(std::get<1>(t));
            point_identifier<3>     p1(std::get<2>(t));
            point_identifier<3>     p2(std::get<3>(t));

            surface_type   tri( { p0, p1, p2 } );

            boundary_surfaces.push_back( tri );

            linecount++;
        }

        std::cout << "Reading triangle: " << linecount;
        std::cout << "/" << lines << std::endl;

        return true;
    }

public:
    netgen_mesh_loader() = default;

    bool read_mesh(const std::string& s)
    {
        std::cout << " *** READING NETGEN MESH ***" << std::endl;
        return netgen_read(s);
    }

    bool populate_mesh(mesh_type& msh)
    {
        auto storage = msh.backend_storage();

        std::cout << "Sorting data...";
        std::cout.flush();

        storage->points = std::move(points);
        storage->nodes = std::move(nodes);

        /* sort edges, make unique and move them in geometry */
        THREAD(edge_thread,
            priv::sort_uniq(edges);
            storage->edges = std::move(edges);
        );

        /* sort triangles, make unique and move them in geometry */
        THREAD(tri_thread,
            priv::sort_uniq(surfaces);
            storage->surfaces = std::move(surfaces);
        );

        /* sort tetrahedra, make unique and move them in geometry */
        THREAD(tet_thread,
            std::sort(volumes.begin(), volumes.end());
            storage->volumes = std::move(volumes);
        );

        /* wait for the threads */
        WAIT_THREAD(edge_thread);
        WAIT_THREAD(tri_thread);
        WAIT_THREAD(tet_thread);

        storage->boundary_surfaces.resize(storage->surfaces.size());
        for (auto& bs : boundary_surfaces)
        {
            auto position = find_element_id(storage->surfaces.begin(),
                                            storage->surfaces.end(), bs);
            if (position.first == false)
            {
                std::cout << "Bad bug at " << __FILE__ << "("
                          << __LINE__ << ")" << std::endl;
                return false;
            }

            storage->boundary_surfaces[position.second] = true;
        }

        std::cout << "done." << std::endl;

        std::cout << "Nodes: " << storage->nodes.size() << std::endl;
        std::cout << "Edges: " << storage->edges.size() << std::endl;
        std::cout << "Faces: " << storage->surfaces.size() << std::endl;
        std::cout << "Volumes: " << storage->nodes.size() << std::endl;

        boundary_surfaces.clear();

        return true;
    }

};






} // namespace disk
