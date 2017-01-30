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

#ifndef _QUADRATURES_HPP_WAS_INCLUDED_
    #error "You must NOT include this file. Include quadratures.hpp"
#endif

#ifndef _QUAD_GENERIC_HPP_
#define _QUAD_GENERIC_HPP_

namespace disk {

template<typename T>
class quadrature<generic_mesh<T,2>, typename generic_mesh<T,2>::cell>
{
    size_t                                          m_order;
    std::vector<std::pair<point<T,2>, T>>           m_quadrature_data;

public:
    typedef generic_mesh<T,2>                       mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef quadrature_point<T,2>                   quadpoint_type;
    typedef typename mesh_type::point_type          point_type;
    typedef T                                       weight_type;

private:

    template<typename PtA>
    std::vector<quadpoint_type>
    integrate_triangle(const mesh_type& msh, const cell_type& cl,
                       const PtA& pts) const
    {
        std::vector<quadpoint_type> ret;

        ret.resize( m_quadrature_data.size() );

        auto col1 = pts[1] - pts[0];
        auto col2 = pts[2] - pts[0];

        /* Compute the area of the sub-triangle */
        auto tm = (col1.x()*col2.y() - col2.x()*col1.y())/2.;

        auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {
            auto point = col1*qd.first.x() + col2*qd.first.y() + pts[0];
            auto weight = qd.second * std::abs(tm);
            return make_qp(point, weight);
        };

        auto retbegin = ret.begin();

        std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                       retbegin, tr);

        return ret;
    }


    template<typename PtA>
    std::vector<quadpoint_type>
    integrate_quad(const mesh_type& msh, const cell_type& cl,
                       const PtA& pts) const
    {
        std::vector<quadpoint_type> ret;

        ret.resize( m_quadrature_data.size() * 2 );
        for (size_t i = 1; i < 3; i++)
        {
            auto pt1 = pts[i];
            auto pt2 = pts[i+1];
            auto col1 = pt1 - pts[0];
            auto col2 = pt2 - pts[0];

            /* Compute the area of the sub-triangle */
            auto tm = (col1.x()*col2.y() - col2.x()*col1.y())/2.;

            auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {
                auto point = col1*qd.first.x() + col2*qd.first.y() + pts[0];
                auto weight = qd.second * std::abs(tm);
                return make_qp(point, weight);
            };

            auto retbegin = ret.begin();
            std::advance(retbegin, m_quadrature_data.size()*(i-1));

            std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                           retbegin, tr);
        }

        return ret;
    }
    template<typename PtA>
    std::vector<quadpoint_type>
    integrate_quad_symmetric(const mesh_type& msh, const cell_type& cl,
                       const PtA& pts) const
    {
        std::vector<std::pair<point<T,1>, T>> m_quadrature_data_sqr;
        m_quadrature_data_sqr = edge_quadrature<T>(m_order);

        auto meas       = measure(msh, cl);

        assert(pts.size() == 4);

        size_t qd_size = m_quadrature_data_sqr.size();
        //std::cout << "qd_size = "<< qd_size << std::endl;
        //std::cout << "measure = "<< meas << std::endl;
        std::vector<std::pair<point<T,2>, T>> qd;
        qd.reserve(qd_size * qd_size);

        for (size_t i = 0; i < qd_size; i++)
        {
            auto px = m_quadrature_data_sqr[i].first;
            auto wx = m_quadrature_data_sqr[i].second;

            for (size_t j = 0; j < qd_size; j++)
            {
                auto py = m_quadrature_data_sqr[j].first;
                auto wy = m_quadrature_data_sqr[j].second;

                auto pt = point<T,2>({px[0], py[0]});
                auto wt = wx * wy;

                qd.push_back( std::make_pair(pt, wt) );
            }
        }

        auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {

            auto xi  = qd.first.x();
            auto eta = qd.first.y();

            auto point =(pts[0]*( 1. - xi) * (1. - eta) +
                         pts[3]*( 1. - xi) * eta  +
                         pts[1]* xi  * (1. - eta) +
                         pts[2]* xi  *  eta );

            auto weight = qd.second * meas;
            return make_qp(point, weight);
        };
        std::vector<quadpoint_type> ret;
        ret.resize(qd_size * qd_size);
        std::transform(qd.begin(), qd.end(), ret.begin(), tr);
        return ret;
    }
//#define OPTIMAL_TRIANGLE_NUMBER

    /* The 'optimal triangle number' version gives almost the same results
     * of the other version. In bigger meshes there are some advantages in
     * assembly time. The problem with this version is that it could generate
     * triangles with a very bad aspect ratio and at the moment I don't know
     * if and how this can affect computations, so I leave it turned off.
     */
#ifdef OPTIMAL_TRIANGLE_NUMBER
    template<typename PtA>
    std::vector<quadpoint_type>
    integrate_other(const mesh_type& msh, const cell_type& cl,
                    const PtA& pts) const
    {
        std::vector<quadpoint_type> ret;

        /* Break the cell in triangles, compute the transformation matrix and
         * map quadrature data in the physical space. Edges of the triangle as
         * column vectors in the transformation matrix. */
        ret.resize( m_quadrature_data.size() * pts.size()-2 );
        for (size_t i = 1; i < pts.size()-1; i++)
        {
            auto pt1 = pts[i];
            auto pt2 = pts[i+1];
            auto col1 = pt1 - pts[0];
            auto col2 = pt2 - pts[0];

            /* Compute the area of the sub-triangle */
            auto tm = (col1.x()*col2.y() - col2.x()*col1.y())/2.;

            auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {
                auto point = col1*qd.first.x() + col2*qd.first.y() + pts[0];
                auto weight = qd.second * std::abs(tm);
                return make_qp(point, weight);
            };

            auto retbegin = ret.begin();
            std::advance(retbegin, m_quadrature_data.size()*(i-1));

            std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                           retbegin, tr);
        }

        return ret;
    }
#else
    template<typename PtA>
    std::vector<quadpoint_type>
    integrate_other(const mesh_type& msh, const cell_type& cl,
                    const PtA& pts) const
    {
        auto c_center   = barycenter(msh, cl);

        std::vector<quadpoint_type> ret;;

        /* Break the cell in triangles, compute the transformation matrix and
         * map quadrature data in the physical space. Edges of the triangle as
         * column vectors in the transformation matrix. */
        ret.resize( m_quadrature_data.size() * pts.size() );
        for (size_t i = 0; i < pts.size(); i++)
        {
            auto pt1 = pts[i];
            auto pt2 = pts[(i+1)%pts.size()];
            auto col1 = pt1 - c_center;
            auto col2 = pt2 - c_center;

            /* Compute the area of the sub-triangle */
            auto tm = (col1.x()*col2.y() - col2.x()*col1.y())/2.;

            auto tr = [&](const std::pair<point<T,2>, T>& qd) -> auto {
                auto point = col1*qd.first.x() + col2*qd.first.y() + c_center;
                auto weight = qd.second * std::abs(tm);
                return make_qp(point, weight);
            };

            auto retbegin = ret.begin();
            std::advance(retbegin, m_quadrature_data.size()*i);

            std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                           retbegin, tr);
        }

        return ret;
    }
#endif
    template<typename N>
    void
    sort_uniq(std::vector<N>& v)
    {
        std::sort(v.begin(), v.end());
        auto uniq_iter = std::unique(v.begin(), v.end());
        v.erase(uniq_iter, v.end());
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

public:
    quadrature()
        : m_order(1)
    {
        m_quadrature_data = triangle_quadrature(1);
    }

    quadrature(size_t order)
        : m_order(order)
    {
        m_quadrature_data = triangle_quadrature(m_order);
    }

    std::vector<quadpoint_type>
    integrate(const mesh_type& msh, const cell_type& cl) const
    {
        auto pts      = points(msh, cl);
        std::vector<point<T,2>> vertices;

        /*WK:  This is ok if there are hanging nodes, otherwise it would be
                unnecessary and costly. One alternative is to bring the variable
                "hanging nodes" and do:

                if (num_pts > 3 & hanging_nodes)
                    vertices   = verify_polygon(msh, cl, pts);
                else
                    vertices   = pts;
                    else
        */


        if(pts.size() > 3) // &  hanging nodes)
        {
            auto pair = is_special_polygon(msh, cl, pts);
            if(pair.first)
                vertices = pair.second;
        }
            vertices = pts;

        switch(vertices.size())
        {
            case 3:
                return integrate_triangle(msh, cl, vertices);

            case 4:
                return integrate_quad_symmetric(msh, cl, vertices);

            default:
                return integrate_other(msh, cl, vertices);
        }

        throw std::logic_error("Shouldn't have arrived here");
    }
};

template<typename T>
class quadrature<generic_mesh<T,2>, typename generic_mesh<T,2>::face>
{
    size_t                                          m_order;
    std::vector<std::pair<point<T,1>, T>>           m_quadrature_data;

public:
    typedef generic_mesh<T,2>                       mesh_type;
    typedef typename mesh_type::face                face_type;
    typedef quadrature_point<T,2>                   quadpoint_type;
    typedef typename mesh_type::point_type          point_type;
    typedef T                                       weight_type;

    quadrature()
        : m_order(1)
    {
        m_quadrature_data = edge_quadrature<T>(1);
    }

    quadrature(size_t order)
        : m_order(order)
    {
        m_quadrature_data = edge_quadrature<T>(m_order);
    }

    std::vector<quadpoint_type>
    integrate(const mesh_type& msh, const face_type& fc) const
    {
        auto meas = measure(msh, fc);
        auto pts = points(msh, fc);
        auto tr = [&](const std::pair<point<T,1>, T>& qd) -> auto {
            auto point = (pts[1] - pts[0])*qd.first.x() + pts[0];
            auto weight = qd.second * meas;
            return make_qp(point, weight);
        };

        std::vector<quadpoint_type> ret(m_quadrature_data.size());
        std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                       ret.begin(), tr);

        return ret;
    }
};

template<typename T>
class quadrature<generic_mesh<T,1>, typename generic_mesh<T,1>::cell>
{
    size_t                                          m_order;
    std::vector<std::pair<point<T,1>, T>>           m_quadrature_data;

public:
    typedef generic_mesh<T,1>                       mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef quadrature_point<T,1>                   quadpoint_type;
    typedef typename mesh_type::point_type          point_type;
    typedef T                                       weight_type;

    quadrature()
        : m_order(1)
    {
        m_quadrature_data = edge_quadrature<T>(1);
    }

    quadrature(size_t order)
        : m_order(order)
    {
        m_quadrature_data = edge_quadrature<T>(m_order);
    }

    std::vector<quadpoint_type>
    integrate(const mesh_type& msh, const cell_type& cl) const
    {
        auto pts        = points(msh, cl);
        auto meas       = measure(msh, cl);

        assert(pts.size() == 2);

        std::vector<quadpoint_type> ret;;
        ret.resize( m_quadrature_data.size() );


        auto tr = [&](const std::pair<point<T,1>, T>& qd) -> auto {
            auto point = (pts[1] - pts[0]) * qd.first.x() + pts[0];
            auto weight = qd.second * meas;
            return make_qp(point, weight);
        };

        auto retbegin = ret.begin();

        std::transform(m_quadrature_data.begin(), m_quadrature_data.end(),
                       retbegin, tr);

        return ret;
    }
};

template<typename T>
class quadrature<generic_mesh<T,1>, typename generic_mesh<T,1>::face>
{

public:
    typedef generic_mesh<T,1>                       mesh_type;
    typedef typename mesh_type::face                face_type;
    typedef quadrature_point<T,1>                   quadpoint_type;
    typedef typename mesh_type::point_type          point_type;
    typedef T                                       weight_type;

    quadrature()
    {}

    quadrature(size_t)
    {}

    std::vector<quadpoint_type>
    integrate(const mesh_type& msh, const face_type& fc) const
    {
        return std::vector<quadpoint_type>( 1, make_qp(barycenter(msh, fc), 1.) );
    }
};

} // namespace disk

#endif /* _QUAD_GENERIC_HPP_ */
