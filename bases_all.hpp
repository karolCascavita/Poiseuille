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

#pragma once

#include <vector>
#include "geometry/geometry.hpp"
#include "common/eigen.hpp"
#include "bases/bases_bones.hpp"
#include "bases/bases_templates.hpp"
#include "bases/monomial_generator.hpp"

#define VERY_HIGH_DEGREE 9999999


namespace disk {

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::cell>
    : public priv::scaled_monomial_basis_base<2>
{
    typedef Mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef typename mesh_type::point_type          point_type;
    typedef priv::scaled_monomial_basis_base<2>     base;

public:
    typedef T                                       function_value_type;
    typedef static_vector<T,2>                      gradient_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    Eigen::Matrix<T, Eigen::Dynamic, 1>
    eval_laplacians(const mesh_type& msh, const cell_type& cl, const point_type& pt,
                        size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h   = diameter(msh, cl);

        auto ep = (pt - bar)/h;
        Eigen::Matrix<T, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());

        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;

            auto px  = iexp_pow(ep.x(), m[0]);
            auto py  = iexp_pow(ep.y(), m[1]);

            auto ddx = (m[0] < 2) ? 0 : (m[0]/h) * ((m[0]-1)/h) * iexp_pow(ep.x(), m[0]-2);
            auto ddy = (m[1] < 2) ? 0 : (m[1]/h) * ((m[1]-1)/h) * iexp_pow(ep.y(), m[1]-2);

            ret(i++) = ddx * py + ddy * px;
        }
        return ret;
    }

    Eigen::Matrix< T, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const cell_type& cl,
                   const point_type& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);

        auto ep = (pt - bar)/h;
        Eigen::Matrix< T, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;
            auto vx = iexp_pow(ep.x(), m[0]);
            auto vy = iexp_pow(ep.y(), m[1]);
            ret(i++) = vx * vy;
        }

        return ret;
    }

    Eigen::Matrix< T, Eigen::Dynamic, 2>
    eval_gradients(const mesh_type& msh, const cell_type& cl,
                   const point_type& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);
        auto ih = 1./h;

        auto ep = (pt - bar)/h;

        Eigen::Matrix< T, Eigen::Dynamic, 2> ret;
        ret.resize( eval_range.size(), 2 );

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;
            gradient_value_type grad;

            auto px = iexp_pow(ep.x(), m[0]);
            auto py = iexp_pow(ep.y(), m[1]);

            auto dx = (m[0] == 0) ? 0 : m[0]*ih*iexp_pow(ep.x(), m[0]-1);
            auto dy = (m[1] == 0) ? 0 : m[1]*ih*iexp_pow(ep.y(), m[1]-1);
            ret(i,0) = dx * py;
            ret(i,1) = dy * px;
            i++;
        }

        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::face>
    : public priv::scaled_monomial_basis_base<1>
{
    typedef Mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::scalar_type         scalar_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::face                face_type;
    typedef priv::scaled_monomial_basis_base<1>     base;

public:
    typedef T                                       function_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    scaled_monomial_scalar_basis(size_t degree, size_t computed_degree)
        : base(degree, computed_degree)
    {}

    Eigen::Matrix<T, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const face_type& fc,
                   const point_type& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto pts = points(msh, fc);
        auto bar = barycenter(msh, fc);
        auto h = (pts[1] - pts[0]).to_vector().norm()/2.0;
        auto v = (bar - pts[0]).to_vector();
        auto t = (pt - bar).to_vector();
        T dot = v.dot(t);
        auto ep = point<T, 1>({dot/(h*h)});

        Eigen::Matrix<T, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;
            auto vx = iexp_pow(ep.x(), m[0]);
            ret(i++) = vx;
        }

        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,1,Storage>, typename Mesh<T,1,Storage>::cell>
    : public priv::scaled_monomial_basis_base<1>
{
    typedef Mesh<T,1,Storage>                       mesh_type;
    typedef typename mesh_type::cell                cell_type;
    typedef typename mesh_type::point_type          point_type;
    typedef priv::scaled_monomial_basis_base<1>     base;

public:
    typedef T                                       function_value_type;
    typedef T                                       gradient_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    Eigen::Matrix<T, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const cell_type& cl,
                   const point_type& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h   = diameter(msh, cl);
        auto ep  = (pt - bar)/h;

        Eigen::Matrix< T, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;
            auto vx = iexp_pow(ep.x(), m[0]);
            ret(i++) = vx ;
        }

        return ret;
    }
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    eval_gradients(const mesh_type& msh, const cell_type& cl,
               const point_type& pt,
               size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h   = diameter(msh, cl);
        auto ep  = (pt - bar)/h;

        Eigen::Matrix<T, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto  m  = *itor;
            auto  dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);
            ret(i++) = dx;
        }

        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,1,Storage>, typename Mesh<T,1,Storage>::face>
    : public priv::scaled_monomial_basis_base<0>
{
    typedef Mesh<T,1,Storage>                           mesh_type;
    typedef typename mesh_type::point_type              point_type;
    typedef typename mesh_type::face                    face_type;
    typedef priv::scaled_monomial_basis_base<0>         base;

public:
    typedef T                                           function_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    Eigen::Matrix<T, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const face_type& fc, const point_type& pt,
                size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const

    {
        return Eigen::Matrix<T, 1, 1>::Ones(1,1);
    }
};


template<template<typename, size_t, typename> class Mesh,
                  typename T, typename Storage>
std::vector<point<T,2>>
make_test_points(const Mesh<T,2,Storage>& msh,
                 const typename Mesh<T,2,Storage>::cell& cl)
{
    std::vector<point<T,2>> test_points;
    auto pts = points(msh, cl);
    auto bar = barycenter(msh, cl);

    test_points.insert(test_points.begin(), pts.begin(), pts.end()); //vertices
    test_points.push_back(bar); //barycenter

    //split in triangles and take barycenter
    for (size_t i = 0; i < pts.size(); i++)
        test_points.push_back( (pts[i] + pts[(i+1) % pts.size()] + bar) / 3. );

    return test_points;
}

template<template<typename, size_t, typename> class Mesh,
                  typename T, typename Storage>
std::vector<point<T,2>>
make_test_points(const Mesh<T,2,Storage>& msh,
                 const typename Mesh<T,2,Storage>::face& fc)
{
    std::vector<point<T,2>> test_points(5);
    auto pts = points(msh, fc);
    assert(pts.size() == 2); /* we are in 2D and we are dealing with a face,
                                which is an edge */

    test_points[0] = pts[0];
    test_points[1] = pts[1];
    test_points[2] = (pts[0] + pts[1]) / 2.;
    test_points[3] = (test_points[2] + pts[0]) / 2.;
    test_points[4] = (test_points[2] + pts[1]) / 2.;

    return test_points;
}



template<template<typename, size_t, typename> class Mesh,
                  typename T, typename Storage>
std::vector<point<T,1>>
make_test_points(const Mesh<T,1,Storage>& msh,
                 const typename Mesh<T,1,Storage>::cell& cl)
{
    std::vector<point<T,1>> test_points;
    auto pts = points(msh, cl);
    auto bar = barycenter(msh, cl);

    test_points.insert(test_points.begin(), pts.begin(), pts.end()); //vertices
    test_points.push_back(bar); //barycenter
    test_points.push_back((pts[0] + bar) / 2.);
    test_points.push_back((pts[1] + bar) / 2.);

    return test_points;
}


template<template<typename, size_t, typename> class Mesh,
                  typename T, typename Storage>
std::vector<point<T,2>>
make_test_points(const Mesh<T,2,Storage>& msh,
                 const typename Mesh<T,2,Storage>::face& fc,
                 const size_t degree)
{
    std::vector<point<T,2>> test_points(degree + 1);
    auto pts = points(msh, fc);
    assert(pts.size() == 2); /* we are in 2D and we are dealing with a face,
                                which is an edge */

    test_points[0] = pts[0];

    point<T,2> h;
    if(degree > 0)
        h = (pts[0] + pts[1]) / degree;

    for(size_t i = 1; i < degree + 1; i++)
        test_points[i] = pts[0] + i * h;

    return test_points;
}
template<template<typename, size_t, typename> class Mesh,
                  typename T, typename Storage>
std::vector<point<T,2>>
make_test_points(const Mesh<T,2,Storage>& msh,
                 const typename Mesh<T,2,Storage>::cell& cl,
                 const size_t degree)
{
    std::vector<point<T,2>> test_points;
    auto pts = points(msh, cl);
    auto bar = barycenter(msh, cl);
    auto vts = msh.get_vertices(cl, pts);

    test_points.insert(test_points.begin(), vts.begin(), vts.end()); //vertices

    //split in triangles and take barycenter
    for (size_t i = 0; i < pts.size(); i++)
    {
        auto h = ( pts[i] + pts[(i+1) % pts.size()] )/ degree ;
        // We use just degree since the numeration starts with vertices and
        // it is not linear as in faces
        for(size_t j = 1; j < degree ; i++)
            test_points.push_back( pts[i]  + j * h );
    }

    if(degree == 3)
        test_points.push_back(bar); //barycenter

    if(degree > 3)
        throw std::invalid_argument("Make_test_points not available for degree > 3");

    return test_points;
}


} // namespace disk
