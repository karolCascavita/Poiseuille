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
struct mesh_parameters
{
    mesh_parameters()
    {
        start_exact = 0;
        mesh_name   = 0;
        num_remesh  = 0;
        marker_name = 2;
        initial_imsh = 0;
        percent     = 0.1;
        recycle     = true;
        diff        = false;
        mark_all    = false; // There is a problem here read comment
        hanging_nodes = true;
        call_mesher   = false;
    }
    // If marl_all ==true, it cannot be used  in circular for iter == 0, that used quadrilaterals
    // Since, we have triangles in the center, that will make hanging_nodes to appear
    // However, since the quadrilaterals are split in triangles and doesn't support
    // hanging nodes, the triangles made by quadrilaterals won-t see the hanging_nodes

    T       percent;
    int     mesh_name;
    int     num_remesh;
    bool    call_mesher;
    int     initial_imsh;
    int     marker_name;
    int     start_exact;
    bool    recycle;
    bool    hanging_nodes;
    bool    diff;
    bool    mark_all;
    std::string     short_mesh_name;
    std::string     directory;
    std::string     summary;
    std::string     summary_old;

    friend std::ostream& operator<<(std::ostream& os, const mesh_parameters<T>& mp) {
        os << "Mesh Parameters: "<<std::endl;
        os << "* percent      : "<< mp.percent<< std::endl;
        os << "* mesh_name    : "<< mp.mesh_name<< std::endl;
        os << "* num_remesh   : "<< mp.num_remesh<< std::endl;
        os << "* call_mesher  : "<< mp.call_mesher<< std::endl;
        os << "* initial_imsh : "<< mp.initial_imsh<< std::endl;
        os << "* marker_name  : "<< mp.marker_name<< std::endl;
        os << "* recycle      : "<< mp.recycle<< std::endl;
        os << "* hanging_nodes: "<< mp.hanging_nodes<< std::endl;
        os << "* diffusion    : "<< mp.diff<< std::endl;
        os << "* mark_all     : "<< mp.mark_all<< std::endl;
        os << "* short_mesh_name: "<< mp.short_mesh_name<< std::endl;
        os << "* directory    : "<< mp.directory<< std::endl;
        os << "* summary      : "<< mp.summary<< std::endl;
        os << "* summary_old  : "<< mp.summary_old<< std::endl;
        return os;
    }
};
template<typename MeshType, typename CellType, typename FaceType>
size_t
face_position(const MeshType& msh, const CellType& cell, const FaceType& face)
{
    auto c_faces = faces(msh, cell);
    size_t   j = 0;
    for(auto& fc : c_faces)
    {
        if(fc == face)
            return j;
        j++;
    }
    throw std::invalid_argument("This is a bug: face not found");
}

template<typename T, size_t DIM, typename Storage>
std::vector<dynamic_vector<T>>
solution_zero_vector(const disk::mesh<T,DIM,Storage>& msh, const size_t degree)
{
    typedef disk::mesh<T,DIM,Storage>              mesh_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, typename mesh_type::cell>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, typename mesh_type::face>    face_basis_type;

    std::vector<dynamic_vector<T>> U(msh.cells_size());
    size_t i = 0;
    for (auto& cl : msh)
    {
        auto fcs = faces(msh, cl);
        auto num_faces       = fcs.size();
        auto num_cell_dofs   = cell_basis_type(degree).size();
        auto num_face_dofs   = face_basis_type(degree).size();
        disk::dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);
        U[i] =  dynamic_vector<T>::Zero(dsr.total_size());
        ++i;
    }
    return U;
};

template<typename U, size_t N>
void
fvca5_save_tuples(std::ofstream& ofs,
                    const std::vector<std::array<U, N>>& tuples)
{
    switch(N)
    {
        case 3 :
            ofs << "triangles"<<std::endl;
            break;
        case 4 :
            ofs << "quadrangles"<<std::endl;
            break;
        case 5 :
            ofs << "pentagons"<<std::endl;
            break;
        case 6 :
            ofs << "hexagons"<<std::endl;
            break;
        case 7 :
            ofs << "heptagons"<<std::endl;
            break;
        case 8 :
            ofs << "octagons"<<std::endl;
            break;
        default:
            std::cout << "ADD MORE POLYGONS" << std::endl;
            break;
    }
    ofs << tuples.size()<<std::endl;
    for (auto& tuple : tuples)
    {
        for (size_t j = 0; j < N; j++)
        {
            ofs << tuple.at(j) + 1 << " ";
        }
        ofs <<std::endl;
    }
    //return true;
};



template <typename T>
void
save_data(const std::vector<T>& vec, const std::string& filename)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file"<<std::endl;

    ofs << vec.size()<<std::endl;

    for(auto& v :  vec)
        ofs << v <<std::endl;
};

template <typename T>
void
save_data(const std::vector<dynamic_vector<T>>& vec, const std::string& filename)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file"<<std::endl;

    ofs << vec.size()<<std::endl;

    for(auto& v :  vec)
    {
        ofs << v.size() << "  ";
        for (size_t i = 0; i < v.size(); i++)
            ofs << v(i) <<" ";

        ofs << std::endl;
    }
};

template <typename T>
void
save_data(const std::vector<dynamic_matrix<T>>& vec, const std::string& filename)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file"<<std::endl;

    ofs << vec.size()<<std::endl;
    for(auto& m :  vec)
    {
        ofs << m.rows()<< " " << m.cols() << "  ";
        for (size_t i = 0; i < m.rows(); i++)
        {
            for (size_t j = 0; j < m.cols();j++)
                ofs << m(i,j) <<" ";
        }
        ofs << std::endl;
    }
};

template <typename MeshType>
void
save_data(const tensors2<MeshType>& vec, const std::string& filename)
{
    typedef typename MeshType::scalar_type scalar_type;
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file"<<std::endl;

    ofs << vec.size()<<std::endl;

    for(auto& tsr :  vec)
    {
        dynamic_matrix<scalar_type> m = tsr.join_all();
        ofs << m.rows()<< " " << m.cols() << "  ";
        for (size_t i = 0; i < m.rows(); i++)
        {
            for (size_t j = 0; j < m.cols();j++)
                ofs << m(i,j) <<" ";
        }
        ofs << std::endl;
    }
};
template <typename T>
void
read_data(std::vector<T>& vec, const std::string& filename)
{

    size_t elements_to_read;
    std::ifstream ifs(filename);

    if (!ifs.is_open())
        std::cout << "Error opening file"<<std::endl;
    std::cout << "Opening file: "<< filename <<std::endl;

    ifs >> elements_to_read;
    std::cout << "Attempting to read " << elements_to_read << " values" << std::endl;

    vec = std::vector<T>(elements_to_read);

    for (size_t i = 0; i < elements_to_read; i++)
        ifs >> vec.at(i);

    return;
};
template <typename T>
void
read_data(std::vector<dynamic_vector<T>>& vec, const std::string& filename)
{

    size_t elements_to_read;
    size_t values_to_read;

    std::ifstream ifs(filename);

    if (!ifs.is_open())
        std::cout << "Error opening file: "<< filename <<std::endl;

    std::cout << "Opening file: "<< filename <<std::endl;
    ifs >> elements_to_read;
    std::cout << "Attempting to read " << elements_to_read << " values" << std::endl;

    vec = std::vector<dynamic_vector<T>>(elements_to_read);

    for (size_t i = 0; i < elements_to_read; i++)
    {
        ifs >> values_to_read;
        dynamic_vector<T>  x(values_to_read);
        //std::cout << " * x.size : "<< x.rows() << "x"<<x.cols() << std::endl;
        //std::cout << " * values_to_read : "<<values_to_read << std::endl;
        for(size_t j = 0; j < values_to_read; j++)
        {
            T val;
            ifs >> val;
            x(j,0) = val;
        }
        vec.at(i) = x;
    }
    return;
};
template<typename T>
void
read_data(std::vector<dynamic_matrix<T>>& vec, const std::string& filename)
{
    std::cout << "INSIDE read_data matrix" << std::endl;

    size_t elements_to_read;
    size_t cols, rows;

    std::ifstream ifs(filename);

    if (!ifs.is_open())
        std::cout << "Error opening file: "<< filename<<std::endl;
    std::cout << "Opening file: "<< filename <<std::endl;

    ifs >> elements_to_read;
    std::cout << "Attempting to read " << elements_to_read << " values" << std::endl;

    vec = std::vector<dynamic_matrix<T>>(elements_to_read);

    for (size_t icell = 0; icell < elements_to_read; icell++)
    {
        ifs >> rows >> cols;
        //std::cout << " * mat.size("<< icell <<") : "<< rows << "x"<< cols << std::endl;
        dynamic_matrix<T>  mat = dynamic_matrix<T>::Zero(rows, cols);

        for(size_t k = 0; k < rows * cols; k++)
        {
            T val;
            ifs >> val;
            auto i = k / cols;
            auto j = k % cols;
            assert(i < rows || j < cols);

            mat(i,j) = val;
        }
        vec.at(icell) = mat;
    }
    return;
};
void
read_data(std::vector<std::pair<size_t, size_t>>& vec,
          std::vector<size_t>& idx_transf_vec,
            const std::string& filename)
{
    std::cout << "INSIDE read_data pair" << std::endl;

    size_t elements_to_read;

    std::ifstream ifs(filename);

    if (!ifs.is_open())
        std::cout << "Error opening file: "<< filename<<std::endl;
    std::cout << "Opening file: "<< filename <<std::endl;

    ifs >> elements_to_read;
    std::cout << "Attempting to read " << elements_to_read << " values" << std::endl;

    vec = std::vector<std::pair<size_t, size_t>>(elements_to_read);
    idx_transf_vec = std::vector<size_t>(elements_to_read);

    for (size_t icell = 0; icell < elements_to_read; icell++)
    {
        size_t value, first, second;
        ifs >> value >>first >> second;
        vec.at(icell) = std::make_pair(first, second);
        idx_transf_vec.at(icell) = value;
    }
};
template<typename TensorsType, typename MeshType>
void
get_from_tensor(std::vector<punctual_tensor<MeshType>>& vec,
                const std::vector<TensorsType> & tsr_vec,
                const std::string& name)
{
    vec = std::vector<punctual_tensor<MeshType>>(tsr_vec.size());
    size_t val = 0;
    if(name == "sigma")
        val = 1;
    else if(name == "gamma")
        val = 2;
    else if(name == "xi_norm")
        val = 3;
    else
        std::cout << "WARNING: name not found in tensor matrix variables." << std::endl;
    //std::cout << "INSIDE get_from_tensor" << std::endl;
    //std::cout << " * vec size: "<< vec.size() << std::endl;
    //std::cout << " * value = "<<val << std::endl;
    size_t i = 0;
    switch (val)
    {
        case 1:
            for(auto& tsr : tsr_vec)
                vec.at(i++) = tsr.sigma;
            return;
        case 2:
            for(auto& tsr : tsr_vec)
                vec.at(i++) = tsr.gamma;
            return;
        case 3:
            for(auto& tsr : tsr_vec)
                vec.at(i++) = tsr.xi_norm;
            return;
        default:
            std::cout << "WARNING: Data was not saved." << std::endl;
            return;
    }
};

template<typename TensorsType>
void
get_from_tensor(std::vector<size_t>& vec,
                const std::vector<TensorsType>& tsr_vec,
                const std::string& name)
{
    vec = std::vector<size_t>(tsr_vec.size());
    if(name != "quad_degree")
    {
        std::cout << "WARNING: The only scalar variable in tensor structure is";
        std::cout << "quad_degree. Review name or function get_from_tensor.";
        std::cout << std::endl;
        return;
    }
    size_t i = 0;
    for(auto& tsr : tsr_vec)
        vec.at(i++) = tsr.sigma.m_quad_degree;
};

template<typename TensorsType, typename T>
void
put_tensor(std::vector<TensorsType>& tsr_vec,
           const std::vector<dynamic_matrix<T>>& vec,
           const std::string& name)
{
    int val;
    if(name == "sigma")
        val = 1;
    else if(name == "gamma")
        val = 2;
    else if(name == "xi_norm")
        val = 3;
    else
        std::cout << "WARNING: name not found in tensor matrix variables." << std::endl;

    size_t i = 0;

    switch(val)
    {
        case 1:
            for(auto mat : vec)
                tsr_vec.at(i++).sigma.split_all(mat, 2);
            return;
        case 2:
            for(auto mat : vec)
                tsr_vec.at(i++).gamma.split_all(mat, 2);
            return;
        case 3:
            for(auto mat : vec)
                tsr_vec.at(i++).xi_norm.split_all(mat, 1);
            return;
        default:
            throw std::invalid_argument("Not known name variable for tensor.");
            return;
    }
};
template<typename TensorsType, typename MeshType, typename T>
void
put_tensor(std::vector<TensorsType>& tsr_vec,
           const MeshType& msh,
           const std::vector<T>& vec,
           const std::string& name)
{
    if(name != "quad_degree")
    {
        std::cout << "WARNING: The only scalar variable in tensor structure is";
        std::cout << "quad_degree. Review name or function get_from_tensor.";
        std::cout << std::endl;
        return;
    }

    size_t i = 0;

    for(auto& cl : msh)
    {
        auto local_quad_degree = vec.at(i);
        tsr_vec.at(i).sigma = punctual_tensor<MeshType>(msh, cl, local_quad_degree);
        tsr_vec.at(i).gamma = punctual_tensor<MeshType>(msh, cl, local_quad_degree);
        tsr_vec.at(i).xi_norm = punctual_tensor<MeshType>(msh, cl, local_quad_degree);
        i++;
    }


    return;
};
#if 0
//WarningK: The first two load_data, should be just one. But I have to
//1 - Leave out index_transf and to be sure that I really dont need it.
//2 - Fix the thing of rmc and rm
//3 - step: imsh for levels_ancestors_vec
//    step: imsh-1 for Uh_Th and tsr since is info on the previous mesh
#endif

template<typename T, typename InputVectorType>
void
 load_data(InputVectorType& vec,
            const mesh_parameters<T>& mp,
            const std::string& name,
            const std::string& R,
            const size_t step)
{
     typedef dynamic_vector<T> vector_type;
     typedef dynamic_matrix<T> matrix_type;

     auto info_other = mp.summary_old + "_" + R + tostr(step) +".txt";
     auto filename   = mp.directory + name + info_other;
     std::cout << "namefile to open:  "<<  filename << std::endl;
     read_data(vec, filename);
 }
template<typename T>
void
load_data(std::vector<std::pair<size_t,size_t>>& levels_ancestors_vec,
                   std::vector<size_t>& index_transf,
                   const mesh_parameters<T>& mp,
                   const std::string& name,
                   const std::string& R,
                   const size_t step)
{
    typedef dynamic_vector<T> vector_type;
    typedef dynamic_matrix<T> matrix_type;

    auto info_other = mp.summary_old + "_" + R + tostr(step) +".txt";
    auto filename   = mp.directory + name + info_other;
    std::cout << "namefile to open:  "<<  filename << std::endl;
    read_data(levels_ancestors_vec, index_transf, filename);
}
#if 0
template< typename T>
void
load_data(std::vector<dynamic_vector<T>>& Uh_Th,
          const mesh_parameters<T>& mp,
           const size_t imsh)
{
    typedef dynamic_vector<T> vector_type;
    typedef dynamic_matrix<T> matrix_type;

    auto info = mp.summary_old + "_R" + tostr(imsh-1) +".txt";

    std::cout << "/* namefile to open:  "<< mp.directory + "/Uh" + info << std::endl;
    read_data(Uh_Th , mp.directory + "/Uh"  + info);

    std::cout << "Uh reading succesful" << std::endl;
}
#endif
template<typename TensorsType, typename MeshType, typename T>
void
load_data(std::vector<TensorsType>& tsr_vec,
        const MeshType& msh,
        const mesh_parameters<T>& mp,
        const std::string& R,
        const size_t step)

{
    typedef dynamic_vector<T> vector_type;
    typedef dynamic_matrix<T> matrix_type;

    std::vector<matrix_type> sigma;
    std::vector<T> quadeg;

    auto info = mp.summary_old + "_" + R + tostr(step-1) +".txt";

    //get_from_tensor(sigma,  tsr_vec, "sigma");
    //get_from_tensor(quadeg, tsr_vec, "quad_degree");

    read_data(sigma , mp.directory + "/Sigma" + info);
    read_data(quadeg, mp.directory + "/QuadDegree" + info);

    assert(sigma.size() == quadeg.size());
    if(sigma.size() != quadeg.size())
        throw std::logic_error("Sizes of loaded data doesn't match");

    tsr_vec = std::vector<TensorsType>(sigma.size());

    //Esto debe ir primero para que ya se tenga el grado y poder hacer split_all;
    put_tensor(tsr_vec, msh, quadeg, "quad_degree");
    put_tensor(tsr_vec, sigma,  "sigma");
};
template<typename T>
struct triangle
{
    std::array<point<T,2>, 3>  points;

    friend std::ostream& operator<<(std::ostream& os, const triangle<T>& t) {
        os << "Triangle: ";
        for (auto& p : t.points)
            os << p << " ";
        return os;
    }
};
template<typename T>
double dot(const point<T,2>& p1, const point<T,2>& p2)
{
    return p1.x()*p2.x() + p1.y()*p2.y();
};
template<typename PointType>
double distance(const PointType& p1, const PointType& p2)
{
    return sqrt(dot(p1, p2));
};
template<typename T>
bool is_inside(const triangle<T>& t, const point<T,2>& pt)
{
    auto pts = t.points;
    auto v0 = pts[1] - pts[0];
    auto v1 = pts[2] - pts[0];
    auto v2 = pt - pts[0];

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

    return (upos && vpos && uvpos);
};
template<typename T>
auto
trilinear_interpolation(const dynamic_matrix<T>& tensor,
                        const triangle<T>& tri , const point<T,2>& ep)
{
    typedef point<T,2>          point_type;
    typedef dynamic_vector<T>   vector_type;
    typedef dynamic_matrix<T>   matrix_type;

    // WARNING: This only works for triangles and linear Lagrange polynomials

    //tri then  3 vertices
    vector_type ret = vector_type::Zero(2);
    vector_type fun = vector_type::Zero(3);

    auto vts  = tri.points;
    auto triarea = [](const std::array<point<T,2>, 3>& pts) -> T
    {
        T acc{};
            for (size_t i = 1; i < 2; i++)
            {
                auto u = (pts[i] - pts[0]).to_vector();
                auto v = (pts[i+1] - pts[0]).to_vector();
                auto n = cross(u, v);
                acc += n.norm() / T(2);
            }
        return acc;
    };

    auto S = triarea(vts);

    for(size_t i = 0; i < 3; i++)
    {
        triangle<T> t;
        t.points[0] = vts.at((i + 1)% 3);
        t.points[1] = vts.at((i + 2)% 3);
        t.points[2] = ep;

        fun(i) = triarea(t.points) / S;
    }

    for(size_t j = 0; j < 2; j++)
    {
        ret(j) = fun.dot(tensor.row(j));
        assert(ret(j) < 1.e10);
    }

    return ret;
};

template<typename T,  typename Storage, typename CellType, typename PunctualTensor>
dynamic_vector<T>
eval_tensor_interpolation(const mesh<T,2, Storage>& msh,
                          const CellType          & cell,
                          const PunctualTensor    & tensor,
                          const point<T,2>        & ep)
{
    typedef dynamic_matrix<T>  matrix_type;
    dynamic_vector<T> ret = dynamic_vector<T>::Zero(2);
    auto is_find = false;
    auto pts   = points( msh, cell);
    auto fcs   = faces(  msh, cell);
    auto cell_barycenter = barycenter(msh, cell);

    for(size_t ifc = 0; ifc < pts.size(); ifc++)
    {
        auto face  = fcs.at(ifc);

        /* Terminar esto
         for(size_t  i = 0; i < 2; i++)
        {
            triangle<T> t1;
            t1.points[0]  =  cell_barycenter;
            t1.points[i + 1]  =  pts.at(ifc);
            t1.points[i + 2 % 2]  =  barycenter(msh, face);
        }
        */
        triangle<T> t1;
        t1.points[0]  =  cell_barycenter;
        t1.points[1]  =  pts.at(ifc);
        t1.points[2]  =  barycenter(msh, face);

        if(is_inside(t1, ep))
        {
            matrix_type tsr_at_t1 = matrix_type::Zero(2,3);
            tsr_at_t1.col(0) = tensor.cell_bar;
            tsr_at_t1.col(1) = tensor.grid_points.col(ifc);
            tsr_at_t1.col(2) = tensor.face_bars.col(ifc);

            ret = trilinear_interpolation(tsr_at_t1, t1, ep);
            is_find = true;
            break;
        }

        triangle<T> t2;
        t2.points[0]  =  cell_barycenter;
        t2.points[1]  =  barycenter(msh, face);
        t2.points[2]  =  pts.at((ifc + 1)% pts.size());

        if(is_inside(t2, ep))
        {
            matrix_type tsr_at_t2 = matrix_type::Zero(2,3);
            tsr_at_t2.col(0) = tensor.cell_bar;
            tsr_at_t2.col(1) = tensor.face_bars.col(ifc);
            tsr_at_t2.col(2) = tensor.grid_points.col((ifc + 1)% pts.size());

            ret = trilinear_interpolation(tsr_at_t2, t2, ep);
            is_find = true;
            break;
        }
    }
    if(!is_find)
        throw std::logic_error("Point not found in this cell.");

    return ret;
}


template<typename T,  typename Storage, typename CellType, typename PunctualTensor,
        typename HangingNodeInfo>
dynamic_vector<T>
eval_tensor_interpolation_tri_P2(const mesh<T,2, Storage>& msh,
                          const CellType        & cell,
                          const PunctualTensor  & tau,
                          const HangingNodeInfo & hgi,
                          const point<T,2>      & ep)
{

    typedef dynamic_vector<T>  vector_type;
    typedef dynamic_vector<T>  matrix_type;

    vector_type ret = vector_type::Zero(tau.grid_points.rows());
    auto pts = points(msh,cell);
    auto vertices_ids = msh.get_vertices_pos(cell);
    auto vertices     = msh.get_vertices(cell, pts);
    assert(vertices.size() == 3);

    std::array<point<T,2>, 3> vert;
    for(size_t i = 0; i < 3; i++)
        vert[i] = vertices[i];

    //Compute coordinates r s t
    auto triarea = [](const std::array<point<T,2>, 3>& vts) -> T
    {
        T acc{};
            for (size_t i = 1; i < 2; i++)
            {
                auto u = (vts[i] - vts[0]).to_vector();
                auto v = (vts[i+1] - vts[0]).to_vector();
                auto n = cross(u, v);
                acc += n.norm() / T(2);
            }
        return acc;
    };

    auto S = triarea(vert);

    std::array< T, 3> coords;

    for(size_t i = 0; i < 3; i++)
    {
        triangle<T> t;
        t.points[0] = vertices.at((i + 1)% 3);
        t.points[1] = vertices.at((i + 2)% 3);
        t.points[2] = ep;

        coords[i] = triarea(t.points) / S;
    }
    std::array< T, 6> phi;

    auto r = coords[0];
    auto s = coords[1];
    auto t = coords[2];

    //Evaluate P2 functions
    phi[0]  =  r * ( T(2) * r - T(1));
    phi[1]  =  s * ( T(2) * s - T(1));
    phi[2]  =  t * ( T(2) * t - T(1));
    phi[3]  =  T(4) * r * s;
    phi[4]  =  T(4) * s * t;
    phi[5]  =  T(4) * r * t;

    //Compute approximation
    auto cell_id = cell.get_id();
    auto num_pts = number_of_faces(msh, cell);

    //std::cout << "tau.grid_points.col = " << tau.grid_pts_cols() << std::endl;

    for(size_t i = 0; i < 3; i++ )
    {
        auto vertex_pos = vertices_ids.at(i);
        vector_type tau_vert = tau.grid_points.col(vertex_pos);
        ret += phi[i] * tau_vert;

        vector_type tau_bar = vector_type::Zero(tau.grid_points.rows());
        if(!hgi.faces_has_hang_nodes.at(vertex_pos))
            tau_bar = tau.face_bars.col(vertex_pos);
        else //if(hgi.faces_has_hang_nodes.at(vertex_pos))
        {
            auto hang_node_pos = hgi.hang_node_pos.at((vertex_pos + 1) % num_pts);
            assert(hang_node_pos != -1);
            tau_bar = tau.grid_points.col(hang_node_pos);
        }
        ret += phi.at(i + 3) * tau_bar;
    }

    return ret;
}





template<typename TensorType, typename T, typename Storage, typename CellType,
        typename HangingNodeInfo>
auto
compute_interpolation(const mesh<T,2,Storage>& msh,
                      const mesh<T,2,Storage>& old_msh,
                      const CellType  & cell,
                      const CellType  & ancestor_cell,
                      const TensorType& tsr,
                      const HangingNodeInfo& hgi)
{
    typedef mesh<T,2,Storage>   mesh_type;
    typedef point<T,2>          point_type;
    typedef dynamic_vector<T>   vector_type;
    typedef dynamic_matrix<T>   matrix_type;

    auto area  = measure(msh, cell);
    if(area < 1.e-10)
        throw std::invalid_argument("Area < 1.e-10. Review how tiny an element could be.");

    auto sigma = tsr.sigma;
    punctual_tensor<mesh_type>  ret(msh, cell, sigma.m_quad_degree);
    ret.Zero(2);
    auto pts_to_eval = ret.tensor_points(msh, cell);
    matrix_type mat  = matrix_type::Zero(2, pts_to_eval.size());

    size_t cont = 0;

    //std::cout << "sigmagrid_points.col = " << sigma.grid_pts_cols() << std::endl;
    for(auto ep : pts_to_eval)
    {
        assert(cont < pts_to_eval.size());
        //mat.col(cont++) = eval_tensor_interpolation(old_msh, ancestor_cell, sigma, ep);
        mat.col(cont++) = eval_tensor_interpolation_tri_P2(old_msh, ancestor_cell, sigma, hgi, ep);
    }

    ret = punctual_tensor<mesh_type>(msh, cell, sigma.m_quad_degree);
    ret.split_all(mat, 2);
    return ret;
};
template<typename TensorType, typename T, typename Storage >
std::vector<TensorType>
sigma_interpolation(const mesh<T,2,Storage>& new_msh,
                    const mesh<T,2,Storage>& old_msh,
                    const std::vector<TensorType>& tsr_vec,
                    const std::vector<size_t>& ancestors,
                    const std::vector<size_t>& cells_marks,
                    const size_t degree)
{
    std::cout << "INSIDE SIGMA INTERPOLATION" << std::endl;
    typedef mesh<T,2,Storage>   mesh_type;
    typedef dynamic_vector<T>   vector_type;
    typedef dynamic_matrix<T>   matrix_type;

    size_t num_new_cells = new_msh.cells_size();
    std::vector<TensorType> new_tsr_vec;
    new_tsr_vec  = disk::tensor_zero_vector(new_msh, degree);

    auto hang_nodes_vec = check4_special_polygons(old_msh);
    for(auto& cell : new_msh)
    {
        auto cell_id        = cell.get_id();
        auto ancestor_id    = ancestors.at(cell_id);
        auto ancestor_cell  = *std::next(old_msh.cells_begin(), ancestor_id);
        auto cell_mark      = cells_marks.at(ancestor_id);

        auto old_tsr      =  tsr_vec.at(ancestor_id);
        auto new_tsr  =  new_tsr_vec.at(cell_id);

        auto hgi = hang_nodes_vec.at(ancestor_id);
        //std::cout << " *  cell_level ("<< cell_id<<"," << ancestor_id  <<") : "<< cell_mark << std::endl;
        //std::cout << "old_tsr.grid_points.col = " << old_tsr.sigma.grid_pts_cols() << std::endl;

        if( cell_mark > 1) //Just cells to refine or its neighbors ( new hanging nodes)
        {
            auto sigma = compute_interpolation(new_msh, old_msh, cell, ancestor_cell, old_tsr, hgi);
            new_tsr_vec.at(cell_id).sigma = sigma;
        }
        else
            new_tsr_vec.at(cell_id).sigma = old_tsr.sigma;
    }
    return new_tsr_vec;
};

template<typename MeshType>
struct hanging_nodes_info
{
    typedef typename  MeshType::point_type::id_type    id_type;
    bool  has_hng_nodes;
    std::vector<id_type>    vertices_ids;
    std::vector<int>        hang_node_pos;
    std::vector<size_t>     faces_has_hang_nodes; //faces with one of its extremes being an hanging node
    hanging_nodes_info(){}

    hanging_nodes_info(const MeshType& msh, const typename MeshType::cell& cell)
    {
        has_hng_nodes = false;
        auto pts_ids  = cell.point_ids();
        auto num_pts  = pts_ids.size();
        vertices_ids  = std::vector<id_type>(num_pts);
        hang_node_pos = std::vector<int>(num_pts);
        faces_has_hang_nodes = std::vector<size_t>(num_pts);
    }

    friend std::ostream& operator<<(std::ostream& os,
                                        const hanging_nodes_info<MeshType>&  hg)
    {
        std::cout << " * has_hng_nodes : "<< hg.has_hng_nodes << std::endl;
        std::cout << " * vertices_ids  : ";
        for(auto v: hg.vertices_ids)
            std::cout << "  "<< v << std::endl;
        std::cout<< std::endl;
        std::cout << " * hang_node_pos  : ";
        for(auto v: hg.hang_node_pos)
            std::cout << "  "<< v << std::endl;
        std::cout<< std::endl;
        std::cout << " * faces_has_hang_nodes : ";
        for(auto v: hg.faces_has_hang_nodes)
            std::cout << "  "<< v << std::endl;
        std::cout<< std::endl;

        return os;
    }
};
template<typename MeshType>
std::vector<hanging_nodes_info<MeshType>>
check4_special_polygons(const MeshType& msh)
{
    auto ret = std::vector<hanging_nodes_info<MeshType>>(msh.cells_size());

    size_t i = 0;

    for(auto cl : msh)
    {
        hanging_nodes_info<MeshType>  hgi(msh, cl);
        if(msh.is_special_cell(cl) ) // &  hanging nodes)
        {
            typedef typename  MeshType::point_type::id_type    id_type;
            auto pts     = points(msh, cl);
            auto pts_ids = cl.point_ids();
            auto num_pts = pts.size();
            hgi.has_hng_nodes = true;


            auto vertices_pos = msh.get_vertices_pos(cl);
            for(size_t col = 0; col < vertices_pos.size(); col++)
            {
                auto pid_itor = std::next(pts_ids.begin(),  vertices_pos.at(col));
                //hgi.vertices_ids.at(col)  =  *pid_itor;
            }
            hgi.vertices_ids = msh.get_vertices_pos(cl);

            //std::cout << "num_pts : "<< num_pts << std::endl;
            //std::cout << "*  pids : " << std::endl;
            //for(auto&  pi:  pts_ids)
            //    std::cout << "  "<< pi;
            //std::cout << std::endl;


            // Search for hanging_nodes
            for(size_t j = 0; j < num_pts; j++)
            {
                //std::cout << "* ----------- j = "<<j<<"-----------*"<< std::endl;
                bool find1 = std::binary_search(vertices_pos.begin(),
                                                vertices_pos.end(), j);
                bool find2 = std::binary_search(vertices_pos.begin(),
                                    vertices_pos.end(), (j + 1)%num_pts);

                hgi.hang_node_pos.at(j) = -1;
                if( !find1 )
                    hgi.hang_node_pos.at(j) = j; //hang_node_id.at(j) = *std::next(pts_ids.begin(), j);


                hgi.faces_has_hang_nodes.at(j) = 0;
                if(!find2 || !find1)
                    hgi.faces_has_hang_nodes.at(j) = 1;

                //std::cout << " ** pi    ("<<find1 <<"): "<< *std::next(pts_ids.begin(),j) << std::endl;
                //std::cout << " ** nxt_pi("<<find2 <<"): "<< *std::next(pts_ids.begin(),(j+1)%num_pts) << std::endl;

                //std::cout << " *********** hgi ("<< cl.get_id()<<") ****************" << std::endl;
                //std::cout << hgi << std::endl;
            }
            //std::cout << "/* special_cell */"<< vts.size() << std::endl;
        }
        ret.at(i) = hgi;
        i++;
    }
    return ret;
}

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

    bool do_refinement;
    std::vector<size_t> cells_marks;
    std::vector<size_t> levels;

    stress_based_mesh(const std::vector<size_t>& levels_vec)
    {
        levels = levels_vec;
    }

    //MArking just below and above
    bool
    marker(const mesh_type& msh,
           const std::vector<tensors<mesh_type>>& tsr_vec,
           const T& yield)
    {
        bool do_refinement = false;
        cells_marks = std::vector<size_t>(msh.cells_size());

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto mat = tsr_vec.at(i).xi_norm;

            for(size_t m = 0; m < mat.size(); m++)
            {
                std::cout << "  "<< m;

                if(std::isnan(mat(0,m)))
                std::cout << "/* WARNING: The norm of the constrain is NaN  */" << std::endl;
                //    throw std::logic_error("The norm of the constrain is NaN");
            }
            if(std::abs(mat.minCoeff() - yield) < 1.e-10
                    || std::abs(mat.minCoeff()) < yield )
            {
                if(std::abs(mat.maxCoeff() - yield) > 1.e-10
                    && std::abs(mat.maxCoeff())     > yield )
                {
                    cells_marks.at(i) = 1;
                    do_refinement = true;
                }
            }
        }
        return do_refinement;
    }
    bool
    marker(const mesh_type& msh,
          const std::vector<tensors<mesh_type>>& tsr_vec,
          const T& yield,
          const T& percent)

    {
        bool do_refinement = false;
        cells_marks = std::vector<size_t>(msh.cells_size());

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto xi_norm = tsr_vec.at(i).xi_norm;

            dynamic_matrix<T> mat = xi_norm.all();

            for(size_t m = 0; m < mat.size(); m++)
            {
                std::cout << "  "<< m;

                if(std::isnan(mat(0,m)))
                    std::cout << "/* WARNING: The norm of the constrain is NaN  */" << std::endl;
                //    throw std::logic_error("The norm of the constrain is NaN");
            }
            //Marking with an interval
            T percent = 0.05;
            for(size_t i = 0; i < mat.size(); i++)
            {
                if(std::abs(mat(i) - (1. - percent) * yield) < 1.e-10
                        || std::abs(mat(i)) > (1. - percent) * yield )
                {
                    if(std::abs(mat(i) - (1. + percent)*yield) < 1.e-10
                        || std::abs(mat(i)) < (1. + percent) * yield )
                    {
                        cells_marks.at(i) = 1;
                        do_refinement = true;
                    }
                }
            }
        }
        return do_refinement;
    }

    template<typename LoaderType>
    auto
    refine(const mesh_type & msh,
           const std::string&  directory)
    {
        mesh_type re_msh;

        auto storage    = msh.backend_storage();
        auto storage_rm = re_msh.backend_storage();
        auto num_cells  = msh.cells_size();
        size_t nds_counter = 0;

        for(auto& cl:msh)
        {
            size_t id  = cl.get_id();
            auto cmkr  = cells_marks.at(id);
            auto lmkr  = (id != 0)?  cells_marks.at(id-1) : 0;
            auto rmkr  = (id != msh.cells_size()-1)?  cells_marks.at(id + 1) : 0;

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

        return re_msh;
    }

};

template<typename T, typename Storage>
class stress_based_mesh<mesh<T,2,Storage>>
{

public:
    typedef mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::scalar_type         scalar_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::cell                cell_type;
    typedef typename mesh_type::face                face_type;
    typedef typename mesh_type::cell_iterator       cell_itor;

    typedef dynamic_matrix<scalar_type>     matrix_type;
    typedef dynamic_vector<scalar_type>     vector_type;

    size_t imsh;
    std::string info;
    mesh_parameters<T> mp;
    mesh_type old_msh;
    plasticity_data<T> pst;
    std::vector<size_t>  cells_marks;
    std::vector<size_t>  levels;
    std::vector<std::pair<size_t,size_t>> level_ancestor;
    std::vector<size_t>  new_levels;
    std::vector<size_t>  ancestors;
    std::vector<vector_type> new_Uh;
    //std::vector<std::pair<bool,typename point<T,2>::id_type>> faces_marks;
    std::vector<std::pair<bool, int>> faces_marks;
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
    std::vector<std::array<ident_impl_t,13>>        m_triadecagons;
    std::vector<std::array<ident_impl_t,14>>        m_tesseradecagons;
    std::vector<std::array<ident_impl_t,15>>        m_pentadecagons;

    std::vector<std::pair<size_t ,size_t>>        l_triangles;
    std::vector<std::pair<size_t ,size_t>>        l_quadrangles;
    std::vector<std::pair<size_t ,size_t>>        l_pentagons;
    std::vector<std::pair<size_t ,size_t>>        l_hexagons;
    std::vector<std::pair<size_t ,size_t>>        l_heptagons;
    std::vector<std::pair<size_t ,size_t>>        l_octagons;
    std::vector<std::pair<size_t ,size_t>>        l_enneagons;
    std::vector<std::pair<size_t ,size_t>>        l_decagons;
    std::vector<std::pair<size_t ,size_t>>        l_hendecagons;
    std::vector<std::pair<size_t ,size_t>>        l_dodecagons;
    std::vector<std::pair<size_t ,size_t>>        l_triadecagons;
    std::vector<std::pair<size_t ,size_t>>        l_tesseradecagons;
    std::vector<std::pair<size_t ,size_t>>        l_pentadecagons;

    std::vector< hanging_nodes_info<mesh_type>>   hang_nodes_vec;

    stress_based_mesh(const mesh_type& msh,
                        const std::vector<size_t>& levels_vec,
                        const plasticity_data<T>&  m_pst,
                        const mesh_parameters<T>&  msh_parameters,
                        const size_t& adaptive_step):
                        pst(m_pst),
                        mp(msh_parameters),
                        imsh(adaptive_step)
    {
        info  =  mp.summary + "_RC" + tostr(imsh);

        std::cout << "Constructor sbm mp.summary :"<< mp.summary << std::endl;

        //check_older_msh(msh);
        cells_marks = std::vector<size_t>(msh.cells_size(),0);
        levels = levels_vec;
        level_ancestor.reserve(msh.cells_size());
        faces_marks.resize(msh.faces_size());

        for(size_t i = 0; i < faces_marks.size(); i++)
            faces_marks.at(i).first = false;

        // Marking all cells was conceived as a test
        if(mp.mark_all)
            cells_marks = std::vector<size_t>(msh.cells_size(),3);
        else
        {
            //This refinement is aimed to start only with triangles.
            //if we want to have more polygons take out this and solve
            //the problem with the refinement_other
            if(imsh == 0)
            {
                for(auto& cell: msh)
                {
                    auto cell_id = cell.get_id();
                    auto pts = points(msh, cell);
                    auto vts = msh.get_vertices(cell , pts);
                    if(vts.size() > 3)
                        cells_marks.at(cell_id) = 3;
                }
            }
        }

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";
        std::cout  << std::endl;

        hang_nodes_vec = check4_special_polygons(msh);

    }

    template< size_t N>
    struct n_gon
    {
       std::array<size_t, N>  p;
       std::array<bool, N>    b;

       friend std::ostream& operator<<(std::ostream& os,  const n_gon<N>& pol)
       {
           os << " Polygon points"<< std::endl;
           for(size_t i = 0; i < N; i++)
                os << pol.p[i] << "  ";
           return os;
       }
    };

    void
    check_levels(const mesh_type& msh, const cell_type & cl)
    {
        auto cl_id   = cl.get_id();
        auto fcs     = faces(msh,cl);

        //std::cout << "cell :"<< cl_id;
        //std::cout << "      level : "<< cl_level << std::endl;
        for(auto& fc : fcs)
        {
            auto cl_level= levels.at(cl_id);
            if(!(msh.is_boundary(fc)))
            {
                size_t ngh_id    = msh.neighbor_id(cl, fc);
                auto   ngh_level = levels.at(ngh_id);

                if( std::abs(int(ngh_level - cl_level)) >= 2 )
                {
                    cell_type search_cl = cl;
                    size_t    search_id = cl_id;

                    if(ngh_level < cl_level)
                    {
                        search_id = ngh_id;
                        search_cl = *std::next(msh.cells_begin(), ngh_id);
                    }

                    ++levels.at(search_id);
                    cells_marks.at(search_id) = 3;
                    std::cout << " ** mark cell : "<< search_id << std::endl;
                    check_levels(msh , search_cl);
                }
            }
        }
        return;
    }
    bool
    test_marker(const mesh_type& msh,
                const T& ratio,
                const size_t  imsh)
    {
        bool do_refinement = false;

        for(auto& cl : msh)
        {
            auto cl_id    = cl.get_id();
            auto pts      = points(msh, cl);
            auto c_faces  = faces(msh, cl);
            std::vector<point_type> bars(pts.size());

            for(size_t i = 0; i < pts.size(); i++)
            {
                auto face  = c_faces.at(i);
                bars.at(i) = barycenter(msh,face);
            }
            pts.insert(pts.end(), bars.begin(), bars.end());

            std::vector<T> r(pts.size());
            for(size_t i = 0; i < pts.size(); i++)
                r.at(i)  = std::sqrt(pts[i].x() * pts[i].x() + pts[i].y() * pts[i].y());
            auto result  = std::minmax_element(r.begin(), r.end());

            std::cout << "cell : "<< cl_id << std::endl;
            std::cout << " * pts : ";
            for(auto& p: pts)
                std::cout << p <<"  ";
            std::cout<< std::endl;
            std::cout << " * r_max: "<< *(result.second) << std::endl;
            std::cout << " * r_min: "<< *(result.first) << std::endl;

            if(*(result.second) >=  ratio )
            {
                if (*(result.first) < ratio )
                {
                    std::cout << " * cell marked" << std::endl;
                    cells_marks.at(cl_id) = 3;
                    if(imsh > 0)
                        ++levels.at(cl_id);
                    do_refinement = true;
                }
            }

        }

        dump_to_matlab(msh, mp.directory +  info + "_wl.m", cells_marks);

        for(auto& cl : msh)
           check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }

    void
    check_for_NaN(const std::vector<tensors<mesh_type>>& tsr_vec)
    {
        for(size_t i = 0; i < tsr_vec.size(); i++)
        {
            auto mat = tsr_vec.at(i).xi_norm.join_all();
            for(size_t m = 0; m < mat.size(); m++)
            {
                if(std::isnan(mat(0,m)))
                    throw std::logic_error("The norm of the constrain is NaN");
            }
        }
        return;
    }

    bool
    marker_xc(const mesh_type& msh,
            const std::vector<tensors<mesh_type>>& tsr_vec,
            const T& yield)
    {
        typedef disk::quadrature<mesh_type, cell_type>      cell_quad_type;
        typedef dynamic_vector<scalar_type>     vector_type;

        std::cout << "INSIDE_XC_MARKER" << std::endl;
        bool do_refinement(false);

        cells_marks = std::vector<size_t>(msh.cells_size());

        check_for_NaN(tsr_vec);

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto  xi_norm = tsr_vec.at(i).xi_norm;
            matrix_type mat1 = xi_norm.cell_quad_pts;
            matrix_type mat2 = xi_norm.grid_points;

            auto mat = xi_norm.all();

            bool in_interval(false);

            //marker XC
            if(std::abs(mat1.minCoeff() - yield) <= 1.e-10 || std::abs(mat1.minCoeff()) < yield )
            {
                if(std::abs(mat1.maxCoeff() - yield) > 1.e-10 && std::abs(mat1.maxCoeff()) > yield )
                    in_interval = true;
            }
            if(std::abs(mat2.minCoeff() - yield) <= 1.e-10 || std::abs(mat2.minCoeff()) < yield )
            {
                if(std::abs(mat2.maxCoeff() - yield) > 1.e-10 && std::abs(mat2.maxCoeff()) > yield )
                    in_interval = true;
            }

            #if 0
            //marker KC
            for(size_t j = 0; j < mat1.size(); j++)
            {
                if(std::abs(mat1.minCoeff() - yield) <= 1.e-10 || std::abs(mat1.minCoeff()) < yield )
                {
                    if(std::abs(mat1(j) - (1. - mp.percent) * yield) < 1.e-10
                                        || std::abs(mat1(j)) > (1. - mp.percent)*yield )
                    {
                        if(std::abs(mat1(j) - (1. + mp.percent) * yield) < 1.e-10
                                        || std::abs(mat1(j)) < (1. + mp.percent)*yield )
                            in_interval = true;
                    }
                }

                if(std::abs(mat1.maxCoeff() - yield) > 1.e-10 && std::abs(mat1.maxCoeff()) > yield )
                {
                    if(std::abs(mat1(j) - (1. - mp.percent) * yield) < 1.e-10
                                        || std::abs(mat1(j)) > (1. - mp.percent)*yield )
                    {
                        if(std::abs(mat1(j) - (1. + mp.percent) * yield) < 1.e-10
                                        || std::abs(mat1(j)) < (1. + mp.percent)*yield )
                            in_interval = true;
                    }
                }
            }

            for(size_t j = 0; j < mat2.size(); j++)
            {
                if(std::abs(mat2.minCoeff() - yield) <= 1.e-10 || std::abs(mat2.minCoeff()) < yield )
                {
                    if(std::abs(mat2(j) - (1. - mp.percent) * yield) < 1.e-10
                                        || std::abs(mat2(j)) > (1. - mp.percent)*yield )
                    {
                        if(std::abs(mat2(j) - (1. + mp.percent) * yield) < 1.e-10
                                        || std::abs(mat2(j)) < (1. + mp.percent)*yield )
                            in_interval = true;
                    }
                }

                if(std::abs(mat2.maxCoeff() - yield) > 1.e-10 && std::abs(mat2.maxCoeff()) > yield )
                {
                    if(std::abs(mat2(j) - (1. - mp.percent) * yield) < 1.e-10
                                        || std::abs(mat2(j)) > (1. - mp.percent)*yield )
                    {
                        if(std::abs(mat2(j) - (1. + mp.percent) * yield) < 1.e-10
                                        || std::abs(mat2(j)) < (1. + mp.percent)*yield )
                            in_interval = true;
                    }
                }
            }
            #endif

            if(in_interval)
            {
                cells_marks.at(i) = 3;
                ++levels.at(i);
                do_refinement = true;
            }
        }
        dump_to_matlab(msh, mp.directory + "/mesh_pr_" +  info + ".m", cells_marks);

        for(auto& cl : msh)
           check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }
    bool
    marker_jb(const mesh_type& msh,
            const std::vector<tensors<mesh_type>>& tsr_vec)
    {
        std::cout << "INSIDE_JB_MARKER" << std::endl;

        check_for_NaN(tsr_vec);

        bool do_refinement(false);
        cells_marks = std::vector<size_t>(msh.cells_size());

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto xi_norm = tsr_vec.at(i).xi_norm;
            matrix_type mat = xi_norm.join_all();

            //Marking with an interval
            bool in_interval = false;
            for(size_t j = 0; j < mat.size(); j++)
            {
                if(std::abs(mat(j) - (1. - mp.percent) * pst.yield) < 1.e-10
                                    || std::abs(mat(j)) > (1. - mp.percent)* pst.yield )
                {
                    if(std::abs(mat(j) - (1. + mp.percent) * pst.yield) < 1.e-10
                                    || std::abs(mat(j)) < (1. + mp.percent)* pst.yield )
                        in_interval = true;
                }
            }

            //marker XC
            if(std::abs(mat.minCoeff() - pst.yield) <= 1.e-10
                                        || std::abs(mat.minCoeff()) < pst.yield )
            {
                if(std::abs(mat.maxCoeff() - pst.yield) > 1.e-10
                                        && std::abs(mat.maxCoeff()) > pst.yield )
                    in_interval = true;
            }

            if(in_interval)
            {
                cells_marks.at(i) = 3;
                ++levels.at(i);
                do_refinement = true;
            }
        }
        dump_to_matlab(msh, mp.directory + "/mesh_pr_" +  info + ".m", cells_marks);
        for(auto& cl : msh)
            check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }

    bool
    marker_m4(const mesh_type& msh,
            const std::vector<tensors<mesh_type>>& tsr_vec)
    {
        std::cout << "INSIDE_M4_MARKER" << std::endl;

        check_for_NaN(tsr_vec);

        bool do_refinement(false);
        cells_marks = std::vector<size_t>(msh.cells_size());

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto xi_norm = tsr_vec.at(i).xi_norm;
            matrix_type mat = xi_norm.join_all();

            //Marking with an interval
            bool in_interval_xc = false;
            bool in_interval_jb = false;
            for(size_t j = 0; j < mat.size(); j++)
            {
                if(std::abs(mat(j) - (1. - mp.percent) * pst.yield) < 1.e-10
                            || std::abs(mat(j)) > (1. - mp.percent)*pst.yield )
                {
                    if(std::abs(mat(j) - (1. + mp.percent) * pst.yield) < 1.e-10
                                    || std::abs(mat(j)) < (1. + mp.percent)*pst.yield )
                        in_interval_jb = true;
                }
            }


            //marker XC
            if(std::abs(mat.minCoeff() -  (1. - mp.percent) * pst.yield) <= 1.e-10
                        || std::abs(mat.minCoeff()) < (1. - mp.percent)*pst.yield )
            {
                if(std::abs(mat.maxCoeff() - (1. + mp.percent) * pst.yield) > 1.e-10
                        && std::abs(mat.maxCoeff()) >(1. + mp.percent)*pst.yield)
                    in_interval_xc = true;
            }

            // State 2
            bool state2 = false;
            if(std::abs(mat.minCoeff() -  (1. - mp.percent)*pst.yield) <= 1.e-10
                        || std::abs(mat.minCoeff()) > (1. - mp.percent)*pst.yield )
            {
                if(std::abs(mat.maxCoeff() - (1. + mp.percent) * pst.yield) > 1.e-10
                        && std::abs(mat.maxCoeff()) <(1. + mp.percent)*pst.yield)
                    state2 = true;
            }

            //State 5
            bool state5 = !state2;

            //state 7 (eventualmente se convierte en el 5, asi que este funciona como y 5-7)
            bool state5_7 = false;
            if(std::abs(mat.minCoeff()) < (1. - mp.percent)*pst.yield )
            {
                if(std::abs(mat.maxCoeff()) >(1. + mp.percent)*pst.yield)
                    state5_7 = true;
            }
            //state3
            bool state3 = false;
            if(std::abs(mat.minCoeff() -  (1. - mp.percent) * pst.yield) <= 1.e-10 ||
                (std::abs(mat.minCoeff()) > (1. - mp.percent) * pst.yield
                    && std::abs(mat.minCoeff()) < (1. + mp.percent)*pst.yield  ))
            {
                if(std::abs(mat.maxCoeff()) >(1. + mp.percent)*pst.yield)
                    state3 = true;
            }
            //state4
            bool state4 = false;
            if(std::abs(mat.maxCoeff() -  (1. + mp.percent) * pst.yield) <= 1.e-10 ||
                (std::abs(mat.maxCoeff()) > (1. - mp.percent) * pst.yield
                    && std::abs(mat.maxCoeff()) < (1. + mp.percent) * pst.yield  ))
            {
                if(std::abs(mat.minCoeff()) <(1. - mp.percent) * pst.yield)
                    state3 = true;
            }

            if(state5_7 || state3 || state4)
            {
                cells_marks.at(i) = 3;
                ++levels.at(i);
                do_refinement = true;
            }
        }
        dump_to_matlab(msh, mp.directory + "/mesh_pr_" +  info + ".m", cells_marks);
        for(auto& cl : msh)
            check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }
#if 0




        // Term T3
        for (size_t fc_id = 0; fc_id < msh.faces_size(); fc_id++)
        {
            //std::cout << "FACE = "<< face_i << std::endl;
            //auto current_face_range = dsr.face_range(face_i);
            T eta_stress(0);
            std::vector<T> Cg(2);
            auto varpi =  0.5;
            auto fc   = *(msh.faces_begin() + fc_id);
            auto fqs  = face_quadrature.integrate(msh, fc);
            auto h_F  = measure(msh, fc);


            auto owner_cells = _vec.at(fc_id);
            std::cout<<"FACE_ID = "<< fc_id<<std::endl;
            for(size_t ii = 0; ii < owner_cells.size(); ii++)
            {
                auto cl_id = owner_cells.at(ii);
                auto cl    = *std::next(msh.cells_begin() , size_t(cl_id));
                auto tsr   =  tsr_vec.at(cl_id);
                auto m_T   =  measure(msh, cl);
                auto vts       =  msh.get_vertices(cl, pts);
                auto h_T       =  diameter(msh, vts);
                auto n     =  normal(msh, cl, fc);

                Cg.at(ii)  = Cp * Cp  +  varpi * (Ct/m_T) * h_T * h_T;

                auto fids_of_cl     = cl.faces_ids();
                auto num_faces      = fids_of_cl.size();
                auto cell_range     = cell_basis.range(0,m_degree);
                auto num_cell_dofs  = cell_range.size();
                auto num_face_dofs  = face_basis.size();
                dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);
                size_t cell_size = dsr.cell_range().size();

                std::cout<<"      CELL_ID :    "<< cl_id<<"   , fc_id:"<< fc_id<<std::endl;
                std::cout<<"      CELL_FCs:" << std::endl;
                for(auto& f : fids_of_cl)
                    std::cout<<"      "<< f;
                std::cout<<std::endl;

                typedef std::pair<size_t, size_t>  myPair;
                std::vector<myPair>                sorted_faces(num_faces);
                for(size_t j = 0 ; j < num_faces; j++)
                {
                    auto face = fids_of_cl.at(j);
                    sorted_faces.at(j) = std::make_pair(face , j);
                }
                std::sort(sorted_faces.begin(),sorted_faces.end());
                //typename face_type::id_type fid(fc_id);
                auto lower = std::lower_bound(sorted_faces.begin(), sorted_faces.end(),
                                    std::make_pair(fc_id, 0), [](myPair lhs, myPair rhs)
                                        -> bool { return lhs.first < rhs.first; });
                auto find   =  std::distance(sorted_faces.begin() , lower);
                auto fc_pos =  sorted_faces.at(find).second;
                //auto fc_pos = 1;
                std::cout<<"      FC_POS :    "<< fc_pos << std::endl;

                if(lower == sorted_faces.end())
                    throw std::invalid_argument("This is a bug: face not found");

                auto current_face_range = dsr.face_range(fc_pos);
                size_t cont = current_face_range.min();

                gradrec_nopre.compute(msh, cl);
                auto uh_TF = Uh_Th.at(cl_id);
                //WK: if ever the polynomial degree change between cells. The Next
                //    "for" will remain true since values in faces are patched. So
                //    u_F|T_1 = u_F|T_2. However, I could have u^1_T1|F1, u^4_T1|F2,
                //    u^2_T1|F3. Then you should pay attention storing the gauss
                //    points for each face, since then number of gpts could change.
                for (auto& qp : fqs)
                {
                    auto c_phi  = cell_basis.eval_functions(msh, cl, qp.point());
                    auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point());
                    auto f_phi  = face_basis.eval_functions(msh, fc, qp.point());
                    auto col_range  =   cell_basis.range(1,m_degree+1);
                    auto row_range  =   dof_range(0,mesh_type::dimension);

                    matrix_type c_dphi_matrix =   make_gradient_matrix(c_dphi);
                    matrix_type c_dphi_taken  =   take(c_dphi_matrix, row_range, col_range);
                    matrix_type c_dphi_rec    =   c_dphi_taken * gradrec_nopre.oper;
                    vector_type dphi_r_uh     =   c_dphi_rec * uh_TF;

                    vector_type tau  = tsr.sigma.col(cont) - pst.alpha * tsr.gamma.col(cont);
                    vector_type str  = tau + pst.alpha * dphi_r_uh;
                    auto jump = make_prod( str , n);
                    //auto proj_jump = make_prod(tau + pst.alpha*dphi_r_uh , n);
                    eta_stress  +=  qp.weight() * varpi * h_F * jump;
                    cont++;
                }
            }

            for(size_t ii = 0; ii < owner_cells.size(); ii++)
                eta.at(ii).first +=  Cg.at(ii) * (eta_stress);
        }

        //WK: try to do this with tansform
        //1st trial: std::transform (vec.begin(), vec.end(), vec.begin(), [](std::pair<double, int>& l)
        //                  {return std::make_pair(std::sqrt(l.first),l.second)});
        for(size_t cl_id = 0; cl_id < msh.cells_size(); cl_id++)
            eta.at(cl_id).first = std::sqrt(eta.at(cl_id).first);

        std::sort(eta.begin(), eta.end());

        std::cout << "ETA : " << std::endl;
        for(auto& e : eta)
            std::cout<< e.first << " "<< e.second<< std::endl;
        size_t num_adapt = int(std::ceil(percent * msh.cells_size()));
        for(size_t i = msh.cells_size() - 1; i >= msh.cells_size() - num_adapt ; i--)
        {
            auto id = eta.at(i).second;
            cells_marks.at(id) = 3;
            ++levels.at(id);
            do_refinement = true;
        }

#endif


    dynamic_matrix<T>
    gradient_lagrange_pol1(const mesh_type& msh, const cell_type& cl,
                           const std::vector<point_type> vts, const point_type& ep )
    {
        // WARNING: This only works for triangles and linear Lagrange polynomials
        dynamic_matrix<T> gradient = dynamic_matrix<T>::Zero(2,3);
        auto S   = measure(msh,cl);
        auto num_vts = vts.size();

        for(size_t i = 0; i < num_vts; i++)
        {
            size_t idx1 = (i + 1)% num_vts;
            size_t idx2 = (i + 2)% num_vts;
            point_type p1 = vts.at(idx1);
            point_type p2 = vts.at(idx2);

            gradient(0, i) = (0.5 / S) * (p1.y() - p2.y()) ;
            gradient(1, i) = (0.5 / S) * (p2.x() - p1.x()) ;
        }
        return gradient;
    }

    dynamic_matrix<T>
    gradient_lagrange_pol1(const mesh_type& msh, const cell_type& cl,
                           const triangle<T>& t, const point_type& ep )
    {
        // WARNING: This only works for triangles and linear Lagrange polynomials
        dynamic_matrix<T> gradient = dynamic_matrix<T>::Zero(2,3);
        auto S   = measure(msh,cl);

        for(size_t i = 0; i < 3; i++)
        {
            size_t idx1 = (i + 1)% 3;
            size_t idx2 = (i + 2)% 3;
            point_type p1 = t.grid_points.at(idx1);
            point_type p2 = t.grid_points.at(idx2);

            gradient(0, i) = (0.5 / S) * (p1.y() - p2.y()) ;
            gradient(1, i) = (0.5 / S) * (p2.x() - p1.x()) ;
        }
        return gradient;
    }


    template<typename TensorType>
    T
    div_tau_P1(const mesh_type& msh, const cell_type& cl,
                      const TensorType tau, const point_type& ep, const size_t& position)
    {
        auto pts = points(msh, cl);
        std::vector<point_type> vts;

        if( msh.is_special_cell(cl) ) // &  hanging nodes)
        {
             vts = msh.get_vertices(cl , pts);
             //std::cout << "/* special_cell */"<< vts.size() << std::endl;
        }
        else
            vts = pts;

        auto num_vts = vts.size();
        if(num_vts != 3)
            assert(num_vts == 3);

        auto vts_ids =  msh.get_vertices_pos(cl);
        matrix_type tau_vec =  dynamic_matrix<T>::Zero(2,3);
        size_t j = 0;

        for(auto& id : vts_ids)
        {
            tau_vec.col(j)  =  tau.col(position + id);
            j++;
        }

        auto grad_matrix = gradient_lagrange_pol1(msh, cl,vts, ep);
        auto grad_tau =  grad_matrix.cwiseProduct(tau_vec);
        auto div_tau  =  grad_tau.sum();

        return div_tau;
    }

    template< typename PunctualTensor>
    T
    div_tau_P2( const mesh_type       & msh,
                const cell_type       & cell,
                const PunctualTensor  & tau,
                const size_t m_degree)
    {
        #if 0
        Warning Quadrature points: we use the points directly from the quadrature
        rule, since we dont want to transform the points to the physical element,
        (as they are stored in the quadrature_data). The reason is that we are
        evaluating the test functions already in the reference element.
        #endif
        auto m_quad_data = triangle_quadrature(2 * m_degree + 2);
        size_t num_test_funct = 6;

        matrix_type dr_phi = matrix_type::Zero(num_test_funct, m_quad_data.size());
        matrix_type ds_phi = matrix_type::Zero(num_test_funct, m_quad_data.size());


        // dr_phi  ds_phi
        for(size_t i = 0; i < m_quad_data.size(); i++)
        {
            auto quad_point = m_quad_data.at(i).first;
            T r = quad_point.x();
            T s = quad_point.y();

            dr_phi(0,i)  =  T(4) * r - T(1);
            dr_phi(1,i)  =  T(0);
            dr_phi(2,i)  =  T(4) * r + T(4) * s - T(3);
            dr_phi(3,i)  =  s;
            dr_phi(4,i)  = -T(4) * s;
            dr_phi(5,i)  = -T(8) * r + T(4);

            ds_phi(0,i)  =  T(0);
            ds_phi(1,i)  =  T(4) * s - T(1);
            ds_phi(2,i)  =  T(4) * r + T(4) * s - T(3);
            ds_phi(3,i)  =  r;
            ds_phi(4,i)  = -T(8) * s + T(4);
            ds_phi(5,i)  = -T(4) * r;
        }

        auto pts_ids  = cell.point_ids();
        auto pts = points(msh, cell);
        auto vertices_ids = msh.get_vertices_pos(cell);
        auto vertices = msh.get_vertices(cell, pts);
        auto num_pts  = number_of_faces(msh, cell);
        matrix_type tau_i = matrix_type::Zero( 2, 6);

        auto cell_id = cell.get_id();
        auto hgi = hang_nodes_vec.at(cell_id);
        for(size_t i = 0; i < 3; i++ )
        {
            auto vertex_pos = vertices_ids.at(i);
            tau_i.col(i) = tau.grid_points.col(vertex_pos);

            if(!hgi.faces_has_hang_nodes.at(vertex_pos))
                tau_i.col(3 + i) = tau.face_bars.col(i);
            else
                tau_i.col(3 + i) = tau.grid_points.col((i + 1)% num_pts);
        }

        // dw/dr  dx/ds   dy/dr  dy/ds
        auto ddr = vertices.at(0) - vertices.at(2);
        auto dds = vertices.at(1) - vertices.at(2);

        // Integral
        vector_type d_tau = vector_type::Zero(2);

        for(size_t i = 0; i < num_test_funct; i++)
        {
            vector_type d_tau_i = vector_type::Zero(2);

            for(size_t j = 0; j < m_quad_data.size(); j ++)
            {
                auto weight = m_quad_data.at(j).second;

                vector_type d_Phi = vector_type::Zero(2);

                d_Phi(0) = ( dr_phi(i,j) * dds.y() - ds_phi(i,j) * dds.y());
                d_Phi(1) = (-dr_phi(i,j) * dds.x() + ds_phi(i,j) * ddr.x());

                d_tau_i +=  weight * tau_i.col(i).cwiseProduct(d_Phi);
            }
            d_tau += d_tau_i;
        }

        return d_tau(0) + d_tau(1);
    }

    template< typename PunctualTensor>
    dynamic_vector<T>
    grad_tau_P2( const mesh_type       & msh,
                 const cell_type       & cell,
                 const PunctualTensor  & tau,
                 const point_type      & ep)
    {
        #if 0
        Warning Quadrature points: Since here we are going to use a point ep,
        already transformed into the physical element, we have to use test functions
        in the physical element. To do so, we evaluate 'r(x,y)' and 's(x,y)'.
        Check the question about the integration in the physical element for
        grad u
        #endif
        size_t num_test_funct = 6;

        auto pts = points(msh, cell);
        auto vts = msh.get_vertices(cell, pts);
        auto num_pts = number_of_faces(msh, cell);
        // dr_phi  ds_phi
        //std::cout << "r : "<< r <<"  ; s : "<< s << std::endl;
        auto x1 = vts[0].x();   auto x2 = vts[1].x();    auto x3 =  vts[2].x();
        auto y1 = vts[0].y();   auto y2 = vts[1].y();    auto y3 =  vts[2].y();

        auto ax  = (x2 - x3) / (x1 - x3) ;
        auto cx  =  x3 /(x1 - x3);
        auto bxy = (y1 - y3) / (x1 - x3);
        auto mxy = (y2 - y3) - ax * (y1 - y3);

        auto M =   ( T(1) / (x1 - x3))  + ( bxy * ax / mxy ) ;
        auto N = - ax / mxy  ;
        auto B = - cx + ( y3 / mxy) * ax - ( x3 / mxy ) * bxy * ax;;

        auto P = - ( T(1) / mxy ) * bxy;
        auto O =   ( T(1) / mxy );
        auto C = - ( y3 / mxy) + (x3 / mxy )* bxy;

        auto U = - P - M;
        auto V = - O - N;
        auto W =   T(1) - B - C;

        matrix_type d_phi = matrix_type::Zero(2,6);

        auto s = (T(1) / mxy) * ( (ep.y() - y3) - (ep.x() - x3) * bxy );
        auto r = - ax * s + ep.x() /(x1 - x3) - cx;
        std::cout << "ep_inside : "<< ep << std::endl;
        std::cout << "r : "<< r <<"  ; s : "<< s << std::endl;
        //dxdPhi
        d_phi(0,0) = T(4) * ( M * M * ep.x() +  M * N * ep.y() +  M * B) - M;
        d_phi(0,1) = T(4) * ( P * P * ep.x() +  P * O * ep.y() +  P * C) - P;
        d_phi(0,2) = T(4) * ( U * U * ep.x() +  U * V * ep.y() +  U * W) - U;
        d_phi(0,3) = T(4) * ( T(2)* M*P*ep.x() + (M*O + N*P)*ep.y() + (M*C + P*B));
        d_phi(0,4) = T(4) * ( T(2)* U*P*ep.x() + (U*O + V*P)*ep.y() + (U*C + P*W));
        d_phi(0,5) = T(4) * ( T(2)* M*U*ep.x() + (U*N + V*M)*ep.y() + (U*B + M*W));

        d_phi(1,0) = T(4) * ( N * N * ep.y() +  M * N * ep.x() +  N * B) - N;
        d_phi(1,1) = T(4) * ( O * O * ep.y() +  P * O * ep.x() +  O * C) - O;
        d_phi(1,2) = T(4) * ( V * V * ep.y() +  U * V * ep.x() +  V * W) - V;
        d_phi(1,3) = T(4) * ( T(2)* N*O*ep.y() + (M*O + N*P)*ep.x() + (N*C + O*B));
        d_phi(1,4) = T(4) * ( T(2)* V*O*ep.y() + (U*O + V*P)*ep.x() + (V*C + O*W));
        d_phi(1,5) = T(4) * ( T(2)* V*N*ep.y() + (U*N + V*M)*ep.x() + (V*B + N*W));

        //dydPhi
        vector_type grad = vector_type::Zero(2);

        auto cell_id = cell.get_id();
        auto hgi = hang_nodes_vec.at(cell_id);
        auto vertices_pos = msh.get_vertices_pos(cell);
        for(size_t i = 0; i < 3; i++)
        {
            auto vertex_pos = vertices_pos.at(i);
            vector_type tau_i = tau.grid_points.col(vertex_pos);
            grad += d_phi.col(i).cwiseProduct(tau_i);

            vector_type tau_bar = vector_type::Zero(2);
            if(!hgi.faces_has_hang_nodes.at(vertex_pos))
                 tau_bar = tau.face_bars.col(vertex_pos);
            else //if(hgi.faces_has_hang_nodes.at(vertex_pos))
            {
                auto hang_node_pos = hgi.hang_node_pos.at((vertex_pos + 1) % num_pts);
                assert(hang_node_pos != -1);
                tau_bar =  tau.grid_points.col(hang_node_pos);
            }
            grad += d_phi.col(i + 3).cwiseProduct(tau_bar);
        }
        return grad;
    }



    template< typename PunctualTensor>
    T
    div_tau_P2( const mesh_type       & msh,
                const cell_type       & cell,
                const PunctualTensor  & tau,
                const point_type      & ep)
    {
        #if 0
        Warning Quadrature points: Since here we are going to use a point ep,
        already transform to the physical element, we have to use test functions
        in the physical element. To do so, we evaluate 'r(x,y)' and 's(x,y)'.
        Check the question about the integration in the physical element for
        grad u
        #endif
        size_t num_test_funct = 6;

        vector_type dr_phi = vector_type::Zero(num_test_funct);
        vector_type ds_phi = vector_type::Zero(num_test_funct);
        auto pts = points(msh, cell);
        auto vts = msh.get_vertices(cell, pts);
        auto num_pts = number_of_faces(msh, cell);
        // dr_phi  ds_phi
        //std::cout << " vts : "<< std::endl;
        //for(auto v :vts)
        //    std::cout << v << std::endl;
        //std::cout << " ep : "<< ep << std::endl;
        auto o = (vts[1].x() - vts[2].x()) / (vts[0].x() - vts[2].x());
        auto b = (ep.x() - vts[2].x()) / (vts[0].x() - vts[2].x());

        auto e = (vts[1].x() - vts[2].x()) / (vts[0].x() - vts[2].x()) ;
        auto d = (ep.x() - vts[2].x()) / (vts[0].x() - vts[2].x()) ;
        auto f = (vts[1].y() - vts[2].y()) - e * (vts[0].y() - vts[2].y());

        T s =  (T(1) / f) * ( (ep.y()   - vts[2].y()) - (vts[0].y() - vts[2].y()) * d);
        T r =  b  - s * o;

        //std::cout << "r : "<< r <<"  ; s : "<< s << std::endl;

        dr_phi(0)  =  T(4) * r - T(1);
        dr_phi(1)  =  T(0);
        dr_phi(2)  =  T(4) * r + T(4) * s - T(3);
        dr_phi(3)  =  s;
        dr_phi(4)  = -T(4) * s;
        dr_phi(5)  = -T(8) * r + T(4);

        ds_phi(0)  =  T(0);
        ds_phi(1)  =  T(4) * s - T(1);
        ds_phi(2)  =  T(4) * r + T(4) * s - T(3);
        ds_phi(3)  =  r;
        ds_phi(4)  = -T(8) * s + T(4);
        ds_phi(5)  = -T(4) * r;

        auto pts_ids  = cell.point_ids();
        auto vertices_pos = msh.get_vertices_pos(cell);
        auto cell_id  = cell.get_id();
        auto hgi   = hang_nodes_vec.at(cell_id);
        matrix_type tau_i = matrix_type::Zero( 2, 6);

        for(size_t i = 0; i < 3; i++ )
        {
            auto vertex_pos = vertices_pos.at(i);
            tau_i.col(i) = tau.grid_points.col(vertex_pos);

            if(!hgi.faces_has_hang_nodes.at(vertex_pos))
                tau_i.col(3 + i) = tau.face_bars.col(vertex_pos);
            else
            {
                auto col = vertices_pos.at(i) + 1;
                tau_i.col(3 + i) = tau.grid_points.col(col % num_pts);
            }
        }

        // dw/dr  dx/ds   dy/dr  dy/ds
        auto ddr = vts.at(0) - vts.at(2);
        auto dds = vts.at(1) - vts.at(2);

        // Integral
        vector_type d_tau = vector_type::Zero(2);

        for(size_t i = 0; i < num_test_funct; i++)
        {
            vector_type d_tau_i = vector_type::Zero(2);


            vector_type d_Phi = vector_type::Zero(2);

            d_Phi(0) = ( dr_phi(i) * dds.y() - ds_phi(i) * dds.y());
            d_Phi(1) = (-dr_phi(i) * dds.x() + ds_phi(i) * ddr.x());

            d_tau_i  +=  tau_i.col(i).cwiseProduct(d_Phi);
        }

        return 0;
    }

    template< typename PunctualTensor>
    T
    div_tau(const mesh_type   & msh,
                              const cell_type       & cell,
                              const PunctualTensor  & tensor,
                              const point_type      & ep)
    {
        typedef dynamic_matrix<T>  matrix_type;
        T div_tau = 0.;
        auto is_find = false;
        auto pts   = points( msh, cell);
        auto fcs   = faces(  msh, cell);
        auto cell_barycenter = barycenter(msh, cell);

        #if 0
        //With barycenter
        for(size_t ifc = 0; ifc < pts.size(); ifc++)
        {
            auto face  = fcs.at(ifc);

            /* Terminar esto
             for(size_t  i = 0; i < 2; i++)
            {
                triangle<T> t1;
                t1.points[0]  =  cell_barycenter;
                t1.points[i + 1]  =  pts.at(ifc);
                t1.points[i + 2 % 2]  =  barycenter(msh, face);
            }
            */


            triangle<T> t1;
            t1.points[0]  =  cell_barycenter;
            t1.points[1]  =  pts.at(ifc);
            t1.points[2]  =  barycenter(msh, face);

            if(is_inside(t1, ep))
            {
                matrix_type tsr_at_t1 = matrix_type::Zero(2,3);
                tsr_at_t1.col(0) = tensor.cell_bar;
                tsr_at_t1.col(1) = tensor.grid_points.col(ifc);
                tsr_at_t1.col(2) = tensor.face_bars.col(ifc);

                matrix_type grad_matrix = gradient_lagrange_pol1(msh, cell, t1, ep);
                matrix_type grad_tau =  grad_matrix.cwiseProduct(tsr_at_t1);

                div_tau =  grad_tau.sum();

                return div_tau;
            }

            triangle<T> t2;
            t2.points[0]  =  cell_barycenter;
            t2.points[1]  =  barycenter(msh, face);
            t2.points[2]  =  pts.at((ifc + 1)% pts.size());

            if(is_inside(t2, ep))
            {
                matrix_type tsr_at_t2 = matrix_type::Zero(2,3);
                tsr_at_t2.col(0) = tensor.cell_bar;
                tsr_at_t2.col(1) = tensor.face_bars.col(ifc);
                tsr_at_t2.col(2) = tensor.grid_points.col((ifc + 1)% pts.size());

                auto grad_matrix = gradient_lagrange_pol1(msh, cell, t2, ep);
                auto grad_tau =  grad_matrix.cwiseProduct(tsr_at_t2);
                div_tau  =  grad_tau.sum();

                return div_tau;
            }

        }
                throw std::logic_error("Point not found in this cell.");
        #endif

        //Lineal
        //#if 0
        auto vts      =  msh.get_vertices(cell, pts);
        auto vts_ids  =  msh.get_vertices_pos(cell);
        auto num_vts = vts.size();
        if(num_vts != 3)
            assert(num_vts == 3);

        triangle<T> t3;
        t3.points[0]  =  vts.at(0);
        t3.points[1]  =  vts.at(1);
        t3.points[2]  =  vts.at(2);

        if(is_inside(t3, ep))
        {
            matrix_type tsr_at_t3 = matrix_type::Zero(2,3);
            tsr_at_t3.col(0) = tensor.grid_points.col(vts_ids.at(0));
            tsr_at_t3.col(1) = tensor.grid_points.col(vts_ids.at(1));
            tsr_at_t3.col(2) = tensor.grid_points.col(vts_ids.at(2));

            auto grad_matrix = gradient_lagrange_pol1(msh, cell, t3, ep);
            auto grad_tau =  grad_matrix.cwiseProduct(tsr_at_t3);
            div_tau  =  grad_tau.sum();

            return div_tau;
        }

        std::cout << "point :" << ep << std::endl;
        std::cout << "vts :" << std::endl;
        for (auto& v: vts)
        {
            std::cout << "* p : "<< v << std::endl;
        }

        std::cout << "pts :" << std::endl;
        for (auto& v: pts)
        {
            std::cout << "* p : "<< v << std::endl;
        }

        throw std::logic_error("Point not found in this cell.");
        //#endif

        #if 0
        // 4 triangles
        for(size_t ifc = 0; ifc < pts.size(); ifc++)
        {
            auto face  = fcs.at(ifc);

            /* Terminar esto
             for(size_t  i = 0; i < 2; i++)
            {
                triangle<T> t1;
                t1.points[0]  =  cell_barycenter;
                t1.points[i + 1]  =  pts.at(ifc);
                t1.points[i + 2 % 2]  =  barycenter(msh, face);
            }
            */


            triangle<T> t1;
            std::cout << "face left : " << (ifc + pts.size() -1) % pts.size() << std::endl;
            std::cout << "face right: " << ifc << std::endl;
            auto ifc_left   = (ifc + pts.size() -1) % pts.size();
            auto ifc_right  = ifc;
            auto face_left  = fcs.at( ifc_left);
            auto face_right = fcs.at( ifc);

            t1.points[0]  =  barycenter(msh, face_left);
            t1.points[1]  =  pts.at(ifc);
            t1.points[2]  =  barycenter(msh,face_right);

            if(is_inside(t1, ep))
            {
                matrix_type tsr_at_t1 = matrix_type::Zero(2,3);
                tsr_at_t1.col(0) = tensor.face_bars.col(ifc_left);
                tsr_at_t1.col(1) = tensor.grid_points.col(ifc);
                tsr_at_t1.col(2) = tensor.face_bars.col(ifc);

                matrix_type grad_matrix = gradient_lagrange_pol1(msh, cell, t1, ep);
                matrix_type grad_tau =  grad_matrix.cwiseProduct(tsr_at_t1);

                div_tau =  grad_tau.sum();
                is_find = true;
                break;
            }
        }

        triangle<T> t2;
        t2.points[0]  =  cell_barycenter;
        t2.points[1]  =  barycenter(msh, face);
        t2.points[2]  =  pts.at((ifc + 1)% pts.size());

        if(is_inside(t2, ep))
        {
            matrix_type tsr_at_t2 = matrix_type::Zero(2,3);
            tsr_at_t2.col(0) = tensor.cell_bar;
            tsr_at_t2.col(1) = tensor.face_bars.col(ifc);
            tsr_at_t2.col(2) = tensor.grid_points.col((ifc + 1)% pts.size());

            auto grad_matrix = gradient_lagrange_pol1(msh, cell, t2, ep);
            auto grad_tau =  grad_matrix.cwiseProduct(tsr_at_t2);
            div_tau  =  grad_tau.sum();

            is_find = true;
            break;
        }

        throw std::logic_error("Point not found in this cell.");
        #endif

    }


    void
    error_to_matlab(const mesh_type& msh,
                    const std::vector<std::pair<T,size_t>>&  eta)
    {
        std::ofstream    mefs(mp.directory + "/estimator_colors.m");

        size_t cont = 0;
        size_t num_intervals = 8;
        size_t num_cells_per_interval = std::ceil(msh.cells_size() / num_intervals);

        std::cout << "numero de intervalos = "<< num_cells_per_interval << std::endl;
        std::cout << "numero de cells per interval "<< num_cells_per_interval << std::endl;

        mefs << "hold on;"<<std::endl;
        mefs << "color_mat = [ "<<std::endl;
        mefs << "    1. 0 0;      \%  red   "<<std::endl;
        mefs << "    0.5 0 0.5;   \%  Magenta"<<std::endl;
        mefs << "    0 0 0.5;     \%  Dark blue"<<std::endl;
        mefs << "    0 0 1;       \%  blue"<<std::endl;
        mefs << "    0 0.8 0.8;   \%  cyan"<<std::endl;
        mefs << "    0 0.4 0;     \%  Dark green"<<std::endl;
        mefs << "    0 1 0;       \%  green"<<std::endl;
        mefs << "    1 1 1;       \%  white"<<std::endl;

        mefs << "    0  0 0 ; "<<std::endl;
        mefs << "    0.2  0.2 0.2; "<<std::endl;
        mefs << "    0.5  0.5 0.5; "<<std::endl;
        mefs << "    1 0.5 1; "<<std::endl;
        mefs << "    1 0.2 1; "<<std::endl;
        mefs << "    0.5 1 1 ;"<<std::endl;
        mefs << "    0.1 1 1] ;"<<std::endl;

    #if 0
        color_mat = [
    1. 0 0;       %   red
    0.5 0 0.5;    %   Magenta
    1 0.8 1       %   pink;
    1 0.2 1       %   rose;
    0 0 0.5;      %   Dark blue
    0 0 1;        %   blue
    0 0.8 0.8;    %   cyan
    0 0.4 0;      %   Dark green
    0.2  0.6 0.4  %   middle green;
    0 1 0;        %   green
    0  0 0        %   black;
    0.5  0.5 0.5  %   gray;
    1 1 1;        %   white
    0.5 0.5 0.2   %   dark ocre
    1 1 0.2       %   yellow;
    0.7 0.7 0.2   %   ocre
    ] ;
    #endif

        size_t fact = 1;
        if ((msh.cells_size() % num_intervals) == 0)
        {
            fact = 0;
            std::cout << "exact division" << std::endl;
        }

        for(size_t i = 0; i < num_intervals + fact; i++)
        {
            size_t num_cells =  num_cells_per_interval;
            if( i == num_intervals - 1 )
                num_cells = msh.cells_size() - i * num_cells_per_interval;

            for(size_t j = 0; j < num_cells; j++ )
            {
                auto cl  = *std::next( msh.cells_begin() ,  eta.at(cont).second);
                auto pts =  points(msh, cl);
                auto b  = barycenter(msh,cl);
                //std::cout << cont <<" "<<  eta.at(cont).second << "  "<< eta.at(cont).first<< std::endl;

                mefs << " coords = [";
                for(auto & p : pts)
                    mefs <<eta.at(cont).second <<" "<< p.x() <<" " << p.y()<<std::endl;
                mefs << " ];"<<std::endl;
                mefs << "fill(coords(:,2), coords(:,3), color_mat("<<i + 1<<",:)); "<<std::endl;
                mefs<< "strName = strtrim(cellstr(num2str("<< eta.at(cont).second<<",'(%d)')));"<<std::endl;
                mefs<< "\%text("<<b.x()<< ","<< b.y() <<",strName,'VerticalAlignment','bottom');"<<std::endl;

                ++cont;

            }
        }
        mefs.close();
        assert(cont == msh.cells_size());
        return;
    }
    void
    error_to_matlab(const mesh_type& msh,
                    const std::vector<std::pair<T,size_t>>&  eta,
                    const std::string& name)
    {
        std::ofstream    mefs(mp.directory + "/estimator_" + name + ".m");
        std::cout << "INSIDE error_to_matlab" << std::endl;
        #if 0
        size_t cont = 0;
        size_t num_intervals = 8;
        size_t num_cells_per_interval = std::ceil(msh.cells_size() / num_intervals);

        std::cout << "numero de intervalos = "<< num_cells_per_interval << std::endl;
        std::cout << "numero de cells per interval "<< num_cells_per_interval << std::endl;

        mefs << "hold on;"<<std::endl;
        mefs << "color_mat = [ "<<std::endl;
        mefs << "    1. 0 0;      \%  red   "<<std::endl;
        mefs << "    0.5 0 0.5;   \%  Magenta"<<std::endl;
        mefs << "    0 0 0.5;     \%  Dark blue"<<std::endl;
        mefs << "    0 0 1;       \%  blue"<<std::endl;
        mefs << "    0 0.8 0.8;   \%  cyan"<<std::endl;
        mefs << "    0 0.4 0;     \%  Dark green"<<std::endl;
        mefs << "    0 1 0;       \%  green"<<std::endl;
        mefs << "    1 1 1;       \%  white"<<std::endl;

        mefs << "    0  0 0 ; "<<std::endl;
        mefs << "    0.2  0.2 0.2; "<<std::endl;
        mefs << "    0.5  0.5 0.5; "<<std::endl;
        mefs << "    1 0.5 1; "<<std::endl;
        mefs << "    1 0.2 1; "<<std::endl;
        mefs << "    0.5 1 1 ;"<<std::endl;
        mefs << "    0.1 1 1] ;"<<std::endl;

        size_t fact = 1;
        if ((msh.cells_size() % num_intervals) == 0)
        {
            fact = 0;
            std::cout << "exact division" << std::endl;
        }

        #endif

        size_t cont = 1;
        T value = eta.at(0).first;
        std::vector<int> cells_per_color;
        for(size_t i = 1; i < eta.size() ; i++)
        {
            auto e = eta.at(i);
            T new_value = e.first;
            if(std::abs((value - new_value)/value)  < 1.e-5)
            {
                cont++;
            }
            else
            {
                cells_per_color.push_back(cont);
                cont = 1;
            }
            value = new_value;
        }
        cells_per_color.push_back(cont);

        mefs << "c = jet("<< cells_per_color.size() <<"); "<< std::endl;
        mefs << "color_mat = flipud(c);"<<std::endl;
        mefs << "hold on;"<< std::endl;
        cont = 0;
        for(size_t i = 0; i < cells_per_color.size(); i++)
        {
            size_t num_cells =  cells_per_color.at(i);
            for(size_t j = 0; j < num_cells; j++ )
            {
                auto cl  = *std::next( msh.cells_begin(),  eta.at(cont).second);
                auto pts =  points(msh, cl);
                auto b   =  barycenter(msh,cl);
                //std::cout << cont <<" "<<  eta.at(cont).second << "  "<< eta.at(cont).first<< std::endl;

                mefs << " coords = [";
                for(auto & p : pts)
                    mefs <<eta.at(cont).second <<" "<< p.x() <<" " << p.y()<<std::endl;
                mefs << " ];"<<std::endl;
                mefs << "hc = fill(coords(:,2), coords(:,3), color_mat("<<i + 1<<",:)) ;"<<std::endl;
                mefs << "set(hc,'EdgeColor','none');"<<std::endl;
                mefs<< "strName = strtrim(cellstr(num2str("<< eta.at(cont).second<<",'(%d)')));"<<std::endl;
                mefs<< "\%text("<<b.x()<< ","<< b.y() <<",strName,'VerticalAlignment','bottom');"<<std::endl;

                mefs<< "strTau = strtrim(cellstr(num2str("<< eta.at(cont).first<<",'(%d)')));"<<std::endl;
                mefs<< "text("<<b.x()<< ","<< b.y() <<",strTau,'VerticalAlignment','bottom');"<<std::endl;

                mefs<<  eta.at(cont).first<<";" <<std::endl;
                ++cont;
            }
        }
        auto vmax = eta.at(0).first;
        auto vmin = eta.at(eta.size()-1).first;
        auto h = (vmax - vmin) / 10.;

        mefs<< "set(gca,'XTick',get(gca,'YTick'));"<<std::endl;
        mefs<< "colormap(fliplr(color_mat));"<<std::endl;;
        mefs<< "hcb = colorbar;"<<std::endl;;
        mefs<< "caxis([3.4577e-14 8.73979e-05]);"<<std::endl;;
        mefs<< "caxis(["<< vmin <<" "<< vmax <<"]);"<<std::endl;
        mefs<< "set(hcb,'YTick',["<<vmin<<":"<< h <<":"<< vmax <<"]);"<<std::endl;

        mefs.close();
        assert(cont == msh.cells_size());
        return;
    }

//#if 0



    template<typename CellBasisType, typename CellQuadType,
             typename FaceBasisType, typename FaceQuadType,
             typename Solution>
    dynamic_vector<T>
    proj_residue(const mesh_type& msh,
                 const cell_type& cell,
                 const dynamic_vector<T>& c_ruh,
                 const punctual_tensor<mesh_type>& c_tau,
                 const Solution& solution,
                 const size_t m_degree)
    {
        std::cout << "INSIDE PROJ_RESIDUE" << std::endl;

        typedef CellBasisType           cell_basis_type;
        typedef FaceBasisType           face_basis_type;
        typedef CellQuadType            cell_quadrature_type;
        typedef FaceQuadType            face_quadrature_type;
        typedef typename mesh_type::point_type               point_type;
        typedef T  scalar_type;
        typedef dynamic_vector<scalar_type>     vector_type;
        typedef dynamic_matrix<scalar_type>     matrix_type;

        cell_basis_type  cb(m_degree + 1);
        auto col_range = cb.range(1,m_degree+1);
        projector_nopre<mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type> projk(m_degree);

        auto cell_id  = cell.get_id();
        auto residue =  [&](const point_type& p, const size_t number) -> scalar_type
        {
            vector_type ddphi     =  cb.eval_laplacians(msh, cell, p);
            vector_type ddphi_zm  =  take(ddphi, col_range);
            T lap_ruh  = ddphi_zm.dot(c_ruh);

            vector_type  grad_tau = grad_tau_P2( msh, cell, c_tau, p);
            auto div_tau = grad_tau.sum();

            auto f = solution.f(p);

            return f + div_tau + lap_ruh;
        };

        vector_type uT = projk.compute_cell(msh, cell, residue);

        return uT;
    }


    template<typename CellBasisType, typename FaceBasisType, typename Solution>
    auto
    error_estimate(T set_error,
            const mesh_type& msh,
            const std::vector<tensors<mesh_type>>& tsr_vec,
            const std::vector<dynamic_vector<T>>&  Uh_Th,
            const size_t& m_degree,
            const Solution& solution,
            const std::vector<matrix_type>& grad_global)

    {
        typedef CellBasisType                   cell_basis_type;
        typedef FaceBasisType                   face_basis_type;

        typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
        typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

        typedef std::vector<matrix_type>                     grad_reconst_type;
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  tensor_matrix;

        //WK: This is ok  whenever k is equal in all cells and then the number of pts gauss are the same.
        //otherwise this should be inside for(auto& cl_id : owner_cells)
        auto quad_degree = tsr_vec[0].sigma.m_quad_degree;

        cell_basis_type                         cell_basis(m_degree + 1);
        face_basis_type                         face_basis(m_degree);
        face_quadrature_type                    face_quadrature(quad_degree);
        cell_quadrature_type                    cell_quadrature(quad_degree);

        std::vector<T> eta_stress(msh.cells_size());
        std::vector<T> eta_residual(msh.cells_size());

        //parameters
        auto Cp = 1. / M_PI;
        auto Ct = 0.77708;

        auto col_range  =   cell_basis.range(1,m_degree+1);
        auto row_range  =   dof_range(0,mesh_type::dimension);

        for(auto& cell : msh)
        {
            auto cell_id   =  cell.get_id();
            auto c_faces   =  faces(msh, cell);
            auto c_tsr     =  tsr_vec.at(cell_id);
            auto meas_T    =  measure( msh, cell);
            auto pts       =  points(msh,cell);
            auto vts       =  msh.get_vertices(cell, pts);
            auto h_T       =  diameter(msh, vts);
            auto num_faces =  c_faces.size();

            vector_type c_uh_TF   =  Uh_Th.at(cell_id);
            matrix_type rec_oper  =  grad_global.at(cell_id);
            vector_type c_ruh     =  rec_oper * c_uh_TF;

            auto c_gamma   =  tsr_vec.at(cell_id).gamma;
            auto c_sigma   =  tsr_vec.at(cell_id).sigma;

            size_t cont_fc  = 0;
            size_t fqs_size = 0;

            //term 4
            auto h_F     =  h_T;

            for( auto& face : c_faces)
            {
                T eta_stress_loc = 0.;

                //********** KWarning:  Review this for hanging nodes *********!!!!
                //auto h_F      =  measure(msh, face);
                //********** KWarning:  Review this for hanging nodes *********!!!!

                auto varpi    =  (msh.is_boundary(face))? 0.: 0.5 ;

                size_t ngh_id;

                if( !msh.is_boundary(face) )
                    ngh_id = msh.neighbor_id(cell, face);
                else
                    ngh_id = cell_id;

                auto neighbor  = *std::next(msh.cells_begin(),  ngh_id);

                auto c_normal  =  normal(msh, cell, face);

                auto c_nFT  =  normal_factor(msh, cell, face);
                auto n_nFT  =  normal_factor(msh, neighbor, face);

                if(c_nFT == n_nFT && !msh.is_boundary(face))
                    throw std::logic_error("Cells sharing face have the same normal. Review normal factor");


                //std::cout << "FACE : "<< fcs_ids.at(borrar++) << std::endl;
                //std::cout << " * cell ("<< cell_id<<") : "<< c_fqs_size << std::endl;
                //std::cout << " * neigh("<< neighbor_id<<") : "<< n_fqs_size << std::endl;


                auto fqs = face_quadrature.integrate(msh, face);
                fqs_size += fqs.size();

                auto n_num_faces = number_of_faces(msh, neighbor);

                matrix_type n_rec_oper =  grad_global.at(ngh_id);
                vector_type nc_uh_TF   =  Uh_Th.at(ngh_id);
                vector_type n_ruh      =  n_rec_oper * nc_uh_TF;

                auto n_gamma   =  tsr_vec.at(ngh_id).gamma;
                auto n_sigma   =  tsr_vec.at(ngh_id).sigma;

                auto c_fqs_size = c_sigma.face_qps_cols();
                auto n_fqs_size = n_sigma.face_qps_cols();

                assert ((c_fqs_size/num_faces) == fqs.size());
                assert ((n_fqs_size/n_num_faces) ==  fqs.size());

                matrix_type c_tau  = c_sigma.face_quad_pts - pst.alpha * c_gamma.face_quad_pts;
                matrix_type n_tau  = n_sigma.face_quad_pts - pst.alpha * n_gamma.face_quad_pts;

                auto c_face_pos  = face_position(msh, cell, face);
                auto n_face_pos  = face_position(msh, neighbor, face);

                auto cq_range = c_sigma.qsr.face_range(c_face_pos);
                auto nq_range = n_sigma.qsr.face_range(n_face_pos);

                auto cqp_pos = cq_range.min();
                auto nqp_pos = nq_range.min();

                //WK: This should change if k different for each face

                for (auto& qp : fqs)
                {
                    auto c_dphi  = cell_basis.eval_gradients(msh, cell, qp.point());
                    //Esto debe ser igual que lo de arriba lo que va a cambiar el gradiente es u_TF en cada elemento
                    auto nc_dphi = cell_basis.eval_gradients(msh, neighbor, qp.point());

                    matrix_type c_dphi_matrix  =   make_gradient_matrix(c_dphi);
                    matrix_type c_dphi_taken   =   take(c_dphi_matrix, row_range, col_range);
                    vector_type c_dphi_ruh     =   c_dphi_taken * c_ruh;

                    matrix_type nc_dphi_matrix =   make_gradient_matrix(nc_dphi);
                    matrix_type nc_dphi_taken  =   take(nc_dphi_matrix, row_range, col_range);
                    vector_type nc_dphi_ruh    =   nc_dphi_taken * n_ruh;

                    //estimator only with plasticity
                    vector_type c_str  = c_tau.col(cqp_pos) + pst.alpha *  c_dphi_ruh;
                    vector_type n_str  = n_tau.col(nqp_pos) + pst.alpha * nc_dphi_ruh;

                    vector_type jump = (c_str * c_nFT)  + (n_str * n_nFT);
                    auto jp_value    = varpi * make_prod( jump , c_normal);
                    eta_stress_loc  += qp.weight() * jp_value * jp_value;
                    cqp_pos++;
                    nqp_pos++;
                }

                auto eT_stress = h_T * std::sqrt( (Ct/meas_T) * h_F * eta_stress_loc);
                eta_stress.at(cell_id) +=  eT_stress;
            }

            T eta_res = 0.;
            T eta_res1(0), eta_res2(0);
            //term 1
            auto cqs = cell_quadrature.integrate(msh, cell);
            auto c_tau =  c_sigma - pst.alpha * c_gamma;
            //auto div_tau  = div_tau_P2( msh, cell, c_tau, m_degree);
            //std::cout << "cell : " << cell_id << std::endl;
            for(auto& qp : cqs)
            {
                vector_type ddphi     =  cell_basis.eval_laplacians(msh, cell, qp.point());
                vector_type ddphi_zm  =  take(ddphi, col_range);

                T    lap_ruh  =  ddphi_zm.dot(c_ruh);

                //matrix_type c_tau_mat = c_tau.join_all();
                //auto div_tau_1  =  div_tau_P1(msh, cell, c_tau_mat, qp.point(), pos);

                //std::cout << "*  div_tau_old : "<< div_tau_1 << std::endl;
                //auto div_tau  = div_tau(msh, cell, c_tau, qp.point());

                //std::cout << "*  div_tau_new : "<< div_tau << std::endl;
                vector_type grad_tau = grad_tau_P2( msh, cell, c_tau, qp.point());
                auto div_tau = grad_tau.sum();
                //auto div_tau_prueba =  div_tau_P2(msh, cell, c_tau, qp.point());

                auto f   = solution.f(qp.point());
                // cuidado aqui porque div_tau calculado con div_tau_P2 entrega la Integral
                //entonces hay qe quitarlo de la integral en el loop de cqs. y se sumaria a eta_res
                //por fuera del loop


                vector_type proj_res_dof = proj_residue<cell_basis_type, cell_quadrature_type,
                                        face_basis_type, face_quadrature_type>
                                        (msh, cell, c_ruh, c_tau, solution, 0);
                auto phi = cell_basis.eval_functions(msh, cell, qp.point());

                T proj_res_fun = proj_res_dof(0) * phi.at(0);

                std::cout << "proj_res_fun ("<< cell_id <<"): "<< proj_res_fun << std::endl;

                eta_res += qp.weight() * (f + div_tau +  pst.alpha * lap_ruh - proj_res_fun)
                            * (f + div_tau +  pst.alpha * lap_ruh - proj_res_fun);

                auto e1 = qp.weight() * (f + div_tau +  pst.alpha * lap_ruh)// - proj_res_fun)
                            * (f + div_tau +  pst.alpha * lap_ruh);// - proj_res_fun);

                auto e2 = qp.weight() * (f + div_tau +  pst.alpha * lap_ruh - proj_res_fun)
                            * (f + div_tau +  pst.alpha * lap_ruh - proj_res_fun);

                eta_res1 += e1;
                eta_res2 += e2;

                std::cout << "e1 = "<< e1 <<"    ; e2 = "<< e2 << std::endl;
                std::cout << "eta_res1 : "<< eta_res1 << "      ; eta_res2 : "<< eta_res2 <<std::endl;
                    #if 0
                    T    lap_ruh  =  ddphi_zm.dot(c_ruh);
                    matrix_type c_tau    =  c_tsr.sigma - pst.alpha * c_tsr.gamma;;
                    auto f   = solution.f(qp.point());
                    auto div_tau  =  divergence_tau_P0(msh, cell, c_tau, qp.point(), pos);
                    //eta_res += qp.weight() * iexp_pow( f, 2.);
                    //eta_res += qp.weight() * iexp_pow( div_tau, 2.);
                    eta_res += iexp_pow( div_tau, 2.);
                    //eta_res += qp.weight() * iexp_pow( lap_ruh, 2.);


                    std::cout << " * c_tau_size : "<< c_tau.rows()<< " x  "<<c_tau.cols()<< std::endl;
                    for(size_t j = 0; j < c_tau.rows(); j++)
                    {
                        for(size_t i = 0; i < c_tau.cols(); i++)
                        {
                            std::cout<< "  " <<c_tau(j , i) ;
                        }
                        std::cout << std::endl;
                    }

                    std::cout << " * div_tau : "<< div_tau << std::endl;
                    std::cout << " * isnan ? : "<< std::isnan(div_tau) << std::endl;
                    #endif

            }

            //std::cout << "      * f_cell : "<< eta_f <<"      lap_u : "<< eta_lap << std::endl;

            eta_residual.at(cell_id) =  Cp  * h_T * std::sqrt(eta_res);
        }

        T eglobal_stress = 0.;
        for(auto& e : eta_stress)
            eglobal_stress += e;
        //ofs << eglobal_stress << "  ";

        T eglobal_res = 0.;
        for(auto& e : eta_residual)
            eglobal_res += e;

        //ofs << eglobal_res << "  ";

        std::cout << "eglobal_res : "<<eglobal_res<<"    ; eglobal_stress : "<< eglobal_stress<< std::endl;


        //WK: this should be vector of vectors if u es a vectorial function,
        //since each tensor is a matrix. Then the jump = [tau + alpha G(u)].n
        // will be a vector and not a scalar.
        #if 0
        for(auto& e: eta)
        {
            std::cout<< e.second<< "  " <<eta_residual.at(e.second) << "  "<<
                                eta_stress.at(e.second) <<"  "<< e.first << std::endl;
        }
        #endif

        set_error = T(0);
        std::vector<std::pair<T,size_t>> eta(msh.cells_size());
        std::vector<std::pair<T,size_t>> etat(msh.cells_size());
        std::vector<std::pair<T,size_t>> etas(msh.cells_size());
        std::vector<std::pair<T,size_t>> etar(msh.cells_size());

        // n_Total = (n_residual + n_jump)_T^2
        for(size_t cl_id = 0; cl_id < msh.cells_size(); cl_id++)
        {
            //Estimator
            auto eT  = iexp_pow(eta_residual.at(cl_id) + eta_stress.at(cl_id), 2.);
            auto eT1 = eta_stress.at(cl_id);
            auto eT2 = eta_residual.at(cl_id);
            set_error += eT;

            eta.at(cl_id).first  = eT;
            eta.at(cl_id).second = cl_id;

            etat.at(cl_id).first  = std::sqrt(eT);
            etat.at(cl_id).second = cl_id;

            etas.at(cl_id).first  = eT1;
            etas.at(cl_id).second = cl_id;

            etar.at(cl_id).first  = eT2;
            etar.at(cl_id).second = cl_id;
        }

        auto comp = [](const std::pair<T, size_t>& a, const std::pair<T, size_t>& b) -> bool {
            return a.first > b.first;
        };

        std::sort(eta.begin(),  eta.end(), comp);
        std::sort(etat.begin(), etat.end(), comp);
        std::sort(etas.begin(), etas.end(), comp);
        std::sort(etar.begin(), etar.end(), comp);

        error_to_matlab(msh, eta,  "tot2" + info);
        error_to_matlab(msh, etat,  "tot" + info);
        error_to_matlab(msh, etas,  "str" + info);
        error_to_matlab(msh, etar,  "res" + info);

        std::ofstream       ofs(mp.directory + "/estimator.txt", std::ios::app);
        if (!ofs.is_open())
            std::cout << "Error opening ofs"<<std::endl;
        ofs<< "No. Adaptation" << imsh <<std::endl;

        disk::projector_nopre<mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type> projk(m_degree);

        for(auto& e : eta)
        {
            auto cell_id   = e.second;
            auto cell  =  *std::next(msh.cells_begin(),cell_id);

            // L2_error
            matrix_type rec_oper =  grad_global.at(cell_id);
            vector_type uTF   =  Uh_Th.at(cell_id);
            vector_type x     =  rec_oper * uTF;

            dynamic_vector<scalar_type> true_dof = projk.compute_cell(msh, cell, solution.sf);
            dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
            dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
            auto err_dof = std::sqrt(diff_dof.dot(projk.cell_mm * diff_dof));

            //L2_error gradient
            T error_dfun(0);
            typedef typename disk::solution<T,2>::gradient_vector   gradient_vector;
            auto cqs = cell_quadrature.integrate(msh, cell);

            for (auto& qp : cqs)
            {
                size_t number = 0;
                #if 0
                //Este cero es solo a causa del test con el prolema difusivo dividido
                // en dos partes. Para tuyau y los otros puede ser cualquier numero.
                // Tener en cuneta esto si se va a llamar el estimador con el proble diff

                if(mp.diff)
                    number = set_cell_number(msh, cell);
                #endif
                gradient_vector dphi_fun  =  solution.df(qp.point(),number);

                auto c_dphi  = cell_basis.eval_gradients(msh, cell, qp.point());
                matrix_type  c_dphi_matrix =  make_gradient_matrix(c_dphi);
                matrix_type  c_dphi_taken  =  take(c_dphi_matrix, row_range, col_range);
                vector_type  c_dphi_ruh    =  c_dphi_taken * x;

                error_dfun  += qp.weight() * (dphi_fun - c_dphi_ruh).dot(dphi_fun - c_dphi_ruh);
            }

            error_dfun = std::sqrt(error_dfun);

            ofs << cell_id <<"  "<< err_dof <<"  " <<error_dfun<<"  ";
            ofs << eta_residual.at(cell_id) <<"  ";
            ofs << eta_stress.at(cell_id)   <<"  "<< std::sqrt(e.first) << std::endl;
        }
        ofs <<std::endl;
        ofs.close();

        return eta;
    }

    template<typename CellBasisType, typename FaceBasisType, typename Solution>
    bool
    marker_ae( const mesh_type& msh,
            const std::vector<tensors<mesh_type>>& tsr_vec,
            const std::vector<dynamic_vector<T>>&  Uh_Th,
            const size_t& m_degree,
            const Solution& solution,
            const std::vector<matrix_type>& grad_global)
    {

        std::cout << "INSIDE_AE_MARKER" << std::endl;

        check_for_NaN(tsr_vec);

        auto set_error = T(0);
        auto eta = error_estimate<CellBasisType, FaceBasisType, Solution>
                (set_error, msh, tsr_vec, Uh_Th, m_degree, solution, grad_global);

        struct group
        {
            typedef typename std::vector<std::pair<T,size_t>>::iterator  itor_type;

            T value;
            itor_type itor_begin;
            itor_type itor_end;

            group(){};

            group(const T val, const itor_type& begin, const itor_type& end)
            {
                value = val;
                itor_begin = begin;
                itor_end   = end;
            }
        };

        auto new_set_error = T(0);
        auto add_group     = false;
        auto do_refinement = false;

        group initial_group(eta[0].first, eta.begin(), eta.begin());
        std::vector<group> eta_group  {initial_group};

        for(auto itor = eta.begin(); itor <= eta.end(); itor++)
        {
            auto e = *itor;
            add_group = true;

            for(auto& eg : eta_group)
            {
                if((std::abs(e.first - eg.value)/eg.value <= 1.e-4))
                {
                    add_group  = false;
                    eg.itor_end  = itor;
                    break;
                }
            }
            if(add_group)
            {
                group new_group(e.first ,itor, itor);
                eta_group.push_back(new_group);
            }
        }

        std::cout << "eta_group size: "<< eta_group.size() << std::endl;

        #if 0
        // Marking by cells
        for(auto& e : eta)
        {
            auto cell_id   = e.second;
            new_set_error += e.first;

            cells_marks.at(cell_id) = 3;
            std::cout << "5.*** cell_id : "<< cell_id << std::endl;
            ++levels.at(cell_id);
            do_refinement  = true;

            std::cout << "cell : "<< cell_id <<std::endl;
            std::cout << " * eta   : "<< e.first << std::endl;
            std::cout << " * error : "<< new_set_error << std::endl;
            if( new_set_error >= mp.percent * set_error )
                break;
        }
        #endif

        cells_marks = std::vector<size_t>(msh.cells_size());

        // Marking by groups
        for(auto eg : eta_group)
        {
            for(auto it = eg.itor_begin; it <= eg.itor_end; it++)
            {
                auto e = *it;
                auto cell_id   = e.second;
                new_set_error += e.first;

                cells_marks.at(cell_id) = 3;
                ++levels.at(cell_id);
                do_refinement  = true;
            }
            if( new_set_error >= mp.percent * set_error )
                break;
        }

        dump_to_matlab( msh, mp.directory + "/mesh_pr_" +  info + ".m", cells_marks);

        if(mp.marker_name == 6 )
        {
            for(size_t i = 0; i < tsr_vec.size();i++)
            {

                dynamic_matrix<T> sigmas = tsr_vec.at(i).sigma.join_all();
                dynamic_vector<T> mat(sigmas.cols());

                for(size_t m = 0; m < sigmas.cols(); m++)
                {
                    dynamic_vector<T> s = sigmas.col(m);
                    mat(m) = s.norm();
                }

                for(size_t m = 0; m < mat.size(); m++)
                {
                    if(std::isnan(mat(0,m)))
                        throw std::logic_error("The norm of the constrain is NaN");
                }

                bool state5_7 = false;
                if(std::abs(mat.minCoeff()) < (1. - 0.05)*pst.yield )
                {
                    if(std::abs(mat.maxCoeff()) >(1. + 0.05)*pst.yield)
                        state5_7 = true;
                }

                if(state5_7)// || state3 || state)
                {
                    cells_marks.at(i) = 3;
                    ++levels.at(i);
                    do_refinement = true;
                }
            }
            dump_to_matlab( msh, mp.directory + "/mesh_agp_" +  info + ".m", cells_marks);
        }

        for(auto& cl : msh)
            check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }
    template<typename CellBasisType, typename FaceBasisType, typename Solution>
    bool
    marker( const mesh_type& msh,
            const std::vector<tensors<mesh_type>>& tsr_vec,
            const plasticity_data<T>& pst,
            const std::vector<dynamic_vector<T>>&  Uh_Th,
            const size_t& m_degree,
            const Solution& solution,
            const std::vector<matrix_type>& grad_global)
    {
        bool ret;
        std::cout << "marker_name : " << mp.marker_name << std::endl;

        switch(mp.marker_name)
        {
            case 2:
                return ret = marker_jb(msh, tsr_vec);
            case 4:
                return ret = marker_m4(msh, tsr_vec);
            case 3:
            case 6:
                return ret = marker_ae<CellBasisType, FaceBasisType, Solution>
                        (msh, tsr_vec, Uh_Th, m_degree, solution, grad_global);
            #if 0
            case 5:
            if(mp.diff)
            {

                // "AE_DIFFUSION"
                do_refinement = sbm.template marker_diffusion<cell_basis_type, cell_quadrature_type,
                                                     face_basis_type, face_quadrature_type,
                                                     gradrec_type, Solution>
                (msh, m_pst, dir_name, info, mp, Uh_Th, m_degree, imsh, solution, grad_global);
                break;
            }
            #endif
            default:
                throw std::invalid_argument("No marker specification.");

        }

    }
//#endif
#if 0

    template<typename CellBasisType, typename CellQuadType,
            typename FaceBasisType, typename FaceQuadType,
             typename GradRecType,  typename Solution,
             typename MarkerParameters>
    bool
    marker(const mesh_type& msh,
            const plasticity_data<T>& pst,
            const std::string& mp.directory,
            const std::string& info,
            const MarkerParameters& mp,
            const std::vector<dynamic_vector<T>>&  Uh_Th,
            const size_t& m_degree,
            const size_t& imsh,
            const Solution& solution,
            const std::vector<matrix_type>& grad_global)
    {
        typedef CellBasisType                   cell_basis_type;
        typedef FaceBasisType                   face_basis_type;
        typedef CellQuadType                    cell_quadrature_type;
        typedef FaceQuadType                    face_quadrature_type;
        typedef GradRecType                     grad_reconst_type;
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  tensor_matrix;

        //WK: This is ok  whenever k is equal in all cells and then the number of pts gauss are the same.
        //otherwise this should be inside the for(auto& cl_id : owner_cells)
        //auto quad_degree = tsr_vec[0].quad_degree;
        auto quad_degree = 2* m_degree + 2;
        cell_basis_type                         cell_basis(m_degree + 1);
        face_basis_type                         face_basis(m_degree);
        face_quadrature_type                    face_quadrature(quad_degree);
        cell_quadrature_type                    cell_quadrature(quad_degree);

        std::cout << "INSIDE_AE_MARKER_DIFFUSION" << std::endl;

        std::ofstream       ofs(mp.directory + "/estimator.txt", std::ios::app);
        if (!ofs.is_open())
            std::cout << "Error opening ofs"<<std::endl;

        bool do_refinement(false);

        cells_marks = std::vector<size_t>(msh.cells_size());

            //Marking with estimators

        //WK: this should be vector of vectors if u es a vectorial function,
        //since each tensor is a matrix. Then the jump = [tau + alpha G(u)].n
        // will be a vector and not a scalar.
        std::vector<T> eta_stress(msh.cells_size());
        std::vector<T> eta_residual(msh.cells_size());

        //parameters
        auto Cp = 1. / M_PI;
        auto Ct = 0.77708;
        auto x0 = 0.5;
        auto col_range  =   cell_basis.range(1,m_degree+1);
        auto row_range  =   dof_range(0,mesh_type::dimension);

        for(auto& cell : msh)
        {
            auto cell_id   =  cell.get_id();
            auto c_faces   =  faces(msh, cell);
            auto meas_T    =  measure( msh, cell);
            auto vts       =  msh.get_vertices(cell, pts);
            auto h_T       =  diameter(msh, vts);
            auto num_faces =  c_faces.size();

            auto pts = points(msh,cell);
            size_t number1 = 0;
            for(auto& p : pts)
            {
                if(p.x() < x0)
                    number1 = 1;
                if(p.x() > x0)
                    number1 = 2;
            }
            if(number1 == 0)
                throw std::invalid_argument("Invalid number domain.");

            vector_type c_uh_TF   =  Uh_Th.at(cell_id);
            matrix_type rec_oper  =  grad_global.at(cell_id);
            vector_type c_ruh     =  rec_oper * c_uh_TF;

            size_t cont_fc  = 0;
            size_t fqs_size = 0;
            //term 4
            for( auto& face : c_faces)
            {
                T eta_stress_loc = 0.;
                auto h_F      =  measure(msh, face);
                auto varpi    =  (msh.is_boundary(face))? 0.: 0.5 ;

                size_t ngh_id;
                if( !msh.is_boundary(face) )
                    ngh_id = msh.neighbor_id(cell, face);
                else
                    ngh_id = cell_id;

                auto neighbor  = *std::next(msh.cells_begin(),  ngh_id);
                auto c_normal  =  normal(msh, cell, face);
                auto n_normal  =  normal(msh, neighbor, face);
                auto cgp_id    =  gauss_points_positions<CellQuadType,FaceQuadType>
                                                    (msh,cell,face, quad_degree);
                auto ngp_id    =  gauss_points_positions<CellQuadType,FaceQuadType>
                                                (msh,neighbor,face, quad_degree);

                auto pts = points(msh, neighbor);
                size_t number2 = 0;
                for(auto& p : pts)
                {
                    if(p.x() < x0)
                    number2 = 1;
                    if(p.x() > x0)
                    number2 = 2;
                }
                if(number2 == 0)
                    throw std::invalid_argument("Invalid number domain.");


                matrix_type n_rec_oper =  grad_global.at(ngh_id);
                vector_type nc_uh_TF   =  Uh_Th.at(ngh_id);
                vector_type n_ruh      =  n_rec_oper * nc_uh_TF;

                //WK: This should change if k different for each face
                auto fqs = face_quadrature.integrate(msh, face);
                fqs_size += fqs.size();

                for (auto& qp : fqs)
                {
                    auto c_phi   = cell_basis.eval_functions(msh, cell, qp.point());
                    auto c_dphi  = cell_basis.eval_gradients(msh, cell, qp.point());
                    //Esto debe ser igual que lo de arriba lo que va a cambiar el gradiente es u_TF en cada elemento
                    auto nc_phi  = cell_basis.eval_functions(msh, neighbor, qp.point());
                    auto nc_dphi = cell_basis.eval_gradients(msh, neighbor, qp.point());

                    matrix_type c_dphi_matrix  =   make_gradient_matrix(c_dphi);
                    matrix_type c_dphi_taken   =   take(c_dphi_matrix, row_range, col_range);
                    vector_type c_dphi_ruh     =   c_dphi_taken * c_ruh;

                    matrix_type nc_dphi_matrix =   make_gradient_matrix(nc_dphi);
                    matrix_type nc_dphi_taken  =   take(nc_dphi_matrix, row_range, col_range);
                    vector_type nc_dphi_ruh    =   nc_dphi_taken * n_ruh;

                    //estimator only with plasticity
                    vector_type c_str  = solution.st(qp.point(),number1) +   c_dphi_ruh;
                    vector_type n_str  = solution.st(qp.point(),number2) +  nc_dphi_ruh;
                    vector_type jump   = c_str - n_str;
                    auto jp_value      = varpi * make_prod( jump , c_normal);
                    //auto proj_jump = make_prod(tau + pst.alpha*dphi_r_uh , n);
                    eta_stress_loc  += qp.weight() * jp_value * jp_value;
                    cgp_id++;
                    ngp_id++;
                }

                auto eT_stress = h_T * std::sqrt( (Ct/meas_T) * h_F * eta_stress_loc);
                eta_stress.at(cell_id) +=  eT_stress;
            }

            T eta_res = 0.;

            //term 1
            auto cqs = cell_quadrature.integrate(msh, cell);
            auto pos = cqs.size() + fqs_size;

            for (auto& qp : cqs)
            {
                    vector_type ddphi     =  cell_basis.eval_laplacians(msh, cell, qp.point());
                    vector_type ddphi_zm  =  take(ddphi, col_range);

                    T    lap_ruh  =  ddphi_zm.dot(c_ruh);
                    auto div_tau  =  solution.dst(qp.point());
                    //auto f   = solution.f(qp.point());
                    std::cout << " * qp  = " << qp.point() << std::endl;
                    std::cout << " * div_tau  = "<< div_tau << std::endl;
                    std::cout << " * lap = "<< lap_ruh << std::endl;
                    eta_res += qp.weight() * iexp_pow(div_tau + lap_ruh, 2.);
                    //eta_res += qp.weight() * iexp_pow(lap_ruh, 2.);

            }
            eta_residual.at(cell_id) =  Cp  * h_T * std::sqrt(eta_res);
        }

        T eglobal_stress = 0.;
        for(auto& e : eta_stress)
            eglobal_stress += e;
        //ofs << eglobal_stress << "  ";

        T eglobal_res = 0.;
        for(auto& e : eta_residual)
            eglobal_res += e;

        //ofs << eglobal_res << "  ";

        std::cout << "eglobal_res : "<<eglobal_res<<"    ; eglobal_stress : "<< eglobal_stress<< std::endl;

        T set_error(0);

        disk::projector_nopre<mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type> projk(m_degree);

        std::vector<std::pair<T,size_t>> eta(msh.cells_size());
        std::vector<std::pair<T,size_t>> etas(msh.cells_size());
        std::vector<std::pair<T,size_t>> etar(msh.cells_size());
        // (n_res + n_stress)_T^2
        for(size_t cl_id = 0; cl_id < msh.cells_size(); cl_id++)
        {
            //Estimator
            auto eT  = iexp_pow(eta_residual.at(cl_id) + eta_stress.at(cl_id), 2.);
            auto eT1 = iexp_pow(eta_stress.at(cl_id), 2.);
            auto eT2 = iexp_pow(eta_residual.at(cl_id), 2.);
            set_error += eT;

            eta.at(cl_id).first  = eT;
            eta.at(cl_id).second = cl_id;

            etas.at(cl_id).first  = eT1;
            etas.at(cl_id).second = cl_id;

            etar.at(cl_id).first  = eT2;
            etar.at(cl_id).second = cl_id;
            //std::cout << cl_id<<"  "<< eta.at(cl_id).first <<"  "<< eta.at(cl_id).second << "  "<< eT2 <<std::endl;
        }
        //ofs << set_error << std::endl;

        auto comp = [](const std::pair<T, size_t>& a, const std::pair<T, size_t>& b) -> bool {
            return a.first > b.first;
        };

        std::sort(eta.begin(), eta.end(), comp);
        std::sort(etas.begin(), etas.end(), comp);
        std::sort(etar.begin(), etar.end(), comp);
        #if 0
        for(auto& e: eta)
        {
            std::cout<< e.second<< "  " <<eta_residual.at(e.second) << "  "<<
                                eta_stress.at(e.second) <<"  "<< e.first << std::endl;
        }
        #endif

        ofs<< "No. Adaptation" << imsh <<std::endl;
        for(auto& e : eta)
        {
            auto cell_id   = e.second;
            auto cell  =  *std::next(msh.cells_begin(),cell_id);
            auto pts   =  points(msh , cell);

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


            // L2_error
            matrix_type rec_oper =  grad_global.at(cell_id);
            vector_type uTF   =  Uh_Th.at(cell_id);
            vector_type x     =  rec_oper * uTF;

            dynamic_vector<scalar_type> true_dof = projk.compute_cell(msh, cell, solution.sf);
            dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
            dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
            auto err_dof = std::sqrt(diff_dof.dot(projk.cell_mm * diff_dof));

            //L2_error gradient
            T error_dfun(0);
            typedef typename disk::solution<T,2>::gradient_vector   gradient_vector;
            auto cqs = cell_quadrature.integrate(msh, cell);

            for (auto& qp : cqs)
            {
                gradient_vector dphi_fun  =  solution.df(qp.point(), number);

                auto c_dphi  = cell_basis.eval_gradients(msh, cell, qp.point());
                matrix_type  c_dphi_matrix =  make_gradient_matrix(c_dphi);
                matrix_type  c_dphi_taken  =  take(c_dphi_matrix, row_range, col_range);
                vector_type  c_dphi_ruh    =  c_dphi_taken * x;

                error_dfun  += qp.weight() * (dphi_fun - c_dphi_ruh).dot(dphi_fun - c_dphi_ruh);
            }

            error_dfun = std::sqrt(error_dfun);

            ofs << cell_id <<"  "<< err_dof <<"  " <<error_dfun<<"  ";
            ofs << eta_residual.at(cell_id) <<"  ";
            ofs << eta_stress.at(cell_id)   <<"  "<< std::sqrt(e.first) << std::endl;

            #if 0
            std::cout << cell_id <<"  "<< err_dof <<"  " <<error_dfun<<"  ";
            std::cout << eta_residual.at(cell_id) <<"  ";
            std::cout << eta_stress.at(cell_id) <<"  "<< std::sqrt(e.first) << std::endl;
            #endif 0
        }
        ofs <<std::endl;
        ofs.close();

        T new_set_error(0);
        for(auto& e : eta)
        {
            auto cell_id = e.second;
            new_set_error += e.first;

            cells_marks.at(cell_id) = 3;
            ++levels.at(cell_id);
            do_refinement  = true;
            if( new_set_error >= mp.percent * set_error )
                break;
        }
        dump_to_matlab( msh, mp.directory + "/mesh_pr_" +  info + ".m", cells_marks);
        error_to_matlab(msh, eta ,  "tot" + info);
        error_to_matlab(msh, etas,  "str" + info);
        error_to_matlab(msh, etar,  "res" + info);


        for(auto& cl : msh)
            check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }
    #endif


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
        auto cmv_temp = cells_marks;

        auto storage = msh.backend_storage();
        for(auto& cl : msh)
        {
            auto cl_id   = cl.get_id();
            auto cl_mark = cells_marks.at(cl_id);

            if(cl_mark == 1)
            {
                auto fcs = faces(msh,cl);
                for(auto& fc:fcs)
                {
                    if(!msh.is_boundary(fc))
                    {
                        auto  ngh_id   = msh.neighbor_id(cl, fc);
                        cmv_temp.at(ngh_id) = 1;
                        auto  neighbor      = *std::next(msh.cells_begin(),size_t(ngh_id));
                        auto  neighbor_fcs  = faces(msh,neighbor);
                        for(auto& nfc : neighbor_fcs)
                        {
                            auto  nfc_id    = msh.lookup(nfc);

                            if(!faces_marks.at(nfc_id).first)
                            {
                                faces_marks.at(nfc_id).first  = true;
                                auto bar = barycenter(msh,nfc);
                                storage ->points.push_back(bar);
                                faces_marks.at(nfc_id).second = storage->points.size() -1;
                                //WK: this should be with an id_type instead of using ints

                                if(!msh.is_boundary(nfc))
                                {
                                    auto  nn_id   = msh.neighbor_id( cl, nfc);
                                    cmv_temp.at(nn_id) = 1;
                                }
                            }
                        }
                    }
                    else
                    {
                        auto  fc_id    = msh.lookup(fc);
                        //std::cout << "inside marker, marking bundary = "<<fc_id << std::endl;
                        faces_marks.at(fc_id).first  = true;
                        auto bar = barycenter(msh,fc);
                        storage ->points.push_back(bar);
                        faces_marks.at(fc_id).second = storage->points.size() -1;
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
        faces_marks.at(fc_id).first  = true;
        auto bar = barycenter(msh,fc);
        storage ->points.push_back(bar);
        faces_marks.at(fc_id).second = storage->points.size() -1;
        std::cout << "     faces_mark.at("<<fc_id <<") = "<< faces_marks.at(fc_id).second << std::endl;
        //std::cout << "inside marker, marking bundary = "<<fc_id << std::endl;
    }

    void
    set_marks_hanging_nodes_4tri(mesh_type & msh, const cell_type & cl)
    {
        auto cl_id   = cl.get_id();
        auto cl_mark = cells_marks.at(cl_id);
        auto pts     = points(msh, cl);
        //To adapt also the neighbors take into account 3 and 2, otherwise just 3
        if(cl_mark == 3)
        //if(cl_mark > 1) /
        {
            bool special_tri = false;
            std::vector<typename point_type::id_type> vertices;
            if(msh.is_special_cell(cl))
            {
                vertices = msh.get_vertices_pos(cl);
                if(vertices.size() == 3)
                    special_tri = true;
            }

            if(special_tri)
            {
                auto fcs = faces(msh, cl);
                auto it0 = pts.begin() + vertices.at(0);
                auto it1 = pts.begin() + vertices.at(1);
                auto it2 = pts.begin() + vertices.at(2);
                auto it3 = pts.end() - 1;

                std::vector<face_type> mark_faces;
                if( it0  + 1 == it1)
                    mark_faces.push_back(fcs.at(vertices.at(0)));
                if( it1  + 1 == it2)
                    mark_faces.push_back(fcs.at(vertices.at(1)));
                if( it2 == it3)
                    mark_faces.push_back(fcs.at(vertices.at(2)));

                // la diferencia con esta parte es que no es recursiva: Asi que debe llamarse de nuevo la funcion en los if anteriores
                for(auto& fc : mark_faces)
                {
                    if(msh.is_boundary(fc))
                       face_mark_hanging_nodes(msh, fc);
                    else
                    {
                        auto  fc_id    = msh.lookup(fc);
                        if(!faces_marks.at(fc_id).first)
                        {

                            face_mark_hanging_nodes(msh, fc);

                            //std::cout << "cell = "<< cl_id <<"    ; fc = "<< fc_id << std::endl;
                            size_t ngh_id = msh.neighbor_id(cl, fc);

                            typename cell_type::id_type id(ngh_id);
                            auto  ngh_mark      = cells_marks.at(id);
                            if(cl_mark > ngh_mark)
                            {
                                cells_marks.at(ngh_id) = cl_mark - 1;
                                auto  ngh  = *std::next(msh.cells_begin(),size_t(ngh_id));
                                set_marks_hanging_nodes(msh, ngh);
                            }

                        }
                    }

                }
            }
            else
            {
                //Use this only if you want to divide  resulting triangles in 4 :
                // This is useful when you want to divide the polygons using barycenter
                // and the face (use also in refine_other_4tri).
                //  Here you are marking the face to add the barycenter for new triangles.

                auto fcs     = faces(msh,cl);
                for(auto & fc: fcs)
                {
                    if(msh.is_boundary(fc))
                        face_mark_hanging_nodes(msh, fc);
                    else
                    {
                        auto  fc_id    = msh.lookup(fc);
                        if(!faces_marks.at(fc_id).first)
                        {
                            face_mark_hanging_nodes(msh, fc);

                            size_t ngh_id = msh.neighbor_id(cl, fc);
                            typename cell_type::id_type id(ngh_id);
                            auto  ngh_mark      = cells_marks.at(id);
                            if(cl_mark > ngh_mark)
                            {
                                cells_marks.at(ngh_id) = cl_mark - 1;
                                auto  ngh  = *std::next(msh.cells_begin(), size_t(ngh_id));
                                set_marks_hanging_nodes(msh, ngh); // shoudn't it be written set_marks_hanging_nodes_4tri ?
                            }

                        }
                    }
                }
            }
        }
        return;

    }
    void
    set_marks_hanging_nodes(mesh_type & msh, const cell_type & cl)
    {
        //std::cout << "INSIDE SET MARKS" << std::endl;

        auto cl_id   = cl.get_id();
        auto cl_mark = cells_marks.at(cl_id);
        auto pts     = points(msh, cl);
        //To adapt also the neighbors take into account 3 and 2, otherwise just 3

        //std::cout << " - CELL: "<< cl_id << std::endl;
        if(cl_mark == 3)
        //if(cl_mark > 1) /
        {
            auto hgi = hang_nodes_vec.at(cl_id);
            bool special_tri = false;

            if( hgi.has_hng_nodes)
            {
                //std::cout << "--inside-cl_mark_hg " << std::endl;

                auto fcs = faces(msh, cl);
                //std::cout << "num_faces: "<< fcs.size() << std::endl;
                size_t fcont = 0;
                // la diferencia con esta parte es que no es recursiva: Asi que debe llamarse de nuevo la funcion en los if anteriores
                for(auto fc : fcs)
                {
                    auto fcs_pts = points(msh, fc);
                    #if 0
                    std::cout << "---fc : "<<fcont << std::endl;
                    std::cout << "---fc_points : "<< fcs_pts[0]<<" "<< fcs_pts[1]<< std::endl;
                    std::cout << "---fc_is_bnd : "<< msh.is_boundary(fc) << std::endl;
                    std::cout << "---fc_has_hg : "<< hgi.faces_has_hang_nodes.at(fcont) << std::endl;
                    #endif
                    if(hgi.faces_has_hang_nodes.at(fcont) == 0)
                    {
                        //std::cout << "----fc without hg" << std::endl;
                        if(msh.is_boundary(fc))
                            face_mark_hanging_nodes(msh, fc);
                        else
                        {
                            auto  fc_id    = msh.lookup(fc);
                            if(!faces_marks.at(fc_id).first)
                            {
                                face_mark_hanging_nodes(msh, fc);

                                //std::cout << "cell = "<< cl_id <<"    ; fc = "<< fc_id << std::endl;
                                size_t ngh_id = msh.neighbor_id(cl, fc);

                                typename cell_type::id_type id(ngh_id);
                                auto  ngh_mark      = cells_marks.at(id);

                                if(cl_mark > ngh_mark)
                                {
                                    cells_marks.at(ngh_id) = cl_mark - 1;
                                    auto  ngh  = *std::next(msh.cells_begin(),size_t(ngh_id));
                                    set_marks_hanging_nodes(msh, ngh);
                                }
                            }
                        }
                    }
                    fcont++;
                }
            }
            if(pts.size() == 3)
            {
                //Use this only if you want to divide triangles in 4 :
                // Here you are marking the face to add the barycenter for new triangles.

                auto fcs     = faces(msh,cl);
                for(auto & fc: fcs)
                {
                    if(msh.is_boundary(fc))
                        face_mark_hanging_nodes(msh, fc);
                    else
                    {
                        auto  fc_id    = msh.lookup(fc);
                        if(!faces_marks.at(fc_id).first)
                        {
                            face_mark_hanging_nodes(msh, fc);

                            size_t ngh_id = msh.neighbor_id(cl, fc);
                            typename cell_type::id_type id(ngh_id);
                            auto  ngh_mark      = cells_marks.at(id);
                            if(cl_mark > ngh_mark)
                            {
                                cells_marks.at(ngh_id) = cl_mark - 1;
                                auto  ngh  = *std::next(msh.cells_begin(), size_t(ngh_id));
                                set_marks_hanging_nodes(msh, ngh);
                            }

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
        auto cl_id   = cl.get_id();
        auto cl_mark = cells_marks.at(cl_id);

        if(cl_mark > 1  &  cl_mark < 4)
        {
            auto fcs = faces(msh,cl);
            for(auto& fc:fcs)
            {
                if(!msh.is_boundary(fc))
                {
                    auto  ngh_id = msh.neighbor_id(cl, fc);
                    auto  ngh_mark    = cells_marks.at(ngh_id);
                    auto  fc_id    = msh.lookup(fc);
                    if(!faces_marks.at(fc_id).first)
                    {
                        faces_marks.at(fc_id).first  = true;
                        auto bar = barycenter(msh,fc);
                        storage ->points.push_back(bar);
                        faces_marks.at(fc_id).second = storage->points.size() -1;
                    }
                    if(ngh_mark < cl_mark)
                    {
                        if(ngh_mark == 1)
                            cells_marks.at(ngh_id) = 2;
                        else
                            cells_marks.at(ngh_id) = cl_mark - 1;
                        auto  ngh = *std::next(msh.cells_begin(), size_t(ngh_id));
                        set_marks(msh,ngh);
                    }

                        #if 0
                        auto  ngh      = *std:next(msh.cells_begin(), size_t(ngh_id));
                        auto  ngh_fcs  = faces(msh,ngh);

                        for(auto& nfc : ngh_fcs)
                        {
                            auto  nfc_id  = msh.lookup(nfc);

                            if(!faces_marks.at(nfc_id).first)
                            {
                                faces_marks.at(nfc_id).first  = true;
                                auto bar = barycenter(msh,nfc);
                                storage ->points.push_back(bar);
                                faces_marks.at(nfc_id).second = storage->points.size() ;
                            }
                        }
                        #endif

                }
                else
                {
                    auto  fc_id    = msh.lookup(fc);
                    faces_marks.at(fc_id).first  = true;
                    auto bar = barycenter(msh,fc);
                    storage ->points.push_back(bar);
                    faces_marks.at(fc_id).second = storage->points.size() -1;
                }
            }
        }

        //else
        //    throw std::logic_error("shouldn't have arrived here");
    }

    void
    face_marker(mesh_type& msh)
    {
        for(auto& cl : msh)
        {
            auto cl_id   = cl.get_id();
            bool cl_mark = (cells_marks.at(cl_id) == 3)? true : false;
            if(cl_mark)
            {
                if(mp.hanging_nodes)
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
          auto fc_mark = faces_marks.at(fcs_ids[i]).first;
          nt.p[j] = t.p[i];
          nt.b[j] = t.b[i];
          if(fc_mark)
           {
               j++;
               nt.p[j] = faces_marks.at(fcs_ids[i]).second;
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
   bool
   new_fc_is_boundary(const mesh_type & msh,
                  const std::vector<ident_impl_t>& pts_ids,
                  const ident_impl_t original_face_id,
                  const size_t& N, const size_t& i)
   {
        auto fo = *std::next(msh.faces_begin(), original_face_id);
        auto fo_pids = fo.point_ids();
        std::sort(fo_pids.begin(), fo_pids.end());

        auto pid1 = pts_ids.at(i);
        auto pid2 = pts_ids.at((i+1)%N);

        bool find1  = std::binary_search(fo_pids.begin(), fo_pids.end(), pid1);
        bool find2  = std::binary_search(fo_pids.begin(), fo_pids.end(), pid2);

        if( find1 || find2 )
            return true;
        return false; //how this could happen?
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
           t.p[i]   = pts_ids.at(i); // Since later in populate_mesh values are subtract  by 1
           t.b[i]   = false;

           auto fid = fcs_ids.at(i);
           auto is_bound = storage->boundary_edges.at(fid);
           if(is_bound)
            t.b[i]   = new_fc_is_boundary(msh, pts_ids, fid, N, i);
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
        if(num_fcs == 13)
        {
            auto t = put_n_gon<13>(msh, cl);
            return  n_gon_base<N, 13>(t, fcs_ids);
        }
        if(num_fcs == 14)
        {
            auto t = put_n_gon<14>(msh, cl);
            return  n_gon_base<N, 14>(t, fcs_ids);
        }
        if(num_fcs == 15)
        {
            auto t = put_n_gon<15>(msh, cl);
            return  n_gon_base<N, 15>(t, fcs_ids);
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
        auto cl_id   = cl.get_id();
        auto cl_mark = cells_marks.at(cl_id);

        for(size_t i = 0; i < num_fcs; i++)
        {
            auto id  = fcs_ids.at(i);
            auto fc_mark = faces_marks.at(id).first;
            t.p[i]   = pts_ids.at(i);
            t.b[i]   = storage->boundary_edges.at(id);
        }

        if(cl_mark > 1 )
        {
            n_gon<3> t0,t1,t2,t3;
            auto bar0pos = faces_marks.at(fcs_ids[0]).second;
            auto bar1pos = faces_marks.at(fcs_ids[1]).second;
            auto bar2pos = faces_marks.at(fcs_ids[2]).second;

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
                auto fc_mark = faces_marks.at(id).first;

                if(fc_mark)
                {
                    auto barpos = faces_marks.at(fcs_ids[i]).second;
                    t0.p[i] = barpos;           t1.p[i] = barpos;
                    t0.b[permut[i]] = false;    t1.b[i] = false;
                }
            }
            store(t0,m_triangles);
            store(t1,m_triangles);
        }
    }


    void
    refine_reference_triangle(const n_gon<3> & t , const std::array<int, 3> & bar_ids,
                             const size_t& cl_level, const size_t& cl_id)
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

        auto pair = std::make_pair(cl_level, cl_id);
        l_triangles.insert(l_triangles.end(), 4, pair);
    }
    void refine_triangle(mesh_type& msh, const cell_type& cl, const size_t& cl_level)
    {
        auto fcs_ids = cl.faces_ids();
        auto bar0pos = faces_marks.at(fcs_ids[0]).second;
        auto bar1pos = faces_marks.at(fcs_ids[1]).second;
        auto bar2pos = faces_marks.at(fcs_ids[2]).second;
        std::array<int, 3> bar_ids = {bar0pos, bar1pos, bar2pos};

        auto t  = put_n_gon<3>(msh, cl);
        refine_reference_triangle( t,  bar_ids, cl_level, cl.get_id());
    }

    void
    refine_other_4tri(mesh_type& msh, const cell_type& cl, const size_t& cl_level)
    {
        // Use this function if you want to divide the resulting triangles in 4
        // Also use the function set_marks_hanging_nodes_4tri instead

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

            tbar_ids.at(0) = faces_marks.at(fcs_ids.at(i)).second;
            tbar_ids.at(1) = fbar_ids.at((i+1)%num_fcs);
            tbar_ids.at(2) = fbar_ids.at(i);

            refine_reference_triangle(nt , tbar_ids, cl_level, cl.get_id());
        }
    }
    void
    refine_other(mesh_type& msh, const cell_type& cl, const size_t& cl_level)
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
            store(nt,m_triangles);
            l_triangles.push_back(std::make_pair(cl_level, cl.get_id()));
        }
    }

    template<typename IdxVector>
    void
    refine_special_other(mesh_type& msh,
                            const cell_type& cl,
                            const IdxVector& vertices,
                            const size_t& cl_level)
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

            //iter it0 = pts.begin() + vertices.at(0);// 0
            //iter it1 = pts.begin() + vertices.at(1);// 3
            //iter it2 = pts.begin() + vertices.at(2);// 4
            //iter it3 = pts.begin() + pts.size();

            //if( it0  + 1 == it1)
            //if( it1  + 1 == it2)
            //if( it2 == it3)
            //std::cout << "INSIDE SPECIAL REFINEMENT" << std::endl;
        size_t i(0), j(0);
        auto num_vertices = vertices.size();
        for(auto & v:vertices)
        {
            int  factor = (i != num_vertices - 1)? 1 : num_fcs;
            auto next_v = vertices.at((i + 1) % num_vertices);

            if(  std::abs(next_v - v)  == factor)
            {
                n_gon<3>   nt;
                nt.p[0] = pts_ids[ v ];
                nt.p[1] = pts_ids[(v + 1) % num_fcs];
                nt.p[2] = cbar_id;

                nt.b[0] = storage->boundary_edges.at(fcs_ids.at(v));
                nt.b[1] = false;
                nt.b[2] = false;
                store(nt,m_triangles);
                l_triangles.push_back(std::make_pair(cl_level, cl.get_id()));

            }
            else
            {
                n_gon<4>   nt;
                nt.p[0] = pts_ids[v];
                nt.p[1] = pts_ids[(v + 1) % num_fcs];
                nt.p[2] = pts_ids[(v + 2) % num_fcs];
                nt.p[3] = cbar_id;

                nt.b[0] = storage->boundary_edges.at(fcs_ids.at(v));
                nt.b[1] = storage->boundary_edges.at(fcs_ids.at((v + 1) % num_fcs));
                nt.b[2] = false;
                nt.b[3] = false;

                store(nt,m_quadrangles);
                l_quadrangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            i++;
        }
    }

    struct point_info
    {
        point_info(){};
        size_t id;
        std::pair<size_t, int> cells;
        size_t point_id;
        size_t face_id;
        bool operator < (const point_info& pi) const
        {
            return (id < pi.id);
        }

        friend std::ostream& operator <<(std::ostream& os, const point_info& p)
        {
            os << "Point info:";
            os << " * id : "<< p.id << std::endl;
            os << " * point_id : "<< p.point_id << std::endl;
            os << " * face_id  : "<< p.face_id << std::endl;
            os << " * cells    :"<<std::endl;
            os << "   "<< p.cells.first << "  "<<p.cells.second  <<std::endl;
            return os;
        }
    };

#if 0
    template<typename IdxVector>
    void
    refine_special_triangle(mesh_type& msh,
                            const cell_type& cl,
                            const IdxVector& vertices,
                            const size_t& cl_level)
    {
        bool do_bar;
        size_t i(0), j(0);
        auto storage  = msh.backend_storage();

        auto pts      = points(msh,cl);
        auto pts_ids  = cl.point_ids();
        auto fcs_ids  = cl.faces_ids();


        //std::vector<std::array<int,5>> owner_cells;

        std::vector<point_info> owner_cells;
        owner_cells.reserve(pts_ids.size());

        //iter it0 = pts.begin() + vertices.at(0);// 0
        //iter it1 = pts.begin() + vertices.at(1);// 3
        //iter it2 = pts.begin() + vertices.at(2);// 4
        //iter it3 = pts.begin() + pts.size();

        //if( it0  + 1 == it1)
        //if( it1  + 1 == it2)
        //if( it2 == it3)
        //std::cout << "INSIDE SPECIAL REFINEMENT" << std::endl;

        for(auto & v:vertices)
        {
            //auto it     = (i != 2)? vertices.at(i + 1) :  pts.size() - 1;
            int  factor = (i != 2)? 1 : 0;
            auto iter_v = pts.begin() + v + factor;
            auto iter_next_vertex = (i != 2)?  pts.begin() + vertices.at(i + 1) : pts.end() - 1;

            if( iter_v  == iter_next_vertex)
            {
                auto fc_id  = fcs_ids.at(v);
                auto bar_id = faces_marks.at(fc_id).second;

                point_info npb, npv;
                npv.id = v + j;                     npv.cells    = std::make_pair(i , -1);
                npb.id = v + j + 1;                 npb.cells    = std::make_pair(i , ( i + 1)%3);
                npv.point_id = pts_ids.at(v);       npv.face_id  = fc_id;
                npb.point_id = bar_id;              npb.face_id  = fc_id;

                //owner_cells.push_back({ v + j, i , -1 , int(pts_ids.at(v)), int(fc_id)});
                //owner_cells.push_back({ v + 1 + j, i , (i +1)%3 , int(bar_id), int(fc_id)});
                owner_cells.push_back(npv);
                owner_cells.push_back(npb);
                j++;
            }
            else
            {
                auto w   = vertices.at((i + 1)% 3);
                auto l_dist  =  0.;
                auto r_dist  =  ( *std::next(pts.begin() , w) -  *std::next(pts.begin(), v)).to_vector().norm();
                auto b_dist  =  0.5 * (r_dist - l_dist) ;

                if( std::abs(b_dist) < 1.e-10)
                    throw std::logic_error("Flat triangle");

                auto r_limit  = (i != 2)? vertices.at(i + 1) :  pts.size();
                for(size_t ip = v; ip < r_limit; ip++)
                {
                    auto fc_id  = fcs_ids.at(ip);
                    auto x      = ( *std::next(pts.begin(), ip) - *std::next(pts.begin() , v) ).to_vector().norm();
                    auto dist   = (2. * x - r_dist - l_dist ) / (r_dist - l_dist);

                    point_info  npx, npf;
                    npx.id  = ip + j;
                    npx.point_id = pts_ids.at(ip);
                    npx.face_id  = fc_id;

                    if( std::abs(dist) < 1.e-10)
                    {
                        npx.cells    = std::make_pair(i , (i + 1) % 3);
                        npf.cells    = std::make_pair((i + 1) % 3 , -1);
                    }
                    else
                    {
                        if(dist < 0.)
                            npx.cells = std::make_pair(i , -1);

                        if(dist > 0.)
                            npx.cells    = std::make_pair((i + 1) % 3 , -1);

                        npf.cells = npx.cells;
                    }
                    owner_cells.push_back(npx);

                    auto fc_mark = faces_marks.at(fc_id).first;
                    if(fc_mark)
                    {
                        npf.id       = ip + j + 1;
                        npf.face_id  = fc_id;
                        npf.point_id = faces_marks.at(fc_id).second;
                        owner_cells.push_back(npf);
                        j++;
                    }
                }

            }

            i++;
        }


        std::sort(owner_cells.begin(), owner_cells.end());
        auto cid  = cl.get_id();
        size_t ii = 0;
        n_gon<3> t3;

        for (auto& np : owner_cells)
        {
                auto is_barycenter = ( np.cells.second != -1)?  true : false;
                if(is_barycenter)
                {
                    t3.p.at(ii) = np.point_id;
                    t3.b.at(ii) = false;
                    ii++;
                }
        }

        if(ii != 3)
            throw std::logic_error(" Triangle in the middle has less/more than 3 faces");

        store(t3,m_triangles);
        l_triangles.push_back(std::make_pair(cl_level,cl.get_id()));

        for(size_t j = 0; j < 3; j++)
        {
            std::vector<size_t> pts_ids;
            std::vector<size_t> fcs_ids;

            pts_ids.reserve(3);
            fcs_ids.reserve(3);

            for (auto i = owner_cells.begin(); i != owner_cells.end(); ++i)
            {
                auto np = *i;
                auto in_triangle = (np.cells.first == j) || (np.cells.second == j)?  true : false;
                if(in_triangle)
                {
                    auto  pt_id = np.point_id;
                    auto  fc_id = np.face_id;
                    pts_ids.push_back(pt_id);
                    fcs_ids.push_back(fc_id);

                }
            }

            auto new_num_fcs = fcs_ids.size();
            if(new_num_fcs == 3)
            {
                auto nt = put_n_gon<3>(msh, fcs_ids, pts_ids);
                store(nt,m_triangles);
                l_triangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 4)
            {
                auto nt = put_n_gon<4>(msh, fcs_ids, pts_ids);
                store(nt,m_quadrangles);
                l_quadrangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 5)
            {
                auto nt = put_n_gon<5>(msh, fcs_ids, pts_ids);
                store(nt,m_pentagons);
                l_pentagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 6)
            {

                auto nt = put_n_gon<6>(msh, fcs_ids, pts_ids);
                store(nt,m_hexagons);
                l_hexagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 7)
            {
                auto nt = put_n_gon<7>(msh, fcs_ids, pts_ids);
                store(nt,m_heptagons);
                l_heptagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 8)
            {
                auto nt = put_n_gon<8>(msh, fcs_ids, pts_ids);
                store(nt,m_octagons);
                l_octagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 9)
            {
                std::cout << "store in 9angles" << std::endl;
                auto nt = put_n_gon<9>(msh, fcs_ids, pts_ids);
                store(nt,m_enneagons);
                l_enneagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 10)
            {
                std::cout << "store in 10angles" << std::endl;
                auto nt = put_n_gon<10>(msh, fcs_ids, pts_ids);
                store(nt,m_decagons);
                l_decagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 11)
            {
                std::cout << "store in 11angles" << std::endl;

                auto nt = put_n_gon<11>(msh, fcs_ids, pts_ids);
                store(nt,m_hendecagons);
                l_hendecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 12)
            {
                std::cout << "store in 12angles" << std::endl;

                auto nt = put_n_gon<12>(msh, fcs_ids, pts_ids);
                store(nt,m_dodecagons);
                l_dodecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }

            if(new_num_fcs == 13)
            {
                std::cout << "store in 13angles" << std::endl;

                auto nt = put_n_gon<13>(msh, fcs_ids, pts_ids);
                store(nt,m_triadecagons);
                l_triadecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 14)
            {
                std::cout << "store in 14angles" << std::endl;

                auto nt = put_n_gon<14>(msh, fcs_ids, pts_ids);
                store(nt,m_tesseradecagons);
                l_tesseradecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 15)
            {
                std::cout << "store in 15angles" << std::endl;

                auto nt = put_n_gon<15>(msh, fcs_ids, pts_ids);
                store(nt,m_pentadecagons);
                l_pentadecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }

            if(new_num_fcs > 15)
            {
                std::cout << "number of faces = "<< new_num_fcs << std::endl;
                throw std::logic_error("number of faces exceeds maximum. Add new array to store it.");
            }
        }
    }
#endif
    template<typename IdxVector>
    void
    refine_special_triangle(mesh_type& msh,
                            const cell_type& cl,
                            const IdxVector& vertices,
                            const size_t& cl_level)
    {
        //std::cout << "INSIDE REFINE SPECIAL TRIANGLE" << std::endl;
        //std::cout << " - CELL : "<< cl.get_id() << std::endl;

        bool do_bar;
        size_t i(0), j(0);
        auto storage  = msh.backend_storage();

        auto pts      = points(msh,cl);
        auto pts_ids  = cl.point_ids();
        auto fcs_ids  = cl.faces_ids();

        std::vector<point_info> owner_cells;
        owner_cells.reserve(pts_ids.size());

        auto cell_id  = cl.get_id();
        auto hgi   = hang_nodes_vec.at(cell_id);

        for(auto & v:vertices)
        {
            if(!hgi.faces_has_hang_nodes.at(v))
            {
                //std::cout << "make bar" << std::endl;
                auto fc_id  = fcs_ids.at(v);
                auto bar_id = faces_marks.at(fc_id).second;

                point_info npb, npv;
                npv.id = v + j;                     npv.cells    = std::make_pair(i , -1);
                npb.id = v + j + 1;                 npb.cells    = std::make_pair(i , ( i + 1)%3);
                npv.point_id = pts_ids.at(v);       npv.face_id  = fc_id;
                npb.point_id = bar_id;              npb.face_id  = fc_id;

                owner_cells.push_back(npv);
                owner_cells.push_back(npb);
                j++;
            }
            else
            {
                auto w   = vertices.at((i + 1)% 3);
                auto l_dist  =  0.;
                auto r_dist  =  ( *std::next(pts.begin() , w) -  *std::next(pts.begin(), v)).to_vector().norm();
                auto b_dist  =  0.5 * (r_dist - l_dist) ;

                if( std::abs(b_dist) < 1.e-8)
                    throw std::logic_error("Flat triangle or angle too small");

                auto fc_id  = fcs_ids.at(v);
                auto next_fc_id  = fcs_ids.at((v + 1)% pts.size());

                auto fc_mark = faces_marks.at(fc_id).first;
                auto next_fc_mark = faces_marks.at(next_fc_id).first;

                //1.
                point_info  npv;
                npv.id  = v + j;
                npv.point_id = pts_ids.at(v);
                npv.face_id  = fc_id;
                npv.cells = std::make_pair(i , -1);
                owner_cells.push_back(npv);

                if(!fc_mark && !next_fc_mark)
                {
                    point_info  npf;

                    //2.
                    npf.id  = v + j + 1;
                    npf.point_id = pts_ids.at((v + 1)%pts.size());
                    npf.face_id  = next_fc_id;
                    npf.cells    = std::make_pair(i , (i + 1) % 3);

                    assert( msh.is_boundary(fc_id) == msh.is_boundary(next_fc_id) );

                    owner_cells.push_back(npf);
                    j++;

                    //std::cout << "  Owner_cells:"<< i << std::endl;
                    //for (auto& np : owner_cells)
                    //    std::cout << np << std::endl;
                }
                else
                {
                    if(fc_mark)
                    {
                        point_info  npf;
                        npf.id       = v + j + 1;
                        npf.face_id  = fc_id;
                        npf.point_id = faces_marks.at(fc_id).second;
                        npf.cells = std::make_pair(i , -1);

                        owner_cells.push_back(npf);
                        j++;
                    }

                    point_info  npb;
                    npb.id  = v + 1 + j;
                    npb.point_id = pts_ids.at((v + 1) %pts.size());
                    npb.face_id  = next_fc_id;
                    npb.cells    = std::make_pair(i , (i + 1) % 3);
                    owner_cells.push_back(npb);

                    if(next_fc_mark)
                    {
                        point_info  npf;
                        npf.id       = v + j + 2;
                        npf.face_id  = next_fc_id;
                        npf.point_id = faces_marks.at(next_fc_id).second;
                        npf.cells    = std::make_pair((i + 1) % 3 , -1);

                        owner_cells.push_back(npf);
                        j++;
                    }
                    //std::cout << "  Owner_cells:"<< i << std::endl;
                    //for (auto& np : owner_cells)
                    //    std::cout << np << std::endl;
                }

            }

            i++;
        }


        std::sort(owner_cells.begin(), owner_cells.end());
        auto cid  = cl.get_id();
        size_t ii = 0;
        n_gon<3> t3;

        for (auto& np : owner_cells)
        {
                auto is_barycenter = ( np.cells.second != -1)?  true : false;
                if(is_barycenter)
                {
                    t3.p.at(ii) = np.point_id;
                    t3.b.at(ii) = false;
                    ii++;
                }
        }

        std::cout << "  Inside triangle: " << std::endl;
        std::cout << "   * "<< t3 << std::endl;
        if(ii != 3)
            throw std::logic_error(" Triangle in the middle has less/more than 3 faces");

        store(t3,m_triangles);
        l_triangles.push_back(std::make_pair(cl_level,cl.get_id()));

        for(size_t j = 0; j < 3; j++)
        {
            std::vector<size_t> pts_ids;
            std::vector<size_t> fcs_ids;

            pts_ids.reserve(3);
            fcs_ids.reserve(3);

            for (auto i = owner_cells.begin(); i != owner_cells.end(); ++i)
            {
                auto np = *i;
                auto in_triangle = (np.cells.first == j) || (np.cells.second == j)?  true : false;
                if(in_triangle)
                {
                    auto  pt_id = np.point_id;
                    auto  fc_id = np.face_id;
                    pts_ids.push_back(pt_id);
                    fcs_ids.push_back(fc_id);

                }
            }

            auto new_num_fcs = fcs_ids.size();
            if(new_num_fcs == 3)
            {
                auto nt = put_n_gon<3>(msh, fcs_ids, pts_ids);

                std::cout << "  Outside triangle: " << std::endl;
                std::cout << "   * "<< nt << std::endl;

                store(nt,m_triangles);
                l_triangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 4)
            {
                auto nt = put_n_gon<4>(msh, fcs_ids, pts_ids);

                std::cout << "  Outside triangle: " << std::endl;
                std::cout << "   * "<< nt << std::endl;

                store(nt,m_quadrangles);
                l_quadrangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 5)
            {
                auto nt = put_n_gon<5>(msh, fcs_ids, pts_ids);

                std::cout << "  Outside triangle: " << std::endl;
                std::cout << "   * "<< nt << std::endl;

                store(nt,m_pentagons);
                l_pentagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }


            if(new_num_fcs > 5)
            {
                std::cout << "number of faces = "<< new_num_fcs << std::endl;
                throw std::logic_error("number of faces exceeds maximum. At least more than 1 hanging node is allowed, triangle cann't have more than 5 faces.");
            }
        }
    }

#if 0
    template<typename IdxVector>
    void
    refine_special_triangle(mesh_type& msh,
                            const cell_type& cl,
                            const IdxVector& vertices,
                            const size_t& cl_level)
    {
        if(vertices.size()>3)
            throw std::logic_error("This is not a triangle. Check refine_single_with_hanging_nodes or number of vertices");
        bool do_bar;
        size_t i(0), j(0);
        auto storage  = msh.backend_storage();

        auto pts      = points(msh,cl);
        auto pts_ids  = cl.point_ids();
        auto fcs_ids  = cl.faces_ids();


        //std::vector<std::array<int,5>> owner_cells;

        std::vector<point_info> owner_cells;
        owner_cells.reserve(pts_ids.size());

        //iter it0 = pts.begin() + vertices.at(0);// 0
        //iter it1 = pts.begin() + vertices.at(1);// 3
        //iter it2 = pts.begin() + vertices.at(2);// 4
        //iter it3 = pts.begin() + pts.size();

        //if( it0  + 1 == it1)
        //if( it1  + 1 == it2)
        //if( it2 == it3)
        //std::cout << "INSIDE SPECIAL REFINEMENT" << std::endl;
        std::cout << "owner_cells :" << std::endl;
        for (auto& np : owner_cells)
            std::cout << np << std::endl;

        std::cout << "POINTS ("<< cl.get_id()<<"):" << std::endl;
        for(auto & p:pts)
        {
            std::cout <<  p.x() << " "<< p.y() << std::endl;
        }
        for(auto & v:vertices)
        {
            std::cout <<" ****** vertex :"<< v <<" *******" << std::endl;
            //auto it     = (i != 2)? vertices.at(i + 1) :  pts.size() - 1;
            int  factor = (i != 2)? 1 : 0;
            auto iter_v = pts_ids.begin() + v + factor;
            auto iter_next_vertex = (i != 2)?  pts_ids.begin() + vertices.at(i + 1) : pts_ids.end() - 1;

            if( iter_v  == iter_next_vertex)
            {
                auto fc_id  = fcs_ids.at(v);
                auto bar_id = faces_marks.at(fc_id).second;

                point_info npb, npv;
                npv.id = v + j;                     npv.cells    = std::make_pair(i , -1);
                npb.id = v + j + 1;                 npb.cells    = std::make_pair(i , ( i + 1)%3);
                npv.point_id = pts_ids.at(v);       npv.face_id  = fc_id;
                npb.point_id = bar_id;              npb.face_id  = fc_id;

                //owner_cells.push_back({ v + j, i , -1 , int(pts_ids.at(v)), int(fc_id)});
                //owner_cells.push_back({ v + 1 + j, i , (i +1)%3 , int(bar_id), int(fc_id)});
                owner_cells.push_back(npv);
                owner_cells.push_back(npb);
                j++;

                std::cout << " ADD BAR  "<< i << std::endl;
                std::cout << "  1-"<<npv<< std::endl;
                std::cout << "  2-"<<npb<< std::endl;

                std::cout << "owner_cells :" << std::endl;
                for (auto& np : owner_cells)
                    std::cout << np << std::endl;

            }
            else
            {
                auto w   = vertices.at((i + 1)% 3);
                auto l_dist  =  0.;
                auto r_dist  =  ( *std::next(pts.begin() , w) -  *std::next(pts.begin(), v)).to_vector().norm();
                auto b_dist  =  0.5 * (r_dist - l_dist) ;

                if( std::abs(b_dist) < 1.e-8)
                    throw std::logic_error("Flat triangle or angle too small");
                auto fc_id  = fcs_ids.at(v);
                auto ffc_id = fcs_ids.at(v+1);

                point_info  npx, npf;

                //1.
                npx.id  = v + j;
                npx.point_id = pts_ids.at(v);
                npx.face_id  = fc_id;
                npx.cells    = std::make_pair(i , -1);

                //2.
                npf.id  = v + j + 1;
                npf.point_id = pts_ids.at(v + 1);
                npf.face_id  = ffc_id;
                npf.cells    = std::make_pair(i , (i + 1) % 3);

                assert( msh.is_boundary(fc_id) == msh.is_boundary(ffc_id) );

                owner_cells.push_back(npx);
                owner_cells.push_back(npf);
                j++;

            }

            i++;


        }

        std::cout << "owner_cells :" << std::endl;
        for (auto& np : owner_cells)
            std::cout << np << std::endl;

        std::sort(owner_cells.begin(), owner_cells.end());

        auto cid  = cl.get_id();
        size_t ii = 0;
        n_gon<3> t3;

        for (auto& np : owner_cells)
        {
                auto is_barycenter = ( np.cells.second != -1)?  true : false;
                if(is_barycenter)
                {
                    t3.p.at(ii) = np.point_id;
                    t3.b.at(ii) = false;
                    ii++;
                }
        }

        if(ii != 3)
            throw std::logic_error(" Triangle in the middle has less/more than 3 faces");

        store(t3,m_triangles);
        l_triangles.push_back(std::make_pair(cl_level,cl.get_id()));

        for(size_t j = 0; j < 3; j++)
        {
            std::vector<size_t> pts_ids;
            std::vector<size_t> fcs_ids;

            pts_ids.reserve(3);
            fcs_ids.reserve(3);

            for (auto i = owner_cells.begin(); i != owner_cells.end(); ++i)
            {
                auto np = *i;
                auto in_triangle = (np.cells.first == j) || (np.cells.second == j)?  true : false;
                if(in_triangle)
                {
                    auto  pt_id = np.point_id;
                    auto  fc_id = np.face_id;
                    pts_ids.push_back(pt_id);
                    fcs_ids.push_back(fc_id);

                }
            }

            auto new_num_fcs = fcs_ids.size();
            if(new_num_fcs == 3)
            {
                auto nt = put_n_gon<3>(msh, fcs_ids, pts_ids);
                store(nt, m_triangles);
                l_triangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 4)
            {
                auto nt = put_n_gon<4>(msh, fcs_ids, pts_ids);
                store(nt,m_quadrangles);
                l_quadrangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
        }
    }
#endif
    void
    refine_single_with_hanging_nodes(mesh_type& msh,
                                     const cell_type& cl,
                                     const size_t& cl_id)
    {
        auto fcs_ids  = cl.faces_ids();
        auto num_fcs  = fcs_ids.size();
        auto pts      = points(msh, cl);
        auto cid      = cl.get_id();
        auto cl_level = levels.at(cid);
        auto cl_mark = cells_marks.at(cid);
        /* Neighbors Adaptation:
        To adapt marked cells(3) and/or their neighbors(2)
        cells and neihgbors => if(cl_mark > 1)
        Only cells          => if(cl_mark == 3)*/
        if(cl_mark == 3)
        {
            if(msh.is_special_cell(cl))
            {
                auto vertices_pos =  msh.get_vertices_pos(cl);

                switch (vertices_pos.size())
                {
                    case 1:
                    case 2:
                        throw std::logic_error("Number of faces cannot be less than 3");
                        break;
                    case 3:
                        refine_special_triangle(msh, cl, vertices_pos,cl_level);
                        break;
                    //case 4:
                    //    refine_quadrangle();
                    //    break;
                    default:
                        refine_special_other(msh, cl, vertices_pos,cl_level); //WK: try to do the same for quadrangles
                        break;
                }
            }
            else
            {
                switch (num_fcs)
                {
                    case 1:
                    case 2:
                        throw std::logic_error("Number of faces cannot be less than 3");
                        break;
                    case 3:
                        refine_triangle(msh, cl, cl_level);
                        break;
                    //case 4:
                    //    refine_quadrangle();
                    //    break;
                    default:
                        refine_other(msh ,cl, cl_level);
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
                    auto fc_mark = faces_marks.at(id).first;
                    if(fc_mark)
                        count++;
                }
            }
            auto new_num_fcs = count + num_fcs;
            if(new_num_fcs == 3)
            {
                auto nt = put_n_gon_with_hanging_node<3>(msh, cl);
                store(nt, m_triangles);
                l_triangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 4)
            {
                auto nt = put_n_gon_with_hanging_node<4>(msh, cl);
                store(nt,m_quadrangles);
                l_quadrangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 5)
            {
                auto nt = put_n_gon_with_hanging_node<5>(msh, cl);
                store(nt,m_pentagons);
                l_pentagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 6)
            {

                auto nt = put_n_gon_with_hanging_node<6>(msh, cl);
                store(nt,m_hexagons);
                l_hexagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 7)
            {
                auto nt = put_n_gon_with_hanging_node<7>(msh, cl);
                store(nt,m_heptagons);
                l_heptagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 8)
            {
                auto nt = put_n_gon_with_hanging_node<8>(msh, cl);
                store(nt,m_octagons);
                l_octagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 9)
            {
                auto nt = put_n_gon_with_hanging_node<9>(msh, cl);
                store(nt,m_enneagons);
                l_enneagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 10)
            {
                auto nt = put_n_gon_with_hanging_node<10>(msh, cl);
                store(nt,m_decagons);
                l_decagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 11)
            {
                auto nt = put_n_gon_with_hanging_node<11>(msh, cl);
                store(nt,m_hendecagons);
                l_hendecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 12)
            {
                auto nt = put_n_gon_with_hanging_node<12>(msh, cl);
                store(nt,m_dodecagons);
                l_dodecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 13)
            {
                auto nt = put_n_gon_with_hanging_node<13>(msh, cl);
                store(nt,m_triadecagons);
                l_triadecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 14)
            {
                auto nt = put_n_gon_with_hanging_node<14>(msh, cl);
                store(nt,m_tesseradecagons);
                l_tesseradecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 15)
            {
                auto nt = put_n_gon_with_hanging_node<15>(msh, cl);
                store(nt,m_pentadecagons);
                l_pentadecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs > 15)
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
    auto
    refine( const mesh_type & old_msh,
            const size_t    & degree)
    {
        mesh_type re_msh, new_mesh, msh;

        msh = old_msh;
        re_msh  = msh;

        auto storage    = msh.backend_storage();
        auto re_storage = re_msh.backend_storage();
        auto num_cells  = msh.cells_size();
        size_t nds_counter = 0;

        face_marker(re_msh);

        auto other_info = mp.summary + "_RC" + tostr(imsh);
        auto marks_filename = mp.directory + "/marks" + other_info + ".txt";
        auto mesh_filename  = mp.directory + "/mesh"  + other_info + ".m";

        dump_to_matlab(msh, mesh_filename, cells_marks);
        save_data(cells_marks, marks_filename);

        for(auto& cl : msh)
        {
            auto cl_id   = cl.get_id();
            auto cl_mark = cells_marks.at(cl_id);

            if(cl_mark > 0)
            {
                if(mp.hanging_nodes)
                    refine_single_with_hanging_nodes(msh, cl, cl_id);
                else
                    refine_single(msh,cl);
            }
            else
            {
                auto cl_level =  levels.at(cl_id);
                auto fcs_ids  =  faces(msh,cl); //WK: There should be direct way to know the number of faces of the cell
                int  num_fcs  =  fcs_ids.size();
                switch(num_fcs)
                {
                    case 3:
                        store(put_n_gon<3>(msh,cl), m_triangles);
                        l_triangles.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 4:
                        store(put_n_gon<4>(msh,cl), m_quadrangles);
                        l_quadrangles.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 5:
                        store(put_n_gon<5>(msh,cl), m_pentagons);
                        l_pentagons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 6:
                        store(put_n_gon<6>(msh,cl), m_hexagons);
                        l_hexagons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 7:
                        store(put_n_gon<7>(msh,cl), m_heptagons);
                        l_heptagons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 8:
                        store(put_n_gon<8>(msh,cl), m_octagons);
                        l_octagons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 9:
                        store(put_n_gon<9>(msh,cl), m_enneagons);
                        l_enneagons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 10:
                        store(put_n_gon<10>(msh,cl), m_decagons);
                        l_decagons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 11:
                        store(put_n_gon<11>(msh,cl), m_hendecagons);
                        l_hendecagons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 12:
                        store(put_n_gon<12>(msh,cl), m_dodecagons);
                        l_dodecagons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 13:
                        store(put_n_gon<13>(msh,cl), m_triadecagons);
                        l_triadecagons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 14:
                        store(put_n_gon<14>(msh,cl), m_tesseradecagons);
                        l_tesseradecagons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 15:
                        store(put_n_gon<15>(msh,cl), m_pentadecagons);
                        l_pentadecagons.push_back(std::make_pair(cl_level, cl.get_id()));
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

        level_ancestor = l_triangles;
        level_ancestor.insert(level_ancestor.end(), l_quadrangles.begin(), l_quadrangles.end());
        level_ancestor.insert(level_ancestor.end(), l_pentagons.begin()  , l_pentagons.end());
        level_ancestor.insert(level_ancestor.end(), l_hexagons.begin()   , l_hexagons.end());
        level_ancestor.insert(level_ancestor.end(), l_heptagons.begin()  , l_heptagons.end());
        level_ancestor.insert(level_ancestor.end(), l_octagons.begin()   , l_octagons.end());
        level_ancestor.insert(level_ancestor.end(), l_enneagons.begin()  , l_enneagons.end());
        level_ancestor.insert(level_ancestor.end(), l_decagons.begin()   , l_decagons.end());
        level_ancestor.insert(level_ancestor.end(), l_hendecagons.begin(), l_hendecagons.end());
        level_ancestor.insert(level_ancestor.end(), l_dodecagons.begin() , l_dodecagons.end());
        level_ancestor.insert(level_ancestor.end(), l_triadecagons.begin()   , l_triadecagons.end());
        level_ancestor.insert(level_ancestor.end(), l_tesseradecagons.begin(), l_tesseradecagons.end());
        level_ancestor.insert(level_ancestor.end(), l_pentadecagons.begin()  , l_pentadecagons.end());

        LoaderType  loader;
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
        loader.m_triadecagons = m_triadecagons;
        loader.m_tesseradecagons  = m_tesseradecagons;
        loader.m_pentadecagons    = m_pentadecagons;

        loader.m_points       = re_storage->points;
        loader.m_boundary_edges = m_boundary_edges;

        auto open1 = save_adapted_mesh(loader);

        auto storage_rm2  = new_mesh.backend_storage();
        loader.populate_mesh(new_mesh);

        auto open2 = save_levels_info(loader, level_ancestor);

        new_levels = std::vector<size_t>(level_ancestor.size());
        ancestors  = std::vector<size_t>(level_ancestor.size());

        for(size_t i = 0; i < level_ancestor.size(); i++)
        {
            size_t idx = loader.m_index_transf.at(i);
            auto pair  = *std::next(level_ancestor.begin(), idx);
            new_levels.at(i) = pair.first;
            ancestors.at(i)  = pair.second;
        }
        std::cout << "levels_size     :" << levels.size()<< std::endl;
        std::cout << "level_ancestor_size:" << level_ancestor.size()<< std::endl;
        std::cout << "new_msh_size    :" << msh.cells_size()<< std::endl;

        for(auto& l : new_levels)
            std::cout << l;
        std::cout<< std::endl;

        return new_mesh;
    }
    template <typename LoaderType>
    bool
    save_adapted_mesh(const LoaderType& loader)
    {
        auto info_other =  mp.summary + "_" + mp.short_mesh_name + ".typ1";
        auto filename   =  mp.directory + "/amesh_" + tostr(imsh) + info_other;

        auto points   = loader.m_points;
        auto edges    = loader.m_edges;
        auto boundary_edges = loader.m_boundary_edges;

        std::ofstream ofs(filename);
        if (!ofs.is_open())
            std::cout << "Error opening file"<<std::endl;

        ofs << "vertices"    << std::endl;
        ofs <<  points.size()<< std::endl;
        for(auto& p : points)
            ofs << p.x() <<"  "<< p.y()<<std::endl;

        fvca5_save_tuples(ofs, loader.m_triangles);
        fvca5_save_tuples(ofs, loader.m_quadrangles);
        fvca5_save_tuples(ofs, loader.m_pentagons);
        fvca5_save_tuples(ofs, loader.m_hexagons);
        fvca5_save_tuples(ofs, loader.m_heptagons);
        fvca5_save_tuples(ofs, loader.m_octagons);
        //fvca5_save_tuples(ofs, loader.agons);

        ofs << "edges of the boundary"<< std::endl;
        ofs << boundary_edges.size()<< std::endl;
        for(auto& be : boundary_edges)
        {
            for(auto c : be)
                ofs << c + 1<< "  ";
            ofs << std::endl;
        }
        ofs << "all edges"<< std::endl;
        ofs << edges.size()<< std::endl;
        for(auto& e : edges)
        {
            for(auto c : e)
                ofs << c + 1 << "  ";
            ofs << std::endl;
        }
    }


    template <typename LoaderType>
    bool
    save_levels_info(const LoaderType& loader,
                     const std::vector<std::pair<size_t,size_t>>& ancestor_info)
    {
        auto other_info = mp.summary + "_R" + tostr(imsh);
        auto filename = mp.directory + "/levels" + other_info + ".txt";

        std::ofstream lfs(filename);
        if (!lfs.is_open())
            std::cout << "Error opening file : "<< filename<<std::endl;

        lfs << ancestor_info.size() <<std::endl;
        for(size_t i = 0; i < ancestor_info.size(); i++)
        {
            auto idx   = loader.m_index_transf.at(i);
            auto pair  = *std::next(ancestor_info.begin(), idx);
            auto ancestor_level = pair.first;
            auto ancestor_id    = pair.second;

            lfs << idx << " "<< ancestor_level << " "<<ancestor_id<< std::endl;
        }
        lfs.close();
    }
};

template<typename MeshType, typename T,
        typename CellBasisType, typename CellQuadType,
        typename FaceBasisType, typename FaceQuadType>
std::vector<dynamic_vector<T>>
proj_rec_solution(const MeshType& new_msh,
                  const MeshType& old_msh,
                  const std::vector<dynamic_matrix<T>>& grad_global,
                  const std::vector<dynamic_vector<T>>& Uh_old,
                  const std::vector<size_t>& ancestors_vec,
                  const size_t m_degree)
{
    std::cout << "INSIDE PROJ_REC_SOLUTION" << std::endl;
    typedef MeshType                mesh_type;
    typedef CellBasisType           cell_basis_type;
    typedef FaceBasisType           face_basis_type;
    typedef CellQuadType            cell_quadrature_type;
    typedef FaceQuadType            face_quadrature_type;
    typedef typename MeshType::point_type     point_type;
    typedef T  scalar_type;
    typedef dynamic_vector<scalar_type>     vector_type;
    typedef dynamic_matrix<scalar_type>     matrix_type;
    auto ret = solution_zero_vector(new_msh, m_degree);
    //std::cout << "ret_size"<< ret.size() << std::endl;
    cell_basis_type         cb(m_degree + 1);
    projector_nopre<mesh_type,
                        cell_basis_type,
                        cell_quadrature_type,
                        face_basis_type,
                        face_quadrature_type> projk(m_degree);

    // Borrar esto luego, es solo para comprobar que los ancestros estén bien
    #if 0
    for(auto& cell : new_msh)
    {
        auto cell_id  = cell.get_id();
        auto b  = barycenter(new_msh, cell);
        std::cout << "strName = strtrim(cellstr(num2str("<< ancestors.at(cell_id) <<",'(%d)')));"<<std::endl;
        std::cout << "text("<<b.x()<< ","<< b.y() <<",strName,'VerticalAlignment','bottom');"<<std::endl;
    }
    #endif
    for(auto& cell : new_msh)
    {
        auto cell_id  = cell.get_id();
        auto ancestor_id    = ancestors_vec.at(cell_id);
        auto ancestor_cell  = *std::next(old_msh.cells_begin(), ancestor_id);

        vector_type   ruh_acst = vector_type::Zero(cb.size());
        vector_type   uTF_acst = Uh_old.at(ancestor_id);
        matrix_type   rec_oper_acst = grad_global.at(ancestor_id);
        ruh_acst.tail(cb.size()-1)  = rec_oper_acst * uTF_acst;
        ruh_acst(0) = uTF_acst(0);

        //typedef std::function<T (const point_type & p, const size_t n)>    function;
        auto rec_fun =  [&](const point_type& p, const size_t number) -> scalar_type
        {
            scalar_type ret(0);

            auto phi  = cb.eval_functions(old_msh, ancestor_cell, p);

            for(size_t i = 0; i < cb.size(); i++)
                ret  += phi.at(i) * ruh_acst(i); //uTF_acst(i);

            return ret;
        };

        vector_type uTF = projk.compute_whole(new_msh, cell, rec_fun);
        ret.at(cell_id) = uTF;
    }

    return ret;
}



}//end Disk
