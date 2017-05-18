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
        mesh_name   = 0;
        num_remesh  = 0;
        marker_name = 2;
        initial_imsh = 0;
        percent     = 0.1;
        recycle     = true;
        diff        = false;
        mark_all    = false;
        hanging_nodes = true;
        call_mesher   = false;
    }
    T       percent;
    int     mesh_name;
    int     num_remesh;
    int     initial_imsh;
    bool    hanging_nodes;
    bool    call_mesher;
    int     marker_name;
    bool    recycle;
    bool    diff;
    bool    mark_all;
    std::string     short_mesh_name;
    std::string     directory;
    std::string     summary;
};

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
template <typename T>
void
read_data(std::vector<T>& vec, const std::string& filename)
{

    size_t elements_to_read;
    std::ifstream ifs(filename);

    if (!ifs.is_open())
        std::cout << "Error opening file"<<std::endl;

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
template<typename TensorsType,  typename T>
void
get_from_tensor(std::vector<dynamic_matrix<T>>& vec,
                const std::vector<TensorsType>& tsr_vec,
                const std::string& name)
{
    std::cout << "INSIDE get_from_tensor" << std::endl;
    vec = std::vector<dynamic_matrix<T>>(tsr_vec.size());
    std::cout << " * vec size: "<< vec.size() << std::endl;
    size_t val = 0;
    if(name == "sigma")
        val = 1;
    else if(name == "gamma")
        val = 2;
    else if(name == "xi_norm")
        val = 3;
    else
        std::cout << "WARNING: name not found in tensor matrix variables." << std::endl;

    std::cout << " * value = "<<val << std::endl;
    size_t i = 0;
    switch (val)
    {
        case 1:
            for(auto& tsr : tsr_vec)
                vec.at(i++) = tsr.siglam;
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

template<typename TensorsType,  typename T>
void
get_from_tensor(std::vector<T>& vec,
                const std::vector<TensorsType>& tsr_vec,
                const std::string& name)
{
    vec = std::vector<T>(tsr_vec.size());
    if(name != "quad_degree")
    {
        std::cout << "WARNING: The only scalar variable in tensor structure is";
        std::cout << "quad_degree. Review name or function get_from_tensor.";
        std::cout << std::endl;
        return;
    }
    size_t i = 0;
    for(auto& tsr : tsr_vec)
        vec.at(i++) = tsr.quad_degree;
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
            for(auto& mat : vec)
                tsr_vec.at(i++).siglam = mat;
            return;
        case 2:
            for(auto& mat : vec)
                tsr_vec.at(i++).gamma  = mat;
            return;
        case 3:
            for(auto& mat : vec)
                tsr_vec.at(i++).xi_norm = mat;
            return;
        default:
            throw std::invalid_argument("Not known name variable for tensor.");
            return;
    }
};
template<typename TensorsType, typename T>
void
put_tensor(std::vector<TensorsType>& tsr_vec,
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

    for(auto& mat : vec)
        tsr_vec.at(i++).quad_degree = mat;

    return;
};
#if 0
//WarningK: The first two load_data, should be just one. But I have to
//1 - Leave out index_transf and to be sure that I really dont need it.
//2 - Fix the thing of rmc and rm
//3 - step: imsh for levels_ancestors_vec
//    step: imsh-1 for Uh_Th and tsr since is info on the previous mesh

template<typename T, typename InputVectorType>
void
 load_data(InputVectorType& vec,
            const mesh_parameters<T>& mp,
            const std::string& name,
            const size_t step)
{
     typedef dynamic_vector<T> vector_type;
     typedef dynamic_matrix<T> matrix_type;

     auto info_other = mp.summary + "_rm" + tostr(step) +".txt";
     std::cout << "namefile to open:  "<< mp.directory + name + info_other << std::endl;
     read_data(vec, mp.directory + name + info_other);
 }
#endif
template<typename T>
void
load_data(std::vector<std::pair<size_t,size_t>>& levels_ancestors_vec,
                   std::vector<size_t>& index_transf,
                   const mesh_parameters<T>& mp,
                   const size_t imsh)
{
    typedef dynamic_vector<T> vector_type;
    typedef dynamic_matrix<T> matrix_type;

    auto info = mp.summary + "_R" + tostr(imsh) +".txt";
    std::cout << "/* namefile to open:  "<< mp.directory + "/levels" + info << std::endl;
    read_data(levels_ancestors_vec, index_transf, mp.directory + "/levels" + info);
    std::cout << "levels_ancestors reading succesful" << std::endl;
}
template< typename T>
void
load_data(std::vector<dynamic_vector<T>>& Uh_Th,
          const mesh_parameters<T>& mp,
           const size_t imsh)
{
    typedef dynamic_vector<T> vector_type;
    typedef dynamic_matrix<T> matrix_type;

    auto info = mp.summary + "_R" + tostr(imsh-1) +".txt";

    std::cout << "/* namefile to open:  "<< mp.directory + "/Uh" + info << std::endl;
    read_data(Uh_Th , mp.directory + "/Uh"  + info);

    std::cout << "Uh reading succesful" << std::endl;
}
template<typename TensorsType, typename T>
void
load_data(std::vector<TensorsType>& tsr_vec,
        const mesh_parameters<T>& mp,
        const size_t imsh)

{
    typedef dynamic_vector<T> vector_type;
    typedef dynamic_matrix<T> matrix_type;

    std::vector<matrix_type> sigma;
    std::vector<T> quadeg;

    auto info = mp.summary + "_R" + tostr(imsh-1) +".txt";

    get_from_tensor(sigma,  tsr_vec, "sigma");
    get_from_tensor(quadeg, tsr_vec, "quad_degree");

    read_data(sigma , mp.directory + "/Sigma" + info);
    read_data(quadeg, mp.directory + "/QuadDegree" + info);

    assert(sigma.size() == quadeg.size());
    if(sigma.size() != quadeg.size())
        throw std::logic_error("Sizes of loaded data doesn't match");

    tsr_vec = std::vector<TensorsType>(sigma.size());

    put_tensor(tsr_vec, sigma,  "sigma");
    put_tensor(tsr_vec, quadeg, "quad_degree");
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
    bool upos  = (std::abs(u) < 1.e-5 * invDenom) || ( u > 0.);
    bool vpos  = (std::abs(v) < 1.e-5 * invDenom) || ( v > 0.);
    bool uvpos = (std::abs(u + v - 1.) < 1.e-10) || (u + v < 1.);

    return (upos && vpos && uvpos);
};
template<typename MeshType, typename CellType>
auto
new_points(const MeshType& msh, const CellType& cell,
           const size_t quad_degree, const size_t num_sig_pts)
{
    typedef MeshType                                    mesh_type;
    typedef typename MeshType::cell                     cell_type;
    typedef typename MeshType::face                     face_type;
    typedef typename MeshType::point_type               point_type;
    typedef disk::quadrature<mesh_type, cell_type>      cell_quad_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quad_type;

    auto cq  = cell_quad_type(quad_degree);
    auto fq  = face_quad_type(quad_degree);

    auto ret = std::vector<point_type>(num_sig_pts);
    //std::cout << "size sig_pts ("<< cell.get_id()<<") : "<< num_sig_pts << std::endl;
    size_t cont = 0;

    auto cqs  = cq.integrate(msh, cell);
    for(auto& cq : cqs)
        ret.at(cont++) = cq.point();

    auto fcs = faces(msh, cell);

    for(auto& fc :fcs)
    {
        auto fqs = fq.integrate(msh, fc);
        for(auto& fq : fqs)
            ret.at(cont++) = fq.point();
    }

    auto pts = points(msh, cell);

    for(auto& p : pts)
        ret.at(cont++) = p;

    for(auto& fc :fcs)
        ret.at(cont++) = barycenter(msh, fc);

    ret.at(cont) = barycenter(msh, cell);

    return ret;
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

template<typename TensorType, typename T, typename Storage, typename CellType>
auto
compute_interpolation(const mesh<T,2,Storage>& msh,
              const CellType  & cell,
              const TensorType& tsr,
              const std::vector<point<T,2>>& pts_to_eval_sigma)
{
    typedef mesh<T,2,Storage>   mesh_type;
    typedef point<T,2>          point_type;
    typedef dynamic_vector<T>   vector_type;
    typedef dynamic_matrix<T>   matrix_type;

    matrix_type sigma = tsr.siglam;

    auto num_sig_pts = pts_to_eval_sigma.size();
    auto num_cols    = sigma.cols();
    matrix_type  ret = matrix_type::Zero(2, num_sig_pts);
    auto pts   = points( msh, cell);
    auto fcs   = faces(  msh, cell);
    auto area  = measure(msh, cell);
    auto cell_barycenter = barycenter(msh, cell);
    auto offset_vertex   = num_cols - 2 * pts.size() - 1;
    auto offset_face_bar = num_cols - pts.size() - 1;
    auto offset_cell_bar = num_cols - 1;

    if(area < 1.e-10)
        throw std::invalid_argument("Area < 1.e-10. Review how tiny an element could be.");

    size_t cont = 0;

    for(auto& ep : pts_to_eval_sigma)
    {
        assert(cont < num_sig_pts);

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

            matrix_type sig_at_vts1 = matrix_type::Zero(2,3);
            sig_at_vts1.col(0) = sigma.col(offset_cell_bar);
            sig_at_vts1.col(1) = sigma.col(offset_vertex   + ifc);
            sig_at_vts1.col(2) = sigma.col(offset_face_bar + ifc);

            if(is_inside(t1, ep))
            {
                ret.col(cont++) = trilinear_interpolation(sig_at_vts1, t1, ep);
                break;
            }


            triangle<T> t2;
            t2.points[0]  =  cell_barycenter;
            t2.points[1]  =  barycenter(msh, face);
            t2.points[2]  =  pts.at((ifc + 1)% pts.size());

            matrix_type sig_at_vts2 = matrix_type::Zero(2,3);
            sig_at_vts2.col(0) = sigma.col(offset_cell_bar);
            sig_at_vts2.col(1) = sigma.col(offset_face_bar + ifc);
            sig_at_vts2.col(2) = sigma.col(offset_vertex   + (ifc + 1)% pts.size());

            if(is_inside(t2, ep))
            {
                ret.col(cont++) = trilinear_interpolation(sig_at_vts2, t2, ep);
                break;
            }


        }
    }
    return ret;
};
template<typename TensorType, typename T, typename Storage >
std::vector<TensorType>
sigma_interpolation(const mesh<T,2,Storage>& new_msh,
                    const mesh<T,2,Storage>& old_msh,
                    const std::vector<TensorType>& tsr_vec,
                    const std::vector<size_t>& ancestors,
                    const size_t degree)
{
    std::cout << "INSIDE SIGMA INTERPOLATION" << std::endl;
    typedef mesh<T,2,Storage>   mesh_type;
    typedef dynamic_vector<T>   vector_type;
    typedef dynamic_matrix<T>   matrix_type;

    size_t num_new_cells = new_msh.cells_size();
    std::vector<TensorType> new_tsr_vec;
    new_tsr_vec  = disk::tensor_zero_vector(new_msh, degree);

    for(auto& cell : new_msh)
    {
        auto cell_id  = cell.get_id();
        auto ancestor_id    = ancestors.at(cell_id);
        auto ancestor_cell  = *std::next(old_msh.cells_begin(), ancestor_id);
        auto tsr      =  tsr_vec.at(ancestor_id);
        auto new_tsr  =  new_tsr_vec.at(cell_id);
        auto num_pts  =  new_tsr.siglam.cols();
        auto qdegree  =  tsr.quad_degree;
        //std::cout << "quad_degree_ancestor ("<< ancestor_id<<"): "<< qdegree << std::endl;
        auto pts_to_eval_sigma  = new_points(new_msh, cell, qdegree, num_pts);
        matrix_type sigma = compute_interpolation(new_msh, cell, tsr, pts_to_eval_sigma);
        new_tsr_vec.at(cell_id).siglam = sigma;
        //std::cout << "quad_degree_cell ("<< cell_id<<")        : "<< new_tsr_vec.at(cell_id).quad_degree << std::endl;
    }
    return new_tsr_vec;
};



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
           const std::vector<tensors<T>>& tsr_vec,
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
          const std::vector<tensors<T>>& tsr_vec,
          const T& yield,
          const T& percent)

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
    void
    refine(mesh_type & msh,
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

        msh = re_msh;
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

    stress_based_mesh(const mesh_type& msh,
                        const std::vector<size_t>& levels_vec,
                        const plasticity_data<T>&  m_pst,
                        const mesh_parameters<T>& msh_parameters,
                        const size_t& adaptive_step):
                        pst(m_pst),
                        mp(msh_parameters),
                        imsh(adaptive_step)
    {
        info  =  mp.summary + "_RC" + tostr(imsh);
        //check_older_msh(msh);
        cells_marks = std::vector<size_t>(msh.cells_size(),0);
        levels = levels_vec;
        level_ancestor.reserve(msh.cells_size());

        // This only to mark all the cells on each adaptation (thought as a test);
        // Comment this for normal adaptations. Also reveiw test adaptation if( imsh > 0)
        // (remember this does a  first adaptation  to ensure only triangles).
        if(mp.mark_all)
        {
            for(auto& cell: msh)
            {
                auto cell_id = cell.get_id();
                cells_marks.at(cell_id) = 3;
            }

            for(auto& b : cells_marks)
                std::cout<<b<<"  ";
            std::cout  << std::endl;

            // Esto creo que es necesqrio incluso para adapt 0, puesto que en ningun otro lado aparece
            // la inicialization con el tamanho de faces_marks. No se como estaba funcionando para imsh ==0
            faces_marks.resize(msh.faces_size());
            for(size_t i = 0; i < faces_marks.size(); i++)
               faces_marks.at(i).first = false;
        }
        else
        {
            //This is only to start only with triangles
            //if we want to have more polygons take out this and solve
            //the problem with the refinement_other
            //#if 0
            if(imsh == 0)
            {
                for(auto& cell: msh)
                {
                    auto cell_id = cell.get_id();
                    auto c_faces = faces(msh, cell);
                    if(c_faces.size() > 3)
                        cells_marks.at(cell_id) = 3;
                }

                for(auto& b : cells_marks)
                    std::cout<<b<<"  ";
                std::cout  << std::endl;

            }
            //#endif
            else
            {
                faces_marks.resize(msh.faces_size());
                for(size_t i = 0; i < faces_marks.size(); i++)
                    faces_marks.at(i).first = false;
            }
            //#endif
        }

    }
    template< size_t N>
    struct n_gon
    {
       std::array<size_t, N>  p;
       std::array<bool, N>    b;
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
                size_t ngh_id    = face_owner_cells_ids(msh, fc, cl);
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
                    check_levels(msh , search_cl);
                }
            }
        }
        return;
    }
    bool
    test_marker(const mesh_type& msh,
            const T& ratio)
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

            if(*(result.second) >=  ratio )
            {
                if (*(result.first) < ratio )
                {
                    cells_marks.at(cl_id) = 3;
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

    bool
    marker_xc(const mesh_type& msh,
            const std::vector<tensors<T>>& tsr_vec,
            const T& yield)
    {
        typedef disk::quadrature<mesh_type, cell_type>      cell_quad_type;
        typedef dynamic_vector<scalar_type>     vector_type;

        std::cout << "INSIDE_XC_MARKER" << std::endl;
        bool do_refinement(false);

        cells_marks = std::vector<size_t>(msh.cells_size());

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto quad_degree = tsr_vec.at(i).quad_degree;
            auto cl  = *std::next(msh.cells_begin() , i);
            auto cq  = cell_quad_type(quad_degree);
            auto cqs_size =  cq.integrate(msh, cl).size();
            auto num_pts  =  points(msh,cl).size();

            vector_type xi_norm = tsr_vec.at(i).xi_norm;
            auto mat1 = xi_norm.head(cqs_size);
            auto mat2 = xi_norm.tail(num_pts);

            bool in_interval(false);

            for(size_t m = 0; m < xi_norm.size(); m++)
            {
                if(std::isnan(xi_norm(m)))
                    throw std::logic_error("The norm of the constrain is NaN");
            }

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
            const std::vector<tensors<T>>& tsr_vec)
    {
        std::cout << "INSIDE_JB_MARKER" << std::endl;

        bool do_refinement(false);
        cells_marks = std::vector<size_t>(msh.cells_size());

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto mat = tsr_vec.at(i).xi_norm;

            for(size_t m = 0; m < mat.size(); m++)
            {
                if(std::isnan(mat(0,m)))
                    throw std::logic_error("The norm of the constrain is NaN");
            }

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
            const std::vector<tensors<T>>& tsr_vec)
    {
        std::cout << "INSIDE_M4_MARKER" << std::endl;

        bool do_refinement(false);
        cells_marks = std::vector<size_t>(msh.cells_size());

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto mat = tsr_vec.at(i).xi_norm;

            for(size_t m = 0; m < mat.size(); m++)
            {
                if(std::isnan(mat(0,m)))
                    throw std::logic_error("The norm of the constrain is NaN");
            }

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
            T eta_stress(0.);
            std::vector<T> Cg(2);
            auto varpi =  0.5;
            auto fc   = *(msh.faces_begin() + fc_id);
            auto fqs  = face_quadrature.integrate(msh, fc);
            auto h_F  = measure(msh, fc);

            auto owner_cells = face_owner_cells_ids(msh, fc);
            std::cout<<"FACE_ID = "<< fc_id<<std::endl;
            for(size_t ii = 0; ii < owner_cells.size(); ii++)
            {
                auto cl_id = owner_cells.at(ii);
                auto cl    = *std::next(msh.cells_begin() , size_t(cl_id));
                auto tsr   =  tsr_vec.at(cl_id);
                auto m_T   =  measure(msh, cl);
                auto h_T   =  diameter(msh, cl);
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

                    vector_type tau  = tsr.siglam.col(cont) - pst.alpha * tsr.gamma.col(cont);
                    vector_type str  = tau + pst.alpha * dphi_r_uh;
                    auto jump = make_prod_stress_n( str , n);
                    //auto proj_jump = make_prod_stress_n(tau + pst.alpha*dphi_r_uh , n);
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


    template<typename CellQuadType, typename FaceQuadType>
    size_t
    gauss_points_positions(const mesh_type& msh, const cell_type& cell,
                            const face_type& face, const size_t& quad_degree)
    {
        auto c_faces = faces(msh, cell);
        size_t   j = 0;
        #if 0
        for(auto& fc : faces)
        {
            sorted_faces.at(j) = std::make_pair(fc , j);
            j++;
        }
        std::sort(sorted_faces.begin(),sorted_faces.end());

        //typename face_type::id_type fid(fc_id);
        auto lower = std::lower_bound(sorted_faces.begin(), sorted_faces.end(),
                            std::make_pair(fc_id, 0), [](myPair lhs, myPair rhs)
                                -> bool { return lhs.first < rhs.first; });
        if(lower == sorted_faces.end())
            throw std::invalid_argument("This is a bug: face not found");

        auto find   =  std::distance(sorted_faces.begin() , lower);
        #endif
        bool find = false;
        for(auto& fc : c_faces)
        {
            if(fc == face)
            {
                find = true;
                break;
            }
            j++;
        }
        if(!find)
            throw std::invalid_argument("This is a bug: face not found");

        auto fc_pos =  j;
        auto cq  = CellQuadType(quad_degree);
        auto fq  = FaceQuadType(quad_degree);
        //WK: This should change if k different for each face
        auto fqs = fq.integrate(msh, face).size();
        auto cqs = cq.integrate(msh, cell).size();

        return cqs + fqs * fc_pos;
    }
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


    template<typename TensorType>
    T
    divergence_tau_P0(const mesh_type& msh, const cell_type& cl,
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

        auto vts_ids =  msh.get_vertices_ids(cl);
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
                mefs << "fill(coords(:,2), coords(:,3), color_mat("<<i + 1<<",:)) "<<std::endl;
                mefs<< "strName = strtrim(cellstr(num2str("<< eta.at(cont).second<<",'(%d)')));"<<std::endl;
                mefs<< "text("<<b.x()<< ","<< b.y() <<",strName,'VerticalAlignment','bottom');"<<std::endl;

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
            if(std::abs(value - new_value) < 1.e-15)
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
                mefs << "fill(coords(:,2), coords(:,3), color_mat("<<i + 1<<",:)) "<<std::endl;
                mefs<< "strName = strtrim(cellstr(num2str("<< eta.at(cont).second<<",'(%d)')));"<<std::endl;
                mefs<< "text("<<b.x()<< ","<< b.y() <<",strName,'VerticalAlignment','bottom');"<<std::endl;
                mefs<<  eta.at(cont).first <<std::endl;
                ++cont;
            }
        }
        auto vmax = eta.at(0).first;
        auto vmin = eta.at(eta.size()-1).first;
        auto h = (vmax - vmin) / 10.;
        mefs<< "h = colorbar;"<<std::endl;
        mefs<< "caxis(["<< vmin <<" "<< vmax <<"])"<<std::endl;
        mefs<< "set(h,'YTick',["<<vmin<<":"<< h <<":"<< vmax <<"])"<<std::endl;

        mefs.close();
        assert(cont == msh.cells_size());
        return;
    }

//#if 0
    template<typename CellBasisType, typename FaceBasisType, typename Solution>
    bool
    marker_ae( const mesh_type& msh,
            const std::vector<tensors<T>>& tsr_vec,
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
        //otherwise this should be inside the for(auto& cl_id : owner_cells)
        auto quad_degree = tsr_vec[0].quad_degree;

        cell_basis_type                         cell_basis(m_degree + 1);
        face_basis_type                         face_basis(m_degree);
        face_quadrature_type                    face_quadrature(quad_degree);
        cell_quadrature_type                    cell_quadrature(quad_degree);

        std::cout << "INSIDE_AE_MARKER" << std::endl;

        std::ofstream       ofs(mp.directory + "/estimator.txt", std::ios::app);
        if (!ofs.is_open())
            std::cout << "Error opening ofs"<<std::endl;

        bool do_refinement(false);

        cells_marks = std::vector<size_t>(msh.cells_size());

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto mat = tsr_vec.at(i).xi_norm;

            for(size_t m = 0; m < mat.size(); m++)
            {
                if(std::isnan(mat(0,m)))
                    throw std::logic_error("The norm of the constrain is NaN");
            }
        }
            //Marking with estimators

        //WK: this should be vector of vectors if u es a vectorial function,
        //since each tensor is a matrix. Then the jump = [tau + alpha G(u)].n
        // will be a vector and not a scalar.
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
            auto h_T       =  diameter(msh, cell);
            auto num_faces =  c_faces.size();

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

                size_t neighbor_id;
                if( !msh.is_boundary(face) )
                    neighbor_id = face_owner_cells_ids(msh, face, cell);
                else
                    neighbor_id = cell_id;

                auto neighbor  = *std::next(msh.cells_begin(),  neighbor_id);
                auto c_normal  =  normal(msh, cell, face);
                auto n_normal  =  normal(msh, neighbor, face);
                auto n_tsr     =  tsr_vec.at(neighbor_id);
                auto cgp_id    =  gauss_points_positions<cell_quadrature_type
                                                         ,face_quadrature_type>
                                                    (msh,cell,face, quad_degree);
                auto ngp_id    =  gauss_points_positions<cell_quadrature_type
                                                         ,face_quadrature_type>
                                                (msh,neighbor,face, quad_degree);
                matrix_type n_rec_oper =  grad_global.at(neighbor_id);
                vector_type nc_uh_TF   =  Uh_Th.at(neighbor_id);
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

                    vector_type c_tau  = c_tsr.siglam.col(cgp_id) - pst.alpha * c_tsr.gamma.col(cgp_id);
                    vector_type n_tau  = n_tsr.siglam.col(ngp_id) - pst.alpha * n_tsr.gamma.col(ngp_id);

                    //estimator only with plasticity
                    vector_type c_str  = c_tau + pst.alpha *  c_dphi_ruh;
                    vector_type n_str  = n_tau + pst.alpha * nc_dphi_ruh;


                    vector_type jump = c_str - n_str;
                    auto jp_value    = varpi * make_prod_stress_n( jump , c_normal);
                    //auto proj_jump = make_prod_stress_n(tau + pst.alpha*dphi_r_uh , n);
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

            //std::cout << "cell : "<< cell_id << std::endl;
            for (auto& qp : cqs)
            {
                    vector_type ddphi     =  cell_basis.eval_laplacians(msh, cell, qp.point());
                    vector_type ddphi_zm  =  take(ddphi, col_range);
                    //#if 0
                    T   lap_ruh =  ddphi_zm.dot(c_ruh);
                    matrix_type c_tau    =  c_tsr.siglam - pst.alpha * c_tsr.gamma;

                    auto div_tau  =  divergence_tau_P0(msh, cell, c_tau, qp.point(), pos);
                    auto f   = solution.f(qp.point());
                    eta_res += qp.weight() * iexp_pow(f + div_tau +  pst.alpha * lap_ruh, 2.);
                    //#endif 0

                    #if 0
                    T    lap_ruh  =  ddphi_zm.dot(c_ruh);
                    matrix_type c_tau    =  c_tsr.siglam - pst.alpha * c_tsr.gamma;;
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
                    #endif 0

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

        T set_error(0.);

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

            // L2_error
            matrix_type rec_oper =  grad_global.at(cell_id);
            vector_type uTF   =  Uh_Th.at(cell_id);
            vector_type x     =  rec_oper * uTF;

            dynamic_vector<scalar_type> true_dof = projk.compute_cell(msh, cell, solution.sf);
            dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
            dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
            auto err_dof = std::sqrt(diff_dof.dot(projk.cell_mm * diff_dof));

            //L2_error gradient
            T error_dfun(0.);
            typedef typename disk::solution<T,2>::gradient_vector   gradient_vector;
            auto cqs = cell_quadrature.integrate(msh, cell);

            for (auto& qp : cqs)
            {
                gradient_vector dphi_fun  =  solution.df(qp.point(),0);
                //Este cero es solo a causa del test con el prolema difusivo dividido en dos partes. Para tuyau y los otros puede ser cualquier nuemero

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

        T new_set_error(0.);
        std::cout << "total error : "<< set_error << std::endl;
        std::cout << "percent     : "<< mp.percent << std::endl;

        for(auto& e : eta)
        {
            auto cell_id   = e.second;
            new_set_error += e.first;

            cells_marks.at(cell_id) = 3;
            ++levels.at(cell_id);
            do_refinement  = true;

            std::cout << "cell : "<< cell_id <<std::endl;
            std::cout << " * eta   : "<< e.first << std::endl;
            std::cout << " * error : "<< new_set_error << std::endl;
            if( new_set_error >= mp.percent * set_error )
                break;
        }
        dump_to_matlab( msh, mp.directory + "/mesh_pr_" +  info + ".m", cells_marks);

        if(mp.marker_name == 6 )
        {
            for(size_t i = 0; i < tsr_vec.size();i++)
            {

                dynamic_matrix<T> sigmas = tsr_vec.at(i).siglam;

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

        error_to_matlab(msh, eta ,  "tot" + info);
        error_to_matlab(msh, etas,  "str" + info);
        error_to_matlab(msh, etar,  "res" + info);


        for(auto& cl : msh)
            check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }
    template<typename CellBasisType, typename FaceBasisType, typename Solution>
    bool
    marker( const mesh_type& msh,
            const std::vector<tensors<T>>& tsr_vec,
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
            auto h_T       =  diameter(msh, cell);
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

                size_t neighbor_id;
                if( !msh.is_boundary(face) )
                    neighbor_id = face_owner_cells_ids(msh, face, cell);
                else
                    neighbor_id = cell_id;

                auto neighbor  = *std::next(msh.cells_begin(),  neighbor_id);
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


                matrix_type n_rec_oper =  grad_global.at(neighbor_id);
                vector_type nc_uh_TF   =  Uh_Th.at(neighbor_id);
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
                    auto jp_value      = varpi * make_prod_stress_n( jump , c_normal);
                    //auto proj_jump = make_prod_stress_n(tau + pst.alpha*dphi_r_uh , n);
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

        T set_error(0.);

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
            T error_dfun(0.);
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

        T new_set_error(0.);
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
                        auto  neighbor_id   = face_owner_cells_ids(msh,fc,cl);
                        cmv_temp.at(neighbor_id) = 1;
                        auto  neighbor      = *std::next(msh.cells_begin(),size_t(neighbor_id));
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
                                    auto  nn_id   = face_owner_cells_ids(msh,nfc,cl);
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
                vertices = msh.get_vertices_ids(cl);
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
                            size_t neighbor_id = face_owner_cells_ids(msh, fc, cl);

                            typename cell_type::id_type id(neighbor_id);
                            auto  ngh_mark      = cells_marks.at(id);
                            if(cl_mark > ngh_mark)
                            {
                                cells_marks.at(neighbor_id) = cl_mark - 1;
                                auto  ngh  = *std::next(msh.cells_begin(),size_t(neighbor_id));
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

                            size_t neighbor_id = face_owner_cells_ids(msh, fc, cl);
                            typename cell_type::id_type id(neighbor_id);
                            auto  ngh_mark      = cells_marks.at(id);
                            if(cl_mark > ngh_mark)
                            {
                                cells_marks.at(neighbor_id) = cl_mark - 1;
                                auto  ngh  = *std::next(msh.cells_begin(), size_t(neighbor_id));
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
    set_marks_hanging_nodes(mesh_type & msh, const cell_type & cl)
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
                vertices = msh.get_vertices_ids(cl);
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
                            size_t neighbor_id = face_owner_cells_ids(msh, fc, cl);

                            typename cell_type::id_type id(neighbor_id);
                            auto  ngh_mark      = cells_marks.at(id);
                            if(cl_mark > ngh_mark)
                            {
                                cells_marks.at(neighbor_id) = cl_mark - 1;
                                auto  ngh  = *std::next(msh.cells_begin(),size_t(neighbor_id));
                                set_marks_hanging_nodes(msh, ngh);
                            }

                        }
                    }

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

                            size_t neighbor_id = face_owner_cells_ids(msh, fc, cl);
                            typename cell_type::id_type id(neighbor_id);
                            auto  ngh_mark      = cells_marks.at(id);
                            if(cl_mark > ngh_mark)
                            {
                                cells_marks.at(neighbor_id) = cl_mark - 1;
                                auto  ngh  = *std::next(msh.cells_begin(), size_t(neighbor_id));
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
                    auto  neighbor_id = face_owner_cells_ids(msh,fc,cl);
                    auto  ngh_mark    = cells_marks.at(neighbor_id);
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
                            cells_marks.at(neighbor_id) = 2;
                        else
                            cells_marks.at(neighbor_id) = cl_mark - 1;
                        auto  ngh = *std::next(msh.cells_begin(), size_t(neighbor_id));
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
        return false;
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

    template<typename IdxVector>
    void
    refine_special_triangle(mesh_type& msh,
                            const cell_type& cl,
                            const IdxVector& vertices,
                            const size_t& cl_level)
    {
        bool do_bar;
        size_t i(0), j(0);
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
        };
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
                auto b_dist  =  (r_dist - l_dist) / 2.;

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
                auto vertices =  msh.get_vertices_ids(cl);

                switch (vertices.size())
                {
                    case 1:
                    case 2:
                        throw std::logic_error("Number of faces cannot be less than 3");
                        break;
                    case 3:
                        refine_special_triangle(msh, cl, vertices,cl_level);
                        break;
                    //case 4:
                    //    refine_quadrangle();
                    //    break;
                    default:
                        refine_special_other(msh, cl, vertices,cl_level); //WK: try to do the same for quadrangles
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
    void refine( mesh_type & msh,
                const std::vector<vector_type>&  Uh_Th,
                const size_t& degree)
    {
        mesh_type re_msh, new_mesh;

        re_msh = msh;

        auto storage    = msh.backend_storage();
        auto re_storage = re_msh.backend_storage();
        auto num_cells  = msh.cells_size();
        size_t nds_counter = 0;

        face_marker(re_msh);
        dump_to_matlab(msh, mp.directory + "/mesh" + mp.summary +".m",cells_marks);

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

        old_msh = msh;
        msh     = new_mesh;
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
        typename FaceBasisType, typename FaceQuadType,
        typename PointType>
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
    typedef PointType               point_type;
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
            scalar_type ret(0.);

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
