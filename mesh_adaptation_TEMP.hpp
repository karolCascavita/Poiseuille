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
#include <list>
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
        load4marker = false;
        binNadapt   = false;
        binning     = false;
        recycle     = true;
        diff        = false;
        mark_all    = false; // There is a problem here read comment
        hanging_nodes = true;
        call_mesher   = false;
        hho = 0;
    }
    // !WARNINGK
    // If mark_all ==true, it cannot be used  in circular for iter == 0, that used quadrilaterals
    // Since, we have triangles in the center, that will make hanging_nodes to appear
    // However, since the quadrilaterals are split in triangles and doesn't support
    // hanging nodes, the triangles made by quadrilaterals won-t see the hanging_nodes
    int     hho;
    T       percent;
    int     mesh_name;
    int     num_remesh;
    bool    call_mesher;
    int     initial_imsh;
    int     marker_name;
    int     start_exact;
    bool    binNadapt;
    bool    load4marker;
    bool    recycle;
    bool    hanging_nodes;
    bool    diff;
    bool    mark_all;
    bool    binning;
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
        os << "* load4marker  : "<< mp.load4marker<< std::endl;
        os << "* binANDadapt  : "<< mp.binNadapt<< std::endl;
        os << "* hanging_nodes: "<< mp.hanging_nodes<< std::endl;
        os << "* binning      : "<< mp.binning<< std::endl;
        os << "* diffusion    : "<< mp.diff<< std::endl;
        os << "* mark_all     : "<< mp.mark_all<< std::endl;
        os << "* short_mesh_name: "<< mp.short_mesh_name<< std::endl;
        os << "* directory    : "<< mp.directory<< std::endl;
        os << "* summary      : "<< mp.summary<< std::endl;
        os << "* summary_old  : "<< mp.summary_old<< std::endl;
        os << "* hho       : "<< mp.hho<< std::endl;

        return os;
    }
};



template<typename T, size_t DIM, typename Storage>
std::vector<dynamic_vector<T>>
solution_zero_vector(const disk::mesh<T,DIM,Storage>& msh, const size_t degree,
                    const bool& change_degree)
{
    typedef disk::mesh<T,DIM,Storage>              mesh_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, typename mesh_type::cell>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, typename mesh_type::face>    face_basis_type;

    std::vector<dynamic_vector<T>> U(msh.cells_size());
    size_t i = 0;

    int degree_cell = degree + change_degree;
    for (auto& cl : msh)
    {
        auto fcs = faces(msh, cl);
        auto num_faces       = fcs.size();
        auto num_cell_dofs   = cell_basis_type(degree_cell).size();
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
        case 9 :
            ofs << "enneagons"<<std::endl;
            break;
        case 10 :
            ofs << "decagons" <<std::endl;
            break;
        case 11 :
            ofs << "hendecagons" <<std::endl;
            break;
        case 12 :
            ofs << "dodecagons" <<std::endl;
            break;
        case 13:
            ofs << "triadecagons" <<std::endl;
            break;
        case 14:
            ofs << "tesseradecagons" <<std::endl;
            break;
        case 15:
            ofs << "pentadecagons"<<std::endl;
            break;
        case 16:
            ofs << "16gons"<<std::endl;
            break;
        case 17:
            ofs << "17gons"<<std::endl;
            break;
        case 18:
            ofs << "18gons"<<std::endl;
            break;
        case 19:
            ofs << "19gons"<<std::endl;
            break;
        case 20:
            ofs << "20gons"<<std::endl;
            break;
        case 21:
            ofs << "21gons"<<std::endl;
            break;
        case 22:
            ofs << "22gons"<<std::endl;
            break;
        case 23:
            ofs << "23gons"<<std::endl;
            break;
        case 24:
            ofs << "24gons"<<std::endl;
            break;
        case 25:
            ofs << "25gons"<<std::endl;
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


template <typename T>
void
save_data(const std::vector<dynamic_vector<T>>& vec,
          const std::string& filename)
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
    ofs.close();
};
template <typename T>
void
save_data(const std::vector<dynamic_vector<T>>& vec,
        std::ofstream & ofs)
{
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
save_data(const std::vector<dynamic_matrix<T>>& vec,
          const std::string& filename)
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
    ofs.close();
};

template <typename T>
void
save_data(const std::vector<dynamic_matrix<T>>& vec,
          std::ofstream & ofs)
{
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
save_data(const std::vector<T>  & vec,
          const std::string     & filename)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file"<<std::endl;

    ofs << vec.size()<<std::endl;

    for(auto& v :  vec)
        ofs << v <<std::endl;
    ofs.close();
};
template< typename MeshType>
void
save_data(const tensors<MeshType>& mpt,
    const std::string& filename)
{
    std::cout << "INSIDE save_data tensor" << std::endl;

    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file: "<< filename<<std::endl;

    typedef typename MeshType::scalar_type   scalar_type;
    std::vector<dynamic_vector<scalar_type>>  sigma = mpt.at_all_cells();
    auto varsigma = mpt.at_all_faces();
    save_data( sigma, ofs);
    save_data( varsigma, ofs);

    ofs.close();
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

template <typename T>
void
read_data(std::vector<dynamic_vector<T>>& vec,
          const std::string& filename)
{
    size_t elements_to_read;
    size_t values_to_read;

    std::ifstream ifs(filename);

    if (!ifs.is_open())
        std::cout << "Error opening file: "<< filename<<std::endl;
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

    ifs.close();
    return;
};
template <typename T>
void
read_data(std::vector<dynamic_vector<T>>& vec,
          std::ifstream& ifs)
{
    size_t elements_to_read;
    size_t values_to_read;

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
read_data(std::vector<dynamic_matrix<T>>& vec,
          const std::string& filename)
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

template<typename T>
void
read_data(std::vector<dynamic_matrix<T>>& vec,
          std::ifstream& ifs)
{
    std::cout << "INSIDE read_data matrix" << std::endl;

    size_t elements_to_read;
    size_t cols, rows;

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

template <typename T>
void
read_data(std::vector<T>& vec,
        const std::string& filename)
{

    size_t elements_to_read;
    std::ifstream ifs(filename);

    if (!ifs.is_open())
        std::cout << "Error opening file: "<< filename<<std::endl;
    std::cout << "Opening file: "<< filename <<std::endl;

    ifs >> elements_to_read;
    std::cout << "Attempting to read " << elements_to_read << " values" << std::endl;

    vec = std::vector<T>(elements_to_read);

    for (size_t i = 0; i < elements_to_read; i++)
        ifs >> vec.at(i);

    ifs.close();

    return;
};
template <typename T>
void
read_data(std::vector<T>& vec,
         std::ifstream& ifs)
{
    size_t elements_to_read;

    ifs >> elements_to_read;
    std::cout << "Attempting to read " << elements_to_read << " values" << std::endl;

    vec = std::vector<T>(elements_to_read);

    for (size_t i = 0; i < elements_to_read; i++)
        ifs >> vec.at(i);

    return;
};
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
template<typename TensorsType, typename MeshType, typename T>
void
load_data(TensorsType& mpt,
        const MeshType& msh,
        const mesh_parameters<T>& mp,
        const std::string& R,
        const size_t step)
{
    std::cout << "INSIDE load_data tensor" << std::endl;

    auto info = mp.summary_old + "_" + R + tostr(step-1) +".txt";


    auto filename =  mp.directory + "/multiplier" + info;

    std::ifstream ifs(filename);
    if (!ifs.is_open())
        std::cout << "Error opening file: "<< filename<<std::endl;
    std::cout << "Opening file: "<< filename <<std::endl;

    std::vector<dynamic_vector<T>> mpt_at_cells =  mpt.at_all_cells();
    std::vector<dynamic_matrix<T>> mpt_at_faces =  mpt.at_all_faces();

    read_data( mpt_at_cells, ifs);
    read_data( mpt_at_faces, ifs);

    assert(mpt.at_all_cells_size() == mpt.at_all_cells_size());
    ifs.close();
    return;

};

template<typename MeshType, typename T,
        typename CellBasisType, typename CellQuadType,
        typename FaceBasisType, typename FaceQuadType>
std::vector<dynamic_vector<T>>
proj_rec_solution(const MeshType& new_msh,
                  const MeshType& old_msh,
                  const std::vector<dynamic_matrix<T>>& reconstruction_opers,
                  const std::vector<dynamic_vector<T>>& Uh_old,
                  const std::vector<size_t>& ancestors_vec,
                  const size_t m_degree,
                  const bool& change_degree)
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
    auto ret = solution_zero_vector(new_msh, m_degree, change_degree);
    //std::cout << "ret_size"<< ret.size() << std::endl;
    cell_basis_type         cb(m_degree + 1);
    projector_nopre<mesh_type,
                        cell_basis_type,
                        cell_quadrature_type,
                        face_basis_type,
                        face_quadrature_type> projk(m_degree, change_degree);

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
        matrix_type   rec_oper_acst = reconstruction_opers.at(ancestor_id);
        ruh_acst.tail(cb.size()-1)  = rec_oper_acst * uTF_acst;
        ruh_acst(0) = uTF_acst(0);

        //typedef std::function<T (const point_type & p, const size_t n)>    function;
        auto rec_fun =  [&](const point_type& p, const size_t number) -> scalar_type
        {
            scalar_type ret(0);

            auto phi  = cb.eval_functions(old_msh, ancestor_cell, p);

            return phi.dot(ruh_acst);
        };

        vector_type uTF = projk.compute_whole(new_msh, cell, rec_fun);
        ret.at(cell_id) = uTF;
    }

    return ret;
}

template<typename MeshType, typename T,
        typename CellBasisType, typename CellQuadType,
        typename FaceBasisType, typename FaceQuadType>
dynamic_vector<T>
lift_face_stress(const MeshType &msh,
                 const typename MeshType::cell & cell,
                 const dynamic_matrix<T> &  varsigma,
                 const size_t m_degree)
{
    typedef CellBasisType           cell_basis_type;
    typedef FaceBasisType           face_basis_type;
    typedef CellQuadType            cell_quadrature_type;
    typedef FaceQuadType            face_quadrature_type;
    typedef dynamic_vector<T>       vector_type;
    typedef dynamic_matrix<T>       matrix_type;

    cell_basis_type         cell_basis(m_degree);
    face_basis_type         face_basis(m_degree);
    cell_quadrature_type    cell_quadrature(2*(m_degree+1));
    face_quadrature_type    face_quadrature(2*m_degree);

    auto num_cell_dofs = cell_basis.range(0, m_degree).size();
    auto DIM = MeshType::dimension;
    auto vec_cbs = num_cell_dofs * DIM;
    matrix_type stiff_mat = matrix_type::Zero(vec_cbs, vec_cbs);

    /* 1 - Lifting of varsgima */
    auto cell_quadpoints = cell_quadrature.integrate(msh, cell);
    for (auto& qp : cell_quadpoints)
    {
        auto c_phi = cell_basis.eval_functions(msh, cell, qp.point());
        matrix_type vec_phi = make_vectorial_matrix(c_phi, DIM);
        stiff_mat += qp.weight() * vec_phi.transpose() * vec_phi;
    }

    /* LHS: take basis functions derivatives from degree 1 to K+1 */
    auto MG_ldlt = stiff_mat.llt();
    auto fcs = faces(msh, cell);
    auto num_faces     = fcs.size();
    auto num_face_dofs = face_basis.size();
    auto BG_row_size   = vec_cbs;

    dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

    vector_type ret = vector_type::Zero(BG_row_size);

    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        auto current_face_range = dsr.face_range(face_i);
        auto fc = fcs[face_i];
        auto n = normal(msh, cell, fc);
        auto face_quadpoints = face_quadrature.integrate(msh, fc);

        matrix_type FG = matrix_type::Zero(BG_row_size, num_face_dofs);

        for (auto& qp : face_quadpoints)
        {
            auto c_phi = cell_basis.eval_functions(msh, cell, qp.point());
            matrix_type vec_phi = make_vectorial_matrix(c_phi, DIM);
            // (  w * n , v_F )
            auto f_phi = face_basis.eval_functions(msh, fc, qp.point());
            FG += qp.weight() * vec_phi.transpose() * n * f_phi.transpose() ;
        }
        auto  pos = face_position(msh, cell, fc);
        vector_type f_varsig  = varsigma.col(pos);
        ret += MG_ldlt.solve(FG) * f_varsig;
    }
    return ret;
}

template<typename Iterator>
std::vector<std::pair<Iterator, typename std::iterator_traits<Iterator>::difference_type>> find_range
(
    Iterator begin,
    Iterator end,
    const typename std::iterator_traits<Iterator>::value_type& val
)
{
    std::vector<std::pair<Iterator, typename std::iterator_traits<Iterator>::difference_type>> res;
    for (Iterator it = std::find(begin, end, val);
         it != end; it = std::find(std::next(it), end, val))
        res.push_back(std::make_pair(it, std::distance(begin, it)));
    return res;
};


template<typename MeshType, typename T,
        typename CellBasisType, typename CellQuadType,
        typename FaceBasisType, typename FaceQuadType, typename TensorsType>
TensorsType
project_stresses_lifted(const MeshType & new_msh,
              const MeshType    & old_msh,
              const TensorsType & mpt,
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

    TensorsType ret;
    ret.template zero_vector<CellBasisType,FaceBasisType>(new_msh, m_degree);
    //std::cout << "ret_size"<< ret.size() << std::endl;
    cell_basis_type         cell_basis(m_degree + 1);
    face_basis_type         face_basis(m_degree);
    cell_quadrature_type    cell_quadrature(2*(m_degree+1));
    face_quadrature_type    face_quadrature(2*m_degree);

    projector_pst<mesh_type,
                        cell_basis_type,
                        cell_quadrature_type,
                        face_basis_type,
                        face_quadrature_type> projk(m_degree);

    for(auto& cell : new_msh)
    {
        auto cell_id  = cell.get_id();
        auto ancestor_id    = ancestors_vec.at(cell_id);
        auto ancestor_cell  = *std::next(old_msh.cells_begin(), ancestor_id);

        vector_type   sigma = mpt.at_cell(old_msh, ancestor_cell);
        ret.save(new_msh, cell, sigma);
    }

    auto DIM = MeshType::dimension;
    auto col_range  =   cell_basis.range(0, m_degree);
    auto row_range  =   dof_range(0, mesh_type::dimension);

    for(auto& ancestor_cell : old_msh)
    {
        /* 1- search child cells */
        auto ancestor_id = ancestor_cell.get_id();
        auto cells_range = find_range(ancestors_vec.begin(), ancestors_vec.end(), ancestor_id);

        for(auto& cr: cells_range)
        {
            int cell_id = cr.second;
            auto cell   = *std::next(new_msh.cells_begin(), cell_id);

            /* 2- lifting */
            matrix_type  old_varsigma =  mpt.at_element_faces(old_msh, ancestor_cell);

            vector_type  lift_dofs = lift_face_stress< mesh_type, T,
                                     cell_basis_type, cell_quadrature_type,
                                     face_basis_type, face_quadrature_type>
                                     (old_msh, ancestor_cell, old_varsigma, m_degree);

            /* 3 - Projection on faces of child cells*/
            //typedef std::function<T (const point_type & p, const size_t n)>    function;
            auto fcs = faces(new_msh, cell);
            auto number = size_t(0);
            for (auto & fc :fcs)
            {
                auto n = normal(new_msh, cell, fc);
                auto lifted_fun =  [&](const point_type& p, const size_t number)-> scalar_type
                {
                    scalar_type ret(0);

                    auto phi  = cell_basis.eval_functions(old_msh, ancestor_cell, p);
                    auto vec_phi = make_vectorial_matrix( phi, DIM);

                    ret =  n.transpose() * vec_phi * lift_dofs;
                    return ret;
                };
                vector_type new_varsigma = projk.compute_face(new_msh, cell,
                                                            fc, lifted_fun);
                ret.save(new_msh, cell, fc, new_varsigma);
            }
        }
    }

    return ret;
};

template<typename Point>
std::vector<Point>
sort_points(const std::vector<Point>& v)
{
    auto ret = v;
    std::stable_sort(ret.begin(), ret.end(), [&](const Point& a, const Point& b)
        {
            return (a.x() < b.x());
        }
    );
    std::stable_sort(ret.begin(), ret.end(), [&](const Point& a, const Point& b)
        {
            if(a.x() == b.x())
                return (a.y() < b.y());
            return false;
        }
    );
    return ret;
};


template<typename MeshType, typename T,
        typename CellBasisType, typename CellQuadType,
        typename FaceBasisType, typename FaceQuadType, typename TensorsType>
TensorsType
project_stresses(const MeshType & new_msh,
                  const MeshType    & old_msh,
                  const TensorsType & mpt,
                  const std::vector<size_t>& ancestors_vec,
                  const size_t m_degree,
                  const bool& change_degree)
{
    typedef MeshType                mesh_type;
    typedef CellBasisType           cell_basis_type;
    typedef FaceBasisType           face_basis_type;
    typedef CellQuadType            cell_quadrature_type;
    typedef FaceQuadType            face_quadrature_type;
    typedef typename MeshType::point_type     point_type;
    typedef T  scalar_type;
    typedef dynamic_vector<scalar_type>     vector_type;
    typedef dynamic_matrix<scalar_type>     matrix_type;

    TensorsType ret;
    ret.template zero_vector<CellBasisType,FaceBasisType>(new_msh, m_degree);
    //std::cout << "ret_size"<< ret.size() << std::endl;
    cell_basis_type         cell_basis(m_degree + 1);
    face_basis_type         face_basis(m_degree);
    cell_quadrature_type    cell_quadrature(2*(m_degree+1));
    face_quadrature_type    face_quadrature(2*m_degree);

    projector_pst<mesh_type,
                        cell_basis_type,
                        cell_quadrature_type,
                        face_basis_type,
                        face_quadrature_type> projk(m_degree, change_degree);

    for(auto& cell : new_msh)
    {
        auto cell_id  = cell.get_id();
        auto ancestor_id    = ancestors_vec.at(cell_id);
        auto ancestor_cell  = *std::next(old_msh.cells_begin(), ancestor_id);

        vector_type   sigma = mpt.at_cell(old_msh, ancestor_cell);
        ret.save(new_msh, cell, sigma);

        auto ancestor_fcs = faces(old_msh, ancestor_cell);
        auto c_fcs = faces(new_msh, cell);
        auto c_pts = points(new_msh, cell);
        for(auto fc : c_fcs)
        {
            auto f_pts = points(new_msh, fc);

            for(auto afc :  ancestor_fcs)
            {
                auto a_pts = points(old_msh, afc);

                auto fpts = sort_points(f_pts);
                auto apts = sort_points(a_pts);
                bool face_found(false);
                if( fpts[0].x() == apts[0].x() && fpts[0].y() == apts[0].y())
                {
                    if( fpts[1].x() == apts[1].x() && fpts[1].y() == apts[1].y())
                        face_found = true;
                }

                if(face_found)
                {
                    auto varsigma =  mpt.at_face(old_msh, ancestor_cell, afc);
                    ret.save(new_msh, cell, fc, varsigma);
                    break;
                }
            }
        }
    }

    return ret;
}

template<typename MeshType>
auto
marks_yield(const MeshType& msh, const tensors<MeshType>& stresses,
                const disk::plasticity_data<typename MeshType::scalar_type>& pst,
                const size_t degree)
{
    typedef MeshType                     mesh_type;
    typedef typename MeshType :: cell             cell_type;
    typedef typename MeshType :: face             face_type;
    typedef typename MeshType::scalar_type        scalar_type;

    typedef disk::quadrature<mesh_type, cell_type>     cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>     face_quadrature_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    typedef dynamic_vector<scalar_type>                 vector_type;
    typedef dynamic_matrix<scalar_type>                 matrix_type;

    std::cout << "INSIDE_YIELD_MARKER " << std::endl;
    std::cout << "degree : "<< degree << std::endl;
    bool do_refinement(false);

    cell_basis_type      cb(degree);
    face_basis_type      fb(degree);

    cell_quadrature_type cq(2 * degree);
    face_quadrature_type fq(2 * degree);

    auto DIM = mesh_type::dimension;
    auto marks = std::vector<size_t>(msh.cells_size());

    for(auto & cl : msh)
    {
        auto fcs = faces(msh, cl);
        auto cell_id = msh.lookup(cl);
        vector_type   sigma    = stresses.at_cell(msh, cl);
        matrix_type   varsigma = stresses.at_element_faces(msh, cl);

        auto cbar = barycenter(msh, cl);
        auto cphi = cb.eval_functions(msh, cl, cbar);
        matrix_type vec_phi = make_vectorial_matrix(cphi, 2);
        vector_type sigma_eval = vec_phi * sigma;

        auto  sigma_norm  = sigma_eval.norm();

        if((sigma_norm >= pst.yield))
            marks.at(cell_id) = 1;
    }
    return marks;
}

template<typename MeshType>
class stress_based_mesh
{};

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
    std::vector<size_t>  binning_marks;
    std::vector<size_t>  new_binning_marks;

    std::vector<size_t>  levels;
    std::vector<std::pair<size_t,size_t>> level_ancestor;
    std::vector<size_t>  new_levels;
    std::vector<size_t>  ancestors;
    std::vector<vector_type> new_Uh;
    //std::vector<std::pair<bool,typename point<T,2>::id_type>> faces_marks;
    std::vector<std::pair<bool, int>> faces_marks;
    std::vector<std::array<ident_impl_t, 4>>        m_edges;
    std::vector<std::array<ident_impl_t, 3>>        m_boundary_edges;

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
    std::vector<std::array<ident_impl_t,16>>        m_16gons;
    std::vector<std::array<ident_impl_t,17>>        m_17gons;
    std::vector<std::array<ident_impl_t,18>>        m_18gons;
    std::vector<std::array<ident_impl_t,19>>        m_19gons;
    std::vector<std::array<ident_impl_t,20>>        m_20gons;
    std::vector<std::array<ident_impl_t,21>>        m_21gons;
    std::vector<std::array<ident_impl_t,22>>        m_22gons;
    std::vector<std::array<ident_impl_t,23>>        m_23gons;
    std::vector<std::array<ident_impl_t,24>>        m_24gons;
    std::vector<std::array<ident_impl_t,25>>        m_25gons;

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
    std::vector<std::pair<size_t ,size_t>>        l_16gons;
    std::vector<std::pair<size_t ,size_t>>        l_17gons;
    std::vector<std::pair<size_t ,size_t>>        l_18gons;
    std::vector<std::pair<size_t ,size_t>>        l_19gons;
    std::vector<std::pair<size_t ,size_t>>        l_20gons;
    std::vector<std::pair<size_t ,size_t>>        l_21gons;
    std::vector<std::pair<size_t ,size_t>>        l_22gons;
    std::vector<std::pair<size_t ,size_t>>        l_23gons;
    std::vector<std::pair<size_t ,size_t>>        l_24gons;
    std::vector<std::pair<size_t ,size_t>>        l_25gons;

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
       std::array<size_t, N>    bid;


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
    test_marker_all(const mesh_type& msh,
                    const size_t  imsh)
    {
        std::cout << " INSIDE MARKER ALL" << std::endl;
        bool do_refinement = false;

        for(auto& cl : msh)
        {
            auto cl_id    = cl.get_id();
            cells_marks.at(cl_id) = 3;
            if(imsh > 0)
                ++levels.at(cl_id);
            do_refinement = true;
        }

        dump_to_matlab(msh, mp.directory +  info + "_wl.m", cells_marks);

        for(auto& cl : msh)
           check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }
    bool
    test_marker_ratio(const mesh_type& msh,
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

            #if 0
            std::cout << "cell : "<< cl_id << std::endl;
            std::cout << " * pts : ";
            for(auto& p: pts)
                std::cout << p <<"  ";
            std::cout<< std::endl;
            std::cout << " * r_max: "<< *(result.second) << std::endl;
            std::cout << " * r_min: "<< *(result.first) << std::endl;
            #endif

            if(*(result.second) >=  ratio )
            {
                if (*(result.first) < ratio )
                {
                    //std::cout << " * cell marked" << std::endl;
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
    check_for_NaN(const tensors<mesh_type>& tsr)
    {
        for(auto itor = tsr.at_cells_begin(); itor != tsr.at_cells_end(); itor++)
        {
            vector_type mat = *itor;

            for(size_t m = 0; m < mat.rows(); m++)
            {
                if(std::isnan(mat(m,0)))
                    throw std::logic_error("The norm of the constrain is NaN");
            }
        }
        for(auto itor = tsr.at_faces_begin(); itor != tsr.at_faces_end(); itor++)
        {
            auto mat = *itor;

            for(size_t m = 0; m < mat.rows(); m++)
                for(size_t n = 0; n < mat.cols(); n++)
                    if(std::isnan(mat(m,n)))
                    {    throw std::logic_error("The norm of the constrain is NaN");}
        }
        return;
    }
    // ******************************  Marker  *********************************
    bool
    marker_prueba(const mesh_type& msh,
                  const tensors<mesh_type>& tsr,
                  const size_t degree)
    {

        typedef disk::quadrature<mesh_type, cell_type>     cell_quadrature_type;
        typedef disk::quadrature<mesh_type, face_type>     face_quadrature_type;

        typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
        typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

        typedef dynamic_vector<scalar_type>                 vector_type;

        std::cout << "INSIDE_MARKER_PRUEBA" << std::endl;
        bool do_refinement(false);

        cell_basis_type      cell_basis(degree + 1);
        cell_quadrature_type cell_quadrature (2 * degree + 2);

        auto DIM = mesh_type::dimension;
        auto one_range = cell_basis.range(1, degree + 1);
        auto row_range = disk::dof_range( 0, DIM);


        cells_marks = std::vector<size_t>(msh.cells_size());

        auto   in_interval = false;
        size_t cell_count  = 0;

        check_for_NaN(tsr);

        for(auto itor = tsr.at_cells_begin(); itor != tsr.at_cells_end(); itor++)
        {
            vector_type sigma = *itor;
            auto cl = *std::next(msh.cells_begin(), cell_count);
            auto cell_quadpoints = cell_quadrature.integrate(msh, cl);
            for (auto& qp : cell_quadpoints)
            {
                auto c_dphi = cell_basis.eval_gradients(msh, cl, qp.point());
                matrix_type c_dphi_0   = make_gradient_matrix(c_dphi, row_range, one_range);
                vector_type sigma_eval = c_dphi_0 * sigma;
                auto sigma_norm = sigma_eval.norm();
                if(std::abs(sigma_norm - pst.yield) <= 1.e-10 || std::abs(sigma_norm< pst.yield ))
                {
                    if(std::abs(sigma_norm - pst.yield) > 1.e-10 && sigma_norm > pst.yield )
                        in_interval = true;
                }
            }

            std::cout << "in_interval :"<< in_interval << std::endl;
            if(in_interval)
            {
                cells_marks.at(cell_count) = 3;
                ++levels.at(cell_count);
                do_refinement = true;
            }
            cell_count++;
        }

        dump_to_matlab(msh, mp.directory + "/mesh_pr_" +  info + ".m", cells_marks);

        for(auto& cl : msh)
           check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }



    auto
    test_marker_yield(const mesh_type& msh, const tensors<mesh_type>& stresses,
                    const disk::plasticity_data<T>& pst, const size_t degree)
    {


        typedef dynamic_vector<scalar_type>                 vector_type;
        typedef dynamic_matrix<scalar_type>                 matrix_type;

        std::cout << "INSIDE_YIELD_MARKER " << std::endl;
        std::cout << "degree : "<< degree << std::endl;
        bool do_refinement(false);

        auto DIM = mesh_type::dimension;
        cells_marks = std::vector<size_t>(msh.cells_size());
        auto initial_marks = marks_yield(msh, stresses, pst, degree);

        dump_to_matlab(msh, mp.directory + "/mesh_initial_marks" +  info + ".m", initial_marks);

        auto node_state_marks = nodes_state(msh, initial_marks);

        for(auto cl : msh)
        {
            auto cid = msh.lookup(cl);
            auto nodes = cl.point_ids();

            bool is_solid(false), is_liquid(false), is_mix(false);

            for(auto n: nodes)
            {
                switch(node_state_marks.at(n))//simplificar esto
                {
                    case 0:
                        is_solid = true;
                        break;
                    case 1:
                        is_liquid = true;
                        break;
                    case 2:
                        is_mix = true;
                        break;
                }
            }

            if(is_solid && is_liquid || is_mix)
            {
                cells_marks.at(cid) = 3;
                ++levels.at(cid);
                do_refinement = true;
            }
        }

        dump_to_matlab(msh, mp.directory + "/refinement_marks" +  info + ".m", cells_marks);

        for(auto& cl : msh)
           check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }

    template<typename VariableType>
    auto
    marker(const mesh_type& msh,
            const VariableType & var,
            const disk::plasticity_data<scalar_type>& pst,
            const size_t m_degree)
    {
        typedef disk::quadrature<mesh_type, cell_type>     cell_quadrature_type;
        typedef disk::quadrature<mesh_type, face_type>     face_quadrature_type;

        typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
        typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

        cell_basis_type         cb(m_degree);
        face_basis_type         fb(m_degree);
        cell_quadrature_type    cq(2 * m_degree);
        face_quadrature_type    fq(2 * m_degree);

        auto DIM = mesh_type::dimension;
        std::vector<size_t> initial_marks(msh.cells_size());
        size_t num_cells_marked = 0;
        for(auto & cl : msh)
        {
            auto fcs = faces(msh, cl);
            auto cell_id = msh.lookup(cl);
            vector_type   sigma    = var.stresses.at_cell(msh, cl);
            matrix_type   varsigma = var.stresses.at_element_faces(msh, cl);

            auto cbar = barycenter(msh, cl);
            auto cphi = cb.eval_functions(msh, cl, cbar, 0, m_degree);
            matrix_type vec_phi = make_vectorial_matrix(cphi, DIM);
            vector_type sigma_eval = vec_phi * sigma;
            auto  sigma_norm  = sigma_eval.norm();

            if( sigma_norm  < pst.yield )
            {
                auto is_flowing(false);
                for(size_t col = 0; col < fcs.size(); col++)
                {
                    auto fc = fcs.at(col);
                    auto fc_varsigma = varsigma.col(col);
                    auto fbar = barycenter(msh, fc);
                    auto fphi = fb.eval_functions(msh, fc, fbar );

                    auto varsig_eval = fphi.dot(fc_varsigma);
                    auto varsig_norm = varsig_eval;

                    if(varsig_norm > pst.yield )
                        is_flowing = true;
                }
                if(!is_flowing)
                    initial_marks.at(cell_id) = 1;
            }
        }

        auto ret = mark_by_node_neighbors(msh, initial_marks);
        return ret;
    }

    template<typename CellBasisType, typename FaceBasisType, typename Solution, typename Variable>
    bool
    marker( const mesh_type& msh,
            const Variable& var,
            const plasticity_data<T>& pst,
            const size_t& m_degree,
            const Solution& solution,
            const std::vector<matrix_type>& grad_global)
    {
        bool ret;
        std::cout << "marker_name : " << mp.marker_name << std::endl;

        switch(mp.marker_name)
        {
            case 1:
                return ret = test_marker_ratio(msh, pst.Bn, 1);
            case 2:
                return ret = test_marker_all(msh, 1);
            case 3:
                return ret = test_marker_yield(msh, var.stresses, pst, m_degree);

                //return ret = marker_prueba(msh, var.stress_Th, m_degree);
            //case 2:
            //    return ret = marker_jb(msh, var.stress_Th);
            //case 4:
            //    return ret = marker_m4(msh, var.stress_Th);
            #if 0
            case 3:
            case 6:
                return ret = marker_ae<CellBasisType, FaceBasisType, Solution>
                        (msh, var.stress_Th, var.velocity_Th, m_degree, solution, grad_global);
            #if 0
            case 5:
            if(mp.diff)
            {

                // "AE_DIFFUSION"
                do_refinement = sbm.template marker_diffusion<cell_basis_type, cell_quadrature_type,
                                                     face_basis_type, face_quadrature_type,
                                                     gradrec_type, Solution>
                (msh, m_pst, dir_name, info, mp, var.velocity_Th, m_degree, imsh, solution, grad_global);
                break;
            }
            #endif
            #endif
            default:
                throw std::invalid_argument("No marker specification.");

        }
    }
    // *****************************  Refinement  ******************************
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
                  m_boundary_edges.push_back({a , b, t.bid.at(i)});
          }
          auto a = t.p.at(0);
          auto b = t.p.at(N-1);
         if(b < a)
             std::swap(a,b);
         m_edges.push_back({a , b});
         if(t.b.at(N -1))
             m_boundary_edges.push_back({a , b, t.bid.at(N -1)});
     }
    // ******************************  Refinement  ******************************
    void
    face_mark_hanging_nodes(mesh_type & msh, const face_type & fc)
    {
        auto storage = msh.backend_storage();
        auto  fc_id  = msh.lookup(fc);
        faces_marks.at(fc_id).first  = true;
        auto bar = barycenter(msh,fc);
        storage ->points.push_back(bar);
        faces_marks.at(fc_id).second = storage->points.size() -1;
        //std::cout << "     faces_mark.at("<<fc_id <<") = "<< faces_marks.at(fc_id).second << std::endl;
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
          nt.bid[j] = t.bid[i];

          if(fc_mark)
           {
               j++;
               nt.p[j] = faces_marks.at(fcs_ids[i]).second;
               nt.b[j] = t.b[i];
               nt.bid[j] = t.bid[i];
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
          auto bnd  = storage->boundary_info.at(id);
          t.bid[i]  = bnd.boundary_id;
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
           t.bid[i] = 10000;

           auto fid = fcs_ids.at(i);
           auto is_bound = storage->boundary_edges.at(fid);
           if(is_bound)
           {
               t.b[i]   = new_fc_is_boundary(msh, pts_ids, fid, N, i);
                if(t.b[i])
                {
                    auto  bnd = storage->boundary_info.at(fid);
                    t.bid[i] = bnd.boundary_id;
                }
                else
                    t.bid[i] = 10000;
           }
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
            auto binfo =  storage->boundary_info.at(id);
            t.bid[i] = binfo.boundary_id;
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

            t0.bid[0] = t.bid[0];   t0.bid[1] = 10000;      t0.bid[2] = t.bid[2];
            t1.bid[0] = t.bid[0];   t1.bid[1] = t.bid[1];   t1.bid[2] = 10000;
            t2.bid[0] = 10000;      t2.bid[1] = t.bid[1];   t2.bid[2] = t.bid[2];
            t3.bid[0] = 10000;      t3.bid[1] = 10000;      t3.bid[2] = 10000;

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

            t0.bid[0] = t.bid[0];   t0.bid[1] = t.bid[1];   t0.bid[2] = t.bid[2];
            t1.bid[0] = t.bid[1];   t1.bid[1] = t.bid[2];   t1.bid[2] = t.bid[0];


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
                    t0.bid[permut[i]] = 10000;      t1.bid[i] = 10000;
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

        t0.bid[0] = t.bid[0];   t0.bid[1] = 10000;      t0.bid[2] = t.bid[2];
        t1.bid[0] = t.bid[0];   t1.bid[1] = t.bid[1];   t1.bid[2] = 10000;
        t2.bid[0] = 10000;      t2.bid[1] = t.bid[1];   t2.bid[2] = t.bid[2];
        t3.bid[0] = 10000;      t3.bid[1] = 10000;      t3.bid[2] = 10000;

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

            auto bnd = storage->boundary_info.at(fcs_ids.at(i));
            nt.bid[0] = bnd.boundary_id;
            nt.bid[1] = 10000;
            nt.bid[2] = 10000;

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

            auto bnd = storage->boundary_info.at(fcs_ids.at(i));
            nt.bid[0] = bnd.boundary_id;
            nt.bid[1] = 10000;
            nt.bid[2] = 10000;

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

                auto bnd = storage->boundary_info.at(fcs_ids.at(v));
                nt.bid[0] = bnd.boundary_id;
                nt.bid[1] = 10000;
                nt.bid[2] = 10000;

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

                auto bnd1 = storage->boundary_info.at(fcs_ids.at(v));
                auto bnd2 = storage->boundary_info.at(fcs_ids.at((v + 1) % num_fcs));

                nt.bid[0] = bnd1.boundary_id;
                nt.bid[1] = bnd2.boundary_id;
                nt.bid[2] = 10000;
                nt.bid[3] = 10000;

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
                    t3.bid.at(ii) = 10000;
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

                std::cout << "  Outside triangle: " << std::endl;
                std::cout << "   * "<< nt << std::endl;

                store(nt,m_triangles);
                l_triangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 4)
            {
                auto nt = put_n_gon<4>(msh, fcs_ids, pts_ids);

                std::cout << "  Outside quadrangles: " << std::endl;
                std::cout << "   * "<< nt << std::endl;

                store(nt,m_quadrangles);
                l_quadrangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 5)
            {
                auto nt = put_n_gon<5>(msh, fcs_ids, pts_ids);

                //std::cout << "  Outside triangle: " << std::endl;
                //std::cout << "   * "<< nt << std::endl;

                store(nt,m_pentagons);
                l_pentagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 6)
            {
                auto nt = put_n_gon<6>(msh, fcs_ids, pts_ids);

                //std::cout << "  Outside triangle: " << std::endl;
                //std::cout << "   * "<< nt << std::endl;

                store(nt,m_hexagons);
                l_hexagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 7)
            {
                auto nt = put_n_gon<7>(msh, fcs_ids, pts_ids);

                //std::cout << "  Outside triangle: " << std::endl;
                //std::cout << "   * "<< nt << std::endl;

                store(nt,m_heptagons);
                l_heptagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }

            if(new_num_fcs == 8)
            {
                auto nt = put_n_gon<8>(msh, fcs_ids, pts_ids);

                //std::cout << "  Outside triangle: " << std::endl;
                //std::cout << "   * "<< nt << std::endl;

                store(nt,m_octagons);
                l_octagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }

            if(new_num_fcs == 9)
            {
                auto nt = put_n_gon<9>(msh, fcs_ids, pts_ids);
                store(nt,m_enneagons);
                l_enneagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }

            if(new_num_fcs == 10)
            {
                auto nt = put_n_gon<10>(msh, fcs_ids, pts_ids);
                store(nt,m_decagons);
                l_decagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 11)
            {
                auto nt = put_n_gon<11>(msh, fcs_ids, pts_ids);
                store(nt,m_hendecagons);
                l_hendecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }

            if(new_num_fcs == 12)
            {
                auto nt = put_n_gon<12>(msh, fcs_ids, pts_ids);
                store(nt,m_dodecagons);
                l_dodecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 13)
            {
                auto nt = put_n_gon<13>(msh, fcs_ids, pts_ids);
                store(nt,m_triadecagons);
                l_triadecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 14)
            {
                auto nt = put_n_gon<14>(msh, fcs_ids, pts_ids);
                store(nt,m_tesseradecagons);
                l_tesseradecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 15)
            {
                auto nt = put_n_gon<15>(msh, fcs_ids, pts_ids);
                store(nt,m_pentadecagons);
                l_pentadecagons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 16)
            {
                auto nt = put_n_gon<16>(msh, fcs_ids, pts_ids);
                store(nt,m_16gons);
                l_16gons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 17)
            {
                auto nt = put_n_gon<17>(msh, fcs_ids, pts_ids);
                store(nt,m_17gons);
                l_17gons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 18)
            {
                auto nt = put_n_gon<18>(msh, fcs_ids, pts_ids);
                store(nt,m_18gons);
                l_18gons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 19)
            {
                auto nt = put_n_gon<19>(msh, fcs_ids, pts_ids);
                store(nt,m_19gons);
                l_19gons.push_back(std::make_pair(cl_level, cl.get_id()));
            }

            if(new_num_fcs == 20)
            {
                auto nt = put_n_gon<20>(msh, fcs_ids, pts_ids);
                store(nt,m_20gons);
                l_20gons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 21)
            {
                auto nt = put_n_gon<21>(msh, fcs_ids, pts_ids);
                store(nt,m_21gons);
                l_21gons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 22)
            {
                auto nt = put_n_gon<22>(msh, fcs_ids, pts_ids);
                store(nt,m_22gons);
                l_22gons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 23)
            {
                auto nt = put_n_gon<23>(msh, fcs_ids, pts_ids);
                store(nt,m_23gons);
                l_23gons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 24)
            {
                auto nt = put_n_gon<24>(msh, fcs_ids, pts_ids);
                store(nt,m_24gons);
                l_24gons.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 25)
            {
                auto nt = put_n_gon<25>(msh, fcs_ids, pts_ids);
                store(nt,m_25gons);
                l_25gons.push_back(std::make_pair(cl_level, cl.get_id()));
            }


            if(new_num_fcs > 25)
            {
                std::cout << "number of faces = "<< new_num_fcs << std::endl;
                throw std::logic_error("number of faces exceeds maximum. At least more than 1 hanging node is allowed, triangle cann't have more than 5 faces.");
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
        std::cout << "1.1 cell_mark ("<< cl_id <<"-"<< cid<< ")= "<< cl_mark << std::endl;
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
                        std::cout << "* refine special_tri" << std::endl;
                        refine_special_triangle(msh, cl, vertices_pos,cl_level);
                        break;
                    //case 4:
                    //    refine_quadrangle();
                    //    break;
                    default:
                        std::cout << "* refine special_other" << std::endl;
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
                        std::cout << "* refine tri" << std::endl;
                        refine_triangle(msh, cl, cl_level);
                        break;
                    //case 4:
                    //    refine_quadrangle();
                    //    break;
                    default:
                        std::cout << "* refine other" << std::endl;
                        refine_other(msh ,cl, cl_level);
                        break;
                }

            }
        }
        else
        {
            std::cout << "* No refinement" << std::endl;

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
                std::cout << " ** add triangle (no ref)" << std::endl;
                std::cout << "  * "<< nt << std::endl;

                l_triangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 4)
            {
                auto nt = put_n_gon_with_hanging_node<4>(msh, cl);
                store(nt, m_quadrangles);
                std::cout << " ** add quadrangle (no ref)" << std::endl;
                std::cout << "  * "<< nt << std::endl;
                l_quadrangles.push_back(std::make_pair(cl_level, cl.get_id()));
            }
            if(new_num_fcs == 5)
            {

                auto nt = put_n_gon_with_hanging_node<5>(msh, cl);
                store(nt,m_pentagons);
                std::cout << " ** add pentagons (no ref)" << std::endl;
                std::cout << "  * "<< nt << std::endl;
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
        fvca5_save_tuples(ofs, loader.m_enneagons);
        fvca5_save_tuples(ofs, loader.m_decagons);
        fvca5_save_tuples(ofs, loader.m_hendecagons);
        fvca5_save_tuples(ofs, loader.m_dodecagons);
        fvca5_save_tuples(ofs, loader.m_triadecagons);
        fvca5_save_tuples(ofs, loader.m_tesseradecagons);
        fvca5_save_tuples(ofs, loader.m_pentadecagons);
        fvca5_save_tuples(ofs, loader.m_16gons);
        fvca5_save_tuples(ofs, loader.m_17gons);
        fvca5_save_tuples(ofs, loader.m_18gons);
        fvca5_save_tuples(ofs, loader.m_19gons);
        fvca5_save_tuples(ofs, loader.m_20gons);
        fvca5_save_tuples(ofs, loader.m_21gons);
        fvca5_save_tuples(ofs, loader.m_22gons);
        fvca5_save_tuples(ofs, loader.m_23gons);
        fvca5_save_tuples(ofs, loader.m_24gons);
        fvca5_save_tuples(ofs, loader.m_25gons);

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

            std::cout << "1. cell_mark ("<< cl_id <<") ="<< cl_mark << std::endl;

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
                        std::cout << " ** add triangle (cell_mark = 0)" << std::endl;
                        l_triangles.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 4:
                        store(put_n_gon<4>(msh,cl), m_quadrangles);
                        std::cout << " ** add quadrangle (cell_mark = 0)" << std::endl;
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
                    case 16:
                        store(put_n_gon<16>(msh,cl), m_16gons);
                        l_16gons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 17:
                        store(put_n_gon<17>(msh,cl), m_17gons);
                        l_17gons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 18:
                        store(put_n_gon<18>(msh,cl), m_18gons);
                        l_18gons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 19:
                        store(put_n_gon<19>(msh,cl), m_19gons);
                        l_19gons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 20:
                        store(put_n_gon<20>(msh,cl), m_20gons);
                        l_20gons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 21:
                        store(put_n_gon<21>(msh,cl), m_21gons);
                        l_21gons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 22:
                        store(put_n_gon<22>(msh,cl), m_22gons);
                        l_22gons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 23:
                        store(put_n_gon<23>(msh,cl), m_23gons);
                        l_23gons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 24:
                        store(put_n_gon<24>(msh,cl), m_24gons);
                        l_24gons.push_back(std::make_pair(cl_level, cl.get_id()));
                        break;
                    case 25:
                        store(put_n_gon<25>(msh,cl), m_25gons);
                        l_25gons.push_back(std::make_pair(cl_level, cl.get_id()));
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
        level_ancestor.insert(level_ancestor.end(), l_16gons.begin()  , l_16gons.end());
        level_ancestor.insert(level_ancestor.end(), l_17gons.begin()  , l_17gons.end());
        level_ancestor.insert(level_ancestor.end(), l_18gons.begin()  , l_18gons.end());
        level_ancestor.insert(level_ancestor.end(), l_19gons.begin()  , l_19gons.end());
        level_ancestor.insert(level_ancestor.end(), l_20gons.begin()  , l_20gons.end());
        level_ancestor.insert(level_ancestor.end(), l_21gons.begin()  , l_21gons.end());
        level_ancestor.insert(level_ancestor.end(), l_22gons.begin()  , l_22gons.end());
        level_ancestor.insert(level_ancestor.end(), l_23gons.begin()  , l_23gons.end());
        level_ancestor.insert(level_ancestor.end(), l_24gons.begin()  , l_24gons.end());
        level_ancestor.insert(level_ancestor.end(), l_25gons.begin()  , l_25gons.end());

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
        loader.m_16gons    = m_16gons;
        loader.m_17gons    = m_17gons;
        loader.m_18gons    = m_18gons;
        loader.m_19gons    = m_19gons;
        loader.m_20gons    = m_20gons;
        loader.m_21gons    = m_21gons;
        loader.m_22gons    = m_22gons;
        loader.m_23gons    = m_23gons;
        loader.m_24gons    = m_24gons;
        loader.m_25gons    = m_25gons;


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

};

template<typename T>
void
print(const std::vector<T> vec, const std::string name)
{
    std::cout << name << std::endl;

    for(auto& v: vec)
        std::cout << " "<< v;
    std::cout << std::endl;
}


struct greater
{
    template<typename T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};

template<typename T>
void
sort_uniq_greater(std::vector<T>& v)
{
    std::sort(v.begin(), v.end(), greater());
    auto uniq_iter = std::unique(v.begin(), v.end());
    v.erase(uniq_iter, v.end());
}

template<typename T>
void
sort_uniq(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    auto uniq_iter = std::unique(v.begin(), v.end());
    v.erase(uniq_iter, v.end());
}

template<typename PointType>
auto
to_angle(const PointType& p, const PointType& o)
{
  return atan2(p.y() - o.y(), p.x() - o.x());
}


template<typename T, typename Storage>
void
sort_cclockwise(const disk::mesh<T,2,Storage>& msh,
        std::vector<size_t>& cells2erase,
        std::vector<typename disk::mesh<T,2,Storage>::face::id_type>& fids_vec,
        std::vector<typename disk::mesh<T,2,Storage>::point_type::id_type>& pids_vec )
{
    typedef disk::mesh<T,2,Storage>     mesh_type;
    typedef std::pair< size_t, typename mesh_type::point_type>       pair_type;
    typedef typename mesh_type::point_type              point_type;
    typedef typename point_type::id_type                point_id_type;
    typedef typename mesh_type::node_type::id_type      node_id_type;
    typedef typename mesh_type::face::id_type           face_id_type;

    sort_uniq(pids_vec);
    sort_uniq(fids_vec);

    auto init_fid = fids_vec.at(0);
    //auto new_fids_vec = init_fid;
    auto init_fc  = *std::next(msh.faces_begin(), init_fid);
    auto new_pids_vec = init_fc.point_ids();

    auto left_id  = new_pids_vec[0];
    auto right_id = new_pids_vec[1];



    for(size_t i = 1; i < pids_vec.size()-1; i++)
    {
        //print(new_pids_vec, "* new_pids_vec : ");

        for(auto& fid : fids_vec )
        {
            auto fc  = *std::next(msh.faces_begin(), fid);
            auto fc_pids  = fc.point_ids();

            if( left_id != fc_pids[0] && right_id == fc_pids[1])
            {
                left_id  = right_id;
                right_id = fc_pids[0];
                new_pids_vec.push_back(right_id);
                break;
            }
            if( left_id != fc_pids[1] && right_id == fc_pids[0])
            {
                left_id  = right_id;
                right_id = fc_pids[1];
                new_pids_vec.push_back(right_id);
                break;
            }
        }
    }
    //print(new_pids_vec, "* new_pids_vec : ");

    auto new_pids_size = new_pids_vec.size();
    assert(new_pids_size == pids_vec.size());

    auto storage = msh.backend_storage();
    std::vector<point_type>  pts(pids_vec.size());

    size_t i = 0;
    for(auto& pid : new_pids_vec)
        pts.at(i++) = storage->points.at(size_t(pid));

    bool cw(false), ccw(false);

    for(size_t i = 0; i < new_pids_size; i++)
    {
        auto p0 = new_pids_vec.at(i);
        auto p1 = new_pids_vec.at((i + 1) % new_pids_size);

        auto is_found =  false;

        for(auto& cl_id: cells2erase)
        {
            auto cl  = *std::next(msh.cells_begin(), cl_id);
            auto cl_pids  = cl.point_ids();
            auto cl_pids_size = cl_pids.size();
            for(size_t j = 0; j < cl_pids_size; j++)
            {
                auto cp0 = cl_pids.at(j);
                auto cp1 = cl_pids.at((j + 1) % cl_pids_size);

                if( p0 == cp0 && p1 == cp1)
                {
                    ccw = true;
                    is_found = true;
                    break;
                }
                if( p0 == cp1 && p1 == cp0)
                {
                    cw = true;
                    is_found = true;
                    break;
                }
            }
            if(is_found)
                continue;
        }
        if(!is_found)
            throw  std::logic_error(" Points not found in any cell points");
    }
    if(cw && ccw)
        throw std::logic_error("Face are ordered in CW and CCW. Check all faces have same direction.");
    if(cw)
        std::reverse(new_pids_vec.begin(),new_pids_vec.end());

    size_t j(0);
    for(auto npid : new_pids_vec)
        pids_vec.at(j++) = point_id_type(npid);
    return;
}

template<typename T, typename Storage>
void
sort_by_polar_angle(const disk::mesh<T,2,Storage>& msh,
        std::vector<typename disk::mesh<T,2,Storage>::point_type::id_type>& vec)
{

    typedef disk::mesh<T,2,Storage>     mesh_type;
    typedef std::pair< size_t, typename mesh_type::point_type>       pair_type;
    typedef typename mesh_type::point_type               point_type;
    typedef typename point_type::id_type                point_id_type;

    sort_uniq(vec);
    auto storage = msh.backend_storage();
    std::vector<point_type>  pts(vec.size());

    size_t i = 0;
    for(auto pid : vec)
        pts.at(i++) = storage->points.at(size_t(pid));

    auto ptsnum = pts.size();

    //Compute barycenter
    point_type  bar = point_type({0.0, 0.0});
    T           area = 0.0;

    for (size_t i = 0; i < ptsnum; i++)
    {
        auto p0 = pts[i];
        auto p1 = pts[(i+1)%ptsnum];

        auto a = p0.x()*p1.y() - p1.x()*p0.y();

        bar = bar + (p0 + p1) * a;
        area += a;
    }

    area *= 0.5;

    auto h = bar/(6.0 * area);

    //sorting point index by polar angle
    std::sort(vec.begin(),vec.end(),
        [&](const point_id_type & va, const point_id_type & vb )
        {
            auto pta = storage->points.at(va);
            auto ptb = storage->points.at(vb);

            auto theta_a = to_angle(pta, h);
            auto theta_b = to_angle(ptb, h);

            return (theta_a < theta_b);
        }
    );
    return;
}

template<typename MeshType>
std::vector<std::list<size_t>>
connectivity_by_face(const MeshType& msh)
{
   std::vector<std::list<size_t>> node_to_face;
   node_to_face.resize( msh.points_size());
   for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
   {
       auto fc = *itor;
       auto pids = fc.point_ids();
       auto fid  = msh.lookup(fc);
       node_to_face.at(pids[0]).push_back(fid);
       node_to_face.at(pids[1]).push_back(fid);
   }
   for(auto& l: node_to_face)
   {
       l.sort();
       l.unique();
   }
   return node_to_face;
}

template<typename MeshType>
std::vector<std::list<size_t>>
connectivity_by_cell(const MeshType& msh)
{
    std::vector<std::list<size_t>> node_to_cell;
    node_to_cell.resize( msh.points_size());
    for (auto & cl : msh)
    {
        auto pids = cl.point_ids();
        auto cid  = msh.lookup(cl);
        for(auto  id : pids )
            node_to_cell.at(id).push_back(cid);
    }
    for(auto& l: node_to_cell)
    {
      l.sort();
      l.unique();
    }
    return node_to_cell;
}
template<typename MeshType>
std::vector<size_t>
nodes_state(const MeshType& msh,
                std::vector<size_t>& initial_marks)
{
    std::vector<size_t> ret(msh.points_size());

    auto mesh_node_to_cell = connectivity_by_cell(msh);

    size_t count = 0;
    for(auto node_to_cell :  mesh_node_to_cell)
    {
        size_t cell_mark = 0;
        for(auto cell_id: node_to_cell)
            cell_mark += initial_marks.at(cell_id);

        if(cell_mark == node_to_cell.size())
            ret.at(count) = 1; // Liquid
        else if(cell_mark > 0)
            ret.at(count) = 2; // Solid - Liquid
        count++;
    }
    return ret;
}

template<typename MeshType>
class binding_meshes
{
    typedef MeshType mesh_type;
    typedef typename mesh_type::scalar_type  scalar_type;

    typedef typename mesh_type::cell         cell_type;
    typedef typename mesh_type::face         face_type;
    typedef typename mesh_type::node_type    node_type;
    typedef typename mesh_type::point_type    point_type;

    typedef typename face_type::id_type      face_id_type;
    typedef typename cell_type::id_type      cell_id_type;
    typedef typename node_type::id_type      node_id_type;
    typedef typename point_type::id_type     point_id_type;

    typedef dynamic_vector<scalar_type>                 vector_type;
    typedef dynamic_matrix<scalar_type>                 matrix_type;

    std::vector<size_t>  m_levels;
    std::vector<size_t>  temp_levels;
    std::vector<size_t>  new_levels;
    std::vector<size_t>  new_ancestors;

    std::vector<size_t>  check_cells;

    struct bin_geometry
    {
        std::vector<size_t> nodes2erase;
        std::vector<size_t> faces2erase;
        std::vector<size_t> cells2erase;
        std::vector<size_t> cellslevels;

        std::vector<face_id_type> faces2bin;
        std::vector<cell_type>   cells2bin;


        bin_geometry(){};

        void print_all()
        {
            print(nodes2erase,"/ * nodes2erase */" );
            print(faces2erase,"/ * faces2erase */" );
            print(cells2erase,"/ * cells2erase */" );
            print(cells2erase,"/ * cellslevels */" );
            print(faces2bin,"/ * faces2bin */" );
        }
        friend bin_geometry operator+(const bin_geometry& bg1, const bin_geometry& bg2)
        {
            bin_geometry ret;
            ret = bg1;
            ret.nodes2erase.insert(ret.nodes2erase.begin(),
                                bg2.nodes2erase.begin(), bg2.nodes2erase.end());
            ret.faces2erase.insert(ret.faces2erase.begin(),
                                bg2.faces2erase.begin(), bg2.faces2erase.end());
            ret.cells2erase.insert(ret.cells2erase.begin(),
                                bg2.cells2erase.begin(), bg2.cells2erase.end());
            ret.cells2bin.insert(ret.cells2bin.begin(),
                                bg2.cells2bin.begin(), bg2.cells2bin.end());
            ret.cellslevels.insert(ret.cellslevels.begin(),
                                bg2.cellslevels.begin(), bg2.cellslevels.end());

            sort_uniq_greater(ret.nodes2erase);
            sort_uniq_greater(ret.faces2erase);
            sort_uniq_greater(ret.cells2erase);
            return ret;
        }

    };
private:
    bin_geometry bgeo_all;


   std::vector<size_t>
   mark_by_node_neighbors(const mesh_type& msh,
                            std::vector<size_t>& initial_marks)
   {
        std::vector<size_t> ret(msh.cells_size());
        std::vector<std::list<size_t>> node_to_face(msh.points_size());

        for(auto& cell : msh)
        {
            auto cl_id = msh.lookup(cell);
            auto cell_mark = initial_marks.at(cl_id);
            if(cell_mark == 1)
            {
                auto fcs =  faces(msh, cell);
                for(auto& fc: fcs)
                {
                    auto pids = fc.point_ids();
                    auto fid  = msh.lookup(fc);
                    node_to_face.at(pids[0]).push_back(fid);
                    node_to_face.at(pids[1]).push_back(fid);
                }
            }
        }
        for(auto& n2f:  node_to_face )
        {
            n2f.sort();
            n2f.unique();
        }
        auto mesh_node_to_face = connectivity_by_face(msh);

        for(auto& cell : msh)
        {
            auto cl_id = msh.lookup(cell);
            auto cell_mark = initial_marks.at(cl_id);
            if(cell_mark == 1)
            {
                bool all_ngh_true = true;
                auto pids = cell.point_ids();
                for(auto pid : pids)
                {
                    auto l1 =  node_to_face.at(pid);
                    auto l2 =  mesh_node_to_face.at(pid);
                    if(! std::equal (l1.begin (), l1.end (), l2.begin ()))  //if( !(l1 == l2) );
                        all_ngh_true = false;
                }
                if(all_ngh_true)
                    ret.at(cl_id) = 1;
            }
        }
        return ret;
   }


    std::vector<size_t>
    mark_by_face_neighbors(const mesh_type& msh, std::vector<size_t>& initial_marks)
    {
        std::vector<size_t> ret(msh.cells_size());

        for(auto& cell : msh)
        {
            auto cl_id = msh.lookup(cell);
            auto cell_mark = initial_marks.at(cl_id);
            if(cell_mark == 1)
            {
                auto all_ngh_true = true;
                auto fcs =  faces(msh, cell);
                for(auto& fc: fcs)
                {
                    auto  ngh_id   = msh.neighbor_id(cell, fc);
                    auto  ngh_mark = initial_marks.at(ngh_id);
                    if(ngh_mark != 1)
                        all_ngh_true = false;
                }
                if(all_ngh_true)
                    ret.at(cl_id) = 1;
            }
        }
        return ret;
    }
public:
    template<typename VariableType>
    auto
    marker(const mesh_type& msh,
            const VariableType & var,
            const disk::plasticity_data<scalar_type>& pst,
            const size_t m_degree)
    {
        typedef disk::quadrature<mesh_type, cell_type>     cell_quadrature_type;
        typedef disk::quadrature<mesh_type, face_type>     face_quadrature_type;

        typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
        typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

        cell_basis_type         cb(m_degree);
        face_basis_type         fb(m_degree);
        cell_quadrature_type    cq(2 * m_degree);
        face_quadrature_type    fq(2 * m_degree);

        auto DIM = mesh_type::dimension;
        std::vector<size_t> initial_marks(msh.cells_size());
        std::vector<size_t> ret(msh.cells_size());

        size_t num_cells_marked = 0;
        for(auto & cl : msh)
        {
            auto fcs = faces(msh, cl);
            auto cell_id = msh.lookup(cl);
            vector_type   sigma    = var.stresses.at_cell(msh, cl);

            auto cbar = barycenter(msh, cl);
            auto cphi = cb.eval_functions(msh, cl, cbar, 0, m_degree);
            matrix_type vec_phi = make_vectorial_matrix(cphi, DIM);
            vector_type sigma_eval = vec_phi * sigma;
            auto  sigma_norm  = sigma_eval.norm();

            if( sigma_norm > pst.yield )
                initial_marks.at(cell_id) = 1;
        }


        dump_to_matlab(msh, "viscoplasticity/circular/Meeting_12_05/k_"+ tostr(m_degree)+ "/BinningAdaptive/borrar_initial_marks.m", initial_marks);
        auto node_state_marks = nodes_state(msh, initial_marks);

        for(auto cl : msh)
        {
            auto cid = msh.lookup(cl);
            auto nodes = cl.point_ids();

            size_t sum_node_marks =  0;

            for(auto n: nodes)
                sum_node_marks += node_state_marks.at(n);

            if( sum_node_marks == 0)
                ret.at(cid) = 3;
        }

        dump_to_matlab(msh, "viscoplasticity/circular/Meeting_12_05/k_"+ tostr(m_degree)+ "/BinningAdaptive/mesh_bin_marks.m", ret);

        return ret;
    }

    void
    search4faces(const mesh_type& msh, const cell_type& cell, bin_geometry & bgeo,
                const std::vector<size_t>& cells_marks)
    {
        auto fcs = faces(msh, cell);
        auto cid =  msh.lookup(cell);

        check_cells.at(cid) = 1;
        bgeo.cells2erase.push_back(cid);

        for(auto& face: fcs)
        {
            auto fid = msh.lookup(face);

            if(!msh.is_boundary(face))
            {
                auto ngh = msh.neighbor(cell, face);
                auto nid = msh.lookup(ngh);
                auto is_checked = check_cells.at(nid) == 1? true: false;
                if(cells_marks.at(nid))
                {
                    check_cells.at(nid) = 1;

                    std::cout << "push_back ("<< fid <<")" << std::endl;
                    bgeo.faces2erase.push_back(fid);

                    if(!is_checked)
                        search4faces(msh, ngh, bgeo, cells_marks);
                }
                else
                    bgeo.faces2bin.push_back(fid);
            }
            else
                bgeo.faces2bin.push_back(fid);
        }
        return;
    }
    auto
    search4nodes(const mesh_type&  msh, bin_geometry& bgeo)
    {
        std::vector<size_t> nodes2erase;

        for(auto& fid : bgeo.faces2erase)
        {
            auto fc   = *std::next(msh.faces_begin(), fid);
            auto pids = fc.point_ids();
            //std::cout << "Face : "<< fid << std::endl;
            for(auto& p : pids)
            {
                auto is_found = false;
                //std::cout << "* Punto: "<< p << std::endl;
                for(auto& poly_fid : bgeo.faces2bin)
                {
                    auto poly_fc  = *std::next(msh.faces_begin(), poly_fid);
                    auto poly_pids = poly_fc.point_ids();

                    if(p == poly_pids[0] || p == poly_pids[1])
                    {
                        is_found = true;
                        break;
                    }
                }
                if(!is_found)
                    nodes2erase.push_back(p);
            }
        }
        sort_uniq_greater(nodes2erase);
        bgeo.nodes2erase = nodes2erase;
        return;
    }
    public:
    binding_meshes(const mesh_type& msh, const std::vector<size_t>& levels,
                                            const std::vector<size_t>& ancestors)
    {
        m_levels  = levels;
        new_ancestors   = ancestors;
        check_cells = std::vector<size_t>(msh.cells_size());
        dump_to_matlab(msh, "VP_initial_mesh.m");
    };
    binding_meshes(const mesh_type& msh)
    {
        temp_levels  = std::vector<size_t>(msh.cells_size());
        new_levels  = std::vector<size_t>(msh.cells_size());
        check_cells = std::vector<size_t>(msh.cells_size());

        dump_to_matlab(msh, "VP_initial_mesh.m");
    };

    void
    find_geometry(const mesh_type& msh, const std::vector<size_t>& cells_marks)
    {
        std::vector<bin_geometry> bgeo_temp;

        for(auto& cell : msh)
        {
            auto cell_id = msh.lookup(cell);
            bool is_checked = check_cells.at(cell_id) == 1? true : false;

            if(!is_checked && cells_marks.at(cell_id))
            {
                check_cells.at(cell_id) = 1;
                bin_geometry bgeo;
                search4faces(msh, cell, bgeo, cells_marks);

                sort_uniq_greater(bgeo.faces2erase);
                sort_uniq_greater(bgeo.faces2bin);
                sort_uniq_greater(bgeo.cells2erase);

                size_t new_cell_level = 0;
                for(auto& cl_id : bgeo.cells2erase)
                    new_cell_level = std::max(m_levels.at(cl_id), new_cell_level);

                search4nodes(msh, bgeo);

                auto f2b = bgeo.faces2bin;

                sort_uniq(f2b);
                bgeo.cells2bin.push_back(cell_type(f2b));
                bgeo.cellslevels.push_back(new_cell_level);

                bgeo.print_all();

                bgeo_temp.push_back(bgeo);

                #if 0
                // Just to check the points to erase are the right ones. Erase this after.
                auto ptid_to_point = [&](const size_t& pi) -> auto {
                    auto itor = msh.points_begin();
                    std::advance(itor, pi);
                    return *itor;
                };
                std::vector<typename mesh_type::point_type> poly_pts(nodes2erase.size());
                std::transform(nodes2erase.begin(), nodes2erase.end(), poly_pts.begin(), ptid_to_point);

                size_t i = 0;
                for(auto& pid: poly_pts)
                    std::cout << nodes2erase.at(i++)<<":  "<< pid.x() << "  "<< pid.y()<< std::endl;
                std::cout << std::endl;
                #endif
            }
        }

        for(auto&  bg_geo : bgeo_temp)
            bgeo_all = bgeo_all + bg_geo;
        return;
    }

    void
    erase_geometry(const mesh_type& msh, mesh_type& new_msh)
    {
        auto storage = msh.backend_storage();
        auto new_storage = new_msh.backend_storage();

        //Erase nodes and points
        for(auto nid: bgeo_all.nodes2erase)
        {
            new_storage->nodes.erase(new_storage ->nodes.begin() + nid);
            new_storage->points.erase(new_storage ->points.begin() + nid);
        }
        //Erase faces
        for(auto fid: bgeo_all.faces2erase)
        {
            new_storage->edges.erase(new_storage->edges.begin() + fid);
            new_storage->boundary_info.erase(new_storage->boundary_info.begin() + fid);
            new_storage->boundary_edges.erase(new_storage->boundary_edges.begin() + fid);
        }
        //Erase cells
        for(auto cid: bgeo_all.cells2erase)
        {
            new_storage->surfaces.erase(new_storage->surfaces.begin() + cid);
            new_storage->special_surfaces.erase(new_storage->special_surfaces.begin() + cid);
            new_levels.erase(new_levels.begin() + cid);
        }

        return;
    }

    void
    copy(const mesh_type& msh, mesh_type& new_msh)
    {
        auto storage = msh.backend_storage();
        auto new_storage = new_msh.backend_storage();

        new_storage->nodes    = storage->nodes;
        new_storage->points   = storage->points;
        new_storage->edges    = storage->edges;
        new_storage->surfaces = storage->surfaces;
        new_storage->boundary_info  = storage->boundary_info;
        new_storage->boundary_edges = storage->boundary_edges;
        new_storage->special_surfaces = storage->special_surfaces;
    }
    std::pair<bool, std::vector<point_id_type>>
    bin_is_special_polygon(const mesh_type& old_msh, const cell_type& cl)
    {
        std::cout << "INSIDE bin_is_special_polygon" << std::endl;

        auto storage = old_msh.backend_storage();
        auto pids  = cl.point_ids();
        auto fids  = cl.faces_ids();

        auto num_pts = pids.size();
        std::vector<static_vector<scalar_type,2>> ns(fids.size());
        size_t i = 0;
        for(size_t i =0; i < pids.size(); i++)
        {
            auto p1 = storage->points.at(pids.at(i));
            auto p2 = storage->points.at(pids.at(i%pids.size()));

            auto v = p2 - p1;
            ns.at(i++) = (point<scalar_type,2>({-v.y(), v.x()})).to_vector();
        }

        size_t hanging_count = 0;
        for (size_t i = 0; i < num_pts; i++)
        {
            auto nl = ns.at(i);
            auto nr = ns.at(( i + 1 )% num_pts);

            if( std::abs(nl.dot(nr) - scalar_type(1)) < 1.e-5)
                hanging_count++;
        }
        if(hanging_count > pids.size())
            throw std::logic_error(" More than one hanging node by face. Check correct counting of hanging nodes.");
        size_t num_vertices = pids.size() - hanging_count;
        if(num_vertices <= 0)
        {
            std::cout << "num_vertices : "<< num_vertices << std::endl;
            std::cout << "num_points   : "<< pids.size() << std::endl;
            std::cout << "hanging_nodes: "<< hanging_count << std::endl;

            throw std::logic_error(" Number of vertices are <= 0. Check correct counting of hanging nodes.");
        }

         //Identify vertices
         std::vector<point_id_type> vertices( num_vertices);

         size_t vcount = 0;
         for(size_t i = 0; i < num_pts; i++)
         {
             auto nl = ns.at(i);
             auto nr = ns.at((i + 1) %num_pts);

            if( std::abs( nr.dot(nl)- scalar_type(1)) >= 1.e-5)
                vertices.at(vcount++) =  point_id_type((i+1)%num_pts);
        }

        if(vcount != num_vertices)
        {
            std::cout << "vcount "<< vcount << std::endl;
            throw  std::logic_error(" Incorrect procedure to find vertices");
        }

        bool has_hang_nodes(false);
        if(vertices.size() != pids.size())
            has_hang_nodes = true;

        return std::make_pair(has_hang_nodes, vertices);
    }

    void
    add_new_cells(const mesh_type& msh, mesh_type& new_msh)
    {
        auto new_storage = new_msh.backend_storage();

        std::cout << "/ * cells2bin */" << std::endl;
        size_t count_level = 0;
        for(auto& cl: bgeo_all.cells2bin)
        {
            auto fids = cl.faces_ids();

            print(fids, " faces" );

            std::vector<typename face_type::id_type> cell_faces(fids.size());
            std::vector<point_id_type> cell_pids(fids.size() * 2);

            size_t i(0), ii(0);
            for(auto fid : fids)
            {
                auto f = *next(msh.faces_begin(), fid);
                auto pids = f.point_ids();

                std::cout << " * face "<< msh.lookup(f)<<":"<< pids[0]<< "  " << pids[1] << std::endl;

                for( auto fp : f.point_ids())
                    cell_pids.at(ii++) = point_id_type(fp);
            }

            //sort_by_polar_angle(msh, cell_pids);
            sort_cclockwise(msh, bgeo_all.cells2erase, fids, cell_pids);

            print(cell_pids, "* cell_pids_after_sort:  :");

            std::vector<face_id_type> new_fids(fids.size());

            size_t count(0);
            for(size_t i = 0; i < fids.size(); i++)
            {
                auto lp = cell_pids.at(i);
                auto rp = cell_pids.at((i+1)%cell_pids.size());

                for(auto fid : fids)
                {
                    std::cout << "face : "<<  fid<< std::endl;
                    auto fc = *next(msh.faces_begin(), fid);
                    auto pids = fc.point_ids();

                    auto left_it  = std::find(pids.begin(), pids.end(), lp);
                    auto right_it = std::find(pids.begin(), pids.end(), rp);

                    if (left_it != pids.end() &&  right_it != pids.end())
                    {
                        new_fids.at(count++) = fid;
                        break;
                    }
                }
            }
            assert(count == fids.size());

            auto cell = cell_type(new_fids);
            cell.set_point_ids(cell_pids.begin(), cell_pids.end());

            print(cell.point_ids(), "cell.pids: ");

            new_storage->surfaces.push_back(cell);
            new_levels.push_back(bgeo_all.cellslevels.at(count_level++));

            auto p = bin_is_special_polygon(msh, cell);
            new_storage->special_surfaces.push_back(p);

        }
        std::cout << "!!! *********** Before Sorting NEW MESH ERASE *********** !!! " << std::endl;
        for(auto cl : new_storage->surfaces)
        {
            auto fids = cl.faces_ids();

            std::cout << "cell :";// << new_msh.lookup(cl)<<":";
            for(auto fid : fids )
                std::cout<<"  "<< fid;

            std::cout << std::endl;
        }
        for(auto cl : new_storage->surfaces)
        {
            auto pids = cl.point_ids();

            std::cout << "cell :";// << new_msh.lookup(cl)<<":";
            for(auto pid : pids )
                std::cout<<"  "<< pid;

            std::cout << std::endl;
        }


        std::vector<int> index(new_storage->surfaces.size(), 0);
        for(int i = 0 ; i != index.size() ; i++)
            index.at(i) = i;


        std::sort(index.begin(), index.end(),[&](const int& a, const int& b)
            {
                cell_type s1  = new_storage->surfaces.at(a);
                cell_type s2  = new_storage->surfaces.at(b);
                return (s1 < s2);
            }
        );

        std::sort(new_storage->surfaces.begin(), new_storage->surfaces.end());

        typedef std::vector<std::pair<bool, std::vector<point_id_type>>> ss_type;

        auto new_special_surfaces = ss_type(new_storage->surfaces.size());

        for(size_t i = 0; i < new_storage->surfaces.size(); i++)
        {
            auto idx = index.at(i);
            new_special_surfaces.at(i) = new_storage->special_surfaces.at(idx);
            new_levels.at(i) = temp_levels.at(idx);
        }

        new_storage->special_surfaces = new_special_surfaces;

        std::cout << "!!! *********** After Sorting NEW MESH ERASE *********** !!! " << std::endl;
        for(auto cl : new_storage->surfaces)
        {
            auto fids = cl.faces_ids();

            std::cout << "cell :";// << new_msh.lookup(cl)<<":";
            for(auto fid : fids )
                std::cout<<"  "<< fid;
            std::cout << std::endl;
        }

        return;
    }

    #if 0
    mesh_type
    bin(const mesh_type& msh)
    {

        //Copy
        mesh_type new_msh;
        auto storage = msh.backend_storage();
        auto new_storage = new_msh.backend_storage();

        copy(msh, new_msh);
        erase_geometry(msh, new_msh);
        add_new_cells(msh, new_msh);

        std::cout << "/* message */" << std::endl;
        std::cout << "!!! _____ OLD MESH _____ !!! " << std::endl;
        for(auto& cl : msh)
        {
            std::cout << "cell : "<< cl.get_id() << " vs "<< msh.lookup(cl)<< std::endl;

            auto fcs =  faces(msh, cl);
            std::vector<typename face_type::id_type> cell_faces(fcs.size());
            size_t j(0), ii(0);
            for(auto f : fcs)
            {
                auto pids = f.point_ids();
                auto msp  = f.get_m_sids_ptrs();
                std::cout << " * face "<< msh.lookup(f)<<":"<< pids[0]<< "  " << pids[1] << std::endl;
                std::cout << " * msp  "<< msh.lookup(f)<<":"<< msp[0]<< "  " << msp[1] << std::endl;

                auto fid = find_element_id(msh.faces_begin(), msh.faces_end(), f);

                if(!fid.first)
                    throw std::invalid_argument("Element not found");
                else
                    cell_faces.at(j++) = fid.second;


            }
            print(cell_faces, "* old_faces :");
        }
        std::cout << "!!! _____ NEW MESH ERASE_____ !!! " << std::endl;
        size_t nn =0;
        std::cout << "Erase Cells" << std::endl;
        for(auto cl : new_storage->surfaces)
        {
            auto msp  = cl.get_m_sids_ptrs();
            auto fids = cl.faces_ids();

            std::cout << "cell :";// << new_msh.lookup(cl)<<":";
            for(auto fid : fids )
                std::cout<<"  "<< fid;
            std::cout << std::endl;
        }
        std::cout << "!!! _____ NEW MESH ERASE_____ !!! " << std::endl;

        for(auto& cl : new_storage->surfaces)
        {
            std::cout << "cell : "<<new_msh.lookup(cl)<< std::endl;  //cl.get_id() << " vs "<< new_msh.lookup(cl)<< std::endl;

            auto fids = cl.faces_ids();
            std::vector<typename face_type::id_type> cell_faces(fids.size());
            size_t j(0), ii(0);

            for(auto fid : fids)
            {
                auto f = *next(msh.faces_begin(), fid);
                auto pids = f.point_ids();
                auto msp  = f.get_m_sids_ptrs();
                std::cout << " * face "<< msh.lookup(f)<<":"<< pids[0]<< "  " << pids[1] << std::endl;
                std::cout << " * msp  "<< msh.lookup(f)<<":"<< msp[0]<< "  " << msp[1] << std::endl;
                auto new_fid = find_element_id(new_storage->edges.begin(),
                                                new_storage->edges.end(), f);

                if(!new_fid.first)
                    throw std::invalid_argument("Element not found");
                else
                    cell_faces.at(j++) = new_fid.second;
            }
            print(cell_faces, "* new_faces :");
        }

        //change numeration
        auto renumbering_faces = [&](const face_type & fc) -> auto
        {
            size_t i(0),j(0);
            auto nids = fc.point_ids();
            std::vector<typename node_type::id_type> face_nodes(nids.size());
            std::vector<typename point_type::id_type> face_points_ids(nids.size());

            for(auto & nid : nids)
            {
                std::cout << "node: "<< nid << std::endl;
                auto new_nid = find_element_id( new_storage ->nodes.begin(),
                                                    new_storage ->nodes.end(),
                                                    node_type(nid));
                if(!new_nid.first)
                    throw std::invalid_argument("Element not found");
                else
                {
                    face_nodes.at(i++) = node_id_type(new_nid.second);
                    face_points_ids.at(j++) = point_id_type(new_nid.second);
                }
            }

            auto face = face_type(face_nodes);
            face.set_point_ids(face_points_ids.begin(), face_points_ids.end());
            return face;
        };

        auto renumbering_cells = [&](const cell_type & cl) -> auto
        {

            //std::cout << "cell : "<< cl.get_id() << " vs "<< new_msh.lookup(cl)<< std::endl;
            auto fids = cl.faces_ids();
            auto pids = cl.point_ids();
            std::vector<typename face_type::id_type> cell_faces(fids.size());
            std::vector<point_id_type> cell_pids(fids.size() * 2);

            size_t i(0), ii(0);
            for(auto fid : fids)
            {
                auto f = *next(msh.faces_begin(), fid);
                auto pids = f.point_ids();
                auto msp  = f.get_m_sids_ptrs();
                std::cout << " * face "<< msh.lookup(f)<<":"<< pids[0]<< "  " << pids[1] << std::endl;
                std::cout << " * msp  "<< msh.lookup(f)<<":"<< msp[0]<< "  " << msp[1] << std::endl;
                auto new_fid = find_element_id(new_storage->edges.begin(),
                                                new_storage->edges.end(), f);
                if(!new_fid.first)
                    throw std::invalid_argument("Element not found");
                else
                    cell_faces.at(i++) = new_fid.second;

                for( auto fp : f.point_ids())
                    cell_pids.at(ii++) = point_id_type(fp);

            }

            //print(cell_faces, "* new_faces :");
            auto cell = cell_type(cell_faces);

            cell.set_point_ids(pids.begin(), pids.end()); //just copy all pids

            return cell;
        };

        std::transform(new_msh.cells_begin(), new_msh.cells_end(),
                            new_msh.cells_begin(), renumbering_cells);

        size_t count_id(0);

        auto change_id = [&](const cell_type & cl) -> auto
        {

            typedef typename cell_type::id_type cell_id_type;
            cell_type ret = cl;
            ret.set_element_id(cell_id_type(count_id++));

            return ret;
        };

        std::transform(new_msh.cells_begin(), new_msh.cells_end(),
                                new_msh.cells_begin(), change_id);

        nn =0;
        std::cout << "!!! ____ Cells transformation ____!!!" << std::endl;
        for(auto & cl : new_storage->surfaces)
        {
            std::cout << "cell  "<< new_msh.lookup(cl)<<":" ;//<< std::endl;

            auto fids = cl.faces_ids();
                //std::cout << "cells "<< n<<":";
            for(auto fid : fids )
                std::cout<<"  "<< fid;
            std::cout << std::endl;
        }

        std::transform(new_msh.faces_begin(), new_msh.faces_end(),
                        new_msh.faces_begin(), renumbering_faces );

        nn =0;
        std::cout << "!!! ____ Faces transformation ____!!!" << std::endl;
        for(auto & f : new_storage->edges)
        {
            auto pids = f.point_ids();
            std::cout << "face "<< nn++<<":"<< pids[0]<< "  " << pids[1] << std::endl;
        }


        auto renumbering_cell_nodes = [&](const cell_type & cl) -> auto
        {
            std::cout << "cell : "<< cl.get_id() << " vs "<< new_msh.lookup(cl)<< std::endl;

            auto pids = cl.point_ids();

            print(pids, "  * pids : ");
            std::vector<point_id_type> new_pids(pids.size());

            size_t i(0);
            for(auto pid : pids)
            {
                std::cout << " * pid : "<< pid << std::endl;
                auto itor = std::lower_bound(new_storage->nodes.begin(),
                                                new_storage->nodes.end(), node_type(pid));

                if (itor != new_storage->nodes.end() && !(node_type(pid) < *itor))
                {
                    auto idx = point_id_type(std::distance(new_storage->nodes.begin(), itor));
                    //auto node = *std::next(new_storage->nodes.begin(), idx);
                    std::cout << "** found in "<< idx << std::endl;

                    new_pids.at(i++)  = idx; //point_id_type(node);
                }
                else
                    throw std::invalid_argument("This is a bug: node not found");
            }

            assert(i == pids.size());
            #if 0
            size_t i(0);
            for(auto fid : fids)
            {
                auto f = *next(new_msh.faces_begin(), fid);
                auto pids = f.point_ids();
                //auto msp  = f.get_m_sids_ptrs();
                std::cout << " * face "<< new_msh.lookup(f)<<":"<< pids[0]<< "  " << pids[1] << std::endl;

                for( auto fp : f.point_ids())
                    cell_pids.at(i++) = point_id_type(fp);
            }
            #endif
            cell_type ret = cl;
            ret.set_point_ids(new_pids.begin(), new_pids.end());

            return ret;
        };

        std::transform(new_msh.cells_begin(), new_msh.cells_end(),
                            new_msh.cells_begin(), renumbering_cell_nodes);

        /* Nodes */
        auto new_nodes_size = new_storage->nodes.size();
        std::cout << "new_nodes_size"<< new_nodes_size << std::endl;
        std::vector<node_type> nodes(new_nodes_size);
        for (size_t i = 0; i < new_nodes_size; i++)
            nodes[i] = node_type(point_identifier<2>(i));

        new_storage->nodes = std::move(nodes);

        nn =0;
        std::cout << "!!! _______ Final Cells _______!!!" << std::endl;

        //for(auto cl: new_msh)
        for(auto & cl : new_storage->surfaces)
        {
            std::cout << "cell : "<< cl.get_id()<< std::endl;//  new_msh.lookup(cl) << std::endl;

            auto fids = cl.faces_ids();
            auto pids = cl.point_ids();

            print(fids, " * faces :");
            print(pids, " * nodes :");

        }
        dump_to_matlab(new_msh, "prueba_new_mesh.m");

        size_t iii(0);
        std::cout << "boundary_edges: " << std::endl;
        //for(auto cl: new_msh)
        for(bool is_boundary : new_storage->boundary_edges)
        {
            std::cout << "face ("<<iii++ <<"): "<< is_boundary << std::endl;//  new_msh.lookup(cl) << std::endl;
        }

        std::cout << "!!!!! Nodes !!!!! " << std::endl;
        //for(auto cl: new_msh)
        for(auto& l : new_storage->nodes)
        {
            auto pids = l.point_ids();
            std::cout << pids[0] << std::endl;//  new_msh.lookup(cl) << std::endl;
        }

        return new_msh;
    }
    #endif


        mesh_type
        bin(const mesh_type& msh)
        {
            //Copy
            mesh_type new_msh;
            auto storage = msh.backend_storage();
            auto new_storage = new_msh.backend_storage();

            copy(msh, new_msh);
            erase_geometry(msh, new_msh);
            add_new_cells(msh, new_msh);

            auto renumbering_cells = [&](const cell_type & cl) -> auto
            {
                auto fids = cl.faces_ids();
                auto pids = cl.point_ids();
                std::vector<typename face_type::id_type> cell_faces(fids.size());

                size_t i(0);
                for(auto fid : fids)
                {
                    auto f = *next(msh.faces_begin(), fid);
                    auto pids = f.point_ids();
                    auto new_fid = find_element_id(new_storage->edges.begin(),
                                                    new_storage->edges.end(), f);
                    std::cout << " * face "<< msh.lookup(f)<<":"<< pids[0]<< "  " << pids[1] << std::endl;

                    if(!new_fid.first)
                        throw std::invalid_argument("Element not found");
                    else
                        cell_faces.at(i++) = new_fid.second;
                }
                print(cell_faces, "* new_faces :");
                print(pids, "* old_pids :");
                auto cell = cell_type(cell_faces);

                cell.set_point_ids(pids.begin(), pids.end()); //just copy all pids

                return cell;
            };

            std::transform(new_msh.cells_begin(), new_msh.cells_end(),
                                new_msh.cells_begin(), renumbering_cells);

            size_t count_id(0);
            auto change_id = [&](const cell_type & cl) -> auto
            {
                typedef typename cell_type::id_type cell_id_type;
                cell_type ret = cl;
                ret.set_element_id(cell_id_type(count_id++));
                return ret;
            };

            std::transform(new_msh.cells_begin(), new_msh.cells_end(),
                                    new_msh.cells_begin(), change_id);

            //change numeration
            auto renumbering_faces = [&](const face_type & fc) -> auto
            {
                size_t i(0),j(0);
                auto nids = fc.point_ids();
                std::vector<typename node_type::id_type> face_nodes(nids.size());
                std::vector<typename point_type::id_type> face_points_ids(nids.size());

                for(auto & nid : nids)
                {
                    auto new_nid = find_element_id( new_storage ->nodes.begin(),
                                                    new_storage ->nodes.end(),
                                                                    node_type(nid));
                    if(!new_nid.first)
                        throw std::invalid_argument("Element not found");
                    else
                    {
                        face_nodes.at(i++) = node_id_type(new_nid.second);
                        face_points_ids.at(j++) = point_id_type(new_nid.second);
                    }
                }

                auto face = face_type(face_nodes);
                face.set_point_ids(face_points_ids.begin(), face_points_ids.end());
                return face;
            };

            std::transform(new_msh.faces_begin(), new_msh.faces_end(),
                            new_msh.faces_begin(), renumbering_faces );

            std::cout << "!!! _______ before cell nodes transformation _______!!!" << std::endl;

            //for(auto cl: new_msh)
            for(auto & cl : new_storage->surfaces)
            {
                std::cout << "cell : "<< cl.get_id()<< std::endl;//  new_msh.lookup(cl) << std::endl;

                auto fids = cl.faces_ids();
                auto pids = cl.point_ids();

                print(fids, " * faces :");
                print(pids, " * nodes :");
            }
            auto renumbering_cell_nodes = [&](const cell_type & cl) -> auto
            {
                std::cout << "cell : "<< cl.get_id() << " vs "<< new_msh.lookup(cl)<< std::endl;

                auto pids = cl.point_ids();

                print(pids, "  * pids : ");

                std::vector<point_id_type> new_pids(pids.size());

                size_t i(0);
                for(auto pid : pids)
                {
                    std::cout << " * pid : "<< pid;

                    auto itor = std::lower_bound(new_storage->nodes.begin(),
                                            new_storage->nodes.end(), node_type(pid));

                    if (itor != new_storage->nodes.end() && !(node_type(pid) < *itor))
                    {
                        auto idx = point_id_type(std::distance(new_storage->nodes.begin(), itor));
                        std::cout << " -> "<< idx << std::endl;
                        new_pids.at(i++)  = idx;
                    }
                    else
                        throw std::invalid_argument("This is a bug: node not found");
                }
                assert(i == pids.size());
                cell_type ret = cl;
                ret.set_point_ids(new_pids.begin(), new_pids.end());

                return ret;
            };

            std::transform(new_msh.cells_begin(), new_msh.cells_end(),
                                new_msh.cells_begin(), renumbering_cell_nodes);

            /* Nodes */
            auto new_nodes_size = new_storage->nodes.size();
            std::vector<node_type> nodes(new_nodes_size);
            for (size_t i = 0; i < new_nodes_size; i++)
                nodes[i] = node_type(point_identifier<2>(i));

            new_storage->nodes = std::move(nodes);

            std::cout << "!!! _______ Final Cells _______!!!" << std::endl;

            //for(auto cl: new_msh)
            for(auto & cl : new_storage->surfaces)
            {
                std::cout << "cell : "<< cl.get_id()<< std::endl;//  new_msh.lookup(cl) << std::endl;

                auto fids = cl.faces_ids();
                auto pids = cl.point_ids();

                print(fids, " * faces :");
                print(pids, " * nodes :");

            }

            size_t iii(0);
            std::cout << "boundary_edges: " << std::endl;

            for(bool is_boundary : new_storage->boundary_edges)
                std::cout << "face ("<<iii++ <<"): "<< is_boundary << std::endl;

            check_older_msh(new_msh);
            dump_to_matlab(new_msh, "VP_new_mesh.m");

            return new_msh;
        }

};
}//disk
