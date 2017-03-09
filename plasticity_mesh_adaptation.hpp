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
           const T yield)
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
          const T yield,
          const T percent)

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
            size_t id  = msh.lookup(cl);
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

    bool m_hanging_nodes;
    std::vector<size_t> cells_marks;
    std::vector<size_t> levels;
    std::vector<size_t> temp_levels;
    std::vector<size_t> new_levels;
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

    std::vector<size_t>        l_triangles;
    std::vector<size_t>        l_quadrangles;
    std::vector<size_t>        l_pentagons;
    std::vector<size_t>        l_hexagons;
    std::vector<size_t>        l_heptagons;
    std::vector<size_t>        l_octagons;
    std::vector<size_t>        l_enneagons;
    std::vector<size_t>        l_decagons;
    std::vector<size_t>        l_hendecagons;
    std::vector<size_t>        l_dodecagons;
    std::vector<size_t>        l_triadecagons;
    std::vector<size_t>        l_tesseradecagons;
    std::vector<size_t>        l_pentadecagons;

    stress_based_mesh(const mesh_type& msh,
                        const std::vector<size_t>& levels_vec,
                        const bool hanging_nodes,
                        const size_t imsh):m_hanging_nodes(hanging_nodes)
    {

        //check_older_msh(msh);
        cells_marks = std::vector<size_t>(msh.cells_size(),0);
        levels = levels_vec;
        temp_levels.reserve(msh.cells_size());

        //This is only to start only with triangles
        //if we want to have more polygons take out this and solve
        //the problem with the refinement_other
        //#if 0
        if(imsh == 0)
        {
            for(auto& cell: msh)
            {
                auto cell_id = msh.lookup(cell);
                auto c_faces = faces(msh, cell);
                if(c_faces.size() == 4)
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
        auto cl_id   = msh.lookup(cl);
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
                        search_cl = *(msh.cells_begin() + ngh_id);
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
            const T ratio,
            const std::string& directory,
            const std::string& info)
    {
        bool do_refinement = false;

        for(auto& cl : msh)
        {
            auto cl_id    = msh.lookup(cl);
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

        dump_to_matlab(msh, directory +  info + "_wl.m", cells_marks);

        for(auto& cl : msh)
           check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }


    bool
    marker(const mesh_type& msh,
            const std::vector<tensors<T>>& tsr_vec,
            const T yield,
            const std::string& directory,
            const std::string& info)
    {
        typedef disk::quadrature<mesh_type, cell_type>      cell_quad_type;
        typedef dynamic_vector<scalar_type>     vector_type;

        T percent = 0.01;
        std::cout << "INSIDE_XC_MARKER" << std::endl;
        bool do_refinement(false);

        cells_marks = std::vector<size_t>(msh.cells_size());

        for(size_t i = 0; i < tsr_vec.size();i++)
        {
            auto quad_degree = tsr_vec.at(i).quad_degree;
            auto cl  = *(msh.cells_begin() + i);
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
                    if(std::abs(mat1(j) - (1. - percent) * yield) < 1.e-10
                                        || std::abs(mat1(j)) > (1. - percent)*yield )
                    {
                        if(std::abs(mat1(j) - (1. + percent) * yield) < 1.e-10
                                        || std::abs(mat1(j)) < (1. + percent)*yield )
                            in_interval = true;
                    }
                }

                if(std::abs(mat1.maxCoeff() - yield) > 1.e-10 && std::abs(mat1.maxCoeff()) > yield )
                {
                    if(std::abs(mat1(j) - (1. - percent) * yield) < 1.e-10
                                        || std::abs(mat1(j)) > (1. - percent)*yield )
                    {
                        if(std::abs(mat1(j) - (1. + percent) * yield) < 1.e-10
                                        || std::abs(mat1(j)) < (1. + percent)*yield )
                            in_interval = true;
                    }
                }
            }

            for(size_t j = 0; j < mat2.size(); j++)
            {
                if(std::abs(mat2.minCoeff() - yield) <= 1.e-10 || std::abs(mat2.minCoeff()) < yield )
                {
                    if(std::abs(mat2(j) - (1. - percent) * yield) < 1.e-10
                                        || std::abs(mat2(j)) > (1. - percent)*yield )
                    {
                        if(std::abs(mat2(j) - (1. + percent) * yield) < 1.e-10
                                        || std::abs(mat2(j)) < (1. + percent)*yield )
                            in_interval = true;
                    }
                }

                if(std::abs(mat2.maxCoeff() - yield) > 1.e-10 && std::abs(mat2.maxCoeff()) > yield )
                {
                    if(std::abs(mat2(j) - (1. - percent) * yield) < 1.e-10
                                        || std::abs(mat2(j)) > (1. - percent)*yield )
                    {
                        if(std::abs(mat2(j) - (1. + percent) * yield) < 1.e-10
                                        || std::abs(mat2(j)) < (1. + percent)*yield )
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
        dump_to_matlab(msh, directory + "/mesh_pr_" +  info + ".m", cells_marks);

        for(auto& cl : msh)
           check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }
    bool
    marker(const mesh_type& msh,
            const std::vector<tensors<T>>& tsr_vec,
            const T yield,
            const std::string& directory,
            const std::string& info,
            const T percent)
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
                if(std::abs(mat(j) - (1. - percent) * yield) < 1.e-10
                                    || std::abs(mat(j)) > (1. - percent)*yield )
                {
                    if(std::abs(mat(j) - (1. + percent) * yield) < 1.e-10
                                    || std::abs(mat(j)) < (1. + percent)*yield )
                        in_interval = true;
                }
            }
            if(in_interval)
            {
                cells_marks.at(i) = 3;
                ++levels.at(i);
                do_refinement = true;
            }
        }
        dump_to_matlab(msh, directory + "/mesh_pr_" +  info + ".m", cells_marks);
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
                auto cl    = *(msh.cells_begin() + size_t(cl_id));
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
                            const face_type& face, const size_t quad_degree)
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
        if(!find && j >= c_faces.size())
            throw std::invalid_argument("This is a bug: face not found");

        auto fc_pos =  j;
        auto cq  = CellQuadType(quad_degree);
        auto fq  = FaceQuadType(quad_degree);
        //WK: This should change if k different for each face
        auto fqs = fq.integrate(msh, face).size();
        auto cqs = cq.integrate(msh, cell).size();

        return cqs + fqs * fc_pos;
    }



    template<typename CellBasisType, typename CellQuadType,
            typename FaceBasisType, typename FaceQuadType,
             typename GradRecType>
    bool
    marker(const mesh_type& msh,
            const std::vector<tensors<T>>& tsr_vec,
            const T yield,
            const std::string& directory,
            const std::string& info,
            const T percent,
            const std::vector<dynamic_vector<T>>&  Uh_Th,
            const plasticity_data<T>& pst,
            const size_t m_degree)
    {
        typedef CellBasisType                   cell_basis_type;
        typedef FaceBasisType                   face_basis_type;
        typedef FaceQuadType                    face_quadrature_type;
        typedef GradRecType                     grad_reconst_type;
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  tensor_matrix;

        //WK: This is ok  whenever k is equal in all cells and then the number of pts gauss are the same.
        //otherwise this should be inside the for(auto& cl_id : owner_cells)
        auto quad_degree = tsr_vec[0].quad_degree;
        grad_reconst_type   gradrec_nopre(m_degree, quad_degree);

        cell_basis_type                         cell_basis(m_degree + 1);
        face_basis_type                         face_basis(m_degree);
        face_quadrature_type                    face_quadrature(quad_degree);

        std::cout << "INSIDE_AE_MARKER" << std::endl;

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
        std::vector<std::pair<T,size_t>> eta(msh.cells_size());
        //parameters
        auto Cp = 1. / M_PI;
        auto Ct = 0.730276;

        for(size_t cl_id = 0; cl_id < msh.cells_size(); cl_id++)
            eta.at(cl_id).second = cl_id;

        for(auto& cell : msh)
        {
            auto c_faces   =  faces(msh, cell);
            auto cell_id   =  msh.lookup(cell);
            auto c_uh_TF   =  Uh_Th.at(cell_id);
            auto c_tsr     =  tsr_vec.at(cell_id);
            auto meas_T    =  measure(msh, cell);
            auto h_T       =  diameter(msh,cell);

            for( auto& face : c_faces)
            {
                T eta_stress = 0.;
                auto h_F      =  measure(msh, face);
                auto varpi    =  msh.is_boundary(face)? 0.: 0.5 ;

                size_t neighbor_id;
                //WK: ARREGLAR ESTO!!! esto es porque enla forntera no hay elemento vecino
                // Sin embargo, debo esta segura si haciendo el vecino el mismo elemento
                // el estomador en la frontera va a ser cero y mas aun que no vaya afectar
                // los otros estimadores
                if(!msh.is_boundary(face))
                    neighbor_id = face_owner_cells_ids(msh, face, cell);
                else
                    neighbor_id = cell_id;

                auto neighbor  = *(msh.cells_begin() +  neighbor_id);
                auto c_normal  =  normal(msh, cell, face);
                auto n_normal  =  normal(msh, neighbor, face);
                auto nc_uh_TF  =  Uh_Th.at(neighbor_id);
                auto n_tsr     =  tsr_vec.at(neighbor_id);
                auto cgp_id    =  gauss_points_positions<CellQuadType,FaceQuadType>
                                                    (msh,cell,face, quad_degree);
                auto ngp_id    = gauss_points_positions<CellQuadType,FaceQuadType>
                                                (msh,neighbor,face, quad_degree);
                //WK: This should change if k different for each face
                auto fqs = face_quadrature.integrate(msh, face);

                for (auto& qp : fqs)
                {
                    auto c_phi   = cell_basis.eval_functions(msh, cell, qp.point());
                    auto c_dphi  = cell_basis.eval_gradients(msh, cell, qp.point());
                    //Esto debe ser igual que lo de arriba lo que va a cambiar el gradiente es u_TF en cada elemento
                    auto nc_phi  = cell_basis.eval_functions(msh, neighbor, qp.point());
                    auto nc_dphi = cell_basis.eval_gradients(msh, neighbor, qp.point());

                    auto col_range  =   cell_basis.range(1,m_degree+1);
                    auto row_range  =   dof_range(0,mesh_type::dimension);

                    matrix_type c_dphi_matrix  =   make_gradient_matrix(c_dphi);
                    matrix_type c_dphi_taken   =   take(c_dphi_matrix, row_range, col_range);
                    matrix_type c_dphi_rec     =   c_dphi_taken * gradrec_nopre.oper;
                    vector_type c_dphi_ruh     =   c_dphi_rec   * c_uh_TF;

                    matrix_type nc_dphi_matrix =   make_gradient_matrix(nc_dphi);
                    matrix_type nc_dphi_taken  =   take(nc_dphi_matrix, row_range, col_range);
                    matrix_type nc_dphi_rec    =   nc_dphi_taken * gradrec_nopre.oper;
                    vector_type nc_dphi_ruh    =   nc_dphi_rec   * nc_uh_TF;

                    vector_type c_tau  = c_tsr.siglam.col(cgp_id) - pst.alpha * c_tsr.gamma.col(cgp_id);
                    vector_type n_tau  = n_tsr.siglam.col(ngp_id) - pst.alpha * n_tsr.gamma.col(ngp_id);
                    vector_type c_str  = c_tau + pst.alpha *  c_dphi_ruh;
                    vector_type n_str  = n_tau + pst.alpha * nc_dphi_ruh;

                    auto c_value = make_prod_stress_n( c_str , c_normal);
                    auto n_value = make_prod_stress_n( n_str , n_normal);

                    auto jump    = (c_value + n_value);
                    //auto proj_jump = make_prod_stress_n(tau + pst.alpha*dphi_r_uh , n);
                    eta_stress  +=  qp.weight() * varpi * h_F * jump * jump;
                    cgp_id++;
                    ngp_id++;
                }

                auto Cg    =  Cp * Cp  +  varpi * (Ct/meas_T) * h_T * h_T;

                eta.at(cell_id).first +=  Cg * (eta_stress);
            }

        }

        std::sort(eta.begin(), eta.end(),std::greater<std::pair<T, size_t>>());

        T set_error(0.), new_set_error(0.);

        for(auto& e : eta)
            set_error += e.first;

        for(auto& e : eta)
        {
            auto cell_id   = e.second;
            new_set_error += e.first;

            cells_marks.at(cell_id) = 3;
            ++levels.at(cell_id);
            do_refinement  = true;
            if( new_set_error >= percent * set_error )
                break;

        }

        dump_to_matlab(msh, directory + "/mesh_pr_" +  info + ".m", cells_marks);
        for(auto& cl : msh)
            check_levels(msh, cl);

        for(auto& b : cells_marks)
            std::cout<<b<<"  ";

        return do_refinement;
    }


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
            auto cl_id   = msh.lookup(cl);
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
                        auto  neighbor      = *(msh.cells_begin() + size_t(neighbor_id));
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
    template<typename PtA>
    std::pair<bool, std::vector<size_t>>
    is_special_polygon(const mesh_type& msh, const cell_type& cl, const PtA& pts)
    {
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
         std::vector<size_t> vertices(ns.size());
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
                 vertices.at(vcount) =  i;
                 vcount++;
             }
         }
        if(vcount != ns.size())
             std::logic_error(" Incorrect procedure to find vertices");

        bool has_hang_nodes(false);
        if(vertices.size() != pts.size())
            has_hang_nodes = true;

        return std::make_pair(has_hang_nodes, vertices);
    }

    void
    set_marks_hanging_nodes_4tri(mesh_type & msh, const cell_type & cl)
    {
        auto cl_id   = msh.lookup(cl);
        auto cl_mark = cells_marks.at(cl_id);
        auto pts     = points(msh, cl);
        //To adapt also the neighbors take into account 3 and 2, otherwise just 3
        if(cl_mark == 3)
        //if(cl_mark > 1) /
        {

            auto e = is_special_polygon(msh, cl, pts);
            if(e.first && (e.second.size() == 3))
            {
                auto vertices = e.second;
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
                                auto  ngh  = *(msh.cells_begin() + size_t(neighbor_id));
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
                                auto  ngh  = *(msh.cells_begin() + size_t(neighbor_id));
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
        auto cl_id   = msh.lookup(cl);
        auto cl_mark = cells_marks.at(cl_id);
        auto pts     = points(msh, cl);
        //To adapt also the neighbors take into account 3 and 2, otherwise just 3
        if(cl_mark == 3)
        //if(cl_mark > 1) /
        {

            auto e = is_special_polygon(msh, cl, pts);
            if(e.first && (e.second.size() == 3))
            {
                auto vertices = e.second;
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
                                auto  ngh  = *(msh.cells_begin() + size_t(neighbor_id));
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
                                auto  ngh  = *(msh.cells_begin() + size_t(neighbor_id));
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
        auto cl_id   = msh.lookup(cl);
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
                        auto  ngh = *(msh.cells_begin() + size_t(neighbor_id));
                        set_marks(msh,ngh);
                    }

                        #if 0
                        auto  ngh      = *(msh.cells_begin() + size_t(ngh_id));
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
            auto cl_id   = msh.lookup(cl);
            bool cl_mark = (cells_marks.at(cl_id) == 3)? true : false;
            if(cl_mark)
            {
                if(m_hanging_nodes)
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
                  const size_t N, const size_t i)
   {
        auto fo = *(msh.faces_begin() + original_face_id);
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
        auto cl_id   = msh.lookup(cl);
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
    refine_reference_triangle(const n_gon<3> & t , const std::array<int, 3> & bar_ids, const size_t cl_level)
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
        l_triangles.insert(l_triangles.end(), 4, cl_level);
    }
    void refine_triangle(mesh_type& msh, const cell_type& cl, const size_t cl_level)
    {
        auto fcs_ids = cl.faces_ids();
        auto bar0pos = faces_marks.at(fcs_ids[0]).second;
        auto bar1pos = faces_marks.at(fcs_ids[1]).second;
        auto bar2pos = faces_marks.at(fcs_ids[2]).second;
        std::array<int, 3> bar_ids = {bar0pos, bar1pos, bar2pos};

        auto t  = put_n_gon<3>(msh, cl);
        refine_reference_triangle( t,  bar_ids, cl_level);
    }

    void
    refine_other_4tri(mesh_type& msh, const cell_type& cl, const size_t cl_level)
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

            refine_reference_triangle(nt , tbar_ids, cl_level);
        }
    }
    void
    refine_other(mesh_type& msh, const cell_type& cl, const size_t cl_level)
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
            l_triangles.push_back(cl_level);
        }
    }

    template<typename IdxVector>
    void
    refine_special_triangle(mesh_type& msh,
                            const cell_type& cl,
                            const IdxVector& vertices,
                            const size_t cl_level)
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
                auto r_dist  =  ( *(pts.begin() + w) -  *(pts.begin() + v)).to_vector().norm();
                auto b_dist  =  (r_dist - l_dist) / 2.;

                if( std::abs(b_dist) < 1.e-10)
                    throw std::logic_error("Flat triangle");

                auto r_limit  = (i != 2)? vertices.at(i + 1) :  pts.size();
                for(size_t ip = v; ip < r_limit; ip++)
                {
                    auto fc_id  = fcs_ids.at(ip);
                    auto x      = ( *(pts.begin() + ip) - *(pts.begin() + v) ).to_vector().norm();
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
        auto cid  = msh.lookup(cl);
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
        l_triangles.push_back(cl_level);

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
                l_triangles.push_back(cl_level);
            }
            if(new_num_fcs == 4)
            {
                auto nt = put_n_gon<4>(msh, fcs_ids, pts_ids);
                store(nt,m_quadrangles);
                l_quadrangles.push_back(cl_level);
            }
            if(new_num_fcs == 5)
            {
                auto nt = put_n_gon<5>(msh, fcs_ids, pts_ids);
                store(nt,m_pentagons);
                l_pentagons.push_back(cl_level);
            }
            if(new_num_fcs == 6)
            {

                auto nt = put_n_gon<6>(msh, fcs_ids, pts_ids);
                store(nt,m_hexagons);
                l_hexagons.push_back(cl_level);
            }
            if(new_num_fcs == 7)
            {
                auto nt = put_n_gon<7>(msh, fcs_ids, pts_ids);
                store(nt,m_heptagons);
                l_heptagons.push_back(cl_level);
            }
            if(new_num_fcs == 8)
            {
                auto nt = put_n_gon<8>(msh, fcs_ids, pts_ids);
                store(nt,m_octagons);
                l_octagons.push_back(cl_level);
            }
            if(new_num_fcs == 9)
            {
                std::cout << "store in 9angles" << std::endl;
                auto nt = put_n_gon<9>(msh, fcs_ids, pts_ids);
                store(nt,m_enneagons);
                l_enneagons.push_back(cl_level);
            }
            if(new_num_fcs == 10)
            {
                std::cout << "store in 10angles" << std::endl;
                auto nt = put_n_gon<10>(msh, fcs_ids, pts_ids);
                store(nt,m_decagons);
                l_decagons.push_back(cl_level);
            }
            if(new_num_fcs == 11)
            {
                std::cout << "store in 11angles" << std::endl;

                auto nt = put_n_gon<11>(msh, fcs_ids, pts_ids);
                store(nt,m_hendecagons);
                l_hendecagons.push_back(cl_level);
            }
            if(new_num_fcs == 12)
            {
                std::cout << "store in 12angles" << std::endl;

                auto nt = put_n_gon<12>(msh, fcs_ids, pts_ids);
                store(nt,m_dodecagons);
                l_dodecagons.push_back(cl_level);
            }

            if(new_num_fcs == 13)
            {
                std::cout << "store in 13angles" << std::endl;

                auto nt = put_n_gon<13>(msh, fcs_ids, pts_ids);
                store(nt,m_triadecagons);
                l_triadecagons.push_back(cl_level);
            }
            if(new_num_fcs == 14)
            {
                std::cout << "store in 14angles" << std::endl;

                auto nt = put_n_gon<14>(msh, fcs_ids, pts_ids);
                store(nt,m_tesseradecagons);
                l_tesseradecagons.push_back(cl_level);
            }
            if(new_num_fcs == 15)
            {
                std::cout << "store in 15angles" << std::endl;

                auto nt = put_n_gon<15>(msh, fcs_ids, pts_ids);
                store(nt,m_pentadecagons);
                l_pentadecagons.push_back(cl_level);
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
                                     const size_t cl_id)
    {
        auto fcs_ids  = cl.faces_ids();
        auto num_fcs  = fcs_ids.size();
        auto pts      = points(msh, cl);
        auto cid      = msh.lookup(cl);
        auto cl_level = levels.at(cid);
        auto cl_mark = cells_marks.at(cid);
        /* Neighbors Adaptation:
        To adapt marked cells(3) and/or their neighbors(2)
        cells and neihgbors => if(cl_mark > 1)
        Only cells          => if(cl_mark == 3)*/

        if(cl_mark == 3)
        {

            auto  e   = is_special_polygon(msh, cl, pts);

            if(e.first)
            {
                auto vertices = e.second;

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
                        std::cout << "/*********it shoudn't be coming here (until now)!!!! *******/" << std::endl;
                        refine_other(msh , cl, cl_level); //WK: try to do the same for quadrangles
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
                store(nt,m_triangles);
                l_triangles.push_back(cl_level);
            }
            if(new_num_fcs == 4)
            {
                auto nt = put_n_gon_with_hanging_node<4>(msh, cl);
                store(nt,m_quadrangles);
                l_quadrangles.push_back(cl_level);
            }
            if(new_num_fcs == 5)
            {
                auto nt = put_n_gon_with_hanging_node<5>(msh, cl);
                store(nt,m_pentagons);
                l_pentagons.push_back(cl_level);
            }
            if(new_num_fcs == 6)
            {

                auto nt = put_n_gon_with_hanging_node<6>(msh, cl);
                store(nt,m_hexagons);
                l_hexagons.push_back(cl_level);
            }
            if(new_num_fcs == 7)
            {
                auto nt = put_n_gon_with_hanging_node<7>(msh, cl);
                store(nt,m_heptagons);
                l_heptagons.push_back(cl_level);
            }
            if(new_num_fcs == 8)
            {
                auto nt = put_n_gon_with_hanging_node<8>(msh, cl);
                store(nt,m_octagons);
                l_octagons.push_back(cl_level);
            }
            if(new_num_fcs == 9)
            {
                auto nt = put_n_gon_with_hanging_node<9>(msh, cl);
                store(nt,m_enneagons);
                l_enneagons.push_back(cl_level);
            }
            if(new_num_fcs == 10)
            {
                auto nt = put_n_gon_with_hanging_node<10>(msh, cl);
                store(nt,m_decagons);
                l_decagons.push_back(cl_level);
            }
            if(new_num_fcs == 11)
            {
                auto nt = put_n_gon_with_hanging_node<11>(msh, cl);
                store(nt,m_hendecagons);
                l_hendecagons.push_back(cl_level);
            }
            if(new_num_fcs == 12)
            {
                auto nt = put_n_gon_with_hanging_node<12>(msh, cl);
                store(nt,m_dodecagons);
                l_dodecagons.push_back(cl_level);
            }
            if(new_num_fcs == 13)
            {
                auto nt = put_n_gon_with_hanging_node<13>(msh, cl);
                store(nt,m_triadecagons);
                l_triadecagons.push_back(cl_level);
            }
            if(new_num_fcs == 14)
            {
                auto nt = put_n_gon_with_hanging_node<14>(msh, cl);
                store(nt,m_tesseradecagons);
                l_tesseradecagons.push_back(cl_level);
            }
            if(new_num_fcs == 15)
            {
                auto nt = put_n_gon_with_hanging_node<15>(msh, cl);
                store(nt,m_pentadecagons);
                l_pentadecagons.push_back(cl_level);
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
                const std::string & directory,
                const std::string & other_info)
    {
        mesh_type re_msh, re_msh_2;

        re_msh = msh;

        auto storage    = msh.backend_storage();
        auto re_storage = re_msh.backend_storage();
        auto num_cells  = msh.cells_size();
        size_t nds_counter = 0;

        face_marker(re_msh);
        dump_to_matlab(msh, directory + "/mesh" + other_info + ".m",cells_marks);

        std::cout << "_refine_hanging_nodes"<<  m_hanging_nodes << std::endl;


        for(auto& cl : msh)
        {
            auto cl_id   = msh.lookup(cl);
            auto cl_mark = cells_marks.at(cl_id);

            if(cl_mark > 0)
            {
                if(m_hanging_nodes)
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
                        l_triangles.push_back(cl_level);
                        break;
                    case 4:
                        store(put_n_gon<4>(msh,cl), m_quadrangles);
                        l_quadrangles.push_back(cl_level);
                        break;
                    case 5:
                        store(put_n_gon<5>(msh,cl), m_pentagons);
                        l_pentagons.push_back(cl_level);
                        break;
                    case 6:
                        store(put_n_gon<6>(msh,cl), m_hexagons);
                        l_hexagons.push_back(cl_level);
                        break;
                    case 7:
                        store(put_n_gon<7>(msh,cl), m_heptagons);
                        l_heptagons.push_back(cl_level);
                        break;
                    case 8:
                        store(put_n_gon<8>(msh,cl), m_octagons);
                        l_octagons.push_back(cl_level);
                        break;
                    case 9:
                        store(put_n_gon<9>(msh,cl), m_enneagons);
                        l_enneagons.push_back(cl_level);
                        break;
                    case 10:
                        store(put_n_gon<10>(msh,cl), m_decagons);
                        l_decagons.push_back(cl_level);
                        break;
                    case 11:
                        store(put_n_gon<11>(msh,cl), m_hendecagons);
                        l_hendecagons.push_back(cl_level);
                        break;
                    case 12:
                        store(put_n_gon<12>(msh,cl), m_dodecagons);
                        l_dodecagons.push_back(cl_level);
                        break;
                    case 13:
                        store(put_n_gon<13>(msh,cl), m_triadecagons);
                        l_triadecagons.push_back(cl_level);
                        break;
                    case 14:
                        store(put_n_gon<14>(msh,cl), m_tesseradecagons);
                        l_tesseradecagons.push_back(cl_level);
                        break;
                    case 15:
                        store(put_n_gon<15>(msh,cl), m_pentadecagons);
                        l_pentadecagons.push_back(cl_level);
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

        temp_levels = l_triangles;
        temp_levels.insert(temp_levels.end(), l_quadrangles.begin(), l_quadrangles.end());
        temp_levels.insert(temp_levels.end(), l_pentagons.begin()  , l_pentagons.end());
        temp_levels.insert(temp_levels.end(), l_hexagons.begin()   , l_hexagons.end());
        temp_levels.insert(temp_levels.end(), l_heptagons.begin()  , l_heptagons.end());
        temp_levels.insert(temp_levels.end(), l_octagons.begin()   , l_octagons.end());
        temp_levels.insert(temp_levels.end(), l_enneagons.begin()  , l_enneagons.end());
        temp_levels.insert(temp_levels.end(), l_decagons.begin()   , l_decagons.end());
        temp_levels.insert(temp_levels.end(), l_hendecagons.begin(), l_hendecagons.end());
        temp_levels.insert(temp_levels.end(), l_dodecagons.begin() , l_dodecagons.end());
        temp_levels.insert(temp_levels.end(), l_triadecagons.begin()   , l_triadecagons.end());
        temp_levels.insert(temp_levels.end(), l_tesseradecagons.begin(), l_tesseradecagons.end());
        temp_levels.insert(temp_levels.end(), l_pentadecagons.begin()  , l_pentadecagons.end());

        LoaderType      loader;
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

        auto open = save_adapted_mesh(loader, directory + "/amesh" + other_info + ".typ1");
        auto storage_rm2  = re_msh_2.backend_storage();
        loader.populate_mesh(re_msh_2);
        msh = re_msh_2;

        new_levels = std::vector<size_t>(temp_levels.size());

        for(size_t i = 0; i < new_levels.size(); i++)
        {
            size_t idx = loader.m_index_transf.at(i);
            new_levels.at(i) = *(temp_levels.begin() + idx);
        }
        std::cout << "levels_size     :" << levels.size()<< std::endl;
        std::cout << "temp_levels_size :" << temp_levels.size()<< std::endl;
        std::cout << "new_msh_size    :" << msh.cells_size()<< std::endl;

        for(auto& l : new_levels)
            std::cout << l;
        std::cout<< std::endl;
    }
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

template <typename LoaderType>
bool
save_adapted_mesh(const LoaderType& loader, const std::string& filename)
{
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
};




template <typename T>
void
save_data(const std::vector<dynamic_vector<T>>& vec, const std::string& filename)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file"<<std::endl;
    for(auto& v :  vec)
    {
        ofs << v.size() << "  ";
        for (size_t i = 0; i < v.size(); i++)
        {
            ofs << v(i) <<" ";
        }
        ofs << std::endl;
    }
};


template <typename T>
bool
read_data(std::vector<dynamic_vector<T>>& vec , const std::string& filename)
{

    size_t elements_to_read;
    size_t values_to_read;

    std::ifstream ifs(filename);

    if (!ifs.is_open())
        std::cout << "Error opening file"<<std::endl;

    ifs >> elements_to_read;
    std::cout << "Attempting to read " << elements_to_read << " values" << std::endl;

    for (size_t i = 0; i < elements_to_read; i++)
    {
        ifs >> values_to_read;
        dynamic_vector<T>  x;
        for(size_t j = 0; j < values_to_read; j++)
        {
            T val;
            ifs >> val;
            x.at(0,j) = val;
        }
        vec.at(i) = x;
    }
};


template<typename TensorsType, typename T>
std::vector<dynamic_vector<T>>
get_tensor(const std::vector<TensorsType>& tsr_vec, const std::string& name)
{

    std::vector<dynamic_vector<T>> vec(tsr_vec.size());
    size_t i = 0;

    int val;
    if(name == "sigma")
        val = 1;
    else if(name == "gamma")
        val = 2;
    else if(name == "xi_norm")
        val = 3;
    else
        std::cout << "WARNING: Not known name for tensor. Data was not saved." << std::endl;

    switch (val)
    {
        case 1:
            for(auto& tsr : tsr_vec)
            {
                vec.at(i) = tsr.siglam;
                i++;
            }
            break;
        case 2:
            for(auto& tsr : tsr_vec)
            {
                vec.at(i) = tsr.gamma;
                i++;
            }
            break;
        case 3:
            std::cout << "Writing xi_norm" << std::endl;
            for(auto& tsr : tsr_vec)
            {
                vec.at(i) = tsr.xi_norm;
                i++;
            }
            break;
        default :
            std::cout << "WARNING: Data was not saved." << std::endl;
            break;
    }
    return vec;
};

}//end Disk
