#include <iostream>
#include <fstream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <map>
#include "loaders/loader.hpp"
#include "hho/hho.hpp"
#include "hho_pst.hpp"
#include "loaders/mesh_adaptation.hpp"
#include "post_processing_all.hpp"
#include "timecounter.h"
#include <thread>

#include <unsupported/Eigen/SparseExtra>


template<typename InputIt, typename OutputIt, typename UnaryOp>
void parallel_transform(InputIt s_first, InputIt s_last,
                        OutputIt d_first, UnaryOp op)
{
    auto processors     = std::thread::hardware_concurrency();
    processors /= 2;  /* Uncomment if you have Hyperthreading */
    auto elements       = std::distance(s_first, s_last);
    auto slice_size     = elements / processors;
    auto remaining      = elements % processors;

#ifdef DEBUG
    std::cout << "Processors : " << processors << std::endl;
    std::cout << "Elements   : " << elements << std::endl;
    std::cout << "Slice size : " << slice_size << std::endl;
    std::cout << "Remaining  : " << remaining << std::endl;
#endif

    std::vector<std::thread> threads;

    auto thread_lambda = [&](InputIt tib, InputIt tie, OutputIt tob) -> auto {
        while (tib != tie)
            *tob++ = op(*tib++);
    };

    for (size_t i = 0; i < processors; i++)
    {
        auto start_ofs      =   i   * slice_size;
        auto end_ofs        = (i+1) * slice_size;
        auto start_itor     = std::next(s_first, start_ofs);
        auto end_itor       = std::next(s_first, end_ofs);
        auto start_oitor    = std::next(d_first, start_ofs);

#ifdef DEBUG
        std::cout << "Thread " << i << ": ";
        std::cout << start_ofs << " " << end_ofs << std::endl;
        thread_lambda(start_itor, end_itor, start_oitor);
#else
        auto thread = std::thread(thread_lambda, start_itor, end_itor, start_oitor);
        threads.push_back( std::move(thread) );
#endif
    }

    auto rem_ofs            = processors * slice_size;
    auto rem_start_itor     = std::next(s_first, rem_ofs);
    auto rem_end_itor       = s_last;
    auto rem_start_oitor    = std::next(d_first, rem_ofs);

    for (auto& thread : threads)
        thread.join();

#ifdef DEBUG
    std::cout << "Remaining: ";
    std::cout << rem_ofs << " " << elements << std::endl;
#endif
    thread_lambda(rem_start_itor, rem_end_itor, rem_start_oitor);
}



//template<template<typename, size_t, typename> class Mesh>
template<typename MeshType>
class plasticity_problem
{};

template< typename T, typename Storage>
class plasticity_problem <disk::mesh<T,2,Storage>>
{
    typedef disk::mesh<T,2,Storage>                    mesh_type;
    typedef typename mesh_type::point_type             point_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;
    typedef typename mesh_type::scalar_type            scalar_type;

    typedef disk::quadrature<mesh_type, cell_type>     cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>     face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    typedef dynamic_vector<scalar_type>                 vector_type;
    typedef dynamic_matrix<scalar_type>                 matrix_type;
    typedef Eigen::SparseMatrix<scalar_type>            sparse_matrix_type;
    typedef disk::gradient_reconstruction_pst
                            <mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type>        gradrec_type;

    typedef disk::diffusion_like_stabilization_pst
                            <mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type>        stab_type;
    typedef disk::diffusion_like_static_condensation_pst
                            <mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type>        statcond_type;
    typedef disk::assembler_pst<mesh_type,
                            face_basis_type,
                            face_quadrature_type>        assembler_type;
    typedef disk::plasticity<scalar_type,
                            mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type>        plast_type;


public:
    size_t                  m_max_iters, m_degree, m_quad_degree, m_msh_name;
    bool                    m_hg_nd;
    std::string             dir_name, info_all;
    std::ofstream           ifs,ofs;
    gradrec_type            gradrec_pst;
    stab_type               stabilization;
    statcond_type           statcond_pst;
    assembler_type          assembly_pst;
    plast_type              plasticity;
    std::vector<size_t>             levels_vec;
    disk::plasticity_data<T>        m_pst;
    std::vector<dynamic_vector<T>>  Uh_Th;
    std::vector<disk::tensors<T>>   tsr_vec;
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>    solver;

    plasticity_problem() = delete;
    plasticity_problem(const mesh_type& msh,
                       const size_t&     degree,
                       const size_t&     max_iters,
                       const disk::mesh_parameters<T>& mp,
                       const std::string&  root,
                       const disk::plasticity_data<T>&  pst)
    : m_degree(degree), m_max_iters(max_iters), m_pst(pst), m_hg_nd(mp.hanging_nodes)
    , m_msh_name(mp.mesh_name)
    {
        auto method_str  =  (pst.method)? tostr(1): tostr(2);
        info_all =  "_n" + tostr(m_msh_name) + "_g"+ tostr(m_hg_nd) +
                         "_m" + method_str  +  "_a" +  tostr(int(m_pst.alpha));
        std::cout << "hanging_nodes? "<< mp.hanging_nodes << std::endl;

        //WK: for p-adaptation quad_degree is not the same for all cells
        levels_vec    = std::vector<size_t>(msh.cells_size());
        Uh_Th         = disk::solution_zero_vector(msh, m_degree);
        tsr_vec       = disk::tensor_zero_vector(msh, m_degree);
        m_quad_degree = tsr_vec.at(0).quad_degree;
        std::cout << "m_quad_degree "<< m_quad_degree << std::endl;
        std::cout << "size tsr    "<<tsr_vec.at(0).xi_norm.rows()<< " x" << tsr_vec.at(0).xi_norm.cols()<< std::endl;
        dir_name      = disk::directory<T, 2>(m_pst.Bn, root);
    }
    #if 0
    void
    pre_computations(const mesh_type & msh,)
    {

    }
    #endif

    struct bilinear_forms
    {
        bilinear_forms() {}
        matrix_type rec_oper;
        matrix_type cell_lhs;
        vector_type cell_rhs;
    };

    std::vector<matrix_type>
    precompute_gradient(const mesh_type& msh,
                        const gradrec_type& gradrec)
    {
        typedef std::pair<matrix_type, matrix_type> gradrec_pair;
        std::vector<matrix_type> vec(msh.cells_size());
        size_t i = 0;
        for (auto& cl : msh)
        {
            gradrec_pair gp   =  gradrec.compute(msh, cl);
            vec.at(i) =  gp.first;
            i++;
        }
        return vec;
    }

    template<typename LoaderType, typename Solution>
    bool
    adaptation_procedure(mesh_type  & msh,
                    const LoaderType& loader,
                    const Solution& solution,
                    const disk::plasticity_data<T>& pst,
                    const disk::mesh_parameters<T>& mp,
                    const size_t imsh)
    {
        std::cout << "INSIDE ADAPTATION_PROCEDURE" << std::endl;



        mesh_type old_msh = msh;

        if(mp.call_mesher) //&& imsh > 0)
        {
            disk::dump_to_matlab(msh, dir_name  +  "/old_mesh.m");

            auto info  =  info_all + "_rmc" + tostr(imsh);
            disk::stress_based_mesh<mesh_type>
                    sbm(msh, levels_vec, pst, mp, m_hg_nd, imsh, dir_name, info);

            // Leave only if(ims > 0) for normal adaptations. Also reveiw test adaptation if( imsh > 0)
            // (remember this does a  first adaptation  to ensure only triangles).
            if(!mp.mark_all)
            {
                if((imsh > 0))
                {
                    std::vector<matrix_type> grad_global =  precompute_gradient(msh, gradrec_pst);

                    auto do_refinement = sbm.template marker<cell_basis_type, face_basis_type, Solution>
                                            (msh, tsr_vec, m_pst, Uh_Th, m_degree, solution, grad_global);

                    if(!do_refinement)
                        return false;

                    sbm.template refine<LoaderType>(msh, Uh_Th, info_all, m_degree);

                    //if(!is_succesful)
                    //    throw std::logic_error("Refinement not done");

                    if(mp.recycle)
                    {
                        Uh_Th = disk::proj_rec_solution< mesh_type, T,
                                        cell_basis_type, cell_quadrature_type,
                                            face_basis_type, face_quadrature_type,
                                                                        point_type>
                        (msh, old_msh,  grad_global, Uh_Th, sbm.ancestors, m_degree);
                    }
                    else
                        Uh_Th = solution_zero_vector(msh, m_degree);
                }
                else
                {
                    sbm.template refine<LoaderType>(msh, Uh_Th, info_all, m_degree);
                    //Normal adaptation
                    std::vector<dynamic_vector<scalar_type>>().swap(Uh_Th);
                    Uh_Th  = solution_zero_vector(msh, m_degree);
                }
            }
            else    //Marking all
            {
                if(imsh == 0)
                {
                    std::vector<dynamic_vector<scalar_type>>().swap(Uh_Th);
                    Uh_Th  = solution_zero_vector(msh, m_degree);
                }
                else
                {
                    //WARNINGK: Revisar si se requiere refine, y la proyection.
                    //Esto deberia estar aqui, de lo contrario no se va a adaptar nada.
                    //ademas tambien deberia preguntarse si se recicla o no la solution y ahi si proj_rec_solution
                    sbm.template refine<LoaderType>(msh, Uh_Th, info_all, m_degree);

                    std::vector<matrix_type> grad_global =  precompute_gradient(msh, gradrec_pst);

                    Uh_Th = disk::proj_rec_solution< mesh_type, T,
                                    cell_basis_type, cell_quadrature_type,
                                        face_basis_type, face_quadrature_type,
                                                                    point_type>
                    (msh, old_msh,  grad_global, Uh_Th, sbm.ancestors, m_degree);
                }
            }

            std::cout << " imsh : "<< imsh << std::endl;
            dump_to_matlab(msh, dir_name + "/mesh" + info_all + "_rm" + tostr(imsh) + ".m", levels_vec, imsh);

            std::cout << "*************BEGIN Refinement No. "<< imsh<<"****************" << std::endl;

            std::vector<size_t>().swap(levels_vec);
            levels_vec = sbm.new_levels;

            if(mp.recycle && imsh > 0)
                tsr_vec = sigma_interpolation(msh, old_msh, tsr_vec, sbm.ancestors, m_degree);
            else
                tsr_vec = disk::tensor_zero_vector(msh, m_degree);



            // Rev_adapt_mesh_level_ancestor
                std::string pfilename = dir_name + "/Rev_adapt_mesh_level_ancestor.txt";
                std::ofstream pfs(pfilename);
                if(!pfs.is_open())
                    std::cout << "Error opening file"<< pfilename << std::endl;

                auto ancestors = sbm.ancestors;
                size_t i = 0;
                for(auto& cl: msh)
                {
                    auto pts = points(msh, cl);
                    for(auto& p : pts )
                        pfs<< p.x()<< " "<< p.y()<< "  ";

                    auto fcs_ids = cl.faces_ids();
                    for(auto& fid : fcs_ids)
                        pfs<< fid<< " ";
                    pfs << ancestors.at(i) << "  "<< levels_vec.at(i)<<std::endl;
                    i++;
                }
                pfs.close();

        }



        disk::save_data(Uh_Th, dir_name + "/Rev_adap_Uh_th.txt");
        std::vector<matrix_type> sigma;
        disk::get_from_tensor(sigma,  tsr_vec, "sigma");
        disk::save_data(Uh_Th, dir_name + "/Rev_adap_tensor.txt");

        return true;
    };
    #if 0
    {
        bool refine = false;

        if(mp.initial_imsh != 0 && imsh == mp.initial_imsh)
        {
            disk::load_previous_data(pp.Uh_Th, pp.tsr_vec, );
            return true;
        }
        else
            refine = pp.test_adaptation(msh, loader, solution, pst, mp, imsh);

        if( !refine && imsh > mp.initial_imsh)
            return false;
    }
    #endif

    template<typename LoaderType>
    void
    load_data_procedure(mesh_type  & old_msh,
                    const disk::plasticity_data<T>& pst,
                    const disk::mesh_parameters<T>& mp,
                    const size_t imsh)
    {
        std::cout << "INSIDE LOAD_DATA_PROCEDURE" << std::endl;
        //1_Loading new mesh
        mesh_type      new_msh;
        LoaderType     new_loader;

        auto info_other   =  info_all + "_" + mp.short_mesh_name + ".typ1";
        auto new_msh_file =  dir_name + "/amesh_" + "rm_"+ tostr(imsh) + info_other;

        if (!new_loader.read_mesh(new_msh_file))
            throw std::logic_error ("Problem loading mesh.");

        new_loader.populate_mesh(new_msh);

        gradrec_pst   = gradrec_type(m_degree);

        //2_Precompute Gradient reconstruction operator
        std::vector<matrix_type> grad_global = precompute_gradient(old_msh, gradrec_pst);

        //3 load Uh_Th, tsr_vec, ancestors_info
        std::vector<std::pair<size_t,size_t>> levels_ancestors;
        std::vector<size_t> index_transf;
        disk::load_previous_data(Uh_Th, tsr_vec, levels_ancestors, index_transf,
                                                    dir_name, info_all, imsh);

        auto num_cells = new_msh.cells_size();
        levels_vec     = std::vector<size_t>(num_cells);
        auto ancestors = std::vector<size_t>(num_cells);

        assert(num_cells == levels_ancestors.size());
        assert(num_cells == Uh_Th.size());
        assert(num_cells == tsr_vec.size());
        //4. Ancestors an levels
        for(size_t i = 0; i < new_msh.cells_size(); i++)
        {
            //auto idx   = index_transf.at(i);
            auto pair  = *std::next(levels_ancestors.begin(), i);
            levels_vec.at(i) = pair.first;
            ancestors.at(i)  = pair.second;
        }


        std::string pfilename = dir_name + "/Rev_load_mesh_level_ancestor.txt";
        std::ofstream pfs(pfilename);
        if(!pfs.is_open())
            std::cout << "Error opening file"<< pfilename << std::endl;

        size_t i = 0;
        for(auto& cl: new_msh)
        {
            auto pts = points(new_msh, cl);
            for(auto& p : pts )
                pfs<< p.x()<< " "<< p.y()<< "  ";

            auto fcs_ids = cl.faces_ids();
            for(auto& fid : fcs_ids)
                pfs<< fid<< " ";
            pfs << ancestors.at(i) << "  "<< levels_vec.at(i)<<std::endl;
            i++;
        }
        pfs.close();

        Uh_Th = disk::proj_rec_solution< mesh_type, T,
                        cell_basis_type, cell_quadrature_type,
                            face_basis_type, face_quadrature_type,
                                                        point_type>
                    (new_msh, old_msh, grad_global, Uh_Th, ancestors, m_degree);

        tsr_vec = sigma_interpolation(new_msh, old_msh, tsr_vec, ancestors, m_degree);


        disk::save_data(Uh_Th, dir_name + "/Rev_load_Uh_th.txt");
        std::vector<matrix_type> sigma;
        disk::get_from_tensor(sigma,  tsr_vec, "sigma");
        disk::save_data(Uh_Th, dir_name + "/Rev_load_tensor.txt");

    }

    template<typename LoaderType, typename Solution>
    bool
    test_adaptation(disk::mesh<T,2,Storage>  & msh,
                    const LoaderType& loader,
                    const Solution  & solution,
                    const disk::plasticity_data<T>& pst,
                    const disk::mesh_parameters<T>& mp,
                    const size_t imsh)
    {
        bool refine;
        auto info  =  info_all + "_rmc" + tostr(imsh);

        disk::stress_based_mesh<mesh_type>
                sbm(msh, levels_vec, pst, mp, m_hg_nd, imsh, dir_name, info);

        if( imsh == mp.initial_imsh)
        {
            if(mp.initial_imsh == 0)
                refine = adaptation_procedure(msh, loader, solution, pst, mp, imsh);
            else
            {
                if(mp.recycle)
                    load_data_procedure<LoaderType>(msh, pst, mp, imsh);
                else
                    throw std::logic_error("No se ha definido nada para este caso.");
            }
            return true;
        }
        else
        {
            refine = adaptation_procedure(msh, loader, solution, pst, mp, imsh);
            return refine;
        }
    }

    void
    output_data( std::string  & er_iters_name,
                   std::string  & er_steps_name,
                   const size_t imsh)
    {
        er_iters_name = dir_name + "/error" + info_all + "_rm" + tostr(imsh) +".dat";
        er_steps_name = dir_name + "/error_by_step_adapt"+ info_all +".dat";

    }
    void
    finalize(const size_t num,
            const mesh_type& msh,
            const scalar_type error_gamma,
            const scalar_type solver_conv,
            const scalar_type error_dof,
            std::ofstream& tfs)
    {
        std::string name;
        switch (num)
        {
            case 1:
                name = "GAMMA    ";
                break;
            case 2:
                name = "Max Iters";
                break;
            case 3:
                name = "Solver   ";
                break;
            default:
                std::cout << "This shouldn't be coming here" << std::endl;
                break;
        }
        auto h_max = mesh_h_max(msh);
        auto h_min = mesh_h_min(msh);

        tfs<<m_degree<<"  " << m_pst.alpha <<"  " <<h_max <<"   "<< error_dof<<std::endl;

        std::cout<< "/***************************************/" << std::endl;
        std::cout <<"/* break by convergence of  " + name + " */"<< std::endl;
        std::cout<< " /***************************************/" << std::endl;
        std::cout << "L2-norm error, dof:         " << error_dof << std::endl;
        std::cout << "l2-norm error, gamma - Du  :" << error_gamma << std::endl;
        std::cout << "l2-norm error, ierations   :" << solver_conv << std::endl;
        std::cout << "Polynomial degree : "<< m_degree <<std::endl;
        std::cout << "Number of cells   : " << msh.cells_size()<<std::endl;
        std::cout << "h_max : "<< h_max  <<std::endl;
        std::cout << "h_min : "<< h_min  <<std::endl;

        return;
    }



    template<typename Solution>
    std::vector<bilinear_forms>
    precompute_bilinear_forms(const mesh_type& msh,
                              const Solution& solution)
    {
        std::vector<bilinear_forms> vec(msh.cells_size());
        typedef std::pair<matrix_type, matrix_type> gradrec_pair;
        size_t i = 0;
        for (auto& cl : msh)
        {
            bilinear_forms v;
            gradrec_pair gradrec    =  gradrec_pst.compute(msh, cl);
            v.rec_oper      =  gradrec.first;

            matrix_type st_form    =  stabilization.compute(msh, cl, v.rec_oper);
            matrix_type at_form    =  gradrec.second;

            v.cell_lhs      =  at_form + st_form;
            v.cell_rhs      =  disk::compute_rhs<cell_basis_type,
                             cell_quadrature_type>(msh, cl, solution.f, m_degree);

            vec.at(i) = v;
            i++;
        }
        return vec;
    }



    void
    make_stiff_matrix(const mesh_type&  msh,
                      const scalar_type factor,
                      const std::vector<bilinear_forms>& biforms,
                      const size_t iter)
    {
        if(iter < 2)
        {
            assembly_pst.initialize_lhs();

            for (auto& cl : msh)
            {
                auto cl_id   = cl.get_id();
                matrix_type loc_mat = biforms.at(cl_id).cell_lhs;
                matrix_type Ac  =  statcond_pst.compute_lhs(msh, cl, factor * loc_mat);
                assembly_pst.assemble_lhs(msh, cl, Ac);
            }
            assembly_pst.impose_boundary_conditions_lhs(msh);
            assembly_pst.finalize();

            solver.analyzePattern(assembly_pst.matrix);
            solver.factorize(assembly_pst.matrix);
            return;
        }
        return;
    }

    template<typename Solution, typename MeshParameters>
    void
    test_plasticity(const mesh_type& msh,
                    const Solution& solution,
                    const MeshParameters& mp,
                    const size_t imsh)
    {

        scalar_type  error_ALG(0.), L2_error(0.);
        scalar_type  kappa_pst  =  m_pst.alpha + m_pst.method * m_pst.mu;
        scalar_type  kappa_diff =  m_pst.mu;

        cell_basis_type         cb(m_degree + 1);
        face_basis_type         fb(m_degree);
        cell_quadrature_type    cq(2 * m_degree + 2);
        cell_quadrature_type    cq_pst(m_quad_degree);

        std::cout << "RECYCLE : "<< mp.recycle << std::endl;
        if(!mp.recycle)
        {
            std::vector<dynamic_vector<scalar_type>>().swap(Uh_Th);
            Uh_Th    = solution_zero_vector(msh, m_degree);
        }

        assembly_pst  = assembler_type (msh, m_degree);
        gradrec_pst   = gradrec_type(m_degree);
        stabilization = stab_type(m_degree);
        statcond_pst  = statcond_type(m_degree);
        plasticity    = plast_type(m_degree, m_quad_degree,m_pst);

        std::string er_iters_name, er_steps_name;
        output_data(er_iters_name, er_steps_name,imsh);
        std::ofstream ifs(er_iters_name);
        std::ofstream tfs(er_steps_name, std::ios::app);

        std::map<std::string, double> timings;
        timecounter tc_detail, tc;

        //#if 0
        if (!ifs.is_open())
            std::cout << "Error opening ifs"<<std::endl;
        if (!tfs.is_open())
            std::cout << "Error opening tfs"<<std::endl;

        auto biforms = precompute_bilinear_forms(msh, solution);

        // Just to verify the projection of former solution in the first adaptation
        //#if 0
        if(imsh == 1)
        {
            std::string  msh_str_tst = tostr(imsh) + tostr(imsh);
            auto er = disk::plasticity_post_processing(msh, solution.sf,solution.df,
                                tsr_vec, Uh_Th, m_pst,dir_name, m_degree,
                                msh_str_tst,info_all, m_hg_nd, m_msh_name  );
        }
        //#endif
        //Iterations
        for(size_t iter = 0; iter < m_max_iters ; iter++)
        {
            std::cout << "/* ____iter "<< imsh<<"-"<<iter<<"_____ */" << std::endl;

            // This is intended when we are recicling the solution of former adaptation step.
            //otherwise, turn off plasticity in the first Iteration
            size_t on_pst;
            if(mp.recycle)
                on_pst = (!(iter == 0 && imsh == 0))? 1 : 0;
            else
                on_pst = (iter != 0);

            if(mp.diff)
            {
                on_pst    = 1; //This is only for test_diffusion;
                kappa_pst = 1.; // WARNING!!!! solo para test diffusion
            }
            std::cout << "on_ plasticity? : "<< on_pst << std::endl;

            scalar_type kappa = (on_pst) * kappa_pst + (1 - on_pst) * kappa_diff;
            make_stiff_matrix(msh, kappa, biforms, iter);

            //#if 0
            //tc.tic();
            assembly_pst.initialize_rhs();

            for(auto& cl: msh)
            {
                auto cell_id  = cl.get_id();
                matrix_type rec_oper = biforms.at(cell_id).rec_oper;
                matrix_type loc_mat  = biforms.at(cell_id).cell_lhs * kappa;
                vector_type cell_rhs = biforms.at(cell_id).cell_rhs;
                auto tsr = tsr_vec.at(cell_id);
                vector_type pst_rhs;
                if(mp.diff)  //Only for diffusion test
                    pst_rhs  = plasticity.compute(msh, cl, rec_oper, Uh_Th.at(cell_id), solution);
                else         //Uncomment this for plasticity
                    pst_rhs = plasticity.compute(msh, cl, rec_oper, Uh_Th.at(cell_id), tsr);

                //std::cout << "pst_rhs : [ "<< pst_rhs.transpose() << " ]"<< std::endl;
                //WK: if diffusion changes with domain, it would be mandatory to include mu here, to do it locally
                vector_type Bc = statcond_pst.compute_rhs(msh, cl, loc_mat, cell_rhs,
                                                scalar_type(on_pst)* pst_rhs);
                assembly_pst.assemble_rhs(msh, cl, Bc);
                //std::cout << "Bc      : [ "<< Bc.transpose() << " ]"<< std::endl;
            }

            #if 0
            tc.tic();
            std::vector<vector_type> rhs_vec(msh.cells_size());
            auto rhs_op = [&](const cell_type& cl) -> auto
            {
                tc_detail.tic();
                auto gradrec  = gradrec_pst.compute(msh, cl);
                auto rec_oper = gradrec.first;
                auto at_form  = gradrec.second;
                timings["Gradient reconstruction"] += tc_detail.to_double();
                tc_detail.toc();

                tc_detail.tic();
                auto st_form  = stabilization.compute(msh, cl, rec_oper);
                timings["Stabilization"] += tc_detail.to_double();
                tc_detail.toc();

                matrix_type loc   =  kappa * (at_form + st_form);
                auto cell_rhs     =  disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, solution.f, m_degree);
                auto cell_id      =  cl.get_id();
                auto tsr          =  tsr_vec.at(cell_id);

                tc_detail.tic();
                auto pst_rhs = plasticity.compute(msh, cl, rec_oper, Uh_Th.at(cell_id), tsr);
                tc_detail.toc();
                timings["Plasticity"] += tc_detail.to_double();

                //WK: if diffusion changes with domain, it would be mandatory to include mu here, to do it locally
                //Bc
                tc_detail.tic();
                vector_type Bc = statcond_pst.compute_rhs(msh, cl, loc, cell_rhs,
                                                scalar_type(on_pst) * pst_rhs);
                tc_detail.toc();
                timings["Static condensation"] += tc_detail.to_double();

                return  Bc;
            };

            parallel_transform(msh.cells_begin(), msh.cells_end(), rhs_vec.begin(), rhs_op);
            assembly_pst.initialize_rhs();
            for (auto& cl : msh)
            {
                auto cl_id  = cl.get_id();
                auto Bc     = rhs_vec.at(cl_id);
                tc_detail.tic();
                assembly_pst.assemble_rhs(msh, cl, Bc);
                tc_detail.toc();
                timings["RHS_Assembly"] += tc_detail.to_double();
            }
            tc.toc();
            timings["Assembly"] += tc.to_double();
            #endif

            assembly_pst.impose_boundary_conditions_rhs(msh, solution.sf);
            /* SOLVE */
            dynamic_vector<scalar_type> X = solver.solve(assembly_pst.rhs);

            //std::cout << "B : [ "<< assembly_pst.rhs << " ]"<< std::endl;

            //_____________________________________________________________________

            scalar_type error_gamma = 0.0;
            scalar_type error_u_Ruh = 0.0;
            scalar_type error_dof   = 0.0;
            scalar_type error_fun   = 0.0;
            scalar_type error_dfun   = 0.0;


            disk::projector_nopre<mesh_type,
                                cell_basis_type,
                                cell_quadrature_type,
                                face_basis_type,
                                face_quadrature_type> projk(m_degree);

            //std::ofstream ofs( dir_name + "/plotnew.dat");
            auto col_range  =   cb.range(1, m_degree + 1);
            auto row_range  =   disk::dof_range(0,mesh_type::dimension);
            size_t  fbs     =   fb.size();

            for (auto& cl : msh)
            {
                auto cell_id   = cl.get_id();
                auto fcs_ids   = cl.faces_ids();
                auto num_faces = fcs_ids.size();

                dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces * fbs);

                for (size_t face_i = 0; face_i < num_faces; face_i++)
                {
                    auto face_id = fcs_ids.at(face_i);
                    dynamic_vector<scalar_type> xF = dynamic_vector<scalar_type>::Zero(fbs);
                    xF = X.block(face_id * fbs, 0, fbs, 1);
                    xFs.block(face_i * fbs, 0, fbs, 1) = xF;
                }

                matrix_type rec_oper = biforms.at(cell_id).rec_oper;
                matrix_type loc_mat  = biforms.at(cell_id).cell_lhs * kappa;
                vector_type cell_rhs = biforms.at(cell_id).cell_rhs;
                vector_type pst_rhs;
                if(mp.diff)  //Only for diffusion test
                    pst_rhs  = plasticity.compute(msh, cl, rec_oper, Uh_Th.at(cell_id), solution);
                else         //Uncomment this for plasticity
                    pst_rhs  = plasticity.compute(msh, cl, rec_oper, Uh_Th.at(cell_id), tsr_vec.at(cell_id));

                dynamic_vector<scalar_type>   x;
                x = statcond_pst.recover(msh, cl, loc_mat, cell_rhs, xFs, on_pst * pst_rhs);
                Uh_Th.at(cell_id)  = x;
                vector_type ruh = rec_oper * x;

                size_t number(1);
                if(mp.diff)
                    number = disk::set_cell_number(msh, cl);

                //____________________________________________________________________
                if( solution.is_exact)
                {
                    auto qps = cq.integrate(msh, cl);
                    auto low_order_range = cb.range(0, m_degree);
                    for (auto& qp : qps )
                    {
                        auto phi  = cb.eval_functions(msh, cl, qp.point());
                        auto dphi = cb.eval_gradients(msh, cl, qp.point());

                        scalar_type pot = 0.0;
                        for (size_t i = 0; i < low_order_range.size(); i++)
                            pot += phi[i] * x(i);

                        vector_type dpot;
                        matrix_type dphi_matrix =   disk::make_gradient_matrix(dphi);
                        matrix_type dphi_taken  =   disk::take(dphi_matrix, row_range, col_range);
                        matrix_type dphi_ruh =   dphi_taken * ruh;
                        dpot = dphi_ruh;

                        //auto potr  = solution.sf(qp.point());
                        auto potr  = solution.sf(qp.point(), number);
                        vector_type dpotr = solution.df(qp.point(), number);
                        scalar_type diff  = 0.0;
                        diff = (pot - potr) * (pot - potr) * qp.weight();
                        //std::cout << pot << " " << potr << " " << qp.weight() << " " << diff << std::endl;

                        scalar_type d_diff = 0.0;
                        d_diff = (dpot - dpotr).dot(dpot - dpotr) * qp.weight();

                        error_fun  += diff;
                        error_dfun += d_diff;
                        auto tp = qp.point();

                        //for (size_t i = 0; i < mesh_type::dimension; i++)
                        //    ofs << tp[i] << " ";
                        //ofs << pot << " " << std::abs(pot - solution.sf(tp)) << std::endl;
                    }

                    dynamic_vector<scalar_type> true_dof = projk.compute_cell(msh, cl, solution.sf);
                    dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
                    dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
                    error_dof += diff_dof.dot(projk.cell_mm * diff_dof);
                }
                //_____________________________________________________________________

                //____________________________________________________________________

                auto qps_pst = cq_pst.integrate(msh, cl);
                size_t col = 0;

                for (auto& qp : qps_pst)
                {
                    auto   dphi =  cb.eval_gradients(msh, cl, qp.point());
                    matrix_type dphi_matrix =   disk::make_gradient_matrix(dphi);
                    matrix_type dphi_taken  =   disk::take(dphi_matrix, row_range, col_range);
                    vector_type dphi_ruh    =   dphi_taken * ruh;

                    vector_type gamma =  tsr_vec.at(cell_id).gamma.col(col);
                    error_gamma  += qp.weight() * (gamma - dphi_ruh).dot(gamma - dphi_ruh);
                    col++;
                }
            }

            error_gamma    = std::sqrt(error_gamma);
            error_dof      = std::sqrt(error_dof);
            scalar_type    solver_conv = std::abs(error_ALG - error_gamma);
            //error_u_Ruh   = std::sqrt(error_u_Ruh);
            error_ALG      = error_gamma;
            ifs<<tostr(iter) <<"  "<<error_gamma<<" "<<error_dof<<" "<<solver_conv<<std::endl;
            std::cout << "l2-norm error, gamma - Du  :" << error_gamma << std::endl;
            std::cout << "l2-norm error, ierations   :" << solver_conv << std::endl;
            std::cout << "L2-norm error, dof:   " << error_dof << std::endl;


            if(error_gamma < 1.e-10)
            {
                finalize(1, msh, error_gamma, solver_conv, error_dof, tfs);
                break;
            }
            if( iter == m_max_iters - 1 )
            {
                finalize(2, msh, error_gamma, solver_conv, error_dof, tfs);
                break;
            }
            if(solver_conv < 1.e-12)
            {
                finalize(3, msh, error_gamma, solver_conv, error_dof, tfs);
                break;
            }
        }

        std::string  msh_str = tostr(imsh);
        auto er = disk::plasticity_post_processing(msh, solution.sf,solution.df,
                                tsr_vec, Uh_Th, m_pst,dir_name, m_degree,
                                msh_str, info_all, m_hg_nd, m_msh_name);
        ifs.close();
        tfs.close();

        //for (auto& t : timings)
        //    std::cout << " * " << t.first << ": " << t.second << " seconds." << std::endl;

        //std::cout << "L2-norm error, dof:         " << error_dof << std::endl;
        return;
    }


};
//#define DEBUG
template<typename Parameters, typename MeshParameters>
bool
read_parameters(const std::string& meshfilename,
                const std::string& directory,
                Parameters& pst,
                MeshParameters& mp,
                size_t& max_iters)
{
    std::string datafname = directory + "/parameters.txt";
    std::ifstream   ifs(datafname);
    std::string     keyword, temp_name, mesh_name;
    int  method_name,hanging_on,refine_on, num_adapts, temp_recycle;

    std::regex  type_regex;
    std::regex base_regex("(.*+)([0-9]+)\\_+([0-9]+)\\.typ1$");
    //}

    std::smatch  match;
    auto  match_found = std::regex_match(meshfilename, match, base_regex);
    if(match_found)
    {
        std::cout<< " inside 1" <<std::endl;
        std::vector<std::string> piece(match.size());
        for (size_t i = 0; i < match.size(); ++i)
        {
            std::ssub_match submatch = match[i];
            piece[i] = submatch.str();
            std::cout << "  submatch " << i << ": " << piece[i] << '\n';
        }
        mp.short_mesh_name = piece[match.size()-2] +"_"+ piece[match.size()-1];
        mesh_name = piece[match.size()-2] + piece[match.size()-1];
        std::cout << "mesh_name : "<< mesh_name << std::endl;

        std::regex type_regex2;
        std::regex base_regex2(directory +"/amesh_rm\\_+([0-9]+)\\_+(.*)");

        std::smatch  match2;
        auto  match_found2 = std::regex_match(piece[0], match2, base_regex2);
        if(match_found2)
        {
            std::cout<< " " <<std::endl;
            std::vector<std::string> piece2(match2.size());
            for (size_t i = 0; i < match2.size(); ++i)
            {
                std::ssub_match submatch2 = match2[i];
                piece2[i] = submatch2.str();
                std::cout << "  submatch " << i << ": " << piece2[i] << '\n';
            }
            mp.initial_imsh = std::stoi(piece2[1])+1;
            std::cout << "Taking "<< mp.initial_imsh<<" adaptive mesh as the initial one." << std::endl;
            if(mp.initial_imsh ==0)
                throw std::invalid_argument("Don't load data from step 0 (rm_0). Use higher adaptive steps.");
        }

    }

    std::cout << "/* message */" << std::endl;
    if (!ifs.is_open())
    {
        std::cout << "Error opening " << datafname << std::endl;
        return false;
    }

    while(std::getline(ifs,keyword))
    {
        if("Name" == keyword)
        {
            ifs >> temp_name;
            pst.hho = std::stoi(temp_name);
        }
        if("Viscosity" == keyword)
            ifs >>pst.mu;
        if ( "Alpha" == keyword)
            ifs >> pst.alpha;
        if ( "Beta" == keyword)
            ifs >> pst.betha;
        if ( "Binhgam number" == keyword)
            ifs >> pst.Bn;
        if ( "Method ALG"== keyword)
            ifs >> method_name;
        if ( "Adaptation"== keyword)
            ifs >> refine_on;
        if ( "Marker"== keyword)
            ifs >> mp.marker_name;
        if ( "Hanging nodes"== keyword)
            ifs >> hanging_on;
        if ( "Number of Adaptations"== keyword)
            ifs >> num_adapts;
        if ( "Percentage"== keyword)
            ifs >> mp.percent;
        if ( "Recycle"== keyword)
            ifs >> temp_recycle;
        if ( "Other"== keyword)
        {
            std::string other;
            ifs >> other;
            temp_name = other + temp_name;
        }
    }
    mesh_name = temp_name + mesh_name + tostr(mp.marker_name);
    mp.mesh_name = std::stoi(mesh_name);

    if (mp.mesh_name < 0)
    {
        std::cout << "Mesh name not specified. Falling back to 0." << std::endl;
        mp.mesh_name = 0;
    }
    if (pst.mu <= 0)
    {
        std::cout << "Viscosity must be bigger than 0. Falling back to 1." << std::endl;
        pst.mu = 1.0;
    }
    if (pst.alpha <= 0)
    {
        std::cout << "Augmentation parameter must be bigger than 0. Falling back to 1." << std::endl;
        pst.alpha = 1.;
    }
    if (pst.betha <= 0)
    {
        std::cout << "Reduction parameter must belong to (0,1). Falling back to 1." << std::endl;
        pst.betha = 1.;
    }
    if (pst.Bn < 0)
    {
        std::cout << "Bingham number must be positive. Falling back to 0. Solving diffusion case." << std::endl;
        pst.Bn    = 0.0;
        max_iters = 1;
    }
    if(pst.Bn == 0.)
        max_iters = 1;
    if( method_name == 1)
        pst.method = true;
    else
    {
        pst.method = false;
        if ( method_name != 2)
            std::cout << "Method must be: Glowinski (1) or Saramito (2). Coming back to Saramito." << std::endl;
    }
    if( refine_on == 1)
    {
        mp.call_mesher   = true;
        if( hanging_on == 1)
            mp.hanging_nodes = true;
    }
    else
    {
        mp.hanging_nodes = false;
        mp.call_mesher   = false;
        if ( refine_on != 0)
            std::cout << "Refinement: Yes (1) or No (0). Coming back to 0." << std::endl;
    }
    if (num_adapts < 0)
    {
        std::cout << "Max number of adaptations must be positive. Falling back to 0." << std::endl;
        num_adapts = 0;
    }

    if(temp_recycle == 1)
        mp.recycle   = true;
    else
    {
        mp.recycle = false;
        if ( refine_on != 0)
            std::cout << "Refinement: Yes (1) or No (0). Falling back to 0." << std::endl;
    }
    if (mp.initial_imsh < 0)
    {
        std::cout << "Initial adaptive step must be positive. Falling back to 0." << std::endl;
        mp.initial_imsh = 0;
    }

    mp.num_remesh = num_adapts * refine_on;


    std::cout << "ALPHA           :  " << pst.alpha    << std::endl;
    std::cout << "Name mesh       :  " << mp.mesh_name << std::endl;
    std::cout << "Adaptation    ? :  " << mp.call_mesher   << std::endl;
    std::cout << "Hanging_nodes ? :  " << mp.hanging_nodes << std::endl;
    std::cout << "Number of adaptations? : " << mp.num_remesh  <<std::endl;
    std::cout << "Estimator Percentage   : " << mp.percent     <<std::endl;
    std::cout << "Initial adaptative step: " << mp.initial_imsh<< std::endl;

    if(mp.num_remesh < mp.initial_imsh)
        throw std::logic_error("'Number of Adaptations' must be >= than the 'Adaptive step'");

    return true;
}
#if 0

residual()
{

}
convergence()
{
    Rnew = residual();
    Rold;
    if( !(Rnew < Rold))
        pst.alpha = 0.5 * pst.alpha;

    if( !pst.alpha < alpha.ref )
        pst.betha = 0.5 * (betha + 1.);
}
#endif
int main (int argc, char** argv )
{

    using RealType = double;

    int     ch;
    int     degree      = 1;
    size_t  elems_1d    = 10;
    size_t  max_iters   = 200000;
    char    *filename   = nullptr;
    RealType Bingham(0.);
    disk::plasticity_data<RealType> pst;
    disk::mesh_parameters<RealType>  mp;

    pst.method = false;

    while ((ch = getopt(argc, argv, "i:k:B:")) != -1 )
    {
        switch(ch)
        {
            case 'i':
            max_iters = atoi(optarg);
            if (max_iters < 0)
            {
                std::cout << "Max number of iterations must be positive. Falling back to 1." << std::endl;
                max_iters = 1;
            }
            break;
            case 'k':
            degree = atoi(optarg);
            if (degree < 0)
            {
                std::cout << "Degree must be positive. Falling back to 1." << std::endl;
                degree = 1;
            }
            case 'B':
                Bingham = atof(optarg);
                if (Bingham < 0)
                {
                    std::cout << "Bingham number must be positive. Falling back to 0. Solving diffusion case." << std::endl;
                    Bingham = 0;
                    max_iters    = 1;
                }
                break;
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);

        }
    }


    argc -= optind;
    argv += optind;

    filename = argv[0];

    if (filename == nullptr)
    {
        std::cout << "no filename specified" << std::endl;
        return 1;
    }

    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<RealType, 2>                 mesh_type;
        typedef disk::fvca5_mesh_loader<RealType, 2>            loader_type;
        typedef typename mesh_type::scalar_type                 scalar_type;
        typedef Eigen::Matrix< RealType, 2, Eigen::Dynamic>     tensor_matrix;
        typedef disk::tuyau<RealType,2>                     solution_type;
        typedef plasticity_problem<mesh_type>                   problem_type;
        typedef typename mesh_type::cell                        cell_type;

        std:: string root = "square/test";
        //std::string root = "square/21_GAUSS_PTS/PER01"
        //std::string root = "square/12_MARKER_4/PER01/5734"
        //std::string root = "square/21_GAUSS_PTS/PER01"
        //std::string root = "square/31_EST_STR/PER01";

        std::string directory = root + "/2D_Bi" + tostr(int(10 * Bingham));

        if(!read_parameters(filename, directory, pst, mp, max_iters))
        {
            std::cout << "Problem loading parameters." << std::endl;
            return 1;
        }
        if(Bingham != pst.Bn)
             throw std::logic_error("Bingham number is not correct in the parameters.txt file");

        pst.yield = 0.5 * pst.Bn * pst.f * pst.Lref;

        solution_type    solution(pst,root);
        mesh_type        msh;
        loader_type      loader;

        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);
        problem_type     pp(msh, degree, max_iters, mp, root, pst);

        if(solution.identity == "diffusion")
            mp.diff = true;


        for( size_t imsh = mp.initial_imsh; imsh < mp.num_remesh + 1; imsh++)
        {

            bool go_on = pp.test_adaptation(msh, loader, solution, pst, mp, imsh);

            if(!go_on)
                break;

            pp.test_plasticity(msh, solution, mp, imsh);

            std::cout << "************* END Refinement No. "<< imsh <<"****************"<<  std::endl;
        }



        return 0;
    }
    return 0;
};
