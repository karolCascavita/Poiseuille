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

// *****************************************************************************
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
template<typename T>
struct  ADMM_parameters
{
    ADMM_parameters()
    {
        alpha_max = 100.;
        alpha_min = 0.0001;
        beta_max  = 0.9999;
        beta_min  = 0.5;
        tol = 1.e-10;
        delta = 0.5;
        MAX_RESETS = 20;
        set_reset_off  = false;
    }
    T alpha_max;
    T alpha_min;
    T beta_max;
    T beta_min;
    T tol;
    T delta;
    size_t MAX_RESETS;
    bool    set_reset_off;
    friend std::ostream& operator <<(std::ostream& os, const ADMM_parameters<T>& ap)
    {
        os << "ADMM parameters: "<<std::endl;
        os << "* alpha_max : "<< ap.alpha_max<< std::endl;
        os << "* alpha_min : "<< ap.alpha_min<< std::endl;
        os << "* beta_max  : "<< ap.beta_max<< std::endl;
        os << "* beta_min  : "<< ap.beta_min<< std::endl;
        os << "* tolerance : "<< ap.tol<< std::endl;
        os << "* delta     : "<< ap.delta<< std::endl;
        os << "* MAX_RESETS: "<< ap.MAX_RESETS<< std::endl;
        return os;
    }
    void
    set_alpha(const T& h_max, const T& alpha)
    {
        alpha_max =std::max(std::pow(h_max,-2.), alpha);
        return;
    }
    void
    reset_all()
    {
        alpha_max = 100.;
        alpha_min = 0.0001;
        beta_max  = 0.9999;
        beta_min  = 0.5;
        tol = 1.e-10;
        delta = 0.5;
        MAX_RESETS = 20;
        // except reset_off since is readed from the file;
    }
};

//#define DEBUG
template<typename Parameters, typename MeshParameters, typename T,
         typename ADMMParameters>
bool
read_parameters(const std::string& meshfilename,
                const std::string& root,
                const T Bingham,
                Parameters      & pst,
                MeshParameters  & mp,
                ADMMParameters  & ap,
                size_t& max_iters,
                bool  & use_ADMM)
{
    std::string directory = root +  tostr(int(10 * Bingham));

    std::string datafname = directory + "/parameters.txt";
    std::ifstream   ifs(datafname);
    std::string     keyword, temp_name, mesh_name;
    int  method_name,hanging_on,refine_on, num_adapts, temp_recycle, temp_ADMM;
    int  temp_reset_off, temp_start_exact;

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
        std::cout << "directory : "<< directory << std::endl;
        std::regex base_regex2(directory +"/amesh+\\_([0-9]+)\\_+(.*)\\_" + mp.short_mesh_name + "\\.typ1$");

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
            std::cout << "Old_mesh : "<< std::stoi(piece2[1])  << std::endl;
            mp.initial_imsh = std::stoi(piece2[1]) + 1;
            mp.summary_old  = "_" + piece2[2];
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
            ifs >> pst.beta;
        if ( "Binhgam number" == keyword)
            ifs >> pst.Bn;
        if ( "Method ALG"== keyword)
            ifs >> method_name;
        if ( "Start with exact solution"== keyword)
            ifs >> temp_start_exact;
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
        if ( "ADMM"   == keyword)
            ifs >> temp_ADMM;
        if ("Set reset off" == keyword)
            ifs >> temp_reset_off;
        if ( "Other"== keyword)
        {
            std::string other;
            ifs >> other;
            temp_name = other + temp_name;
        }
    }
    mesh_name = temp_name + tostr(mp.marker_name) + mesh_name ;
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
    if (pst.beta <= 0)
    {
        std::cout << "Reduction parameter must belong to (0,1). Falling back to 1." << std::endl;
        pst.beta = 1.;
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
    if(!(temp_ADMM == 1 || temp_ADMM ==0))
        std::cout << "To use ADMM solver choose: Yes (1) or No (0). Falling back to 0." << std::endl;
    if(!(temp_reset_off == 1 || temp_reset_off ==0))
        std::cout << "Reset solutions: Yes (1) or No (0). Falling back to 0." << std::endl;
    if(!(temp_start_exact == 1 || temp_start_exact ==0))
        std::cout << "To start from exact solution choose: Yes (1) or No (0). Falling back to 0." << std::endl;

    use_ADMM = (temp_ADMM == 1);

    ap.set_reset_off = (temp_reset_off == 1);
    mp.num_remesh = num_adapts * refine_on;
    mp.start_exact = (temp_start_exact == 1);


    if(mp.num_remesh < mp.initial_imsh)
        throw std::logic_error("'Number of Adaptations' must be >= than the 'Adaptive step'");

    if(Bingham != pst.Bn)
         throw std::logic_error("Bingham number is not correct in the parameters.txt file");

    mp.directory =  directory;

    mp.summary  ="_N" + tostr(mp.mesh_name) + "_A" +  tostr(int(pst.alpha));

    std::cout << "ALPHA           :  " << pst.alpha    << std::endl;
    std::cout << "Name mesh       :  " << mp.mesh_name << std::endl;
    std::cout << "Adaptation    ? :  " << mp.call_mesher   << std::endl;
    std::cout << "Hanging_nodes ? :  " << mp.hanging_nodes << std::endl;
    std::cout << "Number of adaptations? : " << mp.num_remesh  <<std::endl;
    std::cout << "Estimator Percentage   : " << mp.percent     <<std::endl;
    std::cout << "Initial adaptative step: " << mp.initial_imsh<< std::endl;
    std::cout << " *** summary : "<< mp.summary << std::endl;
    return true;
}

// *****************************************************************************

template<typename MeshType>
struct variable
{
        typedef typename MeshType::scalar_type scalar_type;
        variable(){}
        std::vector<dynamic_vector<scalar_type>>   velocity_Th;
        std::vector<disk::tensors<MeshType>>   stress_Th;

        variable& operator=(const variable& other) {
            velocity_Th     = other.velocity_Th;
            stress_Th = other.stress_Th;
            return *this;
        }
        void
        Zero(const MeshType& msh, const size_t degree)
        {
            velocity_Th     = solution_zero_vector(msh, degree);
            stress_Th = disk::tensor_zero_vector(msh, degree);
        }

};
// *****************************************************************************


// *****************************************************************************

template<typename MeshType>
struct variable2
{
    typedef typename MeshType::scalar_type scalar_type;
    variable2(){}
    std::vector<dynamic_vector<scalar_type>>   velocity_Th;
    disk::tensors2<MeshType>                  stress_Th;

    variable2& operator=(const variable2& other) {
        velocity_Th   = other.velocity_Th;
        stress_Th     = other.stress_Th;
        return *this;
    }
    void
    Zero(const MeshType& msh, const size_t degree)
    {
        velocity_Th = solution_zero_vector(msh, degree);
        stress_Th   = disk::tensor_zero_vector(msh, degree);
    }
};
//template<template<typename, size_t, typename> class Mesh>
template<typename MeshType>
class plasticity2_problem
{};
template< typename T, typename Storage>
class plasticity2_problem <disk::mesh<T,2,Storage>>
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
                            face_quadrature_type>        plast2_type;

public:
    size_t                  m_max_iters, m_degree, m_quad_degree;
    gradrec_type            gradrec_pst;
    stab_type               stabilization;
    statcond_type           statcond_pst;
    assembler_type          assembly_pst;
    plast2_type             plasticity2;
    std::vector<size_t>     levels_vec;
    std::vector<size_t>     ancestors_vec;
    std::vector<size_t>     cells_marks;
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>    solver;

    cell_basis_type         cb;
    face_basis_type         fb;
    cell_quadrature_type    cq;
    cell_quadrature_type    cq_pst;

    plasticity2_problem() = delete;
    plasticity2_problem(variable2<mesh_type> var,
                       const mesh_type& msh,
                       const size_t &     degree,
                       const size_t &     max_iters,
                       const disk::mesh_parameters<T>& mp,
                       const std::string&  root)
    : m_degree(degree), m_max_iters(max_iters)
    {
        std::cout << "hanging_nodes? "<< mp.hanging_nodes << std::endl;
        //check_older_msh(msh);

        //WK: for p-adaptation quad_degree is not the same for all cells
        levels_vec    = std::vector<size_t>(msh.cells_size());
        var.Zero(msh, m_degree);
        auto tsr_0    = var.stress_Th.cells.at(0);
        m_quad_degree = tsr_0.quad_degree();

        std::cout << "m_quad_degree "<< m_quad_degree << std::endl;

        cb = cell_basis_type(m_degree + 1);
        fb = face_basis_type(m_degree);
        cq = cell_quadrature_type(2 * m_degree + 2);
        cq_pst = cell_quadrature_type(m_quad_degree);
    }

    struct bilinear_forms
    {
        bilinear_forms() {}
        matrix_type rec_oper;
        matrix_type cell_lhs;
        vector_type cell_rhs;
    };
    std::vector<matrix_type>
    precompute_gradient(const mesh_type& msh)
    {
        typedef std::pair<matrix_type, matrix_type> gradrec_pair;
        std::vector<matrix_type> vec(msh.cells_size());
        gradrec_type gradrec = gradrec_type(m_degree);
        size_t i = 0;
        for(auto& cl : msh)
        {
            gradrec_pair  gp =  gradrec.compute(msh, cl);
            vec.at(i) =  gp.first;
            i++;
        }
        return vec;
    }

    template<typename LoaderType, typename Solution, typename MeshParameters,
                typename PlasticData, typename Variable>
    auto
    mesh_procedure(const mesh_type      & old_msh,
                    const Variable      & var,
                    const LoaderType    & loader,
                    const Solution      & solution,
                    const PlasticData   & pst,
                    const MeshParameters& mp,
                    const size_t imsh)
    {
        std::cout << "INSIDE MESH_PROCEDURE" << std::endl;

        std::cout << mp << std::endl;
        std::cout << pst << std::endl;

        std::pair< bool, mesh_type> ret;
        mesh_type new_msh;

        disk::dump_to_matlab(old_msh, mp.directory  +  "/old_mesh.m");

        auto info  =  mp.summary + "_RC" + tostr(imsh);

        std::cout << "* Size LEVELS  1: "<< levels_vec.size()  << std::endl;

        disk::stress_based_mesh<mesh_type> sbm(old_msh, levels_vec, pst, mp, imsh);

        if(imsh == 0)
            new_msh = sbm.template refine<LoaderType>(old_msh, m_degree);
        else
        {
            if(!mp.mark_all)
            {
                std::vector<matrix_type> grad_global = precompute_gradient(old_msh);

                auto any = sbm.template marker<cell_basis_type, face_basis_type, Solution>
                        (old_msh, var, pst, m_degree, solution, grad_global);

                if(!any)
                    return std::make_pair(false, old_msh);
            }
            new_msh = sbm.template refine<LoaderType>(old_msh, m_degree);
        }
        auto filename =  mp.directory + "/mesh" + mp.summary + "_R" + tostr(imsh) + ".m";
        dump_to_matlab(new_msh, filename, levels_vec, imsh);

        std::vector<size_t>().swap(levels_vec);
        levels_vec    = sbm.new_levels;
        ancestors_vec = sbm.ancestors;
        cells_marks   = sbm.cells_marks;
        std::cout << "Size LEVELS  2: "<< levels_vec.size()  << std::endl;

        return std::make_pair(true, new_msh);
    }

};
// *****************************************************************************



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

    typedef disk::plasticity<scalar_type,
                            mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type>        plast2_type;

public:
    size_t                  m_max_iters, m_degree, m_quad_degree;
    gradrec_type            gradrec_pst;
    stab_type               stabilization;
    statcond_type           statcond_pst;
    assembler_type          assembly_pst;
    plast_type              plasticity;
    plast2_type             plasticity2;
    std::vector<size_t>     levels_vec;
    std::vector<size_t>     ancestors_vec;
    std::vector<size_t>     cells_marks;
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>    solver;

    cell_basis_type         cb;
    face_basis_type         fb;
    cell_quadrature_type    cq;
    cell_quadrature_type    cq_pst;

    plasticity_problem() = delete;
    plasticity_problem(variable<mesh_type> var,
                       const mesh_type& msh,
                       const size_t &     degree,
                       const size_t &     max_iters,
                       const disk::mesh_parameters<T>& mp,
                       const std::string&  root)
    : m_degree(degree), m_max_iters(max_iters)
    {
        std::cout << "hanging_nodes? "<< mp.hanging_nodes << std::endl;
        //check_older_msh(msh);

        //WK: for p-adaptation quad_degree is not the same for all cells
        levels_vec    = std::vector<size_t>(msh.cells_size());
        var.Zero(msh, m_degree);
        m_quad_degree = var.stress_Th.at(0).sigma.m_quad_degree;

        std::cout << "m_quad_degree "<< m_quad_degree << std::endl;

        cb = cell_basis_type(m_degree + 1);
        fb = face_basis_type(m_degree);
        cq = cell_quadrature_type(2 * m_degree + 2);
        cq_pst = cell_quadrature_type(m_quad_degree);
    }

    struct bilinear_forms
    {
        bilinear_forms() {}
        matrix_type rec_oper;
        matrix_type cell_lhs;
        vector_type cell_rhs;
    };

    std::vector<matrix_type>
    precompute_gradient(const mesh_type& msh)
    {
        typedef std::pair<matrix_type, matrix_type> gradrec_pair;
        std::vector<matrix_type> vec(msh.cells_size());
        gradrec_type gradrec = gradrec_type(m_degree);
        size_t i = 0;
        for(auto& cl : msh)
        {
            gradrec_pair  gp =  gradrec.compute(msh, cl);
            vec.at(i) =  gp.first;
            i++;
        }
        return vec;
    }

    template<typename LoaderType, typename Solution, typename MeshParameters,
                typename PlasticData, typename Variable>
    auto
    mesh_procedure(const mesh_type      & old_msh,
                    const Variable      & var,
                    const LoaderType    & loader,
                    const Solution      & solution,
                    const PlasticData   & pst,
                    const MeshParameters& mp,
                    const size_t imsh)
    {
        std::cout << "INSIDE MESH_PROCEDURE" << std::endl;

        std::cout << mp << std::endl;
        std::cout << pst << std::endl;

        std::pair< bool, mesh_type> ret;
        mesh_type new_msh;

        disk::dump_to_matlab(old_msh, mp.directory  +  "/old_mesh.m");

        auto info  =  mp.summary + "_RC" + tostr(imsh);

        std::cout << "* Size LEVELS  1: "<< levels_vec.size()  << std::endl;

        disk::stress_based_mesh<mesh_type> sbm(old_msh, levels_vec, pst, mp, imsh);

        if(imsh == 0)
            new_msh = sbm.template refine<LoaderType>(old_msh, m_degree);
        else
        {
            if(!mp.mark_all)
            {
                std::vector<matrix_type> grad_global = precompute_gradient(old_msh);

                auto any = sbm.template marker<cell_basis_type, face_basis_type, Solution>
                        (old_msh, var.stress_Th, pst, var.velocity_Th, m_degree,
                                                            solution, grad_global);

                if(!any)
                    return std::make_pair(false, old_msh);
            }
            new_msh = sbm.template refine<LoaderType>(old_msh, m_degree);
        }
        auto filename =  mp.directory + "/mesh" + mp.summary + "_R" + tostr(imsh) + ".m";
        dump_to_matlab(new_msh, filename, levels_vec, imsh);

        std::vector<size_t>().swap(levels_vec);
        levels_vec    = sbm.new_levels;
        ancestors_vec = sbm.ancestors;
        cells_marks   = sbm.cells_marks;
        std::cout << "Size LEVELS  2: "<< levels_vec.size()  << std::endl;

        return std::make_pair(true, new_msh);
    }

    template<typename MeshParameters, typename Variable>
    auto
    solution_procedure( const mesh_type  & new_msh,
                        const mesh_type  & old_msh,
                        const Variable   & old_var,
                        const MeshParameters& mp,
                        const size_t imsh)
    {
        std::cout << "INSIDE SOLUTION_PROCEDURE" << std::endl;


        typedef std::vector<size_t> sizet_vector_type;
        Variable var;

        auto num_cells  = new_msh.cells_size();
        auto info  =  mp.summary + "_RC" + tostr(imsh);

        if(imsh > 0 && mp.recycle)
        {
            //2_Precompute Gradient reconstruction operator
            std::vector<matrix_type> grad_global;
            grad_global =  precompute_gradient(old_msh);

            var.velocity_Th = disk::proj_rec_solution< mesh_type, T,
                            cell_basis_type, cell_quadrature_type,
                                face_basis_type, face_quadrature_type>
            (new_msh, old_msh, grad_global, old_var.velocity_Th, ancestors_vec, m_degree);

            var.stress_Th = sigma_interpolation(new_msh, old_msh,
                                            old_var.stress_Th, ancestors_vec,
                                                        cells_marks, m_degree);



            std::vector<disk::punctual_tensor<mesh_type>> sigma(new_msh.cells_size());
            auto quad_degree = var.stress_Th.at(0).sigma.m_quad_degree;
            for(auto& cl: new_msh)
            {
                size_t id    = cl.get_id();
                sigma.at(id) = var.stress_Th.at(id).sigma;
            }
            std::cout << "SAVING SIGMA_INTER_R" << std::endl;
            auto filename = mp.directory + "/SIGMA_INTER_R" + tostr(imsh)+ ".m";
            disk::get_from_tensor(sigma, var.stress_Th, "sigma");
            disk::quiver_matlab(new_msh, filename, quad_degree, sigma);

            return var;
        }
        else
        {
            //std::cout << "* SOLUTION_PROCEDURE ELSE" << std::endl;
            var.Zero(new_msh, m_degree);
            return var;
        }
    }

    template<typename LoaderType, typename Solution, typename PlasticData,
                typename MeshParameters, typename Variable>
    auto
    test_adaptation( Variable & var,
                    const mesh_type   & old_msh,
                    const LoaderType  & loader,
                    const Solution    & solution,
                    const PlasticData & pst,
                    const MeshParameters& mp,
                    const size_t imsh)
    {
        std::cout << "INSIDE TEST_PROCEDURE" << std::endl;

        bool refine;
        mesh_type new_msh;

        if(mp.call_mesher) //&& imsh > 0)
        {
            if( imsh == mp.initial_imsh)
            {
                if(mp.initial_imsh == 0)
                {
                    auto pair = mesh_procedure(old_msh, var, loader, solution,
                                                                pst, mp, imsh);
                    new_msh = pair.second;
                    var  = solution_procedure(new_msh, old_msh, var, mp, imsh);
                }
                else
                {
                        new_msh = load_mesh_procedure<LoaderType>(old_msh, mp, imsh);
                        if(mp.start_exact)
                            var.Zero(new_msh, m_degree);
                        else
                            var  = load_solution_procedure(new_msh, old_msh, mp, imsh);
                        //var = exact_solution_procedure(new_msh, old_msh, mp, imsh, solution.sf);
                        std::cout << "Size LEVELS: "<< levels_vec.size()  << std::endl;
                }
                return std::make_pair(true, new_msh);
            }
            else
            {
                auto pair = mesh_procedure(old_msh, var, loader, solution, pst, mp, imsh);
                new_msh = pair.second;
                var = solution_procedure(new_msh, old_msh, var, mp, imsh);

                return pair;
            }
        }
        else
        {
            // Make the reading of 3 or 4 meshes.
            // sthg like
            // mesh_bs_2_1.typ1
            // mesh_bs_2_2.typ1
            // mesh_bs_2_3.typ1
            // Till now dont let the algorithm come here until not defining this.
            throw std::logic_error("Define procedure for this case or set call_mesher(true) and No. adaptations (0).");
            return std::make_pair(false, old_msh);
        }
    }
    template<typename LoaderType, typename MeshParameters>
    mesh_type
    load_mesh_procedure(const mesh_type     & old_msh,
                        const MeshParameters& mp,
                        const size_t imsh)
    {
        typedef std::vector<size_t> sizet_vector_type;
        std::cout << "INSIDE LOAD_MESH_PROCEDURE ("<< imsh  <<")"<< std::endl;

        //1_Loading new mesh
        mesh_type      new_msh;
        LoaderType     new_loader;

        auto info_other   =  mp.summary_old + "_" + mp.short_mesh_name + ".typ1";
        auto new_msh_file =  mp.directory + "/amesh_" + tostr(imsh) + info_other;

        if (!new_loader.read_mesh(new_msh_file))
            throw std::logic_error ("Problem loading mesh.");

        new_loader.populate_mesh(new_msh);
        auto num_cells  = new_msh.cells_size();

        //if(mp.recycle)
        //{
            ancestors_vec  = sizet_vector_type(num_cells);

            //2 load ancestors and level info
            std::vector<std::pair<size_t,size_t>> levels_ancestors;
            sizet_vector_type   index_transf;
            disk::load_data(levels_ancestors, index_transf, mp, "/levels","R", imsh);
            assert( num_cells  == levels_ancestors.size());

            //3. Ancestors an levels
            levels_vec = sizet_vector_type(num_cells);

            size_t i = 0;
            for(auto& pair : levels_ancestors)
            {
                levels_vec.at(i) = pair.first;
                ancestors_vec.at(i)  = pair.second;
                i++;
            }
        if(mp.recycle)
        {
            //4. cells_marks
            cells_marks = sizet_vector_type(num_cells);
            disk::load_data(cells_marks, mp, "/marks", "RC", imsh);

        }

        //Updating the mesh
        return new_msh;
    }

    template< typename MeshParameters>
    auto
    load_solution_procedure(const mesh_type     & new_msh,
                            const mesh_type     & old_msh,
                            const MeshParameters& mp,
                            const size_t imsh)
    {
        std::cout << "INSIDE LOAD_SOLUTION_PROCEDURE" << std::endl;

        typedef std::vector<size_t>             sizet_vector_type;
        variable<mesh_type>  var, old_var;

        auto info_other =  mp.summary_old + "_" + mp.short_mesh_name + ".typ1";
        auto num_cells  = new_msh.cells_size();
        std::cout << "* recycle = "<< mp.recycle << std::endl;
        if(!mp.recycle)
        {
            std::cout << "Not recycling former solution" << std::endl;
            var.Zero(new_msh, m_degree);
            return var;
        }
        else
        {
            //Load Uh_Th and tensor data
            disk::load_data(old_var.velocity_Th, mp, "/Uh","R", imsh - 1);
            disk::load_data(old_var.stress_Th, old_msh, mp, "R",imsh);
            assert( old_var.velocity_Th.size()  == old_msh.cells_size());
            assert( old_var.stress_Th.size()== old_msh.cells_size());

            //5_Precompute Gradient reconstruction operator
            std::vector<matrix_type> grad_global;
            grad_global = precompute_gradient(old_msh);

            var.velocity_Th = disk::proj_rec_solution< mesh_type, T,
                        cell_basis_type, cell_quadrature_type,
                            face_basis_type, face_quadrature_type>
            (new_msh, old_msh, grad_global, old_var.velocity_Th, ancestors_vec, m_degree);

            var.stress_Th = sigma_interpolation(new_msh, old_msh,
                                         old_var.stress_Th,
                                         ancestors_vec, cells_marks,  m_degree);


            std::vector<disk::punctual_tensor<mesh_type>> sigma(new_msh.cells_size());
            auto quad_degree = var.stress_Th.at(0).sigma.m_quad_degree;
            for(auto& cl: new_msh)
            {
                 size_t id  = cl.get_id();
                 sigma.at(id)   = var.stress_Th.at(id).sigma;
            }
            std::cout << "SAVING SIGMA_INTER_R" << std::endl;
            auto filename = mp.directory + "/SIGMA_INTER_R" + tostr(imsh)+ ".m";
            disk::get_from_tensor(sigma, var.stress_Th, "sigma");
            disk::quiver_matlab(new_msh, filename, quad_degree, sigma);

            return var;
        }
    }

    template<typename Variable>
    T
    ADMM_convergence(const mesh_type& msh,
                     const Variable & new_var,
                     const Variable & old_var,
                     const T& alpha)
    {
        T sigma_norm(0), Guh_norm(0);

        //5_Precompute Gradient reconstruction operator
        std::vector<matrix_type> grad_global = precompute_gradient(msh);

        for(auto& cl : msh)
        {
            //var.get(cl, sigma);

            auto cqs   = cq.integrate(msh, cl);
            auto cl_id = cl.get_id();
            auto new_all_sigma  = new_var.stress_Th.at(cl_id).sigma;
            auto old_all_sigma  = old_var.stress_Th.at(cl_id).sigma;
            matrix_type    new_sigma   = new_all_sigma.cell_quad_pts;
            matrix_type    old_sigma   = old_all_sigma.cell_quad_pts;;
            matrix_type    sigma_diff  = new_sigma - old_sigma;

            matrix_type    rec_oper  = grad_global.at(cl_id);
            vector_type    uh_diff   = new_var.velocity_Th.at(cl_id)
                                        - old_var.velocity_Th.at(cl_id);
            vector_type    ruh_diff  = rec_oper * uh_diff;

            auto qps_pst = cq_pst.integrate(msh, cl);

            size_t cont = 0;
            for(auto& qp : qps_pst)
            {
                vector_type punctual_diff = sigma_diff.col(cont++);
                sigma_norm += qp.weight() * punctual_diff.dot(punctual_diff);

                vector_type Guh_diff = nabla_ruh_at_point(qp.point(), msh, cl, cb,
                                                            ruh_diff, m_degree);

                Guh_norm +=  qp.weight() * Guh_diff.dot(Guh_diff);
            }
        }
        return std::sqrt(sigma_norm + std::pow(alpha,2.) * Guh_norm);
    }

    template<typename PlasticData, typename ADMMParameters>
    bool
    reset_algorithm_parameters(PlasticData& pst,
                    const ADMMParameters& ap,
                    const T& new_residue,
                    const T& old_residue)
    {

        std::cout << "INSIDE RESET ADMM PARAMETERS" << std::endl;
        auto is_contracted = (new_residue <= pst.beta * old_residue);
        //std::cout << ap << std::endl;
        // 7.a
        if(is_contracted || (pst.alpha == ap.alpha_min && pst.beta == ap.beta_max))
        {
            std::cout << "7.A ACCEPT PARAMETERS" << std::endl;
            auto reset = false;
            return reset;
        }
        // 7.b
        else if( !is_contracted  && pst.alpha > ap.alpha_min)
        {
            std::cout << "7.B CHANGE ALPHA" << std::endl;
            pst.alpha = std::max(ap.delta * pst.alpha,  ap.alpha_min);
            auto reset = false;
            return reset;
        }
        // 7.c
        else if(!is_contracted && pst.alpha == ap.alpha_min && pst.beta < ap.beta_max)
        {
            std::cout << "7.C CHANGE ALPHA and BETA and RESET VARIABLES" << std::endl;
            pst.alpha = ap.alpha_max;
            pst.beta  = std::min(0.5 * (pst.beta + 1.), ap.beta_max);
            if(!ap.set_reset_off)
                return true;
            return false;
        }
        else
            throw std::logic_error("ADMM reachs not defined instance.");

    }

    template<typename MeshParameters>
    void
    output_data(  std::ofstream& ifs,
                  std::ofstream& tfs,
                  std::ofstream& exfs,
                  const MeshParameters& mp,
                  const size_t imsh)
    {
        auto er_steps_name = mp.directory + "/error_by_step_adapt"+ mp.summary +".dat";
        auto er_iters_name = mp.directory + "/error" + mp.summary + "_R" +
                                                                   tostr(imsh) +".dat";

       auto execute_name = mp.directory + "/execute" + mp.summary + ".m";
        ifs = std::ofstream(er_iters_name);
        tfs = std::ofstream(er_steps_name, std::ios::app);
        exfs = std::ofstream(execute_name, std::ios::app);

        //#if 0
        if (!ifs.is_open())
            std::cout << "Error opening file :"<<er_iters_name <<std::endl;
        if (!tfs.is_open())
            std::cout << "Error opening file :"<<er_steps_name <<std::endl;
        if (!exfs.is_open())
            std::cout << "Error opening file :"<<execute_name <<std::endl;
        return;
    }

    struct errors
    {
        errors(): dof(0.), fun(0.), dfun(0.), gamma(0.), sigma(0.)
        {}
        T  dof;
        T  fun;
        T  dfun;
        T  gamma;
        T  sigma;

    };
    int factorial(int n)
    {
        return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
    }
    template<typename PlasticData>
    bool
    finalize(const size_t iter,
             const mesh_type    & msh,
             const PlasticData  & pst,
             const errors       & error,
             const scalar_type  & solver_conv,
             std::ofstream      & tfs)
    {
        std::string name;
        bool is_finished = true;
        if(error.gamma < 1.e-10)
            name = "GAMMA    ";
        else if( iter == m_max_iters - 1 )
            name = "MAX ITERS";
        else if(solver_conv < 1.e-14)
            name = "SOLVER   ";
        else
            is_finished = false;

        auto DIM = mesh_type::dimension;
        auto dimension_N_kd = factorial(m_degree + DIM)/ (factorial(m_degree) * factorial(DIM));
        auto dimension_N_kd_1 = factorial(m_degree + DIM -1)/ (factorial(m_degree) * factorial(DIM-1));
        auto F_DOFs = dimension_N_kd * msh.faces_size();
        auto C_DOFs = dimension_N_kd_1 * msh.cells_size();

        disk::punctual_tensor<mesh_type> pt(msh, *msh.cells_begin(), m_degree);
        auto num_cqps = pt.cell_qps_cols();
        auto num_fqps = pt.face_qps_cols();
        size_t tensor_DOFs = 0;
        for(auto& cl: msh)
            tensor_DOFs += num_fqps * number_of_faces(msh, cl);

        tensor_DOFs += num_cqps * msh.cells_size();


        if(is_finished)
        {
            auto h_max = mesh_h_max(msh);
            auto h_min = mesh_h_min(msh);

            tfs<<m_degree<<"  " << pst.alpha <<"  " <<h_max <<"   " << h_min << error.dof;
            tfs<<m_degree<< "  "<< F_DOFs<< "  "<< C_DOFs <<"   " <<  tensor_DOFs <<std::endl;


            #if 0
            std::cout<< "/***************************************/" << std::endl;
            std::cout <<"/* break by convergence of  " + name + " */"<< std::endl;
            std::cout<< " /***************************************/" << std::endl;
            std::cout << "L2-norm error, dof:         " << error.dof << std::endl;
            std::cout << "l2-norm error, gamma - Du  :" << error.gamma << std::endl;
            std::cout << "l2-norm error, ierations   :" << solver_conv << std::endl;
            std::cout << "Polynomial degree : "<< m_degree <<std::endl;
            std::cout << "Number of cells   : " << msh.cells_size()<<std::endl;
            std::cout << "h_max : "<< h_max  <<std::endl;
            std::cout << "h_min : "<< h_min  <<std::endl;
            #endif
            std::cout <<"/* break by convergence of  " + name + " */"<< std::endl;

        }

        return is_finished;
    }



    template<typename Solution>
    std::vector<bilinear_forms>
    precompute_bilinear_forms(const mesh_type& msh,
                              const Solution & solution)
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
    make_stiff_matrix(mesh_type&  msh,
                      const scalar_type & factor,
                      const bool        & is_kappa_diff,
                      const std::vector<bilinear_forms>& biforms)
    {
        //WARNINGK:if kappa changes in each cell it must be included here.
        if(is_kappa_diff)
        {
            std::cout << "REMAKING STIFF MATRIX" << std::endl;
            assembly_pst.initialize_lhs();
            for (auto& cl : msh)
            {
                auto cl_id   = cl.get_id();
                matrix_type loc_mat = biforms.at(cl_id).cell_lhs;
                matrix_type Ac =  statcond_pst.compute_lhs(msh, cl, factor*loc_mat);
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

    template<typename Solution, typename Variable, typename MeshParameters>
    void
    compute_cell_error(errors& e,
                       const mesh_type& msh, const cell_type& cl,
                       const Variable & var,
                       const Variable & old_var,
                       const vector_type& x, const vector_type& ruh,
                       const Solution& solution,
                        const MeshParameters& mp)
    {
        disk::projector_nopre<mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type> projk(m_degree);

        auto cell_id    = cl.get_id();

        if( solution.is_exact)
        {
            auto qps = cq.integrate(msh, cl);
            auto low_order_range = cb.range(0, m_degree);

            for (auto& qp : qps )
            {
                auto phi  = cb.eval_functions(msh, cl, qp.point());
                scalar_type pot = 0.0;
                for (size_t i = 0; i < low_order_range.size(); i++)
                    pot += phi[i] * x(i);

                vector_type dpot = disk::nabla_ruh_at_point(qp.point(), msh, cl,
                                                            cb, ruh, m_degree);


                size_t number = 0;
                if(mp.diff)
                    number = set_cell_number(msh, cl);
                auto potr  = solution.sf(qp.point(), number);
                vector_type dpotr = solution.df(qp.point(), number);
                scalar_type diff  = 0.0;
                diff = (pot - potr) * (pot - potr) * qp.weight();

                scalar_type d_diff = 0.0;
                d_diff = (dpot - dpotr).dot(dpot - dpotr) * qp.weight();

                e.fun  += diff;
                e.dfun += d_diff;
            }

            vector_type  true_dof = projk.compute_cell(msh, cl, solution.sf);
            vector_type  comp_dof = x.block(0,0,true_dof.size(), 1);
            vector_type  diff_dof = (true_dof - comp_dof);
            e.dof += diff_dof.dot(projk.cell_mm * diff_dof);
        }

        auto qps_pst = cq_pst.integrate(msh, cl);
        size_t col = 0;

        for (auto& qp : qps_pst)
        {
            vector_type dphi_ruh  = disk::nabla_ruh_at_point(qp.point(), msh, cl,
                                                            cb, ruh, m_degree);
            auto all_gamma    = var.stress_Th.at(cell_id).gamma;
            auto all_sigma    = var.stress_Th.at(cell_id).sigma;
            auto all_sigma_old= old_var.stress_Th.at(cell_id).sigma;

            vector_type gamma = all_gamma.cell_quad_pts.col(col);
            vector_type sigma = all_sigma.cell_quad_pts.col(col);
            vector_type sigma_old = all_sigma_old.cell_quad_pts.col(col);

            e.gamma  += qp.weight() * (gamma - dphi_ruh).dot(gamma - dphi_ruh);
            e.sigma  += qp.weight() * (sigma - sigma_old).dot(sigma - sigma_old);

            col++;

        }
        return;
    }
    template<typename BiformsType, typename Solution, typename MeshParameters,
             typename PlasticData, typename Variable>
    void
    make_rhs_vector( mesh_type &    msh,
                    const Variable      &  var,
                    const BiformsType   &  biforms,
                    const PlasticData   &  pst,
                    const Solution      &  solution,
                    const MeshParameters&  mp,
                    const scalar_type   &  kappa,
                    const scalar_type   &  on_pst,
                    const int& imsh)
    {
        assembly_pst.initialize_rhs();

        // This is only if mp.start_exact. I cannot put it inside an if,
        // or I couldnot use it later
        disk::projector_nopre<mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type> projk(m_degree);

        for(auto& cl: msh)
        {
            auto cell_id  = cl.get_id();
            matrix_type rec_oper = biforms.at(cell_id).rec_oper;
            matrix_type loc_mat  = biforms.at(cell_id).cell_lhs * kappa;
            vector_type cell_rhs = biforms.at(cell_id).cell_rhs;

            vector_type pst_rhs;
            vector_type u_local = var.velocity_Th.at(cell_id);
            auto tsr = var.stress_Th.at(cell_id);

            if(mp.diff)  //Only for diffusion test
            {
                std::cout << "Do plasticity.compute for diffusion" << std::endl;
                pst_rhs = plasticity.compute(msh, cl, rec_oper, u_local, pst, solution);
            }
            else
            {
                if(imsh == 0 && solution.is_exact && mp.start_exact)
                    u_local = projk.compute_whole(msh, cl, solution.sf);
                pst_rhs = plasticity.compute(msh, cl, rec_oper, u_local, pst, tsr);
                //pst_rhs = statcond_pst.plasticity_rhs(msh, cl, rec_oper, u_local, pst, tsr);
            }

            //WK: if viscosity changes with domain, it would be mandatory to
            //include mu here, to do it locally
            vector_type Bc = statcond_pst.compute_rhs(msh, cl, loc_mat, cell_rhs,
                                                                on_pst* pst_rhs);
            assembly_pst.assemble_rhs(msh, cl, Bc);
        }
        assembly_pst.impose_boundary_conditions_rhs(msh, solution.sf);

        return;
    }

    template<typename BiformsType, typename Solution, typename MeshParameters,
            typename  PlasticData, typename Variable>
    auto
    recover_solution(Variable    & var,
                      const Variable    &  old_var,
                      const mesh_type   &  msh,
                      const vector_type &  X,
                      const BiformsType &  biforms,
                      const PlasticData &  pst,
                      const Solution    &  solution,
                      const MeshParameters& mp,
                      const scalar_type &  kappa,
                      const scalar_type &  on_pst)
    {
        errors  error;
        size_t  fbs = fb.size();

        for (auto& cl : msh)
        {
            auto cell_id   = cl.get_id();
            auto fcs_ids   = cl.faces_ids();
            auto num_faces = number_of_faces(msh, cl);

            vector_type xFs = vector_type::Zero(num_faces * fbs);

            for (size_t face_i = 0; face_i < num_faces; face_i++)
            {
                auto face_id = fcs_ids.at(face_i);
                vector_type xF = vector_type::Zero(fbs);
                xF = X.block(face_id * fbs, 0, fbs, 1);
                xFs.block(face_i * fbs, 0, fbs, 1) = xF;
            }

            matrix_type rec_oper = biforms.at(cell_id).rec_oper;
            matrix_type loc_mat  = biforms.at(cell_id).cell_lhs * kappa;
            vector_type cell_rhs = biforms.at(cell_id).cell_rhs;
            vector_type pst_rhs;
            vector_type u_local  = var.velocity_Th.at(cell_id);
            if(mp.diff)  //Only for diffusion test
            {
                std::cout << "Do plasticity.compute for diffusion" << std::endl;
                pst_rhs = plasticity.compute(msh, cl, rec_oper, u_local, pst, solution);
            }
            else        //normal plasticity
                pst_rhs  = plasticity.compute(msh, cl, rec_oper, u_local, pst, var.stress_Th.at(cell_id));

            vector_type x = statcond_pst.recover(msh, cl, loc_mat, cell_rhs, xFs,
                                                              on_pst * pst_rhs);
            var.velocity_Th.at(cell_id)  = x;
            vector_type ruh = rec_oper * x;

            size_t number(1);
            if(mp.diff)
                number = disk::set_cell_number(msh, cl);

            compute_cell_error(error, msh, cl, var, old_var, x, ruh, solution, mp);
        }
        error.sigma    = std::sqrt(error.sigma);
        error.gamma    = std::sqrt(error.gamma);
        error.dof      = std::sqrt(error.dof);
        error.fun      = std::sqrt(error.fun);
        error.dfun     = std::sqrt(error.dfun);

        return error;
    }
    template< typename MeshParameters, typename PlasticData, typename Solution>
    auto
    compute_kappa(  scalar_type  &old_kappa,
                    bool  & is_kappa_diff,
                    size_t& on_pst,
                    const PlasticData   & pst,
                    const MeshParameters& mp,
                    const Solution      & solution,
                    const size_t iter,
                    const size_t imsh)
    {
        on_pst    = 1;

        if(mp.diff) //test_diffusion
        {
            std::cout << "KAPPA changing for test diffusion" << std::endl;
            if(solution.identity != "diffusion")
                throw std::logic_error("mp.diff doesn't match with problem type.");
            if(iter == 0)
                is_kappa_diff = true;

            scalar_type  kappa_diff =  pst.mu;
            scalar_type  kappa_pst  =  1.;
            #if 0
            std::cout << "on_ plasticity? : "<< on_pst << std::endl;
            std::cout << "kappa_diffusion : "<< kappa_diff << std::endl;
            std::cout << "kappa_plasticity: "<< kappa_pst << std::endl;
            #endif
            return (on_pst) * kappa_pst + (1 - on_pst) * kappa_diff;
        }
        else
        {
            if(iter == 0)
            {
                if( (mp.recycle && imsh == 0) || (!mp.recycle && !mp.start_exact))
                    on_pst =  0;
                std::cout << "on_pst : "<< on_pst << std::endl;
            }
            scalar_type  kappa_diff =  pst.mu;
            scalar_type  kappa_pst  =  pst.alpha + pst.method * pst.mu;
            scalar_type  kappa = (on_pst) * kappa_pst + (1 - on_pst) * kappa_diff;

            if(kappa != old_kappa)
            {
                is_kappa_diff = true;
                std::cout << "on_pst : "<< on_pst << std::endl;
            }
            #if 0
            std::cout << "on_ plasticity? : "<< on_pst << std::endl;
            std::cout << "kappa_diffusion : "<< kappa_diff << std::endl;
            std::cout << "kappa_plasticity: "<< kappa_pst << std::endl;
            #endif
            return kappa;
        }
    }
    template<typename Solution, typename MeshParameters, typename PlasticData,
            typename Variable>
    void
    test_plasticity( mesh_type& msh,
                    Variable        & var,
                    const PlasticData   & pst,
                    const Solution      & solution,
                    const MeshParameters& mp,
                    std::ofstream       & ifs,
                    std::ofstream       & tfs,
                    std::ofstream       & exfs,
                    const size_t imsh)
    {
        std::cout << "****************************** BEGIN Refinement No.";
        std::cout << imsh<<" *********************************" << std::endl;

        scalar_type  old_error_gamma(0), L2_error(0), kappa(0);
        size_t on_pst;

        assembly_pst  = assembler_type (msh, m_degree);
        gradrec_pst   = gradrec_type(m_degree);
        stabilization = stab_type(m_degree);
        statcond_pst  = statcond_type(m_degree);
        plasticity    = plast_type(m_degree, m_quad_degree);
        auto biforms  = precompute_bilinear_forms(msh, solution);

        for(size_t iter = 0; iter < m_max_iters ; iter++)
        {
            std::cout << "/* ____iter "<< imsh<<"-"<<iter<<"_____ */" << std::endl;
            auto old_var  = var;
            auto is_kappa_diff = false;

            kappa = compute_kappa(kappa, is_kappa_diff, on_pst, pst, mp, solution, iter, imsh);
            make_stiff_matrix(msh, kappa, is_kappa_diff, biforms);
            make_rhs_vector(msh, var, biforms, pst, solution, mp, kappa, on_pst, imsh);
            vector_type X = solver.solve(assembly_pst.rhs);

            auto error = recover_solution(var, old_var, msh, X, biforms, pst,
                                                    solution, mp, kappa, on_pst);

            auto solver_conv = std::abs(old_error_gamma - error.gamma);
            old_error_gamma  = error.gamma;

            ifs<< tostr(iter) <<"  "<< error.gamma << "  " << error.sigma << " ";
            ifs<< error.dof   <<"  "<< solver_conv << std::endl;

            std::cout << "L2-norm error, gamma - Du : " << error.gamma << std::endl;
            std::cout << "L2-norm error, sigma      : " << error.sigma << std::endl;
            std::cout << "L2-norm error, ierations  : " << solver_conv << std::endl;
            std::cout << "L2-norm error, dof        : " << error.dof << std::endl;
            if (finalize(iter, msh, pst, error, solver_conv, tfs))
                break;
        }

        auto er = disk::postprocess(msh, solution,  var.stress_Th, var.velocity_Th,
                                    pst, mp, m_degree, imsh);
        auto grad_global = precompute_gradient(msh);
        if(mp.diff)
        {
            disk::jump_to_matlab<mesh_type, T, cell_basis_type, face_basis_type>
        (msh, solution, var.velocity_Th, grad_global, pst, mp, m_degree, imsh);
        }
        else
        {
            disk::jump_to_matlab<mesh_type, T, cell_basis_type, face_basis_type>
            (msh, var.stress_Th, var.velocity_Th, grad_global, pst, mp, m_degree, imsh);
        }

        disk::execute(exfs, mp, pst, imsh);

        std::cout << "****************************** END Refinement No. ";
        std::cout << imsh<<" *********************************" << std::endl;

        return;
    }

    template< typename PlasticData, typename MeshParameters,
                typename Variable, typename ADMMParameters>

    auto test_ADMM( T& old_residue,
                    PlasticData   & pst,
                    const mesh_type       & msh,
                    const Variable        & new_var,
                    const Variable        & old_var,
                    const ADMMParameters & ap,
                    const MeshParameters  & mp)
    {
        auto ret     = std::make_pair(false, false);
        auto residue = ADMM_convergence(msh, new_var, old_var, pst.alpha);


        std::cout << "INSIDE test_ADMM" << std::endl;
        std::cout << "* ALPHA: "<< pst.alpha << std::endl;
        ret.first = (residue < ap.tol);

        if(ret.first)
            std::cout <<"RESIDUE CONVERGENCE" << std::endl;
        else
        {
            ret.second  = reset_algorithm_parameters(pst, ap, residue, old_residue);

            //std::cout << "/***************** ADMM  ******************/" << std::endl;
            std::cout << "* alpha :" << pst.alpha << "      * beta  :" << pst.beta << std::endl;
        }
        old_residue = residue;

        return ret;
    }

    template<typename Solution, typename MeshParameters, typename PlasticData,
            typename Variable, typename ADMMParameters>
    void
    test_plasticity_ADMM( mesh_type     & msh,
                        Variable        & var,
                        PlasticData     & pst,
                        const Solution  & solution,
                        const MeshParameters& mp,
                        ADMMParameters  & ap,
                        std::ofstream   & ifs,
                        std::ofstream   & tfs,
                        std::ofstream   & exfs,
                        const size_t imsh)
    {
        std::cout << "****************************** BEGIN Refinement No.";
        std::cout << imsh<<" *********************************" << std::endl;

        auto hmax =  disk::mesh_h_max(msh);
        ap.reset_all();
        ap.set_alpha(hmax, pst.alpha);
        pst.alpha = ap.alpha_max;
        pst.beta  = ap.beta_min;

        std::cout << "H_MAX     : "<< hmax << std::endl;
        std::cout << "ALPHA_MAX : "<< ap.alpha_max << std::endl;

        T  residue = 1.;
        int cont = 0;
        auto initial_var  = var;
        //************************** ADMM *******************************//

        scalar_type  old_error_gamma(0), L2_error(0), kappa(0);
        size_t on_pst;

        assembly_pst  = assembler_type (msh, m_degree);
        gradrec_pst   = gradrec_type(m_degree);
        stabilization = stab_type(m_degree);
        statcond_pst  = statcond_type(m_degree);
        plasticity    = plast_type(m_degree, m_quad_degree);
        auto biforms  = precompute_bilinear_forms(msh, solution);

        for(size_t iter = 0; iter < m_max_iters ; iter++)
        {
            //************************** ADMM *******************************//
            auto old_var  = var;
            //************************** ADMM *******************************//

            std::cout << "/* ____iter "<< imsh<<"-"<<iter<<"_____ */" << std::endl;

            auto is_kappa_diff = false;
            kappa = compute_kappa(kappa, is_kappa_diff, on_pst, pst, mp, solution, iter, imsh);
            make_stiff_matrix(msh, kappa, is_kappa_diff, biforms);
            make_rhs_vector(msh, var, biforms, pst, solution, mp, kappa, on_pst, imsh);
            vector_type X = solver.solve(assembly_pst.rhs);
            auto    error = recover_solution(var, old_var, msh, X, biforms, pst,
                                                    solution, mp, kappa, on_pst);

            auto solver_conv = std::abs(old_error_gamma - error.gamma);
            old_error_gamma        = error.gamma;

            std::cout << "L2-norm error, gamma - Du : " << error.gamma << std::endl;
            std::cout << "L2-norm error, ierations  : " << solver_conv << std::endl;
            std::cout << "L2-norm error, dof        : " << error.dof << std::endl;


            //************************** ADMM *******************************//
            //VARIABLE ADMM
            auto pair = test_ADMM(residue, pst, msh, var, old_var, ap, mp);
            ifs<< tostr(iter) <<"  "<< error.gamma << "  "<< error.sigma <<"  ";
            ifs<< error.dof   <<"  "<< solver_conv << "  "<< residue     <<"  ";
            ifs<< pst.alpha   <<"  "<< pst.beta    << std::endl;

            // Residue convergence
            if(pair.first)
                break;
            // Restart variables
            if(pair.second)
            {
                std::cout << "¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨ RESTARTS : "<< imsh <<"-"<<iter;
                std::cout << " "<< cont++<<"¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨" << std::endl;
                var    = initial_var;
            }
            //************************** ADMM *******************************//
            if (finalize(iter, msh, pst, error, solver_conv, tfs))
                break;
        }

        auto er = disk::postprocess(msh, solution, var.stress_Th, var.velocity_Th,
                                    pst, mp, m_degree, imsh);
        auto grad_global = precompute_gradient(msh);
        disk::jump_to_matlab<mesh_type, T, cell_basis_type, face_basis_type>
        (msh, var.stress_Th, var.velocity_Th, grad_global, pst, mp, m_degree, imsh);

        disk::execute(exfs, mp, pst, imsh);

        std::cout << "****************************** END Refinement No. ";
        std::cout << imsh<<" *********************************" << std::endl;

        return;
    }



    template<typename Solution, typename MeshParameters, typename PlasticData,
            typename Variable>
    void
    test_plasticity_exact( mesh_type& msh,
                    Variable        & var,
                    const PlasticData   & pst,
                    const Solution      & solution,
                    const MeshParameters& mp,
                    std::ofstream       & ifs,
                    std::ofstream       & tfs,
                    std::ofstream       & exfs,
                    const size_t imsh)
    {
        std::cout << "****************************** BEGIN Refinement No.";
        std::cout << imsh<<" *********************************" << std::endl;

        scalar_type  old_error_gamma(0), L2_error(0);
        scalar_type  kappa(0), old_kappa(0);

        gradrec_pst   = gradrec_type(m_degree);
        plasticity    = plast_type(m_degree, m_quad_degree);

        disk::projector_nopre<mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type> projk(m_degree);

        auto biforms  = precompute_bilinear_forms(msh, solution);

        for(size_t iter = 0; iter < m_max_iters ; iter++)
        {
            std::cout << "/* ____iter "<< imsh<<"-"<<iter<<"_____ */" << std::endl;

            auto old_var  = var;

            //if(iter == 100 && imsh == 5)
            //    pst.alpha = 0.25*pst.alpha;
            scalar_type  kappa_pst  =  pst.alpha + pst.method * pst.mu;
            scalar_type  kappa_diff =  pst.mu;

            size_t on_pst = 1;
            scalar_type  kappa = (on_pst) * kappa_pst;

            errors  error;
            size_t  fbs = fb.size();
            for (auto& cl : msh)
            {
                 auto cell_id   = cl.get_id();

                 matrix_type rec_oper = biforms.at(cell_id).rec_oper;
                 vector_type u_local  = projk.compute_whole(msh, cl, solution.sf);
                 vector_type pst_rhs  = plasticity.compute(msh, cl, rec_oper,
                                    u_local, pst, var.stress_Th.at(cell_id));

                var.velocity_Th.at(cell_id)  = u_local;
                vector_type ruh = rec_oper * u_local;
                compute_cell_error(error, msh, cl, var, old_var, u_local, ruh, solution, mp);
            }
            error.sigma    = std::sqrt(error.sigma);
            error.gamma    = std::sqrt(error.gamma);
            error.dof      = std::sqrt(error.dof);
            error.fun      = std::sqrt(error.fun);
            error.dfun     = std::sqrt(error.dfun);

            auto solver_conv = std::abs(old_error_gamma - error.gamma);
            old_error_gamma        = error.gamma;

            ifs<<tostr(iter) <<"  "<<error.gamma<<" "<<error.sigma<<" "<<error.dof<<" "<<solver_conv<<std::endl;
            #if 0
            std::cout << "on_ plasticity? : "<< on_pst << std::endl;
            std::cout << "kappa_diffusion : "<< kappa_diff << std::endl;
            std::cout << "kappa_plasticity: "<< kappa_pst << std::endl;
            #endif

            std::cout << "L2-norm error, gamma - Du : " << error.gamma << std::endl;
            std::cout << "L2-norm error, sigma      : " << error.sigma << std::endl;
            std::cout << "L2-norm error, ierations  : " << solver_conv << std::endl;
            std::cout << "L2-norm error, dof        : " << error.dof << std::endl;
            if (finalize(iter, msh, pst, error, solver_conv, tfs))
                break;
        }

        auto er = disk::postprocess(msh, solution,  var.stress_Th, var.velocity_Th,
                                    pst, mp, m_degree, imsh);


        auto grad_global = precompute_gradient(msh);
        disk::jump_to_matlab<mesh_type, T, cell_basis_type, face_basis_type>
        (msh, var.stress_Th, var.velocity_Th, grad_global, pst, mp, m_degree, imsh);

        disk::execute(exfs, mp, pst, imsh);

        std::cout << "****************************** END Refinement No. ";
        std::cout << imsh<<" *********************************" << std::endl;

        return;
    }


};



int main (int argc, char** argv )
{

    using RealType = double;

    int     ch;
    int     degree      = 1;
    size_t  elems_1d    = 10;
    size_t  max_iters   = 100000;
    char    *filename   = nullptr;
    RealType Bingham(0);
    disk::plasticity_data<RealType> pst;
    disk::mesh_parameters<RealType>  mp;
    ADMM_parameters<RealType>        ap;
    bool use_ADMM = false;
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
        //typedef disk::diffusion<RealType,2>                solution_type;
        typedef disk::circular_tuyau<RealType,2>                solution_type;
        typedef plasticity_problem<mesh_type>                   problem_type;
        typedef typename mesh_type::cell                        cell_type;

        //std:: string root = "circular/test/";
        //std:: string root = "circular/borrar/Method2/";
        std:: string root = "circular/trash/";
        //std:: string root = "circular/Meeting7/Diffusion/";

        //std::string root = "square/21_GAUSS_PTS/PER01"
        //std::string root = "square/12_MARKER_4/PER01/5734"
        //std::string root = "square/21_GAUSS_PTS/PER01"
        //std::string root = "square/31_EST_STR/PER01";


        if(!read_parameters(filename, root + "2D_Bi", Bingham, pst, mp, ap, max_iters,
                                                                        use_ADMM))
        {
            std::cout << "Problem loading parameters." << std::endl;
            return 1;
        }
        // WARNINGK: Esto deberia estar definido despues de solution, puesto que
        // ahi se modifica pst.f. El problema es que pst.yield ya se debe conocer
        // para circular.
        pst.yield = 0.5 * pst.Bn * pst.f * pst.Lref;

        solution_type    solution(pst, root);
        mesh_type        msh;
        loader_type      loader;

        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);
        variable<mesh_type> var;
        problem_type     pp(var, msh, degree, max_iters, mp, root);

        if(solution.identity == "diffusion")
        {
             mp.diff = true;
             max_iters = 1;
        }


        for( size_t imsh = mp.initial_imsh; imsh < mp.num_remesh + 1; imsh++)
        {

            std::cout <<  pst << std::endl;
            auto ta_pair = pp.test_adaptation(var, msh, loader, solution, pst, mp, imsh);

            if(!ta_pair.first)
                break;
            msh = ta_pair.second;

            std::ofstream ifs, efs, exfs;
            pp.output_data(ifs, efs, exfs, mp, imsh);

            std::cout <<  pst << std::endl;

            if(use_ADMM)
                pp.test_plasticity_ADMM(msh, var, pst, solution, mp, ap, ifs, efs, exfs, imsh);
            else if(mp.start_exact)
                pp.test_plasticity_exact(msh, var, pst, solution, mp, ifs, efs, exfs, imsh);
            else
                pp.test_plasticity(msh, var, pst, solution, mp, ifs, efs, exfs, imsh);

            std::cout <<  pst << std::endl;

            ifs.close();
            efs.close();
            exfs.close();

        }
        return 0;
    }
    return 0;
};
