#include <iostream>
#include <fstream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <map>
#include "loaders/loader.hpp"
#include "hho/hho.hpp"
#include "hho_pst.hpp"
#include "loaders/mesh_adaptation_vpst.hpp"
#include "post_processing_all.hpp"
#include "timecounter.h"
#include <thread>
#include <unsupported/Eigen/SparseExtra>
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
    int  method_name,hanging_on,refine_on, binning_on, num_adapts, temp_recycle, temp_ADMM;
    int  temp_reset_off, temp_start_exact, binNadapt_on;

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
            mp.hho = std::stoi(temp_name);
            if (std::stoi(temp_name) > 1)
            {
                std::cout << "Equal face and cell degrees (0) or higher degree for cells (1). Coming back to the first option." << std::endl;
                mp.hho = 0;
            }

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
        if ( "Binning"== keyword)
            ifs >> binning_on;
        if ( "Binning + adaptation"== keyword)
                ifs >> binNadapt_on;
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
        std::cout << "Refine on!!!!!!!" << std::endl;
        mp.call_mesher   = true;
        if( hanging_on == 1)
            mp.hanging_nodes = true;
        if( binning_on == 1)
        {
            mp.binning  = true;
            if( binNadapt_on == 1)
            {
                std::cout << "Binning and Binning+Adaptation on: Binning chosen." << std::endl;
                mp.binNadapt = false;
            }
        }
        else
        {
            if(binNadapt_on == 1)
                mp.binNadapt  = true;
        }
    }
    else
    {
        mp.hanging_nodes = false;
        mp.call_mesher   = false;
        mp.binning = false;
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
    std::cout << "Binning       ? :  " << mp.binning<< std::endl;
    std::cout << "Hanging_nodes ? :  " << mp.hanging_nodes << std::endl;
    std::cout << "Number of adaptations? : " << mp.num_remesh  <<std::endl;
    std::cout << "Estimator Percentage   : " << mp.percent     <<std::endl;
    std::cout << "Initial adaptative step: " << mp.initial_imsh<< std::endl;
    std::cout << " *** summary : "<< mp.summary << std::endl;
    return true;
}

template<typename MeshType>
struct variable2
{
    typedef typename MeshType::scalar_type scalar_type;
    variable2(){}
    std::vector<dynamic_vector<scalar_type>>  velocity;
    disk::tensors<MeshType>                   stresses;

    variable2& operator=(const variable2& other) {
        velocity   = other.velocity;
        stresses   = other.stresses;
        return *this;
    }
    template<typename CellBasisType, typename FaceBasisType>
    void
    Zero(const MeshType& msh, const size_t degree, const bool& change_degree)
    {
        velocity = solution_zero_vector(msh, degree, change_degree);
        stresses.template zero_vector<CellBasisType,FaceBasisType>(msh, degree);
    }
};

template<typename MeshType,  typename Variable>
class plasticity2_problem
{};
template< typename T, typename Storage>
class plasticity2_problem <disk::mesh<T,2,Storage>, variable2<disk::mesh<T,2,Storage>>>
{
    typedef disk::mesh<T,2,Storage>                    mesh_type;
    typedef disk::tensors<mesh_type>                   tensors_type;
    typedef variable2<disk::mesh<T,2,Storage>>         variable_type;
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
    typedef disk::gradient_reconstruction_nopre
                            <mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type>        gradrec_type;
    typedef disk::diffusion_like_stabilization_nopre
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
    typedef disk::plasticity2<scalar_type,
                            mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type,
                            tensors_type>                plasticity_type;

    cell_basis_type         cb;
    face_basis_type         fb;
    cell_quadrature_type    cq;
    face_quadrature_type    fq;

    gradrec_type            gradrec_pst;
    stab_type               stabilization;
    statcond_type           statcond_pst;
    assembler_type          assembly_pst;
    plasticity_type         plasticity;

    std::vector<size_t>     levels_vec;
    std::vector<size_t>     ancestors_vec;
    std::vector<size_t>     cells_marks;

    size_t                  m_max_iters, m_degree_cell, m_degree;
    bool                    change_degree;

    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>    solver;

public:
    plasticity2_problem() = delete;
    plasticity2_problem(variable_type   &   var,
                       const mesh_type  &   msh,
                       const size_t     &   degree,
                       const size_t     &   max_iters,
                       const disk::mesh_parameters<T>& mp,
                       const std::string&  root)
                       : m_degree(degree), m_degree_cell(degree),
                        m_max_iters(max_iters)
    {
        std::cout << "hanging_nodes? "<< mp.hanging_nodes << std::endl;
        //check_older_msh(msh);

        //WK: for p-adaptation quad_degree is not the same for all cells
        levels_vec    = std::vector<size_t>(msh.cells_size());
        var.template Zero<cell_basis_type, face_basis_type>(msh, m_degree, mp.hho);

        m_degree_cell = degree + mp.hho;

        cb = cell_basis_type(m_degree_cell);
        fb = face_basis_type(m_degree);
        cq = cell_quadrature_type(2 * m_degree_cell + 2);
        fq = face_quadrature_type(2 * m_degree + 2);
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

        if(imsh == mp.initial_imsh)
            exfs = std::ofstream(execute_name);
        else
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

    struct bilinear_forms
    {
        bilinear_forms() {}
        matrix_type cell_lhs;
        vector_type cell_rhs;
    };

    template<typename Solution>
    std::vector<bilinear_forms>
    precompute_bilinear_forms(const mesh_type& msh,
                const Solution & solution,
                const disk::mesh_parameters<T>& mp,
                std::vector<matrix_type>& grec_opers,
                std::vector<matrix_type>& stab_opers)
    {
        stab_opers  = std::vector<matrix_type>(msh.cells_size());
        grec_opers = std::vector<matrix_type>(msh.cells_size());

        std::vector<bilinear_forms> vec(msh.cells_size());
        typedef std::pair<matrix_type, matrix_type> mat_mat_pair;
        size_t i = 0;

        for (auto& cl : msh)
        {

            // reconsruction in the space \nabla P^{k+1}
            gradrec_pst.compute(msh, cl);
            // grad reconsruction in the space [P^{k}]^dim
            gradrec_pst.compute_low_order(msh, cl);

            if(mp.hho)
                stabilization.compute_low_order(msh, cl, gradrec_pst.gradrecoper);
            else
                stabilization.compute(msh, cl, gradrec_pst.oper);

            //save opers
            auto cl_id = msh.lookup(cl);
            grec_opers.at(cl_id)   = gradrec_pst.gradrecoper;
            stab_opers.at(cl_id)   = stabilization.oper;

            matrix_type at_form    =  gradrec_pst.data;
            matrix_type st_form;

            if(mp.hho)
                st_form    =  stabilization.low_order_data;
            else
                st_form    =  stabilization.data;//low_order_data;


            bilinear_forms v;
            v.cell_lhs  =  at_form + st_form;
            v.cell_rhs  =  disk::compute_rhs<cell_basis_type, cell_quadrature_type>
                                            (msh, cl, solution.f, m_degree_cell);

            //std::cout << "m_degree_cell : "<< m_degree_cell << std::endl;
            //std::cout << "v.cell_rhs    : "<< v.cell_rhs.size() << std::endl;
            vec.at(i) = v;
            i++;
        }
        return vec;
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
                //std::cout << "on_pst : "<< on_pst << std::endl;
            }
            scalar_type  kappa_diff =  pst.mu;
            scalar_type  kappa_pst  =  pst.alpha + pst.method * pst.mu;
            scalar_type  kappa = (on_pst) * kappa_pst + (1 - on_pst) * kappa_diff;

            if(kappa != old_kappa)
            {
                is_kappa_diff = true;
                //std::cout << "on_pst : "<< on_pst << std::endl;
            }
            #if 0
            std::cout << "on_ plasticity? : "<< on_pst << std::endl;
            std::cout << "kappa_diffusion : "<< kappa_diff << std::endl;
            std::cout << "kappa_plasticity: "<< kappa_pst << std::endl;
            #endif
            return kappa;
        }
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
                matrix_type Ac =  statcond_pst.compute_lhs(msh, cl, factor * loc_mat);
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

    template<typename BiformsType, typename Solution, typename MeshParameters,
             typename PlasticData>
    void
    make_rhs_vector( mesh_type &    msh,
                    const variable_type &  var,
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
                            face_quadrature_type> projk(m_degree, mp.hho);

        for(auto& cl: msh)
        {
            auto cell_id  = cl.get_id();

            matrix_type loc_mat  = biforms.at(cell_id).cell_lhs * kappa;
            vector_type cell_rhs = biforms.at(cell_id).cell_rhs;

            vector_type pst_rhs;

            vector_type u_local = var.velocity.at(cell_id);
            //std::cout << "u_TF size"<< u_local.size() << std::endl;

            if(mp.diff)  //Only for diffusion test
            {
                std::cout << "Do plasticity.compute for diffusion" << std::endl;
                //REview this for the new plasticity
                //pst_rhs = plasticity.compute(msh, cl, rec_oper, u_local, pst, solution);
            }
            else
            {
                if(imsh == 0 && solution.is_exact && mp.start_exact)
                    u_local = projk.compute_whole(msh, cl, solution.sf);

                pst_rhs = plasticity.compute_integral(msh, cl, var.stresses);
                //pst_rhs = vector_type::Zero(u_local.size());
            }

            //WK: if viscosity changes with domain, it would be mandatory to
            //include mu here, to do it locally
            vector_type Bc = statcond_pst.compute_rhs(msh, cl, loc_mat, cell_rhs,
                                                            on_pst * pst_rhs);
            assembly_pst.assemble_rhs(msh, cl, Bc);
        }
        assembly_pst.impose_boundary_conditions_rhs(msh, solution.sf);

        return;
    }
    struct errors
    {
        errors(): dof(0.), dof_grad(0.), fun(0.), dfun(0.), dof_sigma(0.),
                  gamma(0.), psi(0.),  sigma(0.), varsigma(0.),
                  stab_u(0.), grad_u(0.), residue_all(0.), residue(0.),
                  solver_conv(0.)

        {}
        T  dof  ;  T  dof_grad;     T  fun; T  dfun; T   dof_sigma;
        T  gamma;  T  psi;
        T  sigma;  T  varsigma;
        T  stab_u; T  grad_u;
        T  residue_all;    T  residue;
        T solver_conv;

        void
        set_zero_discrete_norms()
        {
            gamma  = 0.; psi    = 0 ; sigma = 0.      ; varsigma = 0.;
            stab_u = 0.; grad_u = 0.; residue_all = 0.; residue = 0.;
            solver_conv = 0.;
        }
        void
        sqrt_discrete_norms()
        {
            gamma  = std::sqrt(gamma);
            psi    = std::sqrt(psi);

            sigma    = std::sqrt(sigma);
            varsigma = std::sqrt(varsigma);

            stab_u  = std::sqrt(stab_u);
            grad_u  = std::sqrt(grad_u);

            residue_all = std::sqrt(residue_all);
            residue = std::sqrt(residue);
            return;
        }

        void
        set_zero_exact_norms()
        {
            dof  = 0.; dof_grad = 0.; fun = 0.; dfun = 0.; dof_sigma = 0.;
            return;
        }

    };

    int factorial(int n) { return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}

    template<typename Solution, typename MeshParameters>
    void
    compute_exact_error(errors& e,
                       const mesh_type      & msh,
                       const cell_type      & cl,
                       const variable_type  & var,
                       const variable_type  & old_var,
                       const vector_type    & x,
                       const vector_type    & ruh,
                       const Solution       & solution,
                       const MeshParameters & mp)
    {
        disk::projector_nopre<mesh_type,
                            cell_basis_type,
                            cell_quadrature_type,
                            face_basis_type,
                            face_quadrature_type> projk(m_degree, mp.hho);

        auto cell_id    = cl.get_id();
        auto cell_quadpoints = cq.integrate(msh, cl);

        auto mindeg = 0;
        auto maxdeg = m_degree;

        auto row_range  = disk::dof_range( 0, mesh_type::dimension);
        auto zero_range = cb.range( mindeg,  maxdeg);

        auto DIM = mesh_type::dimension;
        auto cbs = zero_range.size();
        auto vec_cbs = DIM * cbs;

        matrix_type mass_mat  = matrix_type::Zero(vec_cbs, vec_cbs);

        /* LHS: take basis functions derivatives from degree 1 to K+1 */

        for (auto& qp : cell_quadpoints)
        {
            auto c_phi = cb.eval_functions(msh, cl, qp.point(), mindeg, maxdeg);
            matrix_type vec_phi = disk::make_vectorial_matrix(c_phi, DIM);
            mass_mat += qp.weight() * vec_phi.transpose() * vec_phi;

            scalar_type pot = 0.0;
            // u_T
            auto u_T =  x.head(zero_range.size());
            pot += c_phi.dot( u_T);

            size_t number = 0;
            if(mp.diff)
                number = set_cell_number(msh, cl);
            auto potr  = solution.sf(qp.point(), number);

            vector_type dpot  = vec_phi * ruh;
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

        vector_type  true_ddof = projk.compute_vec_cell(msh, cl, solution.df);
        vector_type  comp_ddof = ruh;
        vector_type  diff_ddof = (true_ddof - comp_ddof);
        e.dof_grad += diff_ddof.transpose()* mass_mat * diff_ddof;

        vector_type sigma_dofs = var.stresses.at_cell(msh, cl);
        vector_type sig_true_dof = projk.compute_vec_cell(msh, cl, solution.st);

        vector_type diff_sigma = sig_true_dof - sigma_dofs;
        e.dof_sigma += diff_sigma.transpose() * mass_mat * diff_sigma;
        return;
    }

    template<typename PlasticData>
    bool
    finalize(const size_t iter,
             const mesh_type    & msh,
             const PlasticData  & pst,
             const errors       & error,
             std::ofstream      & tfs,
             const variable_type& var)
    {
        std::string name;
        bool is_finished = true;
        //if(error.gamma < 1.e-14)
        //    name = "GAMMA    ";
        if(error.residue_all < 1.e-8) //with multipliers in the faces
        //if(error.residue < 1.e-8)//without multipliers in the faces
            name = "RESIDUE_ALL";
        else if( iter == m_max_iters - 1 )
            name = "MAX ITERS";
        else if(error.solver_conv < 1.e-15)
            name = "SOLVER   ";
        else
            is_finished = false;

        auto DIM = mesh_type::dimension;
        auto dimension_N_kd = factorial(m_degree_cell + DIM)/ (factorial(m_degree_cell) * factorial(DIM));
        auto dimension_N_kd_1 = factorial(m_degree + DIM -1)/ (factorial(m_degree) * factorial(DIM-1));
        auto F_DOFs = msh.faces_size(); //dimension_N_kd * msh.faces_size();
        auto C_DOFs = msh.cells_size();//dimension_N_kd_1 * msh.cells_size();

        if(is_finished)
        {
            auto h_max = mesh_h_max(msh);
            auto h_min = mesh_h_min(msh);

            disk::surface_areas<scalar_type> surfaces;
            auto initial_marks = disk::marks_yield(msh, var, pst, m_degree, surfaces);

            tfs<< m_degree<<"  " << pst.alpha <<"  " <<h_max <<"   " << h_min << "  ";
            tfs<< F_DOFs<< "  "<< C_DOFs << "  ";
            tfs<< error.dof <<"  "<< error.dof_grad << "  "<<  error.dof_sigma << "  ";
            tfs<< surfaces.liquid <<"  "<< surfaces.rigid<< "  ";
            tfs<< surfaces.rigid_zero << "  "<< surfaces.rigid_mov<< std::endl;

            std::cout<< "/***************************************/" << std::endl;
            std::cout <<"/* break by convergence of  " + name + " */"<< std::endl;
            std::cout<< "/***************************************/" << std::endl;
            std::cout << " L2-norm error, dof         : " << error.dof << std::endl;
            std::cout << " l2-norm error, gamma - Du  : " << error.gamma << std::endl;
            std::cout << " l2-norm error, ierations   : " << error.solver_conv << std::endl;
            std::cout << " Polynomial degree          : "<< m_degree <<std::endl;
            std::cout << " Polynomial cell degree     : "<< m_degree_cell <<std::endl;
            std::cout << " Number of cells            : " << msh.cells_size()<<std::endl;
            std::cout << " h_max : "<< h_max  <<std::endl;
            std::cout << " h_min : "<< h_min  <<std::endl;

            std::cout <<"/* break by convergence of  " + name + " */"<< std::endl;

        }

        return is_finished;
    }


    template<typename BiformsType, typename Solution, typename MeshParameters,
            typename  PlasticData>
    void
    recover_solution(variable_type          &  var,
                      const variable_type   &  old_var,
                      const mesh_type       &  msh,
                      const vector_type     &  X,
                      const std::vector<matrix_type>& gradrec_opers,
                      const BiformsType     &  biforms,
                      const PlasticData     &  pst,
                      const Solution        &  solution,
                      const MeshParameters  &  mp,
                      const scalar_type     &  kappa,
                      const scalar_type     &  on_pst,
                      errors&  e)
    {
        size_t  fbs = fb.size();

        e.set_zero_exact_norms();

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

            matrix_type grec_oper = gradrec_opers.at(cell_id);
            matrix_type loc_mat  = biforms.at(cell_id).cell_lhs * kappa;
            vector_type cell_rhs = biforms.at(cell_id).cell_rhs;
            vector_type pst_rhs;
            vector_type u_local  = var.velocity.at(cell_id);

            //std::cout << "grec_oper size"<< grec_oper.rows() <<" x "<< grec_oper.cols() << std::endl;
            //std::cout << "u_TF size"<< u_local.size() << std::endl;


            if(mp.diff)  //Only for diffusion test
            {
                std::cout << "Do plasticity.compute for diffusion" << std::endl;
                //pst_rhs = plasticity.compute(msh, cl, rec_oper, u_local, pst, solution);
            }
            else        //normal plasticity
            {

                pst_rhs = plasticity.compute_integral(msh, cl, var.stresses);
                //pst_rhs = vector_type::Zero(u_local.size());
            }
            vector_type x = statcond_pst.recover(msh, cl, loc_mat, cell_rhs, xFs,
                                                              on_pst * pst_rhs);


            var.velocity.at(cell_id)  = x;
            vector_type Gt_uT = grec_oper * x;

            size_t number(1);
            if(mp.diff)
                number = disk::set_cell_number(msh, cl);
            if(solution.is_exact)
                compute_exact_error(e, msh, cl, var, old_var, x, Gt_uT, solution, mp);
        }
        if(solution.is_exact)
        {
            e.dof      = std::sqrt(e.dof);
            e.dof_grad      = std::sqrt(e.dof_grad);
            e.fun      = std::sqrt(e.fun);
            e.dfun     = std::sqrt(e.dfun);
            e.dof_sigma = std::sqrt(e.dof_sigma);
        }
        return;
    }

    std::vector<matrix_type>
    precompute_gradient(const mesh_type& msh, const std::string& order,
                        const bool& change_degree)
    {
        std::vector<matrix_type> vec(msh.cells_size());
        gradrec_type gradrec = gradrec_type(m_degree, change_degree);
        size_t i = 0;
        if(order == "low")
        {
            for(auto& cl : msh)
            {
                gradrec.compute_low_order(msh, cl);
                matrix_type stg = gradrec.gradrecoper;
                vec.at(i++) = stg;
            }
        }
        if(order == "high")
        {
            for(auto& cl : msh)
            {
                gradrec.compute(msh, cl);
                matrix_type stg = gradrec.oper;
                vec.at(i++) = stg;
            }
        }

        return vec;
    }

    template<typename LoaderType, typename Solution, typename MeshParameters,
                typename PlasticData>
    auto
    mesh_procedure(const mesh_type      & old_msh,
                    const variable_type & var,
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

        std::cout << "******** check_older_msh in mesh_procedure ******** " << std::endl;
        check_older_msh(old_msh);

        if(mp.binning && imsh !=0)
        {
            disk::binding_meshes<mesh_type> bm(old_msh);
            cells_marks = bm.marker(old_msh, var, pst, m_degree);
            bm.find_geometry(old_msh, cells_marks);
            new_msh = bm.bin(old_msh);
            check_older_msh(new_msh);
            dump_to_matlab(old_msh, "VP_marks.m", cells_marks);
        }
        else if(mp.binNadapt && imsh !=0 &&  (imsh > 2 && imsh !=4 ))
        {
            std::cout << "DO REFINEMENT + AGGLOMERATION" << std::endl;
            disk::stress_based_mesh<mesh_type> sbm(old_msh, levels_vec, pst, mp, imsh);
            std::vector<matrix_type> grec_opers = precompute_gradient(old_msh, "low", mp.hho);

            auto any = sbm.test_marker_yield(old_msh, var, pst, m_degree);
            if(!any)
                return std::make_pair(false, old_msh);

            mesh_type temp_msh = sbm.template refine<LoaderType>(old_msh, m_degree);

            disk::dump_to_matlab(temp_msh, mp.directory  +  "/Refinement_mesh" + info + "_R" +tostr(imsh) +".m");

            std::vector<size_t>().swap(levels_vec);
            levels_vec    = sbm.new_levels;
            ancestors_vec = sbm.ancestors;
            auto cells_marks_temp  = sbm.cells_marks;

            check_older_msh(temp_msh);

            disk::binding_meshes<mesh_type> bm(temp_msh, levels_vec, ancestors_vec);
            bm.find_geometry(temp_msh, sbm.new_binning_marks);
            new_msh = bm.bin(temp_msh);
            levels_vec = bm.new_levels;
            disk::dump_to_matlab(new_msh,  mp.directory  +  "/Agglomaration_mesh" + info + "_R" +tostr(imsh) +".m");

            //cells_marks = bm.cells_marks;
            check_older_msh(new_msh);

            //dump_to_matlab(temp_msh, "VP_marks.m", cells_marks);
        }
        else //(adaptation || imsh == 0)
        {

            disk::stress_based_mesh<mesh_type> sbm(old_msh, levels_vec, pst, mp, imsh);

            if(imsh == 0)
            {
                std::cout << "HERREEE - 1" << std::endl;

                new_msh = sbm.template refine<LoaderType>(old_msh, m_degree);

            }
            else
            {
                if(!mp.mark_all)
                {
                    std::cout << "HERREEE -2" << std::endl;

                    std::vector<matrix_type> grec_opers = precompute_gradient(old_msh, "low", mp.hho);

                    auto any = sbm.template marker<cell_basis_type, face_basis_type, Solution>
                            (old_msh, var, pst, m_degree, solution, grec_opers);

                    if(!any)
                        return std::make_pair(false, old_msh);
                }
                new_msh = sbm.template refine<LoaderType>(old_msh, m_degree);
            }
            std::vector<size_t>().swap(levels_vec);
            levels_vec    = sbm.new_levels;
            ancestors_vec = sbm.ancestors;
            cells_marks   = sbm.cells_marks;
            std::cout << "Size LEVELS  2: "<< levels_vec.size()  << std::endl;
        }
        auto filename =  mp.directory + "/mesh" + mp.summary + "_R" + tostr(imsh) + ".m";
        dump_to_matlab(new_msh, filename, levels_vec, imsh);

        return std::make_pair(true, new_msh);
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

        if(mp.recycle)
        {
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
            //if(mp.recycle)
            //{
            //4. cells_marks
            cells_marks = sizet_vector_type(num_cells);
            disk::load_data(cells_marks, mp, "/marks", "RC", imsh);
        }


        //Updating the mesh
        return new_msh;
    }

    template<typename MeshParameters>
    auto
    solution_procedure( const mesh_type  & new_msh,
                        const mesh_type  & old_msh,
                        const variable_type   & old_var,
                        const MeshParameters& mp,
                        const size_t imsh)
    {
        std::cout << "INSIDE SOLUTION_PROCEDURE" << std::endl;


        typedef std::vector<size_t> sizet_vector_type;
        variable_type var;

        auto num_cells  = new_msh.cells_size();
        auto info  =  mp.summary + "_RC" + tostr(imsh);

        //#if 0
        if(imsh > 0 && mp.recycle)
        {
            std::cout << "RECYCLING SOLUTION" << std::endl;

            //2_Precompute Gradient reconstruction operator
            std::vector<matrix_type> grad_global;
            grad_global =  precompute_gradient(old_msh, "high", mp.hho);

            var.velocity = disk::proj_rec_solution< mesh_type, T,
                                cell_basis_type, cell_quadrature_type,
                                face_basis_type, face_quadrature_type>
            (new_msh, old_msh, grad_global, old_var.velocity, ancestors_vec, m_degree,
                                                                mp.hho);

            var.stresses = disk::project_stresses< mesh_type, T,
                                cell_basis_type, cell_quadrature_type,
                                face_basis_type, face_quadrature_type, tensors_type>
            (new_msh, old_msh, old_var.stresses, ancestors_vec, m_degree, mp.hho);


            #if 0
            std::cout << "SAVING SIGMA_INTER_R" << std::endl;
            auto filename = mp.directory + "/SIGMA_INTER_R" + tostr(imsh)+ ".m";
            disk::quiver_matlab(new_msh, filename, quad_degree, sigma);
            #endif

            return var;
        }
        else //#endif
        {
            std::cout << "* SOLUTION_PROCEDURE ELSE" << std::endl;
            var.template Zero<cell_basis_type,face_basis_type>(new_msh, m_degree, mp.hho);
            return var;
        }
    }

    template< typename MeshParameters>
    auto
    load_solution_procedure(const mesh_type     & msh,
                            const MeshParameters& mp,
                            const size_t imsh)
    {
        std::cout << "INSIDE LOAD_SOLUTION_PROCEDURE_1" << std::endl;

        typedef std::vector<size_t>             sizet_vector_type;
        variable_type  var;

        if(!mp.recycle)
        {
            std::cout << "Not recycling former solution" << std::endl;
            var.template Zero<cell_basis_type,face_basis_type>(msh, m_degree, mp.hho);
            return var;
        }
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
        variable_type  var, old_var;

        auto info_other =  mp.summary_old + "_" + mp.short_mesh_name + ".typ1";
        auto num_cells  = new_msh.cells_size();
        std::cout << "* recycle = "<< mp.recycle << std::endl;
        if(!mp.recycle)
        {
            std::cout << "Not recycling former solution" << std::endl;
            var.template Zero<cell_basis_type,face_basis_type>(new_msh, m_degree, mp.hho);
            return var;
        }
        else
            throw std::invalid_argument("Recycling not defined till now");
        #if 0
        else
        {
            //Load Uh_Th and tensor data
            disk::load_data(old_var.velocity, mp, "/Uh","R", imsh - 1);
            disk::load_data(old_var.stresses, old_msh, mp, "R",imsh);
            assert( old_var.velocity.size()  == old_msh.cells_size());
            assert( old_var.stresses.size()== old_msh.cells_size());

            //5_Precompute Gradient reconstruction operator
            std::vector<matrix_type> grad_global;
            grad_global = precompute_gradient(old_msh);

            var.velocity = disk::proj_rec_solution< mesh_type, T,
                        cell_basis_type, cell_quadrature_type,
                            face_basis_type, face_quadrature_type>
            (new_msh, old_msh, grad_global, old_var.velocity, ancestors_vec, m_degree);

            #if 0
            auto filename = mp.directory + "/SIGMA_INTER_R" + tostr(imsh)+ ".m";
            disk::quiver_matlab(new_msh, filename, quad_degree, sigma);
            #endif
            return var;
        }
        #endif
    }

    template<typename LoaderType, typename Solution, typename PlasticData,
                typename MeshParameters>
    auto
    test_adaptation(variable_type     & var,
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
                    if(mp.load4marker)
                    {
                        disk::load_data(var.stresses, old_msh, mp, "R",imsh);
                        auto pair = mesh_procedure(old_msh, var, loader, solution,
                                                                    pst, mp, imsh);
                        new_msh = pair.second;
                        var.template Zero<cell_basis_type,face_basis_type>(new_msh, m_degree, mp.hho);
                    }
                    else
                    {
                        new_msh = load_mesh_procedure<LoaderType>(old_msh, mp, imsh);
                        var.template Zero<cell_basis_type,face_basis_type>(new_msh, m_degree, mp.hho);
                        #if 0
                        new_msh = load_mesh_procedure<LoaderType>(old_msh, mp, imsh);
                        //if(mp.start_exact)
                            var.template Zero<cell_basis_type,face_basis_type>(new_msh, m_degree, mp.hho);
                        else
                            var  = load_solution_procedure(new_msh, old_msh, mp, imsh);
                        //var = exact_solution_procedure(new_msh, old_msh, mp, imsh, solution.sf);
                        std::cout << "Size LEVELS: "<< levels_vec.size()  << std::endl;
                        #endif
                    }

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
            return std::make_pair(true, old_msh);
    }

    auto
    convergence(const mesh_type& msh,
                     const variable_type & new_var,
                     const variable_type & old_var,
                     const std::vector<matrix_type>& grec_opers,
                     const std::vector<matrix_type>& stab_opers,
                     const T& alpha,
                     const size_t iter,
                     errors & e,
                     T      & old_error,
                     std::ofstream & ifs) //, const timecounter& tc)
    {

        e.set_zero_discrete_norms();

        auto mindeg = 0;
        auto maxdeg = m_degree;
        auto zero_range = cb.range( mindeg,  maxdeg);
        auto cell_range = cb.range( 0,  m_degree_cell);

        auto DIM = mesh_type::dimension;
        auto cbs = zero_range.size();
        auto fbs = fb.size();
        auto vec_cbs = DIM * cbs;

        auto cq_pst = cell_quadrature_type(2 * m_degree_cell);
        auto fq_pst = face_quadrature_type(2 * m_degree);

        for(auto& cl : msh)
        {
            //var.get(cl, sigma);

            auto cl_id = cl.get_id();
            vector_type    new_sigma   = new_var.stresses.at_cell(msh, cl);
            vector_type    old_sigma   = old_var.stresses.at_cell(msh, cl);
            vector_type    sigma_diff  = new_sigma - old_sigma;

            vector_type    new_uh   = new_var.velocity.at(cl_id);
            vector_type    old_uh   = old_var.velocity.at(cl_id);
            vector_type    uh_diff  = new_uh - old_uh;
            vector_type    Guh_diff  = grec_opers.at(cl_id) * uh_diff;

            matrix_type  mass_mat  = matrix_type::Zero(vec_cbs, vec_cbs);

            auto cell_quadpoints   = cq_pst.integrate(msh, cl);
            for (auto& qp : cell_quadpoints)
            {
                auto c_phi = cb.eval_functions(msh, cl, qp.point(), mindeg, maxdeg);
                matrix_type vec_phi = disk::make_vectorial_matrix(c_phi, DIM);
                mass_mat  += qp.weight() * vec_phi.transpose() * vec_phi;
            }

            auto gamma =  plasticity.m_decoupled_var.at_cell(msh, cl);
            vector_type   gamma_Guh = gamma - grec_opers.at(cl_id) * new_uh;
            e.gamma   +=  gamma_Guh.transpose() * mass_mat * gamma_Guh;

            e.sigma  +=  sigma_diff.transpose() * mass_mat * sigma_diff;
            e.grad_u +=  Guh_diff.transpose() * mass_mat * Guh_diff;

            matrix_type  f_mass_mat  = matrix_type::Zero(fbs, fbs);
            auto fcs = faces(msh, cl);
            for(auto& fc: fcs)
            {
                auto face_quadpoints   = fq_pst.integrate(msh, fc);
                for (auto& qp : face_quadpoints)
                {
                    auto f_phi = fb.eval_functions(msh, fc, qp.point());
                    f_mass_mat  += qp.weight() * f_phi * f_phi.transpose();
                }
                vector_type  new_varsigma  = new_var.stresses.at_face(msh, cl, fc);
                vector_type  old_varsigma  = old_var.stresses.at_face(msh, cl, fc);
                vector_type  varsigma_diff = new_varsigma - old_varsigma;

                auto pos = face_position(msh, cl, fc);
                auto num_total_dofs = cell_range.size()
                                                + fbs * number_of_faces(msh, cl);

                matrix_type c_stab_opers = stab_opers.at(cl_id);
                matrix_type S_TF = c_stab_opers.block( 0, pos * num_total_dofs,
                                                            fbs, num_total_dofs);

                vector_type Suh_diff  =  S_TF * uh_diff;

                e.varsigma +=  varsigma_diff.transpose() * f_mass_mat * varsigma_diff;
                e.stab_u   +=  Suh_diff.transpose() * f_mass_mat * Suh_diff;

                auto psi =  plasticity.m_decoupled_var.at_face(msh, cl, fc);
                vector_type   psi_Suh = psi - S_TF * new_uh;
                e.psi   +=  psi_Suh.transpose() * f_mass_mat * psi_Suh;
            }
        }

        auto alpha_2 = alpha * alpha;
        e.residue_all = e.sigma + alpha_2 * e.grad_u
                                            + e.varsigma + alpha_2 * e.stab_u;
        e.residue     = e.sigma + alpha_2 * e.grad_u;

        e.sqrt_discrete_norms();

        e.solver_conv = std::abs(old_error - e.residue_all);
        old_error  = e.residue_all;

        ifs << tostr(iter)<<"  "<< e.residue_all  <<"  "<< e.residue <<"  ";
        ifs << e.grad_u   <<"  "<< e.stab_u <<"  ";
        ifs << e.sigma    <<"  "<< e.gamma  <<"  ";
        ifs << e.varsigma <<"  "<< e.psi    << "  ";
        ifs << e.dof      <<"  "<< e.dof_sigma << " ";
        ifs << e.dof_grad <<" ";
        ifs << e.solver_conv<< std::endl; //<<"  "<< tc << std::endl;

        //for (auto& t : timings)
        //    std::cout << " * " << t.first << ": " << t.second << " seconds." << std::endl;


        std::cout << "L2-norm error, residue_all: " << e.residue_all << std::endl;
        std::cout << "L2-norm error, residue    : " << e.residue << std::endl;
        #if 0
        std::cout << "L2-norm error, gradiente  : " << e.grad_u << std::endl;
        std::cout << "L2-norm error, stability  : " << e.stab_u << std::endl;
        std::cout << "L2-norm error, sigma      : " << e.sigma << std::endl;
        std::cout << "L2-norm error, varsigma   : " << e.varsigma << std::endl;
        std::cout << "L2-norm error, ierations  : " << e.solver_conv << std::endl;
        std::cout << "L2-norm error, dof        : " << e.dof << std::endl;
        std::cout << "L2-norm error, dof_sigma  : " << e.dof_sigma << std::endl;
        #endif
        return;
    }

    template<typename Solution, typename MeshParameters, typename PlasticData>
    void
    test_plasticity(mesh_type           & msh,
                    variable_type       & var,
                    const PlasticData   & pst,
                    const Solution      & solution,
                    const MeshParameters& mp,
                    std::ofstream       & ifs,
                    std::ofstream       & tfs,
                    std::ofstream       & exfs,
                    const size_t          imsh,
                    std::map<std::string, double>& total_timings)
    {
        std::cout << "****************************** BEGIN Refinement No.";
        std::cout << imsh<<"- test_plasticity *********************************" << std::endl;

        scalar_type  old_error(0), kappa(0);
        size_t on_pst;

        timecounter tc, tc_pre;
        std::map<std::string, double> step_timings;


        std::vector<matrix_type> stab_opers;
        std::vector<matrix_type> grec_opers;

        assembly_pst  = assembler_type (msh, m_degree);
        gradrec_pst   = gradrec_type(m_degree, mp.hho);
        stabilization = stab_type(m_degree, mp.hho);
        statcond_pst  = statcond_type(m_degree, mp.hho);

        tc_pre.tic();
        auto biforms  = precompute_bilinear_forms(msh, solution, mp, grec_opers, stab_opers);
        plasticity    = plasticity_type(m_degree, pst, grec_opers, stab_opers, mp.hho);
        tc_pre.toc();

        total_timings["Precomputations"] += tc_pre.to_double();
        step_timings["Precomputations"] += tc_pre.to_double();

        tc.tic();

        for(size_t iter = 0; iter < m_max_iters ; iter++)
        {
            std::cout << "/* ____iter "<< imsh<<"-"<<iter<<"_____ */" << std::endl;
            auto old_var  = var;
            auto is_kappa_diff = false;

            kappa = compute_kappa(kappa, is_kappa_diff, on_pst, pst, mp, solution, iter, imsh);

            timecounter tc_detail;

            tc_detail.tic();
            make_stiff_matrix(msh, kappa, is_kappa_diff, biforms);
            tc_detail.toc();
            total_timings["Stiff matrix"] += tc_detail.to_double();
            step_timings["Stiff matrix"] += tc_detail.to_double();

            tc_detail.tic();
            plasticity.compute_decoupling( msh, var.stresses, var.velocity, step_timings);
            tc_detail.toc();
            total_timings["Plasticity"] += tc_detail.to_double();
            step_timings["Plasticity"] += tc_detail.to_double();

            tc_detail.tic();
            make_rhs_vector(msh, var, biforms, pst, solution, mp, kappa, on_pst, imsh);
            tc_detail.toc();
            total_timings["Rhs vector"] += tc_detail.to_double();
            step_timings["Rhs vector"] += tc_detail.to_double();

            tc_detail.tic();
            vector_type X = solver.solve(assembly_pst.rhs);
            tc_detail.toc();
            total_timings["Solver"] += tc_detail.to_double();
            step_timings["Solver"] += tc_detail.to_double();

            tc_detail.tic();
            errors error;
            recover_solution(var, old_var, msh, X, grec_opers, biforms, pst,
                                        solution, mp, kappa, on_pst, error);
            tc_detail.toc();
            total_timings["Recover solution"] += tc_detail.to_double();
            step_timings["Recover solution"] += tc_detail.to_double();

            tc_detail.tic();
            plasticity.update_multiplier(msh, var.stresses, var.velocity, error);
            tc_detail.toc();
            total_timings["Plasticity - multiplier"] += tc_detail.to_double();
            step_timings["Plasticity - multiplier"] += tc_detail.to_double();

            convergence( msh, var, old_var, grec_opers, stab_opers, pst.alpha,
                                            iter, error, old_error, ifs);//, tc);

            if (finalize(iter, msh, pst, error, tfs, var))
                break;
        }
        tc.toc();

        total_timings["Iterations"] += tc_pre.to_double();
        step_timings["Iterations"] += tc_pre.to_double();

        auto rec_opers = precompute_gradient(msh, "high", mp.hho);

        disk::postprocess(msh, var, grec_opers, rec_opers,
                                 mp, pst, m_degree, imsh, solution);
        disk::execute(exfs, mp, pst, imsh, step_timings);
    }
};
// *****************************************************************************


int main (int argc, char** argv )
{

    using RealType = double;

    int     ch;
    int     degree      = 0;
    size_t  elems_1d    = 10;
    size_t  max_iters   = 500000;
    char    *filename   = nullptr;
    RealType Bingham(0), Load4marker(0);
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
                std::cout << "Degree must be positive. Falling back to 0." << std::endl;
                degree = 0;
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
            case 'l':
                Load4marker = atof(optarg);
                if (Load4marker != 1)
                {
                    std::cout << "Specified Load4marker(1) otherwise fall back to 0." << std::endl;
                    mp.load4marker = false;
                }
                else
                    mp.load4marker = true;
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
        typedef typename mesh_type::scalar_type                 scalar_type;
        typedef typename mesh_type::cell                        cell_type;
        typedef disk::fvca5_mesh_loader<RealType, 2>            loader_type;
        typedef Eigen::Matrix< RealType, 2, Eigen::Dynamic>     tensor_matrix;
        //typedef disk::diffusion<RealType,2>                solution_type;
        typedef disk::square_tuyau<RealType,2>                solution_type;
        typedef variable2<mesh_type>                            variable_type;
        typedef plasticity2_problem<mesh_type, variable_type>   problem_type;

        //std:: string root = "circular/Meeting_24_10/k_"+ tostr(degree) + "/";
        std:: string root = "viscoplasticity/square/Meeting_12_11/k_"+ tostr(degree)+ "/adaptive+WM/";
        //std:: string root = "circular/test/";

        if(!read_parameters(filename, root + "2D_Bi", Bingham, pst, mp, ap,
                                                        max_iters, use_ADMM))
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
        variable_type   var;
        problem_type    pp(var, msh, degree, max_iters, mp, root);

        if(solution.identity == "diffusion")
        {
             mp.diff = true;
             max_iters = 1;
        }

        std::map<std::string, double> timings;

        for( size_t imsh = mp.initial_imsh; imsh < mp.num_remesh + 1; imsh++)
        {
            std::cout <<  pst << std::endl;
            std::cout << "/* Start test adaptation */" << std::endl;
            auto ta_pair = pp.test_adaptation(var, msh, loader, solution, pst, mp, imsh);
            std::cout << "/* End test adaptation */" << std::endl;

            if(!ta_pair.first)
                break;
            msh = ta_pair.second;

            std::ofstream ifs, efs, exfs;
            pp.output_data(ifs, efs, exfs, mp, imsh);


            std::cout <<  pst << std::endl;

            //if(use_ADMM)
            //    pp.test_plasticity_ADMM(msh, var, pst, solution, mp, ap, ifs, efs, exfs, imsh);
            //else if(mp.start_exact)
            //    pp.test_plasticity_exact(msh, var, pst, solution, mp, ifs, efs, exfs, imsh);
            //else
            std::cout << "/* Start test plasticity */" << std::endl;
            pp.test_plasticity(msh, var, pst, solution, mp, ifs, efs, exfs, imsh, timings);

            std::cout << "/* End test plasticity */" << std::endl;

            std::cout <<  pst << std::endl;

            ifs.close();
            efs.close();
            exfs.close();

        }
        disk::make_fig(mp, timings );

        return 0;
    }
    return 0;
};

#if 0
1. Test adapatation
This will be simplified so I dont have to think until about refinement and loading procedures
Take care of uncomment #include "loaders/mesh_adaptation.hpp" in this file and fix problems.
2. post_processing_all
To avoid incompatibilities #include "mesh_adaptation_vpst.hpp" has only MeshParameters. It is mandatory
to solve problems with save_data and others when using the complete version of test_adaptation.
#endif
