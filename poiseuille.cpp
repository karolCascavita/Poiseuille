#include <iostream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <map>
#include "loaders/loader.hpp"
#include "hho/hho.hpp"
#include "hho/hho_pst.hpp"
#include "loaders/plasticity_mesh_adaptation.hpp"
#include "visualization/post_processing_all.hpp"

#include <unsupported/Eigen/SparseExtra>


template<typename T, size_t DIM, typename Storage>
std::vector<dynamic_vector<T>>
startup(const disk::mesh<T,DIM,Storage>& msh, const size_t degree,const size_t imsh)
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


template<typename MeshType, typename Solution,typename TensorsType, typename T>//, typename LoaderType>
bool SOLVE(MeshType & msh,   /* handle to the mesh */
           const Solution&     solution,                /* solution of the problem */
           const disk::plasticity_data<T>&  pst,        /* Material and ALG parameters */
           std::vector<dynamic_vector<T>>&  Uh_Th,      /* storage of velocities*/
           std::vector<TensorsType>      &  tsr_vec,    /* storage of tensors in each cell and in each quad point*/
           const std::string             & directory,
           const size_t degree,                         /* degree of the method */
           T & ALG_error,
           T & old_dof_error,
           size_t iter,
           size_t max_iters,
           std::ofstream     & ifs,
           std::ofstream     & tfs)
{
    typedef MeshType mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>     cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>     face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    typedef dynamic_vector<scalar_type>                 vector_type;
    typedef dynamic_matrix<scalar_type>                 matrix_type;

    //WK: for p-adaptation quad_degree is not the same for all cells
    size_t quad_degree = tsr_vec.at(0).quad_degree;

    disk::gradient_reconstruction_nopre<mesh_type,
                                        cell_basis_type,
                                        cell_quadrature_type,
                                        face_basis_type,
                                        face_quadrature_type> gradrec_nopre(degree, quad_degree);

    disk::diffusion_like_stabilization_nopre<mesh_type,
                                             cell_basis_type,
                                             cell_quadrature_type,
                                             face_basis_type,
                                             face_quadrature_type> stab_nopre(degree);

    disk::diffusion_like_static_condensation_pst<mesh_type,
                                                   cell_basis_type,
                                                   cell_quadrature_type,
                                                   face_basis_type,
                                                   face_quadrature_type> statcond_pst(degree);
    disk::assembler_nopre<mesh_type,
                          face_basis_type,
                          face_quadrature_type> assembler_nopre(msh, degree);

    disk::plasticity<scalar_type,
                        mesh_type,
                        cell_basis_type,
                        cell_quadrature_type,
                        face_basis_type,
                        face_quadrature_type>  plasticity(degree,pst);

    cell_basis_type         cb(degree+1);
    cell_quadrature_type    cq(2*degree+2);

    scalar_type factor, plst_switch(1.);

    if(pst.method == true)
        factor = pst.alpha + pst.mu;
    else
        factor = pst.alpha;

    if(iter == 0)
    {
        factor = pst.mu;
        plst_switch = 0.;
    }
    for (auto& cl : msh)
    {
        gradrec_nopre.compute(msh, cl);
        stab_nopre.compute(msh, cl, gradrec_nopre.oper);
        matrix_type   loc =  factor * (gradrec_nopre.data + stab_nopre.data);
        auto cell_rhs     =  disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, solution.f, degree);
        auto cell_id      =  msh.lookup(cl);
        auto tsr          =  tsr_vec[cell_id];
        plasticity.compute(msh, cl, gradrec_nopre.oper, Uh_Th[cell_id], tsr); //WK: if diffusion changes with domain, it would be mandatory to include mu here, to do it locally
        auto scnp         =  statcond_pst.compute(msh, cl, loc, cell_rhs, plst_switch*plasticity.rhs);
        assembler_nopre.assemble(msh, cl, scnp);
    }
    assembler_nopre.impose_boundary_conditions(msh, solution.sf);
    assembler_nopre.finalize();

    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
    solver.analyzePattern(assembler_nopre.matrix);
    solver.factorize(assembler_nopre.matrix);
    dynamic_vector<scalar_type> X = solver.solve(assembler_nopre.rhs);
    face_basis_type face_basis(degree);
    size_t fbs = face_basis.size();

    scalar_type gamma_error = 0.0;
    scalar_type u_Ruh_error = 0.0;
    scalar_type err_fun  = 0.0;
    scalar_type err_dfun = 0.0;
    scalar_type err_dof  = 0.0;

    disk::projector_nopre<mesh_type,
                        cell_basis_type,
                        cell_quadrature_type,
                        face_basis_type,
                        face_quadrature_type> projk(degree);

    std::ofstream ofs( directory + "/plotnew.dat");
    for (auto& cl : msh)
    {
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces*fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            dynamic_vector<scalar_type> xF = dynamic_vector<scalar_type>::Zero(fbs);
            xF = X.block(face_id * fbs, 0, fbs, 1);
            xFs.block(face_i * fbs, 0, fbs, 1) = xF;
        }

        gradrec_nopre.compute(msh, cl);
        stab_nopre.compute(msh, cl, gradrec_nopre.oper);

        dynamic_matrix<scalar_type> loc = factor*(gradrec_nopre.data + stab_nopre.data);
        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, solution.f, degree);
        auto cell_id  = msh.lookup(cl);
        plasticity.compute(msh, cl, gradrec_nopre.oper, Uh_Th[cell_id], tsr_vec[cell_id]);
        dynamic_vector<scalar_type>   x = statcond_pst.recover(msh, cl, loc, cell_rhs, xFs, plst_switch * plasticity.rhs);

        Uh_Th[cell_id]  = x;

        dynamic_vector<scalar_type> rec(cb.size());
        rec.tail(cb.size()-1) = gradrec_nopre.oper * x;
        rec(0) = x(0);

        //_____________________________________________________________________
        cell_basis_type         cell_basis(degree);
        auto qps_1 = cq.integrate(msh, cl);
        for (auto& qp : qps_1)
        {
            auto phi  = cell_basis.eval_functions(msh, cl, qp.point());
            auto dphi = cb.eval_gradients(msh, cl, qp.point());

            scalar_type pot = 0.0;
            for (size_t i = 0; i < cell_basis.size(); i++)
                pot += phi[i] * x(i);

            vector_type dpot;
            auto col_range  =   cb.range(1,degree+1);
            auto row_range  =   disk::dof_range(0,mesh_type::dimension);
            matrix_type dphi_matrix =   disk::make_gradient_matrix(dphi);
            matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
            matrix_type dphi_rec_uh =   dphi_taken * gradrec_nopre.oper*x;
            dpot = dphi_rec_uh;

            auto potr  = solution.sf(qp.point());
            vector_type dpotr = solution.df(qp.point());
            scalar_type diff = 0.0;
            diff = (pot - potr) * (pot - potr) * qp.weight();
            //std::cout << pot << " " << potr << " " << qp.weight() << " " << diff << std::endl;

            scalar_type d_diff = 0.0;
            d_diff = (dpot - dpotr).dot(dpot - dpotr) * qp.weight();

            err_fun  += diff;
            err_dfun += d_diff;
            auto tp = qp.point();
            for (size_t i = 0; i < mesh_type::dimension; i++)
                ofs << tp[i] << " ";
            ofs << pot << " " << std::abs(pot - solution.sf(tp)) << std::endl;
        }

        dynamic_vector<scalar_type> true_dof = projk.compute_cell(msh, cl, solution.sf);
        dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
        dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
        err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
        //_____________________________________________________________________


        auto test_points = make_test_points(msh, cl);
        for (size_t itp = 0; itp < test_points.size(); itp++)
                //for (auto& qp : qps)
        {
            auto tp   = test_points[itp];
            //auto tp = qp.point();
            auto pot  = 0.;
            auto potR = 0.;

            auto phi = cb.eval_functions(msh, cl, tp);

            for (size_t i = 0; i < cb.size(); i++)
                pot  += phi[i] * rec(i);

            //for (size_t i = 0; i < mesh_type::dimension; i++)
                //ofs << tp[i] << " ";
            //ofs << pot << " " << pot << std::endl;
        }

        cell_quadrature_type    cq2(quad_degree);

        auto qps = cq2.integrate(msh, cl);
        size_t col = 0;
        for (auto& qp : qps)
        {
            auto gamma      =   tsr_vec[cell_id].gamma.col(col);
            auto dphi       =   cb.eval_gradients(msh, cl, qp.point());
            auto col_range  =   cb.range(1,degree+1);
            auto row_range  =   disk::dof_range(0,mesh_type::dimension);

            matrix_type dphi_matrix =   disk::make_gradient_matrix(dphi);
            matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
            matrix_type dphi_rec_uh =   dphi_taken * gradrec_nopre.oper*x;

            gamma_error  += ((gamma - dphi_rec_uh).cwiseProduct(gamma - dphi_rec_uh)).sum();
            col++;
        }

    }

    ofs.close();

    scalar_type TOL_solver = 1.e-10;
    scalar_type TOL = 1.e-6;
    u_Ruh_error     = std::sqrt(u_Ruh_error);
    gamma_error     = std::sqrt(gamma_error);
    err_dof         = std::sqrt(err_dof);
    scalar_type stopping_criterio = std::abs(ALG_error - gamma_error);
    scalar_type error_stop_criterion = std::abs(err_dof - old_dof_error);

    old_dof_error = err_dof;
    ALG_error = gamma_error;

    auto diam = disk::diameter(msh);
    disk::errors<T> er;


    ifs<<to_string(iter) <<"  "<<gamma_error<<" "<< std::sqrt(err_dof) <<"  "
       <<std::sqrt(err_dfun)<<"    "<< std::sqrt(err_fun)<<std::endl;


    auto DIM = mesh_type::dimension;


    if( solution.is_exact)
    {
        std::cout << "L2-norm error, dof:   " << std::sqrt(err_dof) << std::endl;
        std::cout << "L2-norm error, fun:   " << std::sqrt(err_fun) << std::endl;
        std::cout << "L2-norm error, dfun:   " << std::sqrt(err_dfun) << std::endl;
        std::cout << std::endl;
        std::cout << "l2-norm error, L2_norm(dof):" << error_stop_criterion << std::endl;
        std::cout << "l2-norm error, gamma - Du  :" << gamma_error << std::endl;
        std::cout << "l2-norm error, ierations   :" << stopping_criterio << std::endl;

        if( error_stop_criterion <= TOL)
        {
            std::cout << std::endl;
            std::cout<< " /***********************************/" << std::endl;
            std::cout << "/* break by convergence of ER.U-UH */" << std::endl;
            std::cout<< " /***********************************/" << std::endl;
            std::cout << std::endl;

            //tfs<<degree<<"  " << pst.alpha <<"  " << diam<<"    "<< er.u_uh<<"   "<< er.Iu_uh<<"   "<< er.Du_Guh<<"   "<< std::sqrt(err_fun)<<"   " << std::sqrt(err_dof)<<"   "<< std::sqrt(err_dfun)<<std::endl;
            tfs<<degree<<"  " << pst.alpha <<"  " << diam<<"   "<< std::sqrt(err_fun)<<"   " << std::sqrt(err_dof)<<"   "<< std::sqrt(err_dfun)<<std::endl;

            return 1;
        }
        if( iter == max_iters - 1 )
            tfs<<degree<<"  " << pst.alpha <<"  " << diam<<"   "<< std::sqrt(err_fun)<<"   " << std::sqrt(err_dof)<<"   "<< std::sqrt(err_dfun)<<std::endl;

    }
    else
    {
        //L2_error just as an indicative that is not going to inf
        #if 0
        std::cout << "L2-norm error:   " << std::sqrt(err_dof) << std::endl;
        std::cout << "l2-norm error, gamma - Du  :" << gamma_error << std::endl;
        std::cout << "l2-norm error, ierations   :" << stopping_criterio << std::endl;
        #endif
        if(gamma_error < TOL)
        {
            std::cout << std::endl;
            std::cout<< " /***********************************/" << std::endl;
            std::cout << "/* break by convergence of  GAMMA  */" << std::endl;
            std::cout<< " /***********************************/" << std::endl;
            std::cout << std::endl;
            //tfs<<degree<<"  " << pst.alpha <<"  " << diam<<"    "<<gamma_error<<"  " << std::sqrt(err_dof) <<"  " << std::sqrt(err_dfun)<<"  "<< std::sqrt(err_fun)<<std::endl;
            tfs<<degree<<"  " << pst.alpha <<"  " << diam<<"   "<< std::sqrt(err_fun)<<"   " << std::sqrt(err_dof)<<"   "<< std::sqrt(err_dfun)<<std::endl;
            return 1;
        }
        if( iter == max_iters - 1 )
            tfs<<degree<<"  " << pst.alpha <<"  " << diam<<"   "<< std::sqrt(err_fun)<<"   " << std::sqrt(err_dof)<<"   "<< std::sqrt(err_dfun)<<std::endl;

    }

    #if 0
    if(stop_criterio < TOL_solver)
    {
        std::cout << "/* break by TOL_solver */"<<TOL_solver << std::endl;
        return 1;
    }
    #endif

    return 0;

};

template<typename Storage, typename TensorsType, typename T, typename LoaderType>
bool
REFINMENT( disk::mesh<T,2,Storage>  & msh,
           std::vector<size_t>      & levels_vec,
           std::string  & er_steps_name,
           std::string  & er_iters_name,
           const std::string  & directory,
           const TensorsType  & tsr_vec,
           const disk::plasticity_data<T>& pst,
           const bool call_mesher,
           const bool hanging_nodes,
           const size_t elems_1d,
           const size_t degree,                         /* degree of the method */
           const size_t imsh,
           const LoaderType& loader)
{
    typedef  disk::mesh<T,2,Storage> mesh_type;
    auto method_str =  (pst.method)? to_string(1): to_string(2);
    auto hanging_on = hanging_nodes? 1 : 0;
    auto info_1 =  "_n" + to_string(elems_1d) + "_g"+ to_string(hanging_on) ;
    auto info_2 =  "_m" + method_str  +  "_a" +  to_string(int(pst.alpha));
    if(call_mesher && imsh > 0)
    {
        bool do_refinement;
        auto info  =  info_1 + "_rmc" + to_string(imsh) + info_2;

        disk::stress_based_mesh<mesh_type>   sbm(msh, levels_vec, tsr_vec, pst.yield, do_refinement, hanging_nodes,directory,info);
        std::cout << "Do Refinement No."<< imsh <<" ?  "<<do_refinement << std::endl;
        if(!do_refinement)
            return false;
        std::cout << "*************BEGIN Refinement No. "<< imsh<<"****************" << std::endl;
        //auto info  =  info_1 + "_rmc" + to_string(imsh) + info_2;
        sbm.template re_populate_mesh<LoaderType>(msh, tsr_vec, pst.yield, directory, info, degree);
        levels_vec = sbm.new_levels;
        tsr_vec    = zero_tensor_vector(msh, degree);
    }

    auto other_info =  info_1 + "_rm" + to_string(imsh) + info_2;
    auto error_info =  info_1 + info_2;
    er_steps_name   =  directory + "/error_by_step_adapt"+ error_info + ".dat";
    er_iters_name   =  directory + "/error" + other_info +  ".dat";

    dump_to_matlab(msh, directory + "/mesh"+ other_info + ".m");

    return true;
};



int main (int argc, char** argv )
{

    using RealType = double;

    char    *filename   = nullptr;
    int     degree      = 1;
    size_t  elems_1d    = 10;
    size_t  max_iters   = 1;
    bool    call_mesher = false;
    bool    hanging_nodes = false;

    disk::plasticity_data<RealType> pst;
    int     ch, method_name(1),refine_on(0), hanging_on(0);

    while ( (ch = getopt(argc, argv, "k:n:h:v:a:B:i:m:g:")) != -1 )
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
                break;

            case 'n':
                elems_1d = atoi(optarg);
                if (elems_1d < 0)
                {
                    std::cout << "Num of elems must be positive. Falling back to 8." << std::endl;
                    elems_1d = 10;
                }
                break;
            case 'v':
                pst.mu = atof(optarg);
                if (pst.mu <= 0)
                {
                    std::cout << "Viscosity must be bigger than 0. Falling back to 1." << std::endl;
                    pst.mu = 1.0;
                }
                break;
            case 'a':
                pst.alpha = atof(optarg);
                if (pst.alpha <= 0)
                {
                    std::cout << "Augmentation parameter must be bigger than 0. Falling back to 1." << std::endl;
                    pst.alpha = 1.;
                }
                break;
            case 'B':
                pst.Bn = atof(optarg);
                if (pst.Bn < 0)
                {
                    std::cout << "Bingham number must be positive. Falling back to 0. Solving diffusion case." << std::endl;
                    pst.Bn = 0.0;
                    max_iters    = 1;
                }
                break;
            case 'm':
                method_name = atoi(optarg);
                if( method_name == 1)
                    pst.method = true;
                else
                {
                    if ( method_name == 2)
                        pst.method = false;
                    else
                    {
                        std::cout << "Method must be: Glowinski (1) or Saramito (2). Coming back to Glowinski." << std::endl;
                        pst.method = true;
                    }
                }
                break;
            case 'h':
                refine_on = atoi(optarg);
                    if( refine_on == 1)
                        call_mesher = true;
                    else
                    {
                        if ( refine_on == 0)
                            call_mesher = false;
                        else
                        {
                            std::cout << "Refinement: Yes (1) or No (0). Coming back to 0." << std::endl;
                            call_mesher = false;
                        }
                    }
                    break;
            case 'g':
                hanging_on = atoi(optarg);
                    if( refine_on == 1)
                        hanging_nodes = true;
                    else
                    {
                        if ( refine_on == 0)
                                    hanging_nodes = false;
                        else
                        {
                            std::cout << "Refinement: Yes (1) or No (0). Coming back to 0." << std::endl;
                            hanging_nodes = false;
                        }
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

    if (argc == 0)
    {
        std::cout << "Running 1D test simulation" << std::endl;

        typedef disk::generic_mesh<RealType, 1>                 mesh_type;
        typedef disk::uniform_mesh_loader<RealType,1>           loader_type;
        typedef typename mesh_type::scalar_type                 scalar_type;
        typedef disk::poiseuille<RealType,1>                    problem_type;

        pst.yield = pst.Bn * pst.f * pst.Lref;

        mesh_type       msh;
        loader_type     loader(0,1,elems_1d);
        problem_type    solution(pst);
        loader.populate_mesh(msh);
        std::vector<size_t> levels_vec(msh.cells_size());

        auto   dir_name   = disk::directory<scalar_type, 1>(pst.Bn, solution.name);
        size_t num_remesh = 10*refine_on;
        auto   tsr_vec    = zero_tensor_vector(msh, degree);

        for( size_t imsh = 0; imsh <= num_remesh; imsh ++)
        {

            if(call_mesher & imsh > 0)
            {
                bool do_refinement;
                disk::stress_based_mesh<mesh_type> sbm(msh,levels_vec,tsr_vec,pst.yield, do_refinement);
                std::cout << "do_refinement ? "<< do_refinement << std::endl;
                if(!do_refinement)
                    break;
                std::cout << "************* Refinement No. " << imsh << " ****************" << std::endl;
                sbm.template re_populate_mesh<loader_type>(msh, tsr_vec, pst.yield, dir_name, degree);
                tsr_vec    = zero_tensor_vector(msh, degree);
            }
            auto msh_str = to_string(imsh);
            auto Uh_Th   = startup(msh,  degree, imsh);
            std::ofstream       ifs(dir_name + "/error_msh" + msh_str + ".dat");
            std::ofstream       tfs(dir_name + "/errors.dat",std::ios::app);
            if (!ifs.is_open())
              std::cout << "Error opening file";
              if (!tfs.is_open())
                std::cout << "Error opening tfs";

            RealType gamma_error(0.),L2_error(0.);
            for(size_t iter = 0; iter < max_iters ; iter++)
            {
                std::cout << "/* ___________________iter "<<iter<<"___________________ */" << std::endl;
                bool brk = SOLVE(msh, solution, pst, Uh_Th, tsr_vec, dir_name, degree, gamma_error, L2_error,iter, max_iters, ifs,tfs);
                if( brk == true)
                    break;
            }
            ifs.close();
            tfs.close();
        }
        return 0;
    }

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
        typedef disk::tuyau<RealType,2>                problem_type;

        pst.yield = 0.5*pst.Bn * pst.f * pst.Lref;

        mesh_type    msh;
        loader_type  loader;
        problem_type solution(pst);

        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);
        std::vector<size_t> levels_vec(msh.cells_size());
        size_t num_remesh = 3* refine_on;

        auto   dir_name   = disk::directory<scalar_type, 2>(pst.Bn, solution.name);
        auto   tsr_vec    = zero_tensor_vector(msh, degree);
        for( size_t imsh = 0; imsh < num_remesh + 1; imsh++)
        {
            std::string er_iters_name, er_steps_name;

            bool do_refinement = REFINMENT(msh, levels_vec, er_steps_name, er_iters_name,
                                                         dir_name , tsr_vec, pst,  call_mesher,
                                                         hanging_nodes, elems_1d, degree, imsh,loader);

            std::ofstream       ifs(er_steps_name);
            std::ofstream       tfs(er_iters_name, std::ios::app);
            if (!ifs.is_open())
                std::cout << "Error opening file"<<std::endl;
            if (!tfs.is_open())
                std::cout << "Error opening tfs";


            RealType gamma_error(0.), L2_error(0.);
            auto Uh_Th   = startup(msh, degree, imsh);

            for(size_t iter = 0; iter < max_iters ; iter++)
            {
                std::cout << "/* ____iter "<<iter<<"_____ */" << std::endl;
                bool brk = SOLVE(msh, solution, pst, Uh_Th, tsr_vec, dir_name, degree, gamma_error, L2_error,iter, max_iters, ifs,tfs);
                if( brk == true)
                    break;
            }
            std::string msh_str =  to_string(imsh);
            auto er = disk::plasticity_post_processing(msh, solution.sf,solution.df,tsr_vec,Uh_Th,pst,dir_name, degree, ifs, msh_str, hanging_on, elems_1d);
            std::cout << "************* END Refinement No. "<< imsh <<"****************"<<  std::endl;

            ifs.close();
            tfs.close();
        }
        return 0;
    }
    return 0;
};
