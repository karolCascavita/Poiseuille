#include <iostream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <map>
#include "loaders/loader.hpp"
#include "hho/hho.hpp"
#include "hho/hho_pst.hpp"
#include "visualization/post_processing_all.hpp"

#include <unsupported/Eigen/SparseExtra>


template<typename T, size_t DIM, typename Storage>
std::vector<dynamic_vector<T>>
startup(const disk::mesh<T,DIM,Storage>& msh, const T Bi, std::string& dir_name,const size_t degree)
{

    typedef disk::mesh<T,DIM,Storage>              mesh_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, typename mesh_type::cell>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, typename mesh_type::face>    face_basis_type;

    T stop_error(0.);

    std::ostringstream      var_Bn, var_DIM;
    var_Bn  << 10.*Bi;
    var_DIM << mesh_type::dimension;
    std::string Bitostr     =   var_Bn.str();
    std::string DIMtostr    =   var_DIM.str();
    dir_name    =   DIMtostr + "D_Bi" + Bitostr;
    std::cout << "dir_name"<< dir_name << std::endl;
    size_t  i =  0;

    std::vector<dynamic_vector<T>> U(msh.cells_size());

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

template<typename MeshType, typename Function, typename Solution, typename Gradient,
                                            typename TensorsType>
bool test_diffusion (MeshType& msh,       /* handle to the mesh */
           const Function&     load,            /* rhs */
           const Solution&     solution,        /* solution of the problem */
           const Gradient&     gradient,
           const size_t        degree,             /* degree of the method */
           const disk::plasticity_data<typename MeshType::scalar_type>&  pst,   /* Material and ALG parameters */
           std::vector<dynamic_vector<typename MeshType::scalar_type>> & Uh_Th,    /* storage of velocities*/
           std::vector<TensorsType>& tsr_vec,
           typename MeshType::scalar_type& stop_error,
           int iter,
           //std::ofstream& its,
           std::ofstream& efs,
           const std::string & directory,
            typename MeshType::scalar_type& gamma_error)   /* storage of tensor sigma  in each cell and in each quad point*/
{

    //typedef disk::mesh<T,DIM,Storage>    mesh_type;
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

    disk::gradient_reconstruction_nopre<mesh_type,
                                        cell_basis_type,
                                        cell_quadrature_type,
                                        face_basis_type,
                                        face_quadrature_type> gradrec_nopre(degree);

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
    disk::assembler_pst<mesh_type,
                          face_basis_type,
                          face_quadrature_type> assembler_pst(msh, degree);

    disk::plasticity<scalar_type,
                        mesh_type,
                        cell_basis_type,
                        cell_quadrature_type,
                        face_basis_type>  plasticity(degree,pst);

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
        matrix_type   loc = factor * (gradrec_nopre.data + stab_nopre.data);
        auto cell_rhs     =  disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, degree);
        auto cell_id      =  msh.lookup(cl);
        auto tsr          =  tsr_vec[cell_id];
        plasticity.compute(msh, cl, gradrec_nopre.oper, Uh_Th[cell_id],tsr); //WK: if diffusion changes with domain, it would be mandatory to include mu here, to do it locally
        auto scnp         =  statcond_pst.compute(msh, cl, loc, cell_rhs, plst_switch*plasticity.rhs);
        assembler_pst.assemble(msh, cl, scnp);
    }
    assembler_pst.impose_boundary_conditions(msh, solution, pst);
    assembler_pst.finalize();

    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
    solver.analyzePattern(assembler_pst.matrix);
    solver.factorize(assembler_pst.matrix);
    dynamic_vector<scalar_type> X = solver.solve(assembler_pst.rhs);
    face_basis_type face_basis(degree);
    size_t fbs = face_basis.size();

    scalar_type diam        = 0.0;
    gamma_error = 0.0;
    scalar_type gamma_error_2 = 0.0;
    scalar_type u_Ruh_error = 0.0;

    std::ofstream ofs( directory + "/plotnew.dat");
    for (auto& cl : msh)
    {
        diam = std::max(diameter(msh, cl), diam);
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
        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, degree);
        auto cell_id  = msh.lookup(cl);
        plasticity.compute(msh, cl, gradrec_nopre.oper, Uh_Th[cell_id], tsr_vec[cell_id]);
        dynamic_vector<scalar_type>   x = statcond_pst.recover(msh, cl, loc, cell_rhs, xFs, plst_switch*plasticity.rhs);

        Uh_Th[cell_id]  = x;

        dynamic_vector<scalar_type> rec(cb.size());
        rec.tail(cb.size()-1) = gradrec_nopre.oper * x;
        rec(0) = x(0);

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

            for (size_t i = 0; i < mesh_type::dimension; i++)
                ofs << tp[i] << " ";
            ofs << pot << " " << pot << std::endl;
        }

        auto qps = cq.integrate(msh, cl);
        size_t col = 0;

        for (auto& qp : qps)
        {
            auto cell_id    =   msh.lookup(cl);
            auto gamma      =   tsr_vec[cell_id].gamma.col(col);
            auto dphi       =   cb.eval_gradients(msh, cl, qp.point());
            auto phi        =   cb.eval_functions(msh, cl, qp.point());
            auto col_range  =   cb.range(1,degree+1);
            auto row_range  =   disk::dof_range(0,mesh_type::dimension);

            matrix_type dphi_matrix =   disk::make_gradient_matrix(dphi);
            matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
            matrix_type dphi_rec_uh =   dphi_taken * gradrec_nopre.oper*x;

            gamma_error_2  += ((gamma - dphi_rec_uh).cwiseProduct(gamma - dphi_rec_uh)).sum();
            gamma_error  += ((gamma - dphi_rec_uh).cwiseProduct(gamma - dphi_rec_uh)).sum();
            col++;
        }

    }

    ofs.close();

    scalar_type TOL_solver = 1.e-12;
    scalar_type TOL = 1.e-12;
    u_Ruh_error = std::sqrt(u_Ruh_error);
    gamma_error_2 = std::sqrt(gamma_error_2);
    gamma_error = std::sqrt(gamma_error);
    scalar_type stop_criterio = std::abs(stop_error - gamma_error);
    std::cout << "l2-norm error, gamma - Du:" << gamma_error << std::endl;
    std::cout << "l2-norm error_2, gamma - Du:" << gamma_error_2 << std::endl;
    std::cout << "l2-norm error, iterations:" << stop_criterio << std::endl;
    //its<< iter << " "<<gamma_error << std::endl;

    if(gamma_error < TOL)
    {
        std::cout << "/* break by convergence of gamma */" << std::endl;
        return 1;
    }
    #if 0
    if(stop_criterio < TOL_solver)
    {
        std::cout << "/* break by TOL_solver */"<<TOL_solver << std::endl;
        return 1;
    }
    #endif
    stop_error = gamma_error;

    disk::plasticity_post_processing(msh,solution,gradient,tsr_vec,Uh_Th,pst, directory, degree, efs,iter);


    return 0;

};






int main (int argc, char** argv )
{

    using RealType = double;

    char    *filename   = nullptr;
    int     degree      = 1;
    int     elems_1d    = 10;
    int     ch, method_name,refine_on;
    size_t  max_iters   = 1;
    bool    do_remesh = false;
    std::string Error_convg  =  "/Error_convergence.dat";
    std::string Error_sol    =  "/Error_solution.dat";

    disk::plasticity_data<RealType> pst;

    while ( (ch = getopt(argc, argv, "k:n:h:v:a:B:i:m:")) != -1 )
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
                        do_remesh = true;
                    else
                    {
                        if ( refine_on == 0)
                            do_remesh = false;
                        else
                        {
                            std::cout << "Refinement: Yes (1) or No (0). Coming back to 0." << std::endl;
                            do_remesh = false;
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

        mesh_type   msh;

        loader_type loader(0,1,elems_1d);
        loader.populate_mesh(msh);

        auto f =  [](const point<RealType, mesh_type::dimension>& p) -> auto { return 1.; };
        auto sf = [](const point<RealType, mesh_type::dimension>& p, const disk::plasticity_data<RealType>& pst)
                    -> auto
                    {
                        RealType ret, fvalue(1.0), xO(0.), xL(pst.Lref);
                        RealType  xc = pst.Lref/2.0 - 0.5*pst.Bn*pst.Lref*fvalue;
                        if( p.x() >= xO & p.x() <  xc)
                            ret = (xc*xc - (p.x() - xc)*(p.x() - xc));
                        if( p.x() >= xc & p.x() <  (xL- xc))
                            ret = xc*xc;
                        if( p.x() <= xL & p.x() >= (xL - xc))
                            ret = (xc*xc - (p.x() - (pst.Lref - xc))*(p.x() - (pst.Lref -xc)));
                        return ret*(0.5*fvalue/pst.mu);
                    };
        auto df = [](const point<RealType, mesh_type::dimension>& p, const disk::plasticity_data<RealType>& pst)
                    -> auto
                    {
                        dynamic_vector<scalar_type> ret   = dynamic_vector<scalar_type>::Zero(1);
                        RealType  fvalue(1.), xO(0.), xL(pst.Lref);
                        RealType  xc = pst.Lref/2.0 - 0.5*pst.Bn*pst.Lref*fvalue;

                        if( p.x() >= xO & p.x() <  xc)
                            ret(0) = - 2.*(p.x() - xc);
                        if( p.x() >= xc & p.x() <  (xL- xc))
                            ret(0) = 0.;
                        if( p.x() <= xL & p.x() >= (xL - xc))
                            ret(0) = - 2.*(p.x() - (pst.Lref - xc));
                        ret = (0.5*fvalue/pst.mu)*ret;
                        return ret;
                    };

        std::string         dir_name;
        std::ofstream       its(dir_name + Error_convg);
        std::ofstream       efs(dir_name + Error_sol);
        if (!its.is_open())
          std::cout << "Error opening file";

        size_t num_remesh = 2*refine_on;
        auto   tsr_vec    = zero_tensor_vector(msh, degree);

        for( size_t imsh = 0; imsh < num_remesh + 1; imsh ++)
        {
            pst.yield = pst.Bn * pst.f * pst.Lref;

            if(do_remesh & imsh > 0)
            {
                disk::stress_based_mesh<mesh_type> sbm(msh,tsr_vec,pst.yield);
                sbm.re_populate_mesh(msh, tsr_vec, pst.yield, degree);
                tsr_vec    = zero_tensor_vector(msh, degree);
                std::cout << "Refine No."<< imsh << std::endl;

            }
            auto Uh_Th = startup(msh, pst.Bn, dir_name, degree);
            RealType stop_error(0.), gamma_error(0.);
            for(int iter = 0; iter < max_iters ; iter++)
            {
                //bool brk = test_diffusion(msh, f, sf, df, degree, pst, Uh_Th, tsr_vec, stop_error,iter,its,efs,dir_name, error_gamma);

                bool brk = test_diffusion(msh, f, sf, df, degree, pst, Uh_Th, tsr_vec, stop_error,iter,efs,dir_name, gamma_error);
                its<< iter << " "<<gamma_error << std::endl;
                if( brk == true)
                    break;
            }
        }

        its.close();
        efs.close();
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

        mesh_type   msh;
        loader_type loader;

        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return 1.;
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p, const disk::plasticity_data<RealType>& pst) -> auto {
            return 0.;
        };
        auto df = [](const point<RealType, mesh_type::dimension>& p, const disk::plasticity_data<RealType>& pst) -> auto {
            return dynamic_vector<scalar_type>::Zero(2);
        };

        std::string         dir_name;
        std::string         dir_error_sol = dir_name + "/Error_solution.dat";
        //std::ofstream       ifs(dir_name + Error_convg);

        std::ofstream ifs("2D_Bi4/Error_convg.dat");
        if (!ifs.is_open())
          std::cout << "Error opening file"<<std::endl;

        std::ofstream       efs(dir_error_sol);
        size_t num_remesh = 1*refine_on;
        auto   tsr_vec    = zero_tensor_vector(msh, degree);
        for( size_t imsh = 0; imsh < num_remesh + 1; imsh ++)
        {
            pst.yield = 0.5*pst.Bn * pst.f * pst.Lref;

            if(do_remesh & imsh > 0)
            {
                std::cout << "************* Refinement No.****************"<< imsh << std::endl;
                std::cout << "============================================"<< imsh << std::endl;

                disk::stress_based_mesh<mesh_type>   sbm(msh,tsr_vec,pst.yield, false);

                sbm.re_populate_mesh<loader_type>(msh, tsr_vec, pst.yield, degree);
                tsr_vec    = zero_tensor_vector(msh, degree);
                std::cout << "============================================"<< imsh << std::endl;

                dump_to_matlab(msh,"sub_mesh_1.m");
            }

            auto Uh_Th = startup(msh, pst.Bn, dir_name, degree);
            RealType stop_error(0.), gamma_error(0.);;
            for(int iter = 0; iter < max_iters ; iter++)
            {
                std::cout << "/* ____iter "<<iter<<"_____ */" << std::endl;
                //bool brk = test_diffusion(msh, f, sf, df, degree, pst, Uh_Th, tsr_vec, stop_error,iter,its,efs,dir_name);

                bool brk = test_diffusion(msh, f, sf, df, degree, pst, Uh_Th, tsr_vec, stop_error,iter,efs,dir_name,gamma_error);
                ifs<< iter << " "<<gamma_error << std::endl;

                if( brk == true)
                    break;
            }
        }
        ifs.close();
        efs.close();
        return 0;
    }
    return 0;
};
