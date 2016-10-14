#include <iostream>
#include <regex>
#include <unistd.h>

#include <map>
#include "loaders/loader.hpp"
#include "hho/hho.hpp"

#include <unsupported/Eigen/SparseExtra>

template<typename MeshType, typename Function, typename Solution, typename PlasticData,
                                            typename StorageTensorQuadMatrixType>
void test_diffusion (MeshType & msh,       /* handle to the mesh */
           const Function& load,            /* rhs */
           const Solution& solution,        /* solution of the problem */
           const size_t degree,             /* degree of the method */
           const PlasticData&  plst_data,   /* Material and ALG parameters */
           std::vector<dynamic_vector<typename MeshType::scalar_type>> & Uh_Th,    /* storage of velocities*/
           StorageTensorQuadMatrixType   & sigma_Th,
           size_t coef_init)   /* storage of tensor sigma  in each cell and in each quad point*/
{
    typedef MeshType    mesh_type;
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

    disk::diffusion_like_static_condensation_nopre<mesh_type,
                                                   cell_basis_type,
                                                   cell_quadrature_type,
                                                   face_basis_type,
                                                   face_quadrature_type> statcond_nopre(degree);

    disk::assembler_nopre<mesh_type,
                          face_basis_type,
                          face_quadrature_type> assembler_nopre(msh, degree);

    disk::plasticity<scalar_type,
                        mesh_type,
                        cell_basis_type,
                        cell_quadrature_type,
                        face_basis_type>  plasticity(degree,
                                                     plst_data.Lref,
                                                     plst_data.Vref,
                                                     plst_data.Bn,
                                                     plst_data.mu,
                                                     plst_data.f,
                                                     plst_data.alpha,
                                                     plst_data.method);

    cell_basis_type         cb(degree+1);
    cell_quadrature_type    cq(2*degree+2);


    char numstr[4];
    std::string ext = ".mtx";
    std::string result;
    std::string name;

    double factor, plst_switch;

    if(coef_init = 0)
    {
        factor = plst_data.mu;
        plst_switch = 0.;
    }
    else
    {
        plst_switch = 1.;
        if(plst_data.method == true)
            factor = plst_data.alpha + plst_data.mu;
        else
        {
            factor = plst_data.alpha;
        }
    }


    for (auto& cl : msh)
    {
        gradrec_nopre.compute(msh, cl);
        stab_nopre.compute(msh, cl, gradrec_nopre.oper);
        matrix_type loc = factor*(gradrec_nopre.data + stab_nopre.data);
        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, degree);

        auto cell_id  = msh.lookup(cl);
        auto local_sigma_storage =  sigma_Th[cell_id];

        plasticity.compute(msh, cl, gradrec_nopre.oper, Uh_Th[cell_id], local_sigma_storage);
        //WK: if diffusion changes with domain, it would be mandatory to include mu here, to do it locally

        auto scnp = statcond_nopre.compute_Bingham(msh, cl, loc, cell_rhs, plst_switch*plasticity.rhs);
        assembler_nopre.assemble(msh, cl, scnp);
    }

    assembler_nopre.impose_boundary_conditions(msh, solution);
    assembler_nopre.finalize();


    //WK: definirlo bien despues, dejarlo como Matteo/
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
    solver.analyzePattern(assembler_nopre.matrix);
    solver.factorize(assembler_nopre.matrix);
    dynamic_vector<scalar_type> X = solver.solve(assembler_nopre.rhs);
    face_basis_type face_basis(degree);
    size_t fbs = face_basis.size();

    scalar_type diam = 0.0;
    scalar_type err_dof_k = 0.0;
    scalar_type err_fun_k = 0.0;
    scalar_type err_dof_kp = 0.0;
    scalar_type err_fun_kp = 0.0;
    std::ofstream ofs("plotnew.dat");
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
        plasticity.compute(msh, cl, gradrec_nopre.oper, Uh_Th[cell_id], sigma_Th[cell_id]);
        dynamic_vector<scalar_type>   x = statcond_nopre.recover_Bingham(msh, cl, loc, cell_rhs, xFs, plst_switch*plasticity.rhs);
        //dynamic_vector<scalar_type>   x = statcond_nopre.recover(msh, cl, loc, cell_rhs, xFs);

        Uh_Th[cell_id]  = x;

        dynamic_vector<scalar_type> rec(cb.size());
        rec.tail(cb.size()-1) = gradrec_nopre.oper * x;
        rec(0) = x(0);

        //rec = (plst_data.Lref*plst_data.Lref /plst_data.mu) * rec;
        auto test_points = make_test_points(msh, cl);
        for (size_t itp = 0; itp < test_points.size(); itp++)
                //for (auto& qp : qps)
        {
            auto tp = test_points[itp];
            //auto tp = qp.point();
            auto pot = 0.;

            auto phi = cb.eval_functions(msh, cl, tp);

            for (size_t i = 0; i < cb.size(); i++)
                pot += phi[i] * rec(i);

            for (size_t i = 0; i < MeshType::dimension; i++)
                ofs << tp[i] << " ";
            ofs << pot << " " << std::abs(pot - solution(tp)) << std::endl;
        }

        auto qps = cq.integrate(msh, cl);

        for (auto& qp : qps)
        {
            auto phi = cb.eval_functions(msh, cl, qp.point());

            scalar_type potk = 0.0;
            for (size_t i = 0; i < cb.range(0,degree).size(); i++)
                potk += phi[i] * rec(i);

            scalar_type potkp = potk;
            for (size_t i = cb.range(0,degree).size(); i < cb.size(); i++)
                potkp += phi[i] * rec(i);

            auto potr = solution(qp.point());

            scalar_type diffk = 0.0;
            scalar_type diffkp = 0.0;
            diffk = (potk - potr) * (potk - potr) * qp.weight();
            diffkp = (potkp - potr) * (potkp - potr) * qp.weight();

            err_fun_k += diffk;
            err_fun_kp += diffkp;
        }

            //err_dof_k += compute_L2_error(dld, degree, solution, x);
            //err_dof_kp += compute_L2_error(dld, degree+1, solution, x);

    }

    ofs.close();

    std::cout << "Mesh diameter: " << diam << std::endl;
    std::cout << "L2-norm error, dof, K:   " << std::sqrt(err_dof_k) << std::endl;
    std::cout << "L2-norm error, fun, K:   " << std::sqrt(err_fun_k) << std::endl;
    std::cout << "L2-norm error, dof, K+1: " << std::sqrt(err_dof_kp) << std::endl;
    std::cout << "L2-norm error, fun, K+1: " << std::sqrt(err_fun_kp) << std::endl;

};



template<typename T>
struct plastic_data
{
    plastic_data(): Lref(1.0), Vref(1.0), Bn(0.1), mu(1.0), alpha(1.), f(1.),method(true)
    {}
    T f;
    T Lref;                 /* Charactetistic length */
    T Vref;                 /* Reference velocity */
    T Bn;                   /* Bingham number */
    T mu;                   /* viscosity */
    T alpha;
    bool   method;

};

template<typename MeshType, typename TensorQuadType, typename CellQuadType>
void
initialization(const MeshType& msh,
                std::vector<TensorQuadType> & sigma_Th,
                std::vector<dynamic_vector<MeshType::scalar_type>> & Uh_Th,
                size_t degree)
{

    typedef MeshType                                   mesh_type;
    typedef TensorQuadType                             tensor_quad_matrix_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, typename mesh_type::cell>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, typename mesh_type::face>    face_basis_type;

    size_t  i =  0;

    for (auto& cl : msh)
    {
        auto fcs = faces(msh, cl);
        auto num_faces       = fcs.size();
        auto num_cell_dofs   = cell_basis_type(degree).size();
        auto num_face_dofs   = face_basis_type(degree).size();
        disk::dofspace_ranges dsr(num_cell_dofs, num_face_dofs, num_faces);

        auto cell_quad      = CellQuadType(2*degree+2);
        // WK: degree of quadrature assumed as for (Dru,Dru)
        // Dejarlo aqui por si camba el grado para cda celda
        auto cell_quadpoints = cell_quad.integrate(msh, cl);

        sigma_Th[i] =  tensor_quad_matrix_type::Zero(MeshType::dimension, cell_quadpoints.size());
        Uh_Th[i]    =  dynamic_vector<typename mesh_type::scalar_type>::Zero(dsr.total_size());
        i++;
    }
};


int main (int argc, char** argv )
{

    using RealType = double;

    char    *filename   = nullptr;
    int     degree      = 1;
    int     elems_1d    = 10;
    int     ch, method_name;
    size_t  max_iters = 1;
    plastic_data<RealType> plst_data;


    while ( (ch = getopt(argc, argv, "k:n:u:l:v:a:B:i:m:")) != -1 )
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
            case 'u':
                plst_data.Vref = atof(optarg);
                if (plst_data.Vref < 0)
                {
                    std::cout << "Velocity norm must be positive. Falling back to 1." << std::endl;
                    plst_data.Vref = 1.0;
                }
                break;
            case 'l':
                plst_data.Lref = atof(optarg);
                if (plst_data.Lref < 0)
                {
                    std::cout << "Length must be positive. Falling back to 1." << std::endl;
                    plst_data.Lref = 1.0;
                }
                break;
            case 'v':
                plst_data.mu = atof(optarg);
                if (plst_data.mu < 0)
                {
                    std::cout << "Viscosity must be positive. Falling back to 1." << std::endl;
                    plst_data.mu = 1.0;
                }
                break;
            case 'a':
                plst_data.alpha = atof(optarg);
                if (plst_data.alpha < 0)
                {
                    std::cout << "Augmentation parameter must be positive. Falling back to 0.1" << std::endl;
                    plst_data.alpha = 0.1;
                }
                break;
            case 'B':
                plst_data.Bn = atof(optarg);
                if (plst_data.Bn < 0)
                {
                    std::cout << "Bingham number must be positive. Falling back to 0." << std::endl;
                    plst_data.Bn = 0.0;
                }
                break;
            case 'm':
                method_name = atoi(optarg);
                if( method_name == 1)
                {
                    plst_data.method = true;
                }
                else
                {
                    if ( method_name == 2)
                    {
                        plst_data.method = false;
                    }
                    else
                    {
                        std::cout << "Method must be: Glowinski (1) or Saramito (2). Coming back to Glowinski." << std::endl;
                        plst_data.method = true;
                    }
                }
                break;
            case 'h':
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
        typedef Eigen::Matrix< RealType, 1, Eigen::Dynamic>     tensor_quad_matrix_type;
        typedef Eigen::Matrix< RealType, 1, 1 >                 tensor_type; //Esto cambia para 2D y 3D
        typedef typename mesh_type::cell                        cell_type;
        typedef disk::quadrature<mesh_type, cell_type>          cell_quad_type;

        mesh_type   msh;
        loader_type loader(0,1,elems_1d);
        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return 1.0;
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return  0.0;
        };

        std::vector<dynamic_vector<RealType>>       Uh_Th( msh.cells_size());
        std::vector<tensor_quad_matrix_type>        sigma_Th( msh.cells_size());
        initialization<mesh_type,tensor_quad_matrix_type,cell_quad_type> (msh, sigma_Th, Uh_Th, degree);

        for(size_t iter = 1; iter < max_iters ; iter++)
            test_diffusion(msh, f, sf, degree, plst_data, Uh_Th, sigma_Th, iter);
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


        typedef disk::generic_mesh<RealType, 2>  mesh_type;
        typedef disk::fvca5_mesh_loader<RealType, 2>            loader_type;
        typedef typename mesh_type::scalar_type                 scalar_type;
        typedef Eigen::Matrix< RealType, 2, Eigen::Dynamic>     tensor_quad_matrix_type;
        typedef Eigen::Matrix< RealType, 2, 1 >                 tensor_type; //Esto cambia para 2D y 3D
        typedef typename mesh_type::cell                        cell_type;
        typedef disk::quadrature<mesh_type, cell_type>          cell_quad_type;


        mesh_type msh;
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

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return 0.;
        };

        std::vector<dynamic_vector<RealType>>    Uh_Th( msh.cells_size());
        std::vector<tensor_quad_matrix_type>     sigma_Th( msh.cells_size());
        initialization<mesh_type,tensor_quad_matrix_type,cell_quad_type> (msh, sigma_Th, Uh_Th, degree);

        for(size_t iter = 0; iter < max_iters ; iter++)
            test_diffusion(msh, f, sf, degree,
                plst_data,
                 Uh_Th,
                  sigma_Th,
                   iter);

        return 0;
    }
    return 0;
};
