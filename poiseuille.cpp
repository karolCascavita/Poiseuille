#include <iostream>
#include <regex>
#include <unistd.h>

#include <map>
#include "loaders/loader.hpp"
#include "hho/hho.hpp"
#include "visualization/post_processing_all.hpp"

#include <unsupported/Eigen/SparseExtra>

template<typename MeshType, typename Function, typename Solution, typename PlasticData,
                                            typename StorageTensorQuadMatrixType>
bool test_diffusion (MeshType& msh,       /* handle to the mesh */
           const Function&     load,            /* rhs */
           const Solution&     solution,        /* solution of the problem */
           const size_t        degree,             /* degree of the method */
           const PlasticData&  plst_data,   /* Material and ALG parameters */
           std::vector<dynamic_vector<typename MeshType::scalar_type>> & Uh_Th,    /* storage of velocities*/
           StorageTensorQuadMatrixType   & xi_norm_Th,    /* storage of velocities*/
           StorageTensorQuadMatrixType   & sigma_Th,
           StorageTensorQuadMatrixType   & gamma_Th,
           typename MeshType::scalar_type& stop_error,
           int iter,
           std::ofstream& its)   /* storage of tensor sigma  in each cell and in each quad point*/
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

    if(iter == 0)
    {
        factor = plst_data.mu;
        plst_switch = 0.;
    }
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
        auto local_sigma_storage  =  sigma_Th[cell_id];
        auto local_gamma_storage  =  gamma_Th[cell_id];
        auto local_xi_norm        =  xi_norm_Th[cell_id];

        plasticity.compute(msh, cl, gradrec_nopre.oper, Uh_Th[cell_id],local_xi_norm, local_sigma_storage, local_gamma_storage);
        //WK: if diffusion changes with domain, it would be mandatory to include mu here, to do it locally

        auto scnp = statcond_nopre.compute_Bingham(msh, cl, loc, cell_rhs, plst_switch*plasticity.rhs);
        //auto scnp = statcond_nopre.compute(msh, cl, loc, cell_rhs);

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

    scalar_type gamma_error = 0.0;


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
        plasticity.compute(msh, cl, gradrec_nopre.oper, Uh_Th[cell_id], xi_norm_Th[cell_id], sigma_Th[cell_id], gamma_Th[cell_id]);
        dynamic_vector<scalar_type>   x = statcond_nopre.recover_Bingham(msh, cl, loc, cell_rhs, xFs, plst_switch*plasticity.rhs);
        //dynamic_vector<scalar_type>   x = statcond_nopre.recover(msh, cl, loc, cell_rhs, xFs);

        Uh_Th[cell_id]  = x;

        dynamic_vector<scalar_type> rec(cb.size());
        rec.tail(cb.size()-1) = gradrec_nopre.oper * x;
        rec(0) = x(0);

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
        size_t col = 0;
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


            auto cell_id    =   msh.lookup(cl);
            auto gamma      =   gamma_Th[cell_id].col(col);
            auto dphi       =   cb.eval_gradients(msh, cl, qp.point());
            auto col_range  =   cb.range(1,degree+1);
            auto row_range  =   disk::dof_range(0,mesh_type::dimension);

            matrix_type dphi_matrix =   disk::make_gradient_matrix(dphi);
            matrix_type dphi_taken  =   take(dphi_matrix, row_range, col_range);
            matrix_type dphi_rec_uh =   dphi_taken * gradrec_nopre.oper*x;

            gamma_error  += ((gamma - dphi_rec_uh).cwiseProduct(gamma - dphi_rec_uh)).sum()   ;
            col++;
        }
            //err_dof_k += compute_L2_error(dld, degree, solution, x);
            //err_dof_kp += compute_L2_error(dld, degree+1, solution, x);

    }

    ofs.close();
    #if 0
    std::cout << "Mesh diameter: " << diam << std::endl;
    std::cout << "L2-norm error, dof, K:   " << std::sqrt(err_dof_k) << std::endl;
    std::cout << "L2-norm error, fun, K:   " << std::sqrt(err_fun_k) << std::endl;
    std::cout << "L2-norm error, dof, K+1: " << std::sqrt(err_dof_kp) << std::endl;
    std::cout << "L2-norm error, fun, K+1: " << std::sqrt(err_fun_kp) << std::endl;
    #endif
    scalar_type TOL_solver = 1.e-13;
    scalar_type TOL = 1.e-4;

    gamma_error =  std::sqrt(gamma_error);
    scalar_type stop_criterio = std::abs(stop_error - gamma_error);
    std::cout << "l2-norm error, gamma - Du:" << gamma_error << std::endl;
    std::cout << "l2-norm error, iterations:" << stop_criterio << std::endl;

    its<< iter << " "<<gamma_error <<" "<<stop_criterio<<std::endl;
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

    return 0;
};



template<typename T>
struct plastic_data
{
    plastic_data(): Lref(1.), Vref(1.), Bn(0.1), mu(1.), alpha(1.), f(1.),method(true)
    {}
    T f;                    //WK: Cuidado porque f deberia ser el valor externo de la fuente.
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
                std::vector<dynamic_vector<MeshType::scalar_type>> & Uh_Th,
                std::vector<TensorQuadType> & xi_norm_Th,
                std::vector<TensorQuadType> & sigma_Th,
                std::vector<TensorQuadType> & gamma_Th,
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

        Uh_Th[i]        =  dynamic_vector<typename mesh_type::scalar_type>::Zero(dsr.total_size());
        xi_norm_Th[i]   =  tensor_quad_matrix_type::Zero(1,cell_quadpoints.size());
        sigma_Th[i]     =  tensor_quad_matrix_type::Zero(MeshType::dimension, cell_quadpoints.size());
        gamma_Th[i]     =  tensor_quad_matrix_type::Zero(MeshType::dimension, cell_quadpoints.size());

        ++i;
    }
};



template<typename T, typename MeshType, typename TensorByQuadsType>
class stress_based_mesh
{
public:

    typedef MeshType                                mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::node_type           node_type;
    typedef typename mesh_type::edge_type           edge_type;
    typedef typename mesh_type::cell_iterator       cell_itor;

    T m_yield;
    stress_based_mesh(const T yield): m_yield(yield)
    {}

    //bool sort_by_x_coordinate( const point_type& p1, const point_type& p2)  { return (p1.x() < p2.x()); }

    struct stress_point
    {
        point_type p;
        bool   b;
    };
    typedef std::vector<stress_point>  stress_point_vector_type;

    bool check_if_to_refine(const TensorByQuadsType& mat)
    {
        std::cout << "maxCoeff   " << std::abs(mat.maxCoeff())<< std::endl;
        std::cout << "minCoeff   " << std::abs(mat.minCoeff())<< std::endl;

        if(std::abs(mat.maxCoeff()) >=  m_yield)
        {
            if (std::abs(mat.minCoeff())< m_yield)
                {
                    std::cout << " refinement true  " << std::endl;
                    return true;
                }
            return false;
        }
        return false;
    }

    template< typename QuadPointType, size_t DIM>
    stress_point_vector_type
    check_stress(const TensorByQuadsType& mat, const std::vector<QuadPointType>& cq_pts)
    {
        stress_point_vector_type spm(cq_pts.size());

        size_t j = 0;
        auto num_qps    = cq_pts.size();

        for(auto itor = cq_pts.begin(); itor != cq_pts.end(); itor++, j++)
        {
                spm.at(j).p = cq_pts[j].point();

                std::cout << "xi_norm = "<< mat(0,j)  << std::endl;
                if(std::abs(mat(0,j)) > m_yield)
                    spm.at(j).b = true;
                else
                    spm.at(j).b = false;

        }
        return spm;
    }

    //template<typename DIM>
    //void refine_cell {};

    //template<>
    void refine_cell(const size_t n, const stress_point_vector_type& spts, point_type& p_left, point_type& p_right, std::vector<point_type>& npts)
    {
        if ( n == 0 )
        {
            if (spts.at(1).b != spts.at(0).b)
            {
                auto new_pt = (spts.at(1).p + spts.at(0).p)/2.;
                npts.push_back(new_pt);
                return;
            }
        }

        if ( n < 0 )
        {
            std::cout << "/* WK: No refinement. Since there are less than two values to compare. */" << std::endl;
            return;
        }
        if ( n > 0 )
        {
            if (spts.at(n+1).b != spts.at(n).b)
            {
                auto new_pt = (spts.at(n+1).p + spts.at(n).p)/2.;
                npts.push_back(new_pt);
                p_right = new_pt;
            }
            refine_cell(n-1, spts, p_left, p_right, npts);
            return;
        }
    }

    template <typename CellQuadType,size_t DIM>
    void re_populate_mesh(const mesh_type & msh, mesh_type & re_msh,const std::vector<TensorByQuadsType>& xi_norm_Th, const int degree)
    {

        auto storage = msh.backend_storage();
        auto storage_rm = re_msh.backend_storage();

        std::vector<size_t> erase_vector;
        typedef typename CellQuadType::quadpoint_type quadpoint_type;
        auto num_cells = msh.cells_size();
        size_t nds_counter = 0;

        //for (cell_itor itor = msh.cells_begin(); itor != msh.cells_begin() + num_cells; ic++)
        //for (size_t ic = 0 ; ic < num_cells; ic++)
        for(auto& cl:msh)
        {
            std::cout << "/*  BEGIN */" << std::endl;

            size_t cell_id = msh.lookup(cl);
            std::cout << "/*  1  reussi */" << std::endl;

            auto xi_norms_storage = xi_norm_Th[cell_id];
            auto cq         = CellQuadType(2*degree+2);
            std::vector< quadpoint_type> cq_points  = cq.integrate(msh, cl);
            auto num_qps    = cq_points.size();

            std::cout << "================== cell "<< cell_id<<"======================" << std::endl;
            auto on_remesh = check_if_to_refine(xi_norms_storage);

            auto cell_pts  = points(msh, cl);
            storage_rm->points.push_back(cell_pts.at(0));

            if( on_remesh )
            {
                stress_point_vector_type spts = check_stress<quadpoint_type,DIM>(xi_norms_storage, cq_points);

                std::vector<point_type> new_pts;

                refine_cell(num_qps - 2, spts,  cell_pts.at(0), cell_pts.at(1), new_pts);
                std::reverse(new_pts.begin(),new_pts.end());
                auto num_npts = new_pts.size();

                std::cout << "BEGIN" << std::endl;
                //std::cout << "num_old_pts"<<storage->points.size() << std::endl;
                std::cout << "old_pts ("<< cell_id<<") =  ["<< cell_pts[0].x() << ","<< cell_pts[1].x()<<"]"<< std::endl;

                /*WK: This is only "usefull" 1D, maybe is better  to remove it*/
                //std::sort(new_pts.begin(),new_pts.end(), sort_by_x_coordinate);

                std::cout << "new_pts = ";
                for(size_t i= 0; i < new_pts.size(); i++)
                    std::cout << new_pts[i].x()<<" ";
                std::cout<<std::endl;

                /* Points */
                for(size_t  i  =  0 ;  i < num_npts ;  i++ )
                    storage_rm->points.push_back(std::move(new_pts.at(i)));

                std::cout << "num_pts after push_back = "<<storage_rm->points.size() << std::endl;


                /* Cells or edges*/
                for(size_t i = 0; i < num_npts + 1; i ++)
                {
                    std::cout << "====== i = "<<i <<"; nds_counter = "<<nds_counter<<"======"<< std::endl;

                    auto n0 = typename node_type::id_type(nds_counter);
                    auto n1 = typename node_type::id_type(nds_counter + 1);
                    auto e = edge_type{{n0, n1}};

                    std::vector<point_identifier<1>> pts(2);
                    pts[0] = point_identifier<1>(nds_counter);
                    pts[1] = point_identifier<1>(nds_counter+1);
                    e.set_point_ids(pts.begin(), pts.end());
                    storage_rm->edges.push_back(e);

                    std::cout << "p ("<< i <<") = ["<< pts[0]<<","<< pts[1]<<"]" << std::endl;
                    std::cout << "====== END i ======"<< std::endl;
                    nds_counter++;
                }
                std::cout << "/* END */" << std::endl;
            }
            else
            {
                auto n0 = typename node_type::id_type(nds_counter);
                auto n1 = typename node_type::id_type(nds_counter + 1);
                auto e = edge_type{{n0, n1}};

                std::vector<point_identifier<1>> pts(2);
                pts[0] = point_identifier<1>(nds_counter);
                pts[1] = point_identifier<1>(nds_counter+1);
                e.set_point_ids(pts.begin(), pts.end());
                storage_rm->edges.push_back(e);
                nds_counter++;
            }

            if( cell_id == msh.cells_size() -1 )
             {
                 storage_rm->points.push_back(cell_pts.at(1));
                 std::cout << "last cell = "<< cell_id << std::endl;
                 std::cout << "p  = ["<< cell_pts.at(1).x() <<"]" << std::endl;

             }
        } //for cells


        /* nodes */
        auto num_rm_pts = storage_rm->points.size();

        std::cout << "num_rm_pts"<< num_rm_pts << std::endl;

        for (size_t i = 0; i < num_rm_pts; i++)
            storage_rm->nodes.push_back(node_type(point_identifier<1>(i)));

        storage_rm->boundary_nodes.resize(num_rm_pts);
        storage_rm->boundary_nodes.at(0) = true;
        storage_rm->boundary_nodes.at(num_rm_pts - 1) = true;

        for(auto& cl:re_msh)
        {
            auto test_p = points(re_msh,cl);
            size_t cell_id = re_msh.lookup(cl);
            std::cout << "points("<< cell_id<<") =  ["<< test_p[0].x() << ","<< test_p[1].x()<<"]"<< std::endl;

        }

    }

};




int main (int argc, char** argv )
{

    using RealType = double;

    char    *filename   = nullptr;
    int     degree      = 1;
    int     elems_1d    = 10;
    int     ch, method_name,refine_on;
    size_t  max_iters   = 1;
    plastic_data<RealType> plst_data;
    bool do_remesh = false;

    while ( (ch = getopt(argc, argv, "k:n:h:u:l:v:a:B:i:m:")) != -1 )
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
                if (plst_data.mu <= 0)
                {
                    std::cout << "Viscosity must be bigger than 0. Falling back to 1." << std::endl;
                    plst_data.mu = 1.0;
                }
                break;
            case 'a':
                plst_data.alpha = atof(optarg);
                if (plst_data.alpha <= 0)
                {
                    std::cout << "Augmentation parameter must be bigger than 0. Falling back to 1." << std::endl;
                    plst_data.alpha = 1.;
                }
                break;
            case 'B':
                plst_data.Bn = atof(optarg);
                if (plst_data.Bn < 0)
                {
                    std::cout << "Bingham number must be positive. Falling back to 0. Solving diffusion case." << std::endl;
                    plst_data.Bn = 0.0;
                    max_iters    = 1;
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


        std::ofstream its("Errors.dat");
        std::vector<dynamic_vector<RealType>>    Uh_Th( msh.cells_size());
        std::vector<tensor_quad_matrix_type>     sigma_Th( msh.cells_size()), gamma_Th( msh.cells_size());
        std::vector<tensor_quad_matrix_type>     xi_norm_Th( msh.cells_size());

        initialization<mesh_type,tensor_quad_matrix_type,cell_quad_type> (msh, Uh_Th, xi_norm_Th, sigma_Th, gamma_Th, degree);

        disk::post_processing<scalar_type, 1, mesh_type> pp(msh,degree);
        scalar_type stop_error(0.);
        char numstr[4];
        std::string result;
        std::string ext = "1D_Bi";

        size_t mm  = (int)(5*plst_data.Bn);
        sprintf(numstr, "%i",mm);
        std::string dir_name = ext + numstr;
std::cout << "dir_name" << dir_name<< std::endl;

        for(int iter = 0; iter < max_iters ; iter++)
            {
                std::cout << "iter = "<<iter << std::endl;
                bool brk = test_diffusion(msh, f, sf, degree, plst_data, Uh_Th, xi_norm_Th,sigma_Th,gamma_Th, stop_error,iter,its);

                /*All this should be inside test_diffusion*/
                pp.vtk_writer(msh,dir_name + "/solution_Th",degree,Uh_Th,iter);
                if (plst_data.method == true)
                    pp.lambda_vtk_writer(msh,dir_name + "/lambda_Th", degree, sigma_Th, iter);
                else
                    pp.lambda_vtk_writer(msh,dir_name +"/sigma_Th", degree, sigma_Th, iter);
                    pp.lambda_vtk_writer(msh,dir_name + "/gamma_Th", degree, gamma_Th, iter);   /* storage of velocities*/
                    pp.lambda_vtk_writer(msh,dir_name + "/xi_norm_Th", degree, xi_norm_Th, iter);   /* storage of velocities*/
                if(brk ==true)
                    {
                        its<< iter<<" "<< stop_error << "0. " <<std::endl;
                        its<<'\n'<<std::endl;
                        break;
                    }
            }

        if(do_remesh)
        {
            auto yield = 0.5*plst_data.Bn * plst_data.f * plst_data.Lref;
            std::cout << "yield"<< yield << std::endl;
            stress_based_mesh<scalar_type, mesh_type, tensor_quad_matrix_type> sbm_1(yield);
            mesh_type re_msh_1;
            sbm_1.re_populate_mesh<cell_quad_type,1>(msh, re_msh_1, xi_norm_Th, degree);

            std::vector<dynamic_vector<RealType>>    Uh_1Th(re_msh_1.cells_size());
            std::vector<tensor_quad_matrix_type>     sigma_1Th(re_msh_1.cells_size()), gamma_1Th(re_msh_1.cells_size());
            std::vector<tensor_quad_matrix_type>     xi_norm_1Th( re_msh_1.cells_size());

            initialization<mesh_type,tensor_quad_matrix_type,cell_quad_type> (re_msh_1, Uh_1Th, xi_norm_1Th, sigma_1Th, gamma_1Th, degree);
            disk::post_processing<scalar_type, 1, mesh_type> pp_1(re_msh_1,degree);
            stop_error = 0.;
            for(int iter = 0; iter < max_iters ; iter++)
            {
                std::cout << "iter = "<<iter << std::endl;
                bool brk = test_diffusion(re_msh_1, f, sf, degree, plst_data, Uh_1Th, xi_norm_1Th,sigma_1Th,gamma_1Th, stop_error,iter,its);
                /*All this should be inside test_diffusion*/
                pp_1.vtk_writer(re_msh_1,dir_name +"/solution_1Th",degree,Uh_1Th,iter);
                if (plst_data.method == true)
                    pp_1.lambda_vtk_writer(re_msh_1,dir_name + "/lambda_1Th", degree, sigma_1Th, iter);
                else
                    pp_1.lambda_vtk_writer(re_msh_1,dir_name + "/sigma_1Th", degree, sigma_1Th, iter);
                    pp_1.lambda_vtk_writer(re_msh_1,dir_name + "/gamma_1Th", degree, gamma_1Th, iter);   /* storage of velocities*/
                    pp_1.lambda_vtk_writer(re_msh_1,dir_name + "/xi_norm_1Th", degree, xi_norm_1Th, iter);   /* storage of velocities*/
                    if(brk ==true)
                        {
                            its<< iter<<" "<< stop_error << "0. " <<std::endl;
                            its<<'\n'<<std::endl;
                            break;
                        }


            }


            stress_based_mesh<scalar_type, mesh_type, tensor_quad_matrix_type> sbm_2(yield);
            mesh_type re_msh_2;
            sbm_2.re_populate_mesh<cell_quad_type,1>(re_msh_1, re_msh_2, xi_norm_1Th, degree);

            std::vector<dynamic_vector<RealType>>    Uh_2Th(re_msh_2.cells_size());
            std::vector<tensor_quad_matrix_type>     sigma_2Th(re_msh_2.cells_size()), gamma_2Th(re_msh_2.cells_size());
            std::vector<tensor_quad_matrix_type>     xi_norm_2Th( re_msh_2.cells_size());

            initialization<mesh_type,tensor_quad_matrix_type,cell_quad_type> (re_msh_2, Uh_2Th, xi_norm_2Th, sigma_2Th, gamma_2Th, degree);
            disk::post_processing<scalar_type, 1, mesh_type> pp_2(re_msh_2,degree);
            stop_error = 0.;
            for(int iter = 0; iter < max_iters ; iter++)
            {
                std::cout << "iter = "<<iter << std::endl;
                bool brk =test_diffusion(re_msh_2, f, sf, degree, plst_data, Uh_2Th, xi_norm_2Th,sigma_2Th,gamma_2Th, stop_error,iter,its);
                /*All this should be inside test_diffusion*/
                pp_2.vtk_writer(re_msh_2,dir_name +"/solution_2Th",degree,Uh_2Th,iter);
                if (plst_data.method == true)
                    pp_2.lambda_vtk_writer(re_msh_2, dir_name +"/lambda_2Th", degree, sigma_2Th, iter);
                else
                    pp_2.lambda_vtk_writer(re_msh_2, dir_name + "/sigma_2Th", degree, sigma_2Th, iter);
                    pp_2.lambda_vtk_writer(re_msh_2, dir_name + "/gamma_2Th", degree, gamma_2Th, iter);   /* storage of velocities*/
                    pp_2.lambda_vtk_writer(re_msh_2, dir_name +"/xi_norm_2Th", degree, xi_norm_2Th, iter);   /* storage of velocities*/
                    if(brk ==true)
                        {
                            its<< iter<<" "<< stop_error << "0. " <<std::endl;
                            its<<'\n'<<std::endl;
                            break;
                        }

            }

            stress_based_mesh<scalar_type, mesh_type, tensor_quad_matrix_type> sbm_3(yield);
            mesh_type re_msh_3;
            sbm_3.re_populate_mesh<cell_quad_type,1>(re_msh_2, re_msh_3, xi_norm_2Th, degree);

            std::vector<dynamic_vector<RealType>>    Uh_3Th(re_msh_3.cells_size());
            std::vector<tensor_quad_matrix_type>     sigma_3Th(re_msh_3.cells_size()), gamma_3Th(re_msh_3.cells_size());
            std::vector<tensor_quad_matrix_type>     xi_norm_3Th( re_msh_3.cells_size());

            initialization<mesh_type,tensor_quad_matrix_type,cell_quad_type> (re_msh_3, Uh_3Th, xi_norm_3Th, sigma_3Th, gamma_3Th, degree);
            disk::post_processing<scalar_type, 1, mesh_type> pp_3(re_msh_3,degree);
            stop_error = 0.;
            for(int iter = 0; iter < max_iters ; iter++)
            {
                std::cout << "iter = "<<iter << std::endl;
                bool brk =test_diffusion(re_msh_3, f, sf, degree, plst_data, Uh_3Th, xi_norm_3Th,sigma_3Th,gamma_3Th, stop_error,iter,its);
                /*All this should be inside test_diffusion*/
                pp_3.vtk_writer(re_msh_3,dir_name +"/solution_3Th",degree,Uh_3Th,iter);
                if (plst_data.method == true)
                    pp_3.lambda_vtk_writer(re_msh_3, dir_name +"/lambda_3Th", degree, sigma_3Th, iter);
                else
                    pp_3.lambda_vtk_writer(re_msh_3, dir_name +"/sigma_3Th", degree, sigma_3Th, iter);
                    pp_3.lambda_vtk_writer(re_msh_3, dir_name + "/gamma_3Th", degree, gamma_3Th, iter);   /* storage of velocities*/
                    pp_3.lambda_vtk_writer(re_msh_3, dir_name + "/xi_norm_3Th", degree, xi_norm_3Th, iter);   /* storage of velocities*/
                    if(brk ==true)
                        {
                            its<< iter<<" "<< stop_error << "0. " <<std::endl;
                            its<<'\n'<<std::endl;
                            break;
                        }

            }
        }
        its.close();

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
        std::vector<tensor_quad_matrix_type>     sigma_Th( msh.cells_size()), gamma_Th( msh.cells_size());
        std::vector<tensor_quad_matrix_type>     xi_norm_Th( msh.cells_size());

        initialization<mesh_type,tensor_quad_matrix_type,cell_quad_type> (msh, Uh_Th, xi_norm_Th, sigma_Th, gamma_Th, degree);

        disk::post_processing<scalar_type, 2, mesh_type> pp(msh,degree);
        scalar_type stop_error(0.);
        char numstr[4];
        std::string result;
        std::string dir = "2D_Bi";


        int mm  = 10*plst_data.Bn;
        sprintf(numstr, "%i",mm);
        std::string dir_name   = dir + numstr;
        std::string Error_name = "/Errors.dat";
        //std::string error_data =
        std::ofstream its(dir_name + Error_name);

        for(int iter = 0; iter < max_iters ; iter++)
        {
            std::cout << "iter = "<<iter << std::endl;
            bool brk = test_diffusion(msh, f, sf, degree, plst_data, Uh_Th, xi_norm_Th, sigma_Th, gamma_Th, stop_error, iter,its);

            pp.vtk_writer(msh,dir_name +"/solution_Th",degree,Uh_Th,iter);

            /*All this should be inside test_diffusion*/
            if (plst_data.method == true)
                pp.lambda_vtk_writer(msh, dir_name +"/lambda_Th", degree, sigma_Th, iter);   /* storage of velocities*/
            else
                pp.lambda_vtk_writer(msh,dir_name + "/sigma_Th", degree, sigma_Th, iter);   /* storage of velocities*/

            pp.lambda_vtk_writer(msh, dir_name +"/gamma", degree, gamma_Th, iter);   /* storage of velocities*/
            pp.lambda_vtk_writer(msh, dir_name +"/xi_norm", degree, xi_norm_Th, iter);   /* storage of velocities*/
            if(brk ==true)
                {
                    its<< iter<<" "<< stop_error << "0. " <<std::endl;
                    its<<'\n'<<std::endl;
                    break;
                }

        }
        its.close();
        return 0;
    }
    return 0;
};
