#include <iostream>
#include <fstream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <map>
#include <list>
#include "loaders/loader.hpp"
#include <thread>
#include <vector>
#include "loaders/strtot.hpp"
#include "hho/hho.hpp"
#include "../viscoplasticity/hho_pst.hpp"

#include "loaders/mesh_adaptation_vpst.hpp"
#include <unsupported/Eigen/SparseExtra>

//#include "post_processing_all.hpp"
template<typename Parameters, typename MeshParameters>
bool
read_parameters(const std::string& meshfilename,
                const std::string& root,
                Parameters      & pst,
                MeshParameters  & mp)
{
    std::string directory = root;

    std::string datafname = root + "/parameters.txt";
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
    if(!(temp_reset_off == 1 || temp_reset_off ==0))
        std::cout << "Reset solutions: Yes (1) or No (0). Falling back to 0." << std::endl;
    if(!(temp_start_exact == 1 || temp_start_exact ==0))
        std::cout << "To start from exact solution choose: Yes (1) or No (0). Falling back to 0." << std::endl;

    mp.num_remesh = num_adapts * refine_on;
    mp.start_exact = (temp_start_exact == 1);


    if(mp.num_remesh < mp.initial_imsh)
        throw std::logic_error("'Number of Adaptations' must be >= than the 'Adaptive step'");

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


int main( int argc, char** argv)
{
    using RealType = double;
    size_t degree = 0;
    char  *filename   = nullptr;

    RealType Bingham(0), Load4marker(0);
    disk::plasticity_data<RealType> pst;
    disk::mesh_parameters<RealType>  mp;
    pst.method = false;

    argc -= optind;
    argv += optind;

    filename = argv[0];

    if (filename == nullptr)
    {
        std::cout << "no filename specified" << std::endl;
        return 1;
    }

    std:: string root = "test_finecoarse";

    if(!read_parameters(filename, root, pst, mp))
    {
        std::cout << "Problem loading parameters." << std::endl;
        return 1;
    }

    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<RealType, 2>         mesh_type;
        typedef disk::fvca5_mesh_loader<RealType, 2>    loader_type;

        mesh_type       msh;
        loader_type     loader;

        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        dump_to_matlab(msh, "Binning_old_mesh.m");

        #if 0
        binding_meshes<mesh_type> bm(msh);
        auto cells_marks = bm.marker(msh);
        bm.find_geometry(msh, cells_marks);
        auto new_msh = bm.bin(msh);
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AFTER BINNING!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        check_older_msh(new_msh);
        dump_to_matlab(new_msh, "Binning_new_mesh.m");
        #endif

        auto levels_vec = std::vector<size_t>(msh.cells_size(),0);

        disk::stress_based_mesh<mesh_type> sbm(msh, levels_vec, pst, mp, 1);

        auto any = sbm.temporal_marker_yield(msh);
        mesh_type temp_msh = sbm.template refine<loader_type>(msh, degree);

        disk::dump_to_matlab(temp_msh, "test_finecoarse/Refinement_mesh.m");

        std::vector<size_t>().swap(levels_vec);
        levels_vec    = sbm.new_levels;
        std::vector<size_t> ancestors_vec = sbm.ancestors;
        #if 0
        auto cells_marks_temp  = sbm.cells_marks;
        #endif

        //#if 0
        disk::binding_meshes<mesh_type> bm(temp_msh, levels_vec, ancestors_vec);
        bm.find_geometry(temp_msh, sbm.new_binning_marks);
        mesh_type new_msh = bm.bin(temp_msh);

        disk::dump_to_matlab(new_msh, "test_finecoarse/Aglomaration_mesh.m");

        //cells_marks = bm.cells_marks;
        check_older_msh(new_msh);
        //#endif
    }
};
