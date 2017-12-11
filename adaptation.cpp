#include <iostream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <map>
#include "loaders/strtot.hpp"
#include "hho/hho.hpp"
#include "hho_pst.hpp"
#include "loaders/loader.hpp"
#include "loaders/mesh_adaptation_vpst.hpp"

#include <iostream>
#include <thread>
#include <vector>

//#define DEBUG

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

int main(int argc, char** argv)
{
    using RealType = double;
    size_t degree = 1;
    char  *filename   = nullptr;
    int   ch, num_remesh(1);
    RealType ratio(1.);
    disk::plasticity_data<RealType> pst;
    disk::mesh_parameters<RealType>  mp;

    while ( (ch = getopt(argc, argv, "r:t:")) != -1 )
    {
        switch(ch)
        {
            case 'r':
                ratio = atof(optarg);
                if (ratio < 0)
                {
                    std::cout << "Ratio must be positive. Falling back to 1." << std::endl;
                    ratio = 1.;
                }
                break;

            case 't':
                num_remesh = atof(optarg);
                if (num_remesh < 0)
                {
                    std::cout << "Max number of adaptations mu  st be positive. Falling back to 1." << std::endl;
                    num_remesh = 1;
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
        typedef typename mesh_type::cell                        cell_type;
        typedef disk::fvca5_mesh_loader<RealType, 2>            loader_type;
        typedef typename mesh_type::scalar_type                 scalar_type;

        mesh_type    msh;
        loader_type  loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);

        auto   dir_name   = "test_adaptation";
        std::vector<size_t> levels_vec(msh.cells_size());
        std::vector<RealType> vec(msh.cells_size());
        mp.mesh_name   = 1;

        for( size_t imsh = 0; imsh < num_remesh + 1; imsh++)
        {
            mp.directory   = "test_adaptation";

            dump_to_matlab(msh, mp.directory + "/new_msh_pr_"+ tostr(mp.mesh_name) + tostr(imsh) + ".m", levels_vec, imsh);

            disk::stress_based_mesh<mesh_type>  sbm(msh, levels_vec, pst, mp, imsh);

            bool do_refinement = sbm.borrar_marker(msh, imsh);
            std::cout << "*************BEGIN Refinement No. "<< imsh<<"****************" << std::endl;
            auto new_msh = sbm.template refine<loader_type>(msh, degree);
            std::vector<size_t>().swap(levels_vec);
            levels_vec = sbm.new_levels;

            dump_to_matlab(new_msh, mp.directory + "/new_msh_"+ tostr(mp.mesh_name) + tostr(imsh) + ".m", levels_vec, imsh);


            msh = new_msh;

            //dump_to_matlab(msh, mp.directory + "/new_msh_pr_"+ tostr(mp.mesh_name) + tostr(imsh) + ".m", levels_vec, imsh);

            //auto op = [&](cell_type& x) -> RealType { return measure(msh, x); };

            //parallel_transform(msh.cells_begin(), msh.cells_end(), vec.begin(), op);

        }

        return 0;
    }
    return 0;
};
