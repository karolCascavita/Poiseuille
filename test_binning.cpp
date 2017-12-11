#include <iostream>
#include <fstream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <map>
#include <list>
#include "loaders/loader.hpp"
//#include "post_processing_all.hpp"
#include <thread>
#include <unsupported/Eigen/SparseExtra>

template<typename T>
void
print(const std::vector<T> vec, const std::string name)
{
    std::cout << name << std::endl;

    for(auto& v: vec)
        std::cout << " "<< v;
    std::cout << std::endl;
}


struct greater
{
    template<typename T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};

template<typename T>
void
sort_uniq_greater(std::vector<T>& v)
{
    std::sort(v.begin(), v.end(), greater());
    auto uniq_iter = std::unique(v.begin(), v.end());
    v.erase(uniq_iter, v.end());
}

template<typename T>
void
sort_uniq(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    auto uniq_iter = std::unique(v.begin(), v.end());
    v.erase(uniq_iter, v.end());
}

struct bin_geometry
{
    std::vector<size_t> nodes2erase;
    std::vector<size_t> faces2erase;
    std::vector<size_t> cells2erase;
    std::vector<size_t> faces2bin;
    std::vector<size_t> cells2bin;

    bin_geometry(){};

    void print_all()
    {
        print(nodes2erase,"/ * nodes2erase */" );
        print(faces2erase,"/ * faces2erase */" );
        print(cells2erase,"/ * cells2erase */" );
        print(faces2bin,"/ * faces2bin */" );
        print(cells2bin,"/ * cells2bin */" );
    }
};

template<typename PointType>
auto
to_angle(const PointType& p, const PointType& o)
{
  return atan2(p.y() - o.y(), p.x() - o.x());
}

template<typename T, typename Storage>
void
sort_by_polar_angle(const disk::mesh<T,2,Storage>& msh,
        std::vector<typename disk::mesh<T,2,Storage>::point_type::id_type>& vec)
{

    typedef disk::mesh<T,2,Storage>     mesh_type;
    typedef std::pair< size_t, typename mesh_type::point_type>       pair_type;
    typedef typename mesh_type::point_type               point_type;
    typedef typename point_type::id_type                point_id_type;

    sort_uniq(vec);
    auto storage = msh.backend_storage();
    std::vector<point_type>  pts(vec.size());

    size_t i = 0;
    for(auto pid : vec)
        pts.at(i++) = storage->points.at(size_t(pid));

    auto ptsnum = pts.size();

    //Compute barycenter
    point_type  bar = point_type({0.0, 0.0});
    T           area = 0.0;

    for (size_t i = 0; i < ptsnum; i++)
    {
        auto p0 = pts[i];
        auto p1 = pts[(i+1)%ptsnum];

        auto a = p0.x()*p1.y() - p1.x()*p0.y();

        bar = bar + (p0 + p1) * a;
        area += a;
    }

    area *= 0.5;

    auto h = bar/(6.0 * area);

    //sorting point index by polar angle
    std::sort(vec.begin(),vec.end(),
        [&](const point_id_type & va, const point_id_type & vb )
        {
            auto pta = storage->points.at(va);
            auto ptb = storage->points.at(vb);

            auto theta_a = to_angle(pta, h);
            auto theta_b = to_angle(ptb, h);

            return (theta_a < theta_b);
        }
    );
    return;
}

template<typename T, typename Storage>
void
sort_cclockwise(const disk::mesh<T,2,Storage>& msh,
        std::vector<size_t>& cells2erase,
        std::vector<typename disk::mesh<T,2,Storage>::face::id_type>& fids_vec,
        std::vector<typename disk::mesh<T,2,Storage>::point_type::id_type>& pids_vec )
{
    typedef disk::mesh<T,2,Storage>     mesh_type;
    typedef std::pair< size_t, typename mesh_type::point_type>       pair_type;
    typedef typename mesh_type::point_type              point_type;
    typedef typename point_type::id_type                point_id_type;
    typedef typename mesh_type::node_type::id_type      node_id_type;
    typedef typename mesh_type::face::id_type           face_id_type;

    sort_uniq(pids_vec);
    sort_uniq(fids_vec);

    auto init_fid = fids_vec.at(0);
    //auto new_fids_vec = init_fid;
    auto init_fc  = *std::next(msh.faces_begin(), init_fid);
    auto new_pids_vec = init_fc.point_ids();

    auto left_id  = new_pids_vec[0];
    auto right_id = new_pids_vec[1];



    for(size_t i = 1; i < pids_vec.size()-1; i++)
    {
        print(new_pids_vec, "* new_pids_vec : ");

        for(auto& fid : fids_vec )
        {
            auto fc  = *std::next(msh.faces_begin(), fid);
            auto fc_pids  = fc.point_ids();

            if( left_id != fc_pids[0] && right_id == fc_pids[1])
            {
                left_id  = right_id;
                right_id = fc_pids[0];
                new_pids_vec.push_back(right_id);
                //new_fids_vec.push_back(fid);
                break;
            }
            if( left_id != fc_pids[1] && right_id == fc_pids[0])
            {
                left_id  = right_id;
                right_id = fc_pids[1];
                new_pids_vec.push_back(right_id);
                //new_fids_vec.push_back(fid);
                break;
            }
        }
    }
    print(new_pids_vec, "* new_pids_vec : ");

    auto new_pids_size = new_pids_vec.size();
    assert(new_pids_size == pids_vec.size());

    auto storage = msh.backend_storage();
    std::vector<point_type>  pts(pids_vec.size());

    std::cout << "pts_size:"<< pts.size()  << std::endl;
    std::cout << "pids_size:"<< pids_vec.size()  << std::endl;
    std::cout << "new_pids_size:"<< new_pids_vec.size()  << std::endl;

    size_t i = 0;
    for(auto& pid : new_pids_vec)
    {
        auto pt = storage->points.at(size_t(pid));
        pts.at(i) = pt;
        std::cout << "pid ("<< i <<"): "<< pid<<"  . p.x() = " << pt.x();
        std::cout << " ; p.y() = "<< pt.y() << std::endl;
        i++;
    }

    #if 0
    bool cw(false), ccw(false);
    for(size_t i = 0; i < pts.size(); i++)
    {
        auto p0 = pts.at(i);
        auto p1 = pts.at((i + 1) % pts.size());
        auto p2 = pts.at((i + 2) % pts.size());

        std::cout << "p0.x() = " << p0.x()<< " ; p0.y() = "<< p0.y() << std::endl;
        std::cout << "p1.x() = " << p1.x()<< " ; p1.y() = "<< p1.y() << std::endl;
        std::cout << "p2.x() = " << p2.x()<< " ; p2.y() = "<< p2.y() << std::endl;


        auto v1 = ( p1 - p0 ).to_vector();
        auto v =  ( p2 - p1 );
        auto v2 = (point<T,2>({v.y(),- v.x()})).to_vector();

        std::cout << " product ("<< i<<") : "<<v1.dot(v2) << std::endl;

        if( v1.dot(v2) > T(0) || std::abs(v1.dot(v2)) < 10.e-10 )
        {
            ccw = true;
            std::cout << "** ccw" << std::endl;
        }
        if( v1.dot(v2) < T(0))
        {
            cw = true;
            std::cout << "** cw" << std::endl;
        }

    }

    if(cw && ccw)
        throw std::logic_error("Face are ordered in CW and CCW. Check all faces have same direction.");
    if(cw)
        std::reverse(new_pids_vec.begin(),new_pids_vec.end());
    #endif

    bool cw(false), ccw(false);

    for(size_t i = 0; i < new_pids_size; i++)
    {
        auto p0 = new_pids_vec.at(i);
        auto p1 = new_pids_vec.at((i + 1) % new_pids_size);

        auto is_found =  false;

        for(auto& cl_id: cells2erase)
        {
            auto cl  = *std::next(msh.cells_begin(), cl_id);
            auto cl_pids  = cl.point_ids();
            auto cl_pids_size = cl_pids.size();
            for(size_t j = 0; j < cl_pids_size; j++)
            {
                auto cp0 = cl_pids.at(j);
                auto cp1 = cl_pids.at((j + 1) % cl_pids_size);

                if( p0 == cp0 && p1 == cp1)
                {
                    ccw = true;
                    is_found = true;
                    break;
                }
                if( p0 == cp1 && p1 == cp0)
                {
                    cw = true;
                    is_found = true;
                    break;
                }
            }
            if(is_found)
                continue;
        }
        if(!is_found)
            throw  std::logic_error(" Points not found in any cell points");
    }
    if(cw && ccw)
        throw std::logic_error("Face are ordered in CW and CCW. Check all faces have same direction.");
    if(cw)
        std::reverse(new_pids_vec.begin(),new_pids_vec.end());

    //std::vector<point_id_type>

    //   pids_vec = new_pids_vec;
    size_t j(0);
    for(auto npid : new_pids_vec)
        pids_vec.at(j++) = point_id_type(npid);
    return;
}


template<typename MeshType>
class binding_meshes
{
    typedef MeshType mesh_type;
    typedef typename mesh_type::scalar_type  scalar_type;

    typedef typename mesh_type::cell         cell_type;
    typedef typename mesh_type::face         face_type;
    typedef typename mesh_type::node_type    node_type;
    typedef typename mesh_type::point_type    point_type;

    typedef typename face_type::id_type      face_id_type;
    typedef typename cell_type::id_type      cell_id_type;
    typedef typename node_type::id_type      node_id_type;
    typedef typename point_type::id_type     point_id_type;

    std::vector<size_t> check_cells;
    struct bin_geometry
    {
        std::vector<size_t> nodes2erase;
        std::vector<size_t> faces2erase;
        std::vector<size_t> cells2erase;
        std::vector<face_id_type> faces2bin;
        std::vector<cell_type>   cells2bin;

        bin_geometry(){};

        void print_all()
        {
            print(nodes2erase,"/ * nodes2erase */" );
            print(faces2erase,"/ * faces2erase */" );
            print(cells2erase,"/ * cells2erase */" );
            print(faces2bin,"/ * faces2bin */" );
        }
        friend bin_geometry operator+(const bin_geometry& bg1, const bin_geometry& bg2)
        {
            bin_geometry ret;
            ret = bg1;
            ret.nodes2erase.insert(ret.nodes2erase.begin(),
                                bg2.nodes2erase.begin(), bg2.nodes2erase.end());
            ret.faces2erase.insert(ret.faces2erase.begin(),
                                bg2.faces2erase.begin(), bg2.faces2erase.end());
            ret.cells2erase.insert(ret.cells2erase.begin(),
                                bg2.cells2erase.begin(), bg2.cells2erase.end());
            ret.cells2bin.insert(ret.cells2bin.begin(),
                                bg2.cells2bin.begin(), bg2.cells2bin.end());

            sort_uniq_greater(ret.nodes2erase);
            sort_uniq_greater(ret.faces2erase);
            sort_uniq_greater(ret.cells2erase);
            return ret;
        }

    };

    bin_geometry bgeo_all;

    auto
    connectivity(const mesh_type& msh)
    {
       std::vector<std::list<size_t>> node_to_face;
       node_to_face.resize( msh.points_size());
       for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
       {
           auto fc = *itor;
           auto pids = fc.point_ids();
           auto fid  = msh.lookup(fc);
           node_to_face.at(pids[0]).push_back(fid);
           node_to_face.at(pids[1]).push_back(fid);
       }
       for(auto& l: node_to_face)
       {
           l.sort();
           l.unique();
       }
       std::cout<<"mesh_node_to_face"<<std::endl;

       size_t i = 0;
       for(auto& n2f:  node_to_face )
       {
           std::cout << i++ << " ";
           for(auto& l : n2f)
               std::cout << l <<" ";
           std::cout << std::endl;
       }

       return node_to_face;
    }
    std::vector<size_t>
    mark_by_node_neighbors(const mesh_type& msh,
                            std::vector<size_t>& initial_marks)
    {
        std::vector<size_t> ret(msh.cells_size());
        std::vector<std::list<size_t>> node_to_face(msh.points_size());

        for(auto& cell : msh)
        {
            auto cl_id = msh.lookup(cell);
            auto cell_mark = initial_marks.at(cl_id);
            if(cell_mark == 1)
            {
                auto fcs =  faces(msh, cell);
                for(auto& fc: fcs)
                {
                    auto pids = fc.point_ids();
                    auto fid  = msh.lookup(fc);
                    node_to_face.at(pids[0]).push_back(fid);
                    node_to_face.at(pids[1]).push_back(fid);
                }
            }
        }

        std::cout<<"node_to_face"<<std::endl;

        size_t i = 0;
        for(auto& n2f:  node_to_face )
        {
            n2f.sort();
            n2f.unique();
            std::cout << i++ << " ";
            for(auto& l : n2f)
                std::cout << l <<" ";
            std::cout << std::endl;
        }

         std::vector<std::list<size_t>>  mesh_node_to_face = connectivity(msh);

        for(auto& cell : msh)
        {
            auto cl_id = msh.lookup(cell);
            auto cell_mark = initial_marks.at(cl_id);
            if(cell_mark == 1)
            {
                bool all_ngh_true = true;
                auto pids = cell.point_ids();
                for(auto pid : pids)
                {
                    auto l1 =  node_to_face.at(pid);
                    auto l2 =  mesh_node_to_face.at(pid);
                    if(! std::equal (l1.begin (), l1.end (), l2.begin ()))  //if( !(l1 == l2) );
                    {
                        std::cout <<" loc_cell: pid("<< pid << ") ";
                        for(auto& l : node_to_face.at(pid))
                            std::cout << l <<" ";
                        std::cout << std::endl;
                        std::cout <<" msh_cell: pid("<< pid << ") ";
                        for(auto& l : mesh_node_to_face.at(pid))
                            std::cout << l <<" ";
                        std::cout << std::endl;

                        all_ngh_true = false;
                    }
                }
                if(all_ngh_true)
                    ret.at(cl_id) = 1;
            }
        }
        return ret;
    }
    auto
    mark_cell(  const mesh_type& msh,
                const cell_type& cell)
    {
        auto id = msh.lookup(cell);
        //check_cells.at(id) = 1;
        #if 0
        if(id < 8 || (id <= 44 && id >= 42) || id ==37)
            return true;
        #endif
        if( (id >= 2 &&  id <= 6) || (id >= 51 && id <= 53) || id == 123|| id == 121 || id == 124 || id == 111)
            return true;

        return false;
    };


    void
    search4faces(const mesh_type& msh, const cell_type& cell, bin_geometry & bgeo,
                const std::vector<size_t>& cells_marks)
    {
        auto fcs = faces(msh, cell);
        auto cid =  msh.lookup(cell);

        check_cells.at(cid) = 1;
        bgeo.cells2erase.push_back(cid);

        std::cout << "cell "<< cid << std::endl;
        std::cout << "* faces : ";
        for(auto& face: fcs)
            std::cout <<  "  " << msh.lookup(face);
        std::cout<< std::endl;

        for(auto& face: fcs)
        {
            auto fid = msh.lookup(face);

            if(!msh.is_boundary(face))
            {
                auto ngh = msh.neighbor(cell, face);
                auto nid = msh.lookup(ngh);
                auto is_checked = check_cells.at(nid) == 1? true: false;
                if(cells_marks.at(nid))
                {
                    check_cells.at(nid) = 1;
                    std::cout << "faces2erase.push_back ("<< fid <<")" << std::endl;

                    bgeo.faces2erase.push_back(fid);

                    if(!is_checked)
                    {
                        std::cout << "keep searching faces with neighbor ("<< nid <<")" << std::endl;

                        search4faces(msh, ngh, bgeo, cells_marks);
                    }
                }
                else
                {
                    std::cout << "faces2bin.push_back ("<< fid <<")" << std::endl;

                    bgeo.faces2bin.push_back(fid);
                }
            }
            else
            {
                std::cout << "boundary: faces2bin.push_back ("<< fid <<")" << std::endl;

                bgeo.faces2bin.push_back(fid);
            }
        }
        return;
    }

    auto
    search4nodes(const mesh_type&  msh, bin_geometry& bgeo)
    {
        std::vector<size_t> nodes2erase;

        for(auto& fid : bgeo.faces2erase)
        {
            auto fc   = *std::next(msh.faces_begin(), fid);
            auto pids = fc.point_ids();

            for(auto& p : pids)
            {
                auto is_found = false;
                for(auto& poly_fid : bgeo.faces2bin)
                {
                    auto poly_fc  = *std::next(msh.faces_begin(), poly_fid);
                    auto poly_pids = poly_fc.point_ids();

                    if(p == poly_pids[0] || p == poly_pids[1])
                    {
                        is_found = true;
                        break;
                    }
                }
                if(!is_found)
                    nodes2erase.push_back(p);
            }
        }
        sort_uniq_greater(nodes2erase);
        bgeo.nodes2erase = nodes2erase;
        return;
    }
public:

    binding_meshes(const mesh_type& msh)
    {
        check_cells = std::vector<size_t>(msh.cells_size());
    };
    auto
    marker(const mesh_type& msh)
    {

        auto DIM = MeshType::dimension;
        std::vector<size_t> initial_marks(msh.cells_size());
        size_t num_cells_marked = 0;
        for(auto & cell : msh)
        {
            if(mark_cell(msh, cell))
            {
                auto cell_id = msh.lookup(cell);
                initial_marks.at(cell_id) = 1;
            }
        }

        //auto ret = mark_by_node_neighbors(msh, initial_marks);
        //return ret;

        return initial_marks;
    }

    void
    find_geometry(const mesh_type& msh,
                         const std::vector<size_t>& cells_marks)
    {
        std::vector<bin_geometry> bgeo_temp;

        for(auto& cell : msh)
        {
            auto cell_id = msh.lookup(cell);
            bool is_checked = check_cells.at(cell_id) == 1? true : false;

            if(!is_checked && cells_marks.at(cell_id))
            {
                check_cells.at(cell_id) = 1;
                bin_geometry bgeo;
                search4faces(msh, cell, bgeo, cells_marks);

                sort_uniq_greater(bgeo.faces2erase);
                sort_uniq_greater(bgeo.faces2bin);
                sort_uniq_greater(bgeo.cells2erase);

                search4nodes(msh, bgeo);

                auto f2b = bgeo.faces2bin;

                sort_uniq(f2b);
                bgeo.cells2bin.push_back(cell_type(f2b));

                bgeo.print_all();

                bgeo_temp.push_back(bgeo);
            }
        }

        for(auto&  bg_geo : bgeo_temp)
            bgeo_all = bgeo_all + bg_geo;
        return;
    }

    void
    erase_geometry(const mesh_type& msh, mesh_type& new_msh)
    {

        auto storage = msh.backend_storage();
        auto new_storage = new_msh.backend_storage();

        //Erase nodes and points
        for(auto nid: bgeo_all.nodes2erase)
        {
            new_storage->nodes.erase(new_storage ->nodes.begin() + nid);
            new_storage->points.erase(new_storage ->points.begin() + nid);
        }
        //Erase faces
        for(auto fid: bgeo_all.faces2erase)
        {
            new_storage->edges.erase(new_storage->edges.begin() + fid);
            new_storage->boundary_info.erase(new_storage->boundary_info.begin() + fid);
            new_storage->boundary_edges.erase(new_storage->boundary_edges.begin() + fid);
        }
        //Erase cells
        for(auto cid: bgeo_all.cells2erase)
        {
            new_storage->surfaces.erase(new_storage->surfaces.begin() + cid);
            new_storage->special_surfaces.erase(new_storage->special_surfaces.begin() + cid);
        }

        return;
    }

    void
    copy(const mesh_type& msh, mesh_type& new_msh)
    {
        auto storage = msh.backend_storage();
        auto new_storage = new_msh.backend_storage();

        new_storage->nodes    = storage->nodes;
        new_storage->points   = storage->points;
        new_storage->edges    = storage->edges;
        new_storage->surfaces = storage->surfaces;
        new_storage->boundary_info  = storage->boundary_info;
        new_storage->boundary_edges = storage->boundary_edges;
        new_storage->special_surfaces = storage->special_surfaces;
    }

    std::pair<bool, std::vector<point_id_type>>
    bin_is_special_polygon(const mesh_type& old_msh, const cell_type& cl)
    {
        std::cout << "INSIDE bin_is_special_polygon" << std::endl;

        auto storage = old_msh.backend_storage();
        auto pids  = cl.point_ids();
        auto fids  = cl.faces_ids();

        auto num_pts = pids.size();
        std::vector<static_vector<scalar_type,2>> ns(fids.size());
        size_t i = 0;
        for(size_t i =0; i < pids.size(); i++)
        {
            auto p1 = storage->points.at(pids.at(i));
            auto p2 = storage->points.at(pids.at(i%pids.size()));

            auto v = p2 - p1;
            ns.at(i++) = (point<scalar_type,2>({-v.y(), v.x()})).to_vector();
        }

        size_t hanging_count = 0;
        for (size_t i = 0; i < num_pts; i++)
        {
            auto nl = ns.at(i);
            auto nr = ns.at(( i + 1 )% num_pts);

            if( std::abs(nl.dot(nr) - scalar_type(1)) < 1.e-5)
                hanging_count++;
        }
        if(hanging_count > pids.size())
            throw std::logic_error(" More than one hanging node by face. Check correct counting of hanging nodes.");
        size_t num_vertices = pids.size() - hanging_count;
        if(num_vertices <= 0)
        {
            std::cout << "num_vertices : "<< num_vertices << std::endl;
            std::cout << "num_points   : "<< pids.size() << std::endl;
            std::cout << "hanging_nodes: "<< hanging_count << std::endl;

            throw std::logic_error(" Number of vertices are <= 0. Check correct counting of hanging nodes.");
        }

         //Identify vertices
         std::vector<point_id_type> vertices( num_vertices);

         size_t vcount = 0;
         for(size_t i = 0; i < num_pts; i++)
         {
             auto nl = ns.at(i);
             auto nr = ns.at((i + 1) %num_pts);

            if( std::abs( nr.dot(nl)- scalar_type(1)) >= 1.e-5)
                vertices.at(vcount++) =  point_id_type((i+1)%num_pts);
        }

        if(vcount != num_vertices)
        {
            std::cout << "vcount "<< vcount << std::endl;
            throw  std::logic_error(" Incorrect procedure to find vertices");
        }

        bool has_hang_nodes(false);
        if(vertices.size() != pids.size())
            has_hang_nodes = true;

        return std::make_pair(has_hang_nodes, vertices);
    }

    void
    add_new_cells(const mesh_type& msh, mesh_type& new_msh)
    {
        auto new_storage = new_msh.backend_storage();

        for(auto& cl: bgeo_all.cells2bin)
        {
            auto fids = cl.faces_ids();

            std::vector<typename face_type::id_type> cell_faces(fids.size());
            std::vector<point_id_type> cell_pids(fids.size() * 2);

            size_t i(0), ii(0);
            for(auto fid : fids)
            {
                auto f = *next(msh.faces_begin(), fid);
                auto pids = f.point_ids();

                for( auto fp : f.point_ids())
                    cell_pids.at(ii++) = point_id_type(fp);
            }

            //sort_by_polar_angle(msh, cell_pids);
            sort_cclockwise(msh, bgeo_all.cells2erase, fids, cell_pids);

            print(cell_pids, "cell_pids_after_sort: ");

            std::vector<face_id_type> new_fids(fids.size());

            size_t count(0);
            for(size_t i = 0; i < fids.size(); i++)
            {
                auto lp = cell_pids.at(i);
                auto rp = cell_pids.at((i+1)%cell_pids.size());

                for(auto fid : fids)
                {
                    std::cout << "face : "<<  fid<< std::endl;
                    auto fc = *next(msh.faces_begin(), fid);
                    auto pids = fc.point_ids();

                    auto left_it  = std::find(pids.begin(), pids.end(), lp);
                    auto right_it = std::find(pids.begin(), pids.end(), rp);

                    if (left_it != pids.end() &&  right_it != pids.end())
                    {
                        new_fids.at(count++) = fid;
                        break;
                    }

                }
            }
            assert(count == fids.size());

            auto cell = cell_type(new_fids);
            cell.set_point_ids(cell_pids.begin(), cell_pids.end());

            print(cell.point_ids(), "cell.pids: ");

            new_storage->surfaces.push_back(cell);

            auto p = bin_is_special_polygon(msh, cell);
            new_storage->special_surfaces.push_back(p);

        }
        std::vector<int> index(new_storage->surfaces.size(), 0);
        for(int i = 0 ; i != index.size() ; i++)
            index.at(i) = i;


        std::sort(index.begin(), index.end(),[&](const int& a, const int& b)
            {
                cell_type s1  = new_storage->surfaces.at(a);
                cell_type s2  = new_storage->surfaces.at(b);
                return (s1 < s2);
            }
        );

        std::sort(new_storage->surfaces.begin(), new_storage->surfaces.end());

        typedef std::vector<std::pair<bool, std::vector<point_id_type>>> ss_type;

        auto new_special_surfaces = ss_type(new_storage->surfaces.size());

        for(size_t i = 0; i < new_storage->surfaces.size(); i++)
        {
            auto idx = index.at(i);
            new_special_surfaces.at(i) = new_storage->special_surfaces.at(idx);
        }

        new_storage->special_surfaces = new_special_surfaces;

        return;
    }

    mesh_type
    bin(const mesh_type& msh)
    {
        //Copy
        mesh_type new_msh;
        auto storage = msh.backend_storage();
        auto new_storage = new_msh.backend_storage();

        copy(msh, new_msh);
        erase_geometry(msh, new_msh);
        add_new_cells(msh, new_msh);

        auto renumbering_cells = [&](const cell_type & cl) -> auto
        {
            auto fids = cl.faces_ids();
            auto pids = cl.point_ids();
            std::vector<typename face_type::id_type> cell_faces(fids.size());

            size_t i(0);
            for(auto fid : fids)
            {
                auto f = *next(msh.faces_begin(), fid);
                auto pids = f.point_ids();
                auto new_fid = find_element_id(new_storage->edges.begin(),
                                                new_storage->edges.end(), f);
                //std::cout << " * face "<< msh.lookup(f)<<":"<< pids[0]<< "  " << pids[1] << std::endl;

                if(!new_fid.first)
                    throw std::invalid_argument("Element not found");
                else
                    cell_faces.at(i++) = new_fid.second;
            }
            //print(cell_faces, "* new_faces :");
            //print(pids, "* old_pids :");
            auto cell = cell_type(cell_faces);

            cell.set_point_ids(pids.begin(), pids.end()); //just copy all pids

            return cell;
        };

        std::transform(new_msh.cells_begin(), new_msh.cells_end(),
                            new_msh.cells_begin(), renumbering_cells);

        size_t count_id(0);
        auto change_id = [&](const cell_type & cl) -> auto
        {
            typedef typename cell_type::id_type cell_id_type;
            cell_type ret = cl;
            ret.set_element_id(cell_id_type(count_id++));
            return ret;
        };

        std::transform(new_msh.cells_begin(), new_msh.cells_end(),
                                new_msh.cells_begin(), change_id);

        //change numeration
        auto renumbering_faces = [&](const face_type & fc) -> auto
        {
            size_t i(0),j(0);
            auto nids = fc.point_ids();
            std::vector<typename node_type::id_type> face_nodes(nids.size());
            std::vector<typename point_type::id_type> face_points_ids(nids.size());

            for(auto & nid : nids)
            {
                auto new_nid = find_element_id( new_storage ->nodes.begin(),
                                                new_storage ->nodes.end(),
                                                                node_type(nid));
                if(!new_nid.first)
                    throw std::invalid_argument("Element not found");
                else
                {
                    face_nodes.at(i++) = node_id_type(new_nid.second);
                    face_points_ids.at(j++) = point_id_type(new_nid.second);
                }
            }

            auto face = face_type(face_nodes);
            face.set_point_ids(face_points_ids.begin(), face_points_ids.end());
            return face;
        };

        std::transform(new_msh.faces_begin(), new_msh.faces_end(),
                        new_msh.faces_begin(), renumbering_faces );

        std::cout << "!!! _______ before cell nodes transformation _______!!!" << std::endl;

        //for(auto cl: new_msh)
        for(auto & cl : new_storage->surfaces)
        {
            std::cout << "cell : "<< cl.get_id()<< std::endl;//  new_msh.lookup(cl) << std::endl;

            auto fids = cl.faces_ids();
            auto pids = cl.point_ids();

            //print(fids, " * faces :");
            //print(pids, " * nodes :");
        }
        auto renumbering_cell_nodes = [&](const cell_type & cl) -> auto
        {
            std::cout << "cell : "<< cl.get_id() << " vs "<< new_msh.lookup(cl)<< std::endl;

            auto pids = cl.point_ids();

            print(pids, "  * pids : ");

            std::vector<point_id_type> new_pids(pids.size());

            size_t i(0);
            for(auto pid : pids)
            {
                std::cout << " * pid : "<< pid;

                auto itor = std::lower_bound(new_storage->nodes.begin(),
                                        new_storage->nodes.end(), node_type(pid));

                if (itor != new_storage->nodes.end() && !(node_type(pid) < *itor))
                {
                    auto idx = point_id_type(std::distance(new_storage->nodes.begin(), itor));
                    std::cout << " -> "<< idx << std::endl;
                    new_pids.at(i++)  = idx;
                }
                else
                    throw std::invalid_argument("This is a bug: node not found");
            }
            assert(i == pids.size());
            cell_type ret = cl;
            ret.set_point_ids(new_pids.begin(), new_pids.end());

            return ret;
        };

        std::transform(new_msh.cells_begin(), new_msh.cells_end(),
                            new_msh.cells_begin(), renumbering_cell_nodes);

        /* Nodes */
        auto new_nodes_size = new_storage->nodes.size();
        std::vector<node_type> nodes(new_nodes_size);
        for (size_t i = 0; i < new_nodes_size; i++)
            nodes[i] = node_type(point_identifier<2>(i));

        new_storage->nodes = std::move(nodes);

        std::cout << "!!! _______ Final Cells _______!!!" << std::endl;

        //for(auto cl: new_msh)
        for(auto & cl : new_storage->surfaces)
        {
            std::cout << "cell : "<< cl.get_id()<< std::endl;//  new_msh.lookup(cl) << std::endl;

            auto fids = cl.faces_ids();
            auto pids = cl.point_ids();

            //print(fids, " * faces :");
            //print(pids, " * nodes :");

        }

        size_t iii(0);
        std::cout << "boundary_edges: " << std::endl;

        for(bool is_boundary : new_storage->boundary_edges)
            std::cout << "face ("<<iii++ <<"): "<< is_boundary << std::endl;

        check_older_msh(new_msh);

        return new_msh;
    }

};



int main (int argc, char** argv )
{

    using RealType = double;
    char    *filename   = nullptr;

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

        dump_to_matlab(msh, "test_binning/Binning_old_mesh.m");
 
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BEFORE BINNING!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;  
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BEFORE BINNING!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;  
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BEFORE BINNING!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;  

        check_older_msh(msh);
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BEFORE BINNING!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;  
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BEFORE BINNING!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;  
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BEFORE BINNING!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;  

        binding_meshes<mesh_type> bm(msh);
        auto cells_marks = bm.marker(msh);
        bm.find_geometry(msh, cells_marks);
        auto new_msh = bm.bin(msh);
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AFTER BINNING!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        check_older_msh(new_msh);
        dump_to_matlab(new_msh, "test_binning/Binning_new_mesh.m");
    }
};
