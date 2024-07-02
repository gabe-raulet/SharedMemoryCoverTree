#include <iostream>
#include <random>
#include <vector>
#include <unistd.h>
#include "pindex.h"
#include "treeindex.h"
#include "mpi_env.h"
#include "ptraits.h"
#include "misc.h"

using namespace std;

using Index = int64_t;
using Comm = MPIEnv::Comm;

using PointTraits = SelectPoint<FPSIZE, PTDIM>::Traits;
using Real = PointTraits::Real;
using Point = PointTraits::Point;
using IndexSet = unordered_set<Index>;
using IndexVector = vector<Index>;
using IndexSetVector = vector<IndexSet>;

template <bool T> struct SelectIndex;
template <> struct SelectIndex<true> { using TIndex = TreeIndex<PointTraits, Index>; };
template <> struct SelectIndex<false> { using TIndex = PointIndex<PointTraits, Index>; };

using TIndex = SelectIndex<static_cast<bool>(TREE)>::TIndex;

void read_options(int argc, char *argv[], char *&fname, double& cutoff, int& iters, double& damping, char *&outprefix, Comm comm);
void read_points_file(vector<Point>& mypoints, const char *fname, Comm comm);
void build_point_index(TIndex& pidx, const vector<Point>& mypoints, double cutoff);
void build_rgraph(TIndex& pidx, double radius, IndexSetVector& rgraph, int iter);
void write_rgraph_file(const IndexSetVector& rgraph, const char *outprefix, int iter, Comm comm);

int main_mpi(int argc, char *argv[]);
int main(int argc, char *argv[])
{
    MPIEnv::initialize(&argc, &argv);
    int err = main_mpi(argc, argv);
    MPIEnv::finalize();
    return err;
}

int main_mpi(int argc, char *argv[])
{
    Comm comm = Comm::comm_world();
    int myrank = comm.rank();
    int nprocs = comm.size();

    auto timer = comm.get_timer();
    timer.start_timer();

    char *fname, *outprefix = NULL;
    double cutoff;
    int iters = 5;
    double damping = 0.5;

    read_options(argc, argv, fname, cutoff, iters, damping, outprefix, comm);

    vector<Point> mypoints;
    vector<IndexSetVector> rgraphs;
    TIndex pidx(comm);

    read_points_file(mypoints, fname, comm);
    build_point_index(pidx, mypoints, cutoff);

    double radius = cutoff;
    for (int i = 0; i < iters; ++i)
    {
        rgraphs.emplace_back();
        IndexSetVector& rgraph = rgraphs.back();
        build_rgraph(pidx, radius, rgraph, i+1);
        radius *= damping;
    }

    if (outprefix)
        for (int iter = 0; const auto& rgraph : rgraphs)
            write_rgraph_file(rgraph, outprefix, ++iter, comm);

    timer.stop_timer();
    if (!myrank) main_msg(argc, argv, timer.get_max_time(), timer.get_sum_time(), nprocs);

    return 0;
}

void read_options(int argc, char *argv[], char *&fname, double& cutoff, int& iters, double& damping, char *&outprefix, Comm comm)
{
    int myrank = comm.rank();
    int nprocs = comm.size();

    auto timer = comm.get_timer();
    timer.start_timer();

    auto usage = [&argv, &cutoff, &iters, &damping] (int err, bool isroot)
    {
        if (isroot)
        {
            fprintf(stderr, "Usage: %s [options] <ptsfname> <cutoff>\n", argv[0]);
            fprintf(stderr, "Options: -n INT   number of iterations [%d]\n", iters);
            fprintf(stderr, "         -d FLOAT damping factor [%.2f]\n", damping);
            fprintf(stderr, "         -o       output graph prefix\n");
            fprintf(stderr, "         -h       help message\n");
        }

        MPIEnv::exit(err);
    };

    int c;
    while ((c = getopt(argc, argv, "n:d:o:h")) >= 0)
    {
        if      (c == 'n') iters = atoi(optarg);
        else if (c == 'd') damping = atof(optarg);
        else if (c == 'o') outprefix = optarg;
        else if (c == 'h') usage(0, !myrank);
    }

    if (argc - optind < 2)
    {
        if (!myrank) fprintf(stderr, "[err::%s] missing argument(s)\n", __func__);
        usage(1, !myrank);
    }

    fname = argv[optind++];
    cutoff = atof(argv[optind]);

    timer.stop_timer();

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] :: [fname='%s',cutoff=%.3f,iters=%d,damping=%.2f,outprefix='%s',tree=%d]\n", timer.get_max_time(), timer.get_avg_time(), __func__, fname, cutoff, iters, damping, outprefix, !!TREE);
    }
}

void read_points_file(vector<Point>& mypoints, const char *fname, Comm comm)
{
    int myrank = comm.rank();
    int nprocs = comm.size();

    auto timer = comm.get_timer();

    timer.start_timer();
    PointTraits::read_from_file(back_inserter(mypoints), fname, comm);
    timer.stop_timer();

    size_t totsize;
    size_t mysize = mypoints.size();

    comm.reduce(mysize, totsize, 0, MPI_SUM);

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] read %lu points from file '%s' [size=%s]\n", timer.get_max_time(), timer.get_avg_time(), __func__, totsize, fname, PrettyFileSize::str(fname).c_str());
    }
}

void build_point_index(TIndex& pidx, const vector<Point>& mypoints, double cutoff)
{
    auto comm = pidx.getcomm();
    int myrank = comm.rank();
    int nprocs = comm.size();

    auto timer = comm.get_timer();

    timer.start_timer();
    pidx.build(mypoints, cutoff);
    timer.stop_timer();

    Index totsize;
    Index mysize = mypoints.size();
    comm.reduce(mysize, totsize, 0, MPI_SUM);

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] brute force indexed %lld points with cutoff %.3f\n", timer.get_max_time(), timer.get_avg_time(), __func__, totsize, cutoff);
    }
}

void build_rgraph(TIndex& pidx, double radius, IndexSetVector& rgraph, int iter)
{
    auto comm = pidx.getcomm();
    int myrank = comm.rank();
    int nprocs = comm.size();

    auto timer = comm.get_timer();

    timer.start_timer();
    Index my_n_edges = pidx.build_rgraph(radius, rgraph);
    timer.stop_timer();

    Index n_edges;
    comm.reduce(my_n_edges, n_edges, 0, MPI_SUM);

    Index totsize;
    Index mysize = rgraph.size();
    comm.reduce(mysize, totsize, 0, MPI_SUM);

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] %d. built %.2f-graph [n_edges=%lld,avg_degree=%.2f]\n", timer.get_max_time(), timer.get_avg_time(), __func__, iter, radius, n_edges, (n_edges+0.0)/totsize);
    }
}

void write_rgraph_file(const IndexSetVector& rgraph, const char *outprefix, int iter, Comm comm)
{
    int myrank = comm.rank();
    int nprocs = comm.size();

    auto timer = comm.get_timer();
    timer.start_timer();

    Index myoffset;
    Index mysize = rgraph.size();
    comm.exscan(mysize, myoffset, MPI_SUM, static_cast<Index>(0));

    stringstream ss;
    IndexVector neighs_v;
    Index my_n_edges = 0, n_edges;
    Index my_n_verts = rgraph.size(), n_verts;

    ss << outprefix << "_" << iter << ".rgraph";
    string fname = ss.str();

    ss.clear();
    ss.str({});

    for (Index u = 0; const auto& neighs : rgraph)
    {
        my_n_edges += neighs.size();
        neighs_v.assign(neighs.cbegin(), neighs.cend());

        for (Index v : neighs_v)
            ss << u+myoffset+1 << " " << v+1 << "\n";

        u++;
    }

    comm.reduce(my_n_edges, n_edges, 0, MPI_SUM);
    comm.reduce(my_n_verts, n_verts, 0, MPI_SUM);

    string sbuf;

    if (!myrank)
    {
        stringstream ss2;
        ss2 << n_verts << " " << n_edges << "\n" << ss.str();
        sbuf = ss2.str();
    }
    else sbuf = ss.str();

    comm.file_write_at_all_once(fname.c_str(), sbuf);

    timer.stop_timer();

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] wrote file '%s' [size=%s]\n", timer.get_max_time(), timer.get_avg_time(), __func__, fname.c_str(), PrettyFileSize::str(fname.c_str()).c_str());
    }
}
