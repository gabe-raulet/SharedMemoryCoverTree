#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include "misc.h"
#include "timer.h"

using namespace std;

/*
 * !!! IMPORTANT !!!
 *
 * this program assumes that vectors are stored with 32-bit floating
 * point values.
 */

static constexpr size_t realsize = sizeof(float);

int main(int argc, char *argv[])
{
    LocalTimer<double> timer;
    timer.start_timer();

    char *fname;
    size_t fsize, rsize, n;
    int dim;
    ifstream is;

    if (argc != 2)
    {
        fprintf(stderr, "usage: %s <ptsfname>\n", argv[0]);
        return 1;
    }

    fname = argv[1];

    is.open(fname, ios::binary | ios::in);
    is.seekg(0, is.end);
    fsize = is.tellg();
    is.seekg(0, is.beg);
    is.read((char*)&dim, sizeof(int));
    is.close();

    rsize = sizeof(int) + (dim * realsize);
    assert(fsize % rsize == 0);
    n = fsize / rsize;

    timer.stop_timer();

    fprintf(stderr, "[time=%.3f,msg::%s] read points file '%s' [num_points=%lu,dim=%d,size=%s]\n", timer.get_elapsed(), __func__, fname, n, dim, PrettyFileSize::str(fname).c_str());

    main_msg(argc, argv, timer.get_elapsed());
    return 0;
}

