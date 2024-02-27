#include "Point.h"
#include "CoverTree.h"
#include "read_args.h"
#include <cassert>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <omp.h>

int main(int argc, char *argv[])
{
    int64_t n;
    int d = 2;
    int seed = -1;
    const char *output;

    if (find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -n INT    number of points [required]\n");
        fprintf(stderr, "         -o STR    output filename [required]\n");
        fprintf(stderr, "         -d INT    point dimension [default: %d]\n", d);
        fprintf(stderr, "         -s INT    seed [default: random]\n");
        fprintf(stderr, "         -h        help message\n");
        return 0;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    output = read_string_arg(argc, argv, "-o", NULL);
    d = read_int_arg(argc, argv, "-d", &d);
    seed = read_int_arg(argc, argv, "-s", &seed);

    std::vector<Point> points = Point::random_points(n, d, seed);
    Point::write_points(points, output);

    return 0;
}

