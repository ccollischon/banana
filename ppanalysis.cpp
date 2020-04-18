// point pattern analysis: compute Voronoi diagram for a set
// of seed points, and the IMTs of the Voronoi cells.
// optionally, with periodic boundary conditions
// 2019 Sebastian Kapfer <sebastian.kapfer@fau.de>
// 2019 Fabian Schaller <physik@fabian-schaller.de>
// FIXME this does not have tests yet
#include "papaya2.hpp"
#include "readarg.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = K::Point_2;
using DT = CGAL::Delaunay_triangulation_2<K>;
using AT = CGAL::Delaunay_triangulation_adaptation_traits_2<DT>;
using AP = CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT>;
using VD = CGAL::Voronoi_diagram_2<DT,AT,AP>;
using Ccb_halfedge_circulator = VD::Ccb_halfedge_circulator;

using namespace papaya2;

static
std::vector<Point_2> parse_plain_ascii_file (string infilename)
{
    std::vector<Point_2>ret;
    std::ifstream ifs(infilename);
    double x, y;
    while (ifs >> x >> y)
        ret.push_back(Point_2(x, y));
    return ret;
}

static
vec_t point_to_vec (const Point_2 &p)
{
    return { p.x(), p.y() };
}

struct VoronoiDiagram : VD
{
    double boxLx = 1., boxLy = 1.;
    bool periodic = false;
    std::vector<Point_2> seeds;

    // check if all the seeds are within the box.
    // this is required for the periodic version.
    // returns the first seed _not_ in the box in *offender if non-null.
    bool seeds_in_box(Point_2 *offender)
    {
        for (auto p : seeds)
        {
            if (! (0 <= p.x() && p.x() < boxLx && 0 <= p.y() && p.y() < boxLy))
            {
                if (offender)
                    *offender = p;
                return false;
            }
        }
        return true;
    }

    void make_voronoi()
    {
        this->clear();

        if (periodic)
        {
            for (auto p : seeds)
                this->insert(p);
            add_padding_seeds(1.);
        }
        else
        {
            for (auto p : seeds)
                this->insert(p);
        }

        if (!this->is_valid())
            throw std::runtime_error("CGAL thinks the Voronoi diagram is invalid.");
    }

    // FIXME make this adaptive
    void add_padding_seeds(double to)
    {
        double bbx_low = -boxLx * to;
        double bbx_high = (1+to) * boxLx;
        double bby_low = -boxLy * to;
        double bby_high = (1+to) * boxLy;
        for (int i = -1; i <= 1; ++i)
        for (int j = -1; j <= 1; ++j)
        {
            if (i==0 && j==0) continue;

            for (auto p : seeds)
            {
                Point_2 pp = { p.x() + i*boxLx, p.y() + j*boxLy };
                if (pp.x() < bbx_low) continue;
                if (pp.x() >= bbx_high) continue;
                if (pp.y() < bby_low) continue;
                if (pp.y() >= bby_high) continue;
                this->insert(pp);
            }
        }
    }
};

int main (int argc, const char **argv)
{
    string infilename, outfilename = "/dev/stdout";
    VoronoiDiagram vd;

    // process command-line arguments
    for (++argv; *argv; ++argv)
    {
        if (string(*argv) == "in")
        {
            infilename = read_arg<string>(argv++);
        }
        else if (string(*argv) == "out")
        {
            outfilename = read_arg<string>(argv++);
        }
        else if (string(*argv) == "boxL")
        {
            vd.boxLx = vd.boxLy = read_arg<double>(argv++);
            vd.periodic = true;
        }
        else if (string(*argv) == "-h")
        {
            std::cerr << "usage:\n (nonperiodic mode)\n\tppanalysis in <input filename> out <output filename>\n";
            std::cerr << " (periodic mode)\n";
            std::cerr << "\tppanalysis in <input filename> out <output filename> boxL <box side length>\n";
            return 0;
        }
        else
        {
            std::cerr << "illegal argument: " << *argv << "\n";
            return 1;
        }
    }

    if (infilename == "")
    {
        std::cerr << "missing argument: input filename\n";
        return 1;
    }

    std::ofstream outfile (outfilename);
    outfile << "# IMTs of Voronoi cells of a point pattern ";
    if (vd.periodic)
        outfile << "with periodic boundary conditions\n# box size: " << vd.boxLx
            << " " << vd.boxLy;
    outfile << "\n# seedx seedy area perimeter q2 q3 q4 q5 q6\n";
    outfile.precision(8);
    if (!outfile)
    {
        std::cerr << "unable to open output file: " << outfilename << "\n";
        return 1;
    }

    vd.seeds = parse_plain_ascii_file(infilename);

    if (vd.periodic)
    {
        Point_2 offender;
        if (!vd.seeds_in_box(&offender))
        {
            std::cerr << "There are points outside of the box: "
                << offender.x() << ", " << offender.y() << "\n";
            return 1;
        }
    }

    vd.make_voronoi();

    int label = 0;
    for (auto fit = vd.faces_begin(); fit != vd.faces_end(); ++fit)
    {
        const Point_2 seed = fit->dual()->point();

        if (seed != vd.seeds.at(label))
        {
            std::cerr << "CGAL decided to reorder our seed points.\n";
            return 1;
        }

        if (!fit->is_unbounded())
        {
            Ccb_halfedge_circulator heh = fit->ccb();
            const Ccb_halfedge_circulator start = heh;

            std::vector<vec_t> vertices;
            vertices.reserve(10);
            do {
                auto p = heh->source()->point();
                vertices.push_back(point_to_vec(p));
            } while (++heh != start);

            auto minkval = imt_polygon(vertices);

            outfile << seed.x() << " " << seed.y() << " " << minkval.area()
                << " " << minkval.perimeter();
            for (auto s : { 2, 3, 4, 5, 6 })
            {
                outfile << " " << minkval.msm(s);
            }
            outfile << "\n";
        }
        else
        {
            outfile << seed.x() << " " << seed.y() << " unbounded\n";
        }

        ++label;
        if (label == (int)vd.seeds.size())
            break;
    }

    return 0;
}
