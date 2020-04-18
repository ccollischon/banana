#ifndef _banana_originalbanana_hpp_

#define _banana_originalbanana_hpp_

#include "papaya2.hpp"
#include "photofunctions.hpp"

template<typename PHOTO>
void originalBanana(string maskfilename, PHOTO infile, string outfilename, double min_thresh, double max_thresh, double num_thresh)
{
    if (maskfilename != "") 
    {
        std::vector<std::vector<double>> boxesToInclude = readTable(settings.objectDIR+maskfilename);
        includeBoxes(infile, boxesToInclude);
    }
    
    Datafile outfile (outfilename);
    outfile.comment ("threshold area perim euler q2 arg2 q3 arg3 q4 arg4 q5 arg5 q6 arg6 q7 arg7 q8 arg8");

    //std::ofstream contours("contours.dat");

    for (auto thresh : logspace (min_thresh, max_thresh, num_thresh, true))
    {
        auto imt = imt_interpolated_marching_squares (infile, thresh);
        outfile << thresh << imt.area() << imt.perimeter() << imt.euler() << imt.msm (2) << std::arg (imt.imt (2)) << imt.msm (3) << std::arg (imt.imt (3)) << imt.msm (4) << std::arg (imt.imt (4)) << imt.msm (5) << std::arg (imt.imt (5)) << imt.msm (6) << std::arg (imt.imt (6)) << imt.msm (7) << std::arg (imt.imt (7)) << imt.msm (8) << std::arg (imt.imt (8)) << std::endl;
        //GnuplottableContour gc(contours);
        //contours << "# contour " << thresh << '\n';
        //trace_isocontour_interpolated_marching_squares (&gc, infile, thresh);
        //contours << '\n';
        //contours << '\n';
    }
    

//if 0
    // dump a contour for demonstration
    // gnuplot it via
    //          plot "contour.dat" w vec
    //GnuplottableContour gc ("contour.dat");
    //trace_isocontour_interpolated_marching_squares (&gc, infile, 1e-3);
//endif
}



#endif
