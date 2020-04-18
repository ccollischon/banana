#pragma once
#include <iostream>
#include <stdlib.h>

// Papaya2
// header-only library for computing 2D irreducible Minkowski tensors
// 2018-2019 Sebastian Kapfer <sebastian.kapfer@fau.de>
// 2019 Fabian Schaller <physik@fabian-schaller.de>

#include "src/tools.hpp"

namespace papaya2
{

// BasicPhoto
//
// This class is an example for the (static) interface for all
// PHOTO types that are intended to be used as image sources.
// should provide the proper member functions, though they
// need not necessarily derive from Photo.
//
// required functions for a PHOTO are:
// operator() for pixel access.
// a typedef "data_t" which gives the type of the value in each
// pixel.
// width() and height() return the width and height, in pixels.
// operator()(-1,y) etc. return a black (intensity 0) padding pixel.
// smaller indices may fail, throw or cause undefined behavior.
// the boundary pixels are not intended to be written to.
// pixel_width and pixel_height are used to convert to physical
// units.
template <typename TYPE>
struct BasicPhoto
{
    using data_t = TYPE;

    void set_coordinates (double x0, double y0, double x1, double y1,
        int width_ /* number of pixels */, int height_ /* dito */)
    {
        origin_x = x0;
        origin_y = y0;
        numx = width_;
        numy = height_;
        lenx = (x1-x0) / width_;
        leny = (y1-y0) / height_;
        data.resize (size_t (numx+2) * size_t (numy+2));
    }

    // helper method to fill the Photo with a discretized
    // (sampled) version of the provided function
    // use set_coordinates first to define the coordinates.
    // the method below uses integration instead of sampling.
    template <typename FUNCTION>
    void sample_function (const FUNCTION &function)
    {
        for (int j = 0; j < numy; ++j)
        {
            for (int i = 0; i < numx; ++i)
            {
                double x = (i+.5) * lenx + origin_x;
                double y = (j+.5) * leny + origin_y;
                at (i, j) = function (x, y);
            }
        }
    }
    
    // trapezoidal rule to numerically integrate the function over
    // the pixel area.
    template <typename FUNCTION>
    void integrate_function (const FUNCTION &function)
    {
        for (int j = 0; j < numy; ++j)
        for (int i = 0; i < numx; ++i)
        {
            const int ss = 5;
            int total_w = 0;
            for (int subj = 0; subj <= ss; ++subj)
            for (int subi = 0; subi <= ss; ++subi)
            {
                double x = (i+subi*1./ss) * lenx + origin_x;
                double y = (j+subj*1./ss) * leny + origin_y;
                int w = (1 + (subi != 0 && subi != ss)) * (1 + (subj != 0 && subj != ss));
                total_w += w;
                at (i, j) += w * function (x, y);
            }
            at (i, j) /= total_w;
        }
    }
    
    BasicPhoto<data_t> divideThresh(const BasicPhoto<data_t>& b, double thresh = 0)
    {
        
        int w = std::min(numx, b.width());
        int h = std::min(numy, b.height());
        for(int i=0; i<w; i++)
        for(int j=0; j<h; j++)
        {
            if (b.at(i,j) > thresh) at(i,j) = at(i,j)/b.at(i,j);
            else at(i,j) = 0;
        }
        return *this;
    }
    
    data_t const max_value()
    {
        return *(max_element(data.begin(), data.end()));
    }
    
    
    
    // read access to the image
    data_t const operator () (int i, int j) const
    {
        //std::cerr << "Jo\n";
        return at (i, j);
    }
    
    data_t &operator () (int i, int j)
    {
        //std::cerr << "Jo\n";
        return at (i, j);
    }
    
    int width () const
    {
        return numx;
    }

    int height () const
    {
        return numy;
    }

    double pixel_width () const
    {
        return lenx;
    }

    double pixel_height () const
    {
        return leny;
    }
    
    vec_t origin() const
    {
        return { origin_x, origin_y };
    }

    vec_t upper_right() const
    {
        return { origin_x + lenx*numx, origin_y + leny*numy };
    }
    
    
    BasicPhoto<data_t> &operator+=(const BasicPhoto<data_t>& b)
    {
        
        int w = std::min(numx, b.width());
        int h = std::min(numy, b.height());
        for(int i=0; i<w; i++)
        for(int j=0; j<h; j++)
        {
            at(i,j) = at(i,j)+b.at(i,j);
        }
        return *this;
    }
    
    BasicPhoto<data_t> &operator+=(const data_t b)
    {
        
        int w = numx;
        int h = numy;
        for(int i=0; i<w; i++)
        for(int j=0; j<h; j++)
        {
            at(i,j) = at(i,j)+b;
        }
        return *this;
    }
    
    BasicPhoto<data_t> &operator-=(const data_t b)
    {
        
        int w = numx;
        int h = numy;
        for(int i=0; i<w; i++)
        for(int j=0; j<h; j++)
        {
            at(i,j) = at(i,j)-b;
        }
        return *this;
    }
    
    BasicPhoto<data_t> &operator-=(const BasicPhoto<data_t>& b)
    {
        
        int w = std::min(numx, b.width());
        int h = std::min(numy, b.height());
        for(int i=0; i<w; i++)
        for(int j=0; j<h; j++)
        {
            at(i,j) = at(i,j)-b.at(i,j);
        }
        return *this;
    }
    
    
    BasicPhoto<data_t> &operator/=(double b)
    {
        for(int i=0; i<numx; i++)
        for(int j=0; j<numy; j++)
        {
            at(i,j) = at(i,j)/b;
        }
        return *this;
    }
    
    BasicPhoto<data_t> &operator/=(const BasicPhoto<data_t>& b)
    {
        
        int w = std::min(numx, b.width());
        int h = std::min(numy, b.height());
        for(int i=0; i<w; i++)
        for(int j=0; j<h; j++)
        {
            if (b.at(i,j) != 0) at(i,j) = at(i,j)/b.at(i,j);
            else at(i,j) = 0;
        }
        return *this;
    }
    
    BasicPhoto<data_t> &operator*=(const BasicPhoto<data_t>& b)
    {
        
        int w = std::min(numx, b.width());
        int h = std::min(numy, b.height());
        for(int i=0; i<w; i++)
        for(int j=0; j<h; j++)
        {
            at(i,j) = at(i,j)*b.at(i,j);
        }
        return *this;
    }
    
    BasicPhoto<data_t> &operator*=(const double b)
    {
        for(int i=0; i<numx; i++)
        for(int j=0; j<numy; j++)
        {
            at(i,j) = at(i,j)*b;
        }
        return *this;
    }
    

protected:
    int numx, numy;        // number of pixels
    double lenx, leny;        // dimensions of each pixel
    double origin_x, origin_y;  // coordinates of lower-left corner of first pixel
    std::vector <data_t> data;

    size_t pixel_index (int i, int j) const
    {
        return size_t (j+1) * size_t (numx+1) + size_t (i) + 1u;
    }

    data_t &at (int i, int j)
    {
        try{
            return data.at(pixel_index (i, j));
        } catch (const std::out_of_range&) {
            std::cerr << "Indices out of range! (i,j)=("<<i<<","<<j<<"), Photo has width "<<width ()<<" and height "<<height()<<std::endl;
            exit (EXIT_FAILURE);
        }
        return data.at(pixel_index (0, 0));
    }

    const data_t &at (int i, int j) const
    {
        try{
            return data.at(pixel_index (i, j));
        } catch (const std::out_of_range&) {
            std::cerr << "Indices out of range! (i,j)=("<<i<<","<<j<<"), Photo has width "<<width ()<<" and height "<<height()<<std::endl;
            exit (EXIT_FAILURE);
        }
        return data.at(pixel_index (0, 0));
    }
};

using Photo = BasicPhoto <double>;

// helper class to define transformed photos.
template <typename PHOTO>
struct PhotoAdapter
{
    const PHOTO &original;
    using data_t = typename PHOTO::data_t;

    PhotoAdapter (const PHOTO &ph)
        : original (ph)
    {
    }

    int width () const
    {
        return original.width ();
    }

    int height () const
    {
        return original.height ();
    }

    double pixel_width () const
    {
        return original.pixel_width ();
    }

    double pixel_height () const
    {
        return original.pixel_height ();
    }
};

// adapter to threshold a photo.
// returns photo >= threshold, i.e. ones and zeros.
// do not use directly, use make_thresholded_view.
template <typename PHOTO, typename THRESHOLD>
struct ThresholdingAdapter : PhotoAdapter <PHOTO>
{
    using PhotoAdapter <PHOTO>::original;
    using data_t = int;
    const THRESHOLD threshold;

    ThresholdingAdapter (const PHOTO &ph, const THRESHOLD &t)
        : PhotoAdapter <PHOTO> (ph), threshold (t)
    {
    }

    int operator() (int x, int y) const
    {
        return int (original(x,y) >= threshold);
    }
};

template <typename PHOTO, typename THRESHOLD>
auto make_thresholded_view (const PHOTO &p, const THRESHOLD &t) -> ThresholdingAdapter <PHOTO, THRESHOLD>
{
    return ThresholdingAdapter <PHOTO, THRESHOLD> (p, t);
}

// result (return type) of IMT computations.
// holds the values of the IMT's and, in addition, the area.
// FIXME make configurable what we actually want to compute
// FIXME proper accumulators
struct MinkowskiAccumulator
{
    static const int MAX_S = 8;

    MinkowskiAccumulator ()
    {
        area_ = peri_ = euler_ = 0.;
        for (int s : range(2, MAX_S+1))
            psi(s) = 0.;
    }

    // add areas.
    friend
    void add_triangle_area (MinkowskiAccumulator *acc, vec_t, vec_t v0, vec_t v1, vec_t v2)
    {
        acc->area_ += kahan_triangle_area(v0, v1, v2);
    }

    template <typename CONTAINER>
    friend
    void add_polygon_area (MinkowskiAccumulator *acc, const CONTAINER &vertices)
    {
        acc->area_ += shoelace_formula(vertices);
    }

    // add a piece of straight-edge contour.
    friend
    void add_contour_segment(MinkowskiAccumulator *acc, vec_t offset, vec_t begin, vec_t end)
    {
        (void)offset;
        double const p = std::hypot(end[0]-begin[0], end[1]-begin[1]);
        if (p == 0.) return;
        complex_t n = { (end[1]-begin[1]) / p, -(end[0]-begin[0]) / p };
        acc->peri_ += p;
        for (int s : range(2, MAX_S+1))
            acc->psi(s) += p * std::pow(n, s);
    }
    
    friend
    void add_corner(MinkowskiAccumulator *acc, double angleOverPi)
    {
        acc-> euler_ += 1./(2.)*angleOverPi;
    }

    // accessors to read data
    double area () const { return area_; }
    double perimeter () const { return peri_; }
    double euler () const { return euler_; }

    double msm (int s) const
    {
        s = std::abs(s);
        if (s >= 2 && s <= MAX_S)
            return std::abs(psi(s)) / peri_;
        else
            throw std::logic_error("msm(" + std::to_string (s) + ") called");
    }

    complex_t imt (int s) const
    {
        if (s == 0)
            return perimeter();
        else if (s < 0)
            return std::conj(imt(-s));
        else if (s >= 2 && s <= MAX_S)
            return psi(s);
        else
            throw std::logic_error("imt(" + std::to_string(s) + ") called");
    }

private:
    static
    double kahan_triangle_area(vec_t v0, vec_t v1, vec_t v2)
    {
        double a = norm(v0 - v1);
        double b = norm(v1 - v2);
        double c = norm(v2 - v0);
        // sort
        if (a < b) std::swap(a, b);
        if (b < c) std::swap(b, c);
        if (a < b) std::swap(a, b);
        double x = (a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c));
        return .25 * std::sqrt(x);
    }

    template <typename CONTAINER>
    static
    double shoelace_formula(const CONTAINER &vertices)
    {
        const int N = vertices.size();
        double twice_area = 0.;
        // shoelace formula for the area
        for (int i = 0; i != N; ++i)
        {
            auto const &begin = vertices.at(i);
            auto const &end   = vertices.at((i+1) % N);
            twice_area += begin[0] * end[1] - end[0] * begin[1];
        }
        return .5 * twice_area;
    }

    double area_, peri_, euler_;
    complex_t &psi(int s) { return int_psi_.at(s-2); }
    complex_t psi(int s) const { return int_psi_.at(s-2); }
    std::array<complex_t, MAX_S-1> int_psi_;
};

// utility class which can be used with marching squares to
// write an isocontour to a file which may be plotted
// by plot "..." w vec
struct GnuplottableContour
{
    GnuplottableContour (const string &filename) : df(filename) {}
    GnuplottableContour (std::ostream &os) : df(os) {}

    friend
    void add_contour_segment (GnuplottableContour *sink, vec_t off, vec_t begin, vec_t end)
    {
        sink->df
            << (off[0] + begin[0])
            << (off[1] + begin[1])
            << (end[0] - begin[0])
            << (end[1] - begin[1]) << std::endl;
    }

private:
    Datafile df;
};

// add a polygonal contour.
// by default, decompose it into edges.
template <typename SINK, typename CONTAINER>
void add_polygon_contour(SINK *sink, const CONTAINER &vertices)
{
    const int N = vertices.size();
    for (int i = 0; i != N; ++i)
    {
        auto const &begin = vertices.at(i);
        auto const &end   = vertices.at((i+1) % N);
        add_contour_segment(sink, {0., 0.}, begin, end);
    }
}

// add a triangular area.
// by default, do nothing.
template <typename SINK>
void add_triangle_area (SINK *, vec_t offset, vec_t v0, vec_t v1, vec_t v2)
{
    (void)offset;
    (void)v0;
    (void)v1;
    (void)v2;
}

// add a polygonal area.
// by default, do nothing.
template <typename SINK, typename CONTAINER>
void add_polygon_area (SINK *sink, const CONTAINER &vertices)
{
    (void)sink;
    (void)vertices;
}

// interpolated marching squares, core routine handling a 2x2 neighborhood
// FIXME this is missing the curvature measures
// FIXME we might want to add extra logic to break the
//       arbitrary choice in cases 6 and 9
template <typename SINK>
void add_interpolated_four_neighborhood (SINK *sink,
    vec_t const &off, vec_t const &pix_diag,
    double ll, double ul, double lr, double ur, double threshold)
{
    unsigned lut_index = (ll>=threshold) * 1 + (ul>=threshold) * 2
                       + (lr>=threshold) * 4 + (ur>=threshold) * 8;

    // interpolate between 0.5 and 1.5.
    auto msq_interp = [] (double left, double threshold, double right) -> double
    {
        return (threshold - left) / (right - left) + .5;
    };

    // vertices, located on the edges of the neighborhood (NOT pixels)
    // i.e. on lines connecting pixel centers
    vec_t lower = { msq_interp(ll, threshold, lr), 0.5 };
    vec_t upper = { msq_interp(ul, threshold, ur), 1.5 };
    vec_t left  = { 0.5, msq_interp(ll, threshold, ul) };
    vec_t right = { 1.5, msq_interp(lr, threshold, ur) };
    vec_t sw    = { 0.5, 0.5 };
    vec_t se    = { 1.5, 0.5 };
    vec_t nw    = { 0.5, 1.5 };
    vec_t ne    = { 1.5, 1.5 };
    // scale to physical coordinates
    lower = elementwise_product(lower, pix_diag);
    upper = elementwise_product(upper, pix_diag);
    left  = elementwise_product(left , pix_diag);
    right = elementwise_product(right, pix_diag);
    sw    = elementwise_product(sw,    pix_diag);
    se    = elementwise_product(se,    pix_diag);
    nw    = elementwise_product(nw,    pix_diag);
    ne    = elementwise_product(ne,    pix_diag);
    switch (lut_index)
    {
    case 0:
        // no area
        // no perimeter
        break;
    case 1:
        add_triangle_area(sink, off, lower, left, sw); // sw corner
        add_contour_segment(sink, off, lower, left);
        add_corner(sink, 0.5);
        break;
    case 2:
        add_triangle_area(sink, off, left, upper, nw); // nw corner
        add_contour_segment(sink, off, left, upper);
        add_corner(sink, 0.5);
        break;
    case 3:
        add_triangle_area(sink, off, lower, upper, nw); // left half
        add_triangle_area(sink, off, nw, sw, lower);
        add_contour_segment(sink, off, lower, upper);
        break;
    case 4:
        add_triangle_area(sink, off, right, lower, se); // se corner
        add_contour_segment(sink, off, right, lower);
        add_corner(sink, 0.5);
        break;
    case 5:
        add_triangle_area(sink, off, right, left, sw); // lower half
        add_triangle_area(sink, off, sw, se, right);
        add_contour_segment(sink, off, right, left);
        break;
    case 6:
        add_triangle_area(sink, off, nw, left, upper); //Diagonal structure, interpreted as separated holes
        add_triangle_area(sink, off, left, lower, upper);
        add_triangle_area(sink, off, upper, lower, right);
        add_triangle_area(sink, off, lower, se, right);
        add_contour_segment(sink, off, left, lower);
        add_contour_segment(sink, off, right, upper);
        add_corner(sink, -1);
        break;
    case 7:
        add_triangle_area(sink, off, upper, nw, sw); // missing ne corner
        add_triangle_area(sink, off, right, upper, sw);
        add_triangle_area(sink, off, se, right, sw);
        add_contour_segment(sink, off, right, upper);
        add_corner(sink, -0.5);
        break;
    case 8:
        add_triangle_area(sink, off, upper, right, ne); // ne corner
        add_contour_segment(sink, off, upper, right);
        add_corner(sink, 0.5);
        break;
    case 9:
        add_triangle_area(sink, off, sw, lower, left); //Diagonal structure, interpreted as separated holes
        add_triangle_area(sink, off, lower, right, left);
        add_triangle_area(sink, off, left, right, upper);
        add_triangle_area(sink, off, upper, right, ne);
        add_contour_segment(sink, off, lower, right);
        add_contour_segment(sink, off, upper, left);
        add_corner(sink, -1);
        break;
    case 10:
        add_triangle_area(sink, off, nw, left, right); // upper half
        add_triangle_area(sink, off, right, ne, nw);
        add_contour_segment(sink, off, left, right);
        break;
    case 11:
        add_triangle_area(sink, off, right, ne, nw);  // missing se corner
        add_triangle_area(sink, off, lower, right, nw);
        add_triangle_area(sink, off, sw, lower, nw);
        add_contour_segment(sink, off, lower, right);
        add_corner(sink, -0.5);
        break;
    case 12:
        add_triangle_area(sink, off, upper, lower, se); // right half
        add_triangle_area(sink, off, se, ne, upper);
        add_contour_segment(sink, off, upper, lower);
        break;
    case 13:
        add_triangle_area(sink, off, ne, upper, se); // missing nw corner
        add_triangle_area(sink, off, upper, left, se);
        add_triangle_area(sink, off, left, sw, se);
        add_contour_segment(sink, off, upper, left);
        add_corner(sink, -0.5);
        break;
    case 14:
        add_triangle_area(sink, off, nw, left, ne); // missing sw corner
        add_triangle_area(sink, off, left, lower, ne);
        add_triangle_area(sink, off, lower, se, ne);
        add_contour_segment(sink, off, left, lower);
        add_corner(sink, -0.5);
        break;
    case 15:
        add_triangle_area(sink, off, sw, se, ne); // full square
        add_triangle_area(sink, off, ne, nw, sw);
        // no perimeter
        break;
    default:
        std::abort();
    }
}


// interpolated marching squares, loop over the whole image
// effectively calls add_area and add_perimeter on "sink" once
// a piece of the contour has been identified.
// processes the whole photo, (w-2)x(h-2) neighborhoods, or
// w x h if padded is set to true.  padding pixels are taken to
// be black, intensity zero.
template <typename SINK, typename PHOTO, typename THRESHOLD>
void trace_isocontour_interpolated_marching_squares (SINK *sink,
    const PHOTO &ph, const THRESHOLD &threshold, bool padded = false)
{
    const int start = -padded;
    const int width = ph.width () + padded;
    const int height = ph.height () + padded;

    const vec_t pixdiag = { ph.pixel_width (), ph.pixel_height () };

    for (int j = start; j < height-1; ++j)
    for (int i = start; i < width-1; ++i)
    {
        vec_t off = { i * ph.pixel_width (), j * ph.pixel_height () };
        add_interpolated_four_neighborhood (sink, off, pixdiag,
            ph(i,j), ph(i,j+1), ph(i+1,j), ph(i+1,j+1), threshold);
    }
}

// convenience wrapper to compute IMT's with interpolated marching squares
template <typename PHOTO>
MinkowskiAccumulator imt_interpolated_marching_squares (
    const PHOTO &ph, double threshold, bool padded = false)
{
    MinkowskiAccumulator acc;
    trace_isocontour_interpolated_marching_squares (&acc, ph,
        threshold, padded);
    return acc;
}

// reproduce (inferior) results of regular marching squares
// by applying interpolated marching squares on a binarized image
template <typename SINK, typename PHOTO, typename THRESHOLD>
void trace_isocontour_regular_marching_squares (SINK *sink,
    const PHOTO &ph, const THRESHOLD &threshold, bool padded = false)
{
    auto thr_ph = make_thresholded_view (ph, threshold);
    trace_isocontour_interpolated_marching_squares (sink, thr_ph, .5, padded);
}

// convenience wrapper to compute IMT's with regular marching squares
template <typename PHOTO>
MinkowskiAccumulator imt_regular_marching_squares (
    const PHOTO &ph, double threshold, bool padded = false)
{
    MinkowskiAccumulator acc;
    trace_isocontour_regular_marching_squares (&acc, ph,
        threshold, padded);
    return acc;
}

// compute the IMTs of a single polygon
template <typename CONTAINER>
MinkowskiAccumulator imt_polygon (const CONTAINER &vertices)
{
    MinkowskiAccumulator acc;
    add_polygon_area(&acc, vertices);
    add_polygon_contour(&acc, vertices);
    return acc;
}

} // namespace papaya2
