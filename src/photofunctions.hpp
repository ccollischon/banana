#ifndef _banana_photofunctions_hpp_
#define _banana_photofunctions_hpp_

#include "src/papaya2.hpp"

using namespace papaya2;

template<typename PHOTO>
void smooth_map_plus (PHOTO& infile, int squaresize)
{
    PHOTO original = infile;
    //infile.set_coordinates(0, 0, (infile.width()-squaresize), (infile.height()-squaresize), (infile.width()-squaresize), (infile.height()-squaresize));
    infile *= 0;
    double rsquared = squaresize*squaresize/4;
    double sqrtarea = pow(rsquared*3.14159,0.5);
    int begin1 = 0,begin2 = 0;
    int end1 = infile.height()-squaresize+1;
    int end2 = infile.width()-squaresize+1;
    
    int stepsize = std::max(squaresize/6,1);
    
    #pragma omp parallel for
    for (int j = begin1; j < end1; j+=stepsize)
    for (int i = begin2; i < end2; i+=stepsize)
    {
        typename PHOTO::data_t count = 0;
        double n = 0;
        
        for(int k=0; k<(squaresize-1); k++)
        for(int l=0; l<(squaresize-1); l++)
        {
            double dist = pow(((squaresize/2.) - (k+1)),2)+pow(((squaresize/2.) - (l+1)),2);
            if( dist < rsquared-1 )
            {
                double mul = rsquared-1-dist;
                count += original(i+k,j+l)*mul;
                n += mul;
            }
        }
        typename PHOTO::data_t newValue = count/n*sqrtarea;
        for(int kstep=0; kstep<stepsize; kstep++)
        for(int lstep=0; lstep<stepsize; lstep++)
        {
            infile(i+kstep,j+lstep) = newValue;
        }
    }
}

template<typename PHOTO>
void smooth_map (PHOTO& infile)
{
    for (int j = 1; j < infile.height()-1; ++j)
    for (int i = 1; i < infile.width()-1; ++i)
    {
        infile(i,j) = (infile(i,j) + infile(i+1,j)/2. + infile(i-1,j)/2. + infile(i,j+1)/2. + infile(i,j-1)/2. + infile(i+1,j+1)/4.+ infile(i+1,j-1)/4. + infile(i-1,j+1)/4. + infile(i-1,j-1)/4. )/4.;
    }
}


template <typename PHOTO>
PHOTO averagemaps(const std::vector<PHOTO>& datamaps)
{
    std::cout << "Averaging ..." << std::endl;
    PHOTO returnPhoto = datamaps.at(0);
    std::cout << "Map 0 \n";
    for(uint i=1; i<datamaps.size(); i++)
    {
        returnPhoto += datamaps.at(i);
        std::cout << "Map " << i << std::endl;
    }
    returnPhoto/=datamaps.size();
    return returnPhoto;
}

template <typename PHOTO>
PHOTO averageAbs(const std::vector<PHOTO>& datamaps)
{
    int w = datamaps.at(0).width();
    int h = datamaps.at(0).height();
    std::cout << "Averaging absolute values..." << std::endl;
    PHOTO returnPhoto = datamaps.at(0);
    for(int k=0; k<w; k++)
    for(int j=0; j<h; j++)
    {
        returnPhoto(k,j) = std::abs(returnPhoto(k,j));
    }
    std::cout << "Map 0" << std::endl;
    
    for(uint i=1; i<datamaps.size(); i++)
    {
        for(int k=0; k<w; k++)
        for(int j=0; j<h; j++)
        {
            returnPhoto(k,j) = returnPhoto(k,j)+std::abs(datamaps.at(i)(k,j));
        }
        std::cout << "Map " << i << std::endl;
    }
    returnPhoto /= datamaps.size();
    return returnPhoto;
}



template<typename PHOTO>
void maskBoxes(PHOTO& infile, std::vector<std::vector<double>>& boxes)
{
    for(uint_fast8_t i=0; i<boxes.size(); i++)
    {
        int centerDEC = (int)boxes.at(i).at(1);
        int centerRA = (int)boxes.at(i).at(0);
        int widthDEC = (int)boxes.at(i).at(3)+1;
        int widthRA = (int)boxes.at(i).at(2)+1;
        
        int begin1 = std::max(0, infile.height()-centerDEC-widthDEC/2); // DEC
        int end1 = std::min(infile.height(),infile.height()-centerDEC+widthDEC/2);
        int begin2 = std::max(0, centerRA-widthRA/2); //RA
        int end2 = std::min(infile.width(),centerRA+widthRA/2);
        
        for (int j = begin1; j < end1; j++)
        for (int k = begin2; k < end2; k++)
        {
            infile(k,j) = 0;
        }
    }
}

template<typename PHOTO>
BasicPhoto<double> getMask(PHOTO& infile, std::vector<std::vector<double>>& boxes)
{
    BasicPhoto<double> mask;
    mask.set_coordinates (0, 0, 1, 1, infile.width(), infile.height());
    for(uint_fast8_t i=0; i<boxes.size(); i++)
    {
        int centerDEC = (int)boxes.at(i).at(1);
        int centerRA = (int)boxes.at(i).at(0);
        int widthDEC = (int)boxes.at(i).at(3);
        int widthRA = (int)boxes.at(i).at(2);
        
        int begin1 = std::max(0, infile.height()-centerDEC-widthDEC/2); // DEC
        int end1 = std::min(infile.height(),infile.height()-centerDEC+widthDEC/2);
        int begin2 = std::max(0, centerRA-widthRA/2); //RA
        int end2 = std::min(infile.width(),centerRA+widthRA/2);
        for (int j = std::max(0,begin1); j < end1; j++)
        for (int k = std::max(0,begin2); k < end2; k++)
        {
            mask(k,j) = 1;
        }
    }
    return mask;
}

template<typename PHOTO>
void includeBoxes(PHOTO& infile, std::vector<std::vector<double>>& boxes)
{
    BasicPhoto<double> mask = getMask(infile, boxes);
    infile *= mask;
}

template<typename PHOTO>
double averageBoxes(PHOTO& infile, std::vector<std::vector<double>>& boxes)
{
    //std::cout << "Getting mask...\n";
    BasicPhoto<double> mask = getMask(infile, boxes);
    
    //std::cout << "Got mask, averaging...\n";
    double value = 0.;
    int count = 0;
    for(int i=0; i<infile.width(); i++)
    for(int j=0; j<infile.height(); j++)
    {
        value += mask(i,j)*infile(i,j);
        count += mask(i,j);
    }
    //std::cout << "Finished averaging...\n";
    return value/count;
}

template<typename PHOTO>
PHOTO cutout(PHOTO& infile, std::vector<std::vector<double>>& boxes)
{
    BasicPhoto<double> mask = getMask(infile, boxes);
    
    int min1 = 99999;
    int max1 = 0;
    int min2 = 99999;
    int max2 = 0;
    for(uint_fast8_t i=0; i<boxes.size(); i++)
    {
        int centerDEC = (int)boxes.at(i).at(1);
        int centerRA = (int)boxes.at(i).at(0);
        int widthDEC = (int)boxes.at(i).at(3);
        int widthRA = (int)boxes.at(i).at(2);
        
        int begin1 = std::max(0, infile.height()-centerDEC-widthDEC/2); // DEC
        int end1 = std::min(infile.height(),infile.height()-centerDEC+widthDEC/2);
        int begin2 = std::max(0, centerRA-widthRA/2); //RA
        int end2 = std::min(infile.width(),centerRA+widthRA/2);
        
        //Find maximum width and height of image
        if(begin1 < min1) min1 = begin1;
        if(begin2 < min2) min2 = begin2;
        if(end1 > max1) max1 = end1;
        if(end2 > max2) max2 = end2;
    }
    
    PHOTO outfile = infile;
    outfile.set_coordinates (0, 0, 1, 1, max2-min2, max1-min1);

    for (int j = min1; j < max1; j++)
    for (int k = min2; k < max2; k++)
    {
        outfile(k-min2,j-min1) = infile(k,j)*mask(k,j);
    }   
    return outfile;
}


typedef BasicPhoto<complex_t> complex_image;

template <typename PHOTO, typename THRESHOLD>
void minkowski_map_interpolated_marching_squares (complex_image *out,
    const PHOTO &ph, const THRESHOLD &threshold, int s, bool padded = false)
{
    const int start = -padded;
    const int width = ph.width () + padded;
    const int height = ph.height () + padded;
    
    const vec_t pixdiag = { ph.pixel_width (), ph.pixel_height () };
    const vec_t half_a_pixdiag = pixdiag / 2;
    vec_t mmap_origin = ph.origin() + half_a_pixdiag;
    if (padded) mmap_origin -= pixdiag;
    vec_t mmap_upperright = ph.upper_right() - half_a_pixdiag;
    if (padded) mmap_upperright += pixdiag;
    out->set_coordinates(mmap_origin[0], mmap_origin[1], mmap_upperright[0], mmap_upperright[1], width-start-1, height-start-1);

    if(s == 0)
    {
        #pragma omp parallel for
        for (int j = start; j < height-1; ++j)
        for (int i = start; i < width-1; ++i)
        {
            vec_t off = { i * ph.pixel_width (), j * ph.pixel_height () };
            MinkowskiAccumulator minkval;
            add_interpolated_four_neighborhood (&minkval, off, pixdiag,
                ph(i,j), ph(i,j+1), ph(i+1,j), ph(i+1,j+1), threshold);
            (*out)(i+padded, j+padded) = minkval.perimeter();
        
        }
    }
    else if(s == 4234)
    {
        #pragma omp parallel for
        for (int j = start; j < height-1; ++j)
        for (int i = start; i < width-1; ++i)
        {
            vec_t off = { i * ph.pixel_width (), j * ph.pixel_height () };
            MinkowskiAccumulator minkval;
            add_interpolated_four_neighborhood (&minkval, off, pixdiag,
                ph(i,j), ph(i,j+1), ph(i+1,j), ph(i+1,j+1), threshold);
            (*out)(i+padded, j+padded) = minkval.area();
        
        }
    }
    else
    {
        #pragma omp parallel for
        for (int j = start; j < height-1; ++j)
        for (int i = start; i < width-1; ++i)
        {
            vec_t off = { i * ph.pixel_width (), j * ph.pixel_height () };
            MinkowskiAccumulator minkval;
            add_interpolated_four_neighborhood (&minkval, off, pixdiag,
                ph(i,j), ph(i,j+1), ph(i+1,j), ph(i+1,j+1), threshold);
            if (s==1) {
                (*out)(i+padded, j+padded) = minkval.euler();
                //std::cout << minkval.euler() << std::endl;
            }
            else (*out)(i+padded, j+padded) = minkval.imt(s);
            
        }
    }
}

template <typename PHOTO, typename THRESHOLD>
void minkowski_map_interpolated_marching_bigsquares (complex_image *out,
    const PHOTO &ph, const THRESHOLD &threshold, int s, int squaresize, bool padded = false)
{
    const int start = -padded;
    const int width = ph.width () + padded;
    const int height = ph.height () + padded;
    
    const vec_t pixdiag = { ph.pixel_width (), ph.pixel_height () };
    const vec_t half_a_pixdiag = pixdiag / 2;
    vec_t mmap_origin = ph.origin() + half_a_pixdiag*0.5*squaresize;
    if (padded) mmap_origin -= pixdiag*0.5*squaresize;
    vec_t mmap_upperright = ph.upper_right() - half_a_pixdiag*0.5*squaresize;
    if (padded) mmap_upperright += pixdiag*0.5*squaresize;
    out->set_coordinates(mmap_origin[0], mmap_origin[1], mmap_upperright[0], mmap_upperright[1], width-start-1, height-start-1);


    for (int j = start; j < (height-squaresize+1); ++j)
    for (int i = start; i < (width-squaresize+1); ++i)
    {
        if(j%100 == 0 && i == 1) std::cout << "line: " << j << std::endl;
        vec_t off = { i * ph.pixel_width (), j * ph.pixel_height () };
        MinkowskiAccumulator minkval;
        
        for(int k=0; k<(squaresize-1); k++)
        for(int l=0; l<(squaresize-1); l++)
        {
            add_interpolated_four_neighborhood (&minkval, off, pixdiag,
                ph(i+k,j+l), ph(i+k,j+l+1), ph(i+k+1,j+l), ph(i+k+1,j+l+1), threshold);
        }
        if (s==0) (*out)(i+padded, j+padded) = minkval.perimeter();
        else if (s==1) (*out)(i+padded, j+padded) = minkval.euler();
        else (*out)(i+padded, j+padded) = minkval.imt(s);
        
    }
}

template <typename PHOTO, typename THRESHOLD>
void minkowski_map_interpolated_marching_bigcircles (complex_image *out,
    const PHOTO &ph, const THRESHOLD &threshold, int s, int squaresize, bool padded = false)
{
    double rsquared = squaresize*squaresize/4;
    const int start = -padded;
    const int width = ph.width () + padded;
    const int height = ph.height () + padded;
    
    const vec_t pixdiag = { ph.pixel_width (), ph.pixel_height () };
    const vec_t half_a_pixdiag = pixdiag / 2;
    vec_t mmap_origin = ph.origin() + half_a_pixdiag*0.5*squaresize;
    if (padded) mmap_origin -= pixdiag*0.5*squaresize;
    vec_t mmap_upperright = ph.upper_right() - half_a_pixdiag*0.5*squaresize;
    if (padded) mmap_upperright += pixdiag*0.5*squaresize;
    out->set_coordinates(mmap_origin[0], mmap_origin[1], mmap_upperright[0], mmap_upperright[1], width-start-1, height-start-1);


    for (int j = start; j < (height-squaresize+1); ++j)
    for (int i = start; i < (width-squaresize+1); ++i)
    {
        if(j%100 == 0 && i == 1) std::cout << "line: " << j << std::endl;
        vec_t off = { i * ph.pixel_width (), j * ph.pixel_height () };
        MinkowskiAccumulator minkval;
        
        for(int k=0; k<(squaresize-1); k++)
        for(int l=0; l<(squaresize-1); l++)
        {
            if( pow(((squaresize/2.) - (k+1)),2)+pow(((squaresize/2.) - (l+1)),2) < rsquared-1 )
                add_interpolated_four_neighborhood (&minkval, off, pixdiag,
                    ph(i+k,j+l), ph(i+k,j+l+1), ph(i+k+1,j+l), ph(i+k+1,j+l+1), threshold);
        }
        if (s==0) (*out)(i+padded, j+padded) = minkval.perimeter();
        else if (s==1) (*out)(i+padded, j+padded) = minkval.euler();
        else (*out)(i+padded, j+padded) = minkval.imt(s);
    }
}




#endif
