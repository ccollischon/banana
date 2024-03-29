#ifndef _banana_fitsfile_hpp_
#define _banana_fitsfile_hpp_

#include <algorithm>
#include <map>

#include "../papaya2.hpp"

extern double crpix_source;

struct FitsFile : papaya2::Photo
{
    std::map<std::string, std::string> WCSdata_;

    FitsFile (const std::string &infilename, int layer = 0, unsigned hdu = 0);

    FitsFile (const std::string &infilename, const std::vector<std::string>& WCSkeynames, int layer = 0, unsigned hdu = 0);

    FitsFile (const std::vector<std::vector<double>> &map, std::map<std::string, std::string> WCSdata); //Create a FitsFile from 2D-vector containing image

    // allow overwriting pixels by the user
    inline double &at (int i, int j) { return papaya2::Photo::at (i, j); }

    // return WCS data in given order
    inline std::map<std::string, std::string> returnWCSdata() const
    {
        return WCSdata_;
    }

    std::string giveKeyvalue(const std::string& name) const;
    
    void setKeyValue(const std::string& name, const std::string& newValue);
    

    FitsFile &operator-=(const FitsFile& b); //subtract two files with equal CDELT at correct sky positions

    FitsFile &operator-=(double b)
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


    inline double atds9pix(double pix1, double pix2) const //returns value at given pixel (ds9), accounts for shift
    {
        double value = 0;
        //Calculate shift between images. Only works for equal CDELT1/2
        double shift = giveShift();//1674.4;

        //std::cout << pix1+shift << " " << pix2-shift << std::endl;

        if( (pix1+shift-1)<width() && (pix1+shift-1)>0 && (height()-pix2+shift)<height() && (height()-pix2+shift)>0 )
            value = operator ()((int)(pix1+shift-1),(int)(height()-pix2+shift));

        return value;
    }

    inline std::vector<int> giveUnshiftedds9Coord(double pix1, double pix2) const //Takes FitsFile coordinates of some shifted minkmap, returns unshifted pixel coordinate in ds9 (to be used for ds9 region files)
    {
        double shift = giveShift();
        std::vector<int> output;
        output.push_back((int)(pix1+shift-1));
        output.push_back((int)(height()-pix2-shift));
        return output;
    }

    inline std::vector<int> giveShiftedCoord(double pix1, double pix2) const // Takes ds9 coordinates of unshifted image and returns pixel coordinates of shifted FitsFile
    {
        double shift = -giveShift();
        std::vector<int> output;
        output.push_back((int)(pix1+shift-1));
        output.push_back((int)(height()-pix2+shift));
        return output;
    }

    inline double giveShift() const //Subtract from x and add to y for going from original to shifted
    {
        double CRPIX1_here = std::stod(giveKeyvalue("CRPIX1"));
        double shift = crpix_source-CRPIX1_here;
        return shift;
    }

    inline void abs()
    {
        for(int i=1; i<numx; i++)
        	for(int j=1; j<numy; j++)
				at(i,j) = std::abs(at(i,j));
    }

};

template<typename PHOTO> //Writes one fits file from photo. Possible to write absolute value, argument, flipped image in any axis
void writeImage(const PHOTO& minkmap, std::string filename,
		const std::map<std::string, std::string>& WCSdata,
		bool absolute = true, bool arg = false, bool flipY = true,
		bool flipX = false, bool highPrec = false);

template<typename PHOTO> //Writes one fits file from vector of photos. Possible to write absolute value, argument, flipped image in any axis
void write3Dimage(const std::vector<PHOTO>& minkmaps, std::string filename,
		const std::map<std::string, std::string>& WCSdata, bool absolute =
		 true, bool arg = false, bool flipY = true,
		bool flipX = false);

void erosionDilation(FitsFile& infile, int diameter);

double sumOverSquares(const FitsFile& image, int smooth);

void takeLogPositive(FitsFile& image);

void binImage(FitsFile& infile, int factor /*= 2*/);

#endif
