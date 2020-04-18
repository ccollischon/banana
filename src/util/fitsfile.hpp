#ifndef _banana_fitsfile_hpp_
#define _banana_fitsfile_hpp_

#include <algorithm>

#include "../papaya2.hpp"

extern double crpix_source;

struct FitsFile : papaya2::Photo
{
    std::vector<std::string> WCSkeynames_here;
    std::vector<std::string> WCSvalues;

    FitsFile (const std::string &infilename, int layer = 0, unsigned hdu = 0);

    FitsFile (const std::string &infilename, std::vector<std::string>& WCSkeynames, int layer = 0, unsigned hdu = 0);

    FitsFile (const std::vector<std::vector<double>> &map, const std::vector<std::vector<std::string>> &WCSdata); //Create a FitsFile from 2D-vector containing image

    // allow overwriting pixels by the user
    inline double &at (int i, int j) { return papaya2::Photo::at (i, j); }

    // return WCS data in given order
    inline std::vector<std::vector<std::string>> returnWCSdata()
    {
        std::vector<std::vector<std::string>> ret;
        ret.push_back(WCSkeynames_here);
        ret.push_back(WCSvalues);
        return ret;
    }

    inline void addnoise () { at (2,2) = 1; }

    std::string giveKeyvalue(std::string name)
    {
        std::string content = "0";
        auto it = std::find(WCSkeynames_here.begin(), WCSkeynames_here.end(), name);
        if (it != WCSkeynames_here.end())
        {
            int pos1 = distance(WCSkeynames_here.begin(), it);
            content = WCSvalues.at(pos1);
            //std::cout << WCSkeynames_here.at(pos1) << " " << content << std::endl;

        }
        return content;
    }

    void setKeyValue(std::string name, std::string newValue)
    {
        auto it = std::find(WCSkeynames_here.begin(), WCSkeynames_here.end(), name);
        if (it != WCSkeynames_here.end())
        {
            int pos1 = distance(WCSkeynames_here.begin(), it);
            WCSvalues.at(pos1) = newValue;
            //std::cout << WCSkeynames_here.at(pos1) << " " << content << std::endl;
        }
    }

    FitsFile &operator-=(FitsFile& b) //subtract two files with equal CDELT at correct sky positions
    {
        std::cout << "Subtracting FITS files..." << std::endl;
        int w = std::min(numx, b.width());
        int h = std::min(numy, b.height());
        int CRPIX1_here = 0, CRPIX1_there = 0;
        std::vector<std::vector<std::string>> WCS_there = b.returnWCSdata();

        //Find positions of CRPIX1 in both images. Only works for equal CDELT1/2
        auto it = std::find(WCSkeynames_here.begin(), WCSkeynames_here.end(), "CRPIX1");
        if (it != WCSkeynames_here.end())
        {
            int pos1 = distance(WCSkeynames_here.begin(), it);
            CRPIX1_here = std::stoi(WCSvalues.at(pos1));
            std::cout << WCSkeynames_here.at(pos1) << " " << CRPIX1_here << std::endl;

        }
        it = std::find(WCS_there.at(0).begin(), WCS_there.at(0).end(), "CRPIX1");
        if (it != WCS_there.at(0).end())
        {
            int pos2 = distance(WCS_there.at(0).begin(), it);
            CRPIX1_there = std::stoi(WCS_there.at(1).at(pos2));
            std::cout << WCS_there.at(0).at(pos2) << " " << CRPIX1_there << std::endl;
        }

        //Calculate shift between images. Only works for equal CDELT1/2
        int shift = CRPIX1_here-CRPIX1_there;
        //Subtract
        for(int i=std::max(shift,0); i<w-std::max(-1*shift,0); i++)
        for(int j=std::max(shift,0); j<h-std::max(-1*shift,0); j++)
        {
            at(i,j) = at(i,j)-b.at((i-shift),(j-shift));
        }
        return *this;
    }

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

    double atCoord(double ra, double dec) //DOES NOT WORK!! TODO FIXME
    //returns pixel value at given sky coordinates in deg
    {
        double value = 0;

        //Calculate distance to ref pixel in deg
        double CRVAL1 = 0., CRVAL2 = 0.;
        CRVAL1 = std::stod(giveKeyvalue("CRVAL1"));
        CRVAL2 = std::stod(giveKeyvalue("CRVAL2"));
        double distdeg1, distdeg2;
        distdeg1 = ra-CRVAL1;
        distdeg2 = dec-CRVAL2;
        std::cout << distdeg1 << " " << distdeg2 << std::endl;

        //...and in pixels
        double CDELT1 = std::stod(giveKeyvalue("CDELT1"));
        double CDELT2 = std::stod(giveKeyvalue("CDELT2"));
        double distpix1, distpix2;
        distpix1 = distdeg1/CDELT1;
        distpix2 = distdeg2/CDELT2;
        std::cout << distpix1 << " " << distpix2 << std::endl;

        //add to ref pixel (pay attention to different y axes here and in ds9!)
        double CRPIX1 = std::stod(giveKeyvalue("CRPIX1"));
        double CRPIX2 = std::stod(giveKeyvalue("CRPIX2"));
        int pix1, pix2;
        pix1 = CRPIX1 + distpix1;
        pix2 = CRPIX2 + distpix2;
        std::cout << pix1 << " " << pix2 << std::endl;

        //take value of pixel here and return
        value = at(pix1, pix2);
        return value;
    }

    double atds9pix(double pix1, double pix2) //returns value at given pixel (ds9), accounts for shift
    {
        double value = 0;
        //Calculate shift between images. Only works for equal CDELT1/2
        double CRPIX1_here = std::stod(giveKeyvalue("CRPIX1"));
        double shift = CRPIX1_here-crpix_source;//1674.4;

        //std::cout << pix1+shift << " " << pix2-shift << std::endl;

        if( (pix1+shift-1)<width() && (pix1+shift-1)>0 && (height()-pix2+shift)<height() && (height()-pix2+shift)>0 )
            value = at((int)(pix1+shift-1),(int)(height()-pix2+shift));

        return value;
    }

    std::vector<int> giveUnshiftedds9Coord(double pix1, double pix2) //Takes FitsFile coordinates of some shifted minkmap, returns pixel coordinate in ds9 Halpha (to be used for ds9 region files)
    {
        double CRPIX1_here = std::stod(giveKeyvalue("CRPIX1"));
        double shift = crpix_source-CRPIX1_here;
        std::vector<int> output;
        output.push_back((int)(pix1+shift-1));
        output.push_back((int)(height()-pix2-shift));
        return output;
    }

    std::vector<int> giveShiftedCoord(double pix1, double pix2) // Takes ds9 coordinates of unshifted image and returns pixel coordinates of shifted FitsFile
    {
        double CRPIX1_here = std::stod(giveKeyvalue("CRPIX1"));
        double shift = -(crpix_source-CRPIX1_here);
        std::vector<int> output;
        output.push_back((int)(pix1+shift-1));
        output.push_back((int)(height()-pix2+shift));
        return output;
    }

    double giveShift() //Subtract from x and add to y for going from original to shifted
    {
        double CRPIX1_here = std::stod(giveKeyvalue("CRPIX1"));
        double shift = crpix_source-CRPIX1_here;
        return shift;
    }

    void abs()
    {
        for(int i=1; i<numx; i++)
        	for(int j=1; j<numy; j++)
				at(i,j) = std::abs(at(i,j));
    }

};

template<typename PHOTO> //Writes one fits file from photo. Possible to write absolute value, argument, flipped image in any axis
void writeImage(PHOTO minkmap, std::string filename,
		const std::vector<std::vector<std::string>>& WCSdata,
		bool absolute = true, bool arg = false, bool flipY = true,
		bool flipX = false, bool highPrec = false);

template<typename PHOTO> //Writes one fits file from vector of photos. Possible to write absolute value, argument, flipped image in any axis
void write3Dimage(const std::vector<PHOTO>& minkmaps, std::string filename,
		const std::vector<std::vector<std::string>>& WCSdata, bool absolute =
		 true, bool arg = false, bool flipY = true,
		bool flipX = false);

void erosionDilation(FitsFile& infile, int diameter);

double sumOverSquares(FitsFile& image, int smooth);

void takeLogPositive(FitsFile& image);

void binImage(FitsFile& infile, int factor /*= 2*/);

#endif
