#ifndef _banana_fitsfile_hpp_

#define _banana_fitsfile_hpp_

#include "papaya2.hpp"
#include <CCfits/CCfits>
#include <algorithm>

extern double crpix_source;

using namespace papaya2;
struct FitsFile : Photo
{
    std::vector<string> WCSkeynames_here;
    std::vector<string> WCSvalues;
    
    FitsFile (const string &infilename, int layer = 0, unsigned hdu = 0)
    {
        if (hdu != 0)
            throw std::runtime_error ("hdu != 0 not implemented");
        CCfits::FITS::setVerboseMode (true);
        CCfits::FITS infile (infilename.c_str (), CCfits::Read, true);
        CCfits::PHDU &phdu = infile.pHDU ();
        
        phdu.readAllKeys ();
        set_coordinates (0, 0, phdu.axis (0), phdu.axis (1), phdu.axis (0), phdu.axis (1));
        //set_coordinates (0, 0, 1, 1, phdu.axis (0), phdu.axis (1));
        std::valarray <double> buffer;
        phdu.read (buffer);
        const int w = width ();
        const int h = height ();
        double max = -1e99;
        double min = 1e99;
        unsigned long offset = (layer*w*h);
        for (int i = 0; i < w; ++i)
        for (int j = 0; j < h; ++j)
        {
            at (i,j) = buffer[((h - j - 1) * w + i) + offset];
            max = fmax (max, at (i,j));
            min = fmin (min, at (i,j));
        }

        std::cerr << "info: loaded " << w << "x" << h
            << " FITS file, max lumin = " << max << ", min = "
            << min << "\n";
    }

    FitsFile (const string &infilename, std::vector<string>& WCSkeynames, int layer = 0, unsigned hdu = 0)
    {
        if (hdu != 0)
            throw std::runtime_error ("hdu != 0 not implemented");
        CCfits::FITS::setVerboseMode (true);
        CCfits::FITS infile (infilename.c_str (), CCfits::Read, true);
        CCfits::PHDU &phdu = infile.pHDU ();
        
        // FIXME extract proper sky coordinates from FITS and insert pixel dimensions here.
        // for now, just map everything on the unit square
        
        std::cout << "Loading WCS data..." << std::endl;
        WCSkeynames_here = WCSkeynames;
        for(auto keyword:WCSkeynames)
        {
            string keyval="";
            try{
                phdu.readKey(keyword,keyval);
            }
            catch(const CCfits::HDU::NoSuchKeyword&){
                keyval = "0.";
            }
            WCSvalues.push_back(keyval);
        }
        //phdu.readKeys(WCSkeynames,WCSvalues);
        
        set_coordinates (0, 0, phdu.axis (0), phdu.axis (1), phdu.axis (0), phdu.axis (1));
        //set_coordinates (0, 0, 1, 1, phdu.axis (0), phdu.axis (1));
        std::valarray <double> buffer;
        phdu.read (buffer);
        const int w = width ();
        const int h = height ();
        double max = -1e99;
        double min = 1e99;
        unsigned long offset = (layer*w*h);
        for (int i = 0; i < w; ++i)
        for (int j = 0; j < h; ++j)
        {
            at (i,j) = buffer[((h - j - 1) * w + i) + offset];
            max = fmax (max, at (i,j));
            min = fmin (min, at (i,j));
        }

        std::cerr << "info: loaded " << w << "x" << h
            << " FITS file, max lumin = " << max << ", min = "
            << min << "\n";
    }
    
    FitsFile (const std::vector<std::vector<double>> &map, const std::vector<std::vector<string>> &WCSdata) //Create a FitsFile from 2D-vector containing image
    {
        std::cout << "Converting map to FitsFile...\n";
        set_coordinates (0, 0, map.size(), map.at(0).size(), map.size(), map.at(0).size()); //Set size
        
        WCSkeynames_here = WCSdata.at(0); //Add Header data
        WCSvalues = WCSdata.at(1);
        
        for(uint i=0; i<map.size(); i++)
        for(uint j=0; j<map.at(0).size(); j++)
        {
            at(i,j) = map.at(i).at(j);
        }
        std::cout << "Converted map to FitsFile!" << std::endl;
    }
        
    // allow overwriting pixels by the user
    double &at (int i, int j)
    {
        return Photo::at (i, j);
    }
    
    // return WCS data in given order
    std::vector<std::vector<string>> returnWCSdata()
    {
        std::vector<std::vector<string>> ret;
        ret.push_back(WCSkeynames_here);
        ret.push_back(WCSvalues);
        return ret;
    }

    void addnoise ()
    {
        at (2,2) = 1;
    }
    
    string giveKeyvalue(string name)
    {
        string content = "0";
        auto it = std::find(WCSkeynames_here.begin(), WCSkeynames_here.end(), name); 
        if (it != WCSkeynames_here.end())
        {
            int pos1 = distance(WCSkeynames_here.begin(), it);
            content = WCSvalues.at(pos1);
            //std::cout << WCSkeynames_here.at(pos1) << " " << content << std::endl;
            
        }
        return content;
    }
    
    void setKeyValue(string name, string newValue)
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
        std::vector<std::vector<string>> WCS_there = b.returnWCSdata();
        
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
        {
//             //std::cout << i << " " << j << std::endl;
            at(i,j) = std::abs(at(i,j));
        }
    }
    
//    void upscale(int factor)
//    {
//        FitsFile thissave = *this;
//        set_coordinates(0,0, factor*thissave.width()*lenx, factor*thissave.height()*leny, factor*thissave.width(), factor*thissave.height());
//         // Go over each line factor times and repeat original value factor times
//        for(int linenr=0; linenr<thissave.width(); linenr++)
//        for(int i1=0; i1<factor; i1++) //1st factor
//        for(int colnr=0; colnr<thissave.height(); colnr++)
//        for(int i2=0; i2<factor; i2++) //2nd factor
//        {
//            //std::cout << colnr << " " << i2 << " " << factor*colnr+i2 << std::endl;
//            at(factor*linenr+i1, factor*colnr+i2) = thissave(linenr,colnr);
//            //std::cout << width() << " " << thissave.width() << std::endl;
//        }
//        
//        //Rescale Coordinate system
//        setKeyValue("CDELT1",std::to_string( std::stod(giveKeyvalue("CDELT1"))/factor ) );
//        setKeyValue("CDELT2",std::to_string( std::stod(giveKeyvalue("CDELT2"))/factor ) );
//        setKeyValue("CRPIX1",std::to_string( factor*(std::stod(giveKeyvalue("CRPIX1"))-1) ) );
//        setKeyValue("CRPIX2",std::to_string( factor*(std::stod(giveKeyvalue("CRPIX2"))-0.5) ) );
//        WCSkeynames_here.push_back("COMMENT");
//        WCSvalues.push_back("File was scaled by factor "+std::to_string(factor));
//    }
};


template <typename PHOTO> //Writes one fits file from photo. Possible to write absolute value, argument, flipped image in any axis
void writeImage(PHOTO minkmap, string filename,const std::vector<std::vector<string>>& WCSdata, bool absolute = true, bool arg = false, bool flipY = true, bool flipX=false, bool highPrec = false)
{
    std::cout << "Writing FITS file " + filename +  " ..." << std::endl;
    long naxis = 2;
    long naxes[2] = { minkmap.width(), minkmap.height() };
    std::unique_ptr<CCfits::FITS> pFits;
    if(!highPrec) pFits.reset(new CCfits::FITS("!" + filename + ".fits" , FLOAT_IMG , naxis , naxes) );
    else pFits.reset(new CCfits::FITS("!" + filename + ".fits" , DOUBLE_IMG , naxis , naxes) );
    
    
    // NOTE: At this point we assume that there is only 1 layer.
    long nelements = std::accumulate(& naxes[0], & naxes[naxis], 1, std::multiplies<long>());
    std::valarray<double> array(nelements);
    
    for(long i=0; i<nelements; i++)
    {
        complex_t val;
        if (flipY && !flipX) val = minkmap(i%minkmap.width() , minkmap.height() - i/minkmap.width()-1);
        else if (!flipY && flipX) val = minkmap(minkmap.width() - i%minkmap.width()-1 , i/minkmap.width());
        else if (flipX && flipY) val = minkmap(minkmap.width() - i%minkmap.width()-1 ,minkmap.height() - i/minkmap.width()-1);
        else val = minkmap(i%minkmap.width() , i/minkmap.width());
        
        if(absolute && !arg) array[i] = 1.*std::abs(val);
        else if (!arg) array[i] = 1.*std::abs(val)*cos(1.*std::arg(val));
        else array[i] = (1.*std::arg(val) + 3.14159)*(std::abs(val)>0.01);
    }
    long fpixel(1);
    
    //Header information
    for(uint i=0; i<WCSdata.at(0).size(); i++)
    {
        //std::cout << i << " " << WCSdata.at(0).at(i) << " " << WCSdata.at(1).at(i) << std::endl;
        if (WCSdata.at(0).at(i) != "COMMENT")
        {
            pFits->pHDU().addKey(WCSdata.at(0).at(i),WCSdata.at(1).at(i)," ");
        }
        else
        {
            pFits->pHDU().writeComment(WCSdata.at(1).at(i));
        }
    }
    
    pFits->pHDU().write(fpixel, nelements, array);
    
    std::cout << "Done! \n";
}

template <typename PHOTO> //Writes one fits file from vector of photos. Possible to write absolute value, argument, flipped image in any axis
void write3Dimage(const std::vector<PHOTO>& minkmaps, string filename,const std::vector<std::vector<string>>& WCSdata, bool absolute = true, bool arg = false, bool flipY = true, bool flipX=false)
{
    std::cout << "Writing 3D FITS file " + filename + " ..." << std::endl;
    
    
    int N = minkmaps.size();
    
    int w = minkmaps.at(0).width();
    int h = minkmaps.at(0).height();
    
    long naxis = 3;
    long naxes[3] = { w, h, N};
    std::unique_ptr<CCfits::FITS> pFits;
    pFits.reset(new CCfits::FITS("!" + filename + ".fits" , FLOAT_IMG , naxis , naxes) );
    
    
    long nelements = std::accumulate(& naxes[0], & naxes[naxis], 1, std::multiplies<long>());
    std::valarray<double> array(nelements);
    
    for(long i=0; i<nelements/N; i++)
    for(int j=0; j<N; j++)
    {
        complex_t val;
        if (flipY && !flipX) val = minkmaps.at(j)(i%w , h - i/w-1);
        else if (!flipY && flipX) val = minkmaps.at(j)(w - i%w -1, i/w);
        else if (flipX && flipY) val = minkmaps.at(j)(w - i%w-1 , h - i/w-1);
        else val = minkmaps.at(j)(i%w , i/w);
        
        //array[i+j*nelements/N]   = 1.*std::abs(minkmaps.at(j)(i%w , h - i/w));
        if(absolute && !arg) array[i+j*nelements/N] = 1.*std::abs(minkmaps.at(j)(i%w , h - i/w));
        else if (!arg) array[i+j*nelements/N] = 1.*std::abs(minkmaps.at(j)(i%w , h - i/w))*cos(1.*std::arg(minkmaps.at(j)(i%w , h - i/w)));
        else array[i+j*nelements/N] = 1.*std::arg(minkmaps.at(j)(i%w , h - i/w));
    }
    long fpixel(1);
    
    
    for(uint i=0; i<WCSdata.at(0).size(); i++)
    {
        if (WCSdata.at(0).at(i) != "COMMENT")
        {
            pFits->pHDU().addKey(WCSdata.at(0).at(i),WCSdata.at(1).at(i)," ");
        }
        else 
        {
            pFits->pHDU().writeComment(WCSdata.at(1).at(i));
        }
    }
    
    pFits->pHDU().write(fpixel, nelements, array);
    
    std::cout << "Done! \n";
}

void erosionDilation(FitsFile& infile, int diameter)
{
    FitsFile original = infile;
    double rsquared = diameter*diameter/4;
    int begin1 = 1,begin2 = 1;
    int end1 = infile.height();
    int end2 = infile.width();
    
    
    //Erosion
    for (int j = begin1; j < end1; j++)
    for (int i = begin2; i < end2; i++)
    {
        double local_minimum = 99999.;
        
        for(int k=-(diameter/2); k<=(diameter/2); k++)
        for(int l=-(diameter/2); l<=(diameter/2); l++)
        {
            double dist = pow(k,2)+pow(l,2);
            if( dist < rsquared )
            {
                double YInImage = std::min(j+l, end1);
                double XInImage = std::min(i+k, end2);
                if(original(XInImage,YInImage)<local_minimum) local_minimum = original(XInImage,YInImage);
            }
        }
        
        infile(i,j) = local_minimum;
    }
    
    //Dilation
    original = infile; //New original is eroded image
    
    for (int j = begin1; j < end1; j++)
    for (int i = begin2; i < end2; i++)
    {
        double local_maximum = 0.;
        
        for(int k=-(diameter/2); k<=(diameter/2); k++)
        for(int l=-(diameter/2); l<=(diameter/2); l++)
        {
            double dist = pow(k,2)+pow(l,2);
            if( dist < rsquared )
            {
                double YInImage = std::min(j+l, end1);
                double XInImage = std::min(i+k, end2);
                if(original(XInImage,YInImage)>local_maximum) local_maximum = original(XInImage,YInImage);
            }
        }
        
        infile(i,j) = local_maximum;
    }
}

double sumOverSquares(FitsFile& image, int smooth)
{ // Adding over all squares (smoothed) in image
    int stepsize = std::max(smooth/6,1);
    
    double sum = 0.;
    for(int i=std::max(stepsize/2,1); i<image.height(); i+=stepsize)
    for(int j=std::max(stepsize/2,1); j<image.width(); j+=stepsize)
    {
        sum += image(i,j);
    }
    return sum;
}

void takeLogPositive(FitsFile& image)
{
    double min=0.;
    for(int i=0;i<image.height();i++)
    for(int j=0;j<image.width();j++)
    if(image(i,j))
    {
        image(i,j) = std::log(image(i,j));
        min = std::min(min,image(i,j));
    }
    for(int i=0;i<image.height();i++)
    for(int j=0;j<image.width();j++)
    if(image(i,j))
    {
        image(i,j) -= min;
    }
}

void binImage(FitsFile& infile, int factor=2)
{
    std::cout << "Binning...\n";
    FitsFile newFile = infile;
    newFile.set_coordinates(0, 0, (int)(infile.width()/factor), (int)(infile.height()/factor), (int)(infile.width()/factor), (int)(infile.height()/factor));
    
    for(int i=1;i<newFile.width();i++)
    for(int j=1;j<newFile.height();j++)
    {
        double sum = 0;
        for(int step1=0;step1<factor;step1++)
        for(int step2=0;step2<factor;step2++)
        {
            sum+=infile(std::min(infile.width(),factor*i+step1), std::min(infile.height(),factor*j+step2) );
        }
        sum/=factor*factor;
        newFile(i,j) = sum;
    }
    infile = newFile;
    string newCRPIX1 = std::to_string(std::stod(infile.giveKeyvalue("CRPIX1"))/factor);
    string newCRPIX2 = std::to_string(std::stod(infile.giveKeyvalue("CRPIX2"))/factor);
    string newCDELT1 = std::to_string(std::stod(infile.giveKeyvalue("CDELT1"))*factor);
    string newCDELT2 = std::to_string(std::stod(infile.giveKeyvalue("CDELT2"))*factor);
    infile.setKeyValue("CRPIX1", newCRPIX1);
    infile.setKeyValue("CRPIX2", newCRPIX2);
    infile.setKeyValue("CDELT1", newCDELT1);
    infile.setKeyValue("CDELT2", newCDELT2);
    
    std::cout << "Binned!\n";
}




#endif
