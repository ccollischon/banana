#include "fitsfile.hpp"

#include <CCfits/CCfits>

using namespace papaya2;

FitsFile::FitsFile(const string &infilename, int layer /*= 0*/, unsigned hdu /*= 0*/)
{
    if (hdu != 0)
        throw std::runtime_error("hdu != 0 not implemented");
    CCfits::FITS::setVerboseMode(true);
    CCfits::FITS infile(infilename.c_str(), CCfits::Read, true);
    CCfits::PHDU &phdu = infile.pHDU();

    phdu.readAllKeys();
    set_coordinates(0, 0, phdu.axis(0), phdu.axis(1), phdu.axis(0),
            phdu.axis(1));
    //set_coordinates (0, 0, 1, 1, phdu.axis (0), phdu.axis (1));
    std::valarray<double> buffer;
    phdu.read(buffer);
    const int w = width();
    const int h = height();
    double max = -1e99;
    double min = 1e99;
    unsigned long offset = (layer * w * h);
    for (int i = 0; i < w; ++i)
        for (int j = 0; j < h; ++j)
        {
            at(i, j) = buffer[((h - j - 1) * w + i) + offset];
            max = fmax(max, at(i, j));
            min = fmin(min, at(i, j));
        }

    std::cerr << "info: loaded " << w << "x" << h << " FITS file, max lumin = "
            << max << ", min = " << min << "\n";
}

FitsFile::FitsFile(const string &infilename, std::vector<string>& WCSkeynames,
        int layer /*= 0*/, unsigned hdu /*= 0*/)
{
    if (hdu != 0)
        throw std::runtime_error("hdu != 0 not implemented");
    CCfits::FITS::setVerboseMode(true);
    CCfits::FITS infile(infilename.c_str(), CCfits::Read, true);
    CCfits::PHDU &phdu = infile.pHDU();

    // FIXME extract proper sky coordinates from FITS and insert pixel dimensions here.
    // for now, just map everything on the unit square

    std::cout << "Loading WCS data..." << std::endl;
    WCSkeynames_here = WCSkeynames;
    for (auto keyword : WCSkeynames)
    {
        string keyval = "";
        try
        {
            phdu.readKey(keyword, keyval);
        } catch (const CCfits::HDU::NoSuchKeyword&)
        {
            if(keyword=="RADECSYS"){ //RADECSYS and RADESYS should be equivalent
                try {
                    phdu.readKey("RADESYS", keyval);
                } catch (const CCfits::HDU::NoSuchKeyword&)
                {
                    std::cerr<<"Warning: neither RADESYS nor RADECSYS could be determined\n";
                    keyval = "0.";
                }
            }
            else{
                keyval = "0.";
            }
        }
        WCSvalues.push_back(keyval);
    }
    //phdu.readKeys(WCSkeynames,WCSvalues);

    set_coordinates(0, 0, phdu.axis(0), phdu.axis(1), phdu.axis(0),
            phdu.axis(1));
    //set_coordinates (0, 0, 1, 1, phdu.axis (0), phdu.axis (1));
    std::valarray<double> buffer;
    phdu.read(buffer);
    const int w = width();
    const int h = height();
    double max = -1e99;
    double min = 1e99;
    unsigned long offset = (layer * w * h);
    for (int i = 0; i < w; ++i)
        for (int j = 0; j < h; ++j)
        {
            at(i, j) = buffer[((h - j - 1) * w + i) + offset];
            max = fmax(max, at(i, j));
            min = fmin(min, at(i, j));
        }

    std::cerr << "info: loaded " << w << "x" << h << " FITS file, max lumin = "
            << max << ", min = " << min << "\n";
}

FitsFile::FitsFile(const std::vector<std::vector<double>> &map,
        const std::vector<std::vector<string>> &WCSdata) //Create a FitsFile from 2D-vector containing image
{
    std::cout << "Converting map to FitsFile...\n";
    set_coordinates(0, 0, map.size(), map.at(0).size(), map.size(),
            map.at(0).size()); //Set size

    WCSkeynames_here = WCSdata.at(0); //Add Header data
    WCSvalues = WCSdata.at(1);

    for (uint i = 0; i < map.size(); i++)
        for (uint j = 0; j < map.at(0).size(); j++)
        {
            at(i, j) = map.at(i).at(j);
        }
    std::cout << "Converted map to FitsFile!" << std::endl;
}

std::string FitsFile::giveKeyvalue(std::string name)
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

void FitsFile::setKeyValue(std::string name, std::string newValue)
{
	auto it = std::find(WCSkeynames_here.begin(), WCSkeynames_here.end(), name);
	if (it != WCSkeynames_here.end())
	{
		int pos1 = distance(WCSkeynames_here.begin(), it);
		WCSvalues.at(pos1) = newValue;
		//std::cout << WCSkeynames_here.at(pos1) << " " << content << std::endl;
	}
}

FitsFile &FitsFile::operator-=(FitsFile& b) //subtract two files with equal CDELT at correct sky positions
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


template<typename PHOTO> //Writes one fits file from photo. Possible to write absolute value, argument, flipped image in any axis
void writeImage(const PHOTO& minkmap, string filename,
        const std::vector<std::vector<string>>& WCSdata, bool absolute /*= true*/,
        bool arg /*= false*/, bool flipY /*= true*/, bool flipX /*= false*/,
        bool highPrec /*= false*/)
{
    std::cout << "Writing FITS file " + filename + " ..." << std::endl;
    long naxis = 2;
    long naxes[2] =
    { minkmap.width(), minkmap.height() };
    std::unique_ptr<CCfits::FITS> pFits;
    if (!highPrec)
        pFits.reset(
                new CCfits::FITS("!" + filename + ".fits", FLOAT_IMG, naxis,
                        naxes));
    else
        pFits.reset(
                new CCfits::FITS("!" + filename + ".fits", DOUBLE_IMG, naxis,
                        naxes));

    // NOTE: At this point we assume that there is only 1 layer.
    long nelements = std::accumulate(&naxes[0], &naxes[naxis], 1,
            std::multiplies<long>());
    std::valarray<double> array(nelements);

    for (long i = 0; i < nelements; i++)
    {
        complex_t val;
        if (flipY && !flipX)
            val = minkmap(i % minkmap.width(),
                    minkmap.height() - i / minkmap.width() - 1);
        else if (!flipY && flipX)
            val = minkmap(minkmap.width() - i % minkmap.width() - 1,
                    i / minkmap.width());
        else if (flipX && flipY)
            val = minkmap(minkmap.width() - i % minkmap.width() - 1,
                    minkmap.height() - i / minkmap.width() - 1);
        else
            val = minkmap(i % minkmap.width(), i / minkmap.width());

        if (absolute && !arg)
            array[i] = 1. * std::abs(val);
        else if (!arg)
            array[i] = 1. * std::abs(val) * cos(1. * std::arg(val));
        else
            array[i] = (1. * std::arg(val) + 3.14159) * (std::abs(val) > 0.01);
    }
    long fpixel(1);

    //Header information
    for (uint i = 0; i < WCSdata.at(0).size(); i++)
    {
        //std::cout << i << " " << WCSdata.at(0).at(i) << " " << WCSdata.at(1).at(i) << std::endl;
        if (WCSdata.at(0).at(i) != "COMMENT")
        {
            //see if field is a number, and convert if it is
            if(WCSdata.at(1).at(i).find_first_not_of("-0123456789") == string::npos) //int
            {
                pFits->pHDU().addKey(WCSdata.at(0).at(i), std::stoi(WCSdata.at(1).at(i)), " ");
            } else if (WCSdata.at(1).at(i).find_first_not_of("-0123456789.") == string::npos) //double
            {
                pFits->pHDU().addKey(WCSdata.at(0).at(i), std::stod(WCSdata.at(1).at(i)), " ");
            } else {
                pFits->pHDU().addKey(WCSdata.at(0).at(i), WCSdata.at(1).at(i), " ");
            }
        }
        else
        {
            pFits->pHDU().writeComment(WCSdata.at(1).at(i));
        }
    }

    pFits->pHDU().write(fpixel, nelements, array);

    std::cout << "Done! \n";
}

template // explicit instanciation
void writeImage(const FitsFile& minkmap, string filename,
        const std::vector<std::vector<string>>& WCSdata, bool absolute /*= true*/,
        bool arg /*= false*/, bool flipY /*= true*/, bool flipX /*= false*/,
        bool highPrec /*= false*/);

template // explicit instanciation
void writeImage(const papaya2::BasicPhoto<std::complex<double> >& minkmap, string filename,
        const std::vector<std::vector<string>>& WCSdata, bool absolute /*= true*/,
        bool arg /*= false*/, bool flipY /*= true*/, bool flipX /*= false*/,
        bool highPrec /*= false*/);

template // explicit instanciation
void writeImage(const papaya2::BasicPhoto<double>& minkmap, string filename,
        const std::vector<std::vector<string>>& WCSdata, bool absolute /*= true*/,
        bool arg /*= false*/, bool flipY /*= true*/, bool flipX /*= false*/,
        bool highPrec /*= false*/);

template<typename PHOTO> //Writes one fits file from vector of photos. Possible to write absolute value, argument, flipped image in any axis
void write3Dimage(const std::vector<PHOTO>& minkmaps, string filename,
        const std::vector<std::vector<string>>& WCSdata,
        bool absolute /*= true*/, bool arg /*= false*/, bool flipY /*= true*/,
        bool flipX /*= false*/)
{
    std::cout << "Writing 3D FITS file " + filename + " ..." << std::endl;

    int N = minkmaps.size();

    int w = minkmaps.at(0).width();
    int h = minkmaps.at(0).height();

    long naxis = 3;
    long naxes[3] =
    { w, h, N };
    std::unique_ptr<CCfits::FITS> pFits;
    pFits.reset(
            new CCfits::FITS("!" + filename + ".fits", FLOAT_IMG, naxis,
                    naxes));

    long nelements = std::accumulate(&naxes[0], &naxes[naxis], 1,
            std::multiplies<long>());
    std::valarray<double> array(nelements);

    for (long i = 0; i < nelements / N; i++)
        for (int j = 0; j < N; j++)
        {
            complex_t val;
            if (flipY && !flipX)
                val = minkmaps.at(j)(i % w, h - i / w - 1);
            else if (!flipY && flipX)
                val = minkmaps.at(j)(w - i % w - 1, i / w);
            else if (flipX && flipY)
                val = minkmaps.at(j)(w - i % w - 1, h - i / w - 1);
            else
                val = minkmaps.at(j)(i % w, i / w);

            //array[i+j*nelements/N]   = 1.*std::abs(minkmaps.at(j)(i%w , h - i/w));
            if (absolute && !arg)
                array[i + j * nelements / N] = 1.
                        * std::abs(minkmaps.at(j)(i % w, h - i / w));
            else if (!arg)
                array[i + j * nelements / N] = 1.
                        * std::abs(minkmaps.at(j)(i % w, h - i / w))
                        * cos(1. * std::arg(minkmaps.at(j)(i % w, h - i / w)));
            else
                array[i + j * nelements / N] = 1.
                        * std::arg(minkmaps.at(j)(i % w, h - i / w));
        }
    long fpixel(1);

    for (uint i = 0; i < WCSdata.at(0).size(); i++)
    {
        if (WCSdata.at(0).at(i) != "COMMENT")
        {
            pFits->pHDU().addKey(WCSdata.at(0).at(i), WCSdata.at(1).at(i), " ");
        }
        else
        {
            pFits->pHDU().writeComment(WCSdata.at(1).at(i));
        }
    }

    pFits->pHDU().write(fpixel, nelements, array);

    std::cout << "Done! \n";
}


template // explicit instanciation
void write3Dimage(const std::vector<papaya2::BasicPhoto<std::complex<double> >>& minkmaps, string filename,
        const std::vector<std::vector<string>>& WCSdata,
        bool absolute /*= true*/, bool arg /*= false*/, bool flipY /*= true*/,
        bool flipX /*= false*/);

void erosionDilation(FitsFile& infile, int diameter)
{
    FitsFile original = infile;
    double rsquared = diameter * diameter / 4;
    int begin1 = 1, begin2 = 1;
    int end1 = infile.height();
    int end2 = infile.width();

    //Erosion
    for (int j = begin1; j < end1; j++)
        for (int i = begin2; i < end2; i++)
        {
            double local_minimum = 99999.;

            for (int k = -(diameter / 2); k <= (diameter / 2); k++)
                for (int l = -(diameter / 2); l <= (diameter / 2); l++)
                {
                    double dist = pow(k, 2) + pow(l, 2);
                    if (dist < rsquared)
                    {
                        double YInImage = std::min(j + l, end1);
                        double XInImage = std::min(i + k, end2);
                        if (original(XInImage, YInImage) < local_minimum)
                            local_minimum = original(XInImage, YInImage);
                    }
                }

            infile(i, j) = local_minimum;
        }

    //Dilation
    original = infile; //New original is eroded image

    for (int j = begin1; j < end1; j++)
        for (int i = begin2; i < end2; i++)
        {
            double local_maximum = 0.;

            for (int k = -(diameter / 2); k <= (diameter / 2); k++)
                for (int l = -(diameter / 2); l <= (diameter / 2); l++)
                {
                    double dist = pow(k, 2) + pow(l, 2);
                    if (dist < rsquared)
                    {
                        double YInImage = std::min(j + l, end1);
                        double XInImage = std::min(i + k, end2);
                        if (original(XInImage, YInImage) > local_maximum)
                            local_maximum = original(XInImage, YInImage);
                    }
                }

            infile(i, j) = local_maximum;
        }
}

double sumOverSquares(FitsFile& image, int smooth)
{ // Adding over all squares (smoothed) in image
    int stepsize = std::max(smooth / 6, 1);

    double sum = 0.;
    for (int i = std::max(stepsize / 2, 1); i < image.height(); i += stepsize)
        for (int j = std::max(stepsize / 2, 1); j < image.width(); j +=
                stepsize)
        {
            sum += image(i, j);
        }
    return sum;
}

void takeLogPositive(FitsFile& image)
{
    double min = 0.;
    for (int i = 0; i < image.height(); i++)
        for (int j = 0; j < image.width(); j++)
            if (image(i, j))
            {
                image(i, j) = std::log(image(i, j));
                min = std::min(min, image(i, j));
            }
    for (int i = 0; i < image.height(); i++)
        for (int j = 0; j < image.width(); j++)
            if (image(i, j))
            {
                image(i, j) -= min;
            }
}

void binImage(FitsFile& infile, int factor = 2)
{
    std::cout << "Binning...\n";
    FitsFile newFile = infile;
    newFile.set_coordinates(0, 0, (int) (infile.width() / factor),
            (int) (infile.height() / factor), (int) (infile.width() / factor),
            (int) (infile.height() / factor));

    for (int i = 1; i < newFile.width(); i++)
        for (int j = 1; j < newFile.height(); j++)
        {
            double sum = 0;
            for (int step1 = 0; step1 < factor; step1++)
                for (int step2 = 0; step2 < factor; step2++)
                {
                    sum += infile(std::min(infile.width(), factor * i + step1),
                            std::min(infile.height(), factor * j + step2));
                }
            sum /= factor * factor;
            newFile(i, j) = sum;
        }
    infile = newFile;
    string newCRPIX1 = std::to_string(
            std::stod(infile.giveKeyvalue("CRPIX1")) / factor);
    string newCRPIX2 = std::to_string(
            std::stod(infile.giveKeyvalue("CRPIX2")) / factor);
    string newCDELT1 = std::to_string(
            std::stod(infile.giveKeyvalue("CDELT1")) * factor);
    string newCDELT2 = std::to_string(
            std::stod(infile.giveKeyvalue("CDELT2")) * factor);
    infile.setKeyValue("CRPIX1", newCRPIX1);
    infile.setKeyValue("CRPIX2", newCRPIX2);
    infile.setKeyValue("CDELT1", newCDELT1);
    infile.setKeyValue("CDELT2", newCDELT2);

    std::cout << "Binned!\n";
}

