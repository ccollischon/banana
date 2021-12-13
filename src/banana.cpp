// banana: trivial program which load a FITS file, and computes
// IMT's for a number of interpolated marching squares thresholds
// using the preliminary Papaya2 library.
// 2018 Sebastian Kapfer <sebastian.kapfer@fau.de>
// 2018 Jenny Wagner <j.wagner@uni-heidelberg.de>
//
// 2020: no longer trivial, large extension to functionality
// Minkmaps, line density maps, bubble detection, pointspread, histograms
// Caroline Collischon <caroline.collischon@fau.de>

#include <iostream>
#include <stdlib.h>
#include <string>
struct {
    std::string objectDIR;
    std::string pngDIR;
    std::string resultDIR;
    std::string linedensDIR;
} settings;


#include "src/papaya2.hpp"
#include "src/readarg.hpp"
#include "src/photofunctions.hpp"
#include "util/fitsfile.hpp"
#include "src/readdata.hpp"
#include "util/savepng.hpp"
#include "src/originalBanana.hpp"
#include "src/hedgehog.hpp"
#include "src/peach.hpp"
#include "src/sugarsnappea.hpp"
#include <valarray>
#include <sstream>

#include <CCfits/CCfits>

using namespace papaya2;

//global
double crpix_source;

void convertRegions(string boxname, string namename)
{
    string prefix = "DEM_ext";
    std::ofstream outfile_names("DEM_ext_objects");
    std::vector<std::vector<double>> boxes = readTable(boxname);
    std::vector<string> names = readLines(namename);
    for(uint i=0; i<boxes.size();i++)
    {
        string filename;
        int fileA = std::stoi(names.at(i).substr(7,3));
        char fileB = names.at(i).at(11);
        if (fileB ==' ')
            fileB = '_';
        filename = settings.objectDIR + prefix + std::to_string(fileA) + "_" + fileB;
        //filename = "./objects/"+ prefix + std::to_string(i);
        Datafile outfile(filename);
        
        outfile << boxes.at(i).at(0) << boxes.at(i).at(1) << boxes.at(i).at(2) << boxes.at(i).at(3) << std::endl;
        outfile_names << prefix << std::to_string(fileA) + "_" + fileB << std::endl;
    }
    outfile_names.close();
}

void imageToPointPattern(FitsFile& image, int smooth, std::string name, int ENpoints=5000, double thresh=0.)
{ // Converts a (smoothed) image to a R-readable point pattern by placing n points at every square where n corresponds to the (int)pixel value of that point
    //ENpoints = expected number of points in final image
    //double CRPIX1_here = std::stod(image.giveKeyvalue("CRPIX1")); // I tried using the theoretical shift, but it didn't work
    //double shift = crpix_source-CRPIX1_here;
    int stepsize = std::max(smooth/6,1);
    
    std::cout << "Converting image to pattern, saving in file " << name << std::endl;
    
    double p = ENpoints/sumOverSquares(image,smooth); //Probability to keep point    
    //std::cout << "2*shift: " << 2*shift << " stepsize: " << stepsize << std::endl;
    std::ofstream ofs(name);
    for(int i=std::max(stepsize/2,1); i<image.height(); i+=stepsize)
    for(int j=std::max(stepsize/2,1); j<image.width(); j+=stepsize)
    {
        double valueHere = image(i,j);
        if(valueHere > thresh)
        {
            //std::cout << valueHere << std::endl;
            std::vector<int> unshiftedCoords = image.giveUnshiftedds9Coord(i,j);
            for(int density=0;density<valueHere;density++)
            {
                double addToRa = randInInterval(-stepsize/2, stepsize/2);
                double addToDec = randInInterval(-stepsize/2, stepsize/2);
                double randHere = randInInterval(0.,1.); // Only write number if selected
                if(randHere<p)  ofs << unshiftedCoords.at(0)+addToRa << "\t" << unshiftedCoords.at(1)+addToDec << std::endl;
            }
        }
    }
    ofs.close();
}





void makeMinkmap(std::string infilename, FitsFile& infile, std::string outminmap, int s, bool threeD, int squaresize, bool average, int smooth, bool absolute_avg, double min_thresh, double max_thresh, int num_thresh, bool arg)
{ // make Minkowski map for given parameters
//Get Coordinate (or other) info from infile
    std::vector<std::vector<string>> minkmap_WCSdata = infile.returnWCSdata();

//Collect minkmap-results in these
    complex_image minkmap;
    std::vector<complex_image> minkmaps;

//Shift coordinates by 1/2 pixel each squaresize-1/the right amount for smooth + dividing into smooth/6 pixels
    
    if(infile.giveKeyvalue("CRPIX1")=="0")
    {
        minkmap_WCSdata.at(0).push_back("CRPIX1");
        minkmap_WCSdata.at(1).push_back(std::to_string(crpix_source - 0.5*(std::max(1.*squaresize,1.*smooth-smooth/6)-1)));
        minkmap_WCSdata.at(0).push_back("CRPIX2");
        minkmap_WCSdata.at(1).push_back(std::to_string(crpix_source - 0.5*(std::max(1.*squaresize,1.*smooth-smooth/6)-1)));
    }
    else
    {
        minkmap_WCSdata.at(1).at(4) = std::to_string(stod(minkmap_WCSdata.at(1).at(4)) - 0.5*(std::max(1.*squaresize,1.*smooth-smooth/6)-1));
        minkmap_WCSdata.at(1).at(5) = std::to_string(stod(minkmap_WCSdata.at(1).at(5)) + 0.5*(std::max(1.*squaresize,1.*smooth-smooth/6)-1));
    }
    minkmap_WCSdata.at(0).push_back("COMMENT");
    minkmap_WCSdata.at(1).push_back("Original file: "+infilename);
    
    minkmap_WCSdata.at(0).push_back("HISTORY");
    minkmap_WCSdata.at(1).push_back("Edited with banana.cpp using papaya2.hpp");


    std::string sstring;
    if(s==1) sstring = "_euler";
    else if (s==4234) sstring = "_area";
    else sstring = "_s="+ std::to_string(s);

    outminmap = outminmap + sstring;
    if(arg) outminmap += "_arg";
    
/***************************************/
//Start actual processing
/***************************************/

    if (threeD && squaresize == 2 && !average)
    {
        minkmap_WCSdata.at(0).push_back("COMMENT");
        if (s==0) minkmap_WCSdata.at(1).push_back("Minkowski map for perimeter, smoothcount "+std::to_string(smooth));
        else if (s==1) minkmap_WCSdata.at(1).push_back("Minkowski map for Euler characteristic, smoothcount "+std::to_string(smooth));
        else if (s==4234) minkmap_WCSdata.at(1).push_back("Minkowski map for area, smoothcount "+std::to_string(smooth));
        else minkmap_WCSdata.at(1).push_back("Minkowski map for s = "+std::to_string(s)+", smoothcount "+std::to_string(smooth));
        int i = 1;
        for (auto thresh : logspace (min_thresh, max_thresh, num_thresh, true))
        {
            complex_image minkmap2;
            minkowski_map_interpolated_marching_squares(&minkmap2, infile, thresh, s);
            
            if(smooth!=0) smooth_map_plus(minkmap2,smooth);
            
            minkmaps.push_back(minkmap2);
            minkmap_WCSdata.at(0).push_back("COMMENT");
            minkmap_WCSdata.at(1).push_back("Threshold for layer "+std::to_string(i)+": "+std::to_string(thresh));
            std::cout << "Threshold for layer "+std::to_string(i)+": "+std::to_string(thresh) << std::endl;
            i++;
        }
        char buffer1 [40];
        int n = sprintf(buffer1,"%g_%g_%d",min_thresh,max_thresh,num_thresh);
        n++;
        if (smooth==0) write3Dimage(minkmaps, outminmap+"_3d"+"_thresh_"+buffer1, minkmap_WCSdata, s, arg);
        else write3Dimage(minkmaps, outminmap+"_smooth="+std::to_string(smooth)+"_3d"+"_thresh_"+buffer1, minkmap_WCSdata, s, arg);
    }
    else if (threeD && squaresize != 2 && smooth == 0 && !average)
    {
        minkmap_WCSdata.at(0).push_back("COMMENT");
        if (s==0) minkmap_WCSdata.at(1).push_back("Minkowski map for perimeter, squaresize "+std::to_string(squaresize));
        else if (s==1) minkmap_WCSdata.at(1).push_back("Minkowski map for Euler characteristic, squaresize "+std::to_string(squaresize));
        else if (s==4234) minkmap_WCSdata.at(1).push_back("Minkowski map for area, smoothcount "+std::to_string(squaresize));
        else minkmap_WCSdata.at(1).push_back("Minkowski map for s = "+std::to_string(s)+", squaresize "+std::to_string(squaresize));
        int i = 1;
        for (auto thresh : logspace (min_thresh, max_thresh, num_thresh, true))
        {
            complex_image minkmap2;
            minkowski_map_interpolated_marching_bigsquares(&minkmap2, infile, thresh, s, squaresize);
            minkmaps.push_back(minkmap2);
            minkmap_WCSdata.at(0).push_back("COMMENT");
            minkmap_WCSdata.at(1).push_back("Threshold for layer "+std::to_string(i)+": "+std::to_string(thresh));
            std::cout << "Threshold for layer "+std::to_string(i)+": "+std::to_string(thresh) << std::endl;
            i++;
        }
        char buffer1 [40];
        int n = sprintf(buffer1,"%g_%g_%d",min_thresh,max_thresh,num_thresh);
        n++;
        
        complex_image averagemap;
        if(average)
        {
            averagemap = averagemaps(minkmaps);
            outminmap = outminmap+"_avg";
            writeImage(averagemap, outminmap+"_squaresize="+std::to_string(smooth)+"_thresh_"+buffer1, minkmap_WCSdata, s, arg);
        }
        else write3Dimage(minkmaps, outminmap+"_squaresize="+std::to_string(squaresize)+"_3d"+"_thresh_"+buffer1, minkmap_WCSdata, s, arg);
    }
    else if (!threeD && squaresize != 2)
    {
        
        minkowski_map_interpolated_marching_bigcircles(&minkmap, infile, max_thresh, s, squaresize);
        
        minkmap_WCSdata.at(0).push_back("COMMENT");
        if (s==0) minkmap_WCSdata.at(1).push_back("Minkowski map for perimeter, squaresize "+std::to_string(squaresize));
        else if (s==1) minkmap_WCSdata.at(1).push_back("Minkowski map for Euler characteristic, squaresize "+std::to_string(squaresize));
        else if (s==4234) minkmap_WCSdata.at(1).push_back("Minkowski map for area, smoothcount "+std::to_string(squaresize));
        else minkmap_WCSdata.at(1).push_back("Minkowski map for s = "+std::to_string(s)+", squaresize "+std::to_string(squaresize));
        minkmap_WCSdata.at(0).push_back("COMMENT");
        minkmap_WCSdata.at(1).push_back("Threshold: "+std::to_string(max_thresh));
        char buffer1 [15];
        int n = sprintf(buffer1,"%g",max_thresh);
        n++;
        writeImage(minkmap, outminmap+"_squaresize="+std::to_string(squaresize)+"_thresh_"+buffer1, minkmap_WCSdata, s, arg);
    }
    else if (threeD && average)
    {
        minkmap_WCSdata.at(0).push_back("COMMENT");
        if (s==0) minkmap_WCSdata.at(1).push_back("Averaged Minkowski map for perimeter, smoothcount "+std::to_string(smooth));
        else if (s==1) minkmap_WCSdata.at(1).push_back("Minkowski map for Euler characteristic, smoothcount "+std::to_string(smooth));
        else if (s==4234) minkmap_WCSdata.at(1).push_back("Minkowski map for area, smoothcount "+std::to_string(smooth));
        else minkmap_WCSdata.at(1).push_back("Averaged Minkowski map for s = "+std::to_string(s)+", smoothcount "+std::to_string(smooth));

        complex_image averagemap;
        if(absolute_avg)
        {
            int i = 1;
            for (auto thresh : logspace (min_thresh, max_thresh, num_thresh, true))
            {
                complex_image minkmap2;
                minkowski_map_interpolated_marching_squares(&minkmap2, infile, thresh, s);
                
                if(smooth!=0)
                {
                    std::cout << "Smoothing layer " << i << " ...\n";
                    smooth_map_plus(minkmap2,smooth);
                }
                i==1 ? averagemap=minkmap2.abs() : averagemap+=minkmap2.abs();
                minkmap_WCSdata.at(0).push_back("COMMENT");
                minkmap_WCSdata.at(1).push_back("Threshold for layer "+std::to_string(i)+": "+std::to_string(thresh));
                std::cout << "Threshold for layer "+std::to_string(i)+": "+std::to_string(thresh) << std::endl;
                i++;
            }
            averagemap/=num_thresh;
            outminmap = outminmap+"_absavg";
        }
        else
        {
            int i = 1;
            for (auto thresh : logspace (min_thresh, max_thresh, num_thresh, true))
            {
                complex_image minkmap2;
                minkowski_map_interpolated_marching_squares(&minkmap2, infile, thresh, s);
                
                i==1 ? averagemap=minkmap2 : averagemap+=minkmap2;
                minkmap_WCSdata.at(0).push_back("COMMENT");
                minkmap_WCSdata.at(1).push_back("Threshold for layer "+std::to_string(i)+": "+std::to_string(thresh));
                std::cout << "Threshold for layer "+std::to_string(i)+": "+std::to_string(thresh) << std::endl;
                i++;
            }
            averagemap/=num_thresh;
            if(smooth!=0) smooth_map_plus(averagemap,smooth);
            outminmap = outminmap+"_avg";
        }
        
        char buffer1 [40];
        int n = sprintf(buffer1,"%g_%g_%d",min_thresh,max_thresh,num_thresh);
        n++;
        if (smooth==0) writeImage(averagemap, outminmap+"_thresh_"+buffer1, minkmap_WCSdata, s, arg);
        else writeImage(averagemap, outminmap+"_smooth="+std::to_string(smooth)+"_thresh_"+buffer1, minkmap_WCSdata, s, arg);
    }
    else
    {
        minkowski_map_interpolated_marching_squares(&minkmap, infile, max_thresh, s);
        if(smooth!=0) smooth_map_plus(minkmap,smooth);
        
        minkmap_WCSdata.at(0).push_back("COMMENT");
        if (s==0) minkmap_WCSdata.at(1).push_back("Minkowski map for perimeter, smoothcount "+std::to_string(smooth));
        else if (s==1) minkmap_WCSdata.at(1).push_back("Minkowski map for Euler characteristic, smoothcount "+std::to_string(smooth));
        else if (s==4234) minkmap_WCSdata.at(1).push_back("Minkowski map for area, smoothcount "+std::to_string(smooth));
        else minkmap_WCSdata.at(1).push_back("Minkowski map for s = "+std::to_string(s)+", smoothcount "+std::to_string(smooth));
        minkmap_WCSdata.at(0).push_back("COMMENT");
        minkmap_WCSdata.at(1).push_back("Threshold: "+std::to_string(max_thresh));
        char buffer1 [15];
        int n = sprintf(buffer1,"%g",max_thresh);
        n++;
        if (smooth==0)   writeImage(minkmap, outminmap+"_thresh_"+buffer1, minkmap_WCSdata, s, arg);
        else writeImage(minkmap, outminmap+"_smooth="+std::to_string(smooth)+"_thresh_"+buffer1, minkmap_WCSdata, s, arg);
    }
}



template<typename PHOTO>
void makeRegionImages(PHOTO& infileR, PHOTO& infileG, PHOTO& infileB, string objectname, int erd, bool monochrome, string wavelength)
{ //Create png images of objects given in objectname from photos
    std::vector<std::vector<double>> boxesToInclude = readTable(settings.objectDIR+objectname);
    if(erd) objectname = "erd_" + objectname;
    if(monochrome)
    {
        PHOTO outfile = cutout(infileR, boxesToInclude);
        writeMonoPNG(settings.pngDIR+wavelength+"_"+objectname+".png", outfile);
    }
    else
    {
        PHOTO outfileR = cutout(infileR, boxesToInclude);
        PHOTO outfileG = cutout(infileG, boxesToInclude);
        PHOTO outfileB = cutout(infileB, boxesToInclude);
    
        writeColorPNG(settings.pngDIR+"rgb_"+objectname+".png", outfileR, outfileG, outfileB);
    }
}

void testBanana(string filename)
{
    unsigned width,height;
    string outname = "test_";
    //outname += "smooth3_";
    outname += filename+".dat";
    filename = "./testdata/" + filename;
    outname = "./testdata/" + outname;
    
    std::vector<unsigned char> imageData = decodeOneStep(filename.c_str(),width,height);
    BasicPhoto<double> image;
    PNGtoBasicPhoto(imageData, width, height, image, image, image);
    
    originalBanana("", image, outname, 5, 255, 15);
}



int main (int argc, const char **argv)
{
    string infilename, outfilename, outminmap, histfile, histappendix, outlinedens, outpointspread;
    string maskfilename, combination;
    string boxesToExcludeName = "";
    int smooth = 0, erd = 0;
    bool threeD = false, average = false, absolute_avg = false, make_minkmap = true, makeBanana = false, make_peach = false, make_hist = false, makePNGs = false, arg = false, make_hedgehog = false, make_bubbles = false,make_pointspread = false, make_patternFromImage = false, monochrome = true, combine_regions = false;
    bool filewithoutWCS = false;
    int num_thresh = 9;
    int s = 2;
    double min_thresh = .1;
    double max_thresh = 40.;
    double line_scale = 0.3;
    double line_thresh = 20;
    crpix_source = 0.;
    int squaresize = 2;
    string wavelength = "unspecified", linesorminkmap = "linedensity", greenfile, bluefile;
    std::vector<int> smooths = {10,20,30,40,50,60,100,150,200,250,300,450};
    
    std::vector<string> WCSkeynames = {"CTYPE1","CTYPE2","RADECSYS","EQUINOX","CRPIX1","CRPIX2","CRVAL1","CRVAL2","CDELT1","CDELT2"};

    // process command-line arguments
    for (++argv; *argv; ++argv)
    {
        if (string (*argv) == "in")
        {
            infilename = read_arg <string> (argv++);
        }
        else if (string (*argv) == "mask")
        {
            maskfilename = read_arg <string> (argv++);
        }
        else if (string (*argv) == "combination")
        {
            combination = read_arg <string> (argv++);
        }
        else if (string (*argv) == "stars")
        {
            histfile = read_arg <string> (argv++);
        }
        else if (string (*argv) == "histapp")
        {
            histappendix = read_arg <string> (argv++);
        }
        else if (string (*argv) == "boxesToExclude")
        {
            boxesToExcludeName = read_arg <string> (argv++);
        }
        else if (string (*argv) == "filePrefix")
        {
            wavelength = read_arg <string> (argv++);
        }
        else if (string (*argv) == "mint")
        {
            min_thresh = read_arg <double> (argv++);
        }
        else if (string (*argv) == "maxt")
        {
            max_thresh = read_arg <double> (argv++);
        }
        else if (string (*argv) == "lineThresh")
        {
            line_thresh = read_arg <double> (argv++);
        }
        else if (string (*argv) == "lineScale")
        {
            line_scale = read_arg <double> (argv++);
        }
        else if (string (*argv) == "numt")
        {
            num_thresh = read_arg <unsigned long> (argv++);
        }
        else if (string (*argv) == "outascii")
        {
            outfilename = read_arg <string> (argv++);
            settings.resultDIR = outfilename;
        }
        else if (string (*argv) == "outminmap")
        {
            outminmap = read_arg <string> (argv++);
        }
        else if (string (*argv) == "outlinedens")
        {
            outlinedens = read_arg <string> (argv++);
            settings.linedensDIR = outlinedens;
        }
        else if (string (*argv) == "outpointspread")
        {
            outpointspread = read_arg <string> (argv++);
        }
        else if (string (*argv) == "outPNGs")
        {
            settings.pngDIR = read_arg <string> (argv++);
        }
        else if (string (*argv) == "objectDIR")
        {
            settings.objectDIR = read_arg <string> (argv++);
        }
        else if (string (*argv) == "s")
        {
            s = read_arg <double> (argv++);
        }
        else if (string (*argv) == "area")
        {
            s = 4234; //1337 for AREA
        }
        else if (string (*argv) == "euler")
        {
            s = 1;
        }
        else if (string (*argv) == "smooth")
        {
            smooth = read_arg <double> (argv++);
        }
        else if (string (*argv) == "erd")
        {
            erd = read_arg <double> (argv++);
        }
//        else if (string (*argv) == "CRPIX1")
//        {
//            crpix_source = read_arg <double> (argv++);
//        }
        else if (string (*argv) == "3d")
        {
            if (read_arg <string> (argv++) == "true") threeD = true;
        }
        else if (string (*argv) == "avg")
        {
            if (read_arg <string> (argv++) == "true") average = true;
        }
        else if (string (*argv) == "abs_avg")
        {
            if (read_arg <string> (argv++) == "true") absolute_avg = true;
        }
        else if (string (*argv) == "arg")
        {
            if (read_arg <string> (argv++) == "true") arg = true;
        }
        else if (string (*argv) == "makeBanana")
        {
            if (read_arg <string> (argv++) == "true") {
                makeBanana = true;
                make_minkmap = false;
            }
        }
        else if (string (*argv) == "makePeach")
        {
            if (read_arg <string> (argv++) == "true") {
                make_peach = true;
                make_minkmap = false;
            }
        }
        else if (string (*argv) == "makeHist")
        {
            if (read_arg <string> (argv++) == "true") {
                make_hist = true;
                make_minkmap = false;
            }
        }
        else if (string (*argv) == "makePNGs")
        {
            if (read_arg <string> (argv++) == "true") {
                makePNGs = true;
                make_minkmap = false;
            }
        }
        else if (string (*argv) == "combineRegions")
        {
            if (read_arg <string> (argv++) == "true") {
                combine_regions = true;
                make_minkmap = false;
            }
        }
        else if (string (*argv) == "monochrome")
        {
            if (read_arg <string> (argv++) == "true") {
                monochrome = true;
            }
            else monochrome = false;
        }
        else if (string (*argv) == "greenfile")
        {
            greenfile = read_arg <string> (argv++);
        }
        else if (string (*argv) == "bluefile")
        {
            bluefile = read_arg <string> (argv++);
        }
        else if (string (*argv) == "makeHedgehog")
        {
            if (read_arg <string> (argv++) == "true") {
                make_hedgehog = true;
                make_minkmap = false;
            }
        }
        else if (string (*argv) == "makeBubbles")
        {
            if (read_arg <string> (argv++) == "true") {
                make_bubbles = true;
                make_minkmap = false;
            }
        }
        else if (string (*argv) == "makePointspread")
        {
            if (read_arg <string> (argv++) == "true") {
                make_pointspread = true;
                make_minkmap = false;
            }
        }
        else if (string (*argv) == "makePattern")
        {
            string nextcommand = read_arg <string> (argv++);
            if (nextcommand == "linedensity" || nextcommand == "minkmap" || nextcommand.substr(0,6) == "infile") {
                make_patternFromImage = true;
                make_minkmap = false;
                linesorminkmap = nextcommand;
            }
            else
            {
                std::cerr << "illegal argument for pattern generation: " << nextcommand << "\n";
                return 1;
            }
        }
        else if (string (*argv) == "filewithoutWCS")
        {
            filewithoutWCS = true;
        }
        else if (string (*argv) == "readkey")
        {
            WCSkeynames.push_back(read_arg <string> (argv++));
        }
        else if (string (*argv) == "squaresize")
        {
            squaresize = read_arg <double> (argv++);
        }
        else
        {
            std::cerr << "illegal argument: " << *argv << "\n";
            return 1;
        }
    }
    
    if (smooth!=0 && squaresize != 2)
    {
        std::cout << "Change either squaresize or smooth, not both" << std::endl;
        return 1;
    }

    
    
    
    
    
    
    //Read original file
    FitsFile infile (infilename,WCSkeynames);
    if (boxesToExcludeName!="") //Exclude boxes with imaging errors, if file exists
    {
        std::vector<std::vector<double>> boxesToExclude = readTable(boxesToExcludeName);
        maskBoxes(infile, boxesToExclude);
    }
    
    //Add placeholder coordinate data
    if(filewithoutWCS)
    {
        infile.setKeyValue("RADECSYS","FK5");
        infile.setKeyValue("EQUINOX","2000"); 
        infile.setKeyValue("CTYPE1","RA---TAN");
        infile.setKeyValue("CTYPE2","DEC--TAN");
        infile.setKeyValue("CDELT1","-0.001");
        infile.setKeyValue("CDELT2","0.001");
        infile.setKeyValue("CRPIX1","0.");
        infile.setKeyValue("CRPIX2","0.");
        
        std::vector<std::vector<string>> infile_WCSdata = infile.returnWCSdata();
        
        writeImage(infile, infilename, infile_WCSdata);
        infilename += ".fits";
        FitsFile newinfile (infilename,WCSkeynames);
        infile = newinfile;
    }
    crpix_source = std::stod(infile.giveKeyvalue("CRPIX1"));

    
    
   
    if(combine_regions)
    {//read all bubbles in combination and combine, write objectfiles and maskfiles containing objectfiles and ds9 files for each combined bubble and R-readable and ds9-readable files with average positions and sizes of all combined bubbles
        std::vector<std::vector<std::vector<double>>> reg = combineRegions(combination);
        writeRegion(maskfilename, reg);
        writeBubblelist(reg, maskfilename);
    }
    
    if(make_patternFromImage)
    {
        std::string filename;
        std::string filenameWithoutFits;
        double factor = 1;
        if(linesorminkmap=="linedensity")
        {
            //read linedensity file
            std::string path = outlinedens;
            char buffer1 [15];
            int n = sprintf(buffer1,"_%g",line_scale);
            n++;
            filename = path+wavelength+"_lines_smooth="+std::to_string(smooth)+"_thresh"+buffer1+"_scalelength.fits";
            filenameWithoutFits = wavelength+"_lines_smooth="+std::to_string(smooth)+"_thresh"+buffer1+"_scalelength";
        }
        
        if(linesorminkmap=="minkmap")
        {
            //read q_2 minkmaps for this smooth value and one higher value
            std::string buffer2;
            if(wavelength == "halpha") buffer2 = "0.1_40_9";
            else if (wavelength == "sii") buffer2 = "0.045_40_10";
            else if (wavelength == "oiii") buffer2 = "0.1_40_9";
            else if (wavelength == "siioverha") buffer2 = "0.1_40_9";
            else
            {
                char buffer2a [40];
                int n = sprintf(buffer2a,"%g_%g_%d",min_thresh,max_thresh,num_thresh);
                n++;
                buffer2 = buffer2a;
            }
            std::string path = outminmap;
            std::string filenameAabs = wavelength+"_s="+std::to_string(s);
            if(absolute_avg) {
                filenameAabs += "_absavg_smooth=";
            }
            else if (average) {
                filenameAabs += "_avg_smooth=";
            }
            else {
                filenameAabs += "_smooth=";
            }
            
            std::string filenameB = "_thresh_";
            filenameB += buffer2;
            
            if(threeD && !average) std::cout << "No 3d files here! \n";
            
            //Read minkmap A
            filename = path + filenameAabs + std::to_string(smooth) + filenameB + ".fits";
            filenameWithoutFits = filenameAabs + std::to_string(smooth) + filenameB;
            factor = 20;
        }
        if(linesorminkmap.substr(0,6)=="infile")
        {
            filename = infilename;
            filenameWithoutFits = wavelength;
        }
        
        FitsFile image(filename,WCSkeynames);
        
        if(linesorminkmap.substr(0,6)=="infile" && smooth)
        {
            erosionDilation(image,erd);
            
            if(linesorminkmap=="infile_log") 
            {
                takeLogPositive(image);
                filenameWithoutFits += "_log";
                factor = 1./100;
            }
            
            //smooth and get shift in coordinates right
            smooth_map_plus(image,smooth);
            std::vector<std::vector<string>> minkmap_WCSdata = infile.returnWCSdata();
            image.setKeyValue("CRPIX1",std::to_string(stod(image.giveKeyvalue("CRPIX1")) - 0.5*(std::max(1.*squaresize,1.*smooth-smooth/6)-1)));
            image.setKeyValue("CRPIX2",std::to_string(stod(image.giveKeyvalue("CRPIX2")) + 0.5*(std::max(1.*squaresize,1.*smooth-smooth/6)-1)));
            
            filenameWithoutFits += "_smooth_"+std::to_string(smooth);
        }

        image *= factor;
        char bufferthresh [15];
        int n2 = sprintf(bufferthresh,"_%g",line_thresh);
        n2++;
        int ENpoints = 4000;
        if(line_thresh) ENpoints = 1000;
        imageToPointPattern(image, smooth, outfilename + "pattern_"+filenameWithoutFits+bufferthresh, ENpoints, line_thresh*factor);
    }
    
    if(make_hedgehog)
    {
        std::string buffer1;
        char buffer2a [40];
		int n = sprintf(buffer2a,"%g_%g_%d",min_thresh,max_thresh,num_thresh);
		n++;
		buffer1 = buffer2a;
        
        std::string path = outminmap;
        std::string filenameAarg = wavelength+"_s="+std::to_string(s)+"_arg";
        std::string filenameAabs = wavelength+"_s="+std::to_string(s);
        if(absolute_avg) {
            filenameAarg += "_absavg_smooth=";
            filenameAabs += "_absavg_smooth=";
        }
        else if (average) {
            filenameAarg += "_avg_smooth=";
            filenameAabs += "_avg_smooth=";
        }
        else {
            filenameAarg += "_smooth=";
            filenameAabs += "_smooth=";
        }
        
        std::string filenameB = "_thresh_";
        if(threeD) filenameB += buffer1;
        else filenameB += std::to_string(max_thresh);
        
        if(threeD && !average) std::cout << "No 3d files here! \n";
        
        //Read minkmaps
        std::string filename = path + filenameAabs + std::to_string(smooth) + filenameB + ".fits";
        FitsFile* absImage;
        FitsFile* argImage;
        try{
            absImage = new FitsFile(filename, WCSkeynames);
        } catch(const CCfits::FITS::CantOpen&) {
            std::cout << "Minkmap doesn't exist. Make new one...\n";
            makeMinkmap(infilename, infile, outminmap+wavelength, s, threeD, squaresize, average, smooth, absolute_avg, min_thresh, max_thresh, num_thresh, false);
            absImage = new FitsFile(filename, WCSkeynames);
        }
        filename = path + filenameAarg + std::to_string(smooth) + filenameB + ".fits";
        try{
            argImage = new FitsFile(filename, WCSkeynames);
        } catch(const CCfits::FITS::CantOpen&) {
            std::cout << "Minkmap doesn't exist. Make new one...\n";
            makeMinkmap(infilename, infile, outminmap+wavelength, s, threeD, squaresize, average, smooth, absolute_avg, min_thresh, max_thresh, num_thresh, true);
            argImage = new FitsFile(filename, WCSkeynames);
        }
        
        std::vector<std::vector<double>> lines = makeHedgehog(*absImage,*argImage,smooth,line_scale,wavelength);
        makeLinedensity(lines,smooth,line_scale,wavelength,*absImage,min_thresh,max_thresh,num_thresh);
        //bubbles(linedensity, smooth, (int)line_thresh, wavelength, line_scale, absImage, smoothB);
        make_bubbles = true;
        //return 0;
    }
    
    if(make_bubbles)
    {
        //read linedensity file
        std::string path = outlinedens;
        char buffer1 [15];
        int n = sprintf(buffer1,"_%g",line_scale);
        n++;
        std::string filename = path+wavelength+"_lines_smooth="+std::to_string(smooth)+"_thresh"+buffer1+"_scalelength.fits";
        
        FitsFile linedensity(filename,WCSkeynames);
        
        //read q_2 minkmaps for this smooth value and one higher value
        std::string buffer2;
		char buffer2a [40];
		n = sprintf(buffer2a,"%g_%g_%d",min_thresh,max_thresh,num_thresh);
		n++;
		buffer2 = buffer2a;
        
        path = outminmap;
        std::string filenameAabs = wavelength+"_s="+std::to_string(s);
        if(absolute_avg) {
            filenameAabs += "_absavg_smooth=";
        }
        else if (average) {
            filenameAabs += "_avg_smooth=";
        }
        else {
            filenameAabs += "_smooth=";
        }
        
        std::string filenameB = "_thresh_";
        filenameB += buffer2;
        
        if(threeD && !average) std::cout << "No 3d files here! \n";
        
        //Read minkmap A
        std::string filenameS = path + filenameAabs + std::to_string(smooth) + filenameB + ".fits";
        FitsFile smoothA(filenameS, WCSkeynames);
        
        int othersmooth = smooth;
        switch (smooth)
        {
            case 100: othersmooth = 300;
                    break;
            case 80: othersmooth = 250;
                    break;
            case 70: othersmooth = 250;
                    break;
            case 60: othersmooth = 200;
                    break;
            case 50: othersmooth = 150;
                    break;
            case 40: othersmooth = 100;
                    break;
            case 30: othersmooth = 100;
                    break;
            case 20: othersmooth = 60;
                    break;
            case 15: othersmooth = 50;
                    break;
            case 12: othersmooth = 40;
                    break;
            case 10: othersmooth = 30;
                    break;
            default: othersmooth = 3*smooth;
        }
        //Read minkmap B
        filenameS = path + filenameAabs + std::to_string(othersmooth) + filenameB + ".fits";
        FitsFile* smoothB;
        try{
            smoothB = new FitsFile(filenameS, WCSkeynames);
        } catch(const CCfits::FITS::CantOpen&) {
            std::cout << "Minkmap doesn't exist. Make new one...\n";
            makeMinkmap(infilename, infile, outminmap+wavelength, s, threeD, squaresize, average, othersmooth, absolute_avg, min_thresh, max_thresh, num_thresh, false);
            //makeMinkmap(infilename, infile, outminmap+wavelength, s, threeD, squaresize, average,      smooth, absolute_avg, min_thresh, max_thresh, num_thresh, arg);
            smoothB = new FitsFile(filenameS, WCSkeynames);
        }
        
        
        std::cout << "Starting bubbles for file " << filename << std::endl;
        
        
        bubbles(linedensity, smooth, (int)line_thresh, wavelength, line_scale, smoothA, *smoothB);
    }
    
    if(makePNGs)
    {
        if(erd) std::cout << "Warning! Eroding/dilating before making PNGs\n";
        
        FitsFile infile2 = infile, infile3 = infile;
        if(!monochrome) 
        {
            std::cout << "Creating colour PNGs of single objects with R = "+infilename+", G = "+greenfile+", B = "+bluefile+" ...\n";
            FitsFile green (greenfile,WCSkeynames);
            infile2 = green;
            FitsFile blue(bluefile,WCSkeynames);
            infile3 = blue;
        }
        if (erd)
        {
            erosionDilation(infile,erd);
            erosionDilation(infile2,erd);
            erosionDilation(infile3,erd);
        }
        if (maskfilename!="") //Include Boxes with wanted objects
        {
            std::ifstream maskfile (maskfilename);
            std::string objectname;
            // Every line in maskfile should contain one object. Object region files should be saved in settings.objectDIR
            while (std::getline(maskfile, objectname)) 
            {
                if(objectname.at(0) != '#')
                {
                    std::cout << objectname << std::endl;
                    makeRegionImages(infile, infile2, infile3, objectname, erd, monochrome, wavelength);
                }
            }
        }
    }
    
    if(make_pointspread)
    {
        std::ifstream histfiles (histfile);
        std::string objectlist;
        std::vector<std::vector<string>> minkmap_WCSdata = infile.returnWCSdata();
        
        // Every line in histfile should contain one objectlist. Region files should be saved in settings.objectDIR
        while (std::getline(histfiles, objectlist)) 
        {
            if(objectlist.at(0) != '#')
            {
                std::cout << objectlist << std::endl;
                string outname = outpointspread+"pointspread_"+objectlist+"_size_"+std::to_string(smooth)+"_gauss";
                BasicPhoto<double> results = makePointspread(infile, settings.objectDIR+objectlist+".reg", 3.*smooth/2, 1.*smooth/2,gauss);
                writeImage(results, outname, minkmap_WCSdata);
            }
        }
    }
    
    if(make_hist)
    {
        std::ifstream histfiles (histfile);
        std::string objectlist;
        
        
        // Every line in histfile should contain one objectlist. Region files should be saved in settings.objectDIR
        while (std::getline(histfiles, objectlist))
        {
            if(objectlist.at(0) != '#')
            {
                std::cout << objectlist << std::endl;
        
                makepointHist(infile, settings.objectDIR+objectlist+".reg", outfilename+"hist_"+objectlist+"_"+histappendix, min_thresh, max_thresh, num_thresh);
            }
        }
        makeimageHist(infile, outfilename+"hist__"+histappendix, min_thresh, max_thresh, num_thresh);
    }
    
    if(make_peach)
    {
        //Get first and second part of filename right
        char buffer1 [40];
        int n = sprintf(buffer1,"%g_%g_%d",min_thresh,max_thresh,num_thresh);
        n++;
        
        std::string path = outminmap;
        std::string filenameA = wavelength+"_s="+std::to_string(s);
        if(absolute_avg) filenameA += "_absavg_smooth=";
        else if (average) filenameA += "_avg_smooth=";
        else filenameA += "_smooth=";
        
        std::string filenameB = "_thresh_";
        if(threeD) filenameB += buffer1;
        else filenameB += std::to_string(max_thresh);
        
        if(threeD && !average) std::cout << "No 3d files here! \n";
        
        //Read all minkmaps
        std::vector<FitsFile> minkmaps;
        for(uint i=0; i<smooths.size(); i++)
        {
            int smooth = smooths.at(i);
            std::string filename = path + filenameA + std::to_string(smooth) + filenameB + ".fits";
            FitsFile infile(filename, WCSkeynames);
            minkmaps.push_back(infile);
        }
        
        std::cout << "Making peach...\n";
        if (maskfilename!="") //Include Boxes with wanted objects
        {
            std::ifstream maskfile (maskfilename);
            std::string objectname;
            
            
            // Every line in maskfile should contain one object. Object region files should be saved in settings.objectDIR
            while (std::getline(maskfile, objectname)) 
            {
                if(objectname.at(0) != '#')
                {
                    std::cout << objectname << std::endl;
                    makePeach(minkmaps, objectname, smooths, filenameA, filenameB);
                }
            }
        }
        
        else
        {
            makePeach(minkmaps, "", smooths, filenameA, filenameB);
        }
        
    }
    
    if(make_minkmap)
    {
        makeMinkmap(infilename, infile, outminmap+wavelength, s, threeD, squaresize, average, smooth, absolute_avg, min_thresh, max_thresh, num_thresh, arg);
    }
    
    if(makeBanana)
    {
        if(erd)
        {
            std::cout << "Warning! Eroding/dilating before banana!" << std::endl;
            erosionDilation(infile,erd);
        }
        std::cout << "Making Banana...\n";
        if (maskfilename!="") //Include Boxes with wanted objects
        {
            std::ifstream maskfile (maskfilename);
            std::string objectname;
            
            // Every line in maskfile should contain one object. Object region files should be saved in settings.objectDIR
            while (std::getline(maskfile, objectname)) 
            {
                if(objectname.at(0) != '#')
                {
                    std::cout << objectname << std::endl;
                    char buffer1 [40];
                    int n = sprintf(buffer1,"%g_%g_%d",min_thresh,max_thresh,num_thresh);
                    n++;
                    if(erd) originalBanana(objectname,infile, outfilename+"banana_"+wavelength+"_erd_mask_"+objectname+"_thresh_"+buffer1+".dat", min_thresh, max_thresh, num_thresh);
                    else originalBanana(objectname,infile, outfilename+"banana_"+wavelength+"_mask_"+objectname+"_thresh_"+buffer1+".dat", min_thresh, max_thresh, num_thresh);
                    
                }
            }
        }
        
        else
        {
            char buffer1 [40];
            int n = sprintf(buffer1,"%g_%g_%d",min_thresh,max_thresh,num_thresh);
            n++;
            originalBanana("",infile, outfilename+"banana_"+wavelength+"_mask_"+maskfilename+"_thresh_"+buffer1+".dat", min_thresh, max_thresh, num_thresh);
        }
    }


    
    
    return 0;
}
