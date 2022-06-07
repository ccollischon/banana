#ifndef _banana_hedgehog_hpp_

#define _banana_hedgehog_hpp_
// A class for bubble detection related stuff


#include "src/papaya2.hpp"
#include "src/photofunctions.hpp"
#include "util/fitsfile.hpp"

#include <iomanip>

std::vector<std::vector<double>> makeHedgehog(FitsFile& infileAbs, FitsFile& infileArg, int smooth, double thresh, string wavelength)
{
    double length = smooth/2;
    double pi = 3.14159265358979;
    int end1 = infileAbs.height();
    int end2 = infileAbs.width();
    std::vector<std::vector<double>> lines;
    std::cout << "Iterating over " +wavelength+ " file smooth " << smooth << "..." << std::endl; 
    
    char buffer1 [15];
    int n = sprintf(buffer1,"_%g",thresh);
    n++;
    std::ofstream ofs (settings.resultDIR+"angles_"+wavelength+"_smooth="+std::to_string(smooth)+"_thresh"+buffer1+"_scalelength.reg");
    ofs << "# Region file format: DS9 version 4.1 \n global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n image \n";
    
    int lineodd = 0; //helper to write only every 2nd line into region file
    int stepsize = std::max(smooth/6,1);
    //if(smooth<100) stepsize *= 2;
    for (int j = stepsize/2; j < end1; j+=stepsize)
    for (int i = stepsize/2; i < end2; i+=stepsize)
    {
        std::vector<double> line;
        if(infileAbs(i,j)>0.01)
        {
            double angle = std::fmod(0.5*(pi-infileArg(i,j)), pi);
            std::vector<int> coord = infileAbs.giveUnshiftedds9Coord(i,j);
            length = smooth/2*(infileAbs(i,j)/thresh);
            line.push_back(1.*coord.at(0)-length*cos(angle));
            line.push_back(1.*coord.at(1)-length*sin(angle));
            line.push_back(1.*coord.at(0)+length*cos(angle));
            line.push_back(1.*coord.at(1)+length*sin(angle));
            lines.push_back(line);
            
            int color = std::min((int)(200*infileAbs(i,j)/thresh),255);
            if(smooth>=100)  ofs << "line(" << line.at(0) << "," << line.at(1) << "," << line.at(2) << "," << line.at(3) << ") # line=0 0 color=#"<< std::hex<<std::setw(2)<<std::setfill('0') << color << "0000" << std::endl;
            else if(lineodd%2)   ofs << "line(" << line.at(0) << "," << line.at(1) << "," << line.at(2) << "," << line.at(3) << ") # line=0 0 color=#00"<< std::hex<<std::setw(2)<<std::setfill('0') << color << "00" << std::endl;
            lineodd++;
        }
    }
    
    ofs.close();
    
    return lines;
}

void makeLinedensity(const std::vector<std::vector<double>>& lines, int smooth, double thresh, string wavelength, FitsFile& base, double mint, double maxt, int numt)
{
    std::cout<< "Starting line density... \n";
    FitsFile lineImage = base;
    //Set lineImage to zero
    lineImage *= 0;
    
    //Empty table for counting line crossings
    
    int stepsize = std::max(smooth/6,1);
    std::vector<double> zeroline((int)(lineImage.height()/stepsize),0.);
    std::vector<std::vector<double>> values((int)(lineImage.width()/stepsize),zeroline);
    
    int valwidth = values.size();
    int valheight = values.at(0).size();
    //Iterate over every line and find crossings
    #pragma omp parallel for
    for(uint ijk=0; ijk<lines.size(); ijk++)
    {
        auto line = lines.at(ijk);
        std::vector<std::vector<double>> valuesHere((int)(lineImage.width()/stepsize),zeroline);
        std::vector<int> pix1 = base.giveShiftedCoord(line.at(0),line.at(1));
        std::vector<int> pix2 = base.giveShiftedCoord(line.at(2),line.at(3));
        //Find all crossed boxes, set to one
        vec_t pix1v = {(double)(pix1.at(0)), (double)(pix1.at(1))};
        vec_t pix2v = {(double)(pix2.at(0)), (double)(pix2.at(1))};
        vec_t delta = pix2v-pix1v;
        for(double d=0.; d<=1; d+=0.02)
        {
            vec_t newPos = pix1v+delta*d;
            int xPos = newPos[0]/stepsize;
            int yPos = newPos[1]/stepsize;
            if(xPos>=0 && yPos>=0 && (xPos<valwidth) && (yPos<valheight))   valuesHere.at(xPos).at(yPos) = 1;
        }

        //Add to global count
        for(uint i=0; i<valuesHere.size(); i++)
        for(uint j=0; j<valuesHere.at(i).size(); j++)
        {
            values.at(i).at(j)+=valuesHere.at(i).at(j);
        }
    }
    std::cout<< "Filling fits file with values... \n";
    //Fill into minkmap
    #pragma omp parallel for
    for (int j = 0; j < (lineImage.height()/stepsize)-1; j++)
    for (int i = 0; i < (lineImage.width()/stepsize)-1; i++)
    {
        double newValue = values.at(i).at(j);
        for(int kstep=0; kstep<stepsize; kstep++)
        for(int lstep=0; lstep<stepsize; lstep++)
        {
            lineImage(i*stepsize+kstep+1,j*stepsize+lstep+1) = newValue;
        }
    }
    //Write with fitting comments
    char buffer1 [15];
    int n = sprintf(buffer1,"_%g",thresh);
    n++;
    std::vector<std::vector<string>> minkmap_WCSdata = lineImage.returnWCSdata();
    char buffer2a [51];
    n = sprintf(buffer2a,"%d thresholds %g-%g",numt,mint,maxt);
    n++;
    minkmap_WCSdata.at(0).push_back("Comment");
    minkmap_WCSdata.at(1).push_back("Source image: "+wavelength+" smooth "+std::to_string(smooth) +" s 2 minkmap with "+buffer2a);
    minkmap_WCSdata.at(0).push_back("HISTORY");
    minkmap_WCSdata.at(1).push_back("Line densities for above image created with Minkowski absolute value threshold lineScale="+std::to_string(thresh));
    string outname = settings.linedensDIR+wavelength+"_lines_smooth="+std::to_string(smooth)+"_thresh"+buffer1+"_scalelength";
    writeImage(lineImage, outname,minkmap_WCSdata);
    
    //return lineImage;
}

void bubbles(FitsFile& infile, int smooth, int thresh, string wavelength, double line_thresh, FitsFile& smoothA, FitsFile& smoothB)
{
    int width = infile.width();
    int height = infile.height();
    int length = 2*smooth;
    
    char buffer1 [40];
    int n = sprintf(buffer1,"_%g_bubblethresh_%d",line_thresh,thresh);
    n++;
    
    //Sorted: slope of q_2 checked before approving
    std::ofstream ofs (settings.resultDIR+"bubbles_sorted_"+wavelength+"_smooth="+std::to_string(smooth)+"_thresh"+buffer1+".reg");
    std::ofstream ofs_objectnames("bubbles_sorted_regionnames_"+wavelength+"_smooth="+std::to_string(smooth)+"_thresh"+buffer1);
    ofs << "# Region file format: DS9 version 4.1 \n global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n image \n";
    
    std::cout << "Searching for bubbles..." << std::endl;
    int stepsize = std::max(smooth/6,1);
    int objectnumber = 0;
    for (int j = 1; j < (height); j+=stepsize)
    for (int i = 1; i < (width); i+=stepsize)
    {
        double valueHere = infile(i,j);
        if(valueHere>thresh)
        {
            std::cout << "A";
            int inotBelowZero = std::max(0,i-stepsize);
            int jnotBelowZero = std::max(0,j-stepsize);
            //Values of neighbors
            double tl = infile(inotBelowZero,jnotBelowZero), tm = infile(inotBelowZero,j), tr = infile(inotBelowZero,j+stepsize);
            double ml = infile(i,jnotBelowZero), mr = infile(i,j+stepsize);
            double bl = infile(i+stepsize,jnotBelowZero), bm = infile(i+stepsize,j), br = infile(i+stepsize,j+stepsize);
            
            std::cout << "B";
                // Check if value here higher than surroundings, if neighbors of same value: choose the one top-right
            if((tl<=valueHere) && (tm<valueHere) && (tr<valueHere) && (ml<=valueHere) && (mr<valueHere) && (bl<=valueHere) && (bm<=valueHere) && (br<valueHere))
            {
                std::cout << "C";
                //Check if slope of peach graph negative for higher lengthscales
                double q2_avgA, q2_avgB;
                double pix1 = i+stepsize;
                double pix2 = j+stepsize;
                double shiftA = smoothA.giveShift();
                double shiftB = smoothB.giveShift();
                std::vector<int> ds9coord = infile.giveUnshiftedds9Coord(pix1, pix2);
                std::vector<double> boxA = {std::max(0.,1.*ds9coord.at(0)-shiftA), 1.*ds9coord.at(1)+shiftA, 1.*length/4, 1.*length/4};
                std::vector<double> boxB = {std::max(0.,1.*ds9coord.at(0)-shiftB), 1.*ds9coord.at(1)+shiftB, 1.*length/4, 1.*length/4};
                std::vector<std::vector<double>> boxboxA(1,boxA);
                std::vector<std::vector<double>> boxboxB(1,boxB);
                std::cout << "D";
                q2_avgA = averageBoxes(smoothA,boxboxA);
                q2_avgB = averageBoxes(smoothB,boxboxB);
                std::cout << "E";
                
                if(q2_avgA+0.05 > q2_avgB)
                {
                    
                std::cout << "F";
                    std::ofstream ofs_object(settings.objectDIR+"bubble_sorted_"+wavelength+"_smooth="+std::to_string(smooth)+"_thresh"+buffer1+"_nr_"+std::to_string(objectnumber));
                    
                    ofs << "box(" << std::to_string(ds9coord.at(0)) <<","<< std::to_string(ds9coord.at(1)) <<","<< std::to_string(length) <<","<< std::to_string(length) <<",0)"<< std::endl;
                    ofs_object << std::to_string(ds9coord.at(0)) <<" "<< std::to_string(ds9coord.at(1)) <<" "<< std::to_string(length) << " " << std::to_string(length) << std::endl;
                    ofs_objectnames << "bubble_sorted_"<<wavelength<<"_smooth="<<std::to_string(smooth)<<"_thresh"<<buffer1<<"_nr_"<<std::to_string(objectnumber) << std::endl;
                    objectnumber++;
                std::cout << "G\n";
                }
            }    //... or if value same as left or lower middle neighbors
            //else if(tl==valueHere && ml==valueHere && bl==valueHere || bm==valueHere)
            //{
            //    double pix1 = i;
            //    double pix2 = j+stepsize/2;
            //    std::vector<int> ds9coord = infile.giveUnshiftedds9Coord(pix1, pix2);
            //    ofs << "box(" << std::to_string(ds9coord.at(0)) <<","<< std::to_string(ds9coord.at(1)) <<","<< std::to_string(length) <<","<< std::to_string(length) <<",0)"<< std::endl;
            //}
        }
    }
    std::cout << std::endl;
    ofs.close();
}



std::vector<std::vector<std::vector<double>>> combineRegions(string bubbletypes) //Combine all regions from objectlists listed in bubbletypes-named file into actual bubbles
{
    std::ifstream bubblefile (bubbletypes);
    std::string listname;
    std::vector<std::vector<std::vector<double>>> bubbles;
    
    // Every line in bubblefile should contain one objectlist with bubbles rising in size. One objectlist has the structure of a usual maskfile.
    while (std::getline(bubblefile, listname))
    {
        if(listname.at(0) != '#')
        {
            std::vector<std::vector<double>> bubblesThisSize;
            std::ifstream maskfile(listname);
            std::string bubblename;
            while (std::getline(maskfile, bubblename))
            {
                if(bubblename.at(0) != '#')
                {
                    std::vector<std::vector<double>> thisBubble = readTable(settings.objectDIR+bubblename);
                    //bubblesThisSize.push_back(thisBubble.at(0));
                    bubblesThisSize.insert(bubblesThisSize.end(), thisBubble.begin(), thisBubble.end());
                }
            }
            bubbles.push_back(bubblesThisSize);
        }
    }
    
    std::vector<std::vector<std::vector<double>>> combinedBubbles;
    //Have read all data, now combine: start at lowest size, then check if centers of this or next higher size bubbles are contained in this bubble
    //If yes, combine. delete larger bubble from its list, and continue with increased size to higher size.
    for(uint bubbleSizeCounter = 0; bubbleSizeCounter<bubbles.size(); bubbleSizeCounter++) //iterate over all sizes
    {
        std::vector<std::vector<std::vector<double>>> bubblesStartingHere;
        std::cout << bubbleSizeCounter << std::endl;
        
        std::vector<std::vector<double>>::iterator remainingBubbles = bubbles.at(bubbleSizeCounter).begin();
        while(remainingBubbles < bubbles.at(bubbleSizeCounter).end()) //Get all bubbles in this size class
        {
            std::vector<double> coordinatesOfBubble = *remainingBubbles; //Bubble currently working on
            
            std::cout << "First bubble coordinate:" << coordinatesOfBubble.at(0) << std::endl;
            std::vector<std::vector<double>> newBubbleCombined; //List for this bubble
            newBubbleCombined.push_back(coordinatesOfBubble); //Add this one to list
            //
            //Check subsequent bubbles this size and combine if necessary
            //
            double bubblesize = coordinatesOfBubble.at(2)/2;
            std::vector<std::vector<double>>::iterator posThisSize = remainingBubbles+1; //iterate over this size while removing stuff
            while(posThisSize < bubbles.at(bubbleSizeCounter).end())
            {
                std::vector<double> nextbubblehere = *posThisSize;
                if(std::max(std::abs(nextbubblehere.at(0)-coordinatesOfBubble.at(0)),std::abs(nextbubblehere.at(1)-coordinatesOfBubble.at(1))) < bubblesize)
                {
                    newBubbleCombined.push_back(nextbubblehere);
                    posThisSize = bubbles.at(bubbleSizeCounter).erase(posThisSize);
                    std::cout << "Absorbed same sized bubble. New vector size:" << bubbles.at(bubbleSizeCounter).size() << std::endl;
                }
                else posThisSize++;
            }
            
            //
            //See if bubbles in next higher size class close to this one. If none exists, set to false.
            //
            bool notFinished = true;
            uint compareSize = bubbleSizeCounter+1; //increase size
            while(notFinished && compareSize<bubbles.size())
            {
                notFinished = false;
                
                std::vector<std::vector<double>>::iterator it = bubbles.at(compareSize).begin();
                while(it != bubbles.at(compareSize).end())
                {
                    std::vector<double> nextbubble = *it;
                    if(std::max(std::abs(nextbubble.at(0)-coordinatesOfBubble.at(0)),std::abs(nextbubble.at(1)-coordinatesOfBubble.at(1))) < bubblesize)
                    {
                        notFinished = true;
                        bubblesize = nextbubble.at(2)/2;
                        newBubbleCombined.push_back(nextbubble);
                        it = bubbles.at(compareSize).erase(it);
                    std::cout << "Absorbed higher sized bubble" << std::endl;
                    }
                    else it++;
                }
                compareSize++;
            }
            bubblesStartingHere.push_back(newBubbleCombined);
            remainingBubbles++; //step forward in iterator
        }
        
        //Add bubbles this size to list
        combinedBubbles.insert(combinedBubbles.end(), bubblesStartingHere.begin(), bubblesStartingHere.end());
    }
    return combinedBubbles;
}

void writeRegion(string filename, std::vector<std::vector<std::vector<double>>> bubbles)
{ //Write both ds9 region file and objectfile for each bubble and maskfile listing all objectfiles
    std::ofstream maskfile(filename);
    for(uint i=0; i<bubbles.size(); i++)
    {
        std::ofstream ds9file (settings.resultDIR+filename+"_"+std::to_string(i)+".reg");
        ds9file << "# Region file format: DS9 version 4.1 \n global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n image \n";
        std::ofstream objectfile (settings.objectDIR+filename+"_"+std::to_string(i));
        maskfile << filename+"_"+std::to_string(i) << std::endl;
        
        for(auto coord : bubbles.at(i))
        {
            double length = coord.at(2);
            ds9file << "box(" << std::to_string(coord.at(0)) <<","<< std::to_string(coord.at(1)) <<","<< std::to_string(length) <<","<< std::to_string(length) <<",0)"<< std::endl;
            objectfile << std::to_string(coord.at(0)) <<" "<< std::to_string(coord.at(1)) <<" "<< std::to_string(length) << " " << std::to_string(length) << std::endl;
        }
        ds9file.close();
        objectfile.close();
    }
    maskfile.close();
}

#endif
