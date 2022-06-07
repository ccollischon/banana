#ifndef _banana_sugarsnappea_

#define _banana_sugarsnappea_

#include <random>
//std::vector<double> 


double randInInterval(double min_val, double max_val)
{
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_real_distribution<double> distr(min_val, max_val); // define the range

    return distr(eng);
}



// some pointspread functions
double gauss(double x, double s)
{
    double retval = exp(-x*x/(2*s*s));
    retval /= 2*3.14159*s*s*(1-exp(-2));
    return retval;
}
double linDecay(double x, double halfwidth)
{
    double retval = (2*halfwidth-x);
    retval /= 2*3.14159*4/3*halfwidth*halfwidth*halfwidth;
    return retval;
}
double cauchy(double x, double s)
{
    double retval = s/(3.14159*(s*s+x*x));
    retval /= s*log(5);
    return retval;
}
double constant(double x, double s)
{
    x+=s;
    return 1;
}

void makepointHist(FitsFile& infile, string listname, string outname, double min_thresh, double max_thresh, double num_thresh)
{ // Make histogram of brighntess of image at points
    std::vector<std::vector<double>> hist;
    std::vector<double> increments = linspace (min_thresh, max_thresh, num_thresh, true);
    double inc = (max_thresh - min_thresh)/(num_thresh-1);
    hist.push_back(increments);
    std::vector<double> values(increments.size());
    
    //Read coordinates of objects
    std::vector<std::vector<double>> objects = readPixels(listname); //indices of objects must be accessed .at(line).at(column), column gives ra or dec
    
    for(auto coordinates : objects)
    {
        double val = infile.atds9pix(coordinates.at(0), coordinates.at(1));
        if(val>=min_thresh && val <= max_thresh)
            values.at((val-min_thresh)/inc)++;
    }
    hist.push_back(values); //indices of hist must be accessed .at(column).at(line), column gives either threshold or number
    
    Datafile outfile (outname);
    outfile.comment("thresholds    # of objects");
    for(uint i=0; i<num_thresh; i++)
    {
        outfile << hist.at(0).at(i) << hist.at(1).at(i) << std::endl;
        std::cout << hist.at(0).at(i) << " " << hist.at(1).at(i) << std::endl;
    }
}

void makeimageHist(FitsFile& infile, string outname, double min_thresh, double max_thresh, double num_thresh)
{ // make histogram of whole image
    std::vector<std::vector<double>> hist;
    std::vector<double> increments = linspace (min_thresh, max_thresh, num_thresh, true);
    double inc = (max_thresh - min_thresh)/(num_thresh-1);
    hist.push_back(increments);
    std::vector<double> values(increments.size());
    
    for(int i=0; i<infile.width(); i++)
    for(int j=0; j<infile.height(); j++)
    {
        double val = infile.at(i, j);
        if(val>=min_thresh && val <= max_thresh)
            values.at((val-min_thresh)/inc)++;
    }
    hist.push_back(values); //indices of hist must be accessed .at(column).at(line), column gives either threshold or number
    
    Datafile outfile (outname);
    outfile.comment("thresholds    # of pixels");
    for(uint i=0; i<num_thresh; i++)
    {
        outfile << hist.at(0).at(i) << hist.at(1).at(i) << std::endl;
        std::cout << hist.at(0).at(i) << " " << hist.at(1).at(i) << std::endl;
    }
}

template<typename PHOTO>
BasicPhoto<double> makePointspread(PHOTO& infile, std::string listname, double radius, double width, double (*f)(double,double))
{// Spreads points in image to larger blobs with specified radius and function
    BasicPhoto<double> mask;
    mask.set_coordinates (0, 0, infile.width(), infile.height(), infile.width(), infile.height());
    //Read coordinates of objects
    std::vector<std::vector<double>> objects = readPixels(listname); //indices of objects must be accessed .at(line).at(column), column gives ra or dec
    double rsquared = radius*radius;
    int begin1 = 1,begin2 = 1;
    int end1 = infile.width();
    int end2 = infile.height();
    //std::cout << width << std::endl;
    for (auto object : objects)
    {
        double i = object.at(0)-0.5;
        double j = 1.*end2-object.at(1)+1.5;
        for(int k=-(radius); k<=(radius); k++)
        for(int l=-(radius); l<=(radius); l++)
        {
            double dist = pow(k,2)+pow(l,2);
            if( dist <= rsquared && (i+k>begin1) && i+k<end1 && j+l>begin2 && j+l<end2)
            {
                dist = pow(dist,0.5);
                mask(i+k,j+l) += (*f)(dist,width);
            }
        }
    }
    //mask/= (double)objects.size();
    return mask;
}

double getMeanColumn(std::vector<std::vector<double>> bubbleboxes, uint column)
{ //Takes table with one vector for each line and calculates the mean value of columns
    double result = 0.;
    for(auto line:bubbleboxes)
    {
        result+=line.at(column);
    }
    result /= bubbleboxes.size();
    return result;
}

void writeBubblelist(std::vector<std::vector<std::vector<double>>> allBubbles, std::string filename)
{ //Takes complete list of bubbles, combines all info for each bubble, writes table containing only centers (R readable) and ds9-readable file containing centers and average size
    std::ofstream ofs(settings.resultDIR+"R_readable_"+filename);
    std::ofstream ds9file(settings.resultDIR+filename+"_all.reg");
    ds9file << "# Region file format: DS9 version 4.1 \n global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n image \n";
    for(auto bubblecombo:allBubbles)
    {
        double x = getMeanColumn(bubblecombo,0);
        double y = getMeanColumn(bubblecombo,1);
        double size = getMeanColumn(bubblecombo,2);
        ofs << x << "\t" << y << std::endl;
        ds9file << "box(" << x <<","<< y <<","<< size <<","<< size <<",0)"<< std::endl;
    }
    ofs.close();
    ds9file.close();
}


#endif
