#ifndef _banana_peach_hpp_

#define _banana_peach_hpp_

#include "src/papaya2.hpp"
#include "src/photofunctions.hpp"
#include "util/fitsfile.hpp"

void makePeach(std::vector<FitsFile>& infiles, std::string objectname, std::vector<int> smooths, string filenameA, string filenameB)
{
    
    //Read pixel coordinates of object
    std::vector<std::vector<double>> boxesToInclude = readTable(settings.objectDIR+objectname);
    
    Datafile outfile(settings.resultDIR+"peach_"+objectname +"_"+filenameA+std::to_string(smooths.at(0)) + "_" + std::to_string(smooths.back())+filenameB+".dat");
    
    outfile.comment("smooth averagedValue");
    
    for(uint_fast8_t i=0; i<smooths.size(); i++)
    {
        int smooth  = smooths.at(i);
        
        //Shift boxes correctly
        std::vector<std::vector<double>> boxes = boxesToInclude;
        double shift = 0.5*(5*smooth/6 - 1);
        for(uint_fast8_t k=0; k<boxes.size(); k++)
        {
            boxes.at(k).at(0) -= shift;
            boxes.at(k).at(1) += shift;
            //std::cout << boxes.at(k).at(0) << " " << boxes.at(k).at(1) << std::endl;
        }
        
        //Average
        double av = averageBoxes(infiles.at(i), boxes);
        outfile << smooth << av << std::endl;
        //std::vector<std::vector<string>> minkmap_WCSdata = infiles.at(i).returnWCSdata();
        //includeBoxes(infiles.at(i),boxes);
        //writeImage(infiles.at(i),"test"+std::to_string(smooth),minkmap_WCSdata);
    }
}


#endif
