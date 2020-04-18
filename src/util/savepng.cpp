#include "savepng.hpp"
#include "fitsfile.hpp"

using namespace papaya2;

template<typename PHOTO> //This function assumes than a logarithmic scale is the way to go
std::vector<unsigned char> prepareData (PHOTO& inputR, PHOTO& inputG, PHOTO& inputB)
{
    double maxR = std::abs(inputR.max_value());
    double maxG = std::abs(inputG.max_value());
    double maxB = std::abs(inputB.max_value());
    
    double factorR = 255./std::log(maxR+1);
    double factorG = 255./std::log(maxG+1);
    double factorB = 255./std::log(maxB+1);
                                       
    int width = inputR.width();
    int height = inputR.height();
    std::vector<unsigned char> output;
    for(int i=0; i<height; i++)
    for(int j=0; j<width; j++)
    {
        output.push_back((int)(std::min(std::log(inputR(j,i)+1)*factorR,255.)) );
        output.push_back((int)(std::min(std::log(inputG(j,i)+1)*factorG,255.)) );
        output.push_back((int)(std::min(std::log(inputB(j,i)+1)*factorB,255.)) );
        output.push_back(255);
    }
    return output;
}

template // explicit instantiation
std::vector<unsigned char> prepareData (Photo& inputR, Photo& inputG, Photo& inputB);

template // explicit instantiation
std::vector<unsigned char> prepareData (FitsFile& inputR, FitsFile& inputG, FitsFile& inputB);

void PNGtoBasicPhoto(std::vector<unsigned char>& image, unsigned width, unsigned height, papaya2::BasicPhoto<double>& outputR, papaya2::BasicPhoto<double>& outputG, papaya2::BasicPhoto<double>& outputB, bool color/* = false*/)
{
    outputR.set_coordinates(0,0,width,height,width,height);
    if(color)
    {
        outputG.set_coordinates(0,0,width,height,width,height);
        outputB.set_coordinates(0,0,width,height,width,height);
    }
    
    //Read all colors separately
    std::vector<unsigned char> bufferR;
    std::vector<unsigned char> bufferG;
    std::vector<unsigned char> bufferB;
    for(uint i=0; i<image.size(); i+=4)
    {
        bufferR.push_back(image.at(i));
        if(color)
        {
            bufferG.push_back(image.at(i+1));
            bufferB.push_back(image.at(i+2));
        }
    }
    
    //Enter into Photo
    for (uint i = 0; i < width; ++i)
    for (uint j = 0; j < height; ++j)
    {
        outputR(i,height-j) = bufferR.at((height - j - 1) * width + i);
        if(color)
        {
            outputG(i,height-j) = bufferG.at((height - j - 1) * width + i);
            outputB(i,height-j) = bufferB.at((height - j - 1) * width + i);
        }
    }
}
