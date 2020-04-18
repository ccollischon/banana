#ifndef _banana_savepng_hpp_

#define _banana_savepng_hpp_

#include "src/papaya2.hpp"
#include "util/lodepng.h"

void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
    //Encode the image
    unsigned error = lodepng::encode(filename, image, width, height);

    //if there's an error, display it
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

std::vector<unsigned char> decodeOneStep(const char* filename, unsigned &width, unsigned& height)
{
    std::vector<unsigned char> image; //the raw pixels

    //decode
    unsigned error = lodepng::decode(image, width, height, filename);

    //if there's an error, display it
    if(error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    return image;
}

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

void PNGtoBasicPhoto(std::vector<unsigned char>& image, unsigned width, unsigned height, BasicPhoto<double>& outputR, BasicPhoto<double>& outputG, BasicPhoto<double>& outputB, bool color = false)
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


template<typename PHOTO>
void writeColorPNG(std::string filename, PHOTO& inputR, PHOTO& inputG, PHOTO& inputB)
{
    std::vector<unsigned char> imageData = prepareData(inputR, inputG, inputB);
    unsigned width = inputR.width();
    unsigned height = inputR.height();
    filename = filename;
    encodeOneStep(filename.c_str(), imageData, width, height);
}

template<typename PHOTO>
void writeMonoPNG(std::string filename, PHOTO& input)
{
    std::vector<unsigned char> imageData = prepareData(input, input, input);
    unsigned width = input.width();
    unsigned height = input.height();
    encodeOneStep(filename.c_str(), imageData, width, height);
}



void testreadPNG()
{
    unsigned width,height;
    std::vector<unsigned char> imageData = decodeOneStep("testimage.png",width,height);
    BasicPhoto<double> imageR, imageG, imageB;
    PNGtoBasicPhoto(imageData, width, height, imageR, imageG, imageB, true);
    
    writeColorPNG("testoutput.png", imageR, imageG, imageB);
}

#endif
