#ifndef _banana_savepng_hpp_

#define _banana_savepng_hpp_

#include <vector>
#include "../papaya2.hpp"
#include "lodepng.h"


//Example function from lodepng
inline std::vector<unsigned char> decodeOneStep(const char* filename, unsigned &width, unsigned& height)
{
    std::vector<unsigned char> image; //the raw pixels
    //decode
    unsigned error = lodepng::decode(image, width, height, filename);
    //if there's an error, display it
    if(error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    return image;
}

template<typename PHOTO> //This function assumes than a logarithmic scale is the way to go
std::vector<unsigned char> prepareData (PHOTO& inputR, PHOTO& inputG, PHOTO& inputB);

void PNGtoBasicPhoto(std::vector<unsigned char>& image, unsigned width, unsigned height, papaya2::BasicPhoto<double>& outputR, papaya2::BasicPhoto<double>& outputG, papaya2::BasicPhoto<double>& outputB, bool color = false);


template<typename PHOTO>
void writeColorPNG(std::string filename, PHOTO& inputR, PHOTO& inputG, PHOTO& inputB);

template<typename PHOTO>
void writeMonoPNG(std::string filename, PHOTO& input);


/*
inline void testreadPNG()
{
    unsigned width,height;
    std::vector<unsigned char> imageData = decodeOneStep("testimage.png",width,height);
    papaya2::BasicPhoto<double> imageR, imageG, imageB;
    PNGtoBasicPhoto(imageData, width, height, imageR, imageG, imageB, true);
    
    writeColorPNG("testoutput.png", imageR, imageG, imageB);
}
*/

#endif
