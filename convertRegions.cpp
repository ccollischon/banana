#include "readarg.hpp"
#include "readdata.hpp"
#include "papaya2.hpp"
#include <sstream>
#include <vector>
#include <valarray>
using namespace papaya2;

int main (int argc, const char **argv)
{
    std::string boxFileName, nameFileName
    // process command-line arguments
    for (++argv; *argv; ++argv)
    {
        if(std::string (*argv) == "names")
            nameFileName = read_arg <std::string> (argv++);
        else if(std::string (*argv) == "boxes")
            boxFileName = read_arg <std::string> (argv++);
        else
        {
            std::cerr << "illegal argument: " << *argv << "\n";
            return 1;
        }
    }
    
    vector<std::vector<double>> boxes = readTable(boxfilename);
}
