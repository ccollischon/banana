#ifndef _banana_readdata_hpp_

#define _banana_readdata_hpp_

#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>

    //Read table with Center RA, DEC, Size RA, DEC in px (or any table containing only numbers/numbers except characters at beginnig to be searched and discarded)
    //onlyStartingWith means selecting only those lines starting with the given chars, and those chars will not be added to the vector<vector>. Intended for reading Harris & Zaritsky(2009) data
std::vector<std::vector<double>> readTable(std::string filename, std::string onlyStartingWith="")
{
    std::vector<std::vector<double> > numbers;
    std::ifstream infile (filename);
    std::string temp;
    
    uint startingLength = onlyStartingWith.length();

    while (std::getline(infile, temp))
    {
        if(onlyStartingWith!="")
        {
            std::string startingchars = temp.substr(0,startingLength);
            if(startingchars == onlyStartingWith)
            {
                temp = temp.substr(startingLength+1,temp.length()-startingLength);
                std::istringstream buffer(temp);
                std::vector<double> line((std::istream_iterator<double>(buffer)),std::istream_iterator<double>());
                numbers.push_back(line);
            }
        }
        else
        {
            std::istringstream buffer(temp);
            std::vector<double> line((std::istream_iterator<double>(buffer)),std::istream_iterator<double>());
            numbers.push_back(line);
        }
    }
    return numbers;
}

    //Read table with hh:mm:ss,dd:mm:ss in ds9 format, gives result in deg, degnothms to be set when table already in deg
std::vector<std::vector<double>> readCoordinates(std::string filename, bool degnothms = false)
{
    std::vector<std::vector<double> > numbers;
    std::ifstream infile (filename);
    std::string temp;
    
    //Discard first 3 lines in ds9 region file
    std::getline(infile, temp);
    std::getline(infile, temp);
    std::getline(infile, temp);
    
    while (std::getline(infile, temp))
    {
        std::vector<double> line;
        std::size_t found = temp.find(",");
        if(!degnothms)
        {
            line = {0.,0.};
            line.at(0)= 150.*(temp.at(6)-'0') + 15.*(temp.at(7)-'0') + 2.5*(temp.at(9)-'0') + 0.25*(temp.at(10)-'0') + 0.041667*(temp.at(12)-'0') + 0.0041667*(temp.at(13)-'0');
            
            line.at(1)= 10.*(temp.at(2+found)-'0') + 1.*(temp.at(3+found)-'0') + 0.16667*(temp.at(5+found)-'0') + 0.016667*(temp.at(6+found)-'0') + 2.7778e-03*(temp.at(8+found)-'0') + 2.7778e-04*(temp.at(9+found)-'0');
            line.at(1)*=-1.;
            std::cout << temp.at(2+found)-'0' << " " << temp.at(3+found)-'0' << std::endl;
            std::cout << line.at(0) << " " << line.at(1) << std::endl;
        }
        else
        {
            // An object of regex for pattern to be searched 
            //std::regex r("^(-?)(0|([1-9][0-9]*))(\\.[0-9]+)?$"); 
            // flag type for determining the matching behavior 
            // here it is for matches on 'string' objects 
            //std::smatch m; 
            // regex_search() for searching the regex pattern 
            // 'r' in the string 'temp'. 'm' is flag for determining 
            // matching behavior. 
            
            //if (std::regex_search(temp, m, r))
            //{
            //    std::cout << m.str(1) << std::endl;
            //    line.push_back(std::stod(m.str(1)));
            //}
            std::string a = temp.substr(6,9);
            line.push_back(std::stod(a));
            a = temp.substr(found+1,9);
            line.push_back(std::stod(a));
            std::cout << line.at(0) << " " << line.at(1) << std::endl;
        }
        numbers.push_back(line);
    }
    return numbers;
}


    //Read table with rapix,decpix in ds9 format (or any two comma-separated numbers in brackets), gives result in pixels of original lmcdata
std::vector<std::vector<double>> readPixels(std::string filename)
{
    std::vector<std::vector<double> > numbers;
    std::ifstream infile (filename);
    std::string temp;
    
    //Discard first 3 lines in ds9 region file
    std::getline(infile, temp);
    std::getline(infile, temp);
    std::getline(infile, temp);
    
    while (std::getline(infile, temp))
    {
        std::vector<double> line;
        std::size_t klammer1 = temp.find("(");
        std::size_t comma = temp.find(",");
        std::size_t klammer2 = temp.find(")");
        
        int len1 = comma-klammer1;
        int len2 = klammer2-comma;
        
        std::string a = temp.substr(klammer1+1,len1-1);
        line.push_back(std::stod(a));
        a = temp.substr(comma+1,len2-1);
        line.push_back(std::stod(a));
        //std::cout << line.at(0) << " " << line.at(1) << std::endl;
        
        numbers.push_back(line);
    }
    return numbers;
}

std::vector<std::string> readLines(std::string filename)
{
    std::vector<std::string> lines;
    std::ifstream infile (filename);
    std::string temp;

    while (std::getline(infile, temp))
    {
        lines.push_back(temp);
    }
    return lines;
}


#endif
