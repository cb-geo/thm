#include "csv.h"
#include <vector>
void get_parameter(std::string filename, std::vector<double>& coord,std::vector<double>& data)
{

    double v1,v2;
    io::CSVReader<2> in(filename);
    while(in.read_row(v1, v2))
    {
      coord.push_back(v1);
      data.push_back(v2);
    }
      
}
