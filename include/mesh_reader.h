#include <cmath>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <sstream>
#include <vector>


//! Mesh reader class
//! \brief Mesh class for creating nodes and elements from GMSH input
class MeshReader {
 public:

  //! \brief Read keywords
  //! \param[in] file Input file stream
  //! \param[in] keyword Search for keyword
  void read_keyword(std::ifstream& file, const std::string& keyword) {

    bool read_status = false;
    std::string line;
    file.clear();
    file.seekg(0, std::ios::beg);
    while (std::getline(file, line)) {
      if (line != keyword) {
        if (line.find(keyword) != std::string::npos) {
          break;
        };
      } else {
        std::cout << "Read keyword: " << keyword << " successfully\n";
        read_status = true;
        break;
      }
    }
    if (!read_status) {
      std::cerr << "Cannot find keyword: " << keyword << '\n';
    }
  }

  //! \brief Read ids and coordinates of vertices
  //! \param[in] file Input file stream object of msh file
  void read_vertices(std::ifstream& file) {
    read_keyword(file, "$Nodes");

    std::string line;
    std::getline(file, line);
    std::istringstream istream(line);

    // Total number of vertices
    unsigned nvertices;
    istream >> nvertices;
    std::cout << "Total number of vertices = " << nvertices << '\n';

    // Vertex id and coordinates
    unsigned vid = std::numeric_limits<unsigned>::max();
    std::vector<double> coordinates(3, 0);

    // Iterate through all vertices in the mesh file
    for (unsigned i = 0; i < nvertices;) {
      std::getline(file, line);
      std::istringstream istream(line);
      if (line.find('#') == std::string::npos && line != "") {
        // Initialise ids and coordinates
        vid = std::numeric_limits<unsigned>::max();

        // Read ids and coordinates
        istream >> vid;
        for (unsigned j = 0; j < coordinates.size(); ++j)
          istream >> coordinates.at(j);

        // Add vertex coordinates and id to a map
        vertices_[vid] = coordinates;

        // Increament number of vertex on successful read
        ++i;
      } else {
        std::cerr << "Invalid entry for node: " << line << '\n';
      }
    }
  }

  //! \brief Read ids and coordinates of vertices
  //! \param[in] file Input file stream object of msh file
  void read_interior_line(std::ifstream& file) {

    std::string line;
    std::getline(file, line);
    std::istringstream istream(line);

    // Total number of vertices
    unsigned nlines;
    istream >> nlines;
    std::cout << "Total number of interior lines = " << nlines << '\n';

    // line ID
    unsigned lid = std::numeric_limits<unsigned>::max();
    //! Element type
    unsigned element_type = std::numeric_limits<unsigned>::max();
    //! Number of tags
    unsigned ntags = std::numeric_limits<unsigned>::max();
    unsigned tag = std::numeric_limits<unsigned>::max();
    //! Object id
    unsigned object_id = std::numeric_limits<unsigned>::max();
    //! Node id
    unsigned node_id = std::numeric_limits<unsigned>::max();

    for (unsigned i = 0; i < nlines; ++i) {
      std::getline(file, line);
      std::istringstream istream(line);
      if (line.find('#') == std::string::npos && line != "") {
        // Read ids and element type
        istream >> lid;
        istream >> element_type;
        istream >> ntags;
        istream >> object_id;
        // Read element tags
        for (unsigned j = 0; j < ntags - 1; ++j) {
          istream >> tag;
        }

        for (unsigned id = 0; id < 2; ++id) {
          istream >> node_id;
          iline_vert_ids_.push_back(node_id);
        }

      } else {
        //std::cerr << "Invalid entry for node: " << line << '\n';
      }
    }
    // Remove duplicates
    std::sort(iline_vert_ids_.begin(), iline_vert_ids_.end());
    iline_vert_ids_.erase(std::unique(iline_vert_ids_.begin(), iline_vert_ids_.end()), iline_vert_ids_.end());

/*
std::cout<<"vert id on interior line: ";
for(auto id : iline_vert_ids_) std::cout<<id<<", ";
std::cout<<std::endl;
std::cout<<"number of vertices: "<<vertices_.size()<<std::endl;
*/

  }

  //! \brief Check if a msh file exists
  //! \param[in] filename Mesh file name
  void read_mesh(const std::string& mesh_filename, const std::string& line_filename) {
    std::ifstream mesh_file;
    mesh_file.open(mesh_filename.c_str(), std::ios::in);
    if (!mesh_file.is_open())
      throw std::runtime_error("Specified GMSH file does not exist");
    if (mesh_file.good()) {
      read_vertices(mesh_file);
    }
    mesh_file.close();

    std::ifstream line_file;
    line_file.open(line_filename.c_str(), std::ios::in);
    if (!line_file.is_open())
      throw std::runtime_error("Specified txt file does not exist");
    if (line_file.good()) {
      read_interior_line(line_file);
    }
    line_file.close();
  }

  //! Return coords of vertices on interior line
  std::vector<std::vector<double>> interior_line_vertices() {
    std::vector<std::vector<double>> interior_line_vertices;
    for(auto vert_id : iline_vert_ids_) interior_line_vertices.push_back(vertices_.at(vert_id));
    return interior_line_vertices;
  }

 private:
  //! A map of vertex id and coordinates
  std::map<unsigned, std::vector<double>> vertices_;
  //! vector of indices of vector on the interior line
  std::vector<unsigned> iline_vert_ids_;


};

