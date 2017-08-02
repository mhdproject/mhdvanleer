#include "global.h"
class FileWriter {
 public:
  int output(Array3D<zone> grid, Array3D<zone> fx, Array3D<zone> fy, int time, const char *filename);
  void write_general_file(ofstream &gout, const string &str_input_filename) const;
};
