#include "global.h"
class FileWriter {
 public:
  int output(Array3D<zone> grid, int time, const char *filename);
  void write_general_file(ofstream &gout, const string &str_input_filename) const;
};
