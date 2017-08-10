#include "output.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_HDF5
#include "hdf5.h"
#endif

#define RANK   2

int
FileWriter::output(Array3D<zone> grid, int time, const char *filename) {

#ifdef USE_HDF5
  hid_t file, dataset;        /* file and dataset handles */
  hid_t datatype, dataspace;    /* handles */
  hsize_t dimsf[2];        /* dataset dimensions */
  herr_t status;
  double data[nx][ny];
  //array_2d data(boost::extents[nx][ny]);
  const char *names[] =
      {"Density", "Velx", "Vely", "Velz", "Energy", "Bx", "By", "Bz"};
  int ll = 0;
  stringstream hdf5_stream_filename;
  string hdf5_filename;
#endif

  ofstream fout;
  ofstream gout;
  double gammam1 = gammag - 1;
  double rl, ri;
  double px;
  double py;
  double pz;
  double pressure;
  double bx;
  double by;
  double bz;
  double bsquared;
  double et, ul, vl, wl, ke, al;
  int ii = 0;
  int jj = 0;
  int kk = 0;
  char outputdir[50] = "./";
//      char            filename[50] = "out_2d_";
  stringstream s;
  stringstream stream_filename;
  stringstream stream_temp_b;
  string str_file_tag;
  string str_output_filename;
  string str_input_filename;

  s.clear();
  s.width(5);
  s.fill('0');
  s << time;
  s >> str_file_tag;
  stream_filename.clear();
  stream_filename << outputdir << filename << str_file_tag;
  stream_filename >> str_input_filename;

#ifdef USE_HDF5
  hdf5_stream_filename << outputdir << "hdf5_" << filename << str_file_tag <<
                       ".h5";
  hdf5_stream_filename >> hdf5_filename;
  file =
      H5Fcreate(hdf5_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                H5P_DEFAULT);

  for (ll = 0; ll < ne; ll++) {
    dimsf[0] = static_cast<hsize_t>(nx);
    dimsf[1] = (hsize_t) ny;
    dataspace = H5Screate_simple(RANK, dimsf, nullptr);
    /*
     * Define datatype for the data in the file.
     * We will store little endian DOUBLE numbers.
     */
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status = H5Tset_order(datatype, H5T_ORDER_LE);
    assert(status == 0);

    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    dataset = H5Dcreate(file, names[ll], datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for (jj = 0; jj < ny; jj++) {
      for (ii = 0; ii < nx; ii++) {
        data[ii][jj] = grid[ii][jj][kk].array[ll];
      }
    }
    /*
     * Write the data to the dataset using default transfer properties.
     */
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, data);
    assert(status == 0);

    /*
     * Close/release resources.
     */
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
  }

  dimsf[0] = (hsize_t) nx;
  dimsf[1] = (hsize_t) ny;
  dataspace = H5Screate_simple(RANK, dimsf, nullptr);
  /*
   * Define datatype for the data in the file.
   * We will store little endian DOUBLE numbers.
   */
  datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
  status = H5Tset_order(datatype, H5T_ORDER_LE);
  assert(status == 0);
  /*
   * Create a new dataset within the file using defined dataspace and
   * datatype and default dataset creation properties.
   */
  dataset = H5Dcreate(file, "Pressure", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  for (jj = 0; jj < ny; jj++) {
    for (ii = 0; ii < nx; ii++) {

      rl = grid[ii][jj][kk]_MASS;
      px = grid[ii][jj][kk]_MOMX;
      py = grid[ii][jj][kk]_MOMY;
      pz = grid[ii][jj][kk]_MOMZ;
      et = grid[ii][jj][kk]_ENER;
      bx = grid[ii][jj][kk]_B_X;
      by = grid[ii][jj][kk]_B_Y;
      bz = grid[ii][jj][kk]_B_Z;
      ri = 1.0 / rl;
      ul = px * ri;
      vl = py * ri;
      wl = pz * ri;
      ke = 0.5 * rl * (ul * ul + vl * vl + wl * wl);
      bsquared = bx * bx + by * by + bz * bz;
      pressure = et - ke - 0.5 * bsquared;
      pressure = pressure * gammam1;

      data[ii][jj] = pressure;
    }
  }
  /*
   * Write the data to the dataset using default transfer properties.
   */
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, data);
  assert(status == 0);

  /*
   * Close/release resources.
   */
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);

  H5Fclose(file);
#endif /* HDF5 or not  */

  fout.open(str_input_filename.c_str());
  if (fout.fail()) {
    cerr << "unable to open file " << endl;
  }

  // Determine Div B
  Array2D<double> divb(nx, ny);
  double bx1, bx2, by1, by2;
  for (ii = 1; ii < nx - 2; ii++) {
    for (jj = 1; jj < ny - 2; jj++) {

      bx1 = (grid[ii][jj][kk]_B_X + grid[ii - 1][jj][kk]_B_X);
      bx2 = (grid[ii + 1][jj][kk]_B_X + grid[ii][jj][kk]_B_X);
      by1 = (grid[ii][jj][kk]_B_Y + grid[ii][jj - 1][kk]_B_Y);
      by2 = (grid[ii][jj + 1][kk]_B_Y + grid[ii][jj][kk]_B_Y);
      //     bz1 = ( grid[ii  ][jj  ][kk  ]_B_Z + grid[ii  ][jj  ][kk-1]_B_Z );
      //    bz2 = ( grid[ii  ][jj  ][kk+1]_B_Z + grid[ii  ][jj  ][kk  ]_B_Z );
      //divb = (1/delta_x)*(bx2- bx1 + by2 -by1 +bz2 -bz1);
      divb[ii][jj] = (0.5 / delta_x) * (bx2 - bx1 + by2 - by1);

    }
  }
  for (ii = 0; ii < nx; ii++) {
    for (jj = 0; jj < ny; jj++) {
      rl = grid[ii][jj][kk]_MASS;
      px = grid[ii][jj][kk]_MOMX;
      py = grid[ii][jj][kk]_MOMY;
      pz = grid[ii][jj][kk]_MOMZ;
      et = grid[ii][jj][kk]_ENER;
      bx = grid[ii][jj][kk]_B_X;
      by = grid[ii][jj][kk]_B_Y;
      bz = grid[ii][jj][kk]_B_Z;
      ri = 1.0 / rl;
      ul = px * ri;
      vl = py * ri;
      wl = pz * ri;
      ke = 0.5 * rl * (ul * ul + vl * vl + wl * wl);
      bsquared = bx * bx + by * by + bz * bz;
      pressure = et - ke - 0.5 * bsquared;
      pressure = pressure * gammam1;
      al = sqrt(gammag * pressure * ri);

#ifdef DEBUG_BC
      if (ii == 2 && jj == 2 && px != 0)
        {
          cout << px << endl;
          cout << ul << endl;
          cout << "wtf?" << endl;
        }
#endif /* DEBUG_BC */

      fout
          << setiosflags(ios::scientific)
          << " " << (rl)
          << " " << ul
          << " " << vl
          << " " << wl
          << " " << et
          << " " << bx
          << " " << by
          << " " << bz
          << " " << (pressure)
          << " " << (al)
          << " " << divb[ii][jj]
          << endl;
    }

#ifdef TWODIM
    fout << endl;
#endif /* TWODIM */
  }
  fout.close();

  write_general_file(gout, str_input_filename);

  return 0;
}
void FileWriter::write_general_file(ofstream &gout, const string &str_input_filename) const {
  gout.open("gm.general");
  gout << "file = /home/gmurphy/mhdvanleer-0.0.1/" << str_input_filename <<
       endl;
  gout << "grid = " << nx << " x " << ny << endl;
  gout << "format = ascii" << endl;
  gout << "interleaving = field" << endl;
  gout << "majority = row" << endl;
  gout << "field = V_sound, E_tot, Rho, Vel_X, Vel_Y, Pressure" << endl;
  gout << "structure = scalar, scalar, scalar, scalar, scalar, scalar" <<
       endl;
  gout << "type = double, double, double, double, double, double" << endl;
  gout <<
       "dependency = positions, positions, positions, positions, positions, positions"
       << endl;
  gout << "positions = regular, regular, 0, 1, 0, 1" << endl;
  gout << "" << endl;
  gout << "end" << endl;
  gout.close();
}

#ifdef USE_HDF5

#endif /* HDF5 or not  */
