#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>
extern "C" {

#include <fitsio.h>

}

using namespace std;

// c++ MathMorphClose_Cube_v1.cpp -o MathMorphClose_Cube_v1 -O2 -lcfitsio

void MathMorphOpen(float * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag);

void MathMorphClose(float * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag);
 
void MinSmooth(float * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int axis, int round_flag);

void MaxSmooth(float * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int axis, int round_flag);

int Populate_temp_array(float * temp_vals, float * proc_vals, int i, int j, int k, int size_x, int size_y, int size_z, int scale);

float Temp_array_min(float * temp_vals, float * proc_vals, int i, int j, int k, int size_x, int size_y, int size_z, int used);

float Temp_array_max(float * temp_vals, float * proc_vals, int i, int j, int k, int size_x, int size_y, int size_z, int used);

int main(int argc, char* argv[]){

  // define variables and arrays
  fitsfile * inputfile, * outputfile, * outputfile_diff;
  string dummy1, inputfile_name, outputfile_name, outputfile_code;
  stringstream dummy2;
  int sradius, fradius, round_flag, status;
  int file_xsize, file_ysize, file_zsize, x, y, z, i;
  int chunk_xsize, chunk_ysize, chunk_zsize, c, cx, cy, NOc, NOcx, NOcy;
  long fits_read_start[4], fits_read_finish[4], fits_read_inc[4] = {1,1,1,1};
  int chunk_size = 1024 * 1024 * 1024 / 4, chunk_x_overlap, chunk_y_overlap;
  int chunk_x_start, chunk_y_start, chunk_x_finish, chunk_y_finish;
  float * data_vals, * data_vals_diff;

  // get input parameters from command line
  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[1];
  dummy2 >> inputfile_name;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[2];
  dummy2 >> outputfile_code;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[3];
  dummy2 >> sradius;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[4];
  dummy2 >> fradius;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[5];
  dummy2 >> round_flag;

  // open input file
  status = 0;
  cout << "Attempting to open input file: " << inputfile_name << endl;
  fits_open_file(&inputfile,inputfile_name.c_str(),READONLY,&status);
  cout << "File opened, status = " << status << endl;
  if(status > 0){ fits_report_error(stderr,status); return 1; }

  // get size of input file
  fits_read_key(inputfile,TINT,"NAXIS1",&file_xsize,NULL,&status);
  fits_read_key(inputfile,TINT,"NAXIS2",&file_ysize,NULL,&status);
  fits_read_key(inputfile,TINT,"NAXIS3",&file_zsize,NULL,&status);
  if(status > 0){ 
    
    fits_report_error(stderr,status); 
    return 1; 
    
  }
  cout << file_xsize << " x " << file_ysize << " spatial elements, " << file_zsize << " frequency channels in data cube." << endl;

  // check that the smoothing scales of the 2nd iteration aren't too large, adjust if necessary
  if(((2*sradius) + 1) > file_xsize){ 

    sradius = (int) floorf(((float) file_xsize - 1.0)/2.0); 
    cout << "WARNING!!! sradius too large for x-axis, setting to: " << sradius << endl;

  }
  if(((2*sradius) + 1) > file_ysize){ 

    sradius = (int) floorf(((float) file_ysize - 1.0)/2.0);
    cout << "WARNING!!! sradius too large for y-axis, setting to: " << sradius << endl;
  
  }
  if(((2*fradius) + 1) > file_zsize){

    fradius = (int) floorf(((float) file_zsize - 1.0)/2.0);
    cout << "WARNING!!! fradius too large for z-axis, setting to: " << fradius << endl;

  }

  // check if this is a 2-D image, and adjust fradius accordingly
  if(file_zsize == 1){

    fradius = 0;
    cout << "WARNING!!! This is a 2-D image. The fradius is being set to " << fradius << " in response." << endl;

  }

  // calculate overlap between chunks
  chunk_x_overlap = chunk_y_overlap = (2 * sradius);
  if(chunk_x_overlap < 2){ chunk_x_overlap = chunk_y_overlap = 2; }
  cout << "Using an overlap of " << chunk_x_overlap << " along RA & Dec for chunks." << endl;

  // calculate the number of `chunks' that the input file 
  // needs to be split into and display results
  chunk_xsize = chunk_ysize = (int) floorf((sqrtf((floorf((float) (chunk_size / file_zsize))))));
  if(chunk_xsize > file_xsize){ chunk_xsize = file_xsize; }
  if(chunk_ysize > file_ysize){ chunk_ysize = file_ysize; }
  chunk_zsize = file_zsize;
  
  if(chunk_xsize < file_xsize){ NOcx = ((int) ceilf((((float) file_xsize / (float) (chunk_xsize - chunk_x_overlap))))); } else { NOcx = 1; }
  if(chunk_ysize < file_ysize){ NOcy = ((int) ceilf((((float) file_ysize / (float) (chunk_ysize - chunk_y_overlap))))); } else { NOcy = 1; }
  NOc = NOcx * NOcy;

  cout << "Processing file using overlapping " << (chunk_size * 4 / (1024 * 1024)) << "MB `chunks'. This is most efficiently achieved using " << NOc << " chunk/s of size " << chunk_xsize << " (x) by " << chunk_ysize << " (y) by " << chunk_zsize << " (freq.) voxels. This amounts to " << NOcx << " (x) by " << NOcy << " (y) chunks." << endl;

  // create arrays used to iteratively median smooth a chunk
  data_vals = new float[(chunk_xsize * chunk_ysize * chunk_zsize)];
  data_vals_diff = new float[(chunk_xsize * chunk_ysize * chunk_zsize)];

  // open output file
  status = 0;
  outputfile_name = outputfile_code + "_close.fits";
  cout << "Attempting to create output file: " << outputfile_name << " using " << inputfile_name << " as a template." << endl;
  outputfile_name = "!" + outputfile_name;
  fits_create_template(&outputfile,outputfile_name.c_str(),inputfile_name.c_str(),&status);
  cout << "File opened, status = " << status << endl;
  if(status > 0){ fits_report_error(stderr,status); return 1; }

  // check dimensions of output file
  fits_read_key(inputfile,TINT,"NAXIS1",&c,NULL,&status);
  if(c != file_xsize){ 

    cout << "WARNING!!! Size of input & output files RA dimension DOESN'T match." << endl;
    delete [] data_vals;
    delete [] data_vals_diff;
   return 1;

  } else {
    
    cout << "Size of input & output files RA dimension matches." << endl;

  }
  fits_read_key(inputfile,TINT,"NAXIS2",&c,NULL,&status);
  if(c != file_ysize){ 

    cout << "WARNING!!! Size of input & output files Dec dimension DOESN'T match." << endl;
    delete [] data_vals;
    delete [] data_vals_diff;
    return 1;

  } else {
    
    cout << "Size of input & output files Dec dimension matches." << endl;

  }
  fits_read_key(inputfile,TINT,"NAXIS3",&c,NULL,&status);
  if(c != file_zsize){ 

    cout << "WARNING!!! Size of input & output files frequency dimension DOESN'T match." << endl;
    delete [] data_vals;
    delete [] data_vals_diff;
   return 1;

  } else {
    
    cout << "Size of input & output files frequency dimension matches." << endl;

  }

  // open diff output file
  status = 0;
  outputfile_name = outputfile_code + "_cdiff.fits";
  cout << "Attempting to create output file: " << outputfile_name << " using " << inputfile_name << " as a template." << endl;
  outputfile_name = "!" + outputfile_name;
  fits_create_template(&outputfile_diff,outputfile_name.c_str(),inputfile_name.c_str(),&status);
  cout << "File opened, status = " << status << endl;
  if(status > 0){ fits_report_error(stderr,status); return 1; }

  // check dimensions of diff output file
  fits_read_key(inputfile,TINT,"NAXIS1",&c,NULL,&status);
  if(c != file_xsize){ 

    cout << "WARNING!!! Size of input & output files RA dimension DOESN'T match." << endl;
    delete [] data_vals;
    delete [] data_vals_diff;
    return 1;

  } else {
    
    cout << "Size of input & output files RA dimension matches." << endl;

  }
  fits_read_key(inputfile,TINT,"NAXIS2",&c,NULL,&status);
  if(c != file_ysize){ 

    cout << "WARNING!!! Size of input & output files Dec dimension DOESN'T match." << endl;
    delete [] data_vals;
    delete [] data_vals_diff;
    return 1;

  } else {
    
    cout << "Size of input & output files Dec dimension matches." << endl;

  }
  fits_read_key(inputfile,TINT,"NAXIS3",&c,NULL,&status);
  if(c != file_zsize){ 

    cout << "WARNING!!! Size of input & output files frequency dimension DOESN'T match." << endl;
    delete [] data_vals;
    delete [] data_vals_diff;
    return 1;

  } else {
    
    cout << "Size of input & output files frequency dimension matches." << endl;

  }

  // process the chunks
  for(cy = 0; cy < NOcy; cy++){

    for(cx = 0; cx < NOcx; cx++){

      cout << "\nProcessing chunk " << (1 + cx + (cy * NOcx)) << " of " << NOc << " . . . " << endl;
      
      // set parameters for region to be read from memory
      fits_read_start[0] = 1 + (cx * (chunk_xsize - chunk_x_overlap));
      if(fits_read_start[0] < 1){ fits_read_start[0] = 1; }
      fits_read_finish[0] = fits_read_start[0] + chunk_xsize - 1;
      if(fits_read_finish[0] > file_xsize){ fits_read_finish[0] = file_xsize; }

      fits_read_start[1] = 1 + (cy * (chunk_ysize - chunk_y_overlap));
      if(fits_read_start[1] < 1){ fits_read_start[1] = 1; }
      fits_read_finish[1] = fits_read_start[1] + chunk_ysize - 1;
      if(fits_read_finish[1] > file_ysize){ fits_read_finish[1] = file_ysize; }
      
      fits_read_start[2] = 1;
      fits_read_finish[2] = chunk_zsize;
      
      fits_read_start[3] = 1;
      fits_read_finish[3] = 1;
      
      cout << "Covers datacube region: " << fits_read_start[0] << " <= x <= " << fits_read_finish[0] << ", " << fits_read_start[1] << " <= y <= " << fits_read_finish[1] << ", " << fits_read_start[2] << " <= z <= " << fits_read_finish[2] << endl; 

      // read the current chunk into memory
      cout << "Reading chunk of input file into memory . . . " << endl;
      fits_read_subset(inputfile,TFLOAT,fits_read_start,fits_read_finish,fits_read_inc,NULL,data_vals,NULL,&status);
      if(status > 0){ fits_report_error(stderr,status); return 1; }
  
      // copy input data into data_vals_diff
      cout << "Initialising data_vals_diff array . . . " << endl;
      for(i = 0; i < (chunk_xsize * chunk_ysize * chunk_zsize); i++){ data_vals_diff[i] = data_vals[i]; }

      // apply mathematical morphology opening to input data, to filter out signal on small scales
      cout << "Applying mathematical morphology opening . . . " << endl;
      MathMorphClose(data_vals,(fits_read_finish[0] - fits_read_start[0] + 1),(fits_read_finish[1] - fits_read_start[1] + 1),(fits_read_finish[2] - fits_read_start[2] + 1),sradius,fradius,round_flag);

      // update data_vals_diff array for results of mathematical morphology
      cout << "Updating data_vals_diff array for results of mathematical morphology . . . " << endl;
      for(i = 0; i < (chunk_xsize * chunk_ysize * chunk_zsize); i++){ data_vals_diff[i]-=data_vals[i]; }

      // set boundaries of region that hasn't previously been smoothed and/or will be smoothed better by the next datacube

      // x boundaries
      if(cx == 0){

	chunk_x_start = 0;

      } else {

	if(sradius > 1){

	  chunk_x_start = sradius;

	} else {

	  chunk_x_start = 1;

	}

      }
      /*if(cx == (NOcx - 1)){

	chunk_x_finish = fits_read_finish[0];

      } else {

	chunk_x_finish = fits_read_finish[0] - sradius;

	}*/

      // y boundaries
      if(cy == 0){

	chunk_y_start = 0;

      } else {

	if(sradius > 1){
	
	  chunk_y_start = sradius;

	} else {

	  chunk_y_start = 1;

	}

      }
      /*if(cy == (NOcy - 1)){

	chunk_y_finish = fits_read_finish[1];

      } else {

	chunk_y_finish = fits_read_finish[1] - sradius;

	}*/

      // shuffle values in data_vals array to account for overlap
      cout << "Shuffling values in data_vals array to account for overlap between chunks . . . " << endl;
      if((cx != 0) && (cy != 0)){

	for(z = 0; z < chunk_zsize; z++){

	  for(y = 0; y < (fits_read_finish[1] - fits_read_start[1] + 1 - sradius); y++){

	    for(x = 0; x < (fits_read_finish[0] - fits_read_start[0] + 1 - sradius); x++){

	      data_vals[((z * (fits_read_finish[0] - fits_read_start[0] + 1) * (fits_read_finish[1] - fits_read_start[1] + 1)) + (y * (fits_read_finish[0] - fits_read_start[0] + 1)) + x)] = data_vals[((z * (fits_read_finish[0] - fits_read_start[0] + 1) * (fits_read_finish[1] - fits_read_start[1] + 1)) + ((y + sradius) * (fits_read_finish[0] - fits_read_start[0] + 1)) + x + sradius)];

	    }

	  }

	}

      } else if((cx == 0) && (cy != 0)){

	for(z = 0; z < chunk_zsize; z++){
	  
	  for(y = 0; y < (fits_read_finish[0] - fits_read_start[0] + 1 - sradius); y++){
	    
	    for(x = 0; x < (fits_read_finish[0] - fits_read_start[0] + 1); x++){
	      
	      data_vals[((z * (fits_read_finish[0] - fits_read_start[0] + 1) * (fits_read_finish[1] - fits_read_start[1] + 1)) + (y * (fits_read_finish[0] - fits_read_start[0] + 1)) + x)] = data_vals[((z * (fits_read_finish[0] - fits_read_start[0] + 1) * (fits_read_finish[1] - fits_read_start[1] + 1)) + ((y + sradius) * (fits_read_finish[0] - fits_read_start[0] + 1)) + x)];
	      
	    }
	    
	  }
	  
	}

      } else if((cx != 0) && (cy == 0)){

	for(z = 0; z < chunk_zsize; z++){

	  for(y = 0; y < (fits_read_finish[1] - fits_read_start[1] + 1); y++){

	    for(x = 0; x < (fits_read_finish[0] - fits_read_start[0] + 1 - sradius); x++){

	      data_vals[((z * (fits_read_finish[0] - fits_read_start[0] + 1) * (fits_read_finish[1] - fits_read_start[1] + 1)) + (y * (fits_read_finish[0] - fits_read_start[0] + 1)) + x)] = data_vals[((z * (fits_read_finish[0] - fits_read_start[0] + 1) * (fits_read_finish[1] - fits_read_start[1] + 1)) + (y * (fits_read_finish[0] - fits_read_start[0] + 1)) + x + sradius)];

	    }

	  }

	}

      }

      // update fits_read arrays for boundary of non-overlapping region of chunk
      fits_read_start[0] = fits_read_start[0] + chunk_x_start;
      fits_read_start[1] = fits_read_start[1] + chunk_y_start;
      /*fits_read_finish[0] = chunk_x_finish;
	fits_read_finish[1] = chunk_y_finish;*/

      // write the iteratively median smoothed chunk to the output file
      cout << "Writing mathematical morphology filtered chunk covering " << fits_read_start[0] << " <= x <= " << fits_read_finish[0] << ", " << fits_read_start[1] << " <= y <= " << fits_read_finish[1] << ", " << fits_read_start[2] << " <= z <= " << fits_read_finish[2] << " to output file . . . " << endl;
      fits_write_subset(outputfile,TFLOAT,fits_read_start,fits_read_finish,data_vals,&status);
      if(status > 0){ fits_report_error(stderr,status); return 1; } else { cout << "done." << endl; }

      // write the residual chunk to the diff output file
      cout << "Writing residual chunk covering " << fits_read_start[0] << " <= x <= " << fits_read_finish[0] << ", " << fits_read_start[1] << " <= y <= " << fits_read_finish[1] << ", " << fits_read_start[2] << " <= z <= " << fits_read_finish[2] << " to output file . . . " << endl;
      fits_write_subset(outputfile_diff,TFLOAT,fits_read_start,fits_read_finish,data_vals_diff,&status);
      if(status > 0){ fits_report_error(stderr,status); return 1; } else { cout << "done." << endl; }

      // for(cx = 0; cx < NOcx; cx++)
    }

    // for(cy = 0; cy < NOcy; cy++)
  }

  // close the input file
  cout << "Closing fits input file . . . "; cout.flush();
  fits_close_file(inputfile,&status);
  cout << "done, status = " << status << endl;

  // close the output file
  cout << "Closing fits output file . . . "; cout.flush();
  fits_close_file(outputfile,&status);
  cout << "done, status = " << status << endl;

  // close the iff output file
  cout << "Closing fits diff output file . . . "; cout.flush();
  fits_close_file(outputfile_diff,&status);
  cout << "done, status = " << status << endl;

  // free memory
  delete [] data_vals;
  delete [] data_vals_diff;

}

void MathMorphOpen(float * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag){

  void MinSmooth(float * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int axis, int round_flag);
  void MaxSmooth(float * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int axis, int round_flag);
  int axis;

  // determine the longest axis
  axis = 1;
  if((size_y >= size_x) && (size_y > size_z)){ axis = 2; }
  if((size_z >= size_x) && (size_z >= size_y)){ axis = 3; }

  cout << "Smoothing along axis: " << axis << endl;

  // carry out the first smoothing iteration
  cout << "Carrying out 1st smoothing iteration . . . " << endl;
  MinSmooth(data_vals,size_x,size_y,size_z,sradius,fradius,axis,round_flag);

  // carry out the second smoothing iteration
  cout << "Carrying out 2nd smoothing iteration . . . " << endl;
  MaxSmooth(data_vals,size_x,size_y,size_z,sradius,fradius,axis,round_flag);

}

void MathMorphClose(float * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag){

  void MinSmooth(float * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int axis, int round_flag);
  void MaxSmooth(float * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int axis, int round_flag);
  int axis;

  // determine the longest axis
  axis = 1;
  if((size_y >= size_x) && (size_y > size_z)){ axis = 2; }
  if((size_z >= size_x) && (size_z >= size_y)){ axis = 3; }

  cout << "Smoothing along axis: " << axis << endl;

  // carry out the first smoothing iteration
  cout << "Carrying out 1st smoothing iteration . . . " << endl;
  MaxSmooth(data_vals,size_x,size_y,size_z,sradius,fradius,axis,round_flag);

  // carry out the second smoothing iteration
  cout << "Carrying out 2nd smoothing iteration . . . " << endl;
  MinSmooth(data_vals,size_x,size_y,size_z,sradius,fradius,axis,round_flag);

}

void MinSmooth(float * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int axis, int round_flag){

  int Populate_temp_array(float * temp_vals, float * data_vals, int i, int j, int k, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag);
  float Temp_array_min(float * temp_vals, float * data_vals, int i, int j, int k, int size_x, int size_y, int size_z, int used);
  float * buffer_vals, * temp_vals, sep, sep_max, progress;
  int i,j,k,s,bdepth,used,NOvals;
  
  // create temp array to sort
  if(round_flag <= 0){
    
    NOvals = ((2*sradius) + 1)*((2*sradius) + 1)*((2*fradius) + 1);
    temp_vals = new float[NOvals];
    
  } else {
    
    NOvals = 0;
    
    for(k = (-1 * fradius); k < ((-1 * fradius) + sradius); k++){
      
      // calculate projected radius for this channel
      sep_max = acos(((float)(k + fradius - sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; j++){
	
	for(i = (-1 * sradius); i <= sradius; i++){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; i++)
	}
	
	// for(j = (-1 * sradius); j <= sradius; j++)
      }
      
      // for(k = (-1 * fradius); k <= fradius; k++)
    }
    
    for(k = (-1 * fradius) + sradius; k <= (fradius - sradius); k++){
      
      // calculate projected radius for this channel
      sep_max = (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; j++){
	
	for(i = (-1 * sradius); i <= sradius; i++){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; i++)
	}
	
	// for(j = (-1 * sradius); j <= sradius; j++)
      }
      
      // for(k = (-1 * fradius); k <= fradius; k++)
    }
    
    for(k = (fradius - sradius + 1); k <= fradius; k++){
      
      // calculate projected radius for this channel
      sep_max = acos(((float)(k - fradius + sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; j++){
	
	for(i = (-1 * sradius); i <= sradius; i++){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; i++)
	}
	
	// for(j = (-1 * sradius); j <= sradius; j++)
      }
      
      // for(k = (-1 * fradius); k <= fradius; k++)
    }

    temp_vals = new float[NOvals];
    
  }
  cout << "Using temporary array containing " << NOvals << " values." << endl;
  
  // carry out smoothing by shuffling data along longest axis
  // using the buffer array
  switch (axis)
    {
      
    case 1:
      // longest axis is ra
      
      // calculate buffer depth
      bdepth = sradius + 1;

      // create buffer array
      buffer_vals = new float [(size_y * size_z * bdepth)];
      
      // fill initial buffer values
      for(i = 0; i < (bdepth - 1); i++){

	for(k = 0; k < size_z; k++){
	  
	  for(j = 0; j < size_y; j++){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; j++)
	  }
	
	  // for(k = 0; k < size_z; k++)
	}   

	// for(i = 0; i < (bdepth - 1); i++)
      }
	      
      // buffer and shuffle median values
      cout << "0 | |:| | : | |:| | 100% complete" << endl;
      progress = 0.0;
      for(i = (bdepth - 1); i < size_x; i++){
	
	for(k = 0; k < size_z; k++){
	  
	  for(j = 0; j < size_y; j++){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_y) + (j*bdepth) + bdepth - 1)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; j++)
	  }
	  
	  //for(k = 0; k < size_z; k++)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; k++){
	  
	  for(j = 0; j < size_y; j++){

	    data_vals[((k*size_x*size_y) + (j*size_x) + i - bdepth + 1)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth))];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s + 1)];
	      
	    // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; j++)
	  }
	
	  // for(k = 0; k < size_z; k++)
	}

	while(progress <= (((float) (i + 1)) / ((float) size_x))){ cout << "*"; cout.flush(); progress+=0.05; }

	// for(i = (bdepth - 1); i < size_x; i++)
      }
      cout << "* done." << endl;
      
      // write final set of buffer values to array
      for(i = 0; (i < (bdepth - 1)); i++){
	
	for(k = 0; k < size_z; k++){
	  
	  for(j = 0; j < size_y; j++){
	    	      
	    data_vals[((k*size_x*size_y) + (j*size_x) + size_x - bdepth + 1 + i)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)];
	      	    
	    // for(j = 0; j < size_y; j++)
	  }
	  
	  // for(k = 0; k < size_z; k++)
	}

	// (i = 0; (i < (bdepth - 1)); i++)
      }

      // end of case statement
      break;
      
    case 2:
      // longest axis is dec
      
      // calculate buffer depth
      bdepth = sradius + 1;
      
      // create buffer array
      buffer_vals = new float [(size_x * size_z * bdepth)];
      
      // fill initial buffer values	
      for(j = 0; j < (bdepth - 1); j++){

	for(k = 0; k < size_z; k++){
	  
	  for(i = 0; i < size_x; i++){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(i = 0; i < size_x; i++)
	  }
	
	  //for(k = 0; k < size_z; k++)
	}   
	
	// for(j = 0; j < (bdepth - 1); j++)
      }
  
      // buffer and shuffle median values
      cout << "0 | |:| | : | |:| | 100% complete" << endl;
      progress = 0.0;
      for(j = (bdepth - 1); j < size_y; j++){
	
	for(k = 0; k < size_z; k++){
	  
	  for(i = 0; i < size_x; i++){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_x) + ((bdepth - 1)*size_x) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(i = 0; i < size_x; i++)
	  }
	  
	  //for(k = 0; k < size_z; k++)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; k++){
	  
	  for(i = 0; i < size_x; i++){
	    
	    data_vals[((k*size_x*size_y) + ((j - bdepth + 1)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_x) + (s*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + ((s + 1)*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(i = 0; i < size_x; i++)
	  }
	  
	  // for(k = 0; k < size_z; k++)
	}
	
	while(progress <= (((float) (j + 1)) / ((float) size_y))){ cout << "*"; cout.flush(); progress+=0.05; }

	// for(j = (bdepth - 1); j < size_y; j++)
      }
      cout << "* done." << endl;
      
      // write final set of buffer values to array
      for(j = 0; (j < (bdepth - 1)); j++){
	
	for(k = 0; k < size_z; k++){
	  
	  for(i = 0; i < size_x; i++){
	    	      
	    data_vals[((k*size_x*size_y) + ((size_y - bdepth + 1 + j)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)];
	      	    
	    // for(i = 0; i < size_x; i++)
	  }
	  
	  // for(k = 0; k < size_z; k++)
	}

	// for(j = 0; (j < (bdepth - 1)); j++)
      }
      
      // end of case statement
      break;
      
    case 3:
      // longest axis is frequency
      
      // calculate buffer depth
      bdepth = fradius + 1;

      // create buffer array
      buffer_vals = new float [(size_x * size_y * bdepth)];
      
      // fill initial buffer values
      for(k = 0; k < (bdepth - 1); k++){
	
	for(i = 0; i < size_x; i++){
	  
	  for(j = 0; j < size_y; j++){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*size_x*size_y) + (j*size_x) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; j++)
	  }
	  
	  //for(i = 0; i < size_x; i++)
	}   
	
	// for(k = 0; k < (bdepth - 1); k++)
      }  
      
      // buffer and shuffle median values
      cout << "0 | |:| | : | |:| | 100% complete" << endl;
      progress = 0.0;
      for(k = (bdepth - 1); k < size_z; k++){
	
	for(i = 0; i < size_x; i++){
	  
	  for(j = 0; j < size_y; j++){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[(((bdepth - 1)*size_x*size_y) + (j*size_x) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; j++)
	  }
	  
	  //for(i = 0; i < size_x; i++)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(i = 0; i < size_x; i++){
	  
	  for(j = 0; j < size_y; j++){
	    
	    data_vals[(((k - bdepth + 1)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((j*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((s*size_x*size_y) + (j*size_x) + i)] = buffer_vals[(((s + 1)*size_x*size_y) + (j*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; j++)
	  }
	  
	  // for(i = 0; i < size_x; i++)
	}

	while(progress <= (((float) (k + 1)) / ((float) size_z))){ cout << "*"; cout.flush(); progress+=0.05; }
	
	// for(k = (bdepth - 1); k < size_z; k++)
      }  
      cout << "* done." << endl;

      // write final set of buffer values to array
      for(k = 0; (k < (bdepth - 1)); k++){
	
	for(i = 0; i < size_x; i++){
	  
	  for(j = 0; j < size_y; j++){
	    	      
	    data_vals[(((size_z - bdepth + 1 + k)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((k*size_x*size_y) + (j*size_x) + i)];
	      	    
	    // for(j = 0; j < size_y; j++)
	  }
	  
	  // for(i = 0; i < size_x; i++)
	}

	// for(k = 0; (k < (bdepth - 1)); k++)
      }
      
      // end of case statement
      break;
           
      // switch(axis) . . . case
    }
  
  // free up memory 
  delete [] temp_vals;
  delete [] buffer_vals;
  
}

void MaxSmooth(float * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int axis, int round_flag){

  int Populate_temp_array(float * temp_vals, float * data_vals, int i, int j, int k, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag);
  float Temp_array_max(float * temp_vals, float * data_vals, int i, int j, int k, int size_x, int size_y, int size_z, int used);
  float * buffer_vals, * temp_vals, sep, sep_max, progress;
  int i,j,k,s,bdepth,used,NOvals;
  
  // create temp array to sort
  if(round_flag <= 0){
    
    NOvals = ((2*sradius) + 1)*((2*sradius) + 1)*((2*fradius) + 1);
    temp_vals = new float[NOvals];
    
  } else {
    
    NOvals = 0;
    
    for(k = (-1 * fradius); k < ((-1 * fradius) + sradius); k++){
      
      // calculate projected radius for this channel
      sep_max = acos(((float)(k + fradius - sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; j++){
	
	for(i = (-1 * sradius); i <= sradius; i++){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; i++)
	}
	
	// for(j = (-1 * sradius); j <= sradius; j++)
      }
      
      // for(k = (-1 * fradius); k <= fradius; k++)
    }
    
    for(k = (-1 * fradius) + sradius; k <= (fradius - sradius); k++){
      
      // calculate projected radius for this channel
      sep_max = (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; j++){
	
	for(i = (-1 * sradius); i <= sradius; i++){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; i++)
	}
	
	// for(j = (-1 * sradius); j <= sradius; j++)
      }
      
      // for(k = (-1 * fradius); k <= fradius; k++)
    }
    
    for(k = (fradius - sradius + 1); k <= fradius; k++){
      
      // calculate projected radius for this channel
      sep_max = acos(((float)(k - fradius + sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; j++){
	
	for(i = (-1 * sradius); i <= sradius; i++){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; i++)
	}
	
	// for(j = (-1 * sradius); j <= sradius; j++)
      }
      
      // for(k = (-1 * fradius); k <= fradius; k++)
    }

    temp_vals = new float[NOvals];
    
  }
  cout << "Using temporary array containing " << NOvals << " values." << endl;
  
  // carry out smoothing by shuffling data along longest axis
  // using the buffer array
  switch (axis)
    {
      
    case 1:
      // longest axis is ra
      
      // calculate buffer depth
      bdepth = sradius + 1;

      // create buffer array
      buffer_vals = new float [(size_y * size_z * bdepth)];
      
      // fill initial buffer values
      for(i = 0; i < (bdepth - 1); i++){

	for(k = 0; k < size_z; k++){
	  
	  for(j = 0; j < size_y; j++){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; j++)
	  }
	
	  // for(k = 0; k < size_z; k++)
	}   

	// for(i = 0; i < (bdepth - 1); i++)
      }
	      
      // buffer and shuffle median values
      cout << "0 | |:| | : | |:| | 100% complete" << endl;
      progress = 0.0;
      for(i = (bdepth - 1); i < size_x; i++){
	
	for(k = 0; k < size_z; k++){
	  
	  for(j = 0; j < size_y; j++){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_y) + (j*bdepth) + bdepth - 1)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; j++)
	  }
	  
	  //for(k = 0; k < size_z; k++)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; k++){
	  
	  for(j = 0; j < size_y; j++){

	    data_vals[((k*size_x*size_y) + (j*size_x) + i - bdepth + 1)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth))];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s + 1)];
	      
	    // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; j++)
	  }
	
	  // for(k = 0; k < size_z; k++)
	}

	while(progress <= (((float) (i + 1)) / ((float) size_x))){ cout << "*"; cout.flush(); progress+=0.05; }

	// for(i = (bdepth - 1); i < size_x; i++)
      }
      cout << "* done." << endl;
      
      // write final set of buffer values to array
      for(i = 0; (i < (bdepth - 1)); i++){
	
	for(k = 0; k < size_z; k++){
	  
	  for(j = 0; j < size_y; j++){
	    	      
	    data_vals[((k*size_x*size_y) + (j*size_x) + size_x - bdepth + 1 + i)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)];
	      	    
	    // for(j = 0; j < size_y; j++)
	  }
	  
	  // for(k = 0; k < size_z; k++)
	}

	// (i = 0; (i < (bdepth - 1)); i++)
      }

      // end of case statement
      break;
      
    case 2:
      // longest axis is dec
      
      // calculate buffer depth
      bdepth = sradius + 1;
      
      // create buffer array
      buffer_vals = new float [(size_x * size_z * bdepth)];
      
      // fill initial buffer values	
      for(j = 0; j < (bdepth - 1); j++){

	for(k = 0; k < size_z; k++){
	  
	  for(i = 0; i < size_x; i++){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(i = 0; i < size_x; i++)
	  }
	
	  //for(k = 0; k < size_z; k++)
	}   
	
	// for(j = 0; j < (bdepth - 1); j++)
      }
  
      // buffer and shuffle median values
      cout << "0 | |:| | : | |:| | 100% complete" << endl;
      progress = 0.0;
      for(j = (bdepth - 1); j < size_y; j++){
	
	for(k = 0; k < size_z; k++){
	  
	  for(i = 0; i < size_x; i++){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_x) + ((bdepth - 1)*size_x) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(i = 0; i < size_x; i++)
	  }
	  
	  //for(k = 0; k < size_z; k++)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; k++){
	  
	  for(i = 0; i < size_x; i++){
	    
	    data_vals[((k*size_x*size_y) + ((j - bdepth + 1)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_x) + (s*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + ((s + 1)*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(i = 0; i < size_x; i++)
	  }
	  
	  // for(k = 0; k < size_z; k++)
	}
	
	while(progress <= (((float) (j + 1)) / ((float) size_y))){ cout << "*"; cout.flush(); progress+=0.05; }

	// for(j = (bdepth - 1); j < size_y; j++)
      }
      cout << "* done." << endl;
      
      // write final set of buffer values to array
      for(j = 0; (j < (bdepth - 1)); j++){
	
	for(k = 0; k < size_z; k++){
	  
	  for(i = 0; i < size_x; i++){
	    	      
	    data_vals[((k*size_x*size_y) + ((size_y - bdepth + 1 + j)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)];
	      	    
	    // for(i = 0; i < size_x; i++)
	  }
	  
	  // for(k = 0; k < size_z; k++)
	}

	// for(j = 0; (j < (bdepth - 1)); j++)
      }
      
      // end of case statement
      break;
      
    case 3:
      // longest axis is frequency
      
      // calculate buffer depth
      bdepth = fradius + 1;

      // create buffer array
      buffer_vals = new float [(size_x * size_y * bdepth)];
      
      // fill initial buffer values
      for(k = 0; k < (bdepth - 1); k++){
	
	for(i = 0; i < size_x; i++){
	  
	  for(j = 0; j < size_y; j++){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*size_x*size_y) + (j*size_x) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; j++)
	  }
	  
	  //for(i = 0; i < size_x; i++)
	}   
	
	// for(k = 0; k < (bdepth - 1); k++)
      }  
      
      // buffer and shuffle median values
      cout << "0 | |:| | : | |:| | 100% complete" << endl;
      progress = 0.0;
      for(k = (bdepth - 1); k < size_z; k++){
	
	for(i = 0; i < size_x; i++){
	  
	  for(j = 0; j < size_y; j++){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[(((bdepth - 1)*size_x*size_y) + (j*size_x) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; j++)
	  }
	  
	  //for(i = 0; i < size_x; i++)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(i = 0; i < size_x; i++){
	  
	  for(j = 0; j < size_y; j++){
	    
	    data_vals[(((k - bdepth + 1)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((j*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((s*size_x*size_y) + (j*size_x) + i)] = buffer_vals[(((s + 1)*size_x*size_y) + (j*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; j++)
	  }
	  
	  // for(i = 0; i < size_x; i++)
	}

	while(progress <= (((float) (k + 1)) / ((float) size_z))){ cout << "*"; cout.flush(); progress+=0.05; }
	
	// for(k = (bdepth - 1); k < size_z; k++)
      }  
      cout << "* done." << endl;

      // write final set of buffer values to array
      for(k = 0; (k < (bdepth - 1)); k++){
	
	for(i = 0; i < size_x; i++){
	  
	  for(j = 0; j < size_y; j++){
	    	      
	    data_vals[(((size_z - bdepth + 1 + k)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((k*size_x*size_y) + (j*size_x) + i)];
	      	    
	    // for(j = 0; j < size_y; j++)
	  }
	  
	  // for(i = 0; i < size_x; i++)
	}

	// for(k = 0; (k < (bdepth - 1)); k++)
      }
      
      // end of case statement
      break;
           
      // switch(axis) . . . case
    }
  
  // free up memory 
  delete [] temp_vals;
  delete [] buffer_vals;
  
}

int Populate_temp_array(float * temp_vals, float * data_vals, int i, int j, int k, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag){

  int used = 0,x,y,z;
  int x_start, y_start, z_start;
  float sep, sep_max;

  // generate boundaries of smoothing box
  x_start = i - sradius;
  if(x_start < 0){ x_start = 0; }
  y_start = j - sradius;
  if(y_start < 0){ y_start = 0; }
  z_start = k - fradius;
  if(z_start < 0){ z_start = 0; }

  // populate smoothing box depending on the value of round_flag
  if(round_flag <= 0){

    // rectangular smoothing box

    for(z = z_start; ((z < size_z) && (z <= (k + fradius))); z++){

      for(y = y_start; ((y < size_y) && (y <= (j + sradius))); y++){

	for(x = x_start; ((x < size_x) && (x <= (i + sradius))); x++){

	  temp_vals[used] = data_vals[((z * size_x * size_y) + (y * size_x) + x)];
	  used++;
	  
	}

      }

    }

  } else {

    // rounded cylindrical smoothing box

    for(z = z_start; ((z < (z_start + sradius)) && (z < size_z)); z++){

      // calculate projected radius for this channel
      sep_max = acos(((float)(z - z_start - sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;

      for(y = y_start; ((y <= (j + sradius)) && (y < size_y)); y++){

	for(x = x_start; ((x <= (i + sradius)) && (x < size_x)); x++){
    
	  sep = sqrtf(((float)((x - i) * (x - i)) + (float)((y - j) * (y - j))));
	  if(sep <= sep_max){ 

	    temp_vals[used] = data_vals[((z * size_x * size_y) + (y * size_x) + x)];
	    used++; 

	  }

	  // for(i = (-1 * sradius); i <= sradius; i++)
	}
	  
	// for(j = (-1 * sradius); j <= sradius; j++)
      }

      // for(k = (-1 * fradius); k <= fradius; k++)
    }

    for(z = z_start + sradius; ((z <= (k + fradius - sradius)) && (z < size_z)); z++){

      // calculate projected radius for this channel
      sep_max = (float) sradius;

      for(y = y_start; ((y <= (j + sradius)) && (y < size_y)); y++){

	for(x = x_start; ((x <= (i + sradius)) && (x < size_x)); x++){
    
	  sep = sqrtf(((float)((x - i) * (x - i)) + (float)((y - j) * (y - j))));
	  if(sep <= sep_max){ 

	    temp_vals[used] = data_vals[((z * size_x * size_y) + (y * size_x) + x)];
	    used++; 
	    
	  }

	  // for(i = (-1 * sradius); i <= sradius; i++)
	}
	  
	// for(j = (-1 * sradius); j <= sradius; j++)
      }

      // for(k = (-1 * fradius); k <= fradius; k++)
    }

    for(z = k + fradius - sradius + 1; ((z <= (k + fradius)) && (z < size_z)); z++){

      // calculate projected radius for this channel
      sep_max = acos(((float)(z - k - fradius + sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;

      for(y = y_start; ((y <= (j + sradius)) && (y < size_y)); y++){

	for(x = x_start; ((x <= (i + sradius)) && (x < size_x)); x++){
    
	  sep = sqrtf(((float)((x - i) * (x - i)) + (float)((y - j) * (y - j))));
	  if(sep <= sep_max){

	    temp_vals[used] = data_vals[((z * size_x * size_y) + (y * size_x) + x)];
	    used++; 
	    
	  }

	  // for(i = (-1 * sradius); i <= sradius; i++)
	}
	  
	// for(j = (-1 * sradius); j <= sradius; j++)
      }

      // for(k = (-1 * fradius); k <= fradius; k++)
    }

  }

  return used;

}

float Temp_array_min(float * temp_vals, float * data_vals, int i, int j, int k, int size_x, int size_y, int size_z, int used){

  float min;
  int u;
  
  min = data_vals[((k*size_x*size_y) + (j*size_x) + i)];

  // find minimum value
  for(u = 0; u < used; u++){
    
    if(temp_vals[u] <= min){ min = temp_vals[u]; }
    
  } 

  return min;

}


float Temp_array_max(float * temp_vals, float * data_vals, int i, int j, int k, int size_x, int size_y, int size_z, int used){

  float max;
  int u;
  
  max = data_vals[((k*size_x*size_y) + (j*size_x) + i)];

  // find minimum value
  for(u = 0; u < used; u++){
    
    if(temp_vals[u] >= max){ max = temp_vals[u]; }
    
  } 

  return max;

}




