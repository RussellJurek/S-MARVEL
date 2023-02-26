#include <iostream>
#include <cstdio>
#include <cstring>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <GaussFit_OMP.h>

extern "C" {

#include <cpgplot.h>
#include <tiffio.h>

}
#include<RJJ_ObjGen.h>
#include<RJJ_ObjGen_Plots.h>

using namespace std;

int main(int argc, char* argv[]){

  // define variables

  // a. generating the list of *.tif & *.tiff files to process
  stringstream tiff_file_list;
  string tiff_file_list_string;
  vector<string> tiff_files;
  int tf;
  unsigned int tiff_processed_count = 0;
  vector<string>::iterator tiff_file;
  FILE * tiff_file_system_call;
  char tiff_file_system_call_buffer[1000];

  // b. reading and writing *.tif & *.tiff files into memory
  TIFF * tiff_input, * tiff_output;
  uint32 x, y, size_x, size_y, npixels, tiff_row, * raster, tiff_bitdepth, tiff_spp, tiff_config, tiff_sample_format, tiff_colour_format;
  unsigned char * tiff_image_data;
  tdata_t scan_data;
  tsize_t tiff_linebytes, scanline_size;
  unsigned char * tiff_buffer;
  unsigned short int * scan_values, sample_val;

  // c. processing *.tif & *.tiff files
  string dummy1;
  stringstream dummy2;
  int mm_scale = 5, fft_size_x, fft_size_y, fft_start_x, fft_start_y;
  float * img_vals_1, * img_vals_2, img_min, img_max;
  double * img_vals_3, * img_vals_3_speq, fft_median, fft_sigma, fft_val, fft_angle, fft_max, fft_rel_clip;
  float fft_x_r2,fft_y_r2;
  int fft_x_min,fft_x_mid,fft_x_max,fft_y_min,fft_y_mid,fft_y_max;

  // d. generating histogram
  float hist_step, * x_vals, * y_vals, hist_max_val, hist_mean, hist_sigma;
  int Nhist = 100,hist_max_pos;
  
  // e. intensity thresholding
  float sigma_clip = 2.5, abs_clip;
  int * mask;

  // f. catalogue of sources
  vector<object_props *> detections;
  vector<int> obj_ids, check_obj_ids;
  int NOobj, obj_limit = 1000, obj_batch_limit = 100000, obj_batch, * xyz_order;
  int NO_check_obj_ids, NO_obj_ids, cat_mode;
  fstream outputfile;
  string output_code;
  size_t * data_metric;
  
  // g. variables used to scale velocity field postage stamps
  float cdelt3 = 1.0, crpix3 = 0.0, crval3 = 0.0, restfreq = 1.0;
  int ctype3 = 1;

  // h. gaussian fitting to positive side of pixel intensity histogram
  float * fit_params, ** fit_covar;
  int best_NOp;

  // i. calculate mean and standard deviation temp variables
  double val_sum, val2_sum;

  // j. dot thresholding variables
  int dot_x_min, dot_y_min, dot_v_min, dot_x_merge, dot_y_merge;
  float dot_rel_sum_max;

  // k. generate plots
  int plots_flag = 0;

  float cfd, cfd2, fft_cfd_25, fft_cfd_75;

  int i,j;
  float angle;

  // i. pgplot variables
  float tr[6] = {0.0,1.0,0.0,0.0,0.0,1.0};

  // create Gaussian fitting arrays
  fit_params = new float[7];
  fit_covar = new float * [3];
  for(i = 0; i < 3; ++i){ fit_covar[i] = new float[3]; }  

  // seed random number generator
  srand(time(NULL));

  // get the spatial scale of the mathematical morphology operation and the relative clip value of the fourier components
  if(argc == 8 || argc == 9){

    dummy2.str("");
    dummy2.clear();
    dummy2 << argv[1];
    dummy2 >> sigma_clip;

    dummy2.str("");
    dummy2.clear();
    dummy2 << argv[2];
    dummy2 >> dot_x_min;

    dummy2.str("");
    dummy2.clear();
    dummy2 << argv[3];
    dummy2 >> dot_y_min;
	
    dummy2.str("");
    dummy2.clear();
    dummy2 << argv[4];
    dummy2 >> dot_v_min;
	
    dummy2.str("");
    dummy2.clear();
    dummy2 << argv[5];
    dummy2 >> dot_rel_sum_max;
	
    dummy2.str("");
    dummy2.clear();
    dummy2 << argv[6];
    dummy2 >> dot_x_merge;
	
    dummy2.str("");
    dummy2.clear();
    dummy2 << argv[7];
    dummy2 >> dot_y_merge;
	
	if(argc == 9){

		dummy2.str("");
		dummy2.clear();
		dummy2 << argv[8];
		dummy2 >> plots_flag;

	}
	
  } else { 

	cout << "Insufficient number of arguments provided: FindDots <sigma_clip> <dot_x_min> <dot_x_max> <dot_v_min> <dot_rel_sum_max> <dot_x_merge> <dot_y_merge> " << endl;

  }

  // construct a list of all the *.tif & *.tiff files in the directory from which the program was called
  tiff_file_list.str("");
  tiff_file_list.clear();
  tiff_file_system_call = popen("ls *_MMremove.tif *_MMremove.tiff 2>/dev/null","r");
  if(tiff_file_system_call){
    while(!feof(tiff_file_system_call)){
      if(fgets(tiff_file_system_call_buffer,1000,tiff_file_system_call)){
		tiff_file_list << tiff_file_system_call_buffer;
      }
    }
  }
  pclose(tiff_file_system_call);

  while(tiff_file_list >> tiff_file_list_string){

    if((tiff_file_list_string.empty()) || (tiff_file_list_string == "") || (tiff_file_list_string[0] == '*')){ continue; }
    if((tiff_file_list_string.find("tif") < tiff_file_list_string.size()) || (tiff_file_list_string.find("tiff") < tiff_file_list_string.size())){

      tiff_files.push_back(tiff_file_list_string);
    
    }

  }

  // display the list of files found
  cout << "Found " << tiff_files.size() << " tif files to process in this directory." << endl;
  cout << "List of tif files: ";
  tiff_file = tiff_files.begin();
  while(tiff_file != tiff_files.end()){ cout << *tiff_file << " "; ++tiff_file; }
  cout << endl;

  // loop through the list of *.tif and *.tiff files and process them
  tiff_file = tiff_files.begin();
  while(tiff_file != tiff_files.end()){

    // 1. display the file being processed
    tiff_processed_count++;
    cout << "Processing file: " << *tiff_file << ", " << tiff_processed_count << " of " << tiff_files.size() <<  endl;

    // 2. open the input file
    tiff_input = TIFFOpen(tiff_file->c_str(),"r");

    // 3. create 2 buffers to store the image data
    // 3.a retrieve the image width
    TIFFGetField(tiff_input,TIFFTAG_IMAGEWIDTH,&size_x);
    // 3.b retrieve the image height
    TIFFGetField(tiff_input,TIFFTAG_IMAGELENGTH,&size_y);
    // 3.c calculate the number of pixels
    npixels = size_x * size_y;
	// 3.d retrieve the scanline size in Bytes
	scanline_size = TIFFScanlineSize(tiff_input);
	cout << "Image has a scanline size of " << scanline_size << " Bytes." << endl;
	// 3.e retrieve the planar configuration 
	TIFFGetField(tiff_input,TIFFTAG_PLANARCONFIG,&tiff_config);
	cout << "Image has planar configuration type, " << tiff_config << endl;
    // 3.f retrieve the bit-depth
    TIFFGetField(tiff_input,TIFFTAG_BITSPERSAMPLE,&tiff_bitdepth);
    cout << "Image has a bit-depth of " << tiff_bitdepth << endl;
    // 3.g retrieve the samples per pixel
    TIFFGetField(tiff_input,TIFFTAG_SAMPLESPERPIXEL,&tiff_spp);
    cout << "Image has " << tiff_spp << " samples per pixel." << endl;
    // 3.h retrieve the samples format
    TIFFGetField(tiff_input,TIFFTAG_SAMPLEFORMAT,&tiff_sample_format);
    cout << "Image pixels have a sample format of " << tiff_sample_format << endl;
	// 3.i retrieve the colouring format
    TIFFGetField(tiff_input,TIFFTAG_PHOTOMETRIC,&tiff_colour_format);
	cout << "Using photometric colur scheme, " << tiff_colour_format << endl;
    // 3.i create arrays 
    img_vals_1 = new float[npixels];
    img_vals_2 = new float[npixels];
    cout << "Created " << size_x << " x " << size_y << " = " << npixels << " element arrays to store 2 copies of image." << endl;
	
    // 4. read the image data into the first buffer
	if(tiff_bitdepth == 16){

		scan_data = _TIFFmalloc(scanline_size);	
		tiff_image_data = static_cast<unsigned char *>(scan_data);

		for(y = 0; y < size_y; ++y){

			TIFFReadScanline(tiff_input,scan_data,y);
			for(x = 0; x < size_x; ++x){ 

				sample_val = (tiff_image_data[2*x + 1] << 8) | (tiff_image_data[2*x]);
				img_vals_1[(y * size_x) + x] = static_cast<float>(sample_val);

			}

		}

		// free up memory
		_TIFFfree(scan_data);

	} else {
    
		// 4.a read in raw blue, red and green channels
		raster = (uint32*) _TIFFmalloc(npixels * sizeof(uint32));
		TIFFReadRGBAImage(tiff_input,size_x,size_y,raster,0);

		// 4.b combine RGB into a grey-scale image
		for(y = 0; y < size_y; y++){
		  
		  for(x = 0; x < size_x; x++){
		
			img_vals_1[(((size_y - 1 - y) * size_x) + x)] = 0.0;
			// process a gray-scale image
			img_vals_1[(((size_y - 1 - y) * size_x) + x)]+=(TIFFGetR(raster[((y * size_x) + x)]));
			img_vals_1[(((size_y - 1 - y) * size_x) + x)]+=(TIFFGetG(raster[((y * size_x) + x)]));
			img_vals_1[(((size_y - 1 - y) * size_x) + x)]+=(TIFFGetB(raster[((y * size_x) + x)]));
			img_vals_1[(((size_y - 1 - y) * size_x) + x)]/=3.0;

			// for(i = 0; i < size_x; i++)
		  }
		  
		  // for(j = 0; j < size_y; j++)
		}

		// free up memory
		_TIFFfree(raster);

	} 
    
    // 5. close the input file
    TIFFClose(tiff_input);
	
    // 6. create a histogram of the image values
    img_min = 0.0;
    img_max = 256.0;
    Nhist = 256;
    hist_step = (img_max - img_min)/static_cast<float>(Nhist);
    cout << "Creating histogram of " << Nhist << " values from " << img_min << " to " << img_max << " with step = " << hist_step << endl;
    img_vals_2 = new float[Nhist];
    for(i = 0; i < Nhist; ++i){ img_vals_2[i] = 0.0; }
    for(y = 0; y < size_y; ++y){
    
      for(x = 0; x < size_x; ++x){

	i = static_cast<int>(floorf((img_vals_1[(y * size_x) + x] - img_min) / hist_step));
	if(i < 0){ i = 0; }
	if(i >= Nhist){ i = Nhist - 1; }
	++img_vals_2[i];

      }
    
    }

    // 7. find the histogram maximum
    hist_max_val = -99.0;
    hist_max_pos = 0;
    for(i = 0; i < Nhist; ++i){

      if(img_vals_2[i] > hist_max_val){

	hist_max_val = img_vals_2[i];
	hist_max_pos = i;

      }
      
    }

    // 8. approximate the standard deviation using the FWHM
    j = hist_max_pos;
    for(i = hist_max_pos + 1; i < Nhist; ++i){

      if(img_vals_2[i] > (0.5*hist_max_val)){ j = i; }

    }
    hist_sigma = 2.0 * (static_cast<float>(j - hist_max_pos) + (img_vals_2[j] - (0.5 * hist_max_val))/(img_vals_2[j] - img_vals_2[j + 1])) / 2.355;

    // 9. convert hist_max_pos to the mean
    hist_mean = img_min + (hist_max_pos*hist_step);
    cout << "Initial estimate of histogram corresponds to gaussian amplitude = " << hist_max_val << ", mean = " << hist_mean << " and std.dev. = " << hist_sigma << endl;
    
	// 10. calculate the standard deviation and mean 
    for(y = 0; y < size_y; ++y){

      for(x = 0; x < size_x; ++x){

		val_sum+=img_vals_1[(y * size_x) + x];
		val2_sum+=(img_vals_1[(y * size_x) + x]*img_vals_1[(y * size_x) + x]);

      }

    }
	val2_sum/=static_cast<double>(size_x * size_y);
	val_sum/=static_cast<double>(size_x * size_y);
	hist_mean = val_sum;
	hist_sigma = val2_sum - (val_sum*val_sum);
	if(isinf(hist_sigma)){
		hist_sigma = (hist_sigma >= 0.0) ? 9E30 : -9E30;
	}
	hist_sigma = isnan(hist_sigma) ? 0.0 : hist_sigma;
	hist_sigma = (hist_sigma > 0.0) ? sqrt(hist_sigma) : 0.0;
	cout << "Brute force mean = " << hist_mean << " and brute force sigma = " << hist_sigma << endl;
    
    // 11. Apply an intensity threshold to the image to generate a binary mask
    abs_clip = hist_mean - (sigma_clip * hist_sigma);
    cout << "Using an absolute threshold of <= " << abs_clip << " for a relative thresold of <= " << sigma_clip << " std. deviations." << endl;
    mask = new int[size_x * size_y];
    j = 0;
    for(y = 0; y < size_y; ++y){

      for(x = 0; x < size_x; ++x){

		if(img_vals_1[(y * size_x) + x] <= abs_clip){

		  mask[(y * size_x) + x] = -1;
		  ++j;

		} else {

		  mask[(y * size_x) + x] = -99;

		}

      }

    }
    cout << "Flagged " << j << " of " << (size_x * size_y) << " pixels as being signal." << endl;

    // 12. Create a catalogue of objects from the mask

    // a. generate the output code
    i = tiff_file->find_first_of('.');
    output_code = tiff_file->substr(0,i);
    cout << "Using output_code = " << output_code << " to generate output files." << endl;

    // b. create initial arrays used for cataloguing and object construction
    InitObjGen(detections,NOobj,obj_limit,obj_ids,check_obj_ids,data_metric,xyz_order);
    xyz_order[0] = 1;
    xyz_order[1] = 2;
    xyz_order[2] = 3;

    // c. create metric for accessing this data chunk in x,y,z order
    CreateMetric(data_metric,xyz_order,size_x,size_y,1);
    cout << "data_metric = " << data_metric[0] << " " << data_metric[1] << " " << data_metric[2] << endl;

    // d. call function to create objects out of array of flagged voxels, and generate sparse representation of objects
    cout << "Creating objects . . . " << endl;
    NOobj = CreateObjects(img_vals_1,mask,size_x,size_y,1,0,0,0,dot_x_merge,dot_y_merge,0,dot_x_min,dot_y_min,0,dot_v_min,-9E30,dot_rel_sum_max*hist_mean,-1,0,detections,obj_ids,check_obj_ids,obj_limit,size_x,size_y,1,1,data_metric,xyz_order);

    // e. apply final size threshold and intensity threshold
    cout << "\nApplying final size threshold and additional intensity threshold . . . " << endl;
    ThresholdObjs(detections,NOobj,obj_limit,2,2,0,4,-9E30,9E30,4);
    
    // f. create catalogue
    outputfile.open((output_code + "_obj.cat").c_str(),ios::out);
    outputfile << "# " << j << " of " << (size_x * size_y) << " pixels are signal." << endl;
	outputfile << "# Using command line arguments: " << sigma_clip << " " << dot_x_min << " " << dot_y_min << " " << dot_v_min << " " << dot_rel_sum_max << " " << dot_x_merge << " " << dot_y_merge << endl;
    PrintCatalogueHeader(outputfile,0);
    i = CreateCatalogue(outputfile,detections,NOobj,obj_limit,-1);
    outputfile.close();

    // g. create global plots
    if(plots_flag == 1){ CreateObjPlots(i,1,output_code,xyz_order,size_x,size_y,1,detections,NOobj,obj_limit,ctype3,crpix3,crval3,cdelt3,restfreq); }

    // 13. free memory 
    if(img_vals_1){ delete [] img_vals_1; }
    if(img_vals_2){ delete [] img_vals_2; }
    if(x_vals){ delete [] x_vals; }
    if(y_vals){ delete [] y_vals; }
    if(mask){ delete [] mask; }
    FreeObjGen(detections,data_metric,xyz_order);

    // 14. advance the file iterator
    ++tiff_file;

    // while(tiff_file != tiff_files.end())
  }
  
  // free up memory
  delete [] fit_params;
  for(i = 0; i < 3; ++i){ delete [] fit_covar[i]; }  
  delete [] fit_covar;


  // int main()
}

