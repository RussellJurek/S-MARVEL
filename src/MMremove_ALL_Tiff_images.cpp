#include <iostream>
#include <cstdio>
#include <cstring>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <bitset>
#include <MathMorph_RJJ.h>
extern "C" {

#include <tiffio.h>

}

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
  int mm_scale = 5, fft_size_x, fft_size_y, fft_start_x, fft_start_y, xyz_order[3], data_metric[3];
  float * img_vals_1, * img_vals_2, img_min, img_max;
  unsigned short int new_value;
  double * img_vals_3, * img_vals_3_speq, fft_median, fft_sigma;

  float cfd, cfd2, fft_cfd_25, fft_cfd_75;

  int i;
  float angle;

  // d. pgplot variables
  float tr[6] = {0.0,1.0,0.0,0.0,0.0,1.0};

  // get the spatial scale of the mathematical morphology operation
  if(argc > 1){

    dummy2.str("");
    dummy2.clear();
    dummy2 << argv[1];
    dummy2 >> mm_scale;

  }
  cout << "Using radius of " << mm_scale << " pixels for circular iterative median smoothing kernel." << endl;

  // construct a list of all the *.tif & *.tiff files in the directory from which the program was called
  tiff_file_list.str("");
  tiff_file_list.clear();
  tiff_file_system_call = popen("ls *.tif *.tiff","r");
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
  tiff_processed_count = 0;
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
    
    // 6. copy the image data to the second buffer
    img_min = 9E30;
    img_max = -9E30;
    for(y = 0; y < size_y; ++y){

      for(x = 0; x < size_x; ++x){

		img_vals_2[(((size_y - 1 - y) * size_x) + x)] = img_vals_1[(((size_y - 1 - y) * size_x) + x)];
		if(img_vals_1[(((size_y - 1 - y) * size_x) + x)] > img_max){ img_max = img_vals_1[(((size_y - 1 - y) * size_x) + x)]; }
		if(img_vals_1[(((size_y - 1 - y) * size_x) + x)] < img_min){ img_min = img_vals_1[(((size_y - 1 - y) * size_x) + x)]; }

      }

    }

    // 10. apply iterative median smoothing to second image buffer
    cout << "Applying iterative median smoothing to image . . . " << endl;
    xyz_order[0] = 1;
    xyz_order[1] = 2;
    xyz_order[2] = 3;
    CreateMetric(data_metric,xyz_order,(int) size_x,(int) size_y,1);
    MedianSmooth(img_vals_2,(int) size_x,(int) size_y,1,1,0,1,0,1);
    MedianSmooth(img_vals_2,(int) size_x,(int) size_y,1,mm_scale,0,1,0,1);

    // 14. subtract the second image buffer from the first image buffer
    img_min = 9E30;
    img_max = -9E30;
    for(y = 0; y < size_y; ++y){

      for(x = 0; x < size_x; ++x){

		img_vals_1[((y * size_x) + x)] = img_vals_1[((y * size_x) + x)] - img_vals_2[((y * size_x) + x)];
		if(img_vals_1[((y * size_x) + x)] > img_max){ img_max = img_vals_1[((y * size_x) + x)]; }
		if(img_vals_1[((y * size_x) + x)] < img_min){ img_min = img_vals_1[((y * size_x) + x)]; }

      }

    }

    // 15.b open the third output file --- _MMremove.tiff
    tf = tiff_file->find_first_of(".");
    cout << "Creating tiff version of the third output file: " << (tiff_file->substr(0,tf) + "_MMremove.tiff") << endl;
	cout << "Scaling image to cover range: " << img_min << " --> " << img_max << endl;
    tiff_output = TIFFOpen((tiff_file->substr(0,tf) + "_MMremove.tiff").c_str(),"w");
    TIFFSetField(tiff_output,TIFFTAG_IMAGEWIDTH,size_x);
    TIFFSetField(tiff_output,TIFFTAG_IMAGELENGTH,size_y);
    TIFFSetField(tiff_output,TIFFTAG_BITSPERSAMPLE,tiff_bitdepth);
    TIFFSetField(tiff_output,TIFFTAG_ORIENTATION,ORIENTATION_TOPLEFT);
    TIFFSetField(tiff_output,TIFFTAG_PLANARCONFIG,PLANARCONFIG_CONTIG);
	TIFFSetField(tiff_output,TIFFTAG_PHOTOMETRIC,tiff_colour_format);

    // 16.b plot the subtraction result in the third output file
    if(tiff_bitdepth == 16){

		TIFFSetField(tiff_output,TIFFTAG_SAMPLESPERPIXEL,1);

		tiff_image_data = new unsigned char[size_x * size_y * 2];
		for(y = 0; y < size_y; ++y){
		  
		  for(x = 0; x < size_x; ++x){

			// generate a gray-scale image
			new_value = ((unsigned short int)floorf(0.5 + (16384.0 * ((img_vals_1[((y * size_x) + x)] - img_min)/(img_max - img_min)))));
			tiff_buffer = ((unsigned char *)&new_value);
			tiff_image_data[(y * size_x * 2) + (x * 2)] = *tiff_buffer;
			tiff_buffer++;
			tiff_image_data[(y * size_x * 2) + (x * 2) + 1] = *tiff_buffer;

		  }

		}

		tiff_linebytes = 2 * size_x;
		TIFFSetField(tiff_output,TIFFTAG_ROWSPERSTRIP,TIFFDefaultStripSize(tiff_output,size_y));

    } else {

		TIFFSetField(tiff_output,TIFFTAG_SAMPLESPERPIXEL,3);

	    tiff_image_data = new unsigned char[size_x * size_y * 3];
	    for(y = 0; y < size_y; ++y){
	      
	      for(x = 0; x < size_x; ++x){

			// generate a gray-scale image
			tiff_image_data[(y * size_x * 3) + (3 * x)] = ((unsigned int)floorf(0.5 + (255.0 * ((img_vals_1[((y * size_x) + x)] - img_min)/(img_max - img_min)))));
			tiff_image_data[(y * size_x * 3) + (3 * x) + 1] = ((unsigned int)floorf(0.5 + (255.0 * ((img_vals_1[((y * size_x) + x)] - img_min)/(img_max - img_min)))));
			tiff_image_data[(y * size_x * 3) + (3 * x) + 2] = ((unsigned int)floorf(0.5 + (255.0 * ((img_vals_1[((y * size_x) + x)] - img_min)/(img_max - img_min)))));
	       
	      }

	    }

		tiff_linebytes = (tiff_bitdepth / 8);
		tiff_linebytes*= (3 * size_x);
		TIFFSetField(tiff_output,TIFFTAG_ROWSPERSTRIP,TIFFDefaultStripSize(tiff_output,size_x*3));

    }
 
    tiff_buffer = NULL;
    tiff_buffer = (unsigned char *)_TIFFmalloc(tiff_linebytes);

    for(tiff_row = 0; tiff_row < size_y; ++tiff_row){
      memcpy(tiff_buffer, &tiff_image_data[tiff_row*tiff_linebytes], tiff_linebytes); 
      if (TIFFWriteScanline(tiff_output,tiff_buffer,tiff_row,0) < 0){ break; }
    }

    // 17.b close the third output file
    delete [] tiff_image_data;
    if(tiff_buffer){ _TIFFfree(tiff_buffer); }
    TIFFClose(tiff_output);

    // 18. release img_vals_2 as it's no longer used
    delete [] img_vals_2;
    img_vals_2 = NULL;

    // 26. free the memory used to store the image vals
    if(img_vals_1){ delete [] img_vals_1; }
    if(img_vals_2){ delete [] img_vals_2; }

    // 27. advance the file iterator
    ++tiff_file;

    // while(tiff_file != tiff_files.end())
  }

  // int main()
}

