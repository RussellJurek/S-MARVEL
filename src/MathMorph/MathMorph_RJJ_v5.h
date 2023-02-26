#include<iostream>
#include<math.h>
#include<omp.h>

#ifndef RJJ_NonLinFilter
#define RJJ_NonLinFilter

using namespace std;

// specify the non linear filters

// a. median
template <typename data_type>
data_type median(data_type * image, size_t * test_region, unsigned int NOt){

  if((NOt % 2) == 0){

    return (0.5 * (image[test_region[(NOt/2) - 1]] + image[test_region[NOt/2]]));
    
  } else {

    return image[test_region[(NOt - 1)/2]];

  }

}

// b. minimum
template <typename data_type>
data_type minimum(data_type * image, size_t * test_region, unsigned int NOt){ return image[test_region[0]]; }

// c. maximum
template <typename data_type>
data_type maximum(data_type * image, size_t * test_region, unsigned int NOt){ return image[test_region[NOt - 1]]; }

// d. specific cumulative frequency distribution value

// e. inter-quartile range


// the heapsort function
template <class T_data, class T_index>
void HeapSort_index(int n, T_data ra[], T_index indices[]){

  int i,ir,j,l;
  T_data rra;
  int rraindex;

  if (n<2) return;
  l = (n >> 1) + 1;
  ir = n;
  for(;;){

    if (l > 1) {
      rra = ra[indices[--l]];
      rraindex = indices[l];
    } else {
      rra=ra[indices[ir]];
      rraindex = indices[ir];
      indices[ir] = indices[1];
      if (--ir == 1) {
	indices[1] = rraindex;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[indices[j]] < ra[indices[j+1]]) j++;
      if (rra < ra[indices[j]]) {
	indices[i] = indices[j];
	i = j;
	j <<= 1;
      } else break;
    }
    indices[i] = rraindex;
  }
}

// specify the non linear filter functions

// a. 1D: spectrum

// b. 2D: image
template <typename data_type>
void NonLinFilter_2D(data_type * image, int NOx, int NOy, bool elliptical, int r_x, int r_y, unsigned short int filter_mode, unsigned short int vb_flag) {

  unsigned short int curr_test_set, next_test_set;
  unsigned int t,NOt,t2,t3,NOt3,os,NOos,NOb;
  size_t ** test_region, d, NOd, * map_p_to_t, buffer_write, buffer_read;
  ssize_t i,j,i2,j2,j_min, * offsets, trailing_pos, leading_pos;
  double r,r_req,prog;
  data_type * buffer;
  data_type (* filter)(data_type *, size_t *, unsigned int);

  // extra variables used to implement parallelisation
  data_type * image_orig;
  int y_min,y_max,tid,NOtid;

  // NOTES:
  // I) elliptical is true for an elliptical mask and false for a rectangular mask
  // II) vb_flag = 0, 1 & 2 for no, minimal and maximal status messages
  
  // 0. Assign the filter function.
  switch(filter_mode){
  case 1:
    filter = &median;
    if(vb_flag > 1){ cout << "Applying median filter." << endl; }
    break;
  case 2:
    filter = &minimum;
    if(vb_flag > 1){ cout << "Applying minimum filter." << endl; }
    break;
  case 3:
    filter = &maximum;
    if(vb_flag > 1){ cout << "Applying maximum filter." << endl; }
    break;
  default:
    filter = &median;
    if(vb_flag > 1){ cout << "ERROR! No matching filter found for filter_mode value. Applying median filter by DEFAULT." << endl; }
    break;
  }

  // 1. Calculate the total size of the image
  NOd = NOx * NOy;

  // 2. Create a copy of the input image to be read from by the various threads, while the original is modified.
  image_orig = new data_type[NOd];
  for(d = 0; d < NOd; ++d){ image_orig[d] = image[d]; }

  // 3. Count the number of pixels that make up the mask --> NOt.
  if(elliptical){

    if(vb_flag > 1){ cout << "Counting the number of pixels in an elliptical mask . . . " << endl; }
    NOt = 0;
    for(j = 0; j < (1 + 2*r_y); ++j){

      for(i = 0; i < (1 + 2*r_x); ++i){
	
	if((((static_cast<double>(i) - static_cast<double>(r_x))*(static_cast<double>(i) - static_cast<double>(r_x))/(static_cast<double>(r_x) * static_cast<double>(r_x))) + 
	    ((static_cast<double>(j) - static_cast<double>(r_y))*(static_cast<double>(j) - static_cast<double>(r_y))/(static_cast<double>(r_y) * static_cast<double>(r_y)))) <= 1.0){

	  ++NOt;

	}
	
      }

    }
    
  } else {

    if(vb_flag > 1){ cout << "Calculating the number of pixels in a rectangular mask . . . " << endl; }
    NOt = (1 + 2*r_x) * (1 + 2*r_y);

  }
  if(vb_flag > 1){ cout << "Mask contains " << NOt << " pixels." << endl; }

  // 4. Count the number of pixels in the leading/trailing edge of the mask --> NOos.
  NOos = 1 + (2*r_y);
  if(vb_flag > 1){ cout << "Calculated that mask leading/trailing edge contains " << NOos << " pixels." << endl; }

  // 5. Calculate the buffer size.
  NOb = NOx * (1 + r_y);
  if(vb_flag > 1){ cout << "Calculated buffer size of " << NOb << "." << endl; }

  // 6. Create an array to hold NOos offsets corresponding to the leading/trailing mask edge.
  offsets = new ssize_t[2*NOos];
  if(vb_flag > 1){ cout << "Allocated offsets array containing " << (2*NOos) << " elements for " << NOos << " pixels." << endl; }
      
  // 7. Populate the offsets array.
  if(vb_flag > 1){ cout << "Populating the offsets array . . . " << endl; }
  NOos = 0;
  if(elliptical){    
    
    // for each j line, find the last pixel that belongs to the mask and store the correspondng offset
    if(vb_flag > 1){ cout << "Finding edge of elliptical mask . . . " << endl; }
    for(j = -1 * r_y; j <= r_y; ++j){
      
      for(i = 0; i <= r_x; ++i){
	
	if(((static_cast<double>(i)*static_cast<double>(i)/(static_cast<double>(r_x) * static_cast<double>(r_x))) + 
	    (static_cast<double>(j)*static_cast<double>(j)/(static_cast<double>(r_y) * static_cast<double>(r_y)))) <= 1.0){
	  
	  offsets[2*NOos] = j;
	    offsets[(2*NOos)+1] = i;
	    
	} else {
	  
	  ++NOos;
	  break;
	  
	}
	
      }
      
      if(i > r_x){ ++NOos; }
      
    }
    
  } else {
    
      // for each j line, store the offset of the last pixel in the mask
    if(vb_flag > 1){ cout << "Finding edge of rectangular mask . . . " << endl; }
    for(j = -r_y; j <= r_y; ++j){
      
      offsets[2*NOos] = j;
      offsets[(2*NOos)+1] = r_x;
      ++NOos;
      
    }
    
  }
  if(vb_flag > 1){ cout << "Done. Stored " << NOos << " offsets." << endl; }
  
  // 8. Retrieve the number of threads available and specify this many to be used
  NOtid = omp_get_max_threads();
  if(vb_flag > 1){ cout << "Using " << NOtid << " threads." << endl; }
  omp_set_num_threads(NOtid);

  // debugging
  /*omp_set_num_threads(2);
  NOtid = 2;
  cout << "Using " << NOtid << " threads for debugging." << endl;*/

  // 9. Start of parallel region
  #pragma omp parallel default(shared) firstprivate(NOt) private(test_region,map_p_to_t,buffer,leading_pos,trailing_pos,i,i2,j,j_min,j2,os,t,t2,t3,NOt3,y_min,y_max,tid,curr_test_set,next_test_set,buffer_read,buffer_write,prog)
  {
    
    // 10. Retrive the thread id
    tid = omp_get_thread_num();

    // 11. Generate the image boundaries for this thread
    y_min = tid * static_cast<int>(ceil(static_cast<double>(NOy) / static_cast<double>(NOtid)));
    y_max = (tid + 1) * static_cast<int>(ceil(static_cast<double>(NOy) / static_cast<double>(NOtid)));
    if(y_max > NOy){ y_max = NOy; }
    
    /*#pragma omp critical
    {
      cout << "thread " << tid << " covers " << y_min << " <= j < " << y_max << endl;
    }
    #pragma omp barrier*/

    // 12. Create an array that holds the mapping of relative positions to test region locations
    map_p_to_t = new size_t[(1 + 2*r_x)*(1 + 2*r_y)];
    for(i = 0; i < ((1 + 2*r_x)*(1 + 2*r_y)); ++i){ map_p_to_t[i] = 0; }
    #pragma omp master
    {
      if(vb_flag > 1){ cout << "Allocated map_p_to_t array containing " << ((1 + 2*r_x)*(1 + 2*r_y)) << " pixels." << endl; }
    }
    
    // debugging
    /*#pragma omp barrier
    #pragma omp master 
    cout << "a" << endl;
    #pragma omp barrier*/

    // 13. Create an array that can hold NOt index values.
    test_region = new size_t * [3];
    test_region[0] = new size_t [NOt];
    test_region[1] = new size_t [NOt];
    test_region[2] = new size_t [1 + (2*r_y)];
    #pragma omp master
    {
      if(vb_flag > 1){ cout << "Allocated test_region array containing " << NOt << " pixels." << endl; }
    }
       
    // debugging
    /*#pragma omp barrier
    #pragma omp master 
    cout << "b" << endl;
    #pragma omp barrier*/

    // 14. Create the buffer array.
    buffer = new data_type[NOb];
    #pragma omp master
    {
      if(vb_flag > 1){ cout << "Allocated buffer array." << endl; }
    }

    // debugging
    /*#pragma omp barrier
    #pragma omp master 
    cout << "c" << endl;
    #pragma omp barrier*/

    // 15. Populate the buffer with its initial values.
    #pragma omp master
    {
      if(vb_flag > 0){ 
	cout << "Populating buffer with its initial values." << endl; 
	cout << "000.00%" << flush;
      }
    }
    curr_test_set = 0;
    next_test_set = 1;
    prog = 0.0;

    // debugging
    /*#pragma omp barrier
    #pragma omp master 
    cout << "d" << endl;
    #pragma omp barrier*/

    // debugging
    /*#pragma omp critical
    {
      if(tid == 0){ cout << y_min << " <= j <= " << (y_min + r_y) << " & j < " << NOy << endl; }
      }*/

    for(j = y_min; (j <= (y_min + r_y)) && (j < NOy); ++j){
      
      // debugging
      /*#pragma omp barrier
      #pragma omp critical
      cout << "tid = " << tid << ", j = " << j << endl;
      #pragma omp barrier*/

      // I) Re-initilise the test_region array for the start of a new line and update buffer with the new value
      NOt = 0;
      j_min = j - r_y;
      j_min = (j_min >= 0) ? j_min : 0;

      // debugging
      /*#pragma omp barrier
      #pragma omp critical
      cout << "tid = " << tid << ", j_min = " << j_min << " for j = " << j << " & r_y = " << r_y << endl;
      #pragma omp barrier*/

      for(j2 = j_min; (j2 <= j + r_y) && (j2 < NOy); ++j2){
	
	for(i = 0; i <= r_x; ++i){
	  
	  if(((static_cast<double>(i)*static_cast<double>(i)/(static_cast<double>(r_x) * static_cast<double>(r_x))) + 
	      (static_cast<double>((j2-j))*static_cast<double>((j2-j))/(static_cast<double>(r_y) * static_cast<double>(r_y)))) <= 1.0){
	    
	    test_region[curr_test_set][NOt] = (j2 * NOx) + i;
	    ++NOt;
	    
	  } else { break; }
	  
	  // for(i = 0; i <= r_x; ++i)
	}
	
	// for(j2 = j_min; (j2 <= j + r_y) && (j2 < NOy); ++j2)
      }

      // debugging
      /*#pragma omp barrier
      #pragma omp critical
      { 
        cout << "tid = " << tid << ", RAW NOt = " << NOt << ": " << flush;
	for(t = 0; t < NOt; ++t){ cout << test_region[curr_test_set][t] << " "; }
	cout << endl;
      }
      #pragma omp barrier*/

      HeapSort_index(NOt,image_orig,test_region[curr_test_set] - 1);

      // debugging
      /*#pragma omp barrier
      #pragma omp critical
      { 
        cout << "tid = " << tid << ", SORTED NOt = " << NOt << ": " << flush;
	for(t = 0; t < NOt; ++t){ cout << test_region[curr_test_set][t] << " "; }
	cout << endl;
      }
      #pragma omp barrier*/

      // advanced version
      buffer[(j - y_min)*NOx] = filter(image_orig,test_region[curr_test_set],NOt);
      // crude version
      //buffer[(j - y_min)*NOx] = filter(image_orig,test_region[curr_test_set],NOt);
     
      // debugging
      /*#pragma omp barrier
      #pragma omp master 
      cout << "d1" << endl;
      #pragma omp barrier*/
      
      // II) Re-initialise the array mapping relative position to test_region index
      for(t = 0; t < NOt; ++t){
    
        j2 = test_region[curr_test_set][t]/NOx;
	i2 = test_region[curr_test_set][t] - (j2 * NOx);
	map_p_to_t[(j2 - j + r_y) * (1 + 2*r_x) + i2 + r_x] = t;
	
      }
            
      // debugging
      /*#pragma omp barrier
      #pragma omp master 
      cout << "d2" << endl;
      #pragma omp barrier*/

      // III) Update the test_region array, then update the buffer with the new value
      for(i = 1; i < NOx; ++i){

        // debugging
        /*#pragma omp barrier
        #pragma omp critical
        cout << "tid = " << tid << ", i = " << i << endl;
	#pragma omp barrier*/

	// debugging
	/*#pragma omp barrier
	#pragma omp critical
	{
          cout << "tid = " << tid << ", offsets (" << NOos << "): " << flush;
	  for(os = 0; os < NOos; ++os){ cout << offsets[2*os] << ", " << offsets[(2*os)+1] << " " << flush; }
	  cout << endl;
	}*/

	// a) Replace each pixel in the previous trailing edge in the test_region array with the mirrored leading edge value if possible.
	//    Replace the previous trailing edge with the max_value_index if the leading edge is out of bounds.
	//    Extend the number of values in the test_region if the previous trailing edge is out of bounds.
	NOt3 = 0;
	for(os = 0; os < NOos; ++os){
	  
	  trailing_pos = ((j + offsets[2*os]) * NOx) + i - offsets[(2*os)+1] - 1;
	  leading_pos = ((j + offsets[2*os]) * NOx) + i + offsets[(2*os)+1];
	  if((j + offsets[2*os]) >= 0 && (j + offsets[2*os]) < NOy && (i - offsets[(2*os)+1] - 1) >= 0 && (i - offsets[(2*os)+1] - 1) < NOx){ 
	    
	    if((i + offsets[(2*os)+1]) >= 0 && (i + offsets[(2*os)+1]) < NOx){ 
	      
	      test_region[curr_test_set][map_p_to_t[(offsets[2*os] + r_y)*(1 + (2*r_x)) - offsets[(2*os)+1] + r_x]] = -1;
	      test_region[2][NOt3++] = leading_pos;
	      
	      // if(leading_pos >= 0 && leading_pos < NOd){
	    } else {
	      
	      test_region[curr_test_set][map_p_to_t[(offsets[2*os] + r_y)*(1 + (2*r_x)) - offsets[(2*os)+1] + r_x]] = -1;
	      
	      // else --- if(leading_pos >= 0 && leading_pos < NOd){
	    }
	    
	    //} else if(leading_pos >= 0 && leading_pos < NOd){
	  } else if((j + offsets[2*os]) >= 0 && (j + offsets[2*os]) < NOy && (i + offsets[(2*os)+1]) >= 0 && (i + offsets[(2*os)+1]) < NOx){
	    
	    test_region[2][NOt3++] = leading_pos;
	    
	    // else --- if(trailing_pos >= 0 && trailing_pos < NOd){ 
	  }
	  
	  // for(os = 0; os < NOos; ++os){
	}

	// debugging
	/*#pragma omp barrier
	#pragma omp critical
	{
           cout << "tid = " << tid << ", RAW NOt3 = " << NOt3 << ": " << flush;
	   for(t3 = 0; t3 < NOt3; ++t3){ cout << test_region[2][t3] << " "; }
	   cout << endl;
        }
	#pragma omp barrier*/
	
	// b. Sort the next_test_set test_region values, merge them with the current set and update the buffer
	HeapSort_index(NOt3,image_orig,test_region[2] - 1);

	// debugging
	/*#pragma omp barrier
	#pragma omp critical
	{
           cout << "tid = " << tid << ", SORTED NOt3 = " << NOt3 << ": " << flush;
           for(t3 = 0; t3 < NOt3; ++t3){ cout << test_region[2][t3] << " "; }
	   cout << endl;
        }
	#pragma omp barrier*/

	for(t = 0, t2 = 0, t3 = 0; t < NOt; ++t){
	  
	  if(test_region[curr_test_set][t] >= NOd){ continue; }
	  
	  if(t3 < NOt3){
	    
	    if(image_orig[test_region[curr_test_set][t]] <= image_orig[test_region[2][t3]]){
	      
	      test_region[next_test_set][t2++] = test_region[curr_test_set][t];
	      
	    } else {
	      
	      test_region[next_test_set][t2++] = test_region[2][t3++];
	      --t;
	      
	    }
	    
	  } else {
	    
	    test_region[next_test_set][t2++] = test_region[curr_test_set][t];
	    
	  }
	  
        }

	// debugging
	/*#pragma omp barrier
	#pragma omp critical
	cout << "tid = " << tid << ", exhausted inital test_region." << endl;
	#pragma omp barrier*/

	while(t3 < NOt3){ test_region[next_test_set][t2++] = test_region[2][t3++]; }

	// debugging
	/*#pragma omp barrier
	#pragma omp critical
	cout << "tid = " << tid << ", exhausted new test_region." << endl;
	#pragma omp barrier*/

	NOt = t2;
	curr_test_set = next_test_set;
	++next_test_set;
	if(next_test_set > 1){ next_test_set = 0; }	

	// advanced version
	buffer[((j - y_min)*NOx) + i] = filter(image_orig,test_region[curr_test_set],NOt);
	// crude version
	//buffer[((j - y_min)*NOx) + i] = filter(image_orig,test_region[curr_test_set],NOt);

	// debugging
	/*#pragma omp barrier
	#pragma omp critical
	{
	  cout << "tid = " << tid << ", NOt = " << NOt << ", curr_test_set = " << curr_test_set << ", next_test_set = " << next_test_set << ": " << flush;
	  for(t = 0; t < NOt; ++t){ cout << test_region[curr_test_set][t] << " "; }
	  cout << endl;
	}
	#pragma omp barrier*/

	// c. Re-create the mapping from relative position to test region index
	//#pragma omp critical
	//{
	  for(t = 0; t < NOt; ++t){

	    // debugging
	    //cout << "tid = " << tid << ", t = " << t << " " << flush;
	    
	    j2 = test_region[curr_test_set][t]/NOx;
	    i2 = test_region[curr_test_set][t] - (j2 * NOx);
	    
	    // debugging
	    //cout << "i2 = " << i2 << ", j2 = " << j2 << " --> " << ((j2 - j + r_y) * (1 + 2*r_x) + i2 - i + r_x) << " " << flush;

	    map_p_to_t[(j2 - j + r_y) * (1 + 2*r_x) + i2 - i + r_x] = t;

	    // debugging
	    //cout << "map_p_to_t = " << map_p_to_t[(j2 - j + r_y) * (1 + 2*r_x) + i2 - i + r_x] << endl;
	    
	  }
	  //}
	
	// debugging
	/*#pragma omp barrier
	#pragma omp critical
	cout << "tid = " << tid << ", finished remapping." << endl;
	#pragma omp barrier*/

	// update progress
	#pragma omp master
	{
	  if(vb_flag > 0){
	    
	    while((((static_cast<double>(j - y_min) * static_cast<double>(NOx)) + static_cast<double>(i) + 1.0) / static_cast<double>(NOx * (1 + r_y))) >= prog){ 
	      
	      prog+=0.001;
	      printf("\r%06.2f%%",(100.0*prog));
	      cout << flush;
	      
	    }
	    
	  }
	}
	
	// for(i = 1; i < NOx; ++i){
      }

      // debugging
      /*#pragma omp barrier
      #pragma omp master 
      cout << "d3" << endl;
      #pragma omp barrier*/
      
      // for(j = 0; j <= r_y; ++j){
    }
    
    // debugging
    /*#pragma omp barrier
    #pragma omp master 
    cout << "e" << endl;
    #pragma omp barrier*/

    #pragma omp master
    {
      if(vb_flag > 0){ 
      
	cout << " Done.\nFinished initialising buffer. \nProcessing input image using buffer . . . " << endl; 
	cout << "000.00%" << flush;
	
      }
    }

    // debugging
    /*#pragma omp critical
    {
      if(tid == 0){ cout << (y_min + 1 + r_y) << " <= j < " << y_max << " & " << NOy << endl; }
      }*/

    // 16. Process the input image using the buffer.
    //buffer_read = 0;
    buffer_write = r_y;
    buffer_read = -1;
    curr_test_set = 0;
    next_test_set = 1;
    prog = 0.0;
    for(j = y_min + 1 + r_y; (j < y_max) && (j < NOy); ++j){
      
      // Update the buffer_read and buffer_write lines.
      ++buffer_read;
      if(buffer_read > r_y){ buffer_read = 0; }

      // Update the input image with the buffer_read line of the buffer.
      for(i = 0; i < NOx; ++i){ image[((j - 1 - r_y) * NOx) + i] = buffer[(buffer_read * NOx) + i]; }
      // crude version
      //for(i = 0; i < NOx; ++i){ image[((j - 1 - r_y) * NOx) + i] = buffer[i]; }

      // Shuffle the buffer contents. ---- ONLY required by crude version
      /*for(j2 = 1; j2 <= r_y; ++j2){
	for(i = 0; i < NOx; ++i){ buffer[((j2 - 1) * NOx) + i] = buffer[(j2 * NOx) + i]; }
      }*/
      
      // debugging
      /*#pragma omp critical
      {
        if(tid == 0){ cout << "j = " << j << " of " << (y_min + 1 + r_y) << " <= j < " << y_max << " & j < " << NOy << endl; }
	}*/

      // Update the buffer_read and buffer_write lines.
      /*buffer_write = buffer_read;
      ++buffer_read;
      if(buffer_read > r_y){ buffer_read = 0; }*/
      
      // Process the current line and write the results to the buffer. 
      
      // I) Re-initilise the test_region array for the start of a new line and update buffer with the new value
      NOt = 0;
      j_min = j - r_y;
      j_min = (j_min >= 0) ? j_min : 0;
      for(j2 = j_min; (j2 <= j + r_y) && (j2 < NOy); ++j2){
	
	for(i = 0; i <= r_x; ++i){
	  
	  if(((static_cast<double>(i)*static_cast<double>(i)/(static_cast<double>(r_x) * static_cast<double>(r_x))) + 
	      (static_cast<double>((j2-j))*static_cast<double>((j2-j))/(static_cast<double>(r_y) * static_cast<double>(r_y)))) <= 1.0){
	    
	    test_region[curr_test_set][NOt] = (j2 * NOx) + i;
	    ++NOt;
	    
	  } else { break; }
	  
	  // for(i = 0; i <= r_x; ++i)
	}
	
	// for(j2 = j_min; (j2 <= j + r_y) && (j2 < NOy); ++j2)
      }
      HeapSort_index(NOt,image_orig,test_region[curr_test_set] - 1);
      // advanced version
      //buffer[buffer_write*NOx] = filter(image_orig,test_region[curr_test_set],NOt);
      buffer[buffer_read*NOx] = filter(image_orig,test_region[curr_test_set],NOt);
      // crude version
      //buffer[r_y*NOx] = filter(image_orig,test_region[curr_test_set],NOt);
      
      // II) Re-initialise the array mapping relative position to test_region index
      for(t = 0; t < NOt; ++t){
	
	j2 = test_region[curr_test_set][t]/NOx;
	i2 = test_region[curr_test_set][t] - (j2 * NOx);
	map_p_to_t[(j2 - j + r_y) * (1 + 2*r_x) + i2 + r_x] = t;
	
      }
      
      // III) Update the test_region array, then update the buffer with the new value
      for(i = 1; i < NOx; ++i){
	
	// a) Replace each pixel in the previous trailing edge in the test_region array with the mirrored leading edge value if possible.
	//    Replace the previous trailing edge with the max_value_index if the leading edge is out of bounds.
	//    Extend the number of values in the test_region if the previous trailing edge is out of bounds.
	NOt3 = 0;
	for(os = 0; os < NOos; ++os){
	  
	  trailing_pos = ((j + offsets[2*os]) * NOx) + i - offsets[(2*os)+1] - 1;
	  leading_pos = ((j + offsets[2*os]) * NOx) + i + offsets[(2*os)+1];
	  if((j + offsets[2*os]) >= 0 && (j + offsets[2*os]) < NOy && (i - offsets[(2*os)+1] - 1) >= 0 && (i - offsets[(2*os)+1] - 1) < NOx){ 
	    
	    if((i + offsets[(2*os)+1]) >= 0 && (i + offsets[(2*os)+1]) < NOx){ 
	      
	      test_region[curr_test_set][map_p_to_t[(offsets[2*os] + r_y)*(1 + (2*r_x)) - offsets[(2*os)+1] + r_x]] = -1;
	      test_region[2][NOt3++] = leading_pos;
	      
	      // if(leading_pos >= 0 && leading_pos < NOd){
	    } else {
	      
	      test_region[curr_test_set][map_p_to_t[(offsets[2*os] + r_y)*(1 + (2*r_x)) - offsets[(2*os)+1] + r_x]] = -1;
	      
	      // else --- if(leading_pos >= 0 && leading_pos < NOd){
	    }
	    
	    //} else if(leading_pos >= 0 && leading_pos < NOd){
	  } else if((j + offsets[2*os]) >= 0 && (j + offsets[2*os]) < NOy && (i + offsets[(2*os)+1]) >= 0 && (i + offsets[(2*os)+1]) < NOx){
	    
	    test_region[2][NOt3++] = leading_pos;
	    
	    // else --- if(trailing_pos >= 0 && trailing_pos < NOd){ 
	  } 
	  
	  // for(os = 0; os < NOos; ++os){
	}
	
	// b. Sort the next_test_set test_region values, merge them with the current set and update the buffer
	HeapSort_index(NOt3,image_orig,test_region[2] - 1);
	for(t = 0, t2 = 0, t3 = 0; t < NOt; ++t){
	  
	  if(test_region[curr_test_set][t] >= NOd){ continue; }
	  
	  if(t3 < NOt3){
	    
	    if(image_orig[test_region[curr_test_set][t]] <= image_orig[test_region[2][t3]]){
	      
	      test_region[next_test_set][t2++] = test_region[curr_test_set][t];
	      
	    } else {
	      
	      test_region[next_test_set][t2++] = test_region[2][t3++];
	      --t;
	      
	    }
	    
	  } else {
	    
	    test_region[next_test_set][t2++] = test_region[curr_test_set][t];
	    
	  }
	  
	}
	while(t3 < NOt3){ test_region[next_test_set][t2++] = test_region[2][t3++]; }
	NOt = t2;
	curr_test_set = next_test_set;
	++next_test_set;
	if(next_test_set > 1){ next_test_set = 0; }
	// advanced version
	//buffer[(buffer_write*NOx) + i] = filter(image_orig,test_region[curr_test_set],NOt);
	buffer[(buffer_read*NOx) + i] = filter(image_orig,test_region[curr_test_set],NOt);
	// crude version
	//buffer[(r_y*NOx) + i] = filter(image_orig,test_region[curr_test_set],NOt);
	
	// c. Re-create the mapping from relative position to test region index
	for(t = 0; t < NOt; ++t){
	  
	  j2 = test_region[curr_test_set][t]/NOx;
	  i2 = test_region[curr_test_set][t] - (j2 * NOx);
	  map_p_to_t[(j2 - j + r_y) * (1 + 2*r_x) + i2 - i + r_x] = t;
	  
	}
		
	// update progress
	#pragma omp master
	{
	  if(vb_flag > 0){
	  
  while((((static_cast<double>(j - y_min) * static_cast<double>(NOx)) + static_cast<double>(i) + 1.0) / static_cast<double>(NOx * (y_max - y_min))) >= prog){ 
	    
	      prog+=0.001;
	      printf("\r%06.2f%%",(100.0*prog));
	      cout << flush;
	      
	    }
	  
	  }
	}
	
	// for(i = 1; i < NOx; ++i){
      }
      
      // for(j = 1 + r_y; j < NOy; ++j){
    }
    #pragma omp master
    {
      if(vb_flag > 0){ cout << " Done." << endl; }
      if(vb_flag > 1){ cout << "Updating image with the final buffer values . . . " << endl; }
    }
    
    // 17. Update the image with the final buffer values
    // advanced version
    if(buffer_read < r_y){
      for(j = 0; j < (r_y - buffer_read); ++j){
	  
	//cout << "j = " << j << endl;
 
	for(i = 0; i < NOx; ++i){
	
	  //cout << "i = " << i << ": updating image[(" << y_max << " - 1 - " << r_y << " + " << j << ") * " << NOx << " + " << i  << " = "  << (((y_max - 1 - r_y + j) * NOx) + i) << "] to buffer[(" << buffer_read << " + 1 + " << j << ") * " << NOx << " + " << i << " = " << (((buffer_read + 1 + j) * NOx) + i) << "] . . . " << flush;
	    
	  image[((y_max - 1 - r_y + j) * NOx) + i] = buffer[((buffer_read + 1 + j) * NOx) + i];
	    
	  //cout << "Done." << endl;
	    
	}
	
      }
    }
    for(j = 0; j <= buffer_read; ++j){

      for(i = 0; i < NOx; ++i){
      
        image[((y_max - 1 - buffer_read + j) * NOx) + i] = buffer[(j * NOx) + i];
      
      }

    }
    // crude version
    /*for(j = 0; j <= r_y; ++j){
      
      for(i = 0; i < NOx; ++i){
	
	image[((y_min + NOy - 1 - r_y + j) * NOx) + i] = buffer[(j * NOx) + i];
	
	// for(i = 0; i < NOx; ++i){
      }
      
      // for(j = 0; j < (1 + r_y); ++j){
      }*/
    
    #pragma omp master
    {
      if(vb_flag > 1){ cout << "Done." << endl; }
      if(vb_flag > 1){ cout << "Freeing up memory . . . " << endl; }
    }

    // 18. free up memory
    delete [] test_region[0];
    delete [] test_region[1];
    delete [] test_region[2];
    delete [] test_region;
    delete [] map_p_to_t;
    delete [] buffer;
    
    // end of parallel section
  }
  
  // 19. free up more memory
  delete [] image_orig;
  delete [] offsets;

}

// c. 3D: datacube


#endif
