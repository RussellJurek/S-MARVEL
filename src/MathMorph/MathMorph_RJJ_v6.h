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
  unsigned int t,NOt,t2,t3,NOt3,os,NOos_x,NOos_y,NOb,NOp;
  size_t ** test_region, * map_p_to_t, * prev_edge_region, d, NOd, buffer_read;
  ssize_t i,j,i2,j2,j_min, ** offsets, trailing_pos, leading_pos;
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

  // 4. Count the number of pixels in the leading/trailing edge of the mask --> NOos_x.
  NOos_x = 1 + (2*r_y);
  NOos_y = 1 + (2*r_x);
  if(vb_flag > 1){ cout << "Calculated that mask leading/trailing edge contains " << NOos_x << " (x direction) & " << NOos_y << " (y direction) pixels." << endl; }

  // 5. Calculate the buffer size.
  NOb = NOx * (1 + r_y);
  if(vb_flag > 1){ cout << "Calculated buffer size of " << NOb << "." << endl; }

  // 6. Create an array to hold NOos_x & NOos_y offsets corresponding to the leading/trailing mask edge.
  offsets = new ssize_t * [2];
  offsets[0] = new ssize_t[2*NOos_x];
  offsets[1] = new ssize_t[2*NOos_y];
  if(vb_flag > 1){ cout << "Allocated offsets[0] array containing " << (2*NOos_x) << " elements for " << NOos_x << " pixels." << endl << "Allocated offsets[1] array containing " << (2*NOos_y) << " elements for " << NOos_y << " pixels." << endl; }
      
  // 7. Populate the offsets array.
  if(vb_flag > 1){ cout << "Populating the offsets array . . . " << endl; }
  NOos_y = NOos_x = 0;
  if(elliptical){    
    
    if(vb_flag > 1){ cout << "Finding edge of elliptical mask . . . " << endl; }

    // for each j line, find the last pixel that belongs to the mask and store the correspondng offset
    for(j = -1 * r_y; j <= r_y; ++j){
      
      for(i = 0; i <= r_x; ++i){
	
	if(((static_cast<double>(i)*static_cast<double>(i)/(static_cast<double>(r_x) * static_cast<double>(r_x))) + 
	    (static_cast<double>(j)*static_cast<double>(j)/(static_cast<double>(r_y) * static_cast<double>(r_y)))) <= 1.0){
	  
	  offsets[0][2*NOos_x] = j;
	  offsets[0][(2*NOos_x)+1] = i;
	  
	} else {
	  
	  ++NOos_x;
	  break;
	  
	}
	
      }
      
      if(i > r_x){ ++NOos_x; }
      
    }
    
    // for each i line, find the last pixel that belongs to the mask and store the correspondng offset
    for(i = -1 * r_x; i <= r_x; ++i){
      
      for(j = 0; j <= r_y; ++j){
	
	if(((static_cast<double>(i)*static_cast<double>(i)/(static_cast<double>(r_x) * static_cast<double>(r_x))) + 
	    (static_cast<double>(j)*static_cast<double>(j)/(static_cast<double>(r_y) * static_cast<double>(r_y)))) <= 1.0){
	  
	  offsets[1][2*NOos_y] = j;
	  offsets[1][(2*NOos_y)+1] = i;
	    
	} else {
	  
	  ++NOos_y;
	  break;
	  
	}
	
      }
      
      if(j > r_y){ ++NOos_y; }
      
    }

  } else {
    
    if(vb_flag > 1){ cout << "Finding edge of rectangular mask . . . " << endl; }

    // for each j line, store the offset of the last pixel in the mask
    for(j = -r_y; j <= r_y; ++j){
      
      offsets[0][2*NOos_x] = j;
      offsets[0][(2*NOos_x)+1] = r_x;
      ++NOos_x;
      
    }
    
    // for each i line, store the offset of the last pixel in the mask
    for(i = -r_x; i <= r_x; ++i){
      
      offsets[1][2*NOos_y] = r_y;
      offsets[1][(2*NOos_y)+1] = i;
      ++NOos_y;
      
    }

  }
  if(vb_flag > 1){ cout << "Done. Stored " << NOos_x << " offsets." << endl; }
  
  // 8. Retrieve the number of threads available and specify this many to be used
  NOtid = omp_get_max_threads();
  if(vb_flag > 1){ cout << "Using " << NOtid << " threads." << endl; }
  omp_set_num_threads(NOtid);
  
  // 9. Start of parallel region
#pragma omp parallel default(shared) firstprivate(NOt) private(test_region,map_p_to_t,buffer,prev_edge_region,leading_pos,trailing_pos,i,i2,j,j_min,j2,os,t,t2,t3,NOt3,NOp,y_min,y_max,tid,curr_test_set,next_test_set,buffer_read,prog)
  {
    
    // 10. Retrive the thread id
    tid = omp_get_thread_num();
    
    // 11. Generate the image boundaries for this thread
    y_min = tid * static_cast<int>(ceil(static_cast<double>(NOy) / static_cast<double>(NOtid)));
    y_max = (tid + 1) * static_cast<int>(ceil(static_cast<double>(NOy) / static_cast<double>(NOtid)));
    if(y_max > NOy){ y_max = NOy; }
    
    // 12. Create an array that holds the mapping of relative positions to test region locations
    map_p_to_t = new size_t[(1 + 2*r_x)*(1 + 2*r_y)];
    for(i = 0; i < ((1 + 2*r_x)*(1 + 2*r_y)); ++i){ map_p_to_t[i] = 0; }
#pragma omp master
    {
      if(vb_flag > 1){ cout << "Allocated map_p_to_t array containing " << ((1 + 2*r_x)*(1 + 2*r_y)) << " pixels." << endl; }
    }
    
    // 13. Create an array that can hold NOt index values for the current and next test regions being filtered.
    test_region = new size_t * [3];
    test_region[0] = new size_t [NOt];
    test_region[1] = new size_t [NOt];
    test_region[2] = new size_t [1 + (2*r_y)];
#pragma omp master
    {
      if(vb_flag > 1){ cout << "Allocated test_region array containing " << NOt << " pixels." << endl; }
    }
    
    // 14. Create an array that can hold NOt index values for the previous edge region at j' = j - 1, i = 0.
    prev_edge_region = new size_t [NOt];
    
    // 14. Create the buffer array.
    buffer = new data_type[NOb];
#pragma omp master
    {
      if(vb_flag > 1){ cout << "Allocated buffer array." << endl; }
    }
    
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
    j = y_min;

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
    buffer[(j - y_min)*NOx] = filter(image_orig,test_region[curr_test_set],NOt);
    NOp = NOt;    
    for(t = 0; t < NOp; ++t){ prev_edge_region[t] = test_region[curr_test_set][t]; }
    
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
      for(os = 0; os < NOos_x; ++os){
	
        trailing_pos = ((j + offsets[0][2*os]) * NOx) + i - offsets[0][(2*os)+1] - 1;
	leading_pos = ((j + offsets[0][2*os]) * NOx) + i + offsets[0][(2*os)+1];
	if((j + offsets[0][2*os]) >= 0 && (j + offsets[0][2*os]) < NOy && (i - offsets[0][(2*os)+1] - 1) >= 0 && (i - offsets[0][(2*os)+1] - 1) < NOx){ 
	  
          if((i + offsets[0][(2*os)+1]) >= 0 && (i + offsets[0][(2*os)+1]) < NOx){ 
	    
            test_region[curr_test_set][map_p_to_t[(offsets[0][2*os] + r_y)*(1 + (2*r_x)) - offsets[0][(2*os)+1] + r_x]] = -1;
	    test_region[2][NOt3++] = leading_pos;
	    
	    // if(leading_pos >= 0 && leading_pos < NOd){
          } else {
	    
            test_region[curr_test_set][map_p_to_t[(offsets[0][2*os] + r_y)*(1 + (2*r_x)) - offsets[0][(2*os)+1] + r_x]] = -1;
	    
	    // else --- if(leading_pos >= 0 && leading_pos < NOd){
          }
	  
	  //} else if(leading_pos >= 0 && leading_pos < NOd){
        } else if((j + offsets[0][2*os]) >= 0 && (j + offsets[0][2*os]) < NOy && (i + offsets[0][(2*os)+1]) >= 0 && (i + offsets[0][(2*os)+1]) < NOx){
	  
          test_region[2][NOt3++] = leading_pos;
	  
	  // else --- if(trailing_pos >= 0 && trailing_pos < NOd){ 
        }
	
	// for(os = 0; os < NOos_x; ++os){
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
      buffer[((j - y_min)*NOx) + i] = filter(image_orig,test_region[curr_test_set],NOt);
      
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
	  
          while((((static_cast<double>(j - y_min) * static_cast<double>(NOx)) + static_cast<double>(i) + 1.0) / static_cast<double>(NOx * (1 + r_y))) >= prog){ 
	    
            prog+=0.001;
	    printf("\r%06.2f%%",(100.0*prog));
	    cout << flush;
	    
          }
	  
        }
      }
      
      // for(i = 1; i < NOx; ++i){
    }

    for(j = y_min + 1; (j <= (y_min + r_y)) && (j < NOy); ++j){

      // I) Replace each pixel in the previous trailing edge in the test_region array with the mirrored leading edge value if possible.
      //    Replace the previous trailing edge with the max_value_index if the leading edge is out of bounds.
      //    Extend the number of values in the test_region if the previous trailing edge is out of bounds.
      NOt = NOp;
      for(t = 0; t < NOt; ++t){
      
        j2 = prev_edge_region[t]/NOx;
        i2 = prev_edge_region[t] - (j2 * NOx);
        map_p_to_t[(j2 - j + 1 + r_y) * (1 + 2*r_x) + i2 + r_x] = t;
        test_region[curr_test_set][t] = prev_edge_region[t];
      
      }
      NOt3 = 0;
      for(os = r_x; os < NOos_y; ++os){
            
        trailing_pos = ((j - offsets[1][2*os] - 1) * NOx) + offsets[1][(2*os)+1];
        leading_pos = ((j + offsets[1][2*os]) * NOx) + offsets[1][(2*os)+1];
        if((j - offsets[1][2*os] - 1) >= 0 && (j - offsets[1][2*os] - 1) < NOy && offsets[1][(2*os)+1] >= 0 && offsets[1][(2*os)+1] < NOx){ 
	  
          if((j + offsets[1][2*os]) >= 0 && (j + offsets[1][2*os]) < NOy){ 
	    
	    test_region[curr_test_set][map_p_to_t[(r_y - offsets[1][2*os])*(1 + (2*r_x)) + offsets[1][(2*os)+1] + r_x]] = -1;
	    test_region[2][NOt3++] = leading_pos;

	    // if(leading_pos >= 0 && leading_pos < NOd){
	  } else {
	    
            test_region[curr_test_set][map_p_to_t[(r_y - offsets[1][2*os])*(1 + (2*r_x)) + offsets[1][(2*os)+1] + r_x]] = -1;

	    // else --- if(leading_pos >= 0 && leading_pos < NOd){
	  }
	  
	  //} else if(leading_pos >= 0 && leading_pos < NOd){
	} else if((j + offsets[1][2*os]) >= 0 && (j + offsets[1][2*os]) < NOy && offsets[1][(2*os)+1] >= 0 && offsets[1][(2*os)+1] < NOx){
	  
	    test_region[2][NOt3++] = leading_pos;

	  // else --- if(trailing_pos >= 0 && trailing_pos < NOd){ 
	}

	// for(os = 0; os < NOos_x; ++os){
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
      buffer[((j - y_min)*NOx)] = filter(image_orig,test_region[curr_test_set],NOt);
      for(t = 0; t < NOt; ++t){ prev_edge_region[t] = test_region[curr_test_set][t]; }
      NOp = NOt;

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
	for(os = 0; os < NOos_x; ++os){
	  
	  trailing_pos = ((j + offsets[0][2*os]) * NOx) + i - offsets[0][(2*os)+1] - 1;
	  leading_pos = ((j + offsets[0][2*os]) * NOx) + i + offsets[0][(2*os)+1];
	  if((j + offsets[0][2*os]) >= 0 && (j + offsets[0][2*os]) < NOy && (i - offsets[0][(2*os)+1] - 1) >= 0 && (i - offsets[0][(2*os)+1] - 1) < NOx){ 
	    
	    if((i + offsets[0][(2*os)+1]) >= 0 && (i + offsets[0][(2*os)+1]) < NOx){ 
	      
	      test_region[curr_test_set][map_p_to_t[(offsets[0][2*os] + r_y)*(1 + (2*r_x)) - offsets[0][(2*os)+1] + r_x]] = -1;
	      test_region[2][NOt3++] = leading_pos;

	      // if(leading_pos >= 0 && leading_pos < NOd){
	    } else {
	      
	      test_region[curr_test_set][map_p_to_t[(offsets[0][2*os] + r_y)*(1 + (2*r_x)) - offsets[0][(2*os)+1] + r_x]] = -1;

	      // else --- if(leading_pos >= 0 && leading_pos < NOd){
	    }
	    
	    //} else if(leading_pos >= 0 && leading_pos < NOd){
	  } else if((j + offsets[0][2*os]) >= 0 && (j + offsets[0][2*os]) < NOy && (i + offsets[0][(2*os)+1]) >= 0 && (i + offsets[0][(2*os)+1]) < NOx){

	    test_region[2][NOt3++] = leading_pos;
	    // else --- if(trailing_pos >= 0 && trailing_pos < NOd){ 
	  }
	  
	  // for(os = 0; os < NOos_x; ++os){
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
	buffer[((j - y_min)*NOx) + i] = filter(image_orig,test_region[curr_test_set],NOt);

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
	    
	    while((((static_cast<double>(j - y_min) * static_cast<double>(NOx)) + static_cast<double>(i) + 1.0) / static_cast<double>(NOx * (1 + r_y))) >= prog){ 
	      
	      prog+=0.001;
	      printf("\r%06.2f%%",(100.0*prog));
	      cout << flush;
	      
	    }
	    
	  }
	}
	
	// for(i = 1; i < NOx; ++i){
      }
      
      // for(j = 0; j <= r_y; ++j){
    }

#pragma omp master
    {
      if(vb_flag > 0){ 
	
	cout << " Done.\nFinished initialising buffer. \nProcessing input image using buffer . . . " << endl; 
	cout << "000.00%" << flush;
	
      }
    }
    
    // 16. Process the input image using the buffer.
    buffer_read = -1;
    curr_test_set = 0;
    next_test_set = 1;
    prog = 0.0;
    for(j = y_min + 1 + r_y; (j < y_max) && (j < NOy); ++j){
      
      // Update the buffer_read line.
      ++buffer_read;
      if(buffer_read > r_y){ buffer_read = 0; }
      
      // Update the input image with the buffer_read line of the buffer.
      for(i = 0; i < NOx; ++i){ image[((j - 1 - r_y) * NOx) + i] = buffer[(buffer_read * NOx) + i]; }      
      
      // Process the current line and write the results to the buffer. 
      
      // I) Replace each pixel in the previous trailing edge in the test_region array with the mirrored leading edge value if possible.
      //    Replace the previous trailing edge with the max_value_index if the leading edge is out of bounds.
      //    Extend the number of values in the test_region if the previous trailing edge is out of bounds.
      NOt = NOp;
      for(t = 0; t < NOt; ++t){
      
        j2 = prev_edge_region[t]/NOx;
        i2 = prev_edge_region[t] - (j2 * NOx);
        map_p_to_t[(j2 - j + 1 + r_y) * (1 + 2*r_x) + i2 + r_x] = t;
        test_region[curr_test_set][t] = prev_edge_region[t];
      
      }
      NOt3 = 0;
      for(os = r_x; os < NOos_y; ++os){
            
        trailing_pos = ((j - offsets[1][2*os] - 1) * NOx) + offsets[1][(2*os)+1];
        leading_pos = ((j + offsets[1][2*os]) * NOx) + offsets[1][(2*os)+1];
        if((j - offsets[1][2*os] - 1) >= 0 && (j - offsets[1][2*os] - 1) < NOy && offsets[1][(2*os)+1] >= 0 && offsets[1][(2*os)+1] < NOx){ 
	  
          if((j + offsets[1][2*os]) >= 0 && (j + offsets[1][2*os]) < NOy){ 
	    
	    test_region[curr_test_set][map_p_to_t[(r_y - offsets[1][2*os])*(1 + (2*r_x)) + offsets[1][(2*os)+1] + r_x]] = -1;
	    test_region[2][NOt3++] = leading_pos;

	    // if(leading_pos >= 0 && leading_pos < NOd){
	  } else {
	    
            test_region[curr_test_set][map_p_to_t[(r_y - offsets[1][2*os])*(1 + (2*r_x)) + offsets[1][(2*os)+1] + r_x]] = -1;

	    // else --- if(leading_pos >= 0 && leading_pos < NOd){
	  }
	  
	  //} else if(leading_pos >= 0 && leading_pos < NOd){
	} else if((j + offsets[1][2*os]) >= 0 && (j + offsets[1][2*os]) < NOy && offsets[1][(2*os)+1] >= 0 && offsets[1][(2*os)+1] < NOx){
	  
	    test_region[2][NOt3++] = leading_pos;

	  // else --- if(trailing_pos >= 0 && trailing_pos < NOd){ 
	}

	// for(os = 0; os < NOos_x; ++os){
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
      buffer[buffer_read*NOx] = filter(image_orig,test_region[curr_test_set],NOt);
      for(t = 0; t < NOt; ++t){ prev_edge_region[t] = test_region[curr_test_set][t]; }
      NOp = NOt;
          
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
	for(os = 0; os < NOos_x; ++os){
	  
	  trailing_pos = ((j + offsets[0][2*os]) * NOx) + i - offsets[0][(2*os)+1] - 1;
	  leading_pos = ((j + offsets[0][2*os]) * NOx) + i + offsets[0][(2*os)+1];
	  if((j + offsets[0][2*os]) >= 0 && (j + offsets[0][2*os]) < NOy && (i - offsets[0][(2*os)+1] - 1) >= 0 && (i - offsets[0][(2*os)+1] - 1) < NOx){ 
	    
	    if((i + offsets[0][(2*os)+1]) >= 0 && (i + offsets[0][(2*os)+1]) < NOx){ 
	      
	      test_region[curr_test_set][map_p_to_t[(offsets[0][2*os] + r_y)*(1 + (2*r_x)) - offsets[0][(2*os)+1] + r_x]] = -1;
	      test_region[2][NOt3++] = leading_pos;
	      
	      // if(leading_pos >= 0 && leading_pos < NOd){
	    } else {
	      
	      test_region[curr_test_set][map_p_to_t[(offsets[0][2*os] + r_y)*(1 + (2*r_x)) - offsets[0][(2*os)+1] + r_x]] = -1;
	      
	      // else --- if(leading_pos >= 0 && leading_pos < NOd){
	    }
	    
	    //} else if(leading_pos >= 0 && leading_pos < NOd){
	  } else if((j + offsets[0][2*os]) >= 0 && (j + offsets[0][2*os]) < NOy && (i + offsets[0][(2*os)+1]) >= 0 && (i + offsets[0][(2*os)+1]) < NOx){
	    
	    test_region[2][NOt3++] = leading_pos;
	    
	    // else --- if(trailing_pos >= 0 && trailing_pos < NOd){ 
	  } 
	  
	  // for(os = 0; os < NOos_x; ++os){
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
	buffer[(buffer_read*NOx) + i] = filter(image_orig,test_region[curr_test_set],NOt);
	
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
    if(buffer_read < r_y){
      for(j = 0; j < (r_y - buffer_read); ++j){
	
	for(i = 0; i < NOx; ++i){
	  
	  image[((y_max - 1 - r_y + j) * NOx) + i] = buffer[((buffer_read + 1 + j) * NOx) + i];
	  
	}
	
      }
    }
    for(j = 0; j <= buffer_read; ++j){
      
      for(i = 0; i < NOx; ++i){
	
        image[((y_max - 1 - buffer_read + j) * NOx) + i] = buffer[(j * NOx) + i];
	
      }
      
    }
    
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
    delete [] prev_edge_region;
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
