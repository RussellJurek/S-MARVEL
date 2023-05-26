#include<iostream>
#include<string>
#include<sstream>
#include<math.h>
#include<vector>

#ifndef MathMorph_RJJ

#define MathMorph_RJJ

using namespace std;
 
void CreateMetric(int data_metric[3], int xyz_order[3], int NOx, int NOy, int NOz);

template <typename data_type>
void HeapSort_voxels(int n, data_type ra[]){

  int i,ir,j,l;
  data_type rra;

  if (n<2) return;
  l = (n >> 1) + 1;
  ir = n;
  for(;;){

    if (l > 1) {
      rra = ra[--l];
    } else {
      rra=ra[ir];
      ra[ir]=ra[1];
      if (--ir == 1) {
	ra[1] = rra;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
	ra[i] = ra[j];
	i = j;
	j <<= 1;
      } else break;
    }
    ra[i] = rra;
  }

}

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

template <typename data_type>
int Populate_temp_array(data_type * temp_vals, data_type * data_vals, int i, int j, int k, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag){

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
  if(round_flag < 0){

    // rectangular smoothing box

    for(z = z_start; ((z < size_z) && (z <= (k + fradius))); ++z){

      for(y = y_start; ((y < size_y) && (y <= (j + sradius))); ++y){

	for(x = x_start; ((x < size_x) && (x <= (i + sradius))); ++x){

	  temp_vals[used] = data_vals[((z * size_x * size_y) + (y * size_x) + x)];
	  used++;
	  
	}

      }

    }

  } else if(round_flag == 0){

    // cylindrical smoothing box

    for(z = z_start; ((z < size_z) && (z <= (k + fradius))); ++z){

      for(y = y_start; ((y < size_y) && (y <= (j + sradius))); ++y){

	for(x = x_start; ((x < size_x) && (x <= (i + sradius))); ++x){

	  sep = sqrtf(((float)((x - i) * (x - i)) + (float)((y - j) * (y - j))));
	  if(sep <= ((float) sradius)){

	    temp_vals[used] = data_vals[((z * size_x * size_y) + (y * size_x) + x)];
	    used++;

	  }
	  
	}

      }

    }


  } else {

    // rounded cylindrical smoothing box

    for(z = z_start; ((z < (z_start + sradius)) && (z<k) && (z < size_z)); ++z){

      // calculate projected radius for this channel
      sep_max = acos(((float)(z - z_start - sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;

      for(y = y_start; ((y <= (j + sradius)) && (y < size_y)); ++y){

	for(x = x_start; ((x <= (i + sradius)) && (x < size_x)); ++x){
    
	  sep = sqrtf(((float)((x - i) * (x - i)) + (float)((y - j) * (y - j))));
	  if(sep <= sep_max){ 

	    temp_vals[used] = data_vals[((z * size_x * size_y) + (y * size_x) + x)];
	    used++; 

	  }

	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	  
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }

      // for(k = (-1 * fradius); k <= fradius; ++k)
    }

    for(z = z_start + sradius; ((z < k) && (z < size_z)); ++z){

      // calculate projected radius for this channel
      sep_max = (float) sradius;

      for(y = y_start; ((y <= (j + sradius)) && (y < size_y)); ++y){

	for(x = x_start; ((x <= (i + sradius)) && (x < size_x)); ++x){
    
	  sep = sqrtf(((float)((x - i) * (x - i)) + (float)((y - j) * (y - j))));
	  if(sep <= sep_max){ 

	    temp_vals[used] = data_vals[((z * size_x * size_y) + (y * size_x) + x)];
	    used++; 
	    
	  }

	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	  
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }

      // for(k = (-1 * fradius); k <= fradius; ++k)
    }


    // calculate projected radius for this channel
    sep_max = (float) sradius;
    
    for(y = y_start; ((y <= (j + sradius)) && (y < size_y)); ++y){
      
      for(x = x_start; ((x <= (i + sradius)) && (x < size_x)); ++x){
	
	sep = sqrtf(((float)((x - i) * (x - i)) + (float)((y - j) * (y - j))));
	if(sep <= sep_max){ 
	  
	  temp_vals[used] = data_vals[((k * size_x * size_y) + (y * size_x) + x)];
	  used++; 
	  
	}
	
	// for(i = (-1 * sradius); i <= sradius; ++i)
      }
      
      // for(j = (-1 * sradius); j <= sradius; ++j)
    }

    for(z = k + 1; ((z <= (k + fradius - sradius)) && (z < size_z)); ++z){

      // calculate projected radius for this channel
      sep_max = (float) sradius;

      for(y = y_start; ((y <= (j + sradius)) && (y < size_y)); ++y){

	for(x = x_start; ((x <= (i + sradius)) && (x < size_x)); ++x){
    
	  sep = sqrtf(((float)((x - i) * (x - i)) + (float)((y - j) * (y - j))));
	  if(sep <= sep_max){ 

	    temp_vals[used] = data_vals[((z * size_x * size_y) + (y * size_x) + x)];
	    used++; 
	    
	  }

	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	  
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }

      // for(k = (-1 * fradius); k <= fradius; ++k)
    }

    for(z = k + fradius - sradius + 1; ((z <= (k + fradius)) && (z>k) && (z < size_z)); ++z){

      // calculate projected radius for this channel
      sep_max = acos(((float)(z - k - fradius + sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;

      for(y = y_start; ((y <= (j + sradius)) && (y < size_y)); ++y){

	for(x = x_start; ((x <= (i + sradius)) && (x < size_x)); ++x){
    
	  sep = sqrtf(((float)((x - i) * (x - i)) + (float)((y - j) * (y - j))));
	  if(sep <= sep_max){

	    temp_vals[used] = data_vals[((z * size_x * size_y) + (y * size_x) + x)];
	    used++; 
	    
	  }

	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	  
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }

      // for(k = (-1 * fradius); k <= fradius; ++k)
    }

  }

  return used;

}

template <typename data_type>
data_type Temp_array_min(data_type * temp_vals, data_type * data_vals, int i, int j, int k, int size_x, int size_y, int size_z, int used){

  data_type min;
  int u;
  
  min = data_vals[((k*size_x*size_y) + (j*size_x) + i)];

  // find minimum value
  for(u = 0; u < used; u++){
    
    if(temp_vals[u] <= min){ min = temp_vals[u]; }
    
  } 

  return min;

}

template <typename data_type>
data_type Temp_array_max(data_type * temp_vals, data_type * data_vals, int i, int j, int k, int size_x, int size_y, int size_z, int used){

  data_type max;
  int u;
  
  max = data_vals[((k*size_x*size_y) + (j*size_x) + i)];

  // find maximum value
  for(u = 0; u < used; u++){
    
    if(temp_vals[u] >= max){ max = temp_vals[u]; }
    
  } 

  return max;

}

template <typename data_type>
data_type Temp_array_median(data_type * temp_vals, data_type * data_vals, int i, int j, int k, int size_x, int size_y, int size_z, int used){
  
  // sort data
  HeapSort_voxels(used,temp_vals - 1);

  // return median value
  if((used % 2) == 0){

    return (0.5 * (temp_vals[(used/2) - 1] + temp_vals[used/2]));

  } else { 

    return temp_vals[(used - 1)/2];

  }

}

template <typename data_type>
void MinSmooth(data_type * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int axis, int round_flag, int vb_flag){

  data_type * buffer_vals, * temp_vals;
  float sep, sep_max, progress;
  int i,j,k,s,bdepth,used,NOvals;
  
  // create temp array to sort
  if(round_flag < 0){
    
    NOvals = ((2*sradius) + 1)*((2*sradius) + 1)*((2*fradius) + 1);
    temp_vals = new data_type[NOvals];
    
  } else if(round_flag == 0){

    NOvals = 0;

    for(j = (-1 * sradius); j <= sradius; ++j){
	
      for(i = (-1 * sradius); i <= sradius; ++i){
	  
	sep = sqrtf(((float)(i * i) + (float)(j * j)));
	if(sep <= ((float) sradius)){ NOvals++; }
	  
	// for(i = (-1 * sradius); i <= sradius; ++i)
      }
	
      // for(j = (-1 * sradius); j <= sradius; ++j)
    }

    NOvals*=((2*fradius) + 1);
    temp_vals = new data_type[NOvals];

  } else {
    
    NOvals = 0;
    
    for(k = (-1 * fradius); (k < ((-1 * fradius) + sradius)) && (k<0); ++k){
      
      // calculate projected radius for this channel
      sep_max = acos(((float)(k + fradius - sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; ++j){
	
	for(i = (-1 * sradius); i <= sradius; ++i){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }
      
      // for(k = (-1 * fradius); k <= fradius; ++k)
    }
    
    for(k = (-1 * fradius) + sradius; (k <= (fradius - sradius)) && (k<0); ++k){
      
      // calculate projected radius for this channel
      sep_max = (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; ++j){
	
	for(i = (-1 * sradius); i <= sradius; ++i){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }
      
      // for(k = (-1 * fradius); k <= fradius; ++k)
    }

    // calculate projected radius for this channel
    sep_max = (float) sradius;
    
    for(j = (-1 * sradius); j <= sradius; ++j){
      
      for(i = (-1 * sradius); i <= sradius; ++i){
	
	sep = sqrtf(((float)(i * i) + (float)(j * j)));
	if(sep <= sep_max){ NOvals++; }
	
	// for(i = (-1 * sradius); i <= sradius; ++i)
      }
      
      // for(j = (-1 * sradius); j <= sradius; ++j)
    }
    
    for(k = 1; k <= (fradius - sradius); ++k){
      
      // calculate projected radius for this channel
      sep_max = (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; ++j){
	
	for(i = (-1 * sradius); i <= sradius; ++i){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }
      
      // for(k = (-1 * fradius); k <= fradius; ++k)
    }
    
    for(k = (fradius - sradius + 1); (k <= fradius) && (k>0); ++k){
      
      // calculate projected radius for this channel
      sep_max = acos(((float)(k - fradius + sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; ++j){
	
	for(i = (-1 * sradius); i <= sradius; ++i){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }
      
      // for(k = (-1 * fradius); k <= fradius; ++k)
    }

    temp_vals = new data_type[NOvals];
    
  }
  if(vb_flag > 0){ cout << "Using temporary array containing " << NOvals << " values." << endl; }
  
  // carry out smoothing by shuffling data along longest axis
  // using the buffer array
  switch (axis) {
      
    case 1:
      // longest axis is ra
      
      // calculate buffer depth
      bdepth = sradius + 1;

      // create buffer array
      buffer_vals = new data_type[(size_y * size_z * bdepth)];
      
      // fill initial buffer values
      for(i = 0; i < (bdepth - 1); ++i){

	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate minimum of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	
	  // for(k = 0; k < size_z; ++k)
	}   

	// for(i = 0; i < (bdepth - 1); ++i)
      }
      cout << "done" << endl;
	      
      // buffer and shuffle minimum values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(i = (bdepth - 1); i < size_x; ++i){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate minimum of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_y) + (j*bdepth) + bdepth - 1)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){

	    data_vals[((k*size_x*size_y) + (j*size_x) + i - bdepth + 1)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth))];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s + 1)];
	      
	    // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	
	  // for(k = 0; k < size_z; ++k)
	}

	if(vb_flag > 0){ while(progress <= (((float) (i + 1)) / ((float) size_x))){ cout << "*"; cout.flush(); progress+=0.05; } }

	// for(i = (bdepth - 1); i < size_x; ++i)
      }
      if(vb_flag > 0){ cout << "* done." << endl; }
      
      // write final set of buffer values to array
      for(i = 0; (i < (bdepth - 1)); ++i){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    	      
	    data_vals[((k*size_x*size_y) + (j*size_x) + size_x - bdepth + 1 + i)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)];
	      	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}

	// (i = 0; (i < (bdepth - 1)); ++i)
      }

      // end of case statement
      break;
      
    case 2:
      // longest axis is dec
      
      // calculate buffer depth
      bdepth = sradius + 1;
      
      // create buffer array
      buffer_vals = new data_type[(size_x * size_z * bdepth)];
      
      // fill initial buffer values	
      for(j = 0; j < (bdepth - 1); ++j){

	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate minimum of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// for(j = 0; j < (bdepth - 1); ++j)
      }
  
      // buffer and shuffle minimum values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(j = (bdepth - 1); j < size_y; ++j){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate minimum of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_x) + ((bdepth - 1)*size_x) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    data_vals[((k*size_x*size_y) + ((j - bdepth + 1)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_x) + (s*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + ((s + 1)*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}
	
	if(vb_flag > 0){ while(progress <= (((float) (j + 1)) / ((float) size_y))){ cout << "*"; cout.flush(); progress+=0.05; } }
	
	// for(j = (bdepth - 1); j < size_y; ++j)
      }
      if(vb_flag > 0){ cout << "* done." << endl; }
      
      // write final set of buffer values to array
      for(j = 0; (j < (bdepth - 1)); ++j){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    	      
	    data_vals[((k*size_x*size_y) + ((size_y - bdepth + 1 + j)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)];
	      	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}

	// for(j = 0; (j < (bdepth - 1)); ++j)
      }
      
      // end of case statement
      break;
      
    case 3:
      // longest axis is frequency
      
      // calculate buffer depth
      bdepth = fradius + 1;

      // create buffer array
      buffer_vals = new data_type[(size_x * size_y * bdepth)];
      
      // fill initial buffer values
      for(k = 0; k < (bdepth - 1); ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate minimum of temp values and update buffer_vals
	    buffer_vals[((k*size_x*size_y) + (j*size_x) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(i = 0; i < size_x; ++i)
	}   
	
	// for(k = 0; k < (bdepth - 1); ++k)
      }  
      
      // buffer and shuffle minimum values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(k = (bdepth - 1); k < size_z; ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate minimum of temp values and update buffer_vals
	    buffer_vals[(((bdepth - 1)*size_x*size_y) + (j*size_x) + i)] = Temp_array_min(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(i = 0; i < size_x; ++i)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    data_vals[(((k - bdepth + 1)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((j*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((s*size_x*size_y) + (j*size_x) + i)] = buffer_vals[(((s + 1)*size_x*size_y) + (j*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(i = 0; i < size_x; ++i)
	}

	if(vb_flag > 0){ while(progress <= (((float) (k + 1)) / ((float) size_z))){ cout << "*"; cout.flush(); progress+=0.05; } }
	
	// for(k = (bdepth - 1); k < size_z; ++k)
      }  
      if(vb_flag > 0){ cout << "* done." << endl; }

      // write final set of buffer values to array
      for(k = 0; (k < (bdepth - 1)); ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    	      
	    data_vals[(((size_z - bdepth + 1 + k)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((k*size_x*size_y) + (j*size_x) + i)];
	      	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(i = 0; i < size_x; ++i)
	}

	// for(k = 0; (k < (bdepth - 1)); ++k)
      }
      
      // end of case statement
      break;
           
      // switch(axis) . . . case
  }
  
  // free up memory 
  delete [] temp_vals;
  delete [] buffer_vals;
  
}

template <typename data_type>
void MaxSmooth(data_type * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int axis, int round_flag, int vb_flag){

  data_type * buffer_vals, * temp_vals;
  float sep, sep_max, progress;
  int i,j,k,s,bdepth,used,NOvals;
  
  // create temp array to sort
  if(round_flag < 0){
    
    NOvals = ((2*sradius) + 1)*((2*sradius) + 1)*((2*fradius) + 1);
    temp_vals = new data_type[NOvals];
    
  } else if(round_flag == 0){

    NOvals = 0;

    for(j = (-1 * sradius); j <= sradius; ++j){
	
      for(i = (-1 * sradius); i <= sradius; ++i){
	  
	sep = sqrtf(((float)(i * i) + (float)(j * j)));
	if(sep <= ((float) sradius)){ NOvals++; }
	  
	// for(i = (-1 * sradius); i <= sradius; ++i)
      }
	
      // for(j = (-1 * sradius); j <= sradius; ++j)
    }

    NOvals*=((2*fradius) + 1);
    temp_vals = new data_type[NOvals];

  } else {
    
    NOvals = 0;
    
    for(k = (-1 * fradius); (k < ((-1 * fradius) + sradius)) && (k<0); ++k){
      
      // calculate projected radius for this channel
      sep_max = acos(((float)(k + fradius - sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; ++j){
	
	for(i = (-1 * sradius); i <= sradius; ++i){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }
      
      // for(k = (-1 * fradius); k <= fradius; ++k)
    }
    
    for(k = (-1 * fradius) + sradius; (k <= (fradius - sradius)) && (k<0); ++k){
      
      // calculate projected radius for this channel
      sep_max = (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; ++j){
	
	for(i = (-1 * sradius); i <= sradius; ++i){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }
      
      // for(k = (-1 * fradius); k <= fradius; ++k)
    }

    // calculate projected radius for this channel
    sep_max = (float) sradius;
    
    for(j = (-1 * sradius); j <= sradius; ++j){
      
      for(i = (-1 * sradius); i <= sradius; ++i){
	
	sep = sqrtf(((float)(i * i) + (float)(j * j)));
	if(sep <= sep_max){ NOvals++; }
	
	// for(i = (-1 * sradius); i <= sradius; ++i)
      }
      
      // for(j = (-1 * sradius); j <= sradius; ++j)
    }
    
    for(k = 1; k <= (fradius - sradius); ++k){
      
      // calculate projected radius for this channel
      sep_max = (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; ++j){
	
	for(i = (-1 * sradius); i <= sradius; ++i){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }
      
      // for(k = (-1 * fradius); k <= fradius; ++k)
    }
    
    for(k = (fradius - sradius + 1); (k <= fradius) && (k>0); ++k){
      
      // calculate projected radius for this channel
      sep_max = acos(((float)(k - fradius + sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; ++j){
	
	for(i = (-1 * sradius); i <= sradius; ++i){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }
      
      // for(k = (-1 * fradius); k <= fradius; ++k)
    }

    temp_vals = new data_type[NOvals];
    
  }
  if(vb_flag > 0){ cout << "Using temporary array containing " << NOvals << " values." << endl; }
  
  // carry out smoothing by shuffling data along longest axis
  // using the buffer array
  switch (axis)
    {
      
    case 1:
      // longest axis is ra
      
      // calculate buffer depth
      bdepth = sradius + 1;

      // create buffer array
      buffer_vals = new data_type[(size_y * size_z * bdepth)];
      
      // fill initial buffer values
      for(i = 0; i < (bdepth - 1); ++i){

	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate maximum of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)] = Temp_array_max(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	
	  // for(k = 0; k < size_z; ++k)
	}   

	// for(i = 0; i < (bdepth - 1); ++i)
      }
	      
      // buffer and shuffle maximum values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(i = (bdepth - 1); i < size_x; ++i){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate maximum of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_y) + (j*bdepth) + bdepth - 1)] = Temp_array_max(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){

	    data_vals[((k*size_x*size_y) + (j*size_x) + i - bdepth + 1)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth))];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s + 1)];
	      
	    // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	
	  // for(k = 0; k < size_z; ++k)
	}

	if(vb_flag > 0){ while(progress <= (((float) (i + 1)) / ((float) size_x))){ cout << "*"; cout.flush(); progress+=0.05; } }

	// for(i = (bdepth - 1); i < size_x; ++i)
      }
      if(vb_flag > 0){ cout << "* done." << endl; }
      
      // write final set of buffer values to array
      for(i = 0; (i < (bdepth - 1)); ++i){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    	      
	    data_vals[((k*size_x*size_y) + (j*size_x) + size_x - bdepth + 1 + i)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)];
	      	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}

	// (i = 0; (i < (bdepth - 1)); ++i)
      }

      // end of case statement
      break;
      
    case 2:
      // longest axis is dec
      
      // calculate buffer depth
      bdepth = sradius + 1;
      
      // create buffer array
      buffer_vals = new data_type[(size_x * size_z * bdepth)];
      
      // fill initial buffer values	
      for(j = 0; j < (bdepth - 1); ++j){

	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate maximum of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)] = Temp_array_max(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// for(j = 0; j < (bdepth - 1); ++j)
      }
  
      // buffer and shuffle maximum values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(j = (bdepth - 1); j < size_y; ++j){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate maximum of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_x) + ((bdepth - 1)*size_x) + i)] = Temp_array_max(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    data_vals[((k*size_x*size_y) + ((j - bdepth + 1)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_x) + (s*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + ((s + 1)*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}
	
	if(vb_flag > 0){ while(progress <= (((float) (j + 1)) / ((float) size_y))){ cout << "*"; cout.flush(); progress+=0.05; } }

	// for(j = (bdepth - 1); j < size_y; ++j)
      }
      if(vb_flag > 0){ cout << "* done." << endl; }
      
      // write final set of buffer values to array
      for(j = 0; (j < (bdepth - 1)); ++j){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    	      
	    data_vals[((k*size_x*size_y) + ((size_y - bdepth + 1 + j)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)];
	      	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}

	// for(j = 0; (j < (bdepth - 1)); ++j)
      }
      
      // end of case statement
      break;
      
    case 3:
      // longest axis is frequency
      
      // calculate buffer depth
      bdepth = fradius + 1;

      // create buffer array
      buffer_vals = new data_type[(size_x * size_y * bdepth)];
      
      // fill initial buffer values
      for(k = 0; k < (bdepth - 1); ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate maximum of temp values and update buffer_vals
	    buffer_vals[((k*size_x*size_y) + (j*size_x) + i)] = Temp_array_max(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(i = 0; i < size_x; ++i)
	}   
	
	// for(k = 0; k < (bdepth - 1); ++k)
      }  
      
      // buffer and shuffle maximum values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(k = (bdepth - 1); k < size_z; ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate maximum of temp values and update buffer_vals
	    buffer_vals[(((bdepth - 1)*size_x*size_y) + (j*size_x) + i)] = Temp_array_max(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(i = 0; i < size_x; ++i)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    data_vals[(((k - bdepth + 1)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((j*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((s*size_x*size_y) + (j*size_x) + i)] = buffer_vals[(((s + 1)*size_x*size_y) + (j*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(i = 0; i < size_x; ++i)
	}

	if(vb_flag > 0){ while(progress <= (((float) (k + 1)) / ((float) size_z))){ cout << "*"; cout.flush(); progress+=0.05; } }
	
	// for(k = (bdepth - 1); k < size_z; ++k)
      }  
      if(vb_flag > 0){ cout << "* done." << endl; }

      // write final set of buffer values to array
      for(k = 0; (k < (bdepth - 1)); ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    	      
	    data_vals[(((size_z - bdepth + 1 + k)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((k*size_x*size_y) + (j*size_x) + i)];
	      	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(i = 0; i < size_x; ++i)
	}

	// for(k = 0; (k < (bdepth - 1)); ++k)
      }
      
      // end of case statement
      break;
           
      // switch(axis) . . . case
    }
  
  // free up memory 
  delete [] temp_vals;
  delete [] buffer_vals;
  
}

template <typename data_type>
void MedianSmooth(data_type * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int axis, int round_flag, int vb_flag){

  data_type * buffer_vals, * temp_vals;
  float sep, sep_max, progress;
  int i,j,k,s,bdepth,used,NOvals;
  
  // create temp array to sort
  if(round_flag < 0){
    
    NOvals = ((2*sradius) + 1)*((2*sradius) + 1)*((2*fradius) + 1);
    temp_vals = new data_type[NOvals];
    
  } else if(round_flag == 0){

    NOvals = 0;

    for(j = (-1 * sradius); j <= sradius; ++j){
	
      for(i = (-1 * sradius); i <= sradius; ++i){
	  
	sep = sqrtf(((float)(i * i) + (float)(j * j)));
	if(sep <= ((float) sradius)){ NOvals++; }
	  
	// for(i = (-1 * sradius); i <= sradius; ++i)
      }
	
      // for(j = (-1 * sradius); j <= sradius; ++j)
    }

    NOvals*=((2*fradius) + 1);
    temp_vals = new data_type[NOvals];

  } else {
    
    NOvals = 0;
    
    for(k = (-1 * fradius); (k < (sradius - fradius)) && (k<0); ++k){
      
      // calculate projected radius for this channel
      sep_max = acos(((float)(k + fradius - sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; ++j){
	
	for(i = (-1 * sradius); i <= sradius; ++i){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }
      
      // for(k = (-1 * fradius); k <= fradius; ++k)
    }
    
    for(k = sradius - fradius; (k <= (fradius - sradius)) && (k<0); ++k){
      
      // calculate projected radius for this channel
      sep_max = (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; ++j){
	
	for(i = (-1 * sradius); i <= sradius; ++i){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }
      
      // for(k = (-1 * fradius); k <= fradius; ++k)
    }

    // calculate projected radius for this channel
    sep_max = (float) sradius;
    
    for(j = (-1 * sradius); j <= sradius; ++j){
      
      for(i = (-1 * sradius); i <= sradius; ++i){
	
	sep = sqrtf(((float)(i * i) + (float)(j * j)));
	if(sep <= sep_max){ NOvals++; }
	
	// for(i = (-1 * sradius); i <= sradius; ++i)
      }
      
      // for(j = (-1 * sradius); j <= sradius; ++j)
    }

    for(k = 1; k <= (fradius - sradius); ++k){
      
      // calculate projected radius for this channel
      sep_max = (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; ++j){
	
	for(i = (-1 * sradius); i <= sradius; ++i){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }
      
      // for(k = (-1 * fradius); k <= fradius; ++k)
    }
    
    for(k = (fradius - sradius + 1); (k <= fradius) && (k>0); ++k){
      
      // calculate projected radius for this channel
      sep_max = acos(((float)(k - fradius + sradius) / (float)(sradius)));
      sep_max = sin(sep_max) * (float) sradius;
      
      for(j = (-1 * sradius); j <= sradius; ++j){
	
	for(i = (-1 * sradius); i <= sradius; ++i){
	  
	  sep = sqrtf(((float)(i * i) + (float)(j * j)));
	  if(sep <= sep_max){ NOvals++; }
	  
	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }
      
      // for(k = (-1 * fradius); k <= fradius; ++k)
    }

    temp_vals = new data_type[NOvals];
    
  }
  if(vb_flag > 0){ cout << "Using temporary array containing " << NOvals << " values." << endl; }
  
  // carry out smoothing by shuffling data along longest axis
  // using the buffer array
  switch (axis) {
      
    case 1:
      // longest axis is ra
      
      // calculate buffer depth
      bdepth = sradius + 1;

      // create buffer array
      buffer_vals = new data_type[(size_y * size_z * bdepth)];
      
      // fill initial buffer values
      for(i = 0; i < (bdepth - 1); ++i){

	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)] = Temp_array_median(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	
	  // for(k = 0; k < size_z; ++k)
	}   

	// for(i = 0; i < (bdepth - 1); ++i)
      }
      cout << "done" << endl;
	      
      // buffer and shuffle median values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(i = (bdepth - 1); i < size_x; ++i){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_y) + (j*bdepth) + bdepth - 1)] = Temp_array_median(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){

	    data_vals[((k*size_x*size_y) + (j*size_x) + i - bdepth + 1)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth))];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s + 1)];
	      
	    // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	
	  // for(k = 0; k < size_z; ++k)
	}

	if(vb_flag > 0){ while(progress <= (((float) (i + 1)) / ((float) size_x))){ cout << "*"; cout.flush(); progress+=0.05; } }

	// for(i = (bdepth - 1); i < size_x; ++i)
      }
      if(vb_flag > 0){ cout << "* done." << endl; }
      
      // write final set of buffer values to array
      for(i = 0; (i < (bdepth - 1)); ++i){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    	      
	    data_vals[((k*size_x*size_y) + (j*size_x) + size_x - bdepth + 1 + i)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)];
	      	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}

	// (i = 0; (i < (bdepth - 1)); ++i)
      }

      // end of case statement
      break;
      
    case 2:
      // longest axis is dec
      
      // calculate buffer depth
      bdepth = sradius + 1;
      
      // create buffer array
      buffer_vals = new data_type[(size_x * size_z * bdepth)];
      
      // fill initial buffer values	
      for(j = 0; j < (bdepth - 1); ++j){

	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)] = Temp_array_median(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// for(j = 0; j < (bdepth - 1); ++j)
      }
  
      // buffer and shuffle median values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(j = (bdepth - 1); j < size_y; ++j){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*bdepth*size_x) + ((bdepth - 1)*size_x) + i)] = Temp_array_median(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    data_vals[((k*size_x*size_y) + ((j - bdepth + 1)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_x) + (s*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + ((s + 1)*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}
	
	if(vb_flag > 0){ while(progress <= (((float) (j + 1)) / ((float) size_y))){ cout << "*"; cout.flush(); progress+=0.05; } }
	
	// for(j = (bdepth - 1); j < size_y; ++j)
      }
      if(vb_flag > 0){ cout << "* done." << endl; }
      
      // write final set of buffer values to array
      for(j = 0; (j < (bdepth - 1)); ++j){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    	      
	    data_vals[((k*size_x*size_y) + ((size_y - bdepth + 1 + j)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)];
	      	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}

	// for(j = 0; (j < (bdepth - 1)); ++j)
      }
      
      // end of case statement
      break;
      
    case 3:
      // longest axis is frequency
      
      // calculate buffer depth
      bdepth = fradius + 1;

      // create buffer array
      buffer_vals = new data_type[(size_x * size_y * bdepth)];
      
      // fill initial buffer values
      for(k = 0; k < (bdepth - 1); ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[((k*size_x*size_y) + (j*size_x) + i)] = Temp_array_median(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(i = 0; i < size_x; ++i)
	}   
	
	// for(k = 0; k < (bdepth - 1); ++k)
      }  
      
      // buffer and shuffle median values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(k = (bdepth - 1); k < size_z; ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // populate array of temp values
	    used = Populate_temp_array(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,sradius,fradius,round_flag);
	    
	    // calculate median of temp values and update buffer_vals
	    buffer_vals[(((bdepth - 1)*size_x*size_y) + (j*size_x) + i)] = Temp_array_median(temp_vals,data_vals,i,j,k,size_x,size_y,size_z,used);
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(i = 0; i < size_x; ++i)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    data_vals[(((k - bdepth + 1)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((j*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((s*size_x*size_y) + (j*size_x) + i)] = buffer_vals[(((s + 1)*size_x*size_y) + (j*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(i = 0; i < size_x; ++i)
	}

	if(vb_flag > 0){ while(progress <= (((float) (k + 1)) / ((float) size_z))){ cout << "*"; cout.flush(); progress+=0.05; } }
	
	// for(k = (bdepth - 1); k < size_z; ++k)
      }  
      if(vb_flag > 0){ cout << "* done." << endl; }

      // write final set of buffer values to array
      for(k = 0; (k < (bdepth - 1)); ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    	      
	    data_vals[(((size_z - bdepth + 1 + k)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((k*size_x*size_y) + (j*size_x) + i)];
	      	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(i = 0; i < size_x; ++i)
	}

	// for(k = 0; (k < (bdepth - 1)); ++k)
      }
      
      // end of case statement
      break;
           
      // switch(axis) . . . case
    }
  
  // free up memory 
  delete [] temp_vals;
  delete [] buffer_vals;
  
}

template <typename data_type>
void MathMorphOpen(data_type * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag, int vb_flag){

  int axis;

  // determine the longest axis
  axis = 1;
  if((size_y >= size_x) && (size_y > size_z)){ axis = 2; }
  if((size_z >= size_x) && (size_z >= size_y)){ axis = 3; }

  if(vb_flag > 0){ cout << "Smoothing along axis: " << axis << endl; }

  // carry out the first smoothing iteration
  if(vb_flag > 0){ cout << "Carrying out 1st smoothing iteration . . . " << endl; }
  MinSmooth(data_vals,size_x,size_y,size_z,sradius,fradius,axis,round_flag,vb_flag);

  // carry out the second smoothing iteration
  if(vb_flag > 0){ cout << "Carrying out 2nd smoothing iteration . . . " << endl; }
  MaxSmooth(data_vals,size_x,size_y,size_z,sradius,fradius,axis,round_flag,vb_flag);

}

template <typename data_type>
void MathMorphClose(data_type * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag, int vb_flag){

  int axis;

  // determine the longest axis
  axis = 1;
  if((size_y >= size_x) && (size_y > size_z)){ axis = 2; }
  if((size_z >= size_x) && (size_z >= size_y)){ axis = 3; }

  if(vb_flag > 0){ cout << "Smoothing along axis: " << axis << endl; }

  // carry out the first smoothing iteration
  if(vb_flag > 0){ cout << "Carrying out 1st smoothing iteration . . . " << endl; }
  MaxSmooth(data_vals,size_x,size_y,size_z,sradius,fradius,axis,round_flag,vb_flag);

  // carry out the second smoothing iteration
  if(vb_flag > 0){ cout << "Carrying out 2nd smoothing iteration . . . " << endl; }
  MinSmooth(data_vals,size_x,size_y,size_z,sradius,fradius,axis,round_flag,vb_flag);

}

template <typename data_type>
void RemHitMiss_4conn(data_type * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, short unsigned int curr_flag_val, int vb_flag){

  data_type * buffer_vals;
  float progress;
  int test_flag,axis,i,j,k,s,bdepth;
  
  // range test sradius and fradius
  sradius = (sradius > 1) ? 1 : sradius;
  sradius = (sradius < 0) ? 0 : sradius;
  fradius = (fradius > 1) ? 1 : fradius;
  fradius = (fradius < 0) ? 0 : fradius;

  // determine the longest axis
  axis = 1;
  if((size_y >= size_x) && (size_y > size_z)){ axis = 2; }
  if((size_z >= size_x) && (size_z >= size_y)){ axis = 3; }

  if(vb_flag > 0){ 
    cout << "Transforming along axis: " << axis << endl;
    cout << "Removing cross hit-and-miss transform . . . " << endl;
  }      

  // carry out smoothing by shuffling data along longest axis
  // using the buffer array
  switch (axis)
    {  
    case 1:
      // longest axis is ra
      
      // calculate buffer depth
      bdepth = sradius + 1;

      // create buffer array
      buffer_vals = new data_type[(size_y * size_z * bdepth)];
      
      // fill initial buffer values
      for(i = 0; i < (bdepth - 1); ++i){

	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // apply cross hit-and-miss transform
	    if(data_vals[((k * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){

	      test_flag = 1;
	      if(sradius > 0){

		if(i > 0){
		  
		  if(data_vals[((k * size_x * size_y) + (j * size_x) + i - 1)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(i < (size_x - 1)){
		  
		  if(data_vals[((k * size_x * size_y) + (j * size_x) + i + 1)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(j > 0){
		  
		  if(data_vals[((k * size_x * size_y) + ((j - 1) * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(j < (size_y - 1)){
		  
		  if(data_vals[((k * size_x * size_y) + ((j + 1) * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}

	      }
	      if(fradius > 0){

		if(k > 0){
		
		  if(data_vals[(((k - 1) * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(k < (size_z - 1)){
		  
		  if(data_vals[(((k + 1) * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}

	      }
	      
	      // update buffer_vals
	      if(test_flag > 0){ 

		buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)] = 0;

	      } else {

		buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)] = data_vals[((k * size_x * size_y) + (j * size_x) + i)];

	      }
	      
	      // if(temp_data_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val)
	    }

	    // for(j = 0; j < size_y; ++j)
	  }
	
	  // for(k = 0; k < size_z; ++k)
	}   

	// for(i = 0; i < (bdepth - 1); ++i)
      }
	      
      // buffer and shuffle values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(i = (bdepth - 1); i < size_x; ++i){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // apply cross hit-and-miss transform
	    if(data_vals[((k * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){

	      test_flag = 1;

	      if(sradius > 0){

		if(i > 0){
		
		  if(data_vals[((k * size_x * size_y) + (j * size_x) + i - 1)] == curr_flag_val){ test_flag = -1; }
		
		}
		if(i < (size_x - 1)){
		  
		  if(data_vals[((k * size_x * size_y) + (j * size_x) + i + 1)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(j > 0){
		  
		  if(data_vals[((k * size_x * size_y) + ((j - 1) * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(j < (size_y - 1)){
		  
		  if(data_vals[((k * size_x * size_y) + ((j + 1) * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}

	      }
	      if(fradius > 0){

		if(k > 0){
		  
		  if(data_vals[(((k - 1) * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(k < (size_z - 1)){
		  
		  if(data_vals[(((k + 1) * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}

	      }
	      
	      // update buffer_vals
	      if(test_flag > 0){ 

		buffer_vals[((k*bdepth*size_y) + (j*bdepth) + bdepth - 1)] = 0;

	      } else {

		buffer_vals[((k*bdepth*size_y) + (j*bdepth) + bdepth - 1)] = data_vals[((k * size_x * size_y) + (j * size_x) + i)];

	      }
	      
	      // if(temp_data_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val)
	    }
	   
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){

	    data_vals[((k*size_x*size_y) + (j*size_x) + i - bdepth + 1)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth))];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s + 1)];
	      
	    // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	
	  // for(k = 0; k < size_z; ++k)
	}

	if(vb_flag > 0){ while(progress <= (((float) (i + 1)) / ((float) size_x))){ cout << "*"; cout.flush(); progress+=0.05; } }

	// for(i = (bdepth - 1); i < size_x; ++i)
      }
      if(vb_flag > 0){ cout << "* done." << endl; }
      
      // write final set of buffer values to array
      for(i = 0; (i < (bdepth - 1)); ++i){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    	      
	    data_vals[((k*size_x*size_y) + (j*size_x) + size_x - bdepth + 1 + i)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)];
	      	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}

	// (i = 0; (i < (bdepth - 1)); ++i)
      }

      // end of case statement
      break;
      
    case 2:
      // longest axis is dec
      
      // calculate buffer depth
      bdepth = sradius + 1;
      
      // create buffer array
      buffer_vals = new data_type[(size_x * size_z * bdepth)];
      
      // fill initial buffer values	
      for(j = 0; j < (bdepth - 1); ++j){

	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    // apply cross hit-and-miss transform
	    if(data_vals[((k * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){

	      test_flag = 1;

	      if(sradius > 0){

		if(i > 0){
		
		  if(data_vals[((k * size_x * size_y) + (j * size_x) + i - 1)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(i < (size_x - 1)){
		  
		  if(data_vals[((k * size_x * size_y) + (j * size_x) + i + 1)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(j > 0){
		
		  if(data_vals[((k * size_x * size_y) + ((j - 1) * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(j < (size_y - 1)){
		
		  if(data_vals[((k * size_x * size_y) + ((j + 1) * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		
		}

	      }
	      if(fradius > 0){
	      
		if(k > 0){
		
		  if(data_vals[(((k - 1) * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		
		}
		if(k < (size_z - 1)){
		
		  if(data_vals[(((k + 1) * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		
		}

	      }
	      
	      // update buffer_vals
	      if(test_flag > 0){ 

		buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)] = 0;

	      } else {

		buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)] = data_vals[((k * size_x * size_y) + (j * size_x) + i)];

	      }
	      
	      // if(temp_data_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val)
	    }

	    // for(i = 0; i < size_x; ++i)
	  }
	
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// for(j = 0; j < (bdepth - 1); ++j)
      }
  
      // buffer and shuffle median values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(j = (bdepth - 1); j < size_y; ++j){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    // apply cross hit-and-miss transform
	    if(data_vals[((k * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){

	      test_flag = 1;
	      if(sradius > 0){

		if(i > 0){
		
		  if(data_vals[((k * size_x * size_y) + (j * size_x) + i - 1)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(i < (size_x - 1)){
		  
		  if(data_vals[((k * size_x * size_y) + (j * size_x) + i + 1)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(j > 0){
		  
		  if(data_vals[((k * size_x * size_y) + ((j - 1) * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(j < (size_y - 1)){
		  
		  if(data_vals[((k * size_x * size_y) + ((j + 1) * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}

	      }
	      if(fradius > 0){

		if(k > 0){
		  
		  if(data_vals[(((k - 1) * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(k < (size_z - 1)){
		  
		  if(data_vals[(((k + 1) * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}
		
	      }

	      // update buffer_vals
	      if(test_flag > 0){ 
		
		buffer_vals[((k*bdepth*size_x) + ((bdepth - 1)*size_x) + i)] = 0;

	      } else {

		buffer_vals[((k*bdepth*size_x) + ((bdepth - 1)*size_x) + i)] = data_vals[((k * size_x * size_y) + (j * size_x) + i)];

	      }
	      
	      // if(temp_data_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val)
	    }
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    data_vals[((k*size_x*size_y) + ((j - bdepth + 1)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_x) + (s*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + ((s + 1)*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}
	
	if(vb_flag > 0){ while(progress <= (((float) (j + 1)) / ((float) size_y))){ cout << "*"; cout.flush(); progress+=0.05; } }

	// for(j = (bdepth - 1); j < size_y; ++j)
      }
      if(vb_flag > 0){ cout << "* done." << endl; }
       
      // write final set of buffer values to array
      for(j = 0; (j < (bdepth - 1)); ++j){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    	      
	    data_vals[((k*size_x*size_y) + ((size_y - bdepth + 1 + j)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)];
	      	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}

	// for(j = 0; (j < (bdepth - 1)); ++j)
      }
      
      // end of case statement
      break;
      
    case 3:
      // longest axis is frequency
      
      // calculate buffer depth
      bdepth = fradius + 1;

      // create buffer array
      buffer_vals = new data_type[(size_x * size_y * bdepth)];
      
      // fill initial buffer values
      for(k = 0; k < (bdepth - 1); ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // apply cross hit-and-miss transform
	    if(data_vals[((k * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){

	      test_flag = 1;

	      if(sradius > 0){

		if(i > 0){
		
		  if(data_vals[((k * size_x * size_y) + (j * size_x) + i - 1)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(i < (size_x - 1)){
		  
		  if(data_vals[((k * size_x * size_y) + (j * size_x) + i + 1)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(j > 0){
		  
		  if(data_vals[((k * size_x * size_y) + ((j - 1) * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(j < (size_y - 1)){
		  
		  if(data_vals[((k * size_x * size_y) + ((j + 1) * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}

	      }
	      if(fradius > 0){

		if(k > 0){
		  
		  if(data_vals[(((k - 1) * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(k < (size_z - 1)){
		  
		  if(data_vals[(((k + 1) * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}

	      }
	      
	      // update buffer_vals
	      if(test_flag > 0){ 

		buffer_vals[((k*size_x*size_y) + (j*size_x) + i)] = 0;

	      } else {

		buffer_vals[((k*size_x*size_y) + (j*size_x) + i)] = data_vals[((k * size_x * size_y) + (j * size_x) + i)];

	      }
	      
	      // if(temp_data_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val)
	    }
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(i = 0; i < size_x; ++i)
	}   
	
	// for(k = 0; k < (bdepth - 1); ++k)
      }  
      
      // buffer and shuffle median values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(k = (bdepth - 1); k < size_z; ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // apply cross hit-and-miss transform
	    if(data_vals[((k * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){

	      test_flag = 1;
	      
	      if(sradius > 0){

		if(i > 0){
		  
		  if(data_vals[((k * size_x * size_y) + (j * size_x) + i - 1)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(i < (size_x - 1)){
		  
		  if(data_vals[((k * size_x * size_y) + (j * size_x) + i + 1)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(j > 0){
		  
		  if(data_vals[((k * size_x * size_y) + ((j - 1) * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(j < (size_y - 1)){
		  
		  if(data_vals[((k * size_x * size_y) + ((j + 1) * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}

	      }
	      if(fradius > 0){

		if(k > 0){
		  
		  if(data_vals[(((k - 1) * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}
		if(k < (size_z - 1)){
		  
		  if(data_vals[(((k + 1) * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){ test_flag = -1; }
		  
		}

	      }

	      // update buffer_vals
	      if(test_flag > 0){ 

		buffer_vals[(((bdepth - 1)*size_x*size_y) + (j*size_x) + i)] = 0;

	      } else {

		buffer_vals[(((bdepth - 1)*size_x*size_y) + (j*size_x) + i)] = data_vals[((k * size_x * size_y) + (j * size_x) + i)];

	      }
	      
	      // if(temp_data_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val)
	    }

	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(i = 0; i < size_x; ++i)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    data_vals[(((k - bdepth + 1)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((j*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((s*size_x*size_y) + (j*size_x) + i)] = buffer_vals[(((s + 1)*size_x*size_y) + (j*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(i = 0; i < size_x; ++i)
	}

	if(vb_flag > 0){ while(progress <= (((float) (k + 1)) / ((float) size_z))){ cout << "*"; cout.flush(); progress+=0.05; } }
	
	// for(k = (bdepth - 1); k < size_z; ++k)
      }  
      if(vb_flag > 0){ cout << "* done." << endl; }

      // write final set of buffer values to array
      for(k = 0; (k < (bdepth - 1)); ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    	      
	    data_vals[(((size_z - bdepth + 1 + k)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((k*size_x*size_y) + (j*size_x) + i)];
	      	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(i = 0; i < size_x; ++i)
	}

	// for(k = 0; (k < (bdepth - 1)); ++k)
      }
      
      // end of case statement
      break;
           
      // switch(axis) . . . case
    }
  
  // free up memory 
  delete [] buffer_vals;

}

template <typename data_type>
void RemHitMiss_8conn(data_type * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, short unsigned int curr_flag_val, int vb_flag){

  data_type * buffer_vals;
  float progress;
  int test_flag,axis,i,j,k,s,bdepth;
  int i2,j2,k2,istart,ifinish,jstart,jfinish,kstart,kfinish;

  // range test sradius and fradius
  sradius = (sradius > 1) ? 1 : sradius;
  sradius = (sradius < 0) ? 0 : sradius;
  fradius = (fradius > 1) ? 1 : fradius;
  fradius = (fradius < 0) ? 0 : fradius;

  // determine the longest axis
  axis = 1;
  if((size_y >= size_x) && (size_y > size_z)){ axis = 2; }
  if((size_z >= size_x) && (size_z >= size_y)){ axis = 3; }

  if(vb_flag > 0){
    cout << "Transforming along axis: " << axis << endl;
    cout << "Removing cross hit-and-miss transform . . . " << endl;
  }      

  // carry out smoothing by shuffling data along longest axis
  // using the buffer array
  switch (axis)
    {  
    case 1:
      // longest axis is ra
      
      // calculate buffer depth
      bdepth = sradius + 1;

      // create buffer array
      buffer_vals = new data_type[(size_y * size_z * bdepth)];
      
      // fill initial buffer values
      for(i = 0; i < (bdepth - 1); ++i){

	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // apply cross hit-and-miss transform
	    if(data_vals[((k * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){

	      test_flag = 1;

	      istart = (i > 0) ? (i - sradius) : i;
	      ifinish = (i < (size_x - 1)) ? (i + sradius) : i;
	      jstart = (j > 0) ? (j - sradius) : j;
	      jfinish = (j < (size_y- 1)) ? (j + sradius) : j;
	      kstart = (k > 0) ? (k - fradius) : k;
	      kfinish = (k < (size_z - 1)) ? (k + fradius) : k;

	      for(k2 = kstart; k2 <= kfinish; k2++){
		for(j2 = jstart; j2 <= jfinish; j2++){
		  for(i2 = istart; i2 <= ifinish; i2++){

		    if((i2 == i) && (j2 == j) && (k2 == k)){ 
		      continue; 
		    } else {
		      if(data_vals[((k2 * size_x * size_y) + (j2 * size_x) + i2)] == curr_flag_val){ test_flag = -1; }
		    }

		  }
		}
	      }
		    
	      // update buffer_vals
	      if(test_flag > 0){ 

		buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)] = 0;

	      } else {

		buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)] = data_vals[((k * size_x * size_y) + (j * size_x) + i)];

	      }
	      
	      // if(temp_data_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val)
	    }

	    // for(j = 0; j < size_y; ++j)
	  }
	
	  // for(k = 0; k < size_z; ++k)
	}   

	// for(i = 0; i < (bdepth - 1); ++i)
      }
	      
      // buffer and shuffle values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(i = (bdepth - 1); i < size_x; ++i){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // apply cross hit-and-miss transform
	    if(data_vals[((k * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){

	      test_flag = 1;

	      istart = (i > 0) ? (i - sradius) : i;
	      ifinish = (i < (size_x - 1)) ? (i + sradius) : i;
	      jstart = (j > 0) ? (j - sradius) : j;
	      jfinish = (j < (size_y- 1)) ? (j + sradius) : j;
	      kstart = (k > 0) ? (k - fradius) : k;
	      kfinish = (k < (size_z - 1)) ? (k + fradius) : k;

	      for(k2 = kstart; k2 <= kfinish; k2++){
		for(j2 = jstart; j2 <= jfinish; j2++){
		  for(i2 = istart; i2 <= ifinish; i2++){

		    if((i2 == i) && (j2 == j) && (k2 == k)){ 
		      continue; 
		    } else {
		      if(data_vals[((k2 * size_x * size_y) + (j2 * size_x) + i2)] == curr_flag_val){ test_flag = -1; }
		    }

		  }
		}
	      }
	      
	      // update buffer_vals
	      if(test_flag > 0){ 

		buffer_vals[((k*bdepth*size_y) + (j*bdepth) + bdepth - 1)] = 0;

	      } else {

		buffer_vals[((k*bdepth*size_y) + (j*bdepth) + bdepth - 1)] = data_vals[((k * size_x * size_y) + (j * size_x) + i)];

	      }
	      
	      // if(temp_data_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val)
	    }
	   
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){

	    data_vals[((k*size_x*size_y) + (j*size_x) + i - bdepth + 1)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth))];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + s + 1)];
	      
	    // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	
	  // for(k = 0; k < size_z; ++k)
	}

	if(vb_flag > 0){ while(progress <= (((float) (i + 1)) / ((float) size_x))){ cout << "*"; cout.flush(); progress+=0.05; } }

	// for(i = (bdepth - 1); i < size_x; ++i)
      }
      if(vb_flag > 0){ cout << "* done." << endl; }
      
      // write final set of buffer values to array
      for(i = 0; (i < (bdepth - 1)); ++i){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(j = 0; j < size_y; ++j){
	    	      
	    data_vals[((k*size_x*size_y) + (j*size_x) + size_x - bdepth + 1 + i)] = buffer_vals[((k*bdepth*size_y) + (j*bdepth) + i)];
	      	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}

	// (i = 0; (i < (bdepth - 1)); ++i)
      }

      // end of case statement
      break;
      
    case 2:
      // longest axis is dec
      
      // calculate buffer depth
      bdepth = sradius + 1;
      
      // create buffer array
      buffer_vals = new data_type[(size_x * size_z * bdepth)];
      
      // fill initial buffer values	
      for(j = 0; j < (bdepth - 1); ++j){

	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    // apply cross hit-and-miss transform
	    if(data_vals[((k * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){

	      test_flag = 1;

	      istart = (i > 0) ? (i - sradius) : i;
	      ifinish = (i < (size_x - 1)) ? (i + sradius) : i;
	      jstart = (j > 0) ? (j - sradius) : j;
	      jfinish = (j < (size_y- 1)) ? (j + sradius) : j;
	      kstart = (k > 0) ? (k - fradius) : k;
	      kfinish = (k < (size_z - 1)) ? (k + fradius) : k;

	      for(k2 = kstart; k2 <= kfinish; k2++){
		for(j2 = jstart; j2 <= jfinish; j2++){
		  for(i2 = istart; i2 <= ifinish; i2++){

		    if((i2 == i) && (j2 == j) && (k2 == k)){ 
		      continue; 
		    } else {
		      if(data_vals[((k2 * size_x * size_y) + (j2 * size_x) + i2)] == curr_flag_val){ test_flag = -1; }
		    }

		  }
		}
	      }
	      
	      // update buffer_vals
	      if(test_flag > 0){ 

		buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)] = 0;

	      } else {

		buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)] = data_vals[((k * size_x * size_y) + (j * size_x) + i)];

	      }
	      
	      // if(temp_data_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val)
	    }

	    // for(i = 0; i < size_x; ++i)
	  }
	
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// for(j = 0; j < (bdepth - 1); ++j)
      }
  
      // buffer and shuffle median values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(j = (bdepth - 1); j < size_y; ++j){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    // apply cross hit-and-miss transform
	    if(data_vals[((k * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){

	      test_flag = 1;

	      istart = (i > 0) ? (i - sradius) : i;
	      ifinish = (i < (size_x - 1)) ? (i + sradius) : i;
	      jstart = (j > 0) ? (j - sradius) : j;
	      jfinish = (j < (size_y- 1)) ? (j + sradius) : j;
	      kstart = (k > 0) ? (k - fradius) : k;
	      kfinish = (k < (size_z - 1)) ? (k + fradius) : k;

	      for(k2 = kstart; k2 <= kfinish; k2++){
		for(j2 = jstart; j2 <= jfinish; j2++){
		  for(i2 = istart; i2 <= ifinish; i2++){

		    if((i2 == i) && (j2 == j) && (k2 == k)){ 
		      continue; 
		    } else {
		      if(data_vals[((k2 * size_x * size_y) + (j2 * size_x) + i2)] == curr_flag_val){ test_flag = -1; }
		    }

		  }
		}
	      }
	      
	      // update buffer_vals
	      if(test_flag > 0){ 

		buffer_vals[((k*bdepth*size_x) + ((bdepth - 1)*size_x) + i)] = 0;

	      } else {

		buffer_vals[((k*bdepth*size_x) + ((bdepth - 1)*size_x) + i)] = data_vals[((k * size_x * size_y) + (j * size_x) + i)];

	      }
	      
	      // if(temp_data_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val)
	    }
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  //for(k = 0; k < size_z; ++k)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    
	    data_vals[((k*size_x*size_y) + ((j - bdepth + 1)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((k*bdepth*size_x) + (s*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + ((s + 1)*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}
	
	if(vb_flag > 0){ while(progress <= (((float) (j + 1)) / ((float) size_y))){ cout << "*"; cout.flush(); progress+=0.05; } }

	// for(j = (bdepth - 1); j < size_y; ++j)
      }
      if(vb_flag > 0){ cout << "* done." << endl; }
      
      // write final set of buffer values to array
      for(j = 0; (j < (bdepth - 1)); ++j){
	
	for(k = 0; k < size_z; ++k){
	  
	  for(i = 0; i < size_x; ++i){
	    	      
	    data_vals[((k*size_x*size_y) + ((size_y - bdepth + 1 + j)*size_x) + i)] = buffer_vals[((k*bdepth*size_x) + (j*size_x) + i)];
	      	    
	    // for(i = 0; i < size_x; ++i)
	  }
	  
	  // for(k = 0; k < size_z; ++k)
	}

	// for(j = 0; (j < (bdepth - 1)); ++j)
      }
      
      // end of case statement
      break;
      
    case 3:
      // longest axis is frequency
      
      // calculate buffer depth
      bdepth = fradius + 1;

      // create buffer array
      buffer_vals = new data_type[(size_x * size_y * bdepth)];
      
      // fill initial buffer values
      for(k = 0; k < (bdepth - 1); ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // apply cross hit-and-miss transform
	    if(data_vals[((k * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){

	      test_flag = 1;
	
	      istart = (i > 0) ? (i - sradius) : i;
	      ifinish = (i < (size_x - 1)) ? (i + sradius) : i;
	      jstart = (j > 0) ? (j - sradius) : j;
	      jfinish = (j < (size_y- 1)) ? (j + sradius) : j;
	      kstart = (k > 0) ? (k - fradius) : k;
	      kfinish = (k < (size_z - 1)) ? (k + fradius) : k;

	      for(k2 = kstart; k2 <= kfinish; k2++){
		for(j2 = jstart; j2 <= jfinish; j2++){
		  for(i2 = istart; i2 <= ifinish; i2++){

		    if((i2 == i) && (j2 == j) && (k2 == k)){ 
		      continue; 
		    } else {
		      if(data_vals[((k2 * size_x * size_y) + (j2 * size_x) + i2)] == curr_flag_val){ test_flag = -1; }
		    }

		  }
		}
	      }

	      // update buffer_vals
	      if(test_flag > 0){ 

		buffer_vals[((k*size_x*size_y) + (j*size_x) + i)] = 0;

	      } else {

		buffer_vals[((k*size_x*size_y) + (j*size_x) + i)] = data_vals[((k * size_x * size_y) + (j * size_x) + i)];

	      }
	      
	      // if(temp_data_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val)
	    }
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(i = 0; i < size_x; ++i)
	}   
	
	// for(k = 0; k < (bdepth - 1); ++k)
      }  
      
      // buffer and shuffle median values
      if(vb_flag > 0){ cout << "0 | |:| | : | |:| | 100% complete" << endl; }
      progress = 0.0;
      for(k = (bdepth - 1); k < size_z; ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    // apply cross hit-and-miss transform
	    if(data_vals[((k * size_x * size_y) + (j * size_x) + i)] == curr_flag_val){

	      test_flag = 1;

	      istart = (i > 0) ? (i - sradius) : i;
	      ifinish = (i < (size_x - 1)) ? (i + sradius) : i;
	      jstart = (j > 0) ? (j - sradius) : j;
	      jfinish = (j < (size_y- 1)) ? (j + sradius) : j;
	      kstart = (k > 0) ? (k - fradius) : k;
	      kfinish = (k < (size_z - 1)) ? (k + fradius) : k;

	      for(k2 = kstart; k2 <= kfinish; k2++){
		for(j2 = jstart; j2 <= jfinish; j2++){
		  for(i2 = istart; i2 <= ifinish; i2++){

		    if((i2 == i) && (j2 == j) && (k2 == k)){ 
		      continue; 
		    } else {
		      if(data_vals[((k2 * size_x * size_y) + (j2 * size_x) + i2)] == curr_flag_val){ test_flag = -1; }
		    }

		  }
		}
	      }

	      // update buffer_vals
	      if(test_flag > 0){ 

		buffer_vals[(((bdepth - 1)*size_x*size_y) + (j*size_x) + i)] = 0;

	      } else {

		buffer_vals[(((bdepth - 1)*size_x*size_y) + (j*size_x) + i)] = data_vals[((k * size_x * size_y) + (j * size_x) + i)];

	      }
	      
	      // if(temp_data_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val)
	    }

	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  //for(i = 0; i < size_x; ++i)
	}   
	
	// write redundant set of buffer vals to proc_vals, and shuffle
	// buffer vals along
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    
	    data_vals[(((k - bdepth + 1)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((j*size_x) + i)];
	    
	    for(s = 0; s < (bdepth - 1); s++){
	      
	      buffer_vals[((s*size_x*size_y) + (j*size_x) + i)] = buffer_vals[(((s + 1)*size_x*size_y) + (j*size_x) + i)];
	      
	      // for(s = 0; s < (bdepth - 1); s++)
	    }
	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(i = 0; i < size_x; ++i)
	}

	if(vb_flag > 0){ while(progress <= (((float) (k + 1)) / ((float) size_z))){ cout << "*"; cout.flush(); progress+=0.05; } }
	
	// for(k = (bdepth - 1); k < size_z; ++k)
      }  
      if(vb_flag > 0){ cout << "* done." << endl; }

      // write final set of buffer values to array
      for(k = 0; (k < (bdepth - 1)); ++k){
	
	for(i = 0; i < size_x; ++i){
	  
	  for(j = 0; j < size_y; ++j){
	    	      
	    data_vals[(((size_z - bdepth + 1 + k)*size_x*size_y) + (j*size_x) + i)] = buffer_vals[((k*size_x*size_y) + (j*size_x) + i)];
	      	    
	    // for(j = 0; j < size_y; ++j)
	  }
	  
	  // for(i = 0; i < size_x; ++i)
	}

	// for(k = 0; (k < (bdepth - 1)); ++k)
      }
      
      // end of case statement
      break;
           
      // switch(axis) . . . case
    }
  
  // free up memory 
  delete [] buffer_vals;

}

template <class T_kernel>
class kernels_s {

 private:
  vector <int> kernels_s_sizes;
  vector<T_kernel> kernels_s_scale;
  vector< vector< vector<T_kernel> > > kernels_s_mask;

 public:
  kernels_s(){ ; }
  ~kernels_s();
  int NOk(){ return (kernels_s_sizes.size() / 2); }
  void AddKernel(int x, int y){

    kernels_s_sizes.push_back(x);
    kernels_s_sizes.push_back(y);
    kernels_s_scale.push_back(0.0);
    kernels_s_mask.push_back( vector< vector<T_kernel> () > () );

  }

  void std(){ ; }

};

template <class T_kernel>
void MakeStdKerns_S(int &NOsf, T_kernel *** &kernels_s, T_kernel * &kernels_s_scale, int ** &kernels_s_sizes){
  
  int i,j,x,y;
  T_kernel kernel_sigma_x,kernel_sigma_y,sum_w,weight;

  NOsf = 10;
  
  if(NOsf > 0){
    
    kernels_s_scale = new T_kernel[NOsf];
    kernels_s_sizes = new int * [NOsf];
    for(i = 0; i < NOsf; ++i){ kernels_s_sizes[i] = new int[2]; }
    kernels_s_sizes[0][0] = kernels_s_sizes[0][1] = 1;
    kernels_s_sizes[1][0] = kernels_s_sizes[1][1] = 1;
    kernels_s_sizes[2][0] = kernels_s_sizes[2][1] = 2;
    kernels_s_sizes[3][0] = kernels_s_sizes[3][1] = 3;
    kernels_s_sizes[4][0] = kernels_s_sizes[4][1] = 4;
    kernels_s_sizes[5][0] = kernels_s_sizes[5][1] = 6;
    kernels_s_sizes[6][0] = kernels_s_sizes[6][1] = 8;
    kernels_s_sizes[7][0] = kernels_s_sizes[7][1] = 11;
    kernels_s_sizes[8][0] = kernels_s_sizes[8][1] = 16;
    kernels_s_sizes[9][0] = kernels_s_sizes[9][1] = 23;
    kernels_s = new T_kernel ** [NOsf];
    for(i = 0; i < NOsf; ++i){
      
      kernels_s[i] = new T_kernel * [(1 + kernels_s_sizes[i][1])];
      for(j = 0; j < (1 + kernels_s_sizes[i][1]); ++j){ kernels_s[i][j] = new T_kernel [(1 + kernels_s_sizes[i][0])]; }
      
    }
    
    // calculate normalised weights
    for(i = 0; i < NOsf; ++i){
      
      kernel_sigma_x = ((T_kernel) kernels_s_sizes[i][0] + 0.5) / 3.0;
      kernel_sigma_y = ((T_kernel) kernels_s_sizes[i][1] + 0.5) / 3.0;
      sum_w = 0.0;
      for(y = 0; y < (1 + kernels_s_sizes[i][1]); ++y){
	
	for(x = 0; x < (1 + kernels_s_sizes[i][0]); ++x){
	  
	  weight = 4.0;
	  if(x == 0){ weight = (y == 0) ? 1.0 : 2.0; }

	  kernels_s[i][y][x] = 0.25 * ((erff(((0.5 + (T_kernel) x)  / (kernel_sigma_x * sqrtf(2.0))))) - (erff((((T_kernel) x - 0.5) / (kernel_sigma_x * sqrtf(2.0)))))) * ((erff(((0.5 + (T_kernel) y) / (kernel_sigma_y * sqrtf(2.0))))) - (erff((((T_kernel) y - 0.5) / (kernel_sigma_y * sqrtf(2.0))))));
	  sum_w+=(weight * kernels_s[i][y][x]);
	  
	}
	
      }      
      
      for(y = 0; y < (1 + kernels_s_sizes[i][1]); ++y){
	
	for(x = 0; x < (1 + kernels_s_sizes[i][0]); ++x){
	  
	  kernels_s[i][y][x]/=sum_w;
	  
	}
	
      }
      
    }
    
    // calculate noise scaling for this kernel
    for(i = 0; i < NOsf; ++i){
      
      kernels_s_scale[i] = 0.0;
      
      for(y = 0; y < (1 + kernels_s_sizes[i][1]); ++y){
	
	for(x = 0; x < (1 + kernels_s_sizes[i][0]); ++x){
	  
	  weight = 4.0;
	  if(x == 0){ weight = (y == 0) ? 1.0 : 2.0; }

	  kernels_s_scale[i]+=(weight * kernels_s[i][y][x]*kernels_s[i][y][x]);
	  
	}
	
      }

      kernels_s_scale[i] = sqrtf(kernels_s_scale[i]);
      
    }
    
    // if(NOsf > 0)
  }
  
}

template <class T_kernel>
void MakeStdKerns_F(int &NOff, T_kernel ** &kernels_f, T_kernel * &kernels_f_scale, int * &kernels_f_sizes){

  int i,x;
  T_kernel sum_w,weight;

  NOff = 6;
  
  kernels_f_scale = new T_kernel[NOff];
  kernels_f_sizes = new int[NOff];
  kernels_f_sizes[0] = 1;
  kernels_f_sizes[1] = 2;
  kernels_f_sizes[2] = 4;
  kernels_f_sizes[3] = 8;
  kernels_f_sizes[4] = 16;
  kernels_f_sizes[5] = 32;
  kernels_f = new T_kernel * [NOff];
  for(i = 0; i < NOff; ++i){
    
    kernels_f[i] = new T_kernel[(1 + kernels_f_sizes[i])];
    
  }
  
  // calculate normalised weights
  for(i = 0; i < NOff; ++i){
    
    sum_w = 0.0;
    for(x = 0; x < (1 + kernels_f_sizes[i]); ++x){
      
      kernels_f[i][x] = 1.0 / ((T_kernel) (1 + (2 * kernels_f_sizes[i])));
      sum_w+=(2.0 * kernels_f[i][x]);
      
    }
    sum_w-=kernels_f[i][0];

    for(x = 0; x < (1 + kernels_f_sizes[i]); ++x){
      
      kernels_f[i][x]/=sum_w;
      
    }
    
  }
    
  // calculate noise scaling for this kernel
  for(i = 0; i < NOff; ++i){
    
    kernels_f_scale[i] = 0.0;
 
    for(x = 0; x < (1 + kernels_f_sizes[i]); ++x){ 
	  
      kernels_f_scale[i]+=(2.0 * kernels_f[i][x]*kernels_f[i][x]);
      
    }

    kernels_f_scale[i]-=(kernels_f[i][0]*kernels_f[i][0]);
    kernels_f_scale[i] = sqrtf(kernels_f_scale[i]);
    
  }
    
}

template <class T_data, class T_kernel, class T_thresh, class T_noise>
  void MultiScale_IThresh(T_data * data_vals, int * flag_vals, short unsigned int * temp_flag_vals, int * sorted_vals, short unsigned int &curr_flag_val, int NOx, int NOy, int NOz, int NOx_max, int NOy_max, int NOz_max, T_noise noise_offset, T_noise noise_level, T_thresh ithresh_1, T_thresh ithresh_2, int thresh_mode, int sradius, int fradius, int flag_val_em, int NOsf, T_kernel *** kernels_s, T_kernel * kernels_s_scale, int ** kernels_s_sizes, int NOff, T_kernel ** kernels_f, T_kernel * kernels_f_scale, int * kernels_f_sizes, int vb_flag){

  void MathMorphOpen(short unsigned int * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag, int vb_flag);
  void RemHitMiss_4conn(short unsigned int * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, short unsigned int curr_flag_val, int vb_flag);
  void RemHitMiss_8conn(short unsigned int * data_vals, int size_x, int size_y, int size_z, int sradius, int fradius, short unsigned int curr_flag_val, int vb_flag);

  double sum_w;
  int x,y,z,i,sf,ff,test_flag;
  int x_pos,x_start,x_finish,y_pos,y_start,y_finish,z_pos,z_start,z_finish;

  // initialise sorted_vals array
  if(vb_flag >= 0){ std::cout << "Initialising sorted_vals array . . . " << std::endl; }
  i = 0;
  for(z = 0; z < NOz; ++z){

    for(y = 0; y < NOy; ++y){

      for(x = 0; x < NOx; ++x){

	sorted_vals[((z * NOx * NOy) + (y * NOx) + x)] = i;
	++i;

      }

    }

  }
  //cout << "Final i value = " << i << endl;

  // sort sorted_vals array according to data_vals
  if(vb_flag >= 0){ std::cout << "Sorting " << (NOx * NOy * NOz) << " values . . . " << std::endl; }
  HeapSort_index((NOx * NOy * NOz),data_vals,(sorted_vals - 1));

  // update the flag_vals array with all individual pixels/voxels above ithresh_1
  if(vb_flag >= 0){ std::cout << "Applying threshold_1 to individual pixels/voxels . . . " << std::endl; }
  for(i = (NOx * NOy * NOz) - 1; i >= 0; --i){

    if(data_vals[sorted_vals[i]] >= ((T_data) noise_offset + ((T_data) ithresh_1 * (T_data) noise_level))){ flag_vals[sorted_vals[i]] = flag_val_em; }

    // for(i = (NOx * NOy * NOz) - 1; i >= 0; i--)
  }

  // for each kernel, populate the temporary array with pixels/voxels above ithresh_2, apply
  // a size/morphology filter and add the remainder to the flag_vals array
  if(vb_flag >= 0){ std::cout << "Applying threshold_2 to smoothed pixels/voxels . . . " << std::endl; }
  // loop through spatial filters
  for(sf = 0; sf < NOsf; ++sf){

    // loop through frequency/spectral filters
    for(ff = 0; ff < NOff; ++ff){

      // populate temp_flag_vals by applying threshold to smoothed data_vals
      if(vb_flag >= 0){ std::cout << "Applying threshold_2 for kernel combination s:" << sf << " and f:" << ff << " . . . " << std::endl; }
      for(i = (NOx * NOy * NOx) - 1; i >= 0; --i){
	
	if(data_vals[sorted_vals[i]] >= ((T_data) noise_offset + ((T_data) ithresh_2 * (T_data) noise_level * (T_data) kernels_s_scale[sf] * (T_data) kernels_f_scale[ff]))){
	  
	  // calculate x,y,z position for this i
	  z_pos = (int) floorf(((float) i / (float) (NOx * NOy)));
	  y_pos = (int) floorf(((float) (i - (z_start * NOx * NOy)) / (float) NOx));
	  x_pos = i - (z_start * NOx * NOy) - (y_start * NOx);
	  
	  // calculate x,y,z range for this kernel about this i
	  z_start = z_pos - kernels_f_sizes[ff];
	  if(z_start < 0){ z_start = 0; }
	  y_start = y_pos - kernels_s_sizes[sf][1];
	  if(y_start < 0){ y_start = 0; }
	  x_start = x_pos - kernels_s_sizes[sf][0];
	  if(x_start < 0){ x_start = 0; }
	  z_finish = 1 + z_pos + kernels_f_sizes[ff];
	  y_finish = 1 + y_pos + kernels_s_sizes[sf][1];
	  x_finish = 1 + x_pos + kernels_s_sizes[sf][0];
	  
	  // calculate weighted sum of pixels/voxels about this position
	  sum_w = 0.0;
	  for(z = z_start; ((z < z_finish) && (z < NOz)); ++z){
	    
	    for(y = y_start; ((y < y_finish) && (y < NOy)); ++y){
	      
	      for(x = x_start; ((x < x_finish) && (x < NOx)); ++x){
		
		sum_w+=((double) data_vals[((z * NOx * NOy) + (y * NOx) + x)] * (double) kernels_s[sf][((int) fabs((y - y_pos)))][((int) fabs((x - x_pos)))] * (double) kernels_f[ff][((int) fabs((z - z_pos)))]);
			
	      }
		  
	    }

	  }

	  // test if the weighted sum exceeds the required noise threshold
	  if(sum_w >= ((double) noise_offset + ((double) ithresh_2 * (double) noise_level * (double) kernels_s_scale[sf] * (double) kernels_f_scale[ff]))){

	    temp_flag_vals[sorted_vals[i]] = curr_flag_val;

	  }

	  // if(data_vals[sorted_vals[i]] >= (ithresh_2 * noise_level * kernels_s_scale[sf]))
	}

	// for(i = (NOx * NOy * NOx) - 1; i >= 0; i--)
      }

      // apply size/morphology filter
      if(vb_flag >= 0){ std::cout << "Applying size/morphology filter for kernel combination s:" << sf <<" and f:" << ff << " . . . " << std::endl; }
      switch(thresh_mode){
      case 0:
	break;
      case 1:
	// removal of 4-connected hit-and-miss transform for solitary pixels/voxels
	RemHitMiss_4conn(temp_flag_vals,NOx,NOy,NOz,1,1,curr_flag_val,vb_flag);
	break;
      case 2:
	// removal of 8-connected hit-and-miss transform for solitary pixels/voxels
	RemHitMiss_8conn(temp_flag_vals,NOx,NOy,NOz,1,1,curr_flag_val,vb_flag);
	break;
      case 3:
	// opening using 4-connected cross kernel
	MathMorphOpen(temp_flag_vals,NOx,NOy,NOz,1,1,1,vb_flag);
	break;
      case 4:
	// opening using 8-connected cross kernel
	MathMorphOpen(temp_flag_vals,NOx,NOy,NOz,1,1,-1,vb_flag);
	break;
      default:
	break;
      }
      
      // update flag_vals array with remaining temp_flag_vals values
      if(vb_flag >= 0){ std::cout << "Updating flag_vals array with the results for kernel combination s:" << sf << " and f:" << ff << " . . . " << std::endl; }
      for(z = 0; z < NOz; ++z){
	
	for(y = 0; y < NOy; ++y){
	  
	  for(x = 0; x < NOx; ++x){
	    
	    if(temp_flag_vals[((z * NOx * NOy) + (y * NOx) + x)] == curr_flag_val){ flag_vals[((z * NOx * NOy) + (y * NOx) + x)] = flag_val_em; }
	    
	  }
	  
	}
	
      }
      
      // increment curr_flag_val
      ++curr_flag_val;
      if(vb_flag >= 0){ std::cout << "Incrementing curr_flag_val to " << curr_flag_val << " . . . " << std::endl; }
      
      // test range of curr_flag_val
      if(curr_flag_val >= 65535){
	
	curr_flag_val = 1;
	for(z = 0; z < NOz_max; ++z){
	  
	  for(y = 0; y < NOy_max; ++y){
	    
	    for(x = 0; x < NOx_max; ++x){
	      
	      temp_flag_vals[((z * NOy_max * NOx_max) + (y * NOx_max) + x)] = 0;
	      
	    }
	    
	  }
	  
	}
	
	// if(curr_flag_val >= 65535)
      }
   
      // for(ff = 0; ff < NOff; ff++)
    }

    // for(sf = 0; sf < NOsf; sf++)
  }

  if(vb_flag >= 0){ std::cout << "Finished finding sources . . . " << std::endl; }

  // void MultiScale_IThresh()
}

#endif

