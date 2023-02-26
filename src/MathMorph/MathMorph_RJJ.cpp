#include<MathMorph_RJJ.h>

using namespace std;

// library compile instructions
//
// c++ MathMorph_RJJ.cpp -c -I/Users/Jurek83/MathMorph -O2
// Libtool -static -o libMathMorph_RJJ.a MathMorph_RJJ.o

void CreateMetric(int data_metric[3], int xyz_order[3], int NOx, int NOy, int NOz){

  // ensure that the x_order, y_order and z_order values are 1,2,3
  if((xyz_order[0] <= xyz_order[1]) && (xyz_order[0] <= xyz_order[2])){
    
    xyz_order[0] = 1;
    if(xyz_order[1] <= xyz_order[2]){ 

      xyz_order[1] = 2; 
      xyz_order[2] = 3;

    } else {

      xyz_order[1] = 3;
      xyz_order[2] = 2;

    }

  } else if((xyz_order[1] <= xyz_order[0]) && (xyz_order[1] <= xyz_order[2])){

    xyz_order[1] = 1;
    if(xyz_order[0] <= xyz_order[2]){ 

      xyz_order[0] = 2; 
      xyz_order[2] = 3;

    } else {

      xyz_order[0] = 3;
      xyz_order[2] = 2;

    }

  } else if((xyz_order[2] <= xyz_order[0]) && (xyz_order[2] <= xyz_order[1])){

    xyz_order[2] = 1;
    if(xyz_order[0] <= xyz_order[1]){ 

      xyz_order[0] = 2; 
      xyz_order[1] = 3;

    } else {

      xyz_order[0] = 3;
      xyz_order[1] = 2;

    }

  }

  // assign values to the data_metric array
  
  // 1. x value
  switch(xyz_order[0]){
  case 1:
    data_metric[0] = 1;
    break;
  case 2:
    if(xyz_order[1] == 1){ data_metric[0] = NOy; } else { data_metric[0] = NOz; }
    break;
  case 3:
    data_metric[0] = NOy * NOz;
    break;
  default:
    data_metric[0] = 1;
    break;
  }
    
  // 2. y value
  switch(xyz_order[1]){
  case 1:
    data_metric[1] = 1;
    break;
  case 2:
    if(xyz_order[0] == 1){ data_metric[1] = NOx; } else { data_metric[1] = NOz; }
    break;
  case 3:
    data_metric[1] = NOx * NOz;
    break;
  default:
    data_metric[1] = 1;
    break;
  }

  // 3. z value
  switch(xyz_order[2]){
  case 1:
    data_metric[2] = 1;
    break;
  case 2:
    if(xyz_order[0] == 1){ data_metric[2] = NOx; } else { data_metric[2] = NOy; }
    break;
  case 3:
    data_metric[2] = NOx * NOy;
    break;
  default:
    data_metric[2] = 1;
    break;
  }

}

void Populate_temp_XtoI_array(size_t * temp_vals, int i, int j, int k, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag){

  int used = 0,x,y,z,i2,j2,k2;
  int x_start, y_start, z_start, xy_stride, x_stride;
  float sep, sep_max;

  xy_stride = (1 + 2 * sradius) * (1 + 2 * sradius);
  x_stride = 1 + 2 * sradius;

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
	  
	  temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = used;
	  ++used;

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

	    temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = used;
	    used++;

	  } else {

	    temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = -1;

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

	    temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = used;
	    used++; 

	  } else {

	    temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = -1;
	    
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

	    temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = used;
	    used++; 
	    
	  } else {

	    temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = -1;
	    
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
	  
	  temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = used;
	  used++; 
	  
	} else {

	  temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = -1;

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

	    temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = used;
	    used++; 
	    
	  } else {

	    temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = -1;
	    
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

	    temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = used;
	    used++; 
	    
	  } else {

	    temp_vals[((z - z_start) * xy_stride) + ((y - y_start) * x_stride) + x - x_start] = -1;

	  }

	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	  
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }

      // for(k = (-1 * fradius); k <= fradius; ++k)
    }

  }

}

void Populate_temp_ItoX_array(size_t * temp_vals, int i, int j, int k, int size_x, int size_y, int size_z, int sradius, int fradius, int round_flag){

  int used = 0,x,y,z,i2,j2,k2;
  int x_start, y_start, z_start, xy_stride, x_stride;
  float sep, sep_max;

  xy_stride = (1 + 2 * sradius) * (1 + 2 * sradius);
  x_stride = 1 + 2 * sradius;

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
	  
	  temp_vals[used] = used;
	  ++used;

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

	    temp_vals[used] = used;
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

	    temp_vals[used] = used;
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

	    temp_vals[used] = used;
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
	  
	  temp_vals[used] = used;
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

	    temp_vals[used] = used;
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

	    temp_vals[used] = used;
	    used++; 
	    
	  } 

	  // for(i = (-1 * sradius); i <= sradius; ++i)
	}
	  
	// for(j = (-1 * sradius); j <= sradius; ++j)
      }

      // for(k = (-1 * fradius); k <= fradius; ++k)
    }

  }

}
