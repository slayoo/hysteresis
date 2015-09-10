#pragma once

#include <stdexcept>                                                            
#include <gsl/gsl_errno.h> 

void gslerr2xcpt(const char* reason, const char* file, int line, int gsl_errno) 
{                                                                                  
  throw std::runtime_error(std::string(file) + ": " + std::string(reason));     
} 

void gslerr_init()
{
  static bool set = false;
  if (!set) gsl_set_error_handler(gslerr2xcpt);
  set = true;
}
