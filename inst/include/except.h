// except.h  Part of the CppADutils package for the R programming language.

// This file is part of CppADutils, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2015 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, you can obtain one at http://mozilla.org/MPL/2.0/.
// See the trustOptim LICENSE file for more information.


#ifndef __MB_EXCEPT
#define __MB_EXCEPT

#include <iostream>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <cstring>

#define ERROR_HANDLER R_Interface_Error_Handler
#include <R_ext/Utils.h>

using Rcpp::Rcout;

class MyException : public std::exception {

public:

  const std::string reason;
  const std::string file;
  const int line;
  
  std::string message;
  
 MyException(const std::string reason_,
	     const std::string file_,
	     const int line_) :
  reason(reason_), file(file_), line(line_)  {    
    
    std::ostringstream oss;
    oss << "\nException thrown from File " << file;
    oss << "  at Line " << line <<".\n";
    oss << "Reason : " << reason << "\n";	
    message = oss.str();
  }
  
  virtual ~MyException() throw() {};
  
  virtual const char* what() const throw() {
    return message.c_str(); 
  }
  
  void print_message() {
    Rcout << message << std::endl;
  }
};

template<typename T>
void R_Interface_Error_Handler(const T & ex) {
  // takes exception object and does R-friendly things to it
  ex.print_message();
  Rf_error("R error\n");
}

static inline void check_interrupt_impl(void* /*dummy*/) {
  R_CheckUserInterrupt();
}

inline bool check_interrupt() {
  return (R_ToplevelExec(check_interrupt_impl, NULL) == FALSE);
}

#endif
