/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/CBout.hxx
    This file is part of ConciBundle, a C/C++ library for convex optimization.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

***************************************************************************** */



#ifndef CBOUT_HXX
#define CBOUT_HXX


/**  @file CBout.hxx
    @brief Header declaring the output class CBout 
    @version 1.0
    @date 2014-07-25
    @author Christoph Helmberg
*/

#include <iostream>
#include <cstdlib>

namespace ConicBundle {

/** @defgroup BasicIOsupport Basic Support for Output Levels

*/
//@{




  /** @brief base class for uniform use of WARNINGS and ERRORS (at some point in time) 
 */

class CBout
{
  /// not output at all if out==0, otherwise use this output stream
  std::ostream *out;


  /// nonnegative level of output, 0 should mean WARNINGS and ERRORS only, 1 should represent normal output, everything above is regarded as more and more detailed information for debugging purposes
  int print_level;  

public:
  /** @brief Specifies the output level (out==NULL: no output at all, 
           out!=NULL and level=0: errors and warnings, 
           level>0 increasingly detailed information)

     @param[in] out  (std::ostream*) 
       direct all output to (*out). If out==NULL, there will be no output at all.

     @param[in] print_level (int) any value <=0 should result in WARNINGS and ERRORS only
  */
  virtual void set_out(std::ostream* out=0,int print_level=1);

  /** @brief Specifies the output level relative to the given CBout class

     @param[in] cb  
       use the same outstream as @a cb but change the print level by @a incr

     @param[in] incr 
      increment the print_level of @a cb by this value 
  */
  virtual void set_cbout(const CBout* cb,int incr=-1);

  ///
  virtual ~CBout();  //implemented in CBout.cxx
  
  /// reset to default settings (out=0,print_level=1)
  void clear_cbout(){out=0;print_level=1;}

  /// calls set_cbout
  CBout(const CBout* cb=0,int incr=-1){set_cbout(cb,incr);}

  ///initialize correspondingly 
  CBout(std::ostream* outp,int pl=1){set_out(outp,pl);}

  ///copy constructor
  CBout(const CBout& cb,int incr=0){set_out(cb.out,cb.print_level+incr);}
  
  /** @brief Returns true if out!=0 and (pl<print_level), pl<0 should be used for WARNINGS and ERRORS only, pl==0 for usual output
  */
  virtual bool cb_out(int pl=-1) const
  {if ((out)&&(pl<print_level)) return true; return false;}

  /** @brief If cb_out() returned true, this returns the output stream, but it will abort if called with out==0.
  */
  std::ostream& get_out() const
  {if (out==0) std::abort(); return *out;}
  
  /** @brief returns the pointer to the output stream
  */
  std::ostream* get_out_ptr() const
  {return out;}
  
  /** @brief returns the print_level
  */
  int get_print_level() const
  {return print_level;}

  /** @brief writes problem data to the given outstream
  */
  virtual int mfile_data(std::ostream& out) const
  { out<<"\n%\n% here something is missing, because mfile_data is not yet implemented\n%\n\n"; return 0;}

  
};


  //@}

}

#endif

