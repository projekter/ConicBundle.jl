/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Tools/BoxPlot.hxx
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



#ifndef CH_TOOLS__BOXPLOT_HXX
#define CH_TOOLS__BOXPLOT_HXX

/**  @file BoxPlot.hxx
    @brief Header declaring the classes CH_Tools::BoxPlot for using pgfplots
    @version 1.0
    @date 2021-07-09
    @author Christoph Helmberg

*/

#include <iostream>
#include <iomanip>
#include "matrix.hxx"

namespace CH_Tools {

/**@defgroup Plot Plot (output of statistical data with pgfplots)
*/


  //@{

  /** @brief for generating LaTeX figures with successive box plots and corresponding tables

      Each call to next_data() with a vector of dimension @a samplesz computes the
      values @a minval,@a lwhisker, @a lquartile, @a median, @a uquartile, @a uwhisker, @a maxval and collects the vector @a outliers. 

      For creating the full LaTeX source for a pdf-file of a single plot that displays the statistics of several vectors use, in this sequence
      - standalone_preamble()   ... documentclass up to begin{document}
      - start_plot()  
      - next_data() ...  generate the data for the first vector
      - add_plot()
      - ...
      - next_data()  ... generate the data for the last vector 
      - add_plot()
      - end_plot()
      - standalone_postamble()  ... end{document}

      There is limited support for adding some options. Once the data has been set in next_data, it can be used for several add_plot() calls to different out streams and some of the data can be output to a table by add_tabularline().
   */

class BoxPlot
{
private:
  CH_Matrix_Classes::Real median;      ///< as set in the last call to next_data()
  CH_Matrix_Classes::Real lquartile;   ///< as set in the last call to next_data()
  CH_Matrix_Classes::Real uquartile;   ///< as set in the last call to next_data()
  CH_Matrix_Classes::Real lwhisker;    ///< as set in the last call to next_data()
  CH_Matrix_Classes::Real uwhisker;     ///< as set in the last call to next_data()
  CH_Matrix_Classes::Real minval;    ///< as set in the last call to next_data()
  CH_Matrix_Classes::Real maxval;     ///< as set in the last call to next_data()
  CH_Matrix_Classes::Matrix outliers;  ///< as set in the last call to next_data()
  CH_Matrix_Classes::Integer samplesz; ///< as set in the last call to next_data()

  bool use_samplesz; ///< if true [=default], the box plot width depends on the relative sample size
  bool use_logscale; ///< if true [default:false], the values of the box plot are output in log scale
  
  bool use_minmax; ///< if true [=default], min and max are plotted as outliers irrespective of whether they are outliers or not, but no outliers are plotted
  
  public:
  /** @name Constructors
   */
  //@{
  BoxPlot();

  ~BoxPlot();
    //@}
    
  /** @name methods
   */
  //@{

  /// if set to true (==default), the samples size is used for the width of the box polots
  void set_use_samplesz(bool us){use_samplesz=us;}
  /// if set to true (default is false), the values are plotted with log scale
  void set_use_logscale(bool ls){use_logscale=ls;}
  /// if set to true (==default), min and max are plotted as outliers irrespective of whether they are outliers or not, but no outliers are plotted
  void set_use_minmax(bool mm){use_minmax=mm;}

  /// to generate LaTeX source that results in a seperate pdf-file of the plot use this to output the corresponding LaTeX preamble 
  int standalone_preamble(std::ostream& out);
  /// if things were started with standalone_preamble() finish the document with calling this
  int standalone_postamble(std::ostream& out);

  /// output the LaTeX commands for starting a new plot
  int start_plot(std::ostream& out, ///< outstream (the file)
		 const char* options=0 ///< as 
		 );
  /// when, after start_plot() all plots have been added by add_plot(), this ouputs the closing part of the plot
  int end_plot(std::ostream& out);

  /// input the next data group to be used in add_plot or add_tabularline; from this data the routine computes min, lwhisker, lquartile, median, uquartile, uwhisker, max, samplesz
  int next_data(const CH_Matrix_Classes::Matrix& datavector);
  
  /// output the box plot to the data of the last call to next_data(); multiple calls without intermediate calls to next_data() all output the same data 
  int add_plot(std::ostream& out, ///< outstream (the file)
	       CH_Matrix_Classes::Real pos, ///< the "draw position" of pgfplots
	       const char* color=0, ///< the color to be used for this box
	       const char* name=0, ///< the name of the data to this box, printed in rotated form
	       CH_Matrix_Classes::Real val=0. /// the name will be printed downwards starting from this value
	       );

  /// output minval & lquartile & median & uquartile & maxval as computed by the latest previous next_data() call. The number of significant digits is given in precision.
  int add_tabularline(std::ostream& out,int prec=6);

  //@}

};

//@}

}

#endif

