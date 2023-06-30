/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Tools/BoxPlot.cxx
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


#include "BoxPlot.hxx"

using namespace CH_Matrix_Classes;

namespace CH_Tools {

  BoxPlot::BoxPlot()
  {
    median=0.;
    lquartile=0.;
    uquartile=0.;
    lwhisker=0.;
    uwhisker=0.;
    outliers.init(0,0,0.);
    samplesz=0;
    minval=max_Real;
    maxval=min_Real;

    use_samplesz=true;
    use_logscale=false;
    use_minmax=true;
  }

  BoxPlot::~BoxPlot()
  {}

  
  int BoxPlot::standalone_preamble(std::ostream& out)
  {
    out<<"\\documentclass{standalone}\n";
    out<<"\\usepackage[utf8]{inputenc}\n";
    out<<"\\usepackage{xcolor}\n";
    out<<"\\usepackage{graphicx}\n";
    out<<"\\usepackage{tikz}\n";
    //out<<"\\usetikzlibrary{decorations,snakes}\n";
    //out<<"\\usetikzlibrary{decorations.pathmorphing}\n";
    //out<<"\\usetikzlibrary{arrows,graphs,external,positioning}\n";
    //out<<"\\usetikzlibrary{calc, decorations.markings, intersections}\n";
    //out<<"\\usetikzlibrary{patterns}\n";
    out<<"\\usepackage{pgfplots} % diagrams\n";
    out<<"\\usepgfplotslibrary{statistics} \n";
    out<<"\\pgfplotsset{compat=newest} % diagrams\n";
    out<<"\\usepackage{filecontents}\n";
    out<<"\\begin{document}\n\n";
    return 0;
  }
  
  int BoxPlot::standalone_postamble(std::ostream& out)
  {
    out<<"\n\n\\end{document}\n";
    return 0;
  }
  
  int BoxPlot::start_plot(std::ostream& out,const char* options)
  {
    out<<"\n\\begin{tikzpicture}\n";
    if (use_logscale)
      out<<"\\begin{semilogyaxis}[\n";
    else 
      out<<"\\begin{axis}[\n";
    out<<"   boxplot/draw direction=y,\n";
    if (use_samplesz)
      out<<"   boxplot/variable width,\n";
    if (options)
      out<<options<<",";
    out<<"]\n";

    return 0;
  }
  
  int BoxPlot::end_plot(std::ostream& out)
  {
    if (use_logscale)
      out<<"\\end{semilogyaxis}\n";
    else 
      out<<"\\end{axis}\n";
    out<<"\\end{tikzpicture}\n\n";

    return 0;
  }

  int BoxPlot::next_data(const Matrix& datavector)
  {
    Indexmatrix sind;
    sortindex(datavector,sind);
    samplesz=datavector.dim();

    //std::cout<<" datavector="<<datavector;
    //std::cout<<" sind="<<transpose(sind);


    if (samplesz==0){
      median=0.;
      lquartile=0.;
      uquartile=0.;
      lwhisker=0.;
      uwhisker=0.;
      minval=max_Real;
      maxval=min_Real;
      outliers.init(0,0,0.);
      return 0;
    }
    if (samplesz==1){
      minval=maxval=median=lquartile=uquartile=lwhisker=uwhisker=datavector(0);
      outliers.init(0,0,0.);
      return 0;
    }
    if (samplesz==2){
      minval=median=lquartile=lwhisker=datavector(sind(0));
      maxval=uquartile=uwhisker=datavector(sind(1));
      outliers.init(0,0,0.);
      return 0;
    }
    if (samplesz==3){
      median=datavector(sind(1));
      minval=lquartile=lwhisker=datavector(sind(0));
      maxval=uquartile=uwhisker=datavector(sind(2));
      outliers.init(0,0,0.);
      return 0;
    }


    minval=datavector(sind(0));
    maxval=datavector(sind(sind.rowdim()-1));
    
    if (samplesz%2==0)
      median=datavector(sind(samplesz/2-1));
    else
      median=.5*(datavector(sind(samplesz/2-1))+datavector(sind(samplesz/2)));
    
    if (samplesz%4==0){
      lquartile=datavector(sind(samplesz/4-1));
      uquartile=datavector(sind(3*samplesz/4-1));
    }
    else {
      lquartile=.5*(datavector(sind(samplesz/4-1))+datavector(sind(samplesz/4)));
      uquartile=.5*(datavector(sind(3*samplesz/4-1))+datavector(sind(3*samplesz/4)));
    }
    Real IQR=uquartile-lquartile;
    Real lwbound=lquartile-1.5*IQR;
    Real uwbound=uquartile+1.5*IQR;
    outliers.init(0,1,0.);
    Integer i=0;
    while((i<sind.rowdim())&&(datavector(sind(i))<lwbound)){
      outliers.concat_below(datavector(sind(i)));
      i++;
    }
    assert(i<sind.rowdim());
    lwhisker=datavector(sind(i));
    i=sind.rowdim();
    while((--i>=0)&&(datavector(sind(i))>uwbound)){
      outliers.concat_below(datavector(sind(i)));
    }
    assert(i>=0);
    uwhisker=datavector(sind(i));
    
    return 0;
  }
  
  int BoxPlot::add_plot(std::ostream& out,
			Real pos,
			const char* color,
			const char* name,
			Real val)
  {
    if (samplesz==0)
      return 0;
    out<<"\\addplot+ [ ";
    if (color)
      out<<color<<",";
    out<<" solid, mark=star, boxplot prepared={\n";
    out<<"draw position="<<pos<<",\n";
    out<<"lower whisker="<<lwhisker<<", lower quartile="<<lquartile<<",\n";
    out<<"median="<<median<<",\n";
    out<<"upper quartile="<<uquartile<<", upper whisker="<<uwhisker;
    if (use_samplesz)
      out<<",\n sample size="<<samplesz;
    out<<"\n}]";
    out<<" coordinates {";
    if (use_minmax){
      out<<" ("<<0<<","<<minval<<")";
      out<<" ("<<0<<","<<maxval<<")";
    }
    else {
      for (Integer i=0;i<outliers.rowdim();i++)
	out<<" ("<<0<<","<<outliers(i)<<")";
    }
    out<<"}";
    if (name) {
      out<<" node[below] at (boxplot box cs: ";
      out<< val<<", 0.5) {\\rotatebox{270}{\\tiny{"<<name<<"}}}";
    }
    out<<";\n";
    return 0;
  }
  
  int BoxPlot::add_tabularline(std::ostream& out,int prec)
  {
    out.precision(prec);
    if (samplesz==0){
      out<<" -- & -- & -- & -- & -- ";
    }
    else if (samplesz<=3){
      out.width(10); out<< lwhisker<<" & -- & ";
      out.width(10); out<< median<<" & -- & ";
      out.width(10); out<< uwhisker<<" ";
    }
    else {
      out.width(10); out<<minval<<" & ";
      out.width(10); out<<lquartile<<" & ";
      out.width(10); out<<median<<" & ";
      out.width(10); out<<uquartile<<" & ";
      out.width(10); out<<maxval<<" ";
    }
    return 0;
  }

}
