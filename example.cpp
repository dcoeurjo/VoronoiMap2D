/*
Copyright (c) 2020 CNRS
David Coeurjolly <david.coeurjolly@liris.cnrs.fr>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIEDi
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <iostream>
#include <vector>
#include <array>
#include <assert.h>
#include <iomanip>
#include <fstream>
#include <cmath>

//#include "VoronoiMap2D.h"

using Index = std::size_t;

//Simple templated Image class
template <typename T>
struct Image{
  Image(const Index width, const Index height): m_width(width), m_height(height){ m_data.resize(width*height);}
  T &operator()(const Index x, const Index y) {return m_data[x+y*m_width];}
  const T &operator()(const Index x, const Index y) const {return m_data[x+y*m_width];}
  T &operator()(const Index idx) {return m_data[idx];}
  const T &operator()(const Index idx) const {return m_data[idx];}
  std::size_t width() const {return m_width;}
  std::size_t height() const {return m_height;}
  Index index(const Index x, const Index y) const {return x+y*m_width;}
  Index getX(const Index pos) const {return pos % m_width;}
  Index getY(const Index pos) const {return pos / m_width;}
  std::size_t m_width;
  std::size_t m_height;
  std::vector<T> m_data;
};

template <typename Image>
void exportSVG(const Image &image, const std::string &filename, const int scale=10)
{
  std::ofstream ofs (filename, std::ofstream::out);
  ofs<<"<svg width=\""<< image.width() * scale<<"\" height=\""<<image.height()*scale<<"\">"<<std::endl;
  for(auto i=0; i < image.m_data.size();++i)
    ofs<<"<rect x=\""<<image.getX(i)*scale<<"\" y=\""<< image.getY(i)*scale<<"\"  width=\""<<scale<<"\"  height=\""<<scale<<"\" style=\"fill:rgb("<< 255-(int)image(i) *12 %256<<","<<255-(int)image(i) *123 %256<<","<<255-(int)image(i) *7 %256<< ")\"  stroke=\"lightgray\"/>"<<std::endl;
  ofs<<"</svg>"<<std::endl;
}
template <typename Image>
void exportSVGDistance(const Image &image, const std::string &filename, const double max=255.0, const int scale=10)
{
  std::ofstream ofs (filename, std::ofstream::out);
  ofs<<"<svg width=\""<< image.width() * scale<<"\" height=\""<<image.height()*scale<<"\">"<<std::endl;
  for(auto i=0; i < image.m_data.size();++i)
  {
    double val  = std::sqrt((image.getX(image(i)) - image.getX(i))*(image.getX(image(i)) - image.getX(i))+
                           (image.getY(image(i)) - image.getY(i))*(image.getY(image(i)) - image.getY(i)));
    
    auto col = std::to_string((int)std::round((val*255/(double)max ))) + ",0,0";
    if (image(i)==i) col="40,40,255";
    ofs<<"<rect x=\""<<image.getX(i)*scale<<"\" y=\""<< image.getY(i)*scale<<"\"  width=\""<<scale<<"\"  height=\""<<scale<<"\" style=\"fill:rgb("<<  col << ")\" />"<<std::endl;
  }
  ofs<<"</svg>"<<std::endl;
}



template<typename T>
Image<Index> computeVoronoiMap(const Image<T> &source, const std::function<bool(const T)>& sitePredicate)
{
  Image<Index> voromap(source.width(), source.height());
  const Index infty =source.width()*source.height() + 1;
  
  auto closest = [&source](const Index x, const Index y,
                           const Index siteA, const Index siteB)
  {
    auto sAx=source.getX(siteA), sAy=source.getY(siteA);
    auto sBx=source.getX(siteB), sBy=source.getY(siteB);
    return (x-sAx)*(x-sAx) + (y-sAy)*(y-sAy) > (x-sBx)*(x-sBx) + (y-sBy)*(y-sBy);
  };
  
  auto hiddenBy = [&source](const Index A, const Index B, const Index C,
                            const Index slabX)
  {
    auto sAx=source.getX(A), sAy=source.getY(A);
    auto sBx=source.getX(B), sBy=source.getY(B);
    auto sCx=source.getX(C), sCy=source.getY(C);
    int a,b,c;
    a = sBy - sAy;
    b = sCy - sBy;
    c = a + b;
    int dA=(sAx-slabX)*(sAx-slabX), dB=(sBx-slabX)*(sBx-slabX), dC=(sCx-slabX)*(sCx-slabX);
    return (c * dB -  b*dA - a*dC - a*b*c) > 0 ;
  };
  
  auto print =[&](const Index pos){if (pos == infty) return std::string("(inf)");
    else return "("+std::to_string(source.getX(pos))+","+std::to_string(source.getY(pos))+")["+std::to_string(pos)+"]";};
  
#pragma omp parallel for schedule(dynamic)
  for(auto y = 0; y < source.height(); ++y)
  {
    Index prev=infty;
    bool firstSite=true;
    for (Index x=0; x < source.width(); ++x)
    {
      voromap(x,y) = std::min(infty,source.index(prev,y));
      if (sitePredicate(source(x,y)))
      {
        Index site =source.index(x,y);
        voromap(x,y) = site;
        if (firstSite)
        {
          firstSite=false;
          for(auto xx=0;xx < x; ++xx)
            voromap(xx,y) = site;
        }
        else //rewrite
          for(Index xx=(prev+x)/2+1 ; xx < x ; ++xx)
            voromap(xx,y) = site;
        prev = x;
      }
    }
  }
  std::cout<<"First Step"<<std::endl;
  for(auto y=0; y < voromap.height(); ++y)
  {
    for(auto x=0; x < voromap.width(); ++x)
      std::cout<<print(voromap(x,y))<<" ";
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
  
#pragma omp parallel for schedule(dynamic)
  for(auto x = 0; x < source.width(); ++x)
  {
    std::vector<Index> sites;
    //Looking for the sites
    for(auto y = 0; y < source.height(); ++y)
      if (voromap(x,y) != infty)
      {
        while (( sites.size() >= 2 ) &&
               (hiddenBy(sites[sites.size()-2], sites[sites.size()-1] , voromap(x,y) ,x)))
        {
          std::cout<<"Poping "<<print(sites[sites.size()-1])<<std::endl;
          sites.pop_back();
        }
        sites.push_back( voromap(x,y) );
      }
    
    std::cout<<"Sites: "<<std::endl;
    for(auto i=0; i < sites.size(); ++i)
      std::cout<<print(sites[i])<<" - ";
    std::cout<<std::endl;
    
    //Rewrite
    Index pos = 0;
    if (sites.size()==1)
      for(auto y = 0; y < source.height(); ++y)
        voromap(x,y) = sites[0];
    else
      for(auto y = 0; y < source.height(); ++y)
      {
        while ((pos < (sites.size()-1)) &&
               (closest( x,y, sites[pos] , sites[pos+1])))
          pos ++;
        voromap(x,y) = sites[pos];
      }
  }
  return voromap;
}


int main()
{
  //init
  Image<double> test(64,64);
  for(auto i = 0 ; i < 64; ++i)
    test( rand() % (64*64) ) = 42.0;
  exportSVG(test,"test.svg");
  
  Image<Index> voro = computeVoronoiMap<double>(test, [](const double val){ return (val==42)? true:false; } );
    
  exportSVG(voro,"result.svg");
  exportSVGDistance(voro,"result-dt.svg",8);

  return 0;
}
