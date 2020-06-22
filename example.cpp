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

//#include "VoronoiMap2D.h"

using Index = std::size_t;

//Simple templated Image class
template <typename T>
struct Image{
  Image(const Index width, const Index height): m_width(width), m_height(height){ m_data.resize(width*height);}
  T &operator()(const Index x,
                const Index y) {return m_data[x+y*m_width];}
  const T &operator()(const Index x,
              const Index y) const {return m_data[x+y*m_width];}
  T &operator()(const Index idx) {return m_data[idx];}
  const T &operator()(const Index idx) const {return m_data[idx];}
  std::size_t width() const {return m_width;}
  std::size_t height() const {return m_width;}
  Index index(const std::size_t x, const std::size_t y) {return x+y*m_width;}
  std::size_t m_width;
  std::size_t m_height;
  std::vector<T> m_data;
};


auto sitePredicate =  [](double val){ return (val==42)? true:false; };

template<typename T>
Image<Index> computeVoronoiMap(const Image<T> &source)
{
  Image<Index> voromap(source.width(), source.height());
  const Index infty =source.width()*source.height() + 1;
  
  auto closest1D=[](const int32_t x, const int32_t siteA, const int32_t siteB){ return std::abs(x-siteB) <= std::abs(x-siteA);};
  
  for(auto i=0 ; i < source.width()* source.height(); ++i)
  {
    if (sitePredicate(source(i)))
      voromap(i) = i;
    else
      voromap(i) = infty;
  }
  
  for(auto y = 0; y < source.height(); ++y)
  {
    int32_t lastsiteabscissa = infty;
    std::vector<Index> sites(source.width());
    unsigned int nbSites=0;
    
    //first scan
    for(auto x = 0; x < source.width(); ++x)
      if (voromap(x,y) != infty)
        sites[nbSites++]= x;
        
    if (nbSites != 0)
    {
      int32_t siteId=0;
      for(auto x = 0; x < source.width(); ++x)
      {
        while ( ( siteId < nbSites - 1 ) &&
               ( closest1D(x, sites[siteId], sites[siteId+1])))
          siteId++;
              
        voromap(x,y) = sites[siteId];
      }
    }
  }
  for(auto x = 0; x < source.width(); ++x)
    for(auto y = 0; y < source.height(); ++y)
    {
      
    }
  return voromap;
}


int main()
{
  //init
  Image<double> test(128,128);
  std::cout<<"Test image: "<<test.m_width<<"x"<<test.m_height<<std::endl;
  std::cout<<"Test read: "<<test(64,64)<<std::endl;
  test(64,64) = 42.0;
  std::cout<<"Test w/r: "<<test(64,64)<<std::endl;
  
  Image<Index> voro = computeVoronoiMap<double>(test);
  
  return 0;
}
