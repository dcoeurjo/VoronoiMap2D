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


auto sitePredicate =  [](double val){ return (val==42)? true:false; };

template<typename T>
Image<Index> computeVoronoiMap(const Image<T> &source)
{
  Image<Index> voromap(source.width(), source.height());
  const Index infty =source.width()*source.height() + 1;
  
  auto closest = [&source](const Index x, const Index y,
                           const Index siteA, const Index siteB)
  {
    std::cout<<"Closest" ;
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
    int a,b, c;
    a = sBy - sAy;
    b = sCy - sBy;
    c = a + b;
    int d2_A=sAx*sAx, d2_B=sBx*sBx, d2_C=sCx*sCx;
    return (c * d2_B -  b*d2_A - a*d2_C - a*b*c) > 0 ;
  };
  
  auto print =[&](const Index pos){if (pos == infty) return std::string("(inf)");
    else return "("+std::to_string(source.getX(pos))+","+std::to_string(source.getY(pos))+") ["+std::to_string(pos)+"]";};
  
  std::cout<<"infty= "<<infty<< " " <<print(infty)<<std::endl;
#pragma omp parallel for schedule(dynamic)
  for(auto y = 0; y < source.height(); ++y)
  {
    Index prev=infty;
    bool firstSite=true;
    for (Index x=0; x < source.width(); ++x)
    {
      voromap(x,y) = source.index(prev,y);
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
  
  
//#pragma omp parallel for schedule(dynamic)
  for(auto x = 0; x < 1; ++x)// source.width(); ++x)
  {
    std::vector<Index> sites(source.height());
    std::size_t nbSites=0;
    //Looking for the sites
    for(auto y = 0; y < source.height(); ++y)
      if (voromap(x,y) != infty)
        sites[nbSites++] = voromap(x,y);
    
    std::cout<<"Sites: "<<std::endl;
    for(auto i=0; i < nbSites; ++i)
      std::cout<<print(sites[i])<<" - ";
    std::cout<<std::endl;
    
    //Prune
    std::vector<bool> remove(nbSites,false);
    Index pos = 2;
    if (nbSites>=2)
    {
      Index A=0, B=1, C=2;
      while (( C < nbSites ) &&
             ( hiddenBy(sites[A],sites[B],sites[C],x)))
      {
        //remove B
        std::cout<<"Hidden by : "<<print(sites[B])<<std::endl;
        remove[B]=true;
        B++;
        C++;
      }
    }
    
    std::cout<<"Sites: "<<std::endl;
      for(auto i=0; i < nbSites; ++i)
        std::cout<<print(sites[i])<< (remove[i]? "X":"ok") <<" - ";
      std::cout<<std::endl;
      
    
    //Rewrite
    pos = 0;
    if (nbSites==1)
      for(auto y = 0; y < source.height(); ++y)
        voromap(x,y) = sites[0];
    else
      for(auto y = 0; y < source.height(); ++y)
      {
        std::cout<<"y= "<<y<<std::endl;
        while ((pos < (nbSites-1)) &&
               remove[pos] &&
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
  Image<double> test(16,4);
  test(4,0)  = 42.0;
  test(12,0) = 42.0;
  test(8,1)  = 42.0;
  test(0,3)  = 42.0;
  test(2,3)  = 42.0;
  Image<Index> voro = computeVoronoiMap<double>(test);

  std::cout<<std::endl;
  std::cout<<std::endl;

  for(auto y=0; y < voro.height(); ++y)
  {
    for(auto x=0; x < voro.width(); ++x)
      std::cout<<voro(x,y)<<" ";
    std::cout<<std::endl;
  }
  return 0;
}
