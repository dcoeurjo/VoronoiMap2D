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
//#include "VoronoiMap2D.h"

template <typename T>
struct Image{
  Image(const std::size_t width, const std::size_t height): m_width(width), m_height(height){ m_data.resize(width*height);}
  T &operator()(const std::size_t x,
                const std::size_t y) {return m_data[x+y*m_width];}
  const T &operator()(const std::size_t x,
              const std::size_t y) const {return m_data[x+y*m_width];}
  std::size_t width(){return m_width;}
  std::size_t height(){return m_width;}
  std::size_t m_width;
  std::size_t m_height;
  std::vector<T> m_data;
};

using Point = std::array<size_t,2>;

int main()
{
  //init
  Image<double> test(128,128);
  std::cout<<"Test image: "<<test.m_width<<"x"<<test.m_height<<std::endl;
  std::cout<<"Test read: "<<test(64,64)<<std::endl;
  test(64,64) = 42.0;
  std::cout<<"Test w/r: "<<test(64,64)<<std::endl;
  
  Image<Point> voromap(128,128);
  voromap(12,12) = {1,2};
  std::cout<<"Test w/r: "<<voromap(12,12)[0]<<","<<voromap(12,12)[1]<<std::endl;
  
  
  return 0;
}
