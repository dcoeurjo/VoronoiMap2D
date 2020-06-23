//////////////////////////////////////////////////////////////////////////////
#include <cstdlib>


template < typename Predicate, typename Image, typename Point>
Image<Point> computeVoronoi()




inline
void
:compute( )
{
  //We copy the image extent
  myLowerBoundCopy = myDomainPtr->lowerBound();
  myUpperBoundCopy = myDomainPtr->upperBound();

  //Point outside the domain
  for ( auto & coord : myInfinity )
    coord = DGtal::NumberTraits< typename Point::Coordinate >::max();

  //Init
  for ( auto const & pt : *myDomainPtr )
    if ( (*myPointPredicatePtr)( pt ))
      myImagePtr->setValue ( pt, myInfinity );
    else
      myImagePtr->setValue ( pt, pt );

  //We process the remaining dimensions
  for ( Dimension dim = 0;  dim< S::dimension ; dim++ )
    computeOtherSteps ( dim );
}

template <typename S, typename P,typename TSep, typename TImage>
inline
void
DGtal::VoronoiMap<S,P, TSep, TImage>::computeOtherSteps ( const Dimension dim ) const
{
#ifdef VERBOSE
  std::string title = "VoronoiMap dimension " +  boost::lexical_cast<std::string>( dim ) ;
  trace.beginBlock ( title );
#endif

  //We setup the subdomain iterator
  //the iterator will scan dimension using the order:
  // {n-1, n-2, ... 1} (we skip the '0' dimension).
  std::vector<Dimension> subdomain;
  subdomain.reserve(S::dimension - 1);
  for ( int k = 0; k < (int)S::dimension ; k++)
    if ( static_cast<Dimension>(((int)S::dimension - 1 - k)) != dim)
      subdomain.push_back( (int)S::dimension - 1 - k );

  Domain localDomain(myLowerBoundCopy, myUpperBoundCopy);

#ifdef WITH_OPENMP
  //Parallel loop
  std::vector<Point> subRangePoints;
  //Starting point precomputation
  for ( auto const & pt : localDomain.subRange( subdomain ) )
    subRangePoints.push_back( pt );

  //We run the 1D problems in //
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < subRangePoints.size(); ++i)
    computeOtherStep1D ( subRangePoints[i], dim);

#else
  //We solve the 1D problems sequentially
  for ( auto const & pt : localDomain.subRange( subdomain ) )
    computeOtherStep1D ( pt, dim);
#endif

#ifdef VERBOSE
  trace.endBlock();
#endif
}

// //////////////////////////////////////////////////////////////////////:
// ////////////////////////// Other Phases
template <typename S,typename P, typename TSep, typename TImage>
void
DGtal::VoronoiMap<S,P,TSep, TImage>::computeOtherStep1D ( const Point &startingPoint,
                                                  const Dimension dim) const
{
  ASSERT(dim < S::dimension);

  // Default starting and ending point for a cycle
  Point startPoint = startingPoint;
  Point endPoint   = startingPoint;
  startPoint[dim]  = myLowerBoundCopy[dim];
  endPoint[dim]    = myUpperBoundCopy[dim];

  // Extent along current dimension.
  const auto extent = myUpperBoundCopy[dim] - myLowerBoundCopy[dim] + 1;

  // Site storage.
  std::vector<Point> Sites;

  // Reserve sites storage.
  // +1 along periodic dimension in order to store two times the site that is on break index.
  Sites.reserve( extent + ( isPeriodic(dim) ? 1 : 0 ) );

  // Pruning the list of sites and defining cycle bounds.
  // In the periodic case, the cycle bounds depend on the so-called break index
  // that defines the start point.
  if ( dim == 0 )
    {
      // For dim = 0, no sites are hidden.
      for ( auto point = startPoint ; point[dim] <= myUpperBoundCopy[dim] ; ++point[dim] )
        {
          const Point psite = myImagePtr->operator()( point );
          if ( psite != myInfinity )
            Sites.push_back( psite );
        }

      // If no sites are found, then there is nothing to do.
      if ( Sites.size() == 0 )
        return;

      // In the periodic case and along the first dimension, the break index
      // is at the first site found.
      if ( isPeriodic(dim) )
        {
          startPoint[dim] = Sites[0][dim];
          endPoint[dim]   = startPoint[dim] + extent - 1;

          // The first site is also the last site (with appropriate shift).
          Sites.push_back( Sites[0] + Point::base(dim, extent) );
        }
    }
  else
    {
      // In the periodic case, the cycle depends on break index
      if ( isPeriodic(dim) )
        {
          // Along other than the first dimension, the break index is at the lowest site found.
          auto minRawDist = DGtal::NumberTraits< typename SeparableMetric::RawValue >::max();

          for ( auto point = startPoint; point[dim] <= myUpperBoundCopy[dim]; ++point[dim] )
            {
              const Point psite = myImagePtr->operator()( point );

              if ( psite != myInfinity )
                {
                  const auto rawDist = myMetricPtr->rawDistance( point, psite );
                  if ( rawDist < minRawDist )
                    {
                      minRawDist = rawDist;
                      startPoint[dim] = point[dim];
                    }
                }
            }

          // If no sites are found, then there is nothing to do.
          if ( minRawDist == DGtal::NumberTraits< typename SeparableMetric::RawValue >::max() )
            return;

          endPoint[dim] = startPoint[dim] + extent - 1;
        }

      // Pruning the list of sites for both periodic and non-periodic cases.
      for( auto point = startPoint ; point[dim] <= myUpperBoundCopy[dim] ; ++point[dim] )
        {
          const Point psite = myImagePtr->operator()(point);

          if ( psite != myInfinity )
            {
              
              while (( Sites.size() >= 2 ) &&
                     ( myMetricPtr->hiddenBy(Sites[Sites.size()-2], Sites[Sites.size()-1] ,
                                             psite, startingPoint, endPoint, dim) ))
                Sites.pop_back();

              Sites.push_back( psite );
            }
        }

      // Pruning the remaining list of sites in the periodic case.
      if ( isPeriodic(dim) )
        {
          auto point = startPoint;
          point[dim] = myLowerBoundCopy[dim];
          for ( ; point[dim] <= endPoint[dim] - extent + 1; ++point[dim] ) // +1 in order to add the break-index site at the cycle's end.
            {
              Point psite = myImagePtr->operator()(point);

              if ( psite != myInfinity )
                {
                  // Site coordinates must be between startPoint and endPoint.
                  psite[dim] += extent;

                  while (( Sites.size() >= 2 ) &&
                         ( myMetricPtr->hiddenBy(Sites[Sites.size()-2], Sites[Sites.size()-1] ,
                                                 psite, startingPoint, endPoint, dim) ))
                    Sites.pop_back();

                  Sites.push_back( psite );
                }
            }
        }
    }

  // No sites found
  if ( Sites.size() == 0 )
    return;

  // Rewriting for both periodic and non-periodic cases.
  std::size_t siteId = 0;
  auto point = startPoint;

  for ( ; point[dim] <= myUpperBoundCopy[dim] ; ++point[dim] )
    {
      while ( ( siteId < Sites.size()-1 ) &&
             ( myMetricPtr->closest(point, Sites[siteId], Sites[siteId+1])
              != DGtal::ClosestFIRST ))
        siteId++;

      myImagePtr->setValue(point, Sites[siteId]);
    }

  // Continuing rewriting in the periodic case.
  if ( isPeriodic(dim) )
    {
      for ( ; point[dim] <= endPoint[dim] ; ++point[dim] )
        {
          while ( ( siteId < Sites.size()-1 ) &&
                 ( myMetricPtr->closest(point, Sites[siteId], Sites[siteId+1])
                  != DGtal::ClosestFIRST ))
            siteId++;

          myImagePtr->setValue(point - Point::base(dim, extent), Sites[siteId] - Point::base(dim, extent) );
        }
    }

}


/**
 * Constructor.
 */
template <typename S,typename P,typename TSep, typename TImage>
inline
DGtal::VoronoiMap<S,P, TSep, TImage>::VoronoiMap( ConstAlias<Domain> aDomain,
                                          ConstAlias<PointPredicate> aPredicate,
                                          ConstAlias<SeparableMetric> aMetric )
     : myDomainPtr(&aDomain)
     , myPointPredicatePtr(&aPredicate)
     , myDomainExtent( aDomain->upperBound() - aDomain->lowerBound() + Point::diagonal(1) )
     , myMetricPtr(&aMetric)
{
  myPeriodicitySpec.fill( false );
  myImagePtr = CountedPtr<OutputImage>( new OutputImage(aDomain) );
  compute();
}

template <typename S,typename P,typename TSep, typename TImage>
inline
DGtal::VoronoiMap<S,P, TSep, TImage>::VoronoiMap( ConstAlias<Domain> aDomain,
                                          ConstAlias<PointPredicate> aPredicate,
                                          ConstAlias<SeparableMetric> aMetric,
                                          PeriodicitySpec const & aPeriodicitySpec )
     : myDomainPtr(&aDomain)
     , myPointPredicatePtr(&aPredicate)
     , myDomainExtent( aDomain->upperBound() - aDomain->lowerBound() + Point::diagonal(1) )
     , myMetricPtr(&aMetric)
     , myPeriodicitySpec(aPeriodicitySpec)
{
  // Finding periodic dimension index.
  for ( Dimension i = 0; i < Space::dimension; ++i )
    if ( isPeriodic(i) )
      myPeriodicityIndex.push_back( i );

  myImagePtr = CountedPtr<OutputImage>( new OutputImage(aDomain) );
  compute();
}

template <typename S,typename P,typename TSep, typename TImage>
inline
typename DGtal::VoronoiMap<S, P, TSep, TImage>::Point
DGtal::VoronoiMap<S, P, TSep, TImage>::projectPoint( Point aPoint  ) const
{
  for ( auto const & dim : myPeriodicityIndex )
    aPoint[ dim ] = projectCoordinate( aPoint[ dim ], dim );

  return aPoint;
}

template <typename S,typename P,typename TSep, typename TImage>
inline
typename DGtal::VoronoiMap<S, P, TSep, TImage>::Point::Coordinate
DGtal::VoronoiMap<S, P, TSep, TImage>::projectCoordinate( typename Point::Coordinate aCoordinate, const Dimension aDim ) const
{
  ASSERT( aCoordinate - myDomainPtr->lowerBound()[aDim] + myDomainExtent[aDim] >= 0 );
  return ( aCoordinate - myDomainPtr->lowerBound()[aDim] + myDomainExtent[aDim] ) % myDomainExtent[aDim] + myDomainPtr->lowerBound()[aDim];
}

template <typename S,typename P,typename TSep, typename TImage>
inline
void
DGtal::VoronoiMap<S,P, TSep, TImage>::selfDisplay ( std::ostream & out ) const
{
  out << "[VoronoiMap] separable metric=" << *myMetricPtr ;
}


// //                                                                           //
// ///////////////////////////////////////////////////////////////////////////////

template <typename S,typename P,typename TSep, typename TImage>
inline
std::ostream&
DGtal::operator<< ( std::ostream & out,
                    const VoronoiMap<S,P,TSep, TImage> & object )
{
  object.selfDisplay( out );
  return out;
}
