// geodesic.c
#include <math.h>
#include <list>

#include "geodesic.h"

using namespace std;

geodesic average(const std::list<geodesic>& g)
{
  geodesic result;
  std::list<geodesic>::const_iterator i = g.begin();
  result.boundary = (*i).boundary;
  for (; i!=g.end(); ++i) {
    result += *i;
  }
  if (g.size()) {
    result /= g.size();
  }
  return result;
}

geodesic std_dev(const std::list<geodesic>& g)
{
  geodesic result;
  std::list<geodesic>::const_iterator i = g.begin();
  result.boundary = (*i).boundary;
  for (; i!=g.end(); ++i) {
    result += (*i) * (*i);
  }
  result /= g.size();
  result -= average(g) * average(g);
  return sqrt(result);
}

geodesic sqrt(const geodesic &g)
{
  geodesic result(g.boundary, sqrt(g.length), sqrt(g.position),
      sqrt(g.angle));
  return result;
}

