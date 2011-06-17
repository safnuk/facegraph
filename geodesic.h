// geodesic.h

/* struct which records the configuration of a geodesic from 
 * boundary to a vertex. The geodesic is perpindicular to the
 * boundary. 
 *      length = distance from boundary to vertex
 *      position = distance from first vertex of the
 *                 boundary cycle to starting point of geodesic
 *      boundary = boundary component
 *      angle = counterclockwise rotation from the geodesic to
 *              the first edge incident to the vertex
 */
struct geodesic {
    int boundary;
    double length;
    double position;
    double angle;
    geodesic(int b=0, double l=0, double p=0, double a=0) :
      boundary(b), length(l), position(p), angle(a) {}
    bool operator<(const geodesic& g) {return length < g.length;}
    bool operator<=(const geodesic& g) {return length <= g.length;}
    bool operator>(const geodesic& g) {return length > g.length;}
    bool operator>=(const geodesic& g) {return length >= g.length;}
    geodesic operator+(const geodesic& g) {
      geodesic result(boundary, length+g.length, position+g.position,
          angle+g.angle);
      return result}
    geodesic operator-(const geodesic& g) {
      geodesic result(boundary, length-g.length, position-g.position,
          angle-g.angle);
      return result}
    geodesic operator*(const geodesic& g) {
      geodesic result(boundary, length*g.length, position*g.position,
          angle*g.angle);
      return result}
    geodesic operator/(const geodesic& g) {
      geodesic result(boundary, length/g.length, position/g.position,
          angle/g.angle);
      return result}
    geodesic operator+=(const geodesic& g) {
      length+=g.length; position+=g.position; angle+=g.angle;
      return *this;
    }
    geodesic operator-=(const geodesic& g) {
      length-=g.length; position-=g.position; angle-=g.angle;
      return *this;
    }
    geodesic operator*=(const geodesic& g) {
      length*=g.length; position*=g.position; angle*=g.angle;
      return *this;
    }
    geodesic operator/=(const geodesic& g) {
      length/=g.length; position/=g.position; angle/=g.angle;
      return *this;
    }
    geodesic operator*=(const double d) {
      length*=d; position*=d; angle*=d;
      return *this;
    }
    geodesic operator/=(const double d) {
      length/=d; position/=d; angle/=d;
      return *this;
    }
    set(int b, double l, double p, double a)
      {boundary=b; length=l; position=p; angle=a;}
};

/* Struct used to encode which incident vertices
 * are closer to the boundary, and which are further.
 *
 * The further vertices are listed by their indices in the
 * array incident_vertices[], while the closer are 
 * listed by direct pointers to the vertices themselves.
 * This is because the further vertices need to know their ordering
 * around the vertex in order to calculate the geodesic, while
 * the closer ones need to be easily removed from the list as they
 * are hit.
 */
typedef struct {
  std::list<int> farther_vertices;
  std::list<vertex*> closer_vertices;
} vertex_config;

geodesic average(const std::list<geodesic>& g);
geodesic std_dev(const std::list<geodesic>& g);
geodesic sqrt(const geodesic &g);

