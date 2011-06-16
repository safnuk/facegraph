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
typedef struct {
    int boundary;
    double length;
    double position;
    double angle;
} geodesic;

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

void set_geodesic(geodesic &g, int boundary, double length, double position, double angle);
