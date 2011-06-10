// cfile_io.h


/* Data structure for a linked list of point coordinates.
 */
typedef struct _point_tmp {
        double x, y, z;
        struct _point_tmp *next;
} _point;

/* Data structure for a linked list of faces (triangles).
 */
typedef struct _face {
        int v[3];
        struct _face *next;
} face;

/* Data structure used to record data from a .mesh file.
 * point_head and face_head are heads of linked lists
 * for vertices and faces. Any filedata object should
 * always be initialized before use (with initialize_filedata), 
 * and deallocated after finished (with deallocate_filedata).  
 */
typedef struct {
        int number_of_points;
        int number_of_triangles;
        _point *point_head;
        face *face_head;
} filedata; 

_point *add_point_node(_point *node, double x, double y, double z);
face *add_face_node(face *node, int v1, int v2, int v3);
_point *parse_new_point(_point *point_node, char *line, int point_count);
face *parse_new_face(face *face_node, char *line, int face_count);
int read(char *filename, filedata *data);
void save_mesh(char *filename, void *data);
void initialize_filedata(filedata *data);
void deallocate_filedata(filedata *data);
void deallocate_point(_point *head);
void deallocate_face(face *head);
void print_filedata(filedata *data);
