/* cfile_io.c
 *
*/

#include <stdlib.h>
#include <stdio.h>

/* Data structure for a linked list of point coordinates.
 */
typedef struct _point {
        float x, y, z;
        struct _point *next;
} point;

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
        point *point_head;
        face *face_head;
} filedata; 

point *add_point_node(point *node, float x, float y, float z);
face *add_face_node(face *node, int v1, int v2, int v3);
point *parse_new_point(point *point_node, char *line, int point_count);
face *parse_new_face(face *face_node, char *line, int face_count);
int read(char *filename, filedata *data);
void initialize_filedata(filedata *data);
void deallocate_filedata(filedata *data);
void deallocate_point(point *head);
void deallocate_face(face *head);
void print_filedata(filedata *data);

int main(int argc, char *argv[]) {
        filedata data;
        if ( argc != 2 ) {
                printf( "usage: %s filename\n", argv[0] );
                return 0;
        }
        initialize_filedata(&data);
        read(argv[1], &data);
        print_filedata(&data);
        deallocate_filedata(&data);
}

/* Read face data from filename, and store it in data.
 * Assumes that data has been properly initialized.
 */
int read(char *filename, filedata* data) 
{
        char line[100];
        FILE *fp;
        int point_count = 0;
        int face_count = 0;
        point *point_node = data->point_head;
        face *face_node = data->face_head;

        fp = fopen(filename, "r");
        if (fp == NULL) {
                printf("Can't open file!\n");
                exit(1);
        }
        while (fgets(line, 100, fp)) {
                if (line[0] == 'V') {
                        point_node = parse_new_point(point_node, 
                                        line, point_count);
                        point_count++;
                } else if (line[0] == 'F') {
                        face_node = parse_new_face(face_node, 
                                        line, face_count);
                        face_count++;
                }
        }
        fclose(fp);
        data->number_of_points = point_count;
        data->number_of_triangles = face_count;
        return 0;
}

/* Allocate a new point and add to tail of linked list.
 */
point *add_point_node(point *node, float x, float y, float z)
{
        point *new_node;
        new_node = (point*)malloc(sizeof(point));
        if (new_node == NULL) {
                printf("Could not allocate memory for new point.\n");
                exit(1);
        }
        node->next = new_node;
        new_node->next = NULL;
        new_node->x = x;
        new_node->y = y;
        new_node->z = z;
        return new_node;
}

/* Allocate a new face and add to tail of linked list.
 */
face *add_face_node(face *node, int v1, int v2, int v3)
{
        face *new_node;
        new_node = (face*)malloc(sizeof(face));
        if (new_node == NULL) {
                printf("Could not allocate memory for new face.\n");
                exit(1);
        }
        node->next = new_node;
        new_node->next = NULL;
        new_node->v[0] = v1;
        new_node->v[1] = v2;
        new_node->v[2] = v3;
        return new_node;
}

/* Adds a new node to linked list with tail point_node.
 * line should be a string of the form
 *      "Vertex index x y z [radius]"
 * The optional radius value is from a circle packing metric.
 */
point *parse_new_point(point* point_node, char *line, int point_count) 
{
        float x, y, z;
        float radius;
        int index;

        sscanf(line, "Vertex %i %f %f %f %f",
                        &index, &x, &y, &z, &radius);
        if (point_count != index) {
                printf("Index mismatch when reading vertex %i", index);
                exit(1);
        }
        return add_point_node(point_node, x, y, z);
}

/* Adds a new node to linked list with tail face_node.
 * line should be a string of the form:
 * "Face index v1 v2 v3"
 * where v1, v2, v3 is a triple of indices to points.
 */
face *parse_new_face(face *face_node, char *line, int face_count)
{
        int v1, v2, v3, index;

        sscanf(line, "Face %i %i %i %i", 
                        &index, &v1, &v2, &v3);
        if (face_count != index) {
                printf("Index mismatch when reading face %i", index);
                exit(1);
        }
        return add_face_node(face_node, v1, v2, v3);
}

void initialize_filedata(filedata *data)
{
        data->point_head = (point*)malloc(sizeof(point));
        data->face_head = (face*)malloc(sizeof(face));
        data->point_head->next = NULL;
        data->face_head->next = NULL;
}

void deallocate_filedata(filedata *data)
{
        deallocate_point(data->point_head);
        deallocate_face(data->face_head);
        data->point_head = NULL;
        data->face_head = NULL;
}

void deallocate_point(point *node)
{
        if (node->next != NULL) {
                deallocate_point(node->next);
        } 
        free(node);
}

void deallocate_face(face *node)
{
        if (node->next != NULL) {
                deallocate_face(node->next);
        } 
        free(node);
}

void print_filedata(filedata *data) 
{
        point *point_node = data->point_head->next;
        face *face_node = data->face_head->next;       
        while (point_node != NULL) {
                printf("V(%f, %f, %f)\n", point_node->x,
                                point_node->y, point_node->z);
                point_node = point_node->next;
        }
        while (face_node != NULL) {
                printf("F(%i, %i, %i)\n", face_node->v[0],
                                face_node->v[1], face_node->v[2]);
                face_node = face_node->next;
        }
        printf("%i vertices, %i faces\n", data->number_of_points,
                        data->number_of_triangles);
}
