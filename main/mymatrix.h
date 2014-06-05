namespace psi{ namespace main{
        
    void init_matrix2d(double**& matrix, int size1, int size2);
    void init_matrix4d(double****& matrix, int size1, int size2, int size3, int size4);

    void free_matrix2d(double**& matrix, int size1, int size2);
    void free_matrix4d(double****& matrix, int size1, int size2, int size3, int size4);

    void print_matrix2d(double**& matrix, int size1, int size2);
    void print_matrix4d(double****& matrix, int size1, int size2, int size3, int size4);

}} //End namespace
        
