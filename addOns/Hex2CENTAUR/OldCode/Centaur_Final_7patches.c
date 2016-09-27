// imposing strict symmetry in geometry and cell numbering
// changing code so that coordinates are generate by reflecting
//  and rotating a small block.
// making code more efficient in memory usage
// changing data structure
/**********************************************************************/
/*	Code to generate TFI grid in 3D domain with spherical hole    */
/*	Author 	: Manoj Kumar Parmar    			      */
/*	Date 	: August 3rd, 2011 	        		      */
/**********************************************************************/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

//#include <iostream>  // I/O
//#include <fstream>   // file I/O
//#include <iomanip>   // format manipulation
//#include <string>

#define PI 3.1415926535897932384626433832792
#define ZERO .000005
#define TOLERANCE 0.050
#define NO_OF_COLUMNS 5
#define XY_PLANE 1
#define YZ_PLANE 2
#define XZ_PLANE 3
#define YZ_Y_PLANE 4
#define YZ_Z_PLANE 5
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define NEXT(x)   (x==3 ? 0 : x+1)
#define PREV(x)	  (x==0 ? 3 : x-1)

typedef struct
{
	double  	x, y, z ;
} vertex ;

typedef struct
{
	int Imax ;
	vertex	Vertices[2] ;
	vertex 	*grid ;
	int *vertexno ;
	int	gridded ;   /** ==1 when TFI is done **/

	double  iclust ;   /** used for clustering **/
	int	iscircle ; /** ==1 when edge is circular section **/
	vertex	centre ;
	double	radius ;
} edge ;

typedef struct
{
	int Imax, Jmax ;
	vertex 	Vertices[4] ;
	edge 	Edges[4] ;
	vertex 	**grid ;
	int **vertexno ;
	int	gridded ; /** ==1 when TFI is done **/
	double  iclust, jclust ;   /** used for clustering **/
} face ;

typedef struct
{
	int Imax, Jmax, Kmax ;
	vertex 	Vertices[8] ;
	edge 	Edges[12] ;
	face	Faces[6] ;
	vertex 	***grid ;
	int ***vertexno ;
	int	gridded ; /** ==1 when TFI is done **/
	double  iclust, jclust, kclust ;   /** used for clustering **/

        /* MULTI-BLOCK FEATURES */
	int	direction[6] ;
	int	shared_face[6] ;
	int	neighbor[6] ;
} block ;

typedef struct
{
	int 	bType ;
	int 	nTris ;
	int	nQuads ;
	int	**Quad2v ;
	char	bName[80] ;
} boundary ;

typedef struct
{
	int     nNodes, nElems, nFaces, nFaces1, nFaces2, nFaces3, nFaces4, nFaces5, nFaces6, nFaces0 ;
	vertex 	*nodes ;
    int     **elems ;
    int     **faces, **faces1, **faces2, **faces3, **faces4, **faces5, **faces6, **faces0 ;
	int     *faceMarkers ;
} tetmesh ;

/* GLOBAL VARIABLE DECLARATIONS */
tetmesh 	Tetmesh ;
block 		*Blocks ;
int		Block_num=6 ; /* there are 6 blocks around sphere */
int             mode=1 ; /* =0 automatic, !=0 manual */
int             partmode=0 ; /* =0 6 blocks, !=0 24 blocks */
int		tetmeshFlag=0 ;
int             stretchmode=0 ; /* =0 default , !=0 user specified spacing */
int             **Cell2v ;
vertex          *Vert2xyz ;
boundary	*Boundaries ;
int		Boundary_num=2 ;
int		nBTris=0, nBQuads=0 ;
int		nCells=0, nVerts=0 ;
int             **Cell2vsingleblock ;
vertex          *Vert2xyzsingleblock ;
boundary	*Boundariessingleblock ;
int		Boundary_num_singleblock=6 ;
int		nBTrissingleblock=0, nBQuadssingleblock=0 ;
int		nCellssingleblock=0, nVertssingleblock=0 ;
double          *LinearMesh ;

/* PARAMETERS FOR DOMAIN */
int     Imax=11, Jmax=11, Kmax=11 ; /* Kmax in radial **/
double	KStretchingRatio=1.0 ;
int     KMeshType = 1 ; /*=0 automatically adjust radial spacing, =1 used specified spacing  */
int     topoFlag = 1 ; /*=0 sector with plane walls, =1 sector in spherical coordinates  */
double 	X,Y,Z,X1,Y1,Z1,YMIN=0.0,YMAX=6.65E-03;
double  D1_D=10.0 ;
double  D=1.0, R=0.5, D1=10.0, R1=5.0 ;
int     nMultiples=4, mMultiples=4 ;

/* DECLARATIONS FOR FUNCTIONS */
void    compute_Kmax(void) ;
double	find_stretching(double, double, int, double) ;
void    generate_3Dmesh(void) ;
void    generate_linearmesh(void) ;
int	generate_TFI_block(int) ;
int	generate_TFI_block_edge(int, int) ;
int	generate_TFI_edge(int, int, int) ;
int	generate_TFI_face(int, int) ;
void    impose_symmetry(void) ;
void   	initialize_flags(void) ;
void	initialize_boundary_patches(void) ;
void	initialize_boundary_patches_singleblock(void) ;
void 	initialize_edges_from_faces(int blockno, int faceno) ;
void 	initialize_from_block_gridpoints(void) ;
void 	initialize_from_block_vertices(void) ;
void 	initialize_from_face_gridpoints(void) ;
void 	initialize_from_face_vertices(void) ;
void   	initialize_meshdata(void) ;
void	number_block_vertices(void) ;
void	number_block_vertices_singleblock(void) ;
void	output_centaur(void) ;
void	output_centaur_singleblock(void) ;
void 	output_centaur_tetmesh(void) ;
void 	output_centaur_tethexmesh(void) ;
void    output_tetgen_singleblock(void) ;
void 	output_single_block(int) ;
void	output_tecplot(void) ;
int	read_inputs(int, char**) ;
void    read_tetmesh(void) ;
void   	reflect_block(int, int, int) ;
void	renumber_tetverts(void) ;
void    scale_block_edgegrid(int, int, double) ;
void    scale_facegrid(int, int, double) ;
void    scale_face_edgegrid(int, int, int, double) ;
void	testcase(void) ;

void main()
{
  read_tetmesh();
  output_centaur_tetmesh();
  //getchar();
}



void read_tetmesh()
{
	int	i,j,marker,dummy ;
	int	nNodes,nDims,nAttrs,nMarkers ;
	int 	nElems,nNpT,nRegions ;
	int     nFaces, count1, count2, count3, count4, count5, count6, count0 ;
	int	v1,v2,v3,v4 ;
	double	x,y,z ;
	char	casename[25],filename[30] ;
	char    line[128] ;
	FILE	*fp ;

	printf("Let's read tetrahedral mesh \n") ;

	sprintf(casename,"Grid_HigherMachN_V8.1") ;

	/* Read nodes ****************************************************/
  	sprintf(filename,"%s.node",casename) ;
	printf("filename= %s \n",filename) ;
      	fp = fopen(filename,"r") ;

	fscanf(fp,"%d %d %d %d",&nNodes,&nDims,&nAttrs,&nMarkers) ;
	printf("%d %d %d %d \n",nNodes,nDims,nAttrs,nMarkers) ;

	Tetmesh.nNodes = nNodes ;

	// Allocate for nodes
	Tetmesh.nodes = (vertex *) malloc(nNodes*sizeof(vertex)) ;

	for(i=0; i<nNodes; i++)
	{
	  fscanf(fp,"%d %lf %lf %lf %d",&dummy,&x,&y,&z,&dummy) ;
//	  printf("%d %lf %lf %lf %d\n",dummy,x,y,z) ;
          Tetmesh.nodes[i].x = x ;
          Tetmesh.nodes[i].y = y ;
          Tetmesh.nodes[i].z = z ;
	}

	fclose(fp) ;

	/* Read elements *************************************************/
  	sprintf(filename,"%s.ele",casename) ;
	printf("filename= %s \n",filename) ;
      	fp = fopen(filename,"r") ;

	fscanf(fp,"%d %d %d",&nElems,&nNpT,&nRegions) ;
	printf("%d %d %d \n",nElems,nNpT,nRegions) ;

	Tetmesh.nElems = nElems ;

	// Allocate for elements
	Tetmesh.elems = (int **) malloc(nElems*sizeof(int *)) ;

	for(i=0; i<nElems; i++)
	{
	  fscanf(fp,"%d %d %d %d %d",&dummy,&v1,&v2,&v3,&v4) ;
//	  printf("%d %d %d %d %d \n",dummy,v1,v2,v3,v4) ;
	  Tetmesh.elems[i] = (int *) malloc(4*sizeof(int)) ;
          Tetmesh.elems[i][0] = v2 ;
          Tetmesh.elems[i][1] = v1 ;
          Tetmesh.elems[i][2] = v3 ;
          Tetmesh.elems[i][3] = v4 ;
	}

	fclose(fp) ;

	/* Read faces ****************************************************/
  	sprintf(filename,"%s.face",casename) ;
	printf("filename= %s \n",filename) ;
      	fp = fopen(filename,"r") ;

	fscanf(fp,"%d %d",&nFaces,&nMarkers) ;
	printf("%d %d \n",nFaces,nMarkers) ;

	Tetmesh.nFaces = nFaces ;
	Tetmesh.nFaces1 = 0 ;
	Tetmesh.nFaces2 = 0 ;
	Tetmesh.nFaces3 = 0 ;
	Tetmesh.nFaces4 = 0 ;
	Tetmesh.nFaces5 = 0 ;
	Tetmesh.nFaces6 = 0 ;
        Tetmesh.nFaces0 = 0 ;

	// Allocate for faces
	Tetmesh.faces = (int **) malloc(nFaces*sizeof(int *)) ;
	Tetmesh.faceMarkers = (int *) malloc(nFaces*sizeof(int)) ;

	for(i=0; i<nFaces; i++)
	{
	  fscanf(fp,"%d %d %d %d %d",&dummy,&v1,&v2,&v3,&marker) ;
//	  printf("%d %d %d %d %d \n",dummy,v1,v2,v3,marker) ;
	  Tetmesh.faces[i] = (int *) malloc(3*sizeof(int)) ;
          Tetmesh.faces[i][0] = v2 ;
          Tetmesh.faces[i][1] = v1 ;
          Tetmesh.faces[i][2] = v3 ;
          Tetmesh.faceMarkers[i] = marker ;
	  if(marker==1)
	  {
	    Tetmesh.nFaces1++ ;
	  }
	  else if(marker==2)
	  {
	    Tetmesh.nFaces2++ ;
	  }
	  else if(marker==3)
      {
	    Tetmesh.nFaces3++ ;
	  }
	  else if(marker==4)
      {
	    Tetmesh.nFaces4++ ;
	  }
	  else if(marker==5)
      {
	    Tetmesh.nFaces5++ ;
	  }
	  else if(marker==6)
      {
	    Tetmesh.nFaces6++ ;
	  }
	  else if(marker==0)
      {
	    Tetmesh.nFaces0++ ;
	  }

}

	fclose(fp) ;

	printf("nFaces1=%d,   nFaces2=%d \n",Tetmesh.nFaces1,Tetmesh.nFaces2) ;
	printf("nFaces3=%d,   nFaces4=%d \n",Tetmesh.nFaces3,Tetmesh.nFaces4) ;
	printf("nFaces5=%d,   nFaces6=%d \n",Tetmesh.nFaces5,Tetmesh.nFaces6) ;
	printf("nFaces0=%d\n",Tetmesh.nFaces0) ;

	// Allocate for faces1,faces2,faces3,faces4,faces5,faces6,faces0
        Tetmesh.faces1 = (int **) malloc(Tetmesh.nFaces1*sizeof(int *)) ;
	Tetmesh.faces2 = (int **) malloc(Tetmesh.nFaces2*sizeof(int *)) ;
	Tetmesh.faces3 = (int **) malloc(Tetmesh.nFaces3*sizeof(int *)) ;
	Tetmesh.faces4 = (int **) malloc(Tetmesh.nFaces4*sizeof(int *)) ;
	Tetmesh.faces5 = (int **) malloc(Tetmesh.nFaces5*sizeof(int *)) ;
	Tetmesh.faces6 = (int **) malloc(Tetmesh.nFaces6*sizeof(int *)) ;
	Tetmesh.faces0 = (int **) malloc(Tetmesh.nFaces0*sizeof(int *)) ;

        count1 = 0 ; count2 = 0 ; count3 = 0;
	count4 = 0 ; count5 = 0 ; count6 = 0 ; count0 = 0;

	for(i=0; i<nFaces; i++)
	{
          marker = Tetmesh.faceMarkers[i] ;
          //printf("Assigning faces to respective boundaries i=%d, marker=%d\n",i,marker) ;

	  if(marker==1)
	  {
	    Tetmesh.faces1[count1] = (int *) malloc(3*sizeof(int)) ;
	    Tetmesh.faces1[count1][0] = Tetmesh.faces[i][0] ;
	    Tetmesh.faces1[count1][1] = Tetmesh.faces[i][1] ;
	    Tetmesh.faces1[count1][2] = Tetmesh.faces[i][2] ;
	    count1++ ;

	  }
	  else if(marker==2)
	  {
	    Tetmesh.faces2[count2] = (int *) malloc(3*sizeof(int)) ;
	    Tetmesh.faces2[count2][0] = Tetmesh.faces[i][0] ;
	    Tetmesh.faces2[count2][1] = Tetmesh.faces[i][1] ;
	    Tetmesh.faces2[count2][2] = Tetmesh.faces[i][2] ;
	    count2++ ;

	  }
	  else if(marker==3)
	  {
	    Tetmesh.faces3[count3] = (int *) malloc(3*sizeof(int)) ;
	    Tetmesh.faces3[count3][0] = Tetmesh.faces[i][0] ;
	    Tetmesh.faces3[count3][1] = Tetmesh.faces[i][1] ;
	    Tetmesh.faces3[count3][2] = Tetmesh.faces[i][2] ;
	    count3++ ;

	  }
	  else if(marker==4)
	  {
	    Tetmesh.faces4[count4] = (int *) malloc(3*sizeof(int)) ;
	    Tetmesh.faces4[count4][0] = Tetmesh.faces[i][0] ;
	    Tetmesh.faces4[count4][1] = Tetmesh.faces[i][1] ;
	    Tetmesh.faces4[count4][2] = Tetmesh.faces[i][2] ;
	    count4++ ;

	  }
	  else if(marker==5)
	  {
	    Tetmesh.faces5[count5] = (int *) malloc(3*sizeof(int)) ;
	    Tetmesh.faces5[count5][0] = Tetmesh.faces[i][0] ;
	    Tetmesh.faces5[count5][1] = Tetmesh.faces[i][1] ;
	    Tetmesh.faces5[count5][2] = Tetmesh.faces[i][2] ;
	    count5++ ;
	  }
	  else if(marker==6)
	  {
	    Tetmesh.faces6[count6] = (int *) malloc(3*sizeof(int)) ;
	    Tetmesh.faces6[count6][0] = Tetmesh.faces[i][0] ;
	    Tetmesh.faces6[count6][1] = Tetmesh.faces[i][1] ;
	    Tetmesh.faces6[count6][2] = Tetmesh.faces[i][2] ;
	    count6++ ;
	  }
	  else if(marker==0)
	  {
	    Tetmesh.faces0[count0] = (int *) malloc(3*sizeof(int)) ;
	    Tetmesh.faces0[count0][0] = Tetmesh.faces[i][0] ;
	    Tetmesh.faces0[count0][1] = Tetmesh.faces[i][1] ;
	    Tetmesh.faces0[count0][2] = Tetmesh.faces[i][2] ;
	    count0++ ;
	  }


	}

		printf("count0=%d  count1=%d  count2=%d  count3=%d  count4=%d  count5=%d  count6=%d\n",count0,count1,count2,count3,count4,count5,count6);
		printf("Face0= %d  Face1= %d  Face2= %d  Face3= %d  Face4= %d  Face5= %d  Face6= %d\n",Tetmesh.nFaces0,Tetmesh.nFaces1,Tetmesh.nFaces2,Tetmesh.nFaces3,Tetmesh.nFaces4,Tetmesh.nFaces5,Tetmesh.nFaces6);

	//if(count1 /= Tetmesh.nFaces1)
	//{
	//	printf("count0=%d  count1=%d  count2=%d  count3=%d  count4=%d  count5=%d  count6=%d\n",count0,count1,count2,count3,count4,count5,count6);
	//	printf("Face0= %d  Face1= %d  Face2= %d  Face3= %d  Face3= %d  Face4= %d  Face5= %d\n",Tetmesh.nFaces0,Tetmesh.nFaces1,Tetmesh.nFaces2,Tetmesh.nFaces3,Tetmesh.nFaces4,Tetmesh.nFaces5,Tetmesh.nFaces6);
 //	  printf("There is problem with face1 counting; count1=%d, nFaces1=%d\n",count1,Tetmesh.nFaces1) ;
	//  //return ;
	//}

	//if(count2 /= Tetmesh.nFaces2)
	//{
	//	printf("count0=%d  count1=%d  count2=%d  count3=%d  count4=%d  count5=%d  count6=%d\n",count0,count1,count2,count3,count4,count5,count6);
	//	printf("Face0= %d  Face1= %d  Face2= %d  Face3= %d  Face3= %d  Face4= %d  Face5= %d\n",Tetmesh.nFaces0,Tetmesh.nFaces1,Tetmesh.nFaces2,Tetmesh.nFaces3,Tetmesh.nFaces4,Tetmesh.nFaces5,Tetmesh.nFaces6);
 //	  printf("There is problem with face2 counting; count2=%d, nFaces2=%d\n",count2,Tetmesh.nFaces2) ;
	//  //return ;
	//}
	//if(count3 /= Tetmesh.nFaces3)
	//{
	//	printf("count0=%d  count1=%d  count2=%d  count3=%d  count4=%d  count5=%d  count6=%d\n",count0,count1,count2,count3,count4,count5,count6);
	//	printf("Face0= %d  Face1= %d  Face2= %d  Face3= %d  Face3= %d  Face4= %d  Face5= %d\n",Tetmesh.nFaces0,Tetmesh.nFaces1,Tetmesh.nFaces2,Tetmesh.nFaces3,Tetmesh.nFaces4,Tetmesh.nFaces5,Tetmesh.nFaces6);
 //	  printf("There is problem with face3 counting; count3=%d, nFaces3=%d\n",count3,Tetmesh.nFaces3) ;
	//  //return ;
	//}
	//if(count4 /= Tetmesh.nFaces4)
	//{
	//	printf("count0=%d  count1=%d  count2=%d  count3=%d  count4=%d  count5=%d  count6=%d\n",count0,count1,count2,count3,count4,count5,count6);
	//	printf("Face0= %d  Face1= %d  Face2= %d  Face3= %d  Face3= %d  Face4= %d  Face5= %d\n",Tetmesh.nFaces0,Tetmesh.nFaces1,Tetmesh.nFaces2,Tetmesh.nFaces3,Tetmesh.nFaces4,Tetmesh.nFaces5,Tetmesh.nFaces6);
 //	  printf("There is problem with face4 counting; count4=%d, nFaces4=%d\n",count4,Tetmesh.nFaces4) ;
	//  //return ;
	//}
	//if(count5 /= Tetmesh.nFaces5)
	//{
	//	printf("count0=%d  count1=%d  count2=%d  count3=%d  count4=%d  count5=%d  count6=%d\n",count0,count1,count2,count3,count4,count5,count6);
	//	printf("Face0= %d  Face1= %d  Face2= %d  Face3= %d  Face3= %d  Face4= %d  Face5= %d\n",Tetmesh.nFaces0,Tetmesh.nFaces1,Tetmesh.nFaces2,Tetmesh.nFaces3,Tetmesh.nFaces4,Tetmesh.nFaces5,Tetmesh.nFaces6);
 //	  printf("There is problem with face5 counting; count5=%d, nFaces5=%d\n",count5,Tetmesh.nFaces5) ;
	//  //return ;
	//}
	//if(count6 /= Tetmesh.nFaces6)
	//{
	//	printf("count0=%d  count1=%d  count2=%d  count3=%d  count4=%d  count5=%d  count6=%d\n",count0,count1,count2,count3,count4,count5,count6);
	//	printf("Face0= %d  Face1= %d  Face2= %d  Face3= %d  Face3= %d  Face4= %d  Face5= %d\n",Tetmesh.nFaces0,Tetmesh.nFaces1,Tetmesh.nFaces2,Tetmesh.nFaces3,Tetmesh.nFaces4,Tetmesh.nFaces5,Tetmesh.nFaces6);
 //	  printf("There is problem with face6 counting; count6=%d, nFaces6=%d\n",count6,Tetmesh.nFaces6) ;
	//  //return ;
	//}
	//if(count0 /= Tetmesh.nFaces0)
	//{
	//	printf("count0=%d  count1=%d  count2=%d  count3=%d  count4=%d  count5=%d  count6=%d\n",count0,count1,count2,count3,count4,count5,count6);
	//	printf("Face0= %d  Face1= %d  Face2= %d  Face3= %d  Face3= %d  Face4= %d  Face5= %d\n",Tetmesh.nFaces0,Tetmesh.nFaces1,Tetmesh.nFaces2,Tetmesh.nFaces3,Tetmesh.nFaces4,Tetmesh.nFaces5,Tetmesh.nFaces6);
 //	  printf("There is problem with face0 counting; count0=%d, nFaces0=%d\n",count0,Tetmesh.nFaces0) ;
	//  //return ;
	//}

  //  for(i=0; i<9999999999; i++)
//	{
//		marker = i;
//	}
}

void output_centaur_tetmesh()
{
	int	i, idum, j, boundno ;
	int	ii, jj, kk, count, quads_count,zero,const1,const2,const3,const4,const5,const6 ;
	int	v1, v2, v3, v4, v5, v6, v7, v8 ;
	char	casename[20], filename[20] ;
	FILE	*fp ;

//*************************************************************
// 	Open file and read title
//*************************************************************

  	printf("Writing tetrahedral CENTAUR grid file...\n") ;
	sprintf(casename,"cylds") ;

  	sprintf(filename,"%s.hyb.asc",casename) ;
  	printf("writing tetrahedral centaur format in \"%s\" \n",filename) ;
      	fp = fopen(filename,"w") ;

	fprintf(fp,"conflu\n") ;

//===================================================================
//	Coordinates
//===================================================================

  	fprintf(fp,"%8d\n",Tetmesh.nNodes) ;

  	if(Tetmesh.nNodes > 0)
	{
  	    printf("Coordinates...\n") ;
            count = 0 ;
            for(i=0; i<Tetmesh.nNodes; i++)
            {
               count++ ;
               if(count > NO_OF_COLUMNS){ fprintf(fp,"\n") ; count = 1 ; }
               fprintf(fp,"%23.16E",Tetmesh.nodes[i].x) ;
               

            }
	    fprintf(fp,"\n") ;
            count = 0 ;
            for(i=0; i<Tetmesh.nNodes; i++)
            {
               count++ ;
               if(count > NO_OF_COLUMNS){ fprintf(fp,"\n") ; count = 1 ; }
               fprintf(fp,"%23.16E",Tetmesh.nodes[i].y) ;
               

            }
	    fprintf(fp,"\n") ;
            count = 0 ;
            for(i=0; i<Tetmesh.nNodes; i++)
            {
               count++ ;
               if(count > NO_OF_COLUMNS){ fprintf(fp,"\n") ; count = 1 ; }
               fprintf(fp,"%23.16E",Tetmesh.nodes[i].z) ;
               

            }
	    fprintf(fp,"\n") ;
  	    idum = 0 ; count = 0 ;
            for(i=0; i<Tetmesh.nNodes; i++)
	    {
		count++ ;
                if(count > 2*NO_OF_COLUMNS){ fprintf(fp,"\n") ; count = 1 ; }
  		fprintf(fp,"%8d",idum) ;
      	    }
	    fprintf(fp,"\n") ;
 	}

//=====================================================================
//	Cell connectivity
//=====================================================================

	/* writing tetrahedrals *********/
	fprintf(fp,"%8d\n",Tetmesh.nElems) ;
  	if(Tetmesh.nElems > 0)
	{
    	    printf("Tetrahedrals...\n") ;
            for(j=0; j<4; j++)
            {
                count = 0 ;
                for(i=0; i<Tetmesh.nElems; i++)
                {
                   count++ ;
               	   if(count > 2*NO_OF_COLUMNS){ fprintf(fp,"\n") ; count = 1 ; }
                   fprintf(fp,"%8d",Tetmesh.elems[i][j]) ;
                }
		fprintf(fp,"\n") ;
            }
  	    idum = 0 ; count = 0 ;
            for(ii=0; ii<Tetmesh.nElems; ii++)
	    {
		count++ ;
                if(count > 2*NO_OF_COLUMNS){ fprintf(fp,"\n") ; count = 1 ; }
  		fprintf(fp,"%8d",idum) ;
      	    }
	    fprintf(fp,"\n") ;
 	}

	/* writing hexahedrals *********/
if(1==2)
{
	fprintf(fp,"%8d\n",nCellssingleblock) ;
  	if(nCellssingleblock > 0)
	{
    	    printf("Hexahedrals...\n") ;
            for(j=0; j<8; j++)
            {
                count = 0 ;
                for(i=0; i<nCellssingleblock; i++)
                {
                   count++ ;
               	   if(count > 2*NO_OF_COLUMNS){ fprintf(fp,"\n") ; count = 1 ; }
                   fprintf(fp,"%8d",Cell2vsingleblock[i][j]) ;
                }
		fprintf(fp,"\n") ;
            }
  	    idum = 0 ; count = 0 ;
            for(ii=0; ii<nCellssingleblock; ii++)
	    {
		count++ ;
                if(count > 2*NO_OF_COLUMNS){ fprintf(fp,"\n") ; count = 1 ; }
  		fprintf(fp,"%8d",idum) ;
      	    }
	    fprintf(fp,"\n") ;
 	}
}
else
{
	fprintf(fp,"%8d\n",0) ;
}

	/* writing prisms *************/
	fprintf(fp,"%8d\n",0) ;

	/* writing pyramids **************/
	fprintf(fp,"%8d\n",0) ;

//==============================================================================
//	Boundary types
//==============================================================================

  	printf("Boundary information...\n") ;

	fprintf(fp,"%8d\n",7) ; // For 7 boundaries

	/* writing Boundary types for all boundaries ***/
	count = 0 ;
	for(i=0; i<7; i++)
	{
		count++ ;
            	if(count <= 2*NO_OF_COLUMNS)
			fprintf(fp,"%8d",400) ;
             	else
                {
                	count = 1 ;
                	fprintf(fp,"\n") ;
			fprintf(fp,"%8d",400) ;
                }
	}
	if(count <= 2*NO_OF_COLUMNS) fprintf(fp,"\n") ;

	/* writing Boundary no of triangles */
	count = 0 ;
        const1 = 0 ;
        const2 = 0 ;
        const3 = 0 ;
        const4 = 0 ;
        const5 = 0 ;
        const6 = 0 ;
        const1 = Tetmesh.nFaces0+Tetmesh.nFaces1;
        const2 = const1+Tetmesh.nFaces2;
        const3 = const2+Tetmesh.nFaces3;
        const4 = const3+Tetmesh.nFaces4;
        const5 = const4+Tetmesh.nFaces5;
        const6 = const5+Tetmesh.nFaces6;
	for(i=0; i<7; i++)
	{
		count++ ;
            	if(count <= 2*NO_OF_COLUMNS)
                {
                if (i==0)
                {
                 fprintf(fp,"%8d",Tetmesh.nFaces0) ;
                }
                else if (i==1)
                {
                fprintf(fp,"%8d",const1) ;
                }
                else if (i==2)
                {
                fprintf(fp,"%8d",const2) ;
                }
                else if (i==3)
                {
                fprintf(fp,"%8d",const3) ;
                }
                else if (i==4)
                {
                fprintf(fp,"%8d",const4) ;
                }
                else if (i==5)
                {
                fprintf(fp,"%8d",const5) ;
                }
                else if (i==6)
                {
                fprintf(fp,"%8d",const6) ;
                }
                else
                {
                count = 1 ;
                fprintf(fp,"\n") ;
	        if (i==0)
                {
                 fprintf(fp,"%8d",Tetmesh.nFaces0) ;
                }
                else if (i==1)
                {
                fprintf(fp,"%8d",const1) ;
                }
                else if (i==2)
                {
                fprintf(fp,"%8d",const2) ;
                }
                else if (i==3)
                {
                fprintf(fp,"%8d",const3) ;
                }
                else if (i==4)
                {
                fprintf(fp,"%8d",const4) ;
                }
                else if (i==5)
                {
                fprintf(fp,"%8d",const5) ;
                }
                else if (i==6)
                {
                fprintf(fp,"%8d",const6) ;
                }
                }
                }
	}
	if(count <= 2*NO_OF_COLUMNS) fprintf(fp,"\n") ;

/* writing Boundary no of Quadrilaterals */
	count = 0 ;
        zero = 0;
	for(i=0; i<7; i++)
	{
		count++ ;
            	if(count <= 2*NO_OF_COLUMNS)
                {
                if (i==0)
                {
                 fprintf(fp,"%8d",zero) ;
                }
                else if (i==1)
                {
                fprintf(fp,"%8d",zero) ;
                }
                else if (i==2)
                {
                fprintf(fp,"%8d",zero) ;
                }
                else if (i==3)
                {
                fprintf(fp,"%8d",zero) ;
                }
                else if (i==4)
                {
                fprintf(fp,"%8d",zero) ;
                }
                else if (i==5)
                {
                fprintf(fp,"%8d",zero) ;
                }
                else if (i==6)
                {
                fprintf(fp,"%8d",zero) ;
                }
                else
                {
                count = 1 ;
                fprintf(fp,"\n") ;
	        if (i==0)
                {
                 fprintf(fp,"%8d",zero) ;
                }
                else if (i==1)
                {
                fprintf(fp,"%8d",zero) ;
                }
                else if (i==2)
                {
                fprintf(fp,"%8d",zero) ;
                }
                else if (i==3)
                {
                fprintf(fp,"%8d",zero) ;
                }
                else if (i==4)
                {
                fprintf(fp,"%8d",zero) ;
                }
                else if (i==5)
                {
                fprintf(fp,"%8d",zero) ;
                }
                else if (i==6)
                {
                fprintf(fp,"%8d",zero) ;
                }
                }
                }
	}
	if(count <= 2*NO_OF_COLUMNS) fprintf(fp,"\n") ;


	/* writing Boundary names */
		fprintf(fp,"%s\n","Inflow") ;
		fprintf(fp,"%s\n","Slip-wall") ;
		fprintf(fp,"%s\n","Slip-wall") ;
		fprintf(fp,"%s\n","Slip-wall") ;
		fprintf(fp,"%s\n","Slip-wall") ;
		fprintf(fp,"%s\n","Slip-wall") ;
		fprintf(fp,"%s\n","Slip-wall") ;
if(1==2)
{
	fprintf(fp,"%8d\n",Boundary_num_singleblock) ;

	/* writing Boundary types for all boundaries ***/
	count = 0 ;
	for(i=0; i<Boundary_num_singleblock; i++)
	{
		count++ ;
            	if(count <= 2*NO_OF_COLUMNS)
			fprintf(fp,"%8d",Boundariessingleblock[i].bType) ;
             	else
                {
                	count = 1 ;
                	fprintf(fp,"\n") ;
			fprintf(fp,"%8d",Boundariessingleblock[i].bType) ;
                }
	}
	if(count <= 2*NO_OF_COLUMNS) fprintf(fp,"\n") ;

	/* writing Boundary no of triangles */
	count = 0 ;
	for(i=0; i<Boundary_num_singleblock; i++)
	{
		count++ ;
            	if(count <= 2*NO_OF_COLUMNS)
			fprintf(fp,"%8d",Boundariessingleblock[i].nTris) ;
             	else
                {
                	count = 1 ;
                	fprintf(fp,"\n") ;
			fprintf(fp,"%8d",Boundariessingleblock[i].nTris) ;
                }
	}
	if(count <= 2*NO_OF_COLUMNS) fprintf(fp,"\n") ;

	/* writing Boundary no of Quadrilaterals */
	count = 0 ; quads_count = 0 ;
	for(i=0; i<Boundary_num_singleblock; i++)
	{
		count++ ;
            	if(count <= 2*NO_OF_COLUMNS)
		{
			quads_count += Boundariessingleblock[i].nQuads ;
			fprintf(fp,"%8d",quads_count) ;
             	}
		else
                {
                	count = 1 ;
                	fprintf(fp,"\n") ;
			quads_count += Boundariessingleblock[i].nQuads ;
			fprintf(fp,"%8d",quads_count) ;
                }
	}
	if(count <= 2*NO_OF_COLUMNS) fprintf(fp,"\n") ;

	/* writing Boundary names */
	for(i=0; i<Boundary_num_singleblock; i++)
		fprintf(fp,"%s\n",Boundariessingleblock[i].bName) ;
}

//==============================================================================
//    	Boundary face connectivity
//==============================================================================

	fprintf(fp,"%8d\n",const6) ;
        for(j=0; j<3; j++)
        {
            count = 0 ;
  	    if(Tetmesh.nFaces0>0)
	    {
    	       /* Tris of each boundary Patch */
	       for(boundno=0; boundno<1; boundno++)
	       {
               /* tris of boundno th boundary */
	           for(i=0; i<Tetmesh.nFaces0; i++)
	           {
		                count++ ;
		        	if(count <= 2*NO_OF_COLUMNS)
				  fprintf(fp,"%8d",Tetmesh.faces0[i][j]) ;
				else
				{
				  count = 1 ;
				  fprintf(fp,"\n") ;
				  fprintf(fp,"%8d",Tetmesh.faces0[i][j]) ;
				}
		   }
		}
             }
            if(Tetmesh.nFaces1>0)
	    {
    	       /* Tris of each boundary Patch */
	       for(boundno=0; boundno<1; boundno++)
	       {
               /* tris of boundno th boundary */
	           for(i=0; i<Tetmesh.nFaces1; i++)
	           {
		                count++ ;
		        	if(count <= 2*NO_OF_COLUMNS)
				  fprintf(fp,"%8d",Tetmesh.faces1[i][j]) ;
				else
				{
				  count = 1 ;
				  fprintf(fp,"\n") ;
				  fprintf(fp,"%8d",Tetmesh.faces1[i][j]) ;
				}
		   }
		}
             }
             if(Tetmesh.nFaces2>0)
	    {
    	       /* Tris of each boundary Patch */
	       for(boundno=0; boundno<1; boundno++)
	       {
               /* tris of boundno th boundary */
	           for(i=0; i<Tetmesh.nFaces2; i++)
	           {
		                count++ ;
		        	if(count <= 2*NO_OF_COLUMNS)
				  fprintf(fp,"%8d",Tetmesh.faces2[i][j]) ;
				else
				{
				  count = 1 ;
				  fprintf(fp,"\n") ;
				  fprintf(fp,"%8d",Tetmesh.faces2[i][j]) ;
				}
		   }
		}
             }
             if(Tetmesh.nFaces3>0)
	    {
    	       /* Tris of each boundary Patch */
	       for(boundno=0; boundno<1; boundno++)
	       {
               /* tris of boundno th boundary */
	           for(i=0; i<Tetmesh.nFaces3; i++)
	           {
		                count++ ;
		        	if(count <= 2*NO_OF_COLUMNS)
				  fprintf(fp,"%8d",Tetmesh.faces3[i][j]) ;
				else
				{
				  count = 1 ;
				  fprintf(fp,"\n") ;
				  fprintf(fp,"%8d",Tetmesh.faces3[i][j]) ;
				}
		   }
		}
             }
             if(Tetmesh.nFaces4>0)
	    {
    	       /* Tris of each boundary Patch */
	       for(boundno=0; boundno<1; boundno++)
	       {
               /* tris of boundno th boundary */
	           for(i=0; i<Tetmesh.nFaces4; i++)
	           {
		                count++ ;
		        	if(count <= 2*NO_OF_COLUMNS)
				  fprintf(fp,"%8d",Tetmesh.faces4[i][j]) ;
				else
				{
				  count = 1 ;
				  fprintf(fp,"\n") ;
				  fprintf(fp,"%8d",Tetmesh.faces4[i][j]) ;
				}
		   }
		}
             }
             if(Tetmesh.nFaces5>0)
	    {
    	       /* Tris of each boundary Patch */
	       for(boundno=0; boundno<1; boundno++)
	       {
               /* tris of boundno th boundary */
	           for(i=0; i<Tetmesh.nFaces5; i++)
	           {
		                count++ ;
		        	if(count <= 2*NO_OF_COLUMNS)
				  fprintf(fp,"%8d",Tetmesh.faces5[i][j]) ;
				else
				{
				  count = 1 ;
				  fprintf(fp,"\n") ;
				  fprintf(fp,"%8d",Tetmesh.faces5[i][j]) ;
				}
		   }
		}
             }
             if(Tetmesh.nFaces6>0)
	    {
    	       /* Tris of each boundary Patch */
	       for(boundno=0; boundno<1; boundno++)
	       {
               /* tris of boundno th boundary */
	           for(i=0; i<Tetmesh.nFaces6; i++)
	           {
		                count++ ;
		        	if(count <= 2*NO_OF_COLUMNS)
				  fprintf(fp,"%8d",Tetmesh.faces6[i][j]) ;
				else
				{
				  count = 1 ;
				  fprintf(fp,"\n") ;
				  fprintf(fp,"%8d",Tetmesh.faces6[i][j]) ;
				}
		   }
		}
             }
             fprintf(fp,"\n") ;	
	}
     

        


if(1==2)
{
	fprintf(fp,"%8d\n",nBQuadssingleblock) ;
  	if(nBQuadssingleblock>0)
	{
    		printf("Boundary quadrilaterals...\n") ;
		/* node j of each quadrilateral in all patches */
		for(j=0; j<4; j++)
		{
			count = 0 ;
			/* Quads of each boundary Patch */
			for(boundno=0; boundno<Boundary_num_singleblock; boundno++)
			{
			  /* quads of boundno th boundary */
			  for(i=0; i<Boundariessingleblock[boundno].nQuads; i++)
			  {
				count++ ;
		        	if(count <= 2*NO_OF_COLUMNS)
				  fprintf(fp,"%8d",Boundariessingleblock[boundno].Quad2v[i][j]) ;
				else
				{
				  count = 1 ;
				  fprintf(fp,"\n") ;
				  fprintf(fp,"%8d",Boundariessingleblock[boundno].Quad2v[i][j]) ;
				}
			  }
			}
        		if(count <= 2*NO_OF_COLUMNS) fprintf(fp,"\n") ;
		}
	}
}

else
{
	fprintf(fp,"%8d\n",0) ;
}

//**************************************************************************
// 	Close file
//**************************************************************************

	fclose(fp) ;

  	printf("Writing tetrahedral CENTAUR grid file done.\n") ;

  	return ;
}

