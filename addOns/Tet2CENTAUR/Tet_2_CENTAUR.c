/**********************************************************************/
/*	Purposse: Code to convert unstructured tetrahedral mesh to    */
/*	          CENTAUR ascii format.                               */
/*	                                                              */
/*	Author 	: Yash Mehta, Christopher Neal			      */
/*	Date 	: July 16th, 2014                                     */
/*	Updated : January 15th, 2015	        		      */
/**********************************************************************/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>


#define NO_OF_COLUMNS 5 //Hardcode for CENTAUR File Format

typedef struct
{
	double  	x, y, z ; // Coordinates of vertex
} vertex ;


typedef struct
{
	int     nNodes, nElems, nFaces,nMarkers;
	vertex 	*nodes ;
   	int     **elems ;
    	int     **faces;
	int     *faceMarkers ;
	int	*faceCounts;
	int 	*boundaryMarkers;
} tetmesh ;

/* GLOBAL VARIABLE DECLARATIONS */
tetmesh 	Tetmesh ;
const char *FLOAT_OUTPUT_FORMAT_SPEC = "%16.9E"; // field width is 26, 16 decimal digits, exponential notation with big E
const char *INT_OUPTUT_FORMAT_SPEC = "%16d"; // Field width is 8, number is integer
/* DECLARATIONS FOR FUNCTIONS */
void 	output_centaur_tetmesh(void) ;
void    read_tetmesh(void) ;


void main()
{
  // Read in tetmesh format
  read_tetmesh();

  // Print tetmesh format in CENTAUR format
  output_centaur_tetmesh();
}



void read_tetmesh()
{

        // Local Variable Declarations
	int	i,j,marker,dummy ;
	int	nNodes,nDims,nAttrs,MarkerFlag ;
	int 	nElems,nNpT,nRegions ;
	int     nFaces ;
	int	v1,v2,v3,v4 ;
	int     alloc_Err ;
	double	x,y,z ;
	char	casename[70],filename[75] ;
	char    line[128] ;
	FILE	*fp ;

	printf("WELCOME TO TETMESH->CENTAUR CONVERTER\n\n");
	printf("The Code reads TETGEN tetrahedral meshes and outputs CENTAUR mesh format. \n\n") ;

        //Store casename provided by user
        printf("What is file casename?  i.e. <casename>.node:");
        scanf("%s",casename);

        // Echo user input
        printf("Selected casename is: %s \n",casename);

	/* Read nodes ****************************************************/
  	sprintf(filename,"%s.node",casename) ;
	printf("Reading File: %s \n\n",filename) ; //Echo filename

      	fp = fopen(filename,"r") ;

	fscanf(fp,"%d %d %d %d",&nNodes,&nDims,&nAttrs,&MarkerFlag) ; // Read first line of .node file

	// Echo data read from .node file
	printf("Node File Data\n");
	printf("# of Nodes: %d \n# of Dimensions: %d\nAttribute Flag: %d \nMarker Flag: %d \n\n",nNodes,nDims,nAttrs,MarkerFlag) ;

	Tetmesh.nNodes = nNodes ; // Store # of nodes into Tetmesh variable

	// Allocate for nodes array(1D array)
	Tetmesh.nodes = (vertex *) malloc(nNodes*sizeof(vertex)) ;

        // Loop through file and read all node data
	for(i=0; i<nNodes; i++)
	{
	
	  if(MarkerFlag == 1 )
	  {
	  	fscanf(fp,"%d %lf %lf %lf %d",&dummy,&x,&y,&z,&dummy) ;
	  	//printf("%lf %lf %lf\n",x,y,z) ; //Echo data from file
          	Tetmesh.nodes[i].x = x ;
          	Tetmesh.nodes[i].y = y ;
          	Tetmesh.nodes[i].z = z ;
	  }
	  else
	  {
		fscanf(fp,"%d %lf %lf %lf",&dummy,&x,&y,&z) ;
                //printf("%lf %lf %lf\n",x,y,z) ; //Echo data from file
                Tetmesh.nodes[i].x = x ;
                Tetmesh.nodes[i].y = y ;
                Tetmesh.nodes[i].z = z ;

	  }


	}

	fclose(fp) ;

	/* Read elements *************************************************/
  	sprintf(filename,"%s.ele",casename) ;
	printf("Reading File: %s \n\n",filename) ; // Echo filename to be read

      	fp = fopen(filename,"r") ; // Open .ele file for reading elements

	fscanf(fp,"%d %d %d",&nElems,&nNpT,&nRegions) ; // Read first line
	printf("Element File Data\n");
	printf("# of Elements: %d\n# of Verts Per Element:  %d\n# of Regions: %d \n\n",nElems,nNpT,nRegions) ; //Echo first line to user

	Tetmesh.nElems = nElems ; // Store number of elements into Tetmesh variable

	// Allocate for elements array ( 2D array)
	Tetmesh.elems = (int **) malloc(nElems*sizeof(int *)) ; // Allocate rows of array
        
        //Loop through and store element vertices
	for(i=0; i<nElems; i++)
	{
	  fscanf(fp,"%d %d %d %d %d",&dummy,&v1,&v2,&v3,&v4) ;
	  //printf("%d %d %d %d %d \n",dummy,v1,v2,v3,v4) ;
	  Tetmesh.elems[i] = (int *) malloc(4*sizeof(int)) ; // Allocate columns of array
          Tetmesh.elems[i][0] = v2 ;
          Tetmesh.elems[i][1] = v1 ;
          Tetmesh.elems[i][2] = v3 ;
          Tetmesh.elems[i][3] = v4 ;
	}

	fclose(fp) ;

	/* Read faces ****************************************************/
  	sprintf(filename,"%s.face",casename) ;
	printf("Reading File: %s \n\n",filename) ; // Echo filename being read

      	fp = fopen(filename,"r") ;

	fscanf(fp,"%d %d",&nFaces,&MarkerFlag) ;//MarkerFlag-->boundary markers flag
	printf("Face File Data\n");
	printf("# of Faces: %d\nBoundary Marker Flag: %d \n",nFaces,MarkerFlag) ; // Echo first line back to screen

	Tetmesh.nFaces = nFaces ; // Store total number of faces in mesh
	

	// Allocate for faces array( 2D array)
	Tetmesh.faces = (int **) malloc(nFaces*sizeof(int *)) ; // Allocate rows of face array
	Tetmesh.faceMarkers = (int *) malloc(nFaces*sizeof(int)) ;// Allocate 1D array for all face markers
	
	Tetmesh.nMarkers = 1;	// Initialize unique boundary marker counter
	Tetmesh.boundaryMarkers = malloc(1*sizeof(int)); // Allocate first element of unique boundary marker array

	for(i=0; i<nFaces; i++)
	{
	  fscanf(fp,"%d %d %d %d %d",&dummy,&v1,&v2,&v3,&marker) ;
	  //printf("%d %d %d %d %d \n",dummy,v1,v2,v3,marker) ;
	  Tetmesh.faces[i] = (int *) malloc(3*sizeof(int)) ; // Allocate columns for row i

          //Order of read matters for Clockwise/Anti-Clockwise ordering
          Tetmesh.faces[i][0] = v2 ;
          Tetmesh.faces[i][1] = v1 ;
          Tetmesh.faces[i][2] = v3 ;

          Tetmesh.faceMarkers[i] = marker ;
	  
	  if(i == 0) // set currentMarker to the first marker that is read from file
	  {
		Tetmesh.boundaryMarkers[i] = marker;
       	  }


	  int newMarker = 1;	//Boolean for tracking if new marker is found. 0=false, 1=true
	  for(j=0;j<Tetmesh.nMarkers;j++) // Compare currently read marker with markers that have been stored
	  {
		if(marker == Tetmesh.boundaryMarkers[j])
		{
			newMarker = 0;
		}

	  }

	  if( newMarker == 1  ) // A new unique boundary marker is found. 
	  {
		Tetmesh.nMarkers=Tetmesh.nMarkers+1;

 		// Dynamically resize boundaryMarkers array to hold newly found boundary marker
		int *temp = realloc(Tetmesh.boundaryMarkers,Tetmesh.nMarkers*sizeof(int));  // Store data in a larger temp array
   		if(temp != NULL) // realloc was successful
		{
			Tetmesh.boundaryMarkers = temp; //Point original array to newly sized array
			Tetmesh.boundaryMarkers[Tetmesh.nMarkers-1] = marker; // Store newly read marker into array
		}
		else
		{
			free(Tetmesh.boundaryMarkers);
			printf("Error allocating memory for Boundary Marker Array!\n");
			return ;
		}

	  }
		

}

	// Display the number of unique markers in .face file
	printf("\nNumber of Unique Boundary Markers: %d \n",Tetmesh.nMarkers);

	// Allocate faceCounts array
	Tetmesh.faceCounts = malloc(Tetmesh.nMarkers*sizeof(int));

	// Initialize faceCounts array to zero. This array will hold the number
	// of faces that have a specific boundary marker.
	for(i=0; i<Tetmesh.nMarkers; i++)
	{
		Tetmesh.faceCounts[i] = 0;
	}
	


	//Loop through faceMarkers array and count the number of faces that have all of the markers
	//that are stored in boundaryMarkers
	for(i=0;i<nFaces; i++) // Loop through every face
	{
		for(j=0;j<Tetmesh.nMarkers; j++) //Loop through all boundary markers
		{

			if(Tetmesh.faceMarkers[i] == Tetmesh.boundaryMarkers[j] ) // see which boundaryMarker the faceMarker matches
			{				
				Tetmesh.faceCounts[j] = Tetmesh.faceCounts[j] + 1; // Increment corresponding faceCount
			}			

		}

	}


	fclose(fp) ;


	//Re-Order the boundaryMarkers array to be in order of ascending boundary marker. Change
	// e.g. [ 1 5 4 6] --> [1 4 5 6]
	for(i=0;i<Tetmesh.nMarkers;i++)
	{
		int small = Tetmesh.boundaryMarkers[i];

		for(j=i;j<Tetmesh.nMarkers;j++)
		{
		
			if(Tetmesh.boundaryMarkers[j] < small)
			{
				//Update small
				small = Tetmesh.boundaryMarkers[j];

				//Store Information about smaller value
				int temp = Tetmesh.boundaryMarkers[j];	//Store smaller value
				int temp2 = Tetmesh.faceCounts[j];	
				
				//Swap boundaryMarkers
				Tetmesh.boundaryMarkers[j] = Tetmesh.boundaryMarkers[i];
				Tetmesh.boundaryMarkers[i] = temp;

				//Swap faceCounts
				Tetmesh.faceCounts[j] = Tetmesh.faceCounts[i];
				Tetmesh.faceCounts[i] = temp2;
			}			

		}


	}

	// Print the number of faces on each boundary that were counted in file
	printf("\nNumber of Faces On Each Boundary\n");	
	int sum = 0;
	for( i=0;i<Tetmesh.nMarkers;i++)
	{
	  	printf("Boundary Marker %d(%d) = %d \n",Tetmesh.boundaryMarkers[i],i+1,Tetmesh.faceCounts[i]);
		sum = sum + Tetmesh.faceCounts[i];  
	}

	//Error Check, Print face counts to screen
	printf("\nTotal Faces Counted:%d\n# Faces In File:%d\n\n",sum,Tetmesh.nFaces);

	printf("Finished Reading Data File\n\n");
}


void output_centaur_tetmesh()
{
	int	i, idum, j, k ;
	int	count,sum ;
	int	bin_or_ascii;
	char	casename[70], filename[70] ;
	FILE	*fp ;

//*************************************************************
// 	Open file and read title
//*************************************************************

	printf("Writing ASCII tetrahedral CENTAUR grid file...\n") ;
	
	/*Store output casename provided by user*/
	printf("What is output casename?  i.e. <casename>.hyb.asc:");
	scanf("%s",casename);
	               
	sprintf(filename,"%s.hyb.asc",casename) ;
	printf("Writing tetrahedral CENTAUR format to: \"%s\" \n",filename) ;
	fp = fopen(filename,"w") ;

	/*Check if the user wants binary or ascii data printed to file*/
	//printf("Write to ASCII or BINARY? : 0-ASCII , 1-BINARY  ");
	//scanf("%d",bin_or_ascii);
	                
/*Begin Printing to File*/

	fprintf(fp,"conflu\n") ;

//===================================================================
//	Print Node Coordinates
//===================================================================

	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.nNodes) ;
	fprintf(fp,"\n");

	printf("Printing Node Coordinates...\n") ;

	count = 0 ;
	for(i=0; i<Tetmesh.nNodes; i++) // Print all x coordinates of nodes
	{
		count++ ;
           	if(count > NO_OF_COLUMNS) // Only allow 5 coordinates per line
	   	{ 
			fprintf(fp,"\n") ;
			count = 1 ; 
	   	}

           	fprintf(fp,FLOAT_OUTPUT_FORMAT_SPEC,Tetmesh.nodes[i].x) ;//Print x coordinate of node i            

        }

	fprintf(fp,"\n") ;

	count = 0 ;
        for(i=0; i<Tetmesh.nNodes; i++) // Print all y coordinates of nodes
        {
        	count++ ;
       	        if(count > NO_OF_COLUMNS) // Only allow 5 coordinates per line 
		{
			fprintf(fp,"\n") ;
			count = 1 ; 
		}
		
		fprintf(fp,FLOAT_OUTPUT_FORMAT_SPEC,Tetmesh.nodes[i].y) ;
               

        }
	
	fprintf(fp,"\n") ;

	count = 0 ;
        for(i=0; i<Tetmesh.nNodes; i++) // Print all z coordinates of nodes
        {
           count++ ;
           	if(count > NO_OF_COLUMNS) // Only allow 5 coordinates per line
		{ 
			fprintf(fp,"\n") ; 
			count = 1 ; 
		}
	
           fprintf(fp,FLOAT_OUTPUT_FORMAT_SPEC,Tetmesh.nodes[i].z) ;
               
        }
	
	fprintf(fp,"\n") ;

  	idum = 0 ; 
	count = 0 ;
        for(i=0; i<Tetmesh.nNodes; i++) // Print a bunch of zeros(per the CENTAUR format)
	{
		count++ ;
                if(count > 2*NO_OF_COLUMNS) // Only allow 10 zeros per line
		{ 
			fprintf(fp,"\n") ; 
			count = 1 ; 
		}

  		fprintf(fp,INT_OUPTUT_FORMAT_SPEC,idum) ;
      	}
	
	fprintf(fp,"\n") ;
	

//=====================================================================
//	Cell connectivity
//=====================================================================

	/* writing tetrahedrals *********/
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.nElems) ;
	fprintf(fp,"\n");	

    	printf("Printing Tetrahedral Cell Connectivity...\n") ;
        
	for(j=0; j<4; j++)
        {
                count = 0 ;
                for(i=0; i<Tetmesh.nElems; i++) // Print nodes that make up elements
                {
                   	count++ ;
               	   	if(count > 2*NO_OF_COLUMNS) // Only print 10 element nodes per line
			{ 
				fprintf(fp,"\n") ; 
				count = 1 ; 
			}

                   	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.elems[i][j]) ;
                }

		fprintf(fp,"\n") ;
        }


  	idum = 0 ; 
	count = 0 ;
        for(i=0; i<Tetmesh.nElems; i++) // Print a bunch of zeros(Per CENTAUR format)
	{
		count++ ;
                if(count > 2*NO_OF_COLUMNS) // Only allow 10 zeros per line
		{ 
			fprintf(fp,"\n") ; 
			count = 1 ; 
		}

  		fprintf(fp,INT_OUPTUT_FORMAT_SPEC,idum) ;
      	}
	
	fprintf(fp,"\n") ;

	/* writing hexahedras *************/
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,0) ;
        fprintf(fp,"\n");

	/* writing prisms *************/
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,0) ;
        fprintf(fp,"\n");

	/* writing pyramids **************/
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,0) ;
	fprintf(fp,"\n");

//==============================================================================
//	Boundary types
//==============================================================================

  	printf("Printing Boundary information...\n") ;

	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.nMarkers) ; //Print number of boundaries to file
	fprintf(fp,"\n");

	printf("Number of Boundaries: %d\n",Tetmesh.nMarkers) ; //Print number of boundaries to screen

	/* writing Boundary types for all boundaries ***/
	count = 0 ;
	for(i=0; i<Tetmesh.nMarkers; i++)
	{
		count++ ;
            	if(count > 2*NO_OF_COLUMNS)
		{
			fprintf(fp,"\n");
			count = 1;
		}

		fprintf(fp,INT_OUPTUT_FORMAT_SPEC,400) ;
	}

	fprintf(fp,"\n") ;

	/* writing Boundary no of triangles */
	count = 0 ;
	
	for(i=0; i<Tetmesh.nMarkers; i++) // Loop to count and print # of faces on each boundary
	{
		count++ ;
            	if(count > 2*NO_OF_COLUMNS) // Limit number of line prints
                {
			fprintf(fp,"\n");
			count=1;
		}
		
		sum=0;	//sum keeps track of the cumulative total number of faces printed so far
		for(j=0;j<=i;j++)
		{
                	sum = sum+Tetmesh.faceCounts[j]; // Add previous face counts to sum
		}

		fprintf(fp,INT_OUPTUT_FORMAT_SPEC,sum) ;
	}

	//Compare final sum to the read in value for the number of faces
	printf("Final Face Count:%d\nFile Face Count: %d\n",sum,Tetmesh.nFaces);
	
	fprintf(fp,"\n") ;

/* writing Boundary no of Quadrilaterals */
	count = 0 ;
	for(i=0; i<Tetmesh.nMarkers; i++)//Print zeros for # of quadrilateral faces
	{
		count++ ;
            	if(count > 2*NO_OF_COLUMNS)
                {
			fprintf(fp,"\n");
			count = 1;
		}

                fprintf(fp,INT_OUPTUT_FORMAT_SPEC,0) ;
	}
       
	fprintf(fp,"\n");


	/* writing Boundary names */
	//Note: Boundary names are arbitrary at this point if code reading this mesh uses
	//      a boundary condition file.
	for(i=0;i<Tetmesh.nMarkers;i++)
	{
		fprintf(fp,"%s\n","Slip-wall");
		
	}

//==============================================================================
//    	Print Boundary Face Connectivity
//==============================================================================

	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,sum) ;
	fprintf(fp,"\n");
        
	for(j=0; j<3; j++) //Loop over all 3 columns of faces array
        {
		count=0;

		for(k=0;k<Tetmesh.nMarkers;k++)
		{
               	
	       		for(i=0; i<Tetmesh.nFaces; i++)
          		{

				if(Tetmesh.faceMarkers[i] == Tetmesh.boundaryMarkers[k])
				{
					count++ ;
		       			if(count > 2*NO_OF_COLUMNS) //Limit # of entries per line
					{
						fprintf(fp,"\n");
						count=1;
					}

					fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.faces[i][j]) ;
				}
	   		}

		}

            	fprintf(fp,"\n") ; //Newline after printing all vertex data
	}

        
	/* writing boundary quadrilaterals *************/
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,0) ;
	fprintf(fp,"\n");

//**************************************************************************
// 	Close file
//**************************************************************************

	fclose(fp) ;
	
//**************************************************************************
//      Print tetmesh object contents
//**************************************************************************

	printf("\nPrinting Tetmesh Data Structure Contents:\n");
	printf("# of Nodes: %d\n# of Elements: %d\n",Tetmesh.nNodes,Tetmesh.nElems);
	printf("# of Faces: %d\n# of Boundaries: %d\n\n",Tetmesh.nFaces,Tetmesh.nMarkers);


  	printf("Writing tetrahedral CENTAUR grid file done.\n") ;

  	return ;
}

