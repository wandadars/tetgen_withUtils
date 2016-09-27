/**********************************************************************/
/*	Purposse: Code to convert unstructured tetrahedral mesh to    */
/*	          CENTAUR ascii format.                               */
/*	                                                              */
/*	Author 	: Christopher Neal				      */
/*	Date 	: July 29th, 2015                                     */
/*	Updated : August 8th, 2015	        		      */
/**********************************************************************/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>


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
	double  *ElemVol;
} tetmesh ;

/* GLOBAL VARIABLE DECLARATIONS */
tetmesh 	Tetmesh ;
const char *FLOAT_OUTPUT_FORMAT_SPEC = "%16.9E"; // field width is 26, 16 decimal digits, exponential notation with big E
const char *INT_OUPTUT_FORMAT_SPEC = "%-16d"; // Field width for integer number

/* DECLARATIONS FOR FUNCTIONS */
void 	output_ugrid_tetmesh(void) ;
void    read_tetmesh(void) ;


void main()
{
  // Read in tetmesh format
  read_tetmesh();

  // Print tetmesh format in UGRID format
  output_ugrid_tetmesh();
}



void read_tetmesh()
{

        // Local Variable Declarations
	int	i,j,marker,dummy ;
	int	nNodes,nDims,nAttrs,MarkerFlag ;
	int 	nElems,nNpT,nRegions ;
	int     nFaces, nodeShift ;
	int	v1,v2,v3,v4 ;
	int     alloc_Err ;
	double	x,y,z ;
	double  *Vol;
	char	casename[100],filename[115] ;
	FILE	*fp ;

	printf("WELCOME TO TETMESH->UGRID CONVERTER\n\n");
	printf("The Code reads TETGEN tetrahedral meshes and outputs UGRID mesh format. \n\n") ;

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
        nodeShift = 0;	//Used for if nodes are numbered from 0(code doesn't like this)
	for(i=0; i<nNodes; i++)
	{
	
		if(MarkerFlag == 1 )
		{
	  		fscanf(fp,"%d %lf %lf %lf %d",&dummy,&x,&y,&z,&dummy) ;
			if( i == 0 && dummy == 0)
			{
				nodeShift = 1;
				printf("Node shifted to start at 1\n");		
			}
	  		//printf("%lf %lf %lf\n",x,y,z) ; //Echo data from file
          		Tetmesh.nodes[i].x = x ;
          		Tetmesh.nodes[i].y = y ;
          		Tetmesh.nodes[i].z = z ;
	  	}
	  	else
	  	{
			fscanf(fp,"%d %lf %lf %lf",&dummy,&x,&y,&z) ;
			if(i == 0 && dummy == 0)
                        {
                                nodeShift = 1;
				printf("Node shifted to start at 1\n");
                        }
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

	  //Usually 1 2 3 4 
          Tetmesh.elems[i][0] = v1 + nodeShift;
          Tetmesh.elems[i][1] = v2 + nodeShift;
          Tetmesh.elems[i][2] = v3 + nodeShift;
          Tetmesh.elems[i][3] = v4 + nodeShift;
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
          Tetmesh.faces[i][0] = v1 + nodeShift;
          Tetmesh.faces[i][1] = v2 + nodeShift;
          Tetmesh.faces[i][2] = v3 + nodeShift;

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
	

	/*Loop through faceMarkers array and count the number of faces that have all of the markers*/
	/*that are stored in boundaryMarkers*/
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



	/*Compute the volume of all elements to see what the minimum volume is*/
	/*Allocate array to hold element areas*/
	Tetmesh.ElemVol = malloc(Tetmesh.nElems*sizeof(double));
	Vol = malloc(Tetmesh.nElems*sizeof(double));

	double MinVol,MaxVol,xP,yP,zP; //For Debugging
        for(i=0;i<Tetmesh.nElems;i++) // Loop over all elements
        {

                /*Initialize centroid location*/
                xP = 0.0 ;
                yP = 0.0 ;
                zP = 0.0 ;

                double x[4], y[4], z[4]; // Array to hold node coordinates

                for(j=0;j<4;j++) // Loop over nodes stores in element array
                {
                  /*Compute X ,Y,Z coordinates of element centroid*/
                  xP = xP +(1.0/4.0)*( Tetmesh.nodes[ Tetmesh.elems[i][j]-1 ].x );
                  yP = yP +(1.0/4.0)*( Tetmesh.nodes[ Tetmesh.elems[i][j]-1 ].y );
                  zP = zP +(1.0/4.0)*( Tetmesh.nodes[ Tetmesh.elems[i][j]-1 ].z );

                  x[j] =  Tetmesh.nodes[ Tetmesh.elems[i][j]-1 ].x;
                  y[j] =  Tetmesh.nodes[ Tetmesh.elems[i][j]-1 ].y;
                  z[j] =  Tetmesh.nodes[ Tetmesh.elems[i][j]-1 ].z;

                }

		/*Compute Element volume using the following Formula*/
                /* V = (1/3)Aface*h ,Aface=area of a face,h=perpendicular distance from face to 4th point */
                double a,b,c,d,s,Aface;
                a = sqrt( (x[0]-x[1])*(x[0]-x[1]) + (y[0]-y[1])*(y[0]-y[1]) + (z[0]-z[1])*(z[0]-z[1]) );
                b = sqrt( (x[0]-x[2])*(x[0]-x[2]) + (y[0]-y[2])*(y[0]-y[2]) + (z[0]-z[2])*(z[0]-z[2]) );
                c = sqrt( (x[1]-x[2])*(x[1]-x[2]) + (y[1]-y[2])*(y[1]-y[2]) + (z[1]-z[2])*(z[1]-z[2]) );
                s = 0.5*(a+b+c);

                Aface = sqrt( s*(s-a)*(s-b)*(s-c) ); //Area of a face from Heron's Formula


                /*Create vector that is perpendicular to the plane formed by the 3 points above for Aface*/
                double A, B, C, D; // Define normal vector components for plane
                A = (y[1]-y[0])*(z[2]-z[0])-(y[2]-y[0])*(z[1]-z[0]);
                B = -1*( (x[1]-x[0])*(z[2]-z[0])-(x[2]-x[0])*(z[1]-z[0])  );
                C = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
                D = A*x[0] + B*y[0] + C*z[0];

                double h; //Perpendicular distance from plane to point
                h = fabs( A*x[3] + B*y[3] + C*z[3] - D )/sqrt(A*A + B*B + C*C);

		/*
		printf("Element(%d), h = %10.2E \t Aface= %10.2E \n",i+1,h,Aface);		
		*/

                Vol[i] = (1.0/3.0)*Aface*h; // Store element volume

		Tetmesh.ElemVol[i] = Vol[i]; // Store element volume;


                /*Compare computed values to find max and min of r and Vol*/
                if( i == 1 ) // first iteration assume values of max and min are same
                {
                        MinVol = Vol[i];
                        MaxVol = Vol[i];
                }
		else // update max and min values
                {
                        if(Vol[i] < MinVol)
                        {
                                MinVol = Vol[i];
                        }

                        if(Vol[i] > MaxVol)
                        {
                                MaxVol = Vol[i];
                        }

                }

        }

	/*Print Out Min and Max element volume*/
        printf("Maximum element volume: %10.2E \n",MaxVol);
        printf("Minimum element volume: %10.2E \n",MinVol);



	/*Re-Order the boundaryMarkers array to be in order of ascending boundary marker. Change*/
	/* e.g. [ 1 5 4 6] --> [1 4 5 6] */
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


void output_ugrid_tetmesh()
{
	int	i, j, k ;
	int	print_tags;
	char	casename[20], filename[20] ;
	char 	TagOutputFmt[]= "%d\t%s%-d\t\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d";
	FILE	*fp ;

//*************************************************************
// 	Open file and read title
//*************************************************************

	printf("Writing ASCII tetrahedral UGRID grid file...\n\n") ;
	
	/*Store output casename provided by user*/
	printf("What is output casename?  i.e. <casename>.ugrid:");
	scanf("%s",casename);

	 /*Store choice for whether to output a .tags file*/
        printf("\n\nPrint additional <>.tags file?  0-No   1-Yes:");
        scanf("%d",&print_tags);
	               
	sprintf(filename,"%s.ugrid",casename) ;
	printf("\nWriting tetrahedral UGRID format to: \"%s\" \n\n",filename) ;
	fp = fopen(filename,"w") ;

	/*Check if the user wants binary or ascii data printed to file*/
	//printf("Write to ASCII or BINARY? : 0-ASCII , 1-BINARY  ");
	//scanf("%d",&bin_or_ascii);
	

                
/*Begin Printing to File*/


//===================================================================
//	Print Node Coordinates
//===================================================================

	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.nNodes);
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.nFaces);
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC, 0);
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.nElems);
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,0);
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,0);
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,0);

	fprintf(fp,"\n");

	printf("Writing Node Coordinates...\n") ;

	for(i=0; i<Tetmesh.nNodes; i++) // Print all coordinates of nodes
	{

           	fprintf(fp,FLOAT_OUTPUT_FORMAT_SPEC,Tetmesh.nodes[i].x) ;//Print x coordinate of node i
		fprintf(fp,FLOAT_OUTPUT_FORMAT_SPEC,Tetmesh.nodes[i].y) ;//Print y coordinate of node i
		fprintf(fp,FLOAT_OUTPUT_FORMAT_SPEC,Tetmesh.nodes[i].z) ;//Print z coordinate of node i
		
		fprintf(fp,"\n"); //Move down one line after printing x,y,z coordinates of node            

        }



/*==============================================================================*/
/*      Print Boundary Face Connectivity                                        */
/*==============================================================================*/
        for(i=0; i<Tetmesh.nFaces; i++) //Loop over all rows of faces array
        {
                
        	for(j=0; j<3; j++)
                {
                        fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.faces[i][j]) ;
                }

                fprintf(fp,"\n") ; //Newline after printing facec onnectivity data
        }

/*==============================================================================*/
/*      Print  Face Boundary Markers		                                */
/*==============================================================================*/
	for(i=0; i<Tetmesh.nFaces; i++) //Loop over all elements of face marker array
	{
		fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.faceMarkers[i]);
		fprintf(fp,"\n");
	}

//=====================================================================
//	Cell connectivity
//=====================================================================

	/* writing tetrahedral connectivities *********/

    	printf("Writing Tetrahedral Cell Connectivity...\n") ;
        
	for(i=0; i<Tetmesh.nElems; i++)
        {
                for(j=0; j<4; j++) // Print nodes that make up elements
                {

                   	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.elems[i][j]) ;
                }

		fprintf(fp,"\n") ;
        }


/****************************************************************************/
/* 	Close UGRID File						    */
/****************************************************************************/

	fclose(fp) ;

/****************************************************************************/
/*      Write .tags filee                                                   */
/****************************************************************************/
	if(print_tags == 1)
	{
		sprintf(filename,"%s.tags",casename) ;
        	printf("Writing tags file for UGRID mesh: \"%s\" \n",filename) ;
	        fp = fopen(filename,"w") ;

		fprintf(fp,"# Generated by SolidMesh version 5.10.0\n");
		fprintf(fp,"#\t\t\t\t\t\t\t\t\t\tTrans\tInitial\tB-L\t# of\n") ;
		fprintf(fp,"#ID\tGroup\t\t\tVisc\tRecon\tRebuild\tFixed\tSource\tTrans\tDelete\tSpacing\tThcknss\tLayers") ;
		fprintf(fp,"\n") ;


		for(i=0;i<Tetmesh.nMarkers;i++)
		{
			char tempName[]="BC_";
			fprintf(fp,TagOutputFmt,Tetmesh.boundaryMarkers[i],"BC_",Tetmesh.boundaryMarkers[i],0,0,0,0,0,0,0,0,0,0);
			fprintf(fp,"\n");

		}


		fclose(fp);

	}




	
/****************************************************************************/
/*      Print tetmesh object contents					    */
/****************************************************************************/

	printf("\nPrinting Tetmesh Data Structure Contents:\n");
	printf("# of Nodes: %d\n# of Elements: %d\n",Tetmesh.nNodes,Tetmesh.nElems);
	printf("# of Faces: %d\n# of Boundaries: %d\n\n",Tetmesh.nFaces,Tetmesh.nMarkers);


  	printf("Writing tetrahedral UGRID grid file done.\n") ;

  	return ;
}

