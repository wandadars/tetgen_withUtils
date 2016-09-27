/**********************************************************************/
/*	Purposse: Code to take an unstructured tetrahedral mesh and   */
/*	          convert to a .surf surface mesh format for AFLR3D   */
/*	          input.                                              */
/*                                                                    */
/*	Author 	: Christopher Neal				      */
/*	Date 	: November 20, 2015                                   */
/*	Updated : March 04, 2016	        		      */
/**********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include<math.h>
#include<omp.h>

typedef struct
{
	double  	x, y, z ; // Coordinates of vertex
} vertex ;


typedef struct
{
	int     nNodes, nFaces,nFaceNodes,nMarkers;
	vertex 	*nodes ;
	vertex  *FaceNodes;
    	int     **faces;
	int     *faceMarkers ;
	int	*faceCounts;
	int 	*boundaryMarkers;
} tetmesh ;

/* GLOBAL VARIABLE DECLARATIONS */
tetmesh 	Tetmesh ;
const char *FLOAT_OUTPUT_FORMAT_SPEC = "%16.9E"; // field width is 26, 16 decimal digits, exponential notation with big E
const char *INT_OUPTUT_FORMAT_SPEC = "%-16d"; // Field width for integer number
int NUM_OF_OPENMP_THREADS=1;		//Global variable used for openmp parallelization

/* DECLARATIONS FOR FUNCTIONS */
void 	output_surf_tetmesh(void) ;
void    read_tetmesh(void) ;


void main()
{

  // Read in tetmesh format
  read_tetmesh();

  // Print tetmesh format in SURF format
  output_surf_tetmesh();
}



void read_tetmesh()
{

        // Local Variable Declarations
	int	i,j,k,marker,dummy ;
	int	nNodes,nDims,nAttrs,MarkerFlag ;
	int 	nElems,nNpT,nRegions ;
	int     nFaces, nodeShift ;
	int	v1,v2,v3,v4 ;
	int     alloc_Err ;
	int     tid,nthreads;  //OpenMP variables
	double	x,y,z ;
	double  *Vol;
	char	casename[100],filename[115] ;
	FILE	*fp ;

	printf("WELCOME TO TETMESH->SURF CONVERTER\n\n");
	printf("The Code reads TETGEN tetrahedral meshes and outputs triangular surface mesh, SURF, mesh format. \n\n") ;

	printf("Number of Threads to Use(OpenMP Parallelization)?");
	scanf("%d",&NUM_OF_OPENMP_THREADS);			
		
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


	/* Read faces ****************************************************/
  	sprintf(filename,"%s.face",casename) ;
	printf("Reading File: %s \n\n",filename) ; // Echo filename being read

      	fp = fopen(filename,"r") ;

	fscanf(fp,"%d %d",&nFaces,&MarkerFlag) ;//MarkerFlag-->boundary markers flag
	printf("Face File Data\n");
	printf("# of Faces: %d\nBoundary Marker Flag: %d \n",nFaces,MarkerFlag) ; // Echo first line back to screen

	Tetmesh.nFaces = nFaces ; // Store total number of faces in mesh
	

	// Allocate for faces array( 2D array)
	printf("Allocating memory for storing boundary face data \n");
	printf("Estimated memory usage for face array is: %10.6f\t Megabytes \n",(float)sizeof(int)*3*nFaces/(1.0e6));
	printf("Estimated memory usage for storing all face markers(not just unique) is: %10.6f\tMegabytes \n",(float)sizeof(int)*nFaces/(1.0e6));

	int **FaceData = (int **) malloc(nFaces*sizeof(int *));
	Tetmesh.faceMarkers = (int *) malloc(nFaces*sizeof(int)) ;// Allocate 1D array for all face markers
	
	Tetmesh.nMarkers = 1;	// Initialize unique boundary marker counter
	Tetmesh.boundaryMarkers = malloc(1*sizeof(int)); // Allocate first element of unique boundary marker array

	/* Allocate structures for storing unique face node data */
	int NumFaceNodes = 1;   // Initialize unique face nodes counter
        int *UniqueFaceNodes = malloc(1*sizeof(int)); // Allocate first element of unique face nodes array

	printf("Extracting Unique Boundary Markers and Unique Face Nodes from Volume Mesh \n");
	for(i=0; i<Tetmesh.nFaces; i++)
	{

	  fscanf(fp,"%d %d %d %d %d",&dummy,&v1,&v2,&v3,&marker) ;
	  //printf("%d %d %d %d %d \n",dummy,v1,v2,v3,marker) ;
	  FaceData[i] =(int *) malloc(3*sizeof(int)) ; // Allocate columns for row i
	//Tetmesh.faces[i] = (int *) malloc(3*sizeof(int)) ; // Allocate columns for row i

          //Order of read matters for Clockwise/Anti-Clockwise ordering
          FaceData[i][0] = v1 + nodeShift;
          FaceData[i][1] = v2 + nodeShift;
          FaceData[i][2] = v3 + nodeShift;

          Tetmesh.faceMarkers[i] = marker ;
	 

	  /*  Store the unique Boundary Markers and Face Nodes that are Read from the File  */    
	  if(i == 0) // set current Marker and face node to the first entries that are read from file
	  {
		Tetmesh.boundaryMarkers[i] = marker;
		UniqueFaceNodes[i] = FaceData[i][0];
       	  }

	  /* Store unique boundary Markers */
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
	

       	  /*Store unique Face nodes */
          for( j=0; j<3;j++) //Loop over all nodes for this face
          {

             int newNode = 1;   //Boolean for tracking if new node is found. 0=false, 1=true
             for(k=0;k<NumFaceNodes;k++) // Compare currently read marker with markers that have been stored
             {
                if(FaceData[i][j] == UniqueFaceNodes[k])
                {
                        newNode = 0; // We have found this node before
                }

             }

             if( newNode == 1  ) // A new unique node is found. 
             {
                NumFaceNodes = NumFaceNodes + 1;

                // Dynamically resize UniqueFaceNodes array to hold newly found node
                int *temp2 = realloc(UniqueFaceNodes,NumFaceNodes*sizeof(int));  // Store data in a larger temp array
                if(temp2 != NULL) // realloc was successful
                {
                        UniqueFaceNodes = temp2; //Point original array to newly sized array
                        UniqueFaceNodes[NumFaceNodes-1] = FaceData[i][j]; // Store newly read marker into array
                }
                else
                {
                        free(UniqueFaceNodes);
                        printf("Error allocating memory for UniqueFaceNodes Array!\n");
                        return ;
                }

             }
          }



	}

	// Display the number of unique markers in .face file
	printf("Number of Unique Boundary Markers: %d \n",Tetmesh.nMarkers);


	/*Sort the UniqueFaceNodes array into ascending order*/
	printf("Sorting Unique Face Nodes Array into Ascending Order \n");
	int a;
	for (i=0; i < NumFaceNodes; ++i)
        {
           for (j=i+1; j < NumFaceNodes; ++j)
           {
               if (UniqueFaceNodes[i] > UniqueFaceNodes[j])
               {
                  a =  UniqueFaceNodes[i];
                  UniqueFaceNodes[i] = UniqueFaceNodes[j];
                  UniqueFaceNodes[j] = a;
               }
           }
        }


	/* The code now has a list of all of the unique nodes. We now create a temporary FaceNodes array holding only the nodes for the faces*/
        /* and then update the Tetmesh.nodes array to have only the face node data */
	printf("Total Boundary Face Nodes:%d\n\n",NumFaceNodes);
	
	//for(i=0;i<NumFaceNodes;i++)
	//{
	//	 printf("Unique Boundary Node #%d \t is:%d\n\n",i+1,UniqueFaceNodes[i]);
	//}


	/*Loop over all nodes in the nodes array for the volume mesh and store only ones that are on a face*/
	vertex *FaceNodes =(vertex *) malloc(NumFaceNodes*sizeof(vertex));  // Store data in a temp array
	for(i=0;i<NumFaceNodes;i++)
	{
		FaceNodes[i].x = Tetmesh.nodes[UniqueFaceNodes[i]-1].x;
		FaceNodes[i].y = Tetmesh.nodes[UniqueFaceNodes[i]-1].y;
		FaceNodes[i].z = Tetmesh.nodes[UniqueFaceNodes[i]-1].z;
	}


	/*Update the faces array to have the re-mapped connectivity that includes only nodes on the faces */
	
	/*
	for(i=0;i<Tetmesh.nFaces;i++)
	{
		for(j=0;j<3;j++)
		{
			 printf("%d \t ",Tetmesh.faces[i][j]);
		}
		printf("\n");
	}
	*/

	printf("Re-Mapping face node numbers. This may take a while. \n");

	//Debugging
	/*
	for(i=0;i<nFaces;i++)
	{
	  for(j=0;j<3;j++)
	  {
	    printf("%d\t",FaceData[i][j]);
	  }
	  printf("\n");
	}
	*/

	omp_set_num_threads(NUM_OF_OPENMP_THREADS);
	#pragma omp parallel default(none) shared(FaceData) firstprivate(nthreads,i,j,k,nFaces,NumFaceNodes,UniqueFaceNodes,tid)
	{
		tid = omp_get_thread_num();
		if (tid == 0)
    		{
    			nthreads = omp_get_num_threads();
   			printf("Number of threads = %d\n", nthreads);
        	}

		#pragma omp for
		for(i=0;i<nFaces;i++)
		{
	   	    for(j=0;j<3;j++)
	   	    {
	      	       	for(k=0;k < NumFaceNodes;k++)
	      		{
	         	    if(FaceData[i][j] == UniqueFaceNodes[k])
		 	    {
			        FaceData[i][j] = k+1;
				break;
		 	    }
	      		}
	   	    }	
		}


	}


	//Debugging
	/*
        for(i=0;i<nFaces;i++)
        { 
          for(j=0;j<3;j++)
          {
            printf("%d\t",FaceData[i][j]);
          }
          printf("\n");
        }
	*/

	int CheckSum = 0;
	for(i=0;i<nFaces;i++)
        {   
            for(j=0;j<3;j++)
            {   
                CheckSum = CheckSum + FaceData[i][j];
            }
        }



	printf("Adding face data to Tetmesh struct variable\n");
	Tetmesh.faces = (int **) malloc(nFaces*sizeof(int *)) ; // Allocate rows of face array
	for(i=0;i<nFaces;i++)
	{
	  Tetmesh.faces[i] = (int *) malloc(3*sizeof(int)) ; // Allocate columns for row i
	  for(j=0;j<3;j++)
	  {
		Tetmesh.faces[i][j] = FaceData[i][j];
	  }
	}
	free(FaceData);



	printf("Sum of Face Array Node Data for Error Checking:%d\n",CheckSum);
	
	/*Re-allocate the nodes array and use the face nodes data only*/
	printf("Storing Face Node Data \n");
	free(Tetmesh.nodes);
	Tetmesh.nodes = (vertex *) malloc(NumFaceNodes*sizeof(vertex)) ;
	for(i=0;i<NumFaceNodes;i++)
	{
		Tetmesh.nodes[i] = FaceNodes[i];
	}

	Tetmesh.nNodes = NumFaceNodes;





	/* Allocate faceCounts array */
	Tetmesh.faceCounts = malloc(Tetmesh.nMarkers*sizeof(int));

	// Initialize faceCounts array to zero. This array will hold the number
	// of faces that have a specific boundary marker.
	for(i=0; i<Tetmesh.nMarkers; i++)
	{
		Tetmesh.faceCounts[i] = 0;
	}
	

	/*Loop through faceMarkers array and count the number of faces that have all of the markers*/
	/*that are stored in boundaryMarkers*/
	printf("Counting number of faces associated with each boundary marker \n");
	for(i=0;i<Tetmesh.nFaces; i++) // Loop through every face
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



	/*Re-Order the boundaryMarkers array to be in order of ascending boundary marker. Change*/
	/* e.g. [ 1 5 4 6] --> [1 4 5 6] */
	printf("Re-Ordering boundary markers to be in ascending order \n");
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
				int temp3 = Tetmesh.boundaryMarkers[j];	//Store smaller value
				int temp4 = Tetmesh.faceCounts[j];	
				
				//Swap boundaryMarkers
				Tetmesh.boundaryMarkers[j] = Tetmesh.boundaryMarkers[i];
				Tetmesh.boundaryMarkers[i] = temp3;

				//Swap faceCounts
				Tetmesh.faceCounts[j] = Tetmesh.faceCounts[i];
				Tetmesh.faceCounts[i] = temp4;
			}			
		}
	}

	// Print the number of faces on each boundary that were counted in file
	printf("Number of Faces On Each Boundary\n");	
	int sum = 0;
	for( i=0;i<Tetmesh.nMarkers;i++)
	{
	  	printf("Boundary Marker %d(%d) = %d \n",Tetmesh.boundaryMarkers[i],i+1,Tetmesh.faceCounts[i]);
		sum = sum + Tetmesh.faceCounts[i];  
	}

	//Error Check, Print face counts to screen
	printf("Total Faces Counted:%d\n# Faces In File:%d\n\n",sum,Tetmesh.nFaces);

	printf("Finished Reading Data File\n\n");
}


void output_surf_tetmesh()
{
	int	i, j, k ;
	int	print_tags;
	char	casename[20], filename[20] ;
	char 	TagOutputFmt[]= "%d\t%s%-d\t\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d";
	FILE	*fp ;

//*************************************************************
// 	Open file and read title
//*************************************************************

	printf("Writing ASCII triangular SURF grid file...\n\n") ;
	
	/*Store output casename provided by user*/
	printf("What is output casename?  i.e. <casename>.ugrid:");
	scanf("%s",casename);

	 /*Store choice for whether to output a .tags file*/
        printf("\n\nPrint additional <>.tags file?  0-No   1-Yes:");
        scanf("%d",&print_tags);
	               
	sprintf(filename,"%s.surf",casename) ;
	printf("\nWriting triangular SURF format to: \"%s\" \n\n",filename) ;
	fp = fopen(filename,"w") ;

                
/*Begin Printing to File*/


//===================================================================
//	Print Node Coordinates
//===================================================================

	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.nFaces);
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC, 0);
	fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.nNodes);

	fprintf(fp,"\n");

	printf("Writing Node Coordinates...\n") ;

	for(i=0; i<Tetmesh.nNodes; i++) // Print all coordinates of nodes
	{

           	fprintf(fp,FLOAT_OUTPUT_FORMAT_SPEC,Tetmesh.nodes[i].x) ;//Print x coordinate of node i
		fprintf(fp,FLOAT_OUTPUT_FORMAT_SPEC,Tetmesh.nodes[i].y) ;//Print y coordinate of node i
		fprintf(fp,FLOAT_OUTPUT_FORMAT_SPEC,Tetmesh.nodes[i].z) ;//Print z coordinate of node i
		fprintf(fp,FLOAT_OUTPUT_FORMAT_SPEC, 0.0) ;			 //Print Initial Normal Spacing of node i

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
	
		fprintf(fp,INT_OUPTUT_FORMAT_SPEC,Tetmesh.faceMarkers[i]);
		fprintf(fp,INT_OUPTUT_FORMAT_SPEC,0);
		fprintf(fp,INT_OUPTUT_FORMAT_SPEC,-1);

                fprintf(fp,"\n") ; //Newline after printing face onnectivity data
        }


/****************************************************************************/
/* 	Close SURF File 						    */
/****************************************************************************/

	fclose(fp) ;

/****************************************************************************/
/*      Write .tags filee                                                   */
/****************************************************************************/
	if(print_tags == 1)
	{
		sprintf(filename,"%s.tags",casename) ;
        	printf("Writing tags file for SURF mesh: \"%s\" \n",filename) ;
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
	printf("# of Nodes: %d\n",Tetmesh.nNodes);
	printf("# of Faces: %d\n# of Boundaries: %d\n\n",Tetmesh.nFaces,Tetmesh.nMarkers);


  	printf("Writing tetrahedral SURF grid file done.\n") ;

  	return ;
}

