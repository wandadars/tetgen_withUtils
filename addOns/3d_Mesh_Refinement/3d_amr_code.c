/**********************************************************************/
/*	Purpose: Code to read an unstructured 3D tetrahedral mesh and */
/*	         and create a .vol file detailing specific volumes for*/
/*	         each tetrahedral element for mesh refinement.        */
/*	         	                        		      */
/*	         	                                              */
/*	 Notes: Code reads 2 files: .node and .ele to produce .vol    */
/*	 	file. This code is specifically designed to work only */
/*	 	with the 3D meshes that Tetgen is used to create      */
/*	                                                              */
/*	Author 	: Christopher Neal				      */
/*	Date 	: September 22, 2014 	        		      */
/*                                                                    */
/*      Edited: October, 24th, 2014                                   */
/*                                                                    */
/*                                                                    */
/**********************************************************************/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct
{
	double  	x, y, z; // Coordinates of vertex
} vertex ;

typedef struct
{
	int     nNodes, nElems, nNodeMarkers;
	vertex 	*nodes ;
	int	*nodeNums;
	int     *nodeMarkers;
	int     *nodeBoundaryMarkers;
   	int     **elems ;
	double  *ElemVol;
} tetmesh ;


/* GLOBAL VARIABLE DECLARATIONS */
tetmesh 	Tetmesh ;
char *FLOAT_OUTPUT_FORMAT_SPEC = "%23.16E"; // field width is 26, 16 decimal digits, exponential notation with big E
char *INT_OUPTUT_FORMAT_SPEC = "%8d"; // Field width is 8, number is integer
char    casename[25];

/* DECLARATIONS FOR FUNCTIONS */
void 	output_tetrahedral_volume_file(void) ;
void    read_tetmesh(void) ;
void	cramer_rule(double [3][4], double [3]);
double  tet_vol_calc(double[3][4]);

void main()
{
	
   printf("\n\n----------WELCOME TO TEGEN ADAPTIVE MESH REFINEMENT PROGRAM----------\n\n");
   printf("The Code reads in data from a 3D tetrahedral mesh and produces an adaptive mesh\n") ;
   printf("Run this code in the directory where the Tetgen mesh files are located.\n") ;

  // Read in tetmesh format
  read_tetmesh();
  
  // Print hybrid mesh format in CENTAUR format
  output_tetrahedral_volume_file();
}




void read_tetmesh()
{

        // Local Variable Declarations
	int	i,j,marker,dummy ;
	int	nNodes,nDims,nAttrs,MarkerFlag ;
	int 	nElems,nNpT,nAttributes ;
	int	debug = 0; // Debugging flag
	int	v1,v2,v3,v4 ;
	double	x,y,z ;
	char	filename[30] ;
	FILE	*fp ;

	printf("----------STARTING TETGEN MESH READER----------\n\n");
	printf("The Code reads tetrahedral meshes created by the program Tetgen. \n\n") ;

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

	//Allocate for node numbering array(1D array)
	Tetmesh.nodeNums = malloc(nNodes*sizeof(int));

	//Allocate for nodeMarkers array(1D array)
	Tetmesh.nodeMarkers = malloc(nNodes*sizeof(int));
        
	// Loop through file and read all node data
	for(i=0; i<nNodes; i++)
	{
	
	  if(MarkerFlag == 1 ) // additional marker line is present
	  {
	  	fscanf(fp,"%d %lf %lf %lf %d",&dummy,&x,&y,&z,&dummy) ;
	  	//printf("%lf %lf \n",x,y) ; //Echo data from file
          	Tetmesh.nodes[i].x = x ;
          	Tetmesh.nodes[i].y = y ;
		Tetmesh.nodes[i].z = z ;
		Tetmesh.nodeMarkers[i] = dummy;
		Tetmesh.nodeNums[i] = i+1;
	  }
	  else
	  {
		fscanf(fp,"%d %lf %lf %lf",&dummy,&x,&y,&z) ;
                //printf("%lf %lf \n",x,y) ; //Echo data from file
                Tetmesh.nodes[i].x = x ;
                Tetmesh.nodes[i].y = y ;
		Tetmesh.nodes[i].z = z ;
		Tetmesh.nodeNums[i] = i+1;

	  }


	}

	fclose(fp) ;

	/* Read elements *************************************************/
  	sprintf(filename,"%s.ele",casename) ;
	printf("Reading File: %s \n\n",filename) ; // Echo filename to be read

      	fp = fopen(filename,"r") ; // Open .ele file for reading elements

	fscanf(fp,"%d %d %d",&nElems,&nNpT,&nAttributes) ; // Read first line
	printf("Element File Data\n");
	printf("# of Elements: %d\n# of Verts Per Element:  %d\n# of Attributes: %d \n\n",nElems,nNpT,nAttributes) ; //Echo first line to user

	Tetmesh.nElems = nElems ; // Store number of elements into Tetmesh variable

	// Allocate for elements array ( 2D array)
	Tetmesh.elems = (int **) malloc(nElems*sizeof(int *)) ; // Allocate rows of array
        
        //Loop through and store element vertices
	for(i=0; i<nElems; i++)
	{
	  fscanf(fp,"%d %d %d %d %d",&dummy,&v1,&v2,&v3,&v4) ;
	  //printf("%d %d %d %d \n",dummy,v1,v2,v3,v4) ;
	  Tetmesh.elems[i] = (int *) malloc(nNpT*sizeof(int)) ; // Allocate columns of array

          //Check this section if grid ouput is not correct(adjust clockwise/counter-clockwise read)
	  Tetmesh.elems[i][0] = v2 ; 
          Tetmesh.elems[i][1] = v1 ;
          Tetmesh.elems[i][2] = v3 ;
	  Tetmesh.elems[i][3] = v4 ;
	}

	fclose(fp) ;

	printf("Finished Reading Tetgen Mesh .node and .ele Data Files\n\n");
	
	//DEBUGGING PURPOSES
	if(debug == 1)//print entire contents of tetmesh data structure
        {

                printf("\n\nPrinting Contents of Tetmesh Data Structure");
                printf("\nNumber of Nodes: %d",Tetmesh.nNodes);
                printf("\nNumber of Elements: %d",Tetmesh.nElems);


                printf("\n\n------------Printing Node Data------------");
                printf("\nNode # \t X Coord \t Y Coord \t Z Coord");
                for(i=0;i<Tetmesh.nNodes;i++)
                {
                        printf("\n %d %4.2e %4.2e %4.2e",Tetmesh.nodeNums[i],Tetmesh.nodes[i].x,Tetmesh.nodes[i].y,Tetmesh.nodes[i].z);

                }


        }


}





void output_tetrahedral_volume_file(void)
{
	int	i, j,k;
	int 	AMR_Flag;
	double  xc, yc, zc, xP, yP,zP, r, R0, V0, Vol_Percent,Vmin_Factor,Vmax_inner;
	double  Vmin, Vmax, rmin,rmax;
	double  *Vol, *radial_pos;//area and radial location of elements
	char	filename[30];
	FILE	*fp ;

//*************************************************************
// 	Open file
//*************************************************************

	//Read Necessary Parameters from Input File
	printf("\n%s\n","Reading Data from input files: amr_inputs.dat");
	printf("%s\n","Reading X,Y,Z coordinates of refinement Center and the AMR Flag");
        fp = fopen("amr_inputs.dat","r") ;
	
	fscanf(fp,"%lf %lf %lf %d",&xc,&yc,&zc,&AMR_Flag);
	
	// Read additional Parameters
	if(AMR_Flag == 1)//linear element scaling, fixed bounding volumes
	{
		printf("\n%s\n","Linear Element Scaling with Minimum Volume Computed from Mesh Selected");
	}
	else if(AMR_Flag == 2)//quadratic element scaling, fixed bounding volumes 
	{
		printf("\n%s\n","Quadratic Element Scaling with Fixed Minimum Volume Selected");

		printf("\n%s\n%s\n%s","Reading:"," Refinement Radius","Volume Percent Threshold");
		fscanf(fp,"%lf %lf",&R0,&Vol_Percent);
	}
	else if(AMR_Flag == 3)// Prescribed Minimum Volume With Linear Variation
	{
		printf("\n%s\n","Linear Element Scaling with Minimum Volume Scaling Selected");

		printf("\n%s\n%s\n%s","Reading:","Minimum Volume Factor(0 to 1)");
		fscanf(fp,"%lf",&Vmin_Factor);
	}
	else if(AMR_Flag == 4)//Qaudratic element scaling, prescribed minimum Volume
	{
		printf("\n%s\n","Quadratic Element Scaling with Minimum Volume Scaling Selected");

		printf("\n%s\n%s\n%s","Reading:"," Refinement Radius"," Volume Percent Threshold");
                fscanf(fp,"%lf %lf",&R0,&Vol_Percent);

		printf("\n%s\n"," Minimum Volume Factor(0 to 1)");
                fscanf(fp,"%lf",&Vmin_Factor);
	}
	else if(AMR_Flag == 5)//Quadratic element scaling, constant inner region element size
	{
	
		printf("\n%s\n","Quadratic Element Scaling with Uniform Mesh inside Refinement Radius Selected");

                printf("\n%s \n%s","Reading:"," Refinement Radius");
                fscanf(fp,"%lf",&R0);

                printf("\n%s\n"," Maximum Inner Volume ");
                fscanf(fp,"%lf",&Vmax_inner);

	}

	//Echo User Inputs Back
	printf("\n%s\n","User Input Is:");
	printf("%s %lf\n","X-Center Coord: ",xc);
	printf("%s %lf\n","Y-Center Coord: ",yc);
	printf("%s %lf\n","Z-Center Coord: ",zc);
	printf("%s %d\n","Element Scaling Flag: ",AMR_Flag);

	if(AMR_Flag == 2)
	{
	   printf("%s %lf\n","Refinement Radius: ",R0);
 	   printf("%s %lf\n","Percent Volume Inside Refinement Radius: ",Vol_Percent);
	}
	else if(AMR_Flag == 3)
	{
	   printf("%s %lf\n","Minimum Volume Factor(0 to 1): ",Vmin_Factor);
	}
	else if(AMR_Flag == 4)
	{
	   printf("%s %lf\n","Refinement Radius: ",R0);
           printf("%s %lf\n","Percent Volume Inside Refinement Radius: ",Vol_Percent);
	   printf("%s %lf\n","Minimum Volume Factor(0 to 1): ",Vmin_Factor);
	}
	else if(AMR_Flag == 5)
	{
	   printf("%s %lf\n","Refinement Radius: ",R0);
           printf("%s %lf\n","Maximum Volume Constraint Inside Refinement Radius: ",Vmax_inner);

	}

	fclose(fp); // Close input file
	
	printf("\n\nComputing Tetrahedral area constraints for .vol file...\n") ;

	sprintf(filename,"%s.vol",casename) ;

	printf("Data will be written to: \"%s\" \n",filename) ;
	fp = fopen(filename,"w") ;


//===================================================================
//	Compute Element Areas and Print to File
//===================================================================

	//Allocate array to hold element areas
	Tetmesh.ElemVol = malloc(Tetmesh.nElems*sizeof(double));

	Vol = malloc(Tetmesh.nElems*sizeof(double));
	radial_pos = malloc(Tetmesh.nElems*sizeof(double));

	printf("%s\n","Allocation of Arrays Complete");

	double MinVol, Minh, MinLoc, MaxVol, Maxh, MaxLoc; //For Debugging
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
		
		Vol[i] = (1.0/3.0)*Aface*h; // Store element volume
		

		/************************DEBUGGING**************************************/
			
//		if( Vol[i] > 1E6 )
//                {

//			printf("\n\n%s \n%s %lf \n%s %lf \n%s %lf \n%s %lf","Plane Parameters..."," A: ",A," B: ",B,"C: ",C,"D: ",D);
			
//			printf("\n\n%s \n%s %lf \n%s %lf \n%s %lf ","4th Point Coordinates..."," X: ",x[3]," Y: ",y[3],"Z: ",z[3]);
			
			/*Print Values of h*/
//                       printf("\n\n%s %10.3E %s %10.3E","Face Area:  ",Aface," Perpendicular Distance: ", h);

                        /*Print Element Volumes*/
//                        printf("\n%s %d %s %10.3E","Element ",i+1," Volume",Vol[i]);

//                }


		/**********************************************************************/
		Tetmesh.ElemVol[i] = Vol[i]; // Store element volume;

		/*Compute element distance from center provided by user*/
		radial_pos[i] = sqrt( (xc-xP)*(xc-xP)  + (yc-yP)*(yc-yP) + (zc-zP)*(zc-zP) );
	
		/*Compare computed values to find max and min of r and Vol*/
		if( i == 1 ) // first iteration assume values of max and min are same
		{
			Vmin = Vol[i];
			Vmax = Vol[i];
			rmin = radial_pos[i];
			rmax = radial_pos[i];
		}
		else // update max and min values
		{
			if(Vol[i] < Vmin)
			{
				Vmin = Vol[i];
			}

			if(Vol[i] > Vmax)
			{
				Vmax = Vol[i];
			}
			
			if(radial_pos[i] < rmin)
			{
				rmin = radial_pos[i];
			}
			
			if(radial_pos[i] > rmax)
			{
				rmax = radial_pos[i];
			}

		}

	}


	//Define fitting constants
	double C1,C2,C3,C4,C5,C6;

	if(AMR_Flag == 1) //Linear Variation
	{
	   
	   //Echo curve fit parameters
           printf("\n%s %4.2E \n%s %4.2E \n","Mininum Volume: ",Vmin,"Maximum Volume: ",Vmax);
           printf("%s %f \n%s %f\n","Minimum Radius: ",rmin,"Maximum Radius: ", rmax);
           
           //Linear Variation
           C1 = (Vmax-Vmin)/(rmax-rmin);
           C2 = Vmin - C1*rmin;
           
           printf("\n%s\n","Coefficients to linear volume sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n","A(r) = C1r +C2","C1 = ",C1,"C2 = ",C2);

	   printf("\n\n%s","Predicted Volumes Using Sizing Function(Check to see if it matches input)");
           printf("\n%s %4.2E","Max Volume: ",C1*rmax + C2);
           printf("\n%s %4.2E\n","Min Volume: ",C1*rmin+C2);

                                  
	}
	else if(AMR_Flag == 2) // Quadratic variation
	{

	   /*Compute intermediate volume*/
	   V0 = Vmin + (Vmax-Vmin)*(Vol_Percent/100);
	          		
	   /*Echo curve fit parameters*/
	   printf("\n%s %4.2E \n%s %4.2E \n%s %4.2E \n","Mininum Volume: ",Vmin,"Critical Volume: ",V0,"Maximum Volume: ",Vmax);
	   printf("%s %f \n%s %f\n","Minimum Radius: ",rmin,"Maximum Radius: ", rmax);
	   

	   /*Quadratic variation with coefficients solved by CRAMER'S RULE*/
	   printf("\n\n%s\n","Solving the following Equations to Scale Element Volumes");
	   printf("%s\n%s\n%s\n%s\n","OUTER SOLUTION","2*C1*rmax + C2 = 0","C1*rmax^2 + C2*rmax + C3 = Vmax","C1*R0^2 + C2*R0 + C3 = V0");

	   /*Quadratic variation with coefficients solved by CRAMER'S RULE*/
	   printf("\n%s\n%s\n%s\n%s\n","INNER SOLUTION","2*C1*rmin + C2 = 0","C1*rmin^2 + C2*rmin + C3 = Vmin","C1*R0^2 + C2*R0 + C3 = V0");
	  

           double A[3][4];
	   double x[3];
		
	   /*Construct the augmented matrix for outer solution coefficients*/
	   A[0][0] = 2*rmax;
	   A[0][1] = 1;
	   A[0][2] = 0;
	   A[0][3] = 0;

	   A[1][0] = rmax*rmax;
           A[1][1] = rmax;
           A[1][2] = 1;
           A[1][3] = Vmax;

	   A[2][0] = R0*R0;
	   A[2][1] = R0;
	   A[2][2] = 1;
	   A[2][3] = V0;

	   /*Call CRAMER'S RULE solver to get coefficients*/
	   cramer_rule(A,x); // x has the coefficients
	   C1 = x[0];
	   C2 = x[1];
	   C3 = x[2];
	   

	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
           printf("\n\n%s\n","Coefficients for OUTER SOLUTION...");
           for(i=0;i<3;i++) // Print Matrix
           {
                printf("[%4.2E]C1 + [%4.2E]C2 + [%4.2E]C3 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
           }



	   /*Construct Augmented Matrix for inner solution coefficients*/

	   A[0][0] = 2*rmin;
           A[0][1] = 1;
           A[0][2] = 0;
           A[0][3] = 0;

           A[1][0] = rmin*rmin;
           A[1][1] = rmin;
           A[1][2] = 1;
           A[1][3] = Vmin;

           A[2][0] = R0*R0;
           A[2][1] = R0;
           A[2][2] = 1;
           A[2][3] = V0;	   


	   /*Call CRAMER'S RULE solver to get coefficients*/
	   cramer_rule(A,x); // x has the coefficients
 	   
	   C4 = x[0];
	   C5 = x[1];
	   C6 = x[2];

	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
	   printf("\n\n%s\n","Coefficients for INNER SOLUTION...");
	   for(i=0;i<3;i++) // Print Matrix
	   { 
		printf("[%4.2E]C4 + [%4.2E]C5 + [%4.2E]C6 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
	   }
		

	   printf("\n%s\n","Coefficients to the OUTER quadratic volume sizing function ");
 	   printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","V(r) = C1r^2 +C2r+C3","C1 = ",C1,"C2 = ",C2,"C3 = ",C3);

           printf("\n\n%s\n","Coefficients to the INNER quadratic volume sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","V(r) = C4r^2 +C5r+C6","C4 = ",C4,"C5 = ",C5,"C6 = ",C6);

	
	   printf("\n\n%s","Predicted Volumes Using Sizing Function(Check to see if it matches input)");
	   printf("\n%s","OUTER volume sizing function predictions");
	   printf("\n%s %4.2E","Max Volume: ",C1*rmax*rmax+C2*rmax+C3);
	   printf("\n%s %4.2E\n","Critical Volume: ",C1*R0*R0+C2*R0+C3);

           printf("\n\n%s","INNER volume sizing function predictions");
           printf("\n%s %4.2E","Min Volume: ",C4*rmin*rmin+C5*rmin+C6);
           printf("\n%s %4.2E\n","Critical Volume: ",C4*R0*R0+C5*R0+C6);


	}
	else if(AMR_Flag == 3) // Prescribed Minimum Volume With Linear Variation
	{
	
	   /*Update minimum volume using user-defined factor*/
	   Vmin = Vmin*(1-Vmin_Factor/100);
	
	   /*Echo curve fit parameters*/
	   printf("\n%s %4.2E \n%s %4.2E \n","Mininum Volume: ",Vmin,"Maximum Volume: ",Vmax);
           printf("%s %f \n%s %f\n","Minimum Radius: ",rmin,"Maximum Radius: ", rmax);

	   /*Linear Variation*/
	   C1 = (Vmax-Vmin)/(rmax-rmin);
           C2 = Vmin - C1*rmin;

           printf("\n%s\n","Coefficients to linear volume sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n","V(r) = C1r +C2","C1 = ",C1,"C2 = ",C2);

		

	}
	else if(AMR_Flag == 4)
	{
	

	   Vmin = Vmin*Vmin_Factor; // Update Minimum Volume

	   /*Compute intermediate area*/
	   V0 = Vmin + (Vmax-Vmin)*(Vol_Percent/100);
	   
	   /*Echo curve fit parameters*/
	   printf("\n%s %4.2E \n%s %4.2E \n%s %4.2E \n","Mininum Volume: ",Vmin,"Critical Volume: ",V0,"Maximum Volume: ",Vmax);
	   printf("%s %f \n%s %f\n","Minimum Radius: ",rmin,"Maximum Radius: ", rmax);
	   
	   /*Quadratic variation with coefficients solved by CRAMER'S RULE*/
	   printf("\n\n%s\n","Solving the following Equations to Scale Element Volume");
	   printf("%s\n%s\n%s\n%s\n","OUTER SOLUTION","2*C1*rmax + C2 = 0","C1*rmax^2 + C2*rmax + C3 = Vmax","C1*R0^2 + C2*R0 + C3 = V0");
	   
	   /*Quadratic variation with coefficients solved by CRAMER'S RULE*/
	   printf("\n%s\n%s\n%s\n%s\n","INNER SOLUTION","2*C1*rmin + C2 = 0","C1*rmin^2 + C2*rmin + C3 = Vmin","C1*R0^2 + C2*R0 + C3 = V0");
	   
	   
	   double A[3][4];
	   double x[3];
	   
	   /*Construct the augmented matrix for outer solution coefficients*/
	   A[0][0] = 2*rmax;
	   A[0][1] = 1;
	   A[0][2] = 0;
	   A[0][3] = 0;
	   
	   A[1][0] = rmax*rmax;
	   A[1][1] = rmax;
	   A[1][2] = 1;
	   A[1][3] = Vmax;
	   
	   A[2][0] = R0*R0;
	   A[2][1] = R0;
	   A[2][2] = 1;
	   A[2][3] = V0;
	   
	   /*Call CRAMER'S RULE solver to get coefficients*/
	   cramer_rule(A,x); // x has the coefficients
	   C1 = x[0];
	   C2 = x[1];
	   C3 = x[2];                                                                                                                                                                                    

  	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
           printf("\n\n%s\n","Coefficients for OUTER SOLUTION...");
           for(i=0;i<3;i++) // Print Matrix
           {
                printf("[%4.2E]C1 + [%4.2E]C2 + [%4.2E]C3 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
           }



           /*Construct Augmented Matrix for inner solution coefficients*/

           A[0][0] = 2*rmin;
           A[0][1] = 1;
           A[0][2] = 0;
           A[0][3] = 0;

           A[1][0] = rmin*rmin;
           A[1][1] = rmin;
           A[1][2] = 1;
           A[1][3] = Vmin;

           A[2][0] = R0*R0;
           A[2][1] = R0;
           A[2][2] = 1;
           A[2][3] = V0;


           /*Call CRAMER'S RULE solver to get coefficients*/
           cramer_rule(A,x); // x has the coefficients

           C4 = x[0];
           C5 = x[1];
           C6 = x[2];

	   
	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
           printf("\n\n%s\n","Coefficients for INNER SOLUTION...");
           for(i=0;i<3;i++) // Print Matrix
           {
                printf("[%4.2E]C4 + [%4.2E]C5 + [%4.2E]C6 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
           }


           printf("\n%s\n","Coefficients to the OUTER quadratic volume sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","V(r) = C1r^2 +C2r+C3","C1 = ",C1,"C2 = ",C2,"C3 = ",C3);

           printf("\n\n%s\n","Coefficients to the INNER quadratic volume sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","V(r) = C4r^2 +C5r+C6","C4 = ",C4,"C5 = ",C5,"C6 = ",C6);


           printf("\n\n%s","Predicted Volumes Using Sizing Function(Check to see if it matches input)");
           printf("\n%s","OUTER area sizing function predictions");
           printf("\n%s %4.2E","Max Volume: ",C1*rmax*rmax+C2*rmax+C3);
           printf("\n%s %4.2E\n","Critical Volume: ",C1*R0*R0+C2*R0+C3);

           printf("\n\n%s","INNER volume sizing function predictions");
           printf("\n%s %4.2E","Min Volume: ",C4*rmin*rmin+C5*rmin+C6);
           printf("\n%s %4.2E\n","Critical Volume: ",C4*R0*R0+C5*R0+C6);

 
	}
	else if(AMR_Flag ==5) // Fixed inner region cell volume.
	{
	
	   Vmin = Vmax_inner; // Update Minimum Volume

           /*Intermediate volume is dictacted by minimum volume now*/
           V0 = Vmin;

           /*Echo curve fit parameters*/
           printf("\n%s %4.2E \n%s %4.2E \n%s %4.2E \n","Mininum Volume: ",Vmin,"Critical Volume: ",V0,"Maximum Volume: ",Vmax);
           printf("%s %f \n%s %f\n","Minimum Radius: ",rmin,"Maximum Radius: ", rmax);

           /*Quadratic variation with coefficients solved by CRAMER'S RULE*/
           printf("\n\n%s\n","Solving the following Equations to Scale Element Volume");
           printf("%s\n%s\n%s\n%s\n","OUTER SOLUTION","2*C1*rmax + C2 = 0","C1*rmax^2 + C2*rmax + C3 = Vmax","C1*R0^2 + C2*R0 + C3 = A0");

           /*Constant Variation Governed By User Input*/
           printf("\n%s\n%s","INNER SOLUTION","V=Vmin (constant)");


           double A[3][4];
           double x[3];

           /*Construct the augmented matrix for outer solution coefficients*/
           A[0][0] = 2*rmax;
           A[0][1] = 1;
           A[0][2] = 0;
           A[0][3] = 0;

           A[1][0] = rmax*rmax;
           A[1][1] = rmax;
           A[1][2] = 1;
           A[1][3] = Vmax;

           A[2][0] = R0*R0;
           A[2][1] = R0;
           A[2][2] = 1;
           A[2][3] = V0;

           /*Call CRAMER'S RULE solver to get coefficients*/
           cramer_rule(A,x); // x has the coefficients
           C1 = x[0];
           C2 = x[1];
           C3 = x[2];   
	 
	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
           printf("\n\n%s\n","Coefficients for OUTER SOLUTION...");
           for(i=0;i<3;i++) // Print Matrix
           {
                printf("[%4.2E]C1 + [%4.2E]C2 + [%4.2E]C3 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
           }


	   /*Construct Augmented Matrix for inner solution coefficients*/

           A[0][0] = 2*rmin;
           A[0][1] = 1;
           A[0][2] = 0;
           A[0][3] = 0;

           A[1][0] = rmin*rmin;
           A[1][1] = rmin;
           A[1][2] = 1;
           A[1][3] = Vmin;

           A[2][0] = R0*R0;
           A[2][1] = R0;
           A[2][2] = 1;
           A[2][3] = Vmin;


           /*Call CRAMER'S RULE solver to get coefficients*/
           cramer_rule(A,x); // x has the coefficients

           C4 = x[0];
           C5 = x[1];
           C6 = x[2];

	   
	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
           printf("\n\n%s\n","Coefficients for INNER SOLUTION...");
           for(i=0;i<3;i++) // Print Matrix
           {
                printf("[%4.2E]C4 + [%4.2E]C5 + [%4.2E]C6 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
           }



           printf("\n%s\n","Coefficients to the OUTER quadratic volume sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","V(r) = C1r^2 +C2r+C3","C1 = ",C1,"C2 = ",C2,"C3 = ",C3);

           printf("\n\n%s\n","Coefficients to the INNER quadratic volume sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","V(r) = C4r^2 +C5r+C6","C4 = ",C4,"C5 = ",C5,"C6 = ",C6);


           printf("\n\n%s","Predicted Volumes Using Sizing Function(Check to see if it matches input)");
           printf("\n%s","OUTER volume sizing function predictions");
           printf("\n%s %4.2E","Max Volume: ",C1*rmax*rmax+C2*rmax+C3);
           printf("\n%s %4.2E\n","Critical Volume: ",C1*R0*R0+C2*R0+C3);

           printf("\n\n%s","INNER volume sizing function predictions");
           printf("\n%s %4.2E","Min Volume: ",C4*rmin*rmin+C5*rmin+C6);
           printf("\n%s %4.2E\n","Critical Volume: ",C4*R0*R0+C5*R0+C6);



	}



	/*************UPDATE ELEMENT VOLUMES***************************************/
	printf("\n%s","Creating List of Updated Element Volumes to .vol file...");
	for(i = 0; i<Tetmesh.nElems;i++) // Scale all volumes between min and max according to r
	{

		if(AMR_Flag == 1) // Linear+Unaltered Min Volume
		{
		  Vol[i] = C1*radial_pos[i] + C2;
		}
		else if (AMR_Flag == 2)//Quadratic+Unaltered Min Volume
		{
		
		  if(radial_pos[i] <= R0) //Scale Inner
		  {
		    Vol[i] = C4*radial_pos[i]*radial_pos[i] + C5*radial_pos[i] + C6;
		  }
		  else // Scale Outer
		  {
		    Vol[i] = C1*radial_pos[i]*radial_pos[i] + C2*radial_pos[i] + C3;
		  }

		}
		else if (AMR_Flag == 3)//Linear+Altered Min Vol
		{
		   Vol[i] = C1*radial_pos[i] + C2;
		
		}
		else if(AMR_Flag == 4) //Quadratic+Altered Min Vol
		{

             	   if(radial_pos[i] <= R0) // Scale Inner
                   {
                     Vol[i] = C4*radial_pos[i]*radial_pos[i] + C5*radial_pos[i] + C6;
                   }
                   else // Scale Outer
                   {
                     Vol[i] = C1*radial_pos[i]*radial_pos[i] + C2*radial_pos[i] + C3;
                   }

		}
		else if(AMR_Flag == 5) // Quadratic+Uniform inner
		{
		
		  if(radial_pos[i] <= R0) // Scale Inner
                   {
                     Vol[i] = C4*radial_pos[i]*radial_pos[i] + C5*radial_pos[i] + C6;
                   }
                   else // Scale Outer
                   {
                     Vol[i] = C1*radial_pos[i]*radial_pos[i] + C2*radial_pos[i] + C3;
                   }


		}

	}

	fprintf(fp,"%8d\n",Tetmesh.nElems);
	
	printf("%s\n","Printing to File Now");

	//Print output to file
	for(i=0;i<Tetmesh.nElems;i++)
	{
	   fprintf(fp,"%8d %23.16E",i+1,Vol[i]);
	   fprintf(fp,"\n");
	}

//**************************************************************************
// 	Close file
//**************************************************************************

	fclose(fp) ;

//**************************************************************************
//      Print trimesh object contents
//**************************************************************************

  	printf("Writing Adaptive Meshed Tetgen .vol file done.\n") ;

  	return ;
}


void    cramer_rule(double A[3][4], double x[3])
{

	/* Local variable declarations*/
	double D,D1,D2,D3;
	double first, second, third;;

	// The Ds are determinants, and first, second, and third are the
	// three parts of the 3x3 determinant

	first=    A[0][0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2]);
	second =  -1*A[0][1]*(A[1][0]*A[2][2]-A[2][0]*A[1][2]);
	third =   A[0][2]*(A[1][0]*A[2][1]-A[2][0]*A[1][1]);
	
	D = first + second + third;


	first=    A[0][3]*(A[1][1]*A[2][2]-A[2][1]*A[1][2]);
        second =  -1*A[0][1]*(A[1][3]*A[2][2]-A[2][3]*A[1][2]);
        third =   A[0][2]*(A[1][3]*A[2][1]-A[2][3]*A[1][1]);

        D1 = first + second + third;


	first=    A[0][0]*(A[1][3]*A[2][2]-A[2][3]*A[1][2]);
        second =  -1*A[0][3]*(A[1][0]*A[2][2]-A[2][0]*A[1][2]);
        third =   A[0][2]*(A[1][0]*A[2][3]-A[2][0]*A[1][3]);

        D2 = first + second + third;


	first=    A[0][0]*(A[1][1]*A[2][3]-A[2][1]*A[1][3]);
        second =  -1*A[0][1]*(A[1][0]*A[2][3]-A[2][0]*A[1][3]);
        third =   A[0][3]*(A[1][0]*A[2][1]-A[2][0]*A[1][1]);

        D3 = first + second + third;

	
	//Compute solutions using Cramer's Rule
	x[0] = D1/D;
	x[1] = D2/D;
	x[2] = D3/D;

}
