
/* 	RNA Secondary Structure Prediction, Using the Algorithm of Zuker
		C++ version by David Mathews, copyright 1996, 1997, 1998, 1999

		Programmed for Isis Pharmaceuticals and the Turner Lab,
      	Chemistry Department, University of Rochester
*/

#include "stdafx.h"
#include <stdio.h>
# include <string.h>
# include <iostream.h>
# include <fstream.h>
              
# include <math.h>
#include <stdlib.h>

#include "platform.cpp"
#include "algorithm.h"
//platform.cpp contains information specific to the
			//	machine




#define maxfil 100    //maximum length of file names
#define maxtloop 100 //maximum tetraloops allowed (info read from tloop)
#define maxstructures 1010 //maximum number of structures in ct file
#define maxbases 3000   //maximum number of bases in a structure
#define ctheaderlength 125 //maximum length of string containing info on sequence
#define ga_bonus -10 //the value of a bonus for the "almost coaxial stacking"
								//case in efn2
#define amax 400 //this is a maximum line length for void linout (below)
#define col 80  //this is the number of columns in an output file
#define numlen 8  //maximum digits in a number
#define maxnopair 600 //maximum number of bases that can be forced single

#include "structure.cpp"//uncommented out for sgis
#include "algorithm.h"


#define debugmode false //a flag to turn on/off debugging features


//********************************functions:

/*	Function efn2

	The energy calculator of Zuker
   calculates the free energy of each structural conformation in a structure

   structures cannot have pseudoknots

   structnum indicates which structure to calculate the free energy
   	the default, 0, indicates "all" structures

*/




/*  Function opendat

		Function opens data files to read thermodynamic data
*/

int opendat (char *loop2,char *stackf,char *tstackh,char *tstacki,
		char *tloop,char *miscloop, char *danglef, char *int22, char *int21,
      char *coax,char *tstackcoax,char *coaxstack,
      char *tstack,char *tstackm, char *triloop, char *int11, datatable* data)
{
char lineoftext[100],base[110];
int count,i,j,k,l, m, a, b, c, d,e,f,g;
float temp;
FILE *check;

//eparam[1] is a basepair bonus
//eparam[2] is a bulge loop bonus
//eparam[3] is an interior loop bonus

 ifstream ml1;
 ifstream lo1;
 ifstream st1;
 ifstream th1;
 ifstream ti1;
 ifstream tl1;
 ifstream da1;
 ifstream in1;
 ifstream in2;
 ifstream tri;
 ifstream co1;
 ifstream co2;
 ifstream co3;
 ifstream st2;
 ifstream tsm;
 ifstream i11;

 //ofstream out("check.out");


//check that all the files exist with a C i/o function
if ((check = fopen(miscloop, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(loop2, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(stackf, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tstackh, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tstacki, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tloop, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(danglef, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(int22, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(int21, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(triloop, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(coax, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tstackcoax, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(coaxstack, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tstack, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tstackm, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(int11,"r")) == NULL) return 0;

fclose(check);

//open the files using the C++ method for reading
 ml1.open(miscloop);
 lo1.open(loop2);
 st1.open(stackf);
 th1.open(tstackh);
 ti1.open(tstacki);
 tl1.open(tloop);
 da1.open(danglef);
 in1.open(int22);
 in2.open(int21);
 tri.open(triloop);
 co1.open(coax);
 co2.open(tstackcoax);
 co3.open(coaxstack);
 st2.open(tstack);
 tsm.open(tstackm);
 i11.open(int11);

/* Read information from miscloop */

// the key sequence "-->" now indicates a record


ml1>>lineoftext;
 while(strcmp(lineoftext,"-->")) ml1>>lineoftext;


 ml1 >> (data->prelog);

 data->prelog = (data->prelog)*10;


 ml1>>lineoftext;
 while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 ml1 >> temp;
 data->maxpen = int (temp*10 + .5);


 ml1>>lineoftext;
 while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 for (count=1;count<= 4;count ++)
 { ml1 >> temp;
 (data->poppen[count])= (int) (temp*10 + .5);



 } 										//this reads float values, converts
											// 	them int and assigns them into
											//		array poppen
 ml1>>lineoftext;
 while(strcmp(lineoftext,"-->")) ml1>>lineoftext;
 

 data->eparam[1] = 0;						 // assign some variables that are
 data->eparam[2] = 0;						 //	"hard-wired" into code
 data->eparam[3] = 0;
 data->eparam[4] = 0;
 ml1 >> temp;
 data->eparam[5] = short (floor (temp*10+.5));  //constant multi-loop penalty



 ml1 >> temp;
 data->eparam[6] = short (floor (temp*10+.5));



 data->eparam[7] = 30;
 data->eparam[8] = 30;
 data->eparam[9] = -500;
 ml1 >> temp;
 data->eparam[10] = short (floor (temp*10+.5));
 ml1>>lineoftext;
 while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 ml1 >> temp;

 if (ml1.peek()==EOF) {
  	//these are old energy rules -- treat the other constants properly
   data->efn2a = data->eparam[5];
   data->efn2b = data->eparam[6];
   data->efn2c = data->eparam[10];
   data->auend=0;
   data->gubonus=0;
   data->cslope = 0;
   data->cint=0;
   data->c3=0;
   data->init=0;
   data->gail = 0;


 }

 else {
 	data->efn2a = short (floor (temp*10+.5));  //constant multi-loop penalty for efn2

	 ml1 >> temp;
 	data->efn2b= short (floor(temp*10+.5));

 	ml1 >> temp;
 	data->efn2c= short (floor(temp*10+.5));

 	//now read the terminal AU penalty:
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;
   ml1>> temp;

 	data->auend = short (floor (temp*10+.5));

 	//now read the GGG hairpin bonus:
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 	ml1 >> temp;
 	data->gubonus = short (floor (temp*10+.5));

	//now read the poly c hairpin penalty slope:
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 	ml1 >> temp;
 	data->cslope = short (floor (temp*10+.5));

 	//now read the poly c hairpin penalty intercept:
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 	ml1 >> temp;
 	data->cint = short (floor (temp*10+.5));

 	//now read the poly c penalty for a loop of 3:
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->"))ml1>>lineoftext;

 	ml1 >> temp;
 	data->c3 = short (floor (temp*10+.5));
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;
 	
   ml1 >> temp;
 	data->init = short (floor (temp*10+.5));

 	//now read the GAIL rule indicator
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 	ml1>> temp;
 	data->gail = short (floor (temp+.5));
 }


 /*	read info from dangle */
 //add to dangle the case where X (represented as 0) is looked up
for (l = 1;l <=2; l++){
	for (i = 0;i <=5; i++){
		if ((i!=0)&&(i!=5)) for (count=1;count <=60;count++) da1 >> lineoftext;
		for (j=0;j<=5; j++) {
			for (k=0;k<=5; k++) {
				if ((i==0)||(j==0)||(k==0)) {
				    data->dangle[i][j][k][l] = 0;
				}
            else if ((i==5)||(j==5)||(k==5)) {
             	data->dangle[i][j][k][l] = 0;
            }
				else {
				    da1 >> lineoftext;
                //cout << lineoftext<<"\n";
				    if (strcmp(lineoftext,".")){
					data->dangle[i][j][k][l] = short (floor (10*(atof(lineoftext))+.5));
				    }
				    else data->dangle[i][j][k][l] = infinity;
				}
            //cout <<"dangle = "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<data->dangle[i][j][k][l]<<"\n";
			}
         //cin >> m;
		}
	}
}


/*	read info from loop for internal loops, hairpin loops, and bulge loops */

for (count = 1; count <=26; count++) lo1 >> lineoftext; //get past text in file
for (i=1;i <= 30; i++) {
	lo1 >> lineoftext;//get past the size column in table
	lo1 >> lineoftext;
					if (strcmp(lineoftext,".")){
					data->inter[i] = short (floor (10*(atof(lineoftext))+.5));
					}
					else data->inter[i] = infinity;

               //cout <<"inter = "<<data->inter[i]<<"\n";
	lo1 >> lineoftext;
					if (strcmp(lineoftext,"."))
					data->bulge[i] = short (floor(10*(atof(lineoftext))+.5));
					else data->bulge[i] = infinity;

               //cout <<"bulge = "<<data->bulge[i]<<"\n";
	lo1 >> lineoftext;
					if (strcmp(lineoftext,".")){
					data->hairpin[i] = short (floor(10*(atof(lineoftext))+.5));
					}
					else data->hairpin[i] = infinity;

               //cout <<"hair = "<<data->hairpin[i]<<"\n";
}


/* Read info from stack */
//add to the stack table the case where X (represented as 0) is looked up:

for (count=1;count<=42;count++) st1 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) st1 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {
				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->stack[i][j][k][l]=0;
				}
            else if ((i==5)||(j==5)||(k==5)||(l==5)) {
             	data->stack[i][j][k][l] = infinity;
            }
				else {
					st1 >> lineoftext;
					if (strcmp(lineoftext,".")){
						data->stack[i][j][k][l] =short (floor(10*(atof(lineoftext))+.5));
					}
					else data->stack[i][j][k][l] = infinity;
				}

			}

		}
	}
}





/* Read info from tstackh */
//add to the tstackh table the case where X (represented as 0) is looked up:

for (count=1;count<=46;count++) th1 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) th1 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {
				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->tstkh[i][j][k][l]=0;
				}
            else if ((i==5)||(j==5)||(k==5)||(l==5)) {
             	data->tstkh[i][j][k][l] = infinity;
            }
				else {
				    th1 >> lineoftext;
				    if (strcmp(lineoftext,".")){
					    data->tstkh[i][j][k][l] =short (floor(10*(atof(lineoftext))+.5));
				    }
				    else data->tstkh[i][j][k][l] = infinity;
				}
            //cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->tstkh[i][j][k][l]<<"\n";
			}
         //cin >>m;
		}
	}
}

/* Read info from tstacki */
//add to the tstacki table the case where X (represented as 0) is looked up:

for (count=1;count<=46;count++) ti1 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) ti1 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {

				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->tstki[i][j][k][l]=0;
				}
            else if ((i==5)||(j==5)||(k==5)||(l==5)) {
             	data->tstki[i][j][k][l] = infinity;
            }
				//else if ((k==5)||(l==5)) {
				//include "5", linker for intermolecular for case of flush ends
				//	data->tstki[i][j][k][l]=0;
				//}
				else {
				    ti1 >> lineoftext;
                //cout <<lineoftext<<"\n";

				    if (strcmp(lineoftext,".")){
					data->tstki[i][j][k][l] =short (floor (10*(atof(lineoftext))+.5));
				    }
				    else data->tstki[i][j][k][l] = infinity;
				}
            //cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->tstki[i][j][k][l]<<"\n";
			}
         //cin >> m;
		}
	}
}

/*	Read info from tloops */
for (count=1;count<=3;count++)	tl1 >> lineoftext;//get past text in file
data->numoftloops=0;
tl1>>lineoftext;


for (count=1;count<=maxtloop&&!tl1.eof();count++){
	//cout << lineoftext;

	(data->numoftloops)++;
	strcpy(base,lineoftext);
	strcpy(base+1,"\0");
	data->tloop[data->numoftloops][0] = tonumi(base);
//cout << base << "\n";
//cout << data->tloop[data->numoftloops][0] << "\n";

	strcpy(base,lineoftext+1);
	strcpy(base+1,"\0");
	data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
		5*tonumi(base);
//cout << base << "\n";
//cout << data->tloop[data->numoftloops][0] << "\n";
	strcpy(base,lineoftext+2);
	strcpy(base+1,"\0");
	data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
		25*tonumi(base);
//cout << base << "\n";
//cout << data->tloop[data->numoftloops][0] << "\n";
	strcpy(base,lineoftext+3);
	strcpy(base+1,"\0");
	data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
		125*tonumi(base);
//cout << base << "\n";
//cout << data->tloop[data->numoftloops][0] << "\n";
	strcpy(base,lineoftext+4);
	strcpy(base+1,"\0");
	data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
		625*tonumi(base);
//cout << base << "\n";
//cout << data->tloop[data->numoftloops][0] << "\n";
	strcpy(base,lineoftext+5);
	strcpy(base+1,"\0");
	data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
		3125*tonumi(base);
//cout << base << "\n";
//cout << data->tloop[data->numoftloops][0] << "\n";
	tl1 >> temp;
	data->tloop[data->numoftloops][1] = short (floor (10*temp+0.5));


	//cout << "key = "<<data->tloop[data->numoftloops][0]<<"\n";
	//cout << "bonus = "<<data->tloop[data->numoftloops][1]<<"\n";
//	cin >> j;

   tl1 >> lineoftext;
}


//Read the 2x2 internal loops
//key iloop22[a][b][c][d][j][l][k][m] =
//a j l b
//c k m d


for (count=1;count<=340;count++) in1 >> lineoftext;//get past text in file

for (i=1;i<=36;i++) {//read each of 21 tables
	for (j=1;j<=39;j++) in1 >> lineoftext;//get past text in file
	strcpy(base,lineoftext);
	strcpy(base+1, "\0");
	a = tonumi(base);
	for (j=1;j<=3;j++) in1 >> lineoftext;
	strcpy(base, lineoftext);
	strcpy(base+1, "\0");
	b = tonumi(base);
	in1>>lineoftext;
	strcpy(base, lineoftext);
	strcpy(base+1, "\0");
	c = tonumi(base);
	for (j=1;j<=3;j++) in1 >> lineoftext;
	strcpy(base, lineoftext);
	strcpy(base+1, "\0");
	d = tonumi(base);
	for (j=1;j<=3;j++) in1 >> lineoftext;//get past text in file
	for (j=1;j<=4;j++) {
	    for (k=1;k<=4;k++) {
		for (l=1;l<=4;l++) {
		    for (m=1;m<=4;m++) {
			in1 >> temp;
			data->iloop22[a][b][c][d][j][l][k][m] = short (floor(10*temp+0.5));

         //no longer need to store the reverse order at same time because
         //the tables contain redundancy:
			//data->iloop22[d][c][b][a][m][k][l][j] = floor(10*temp+0.5);


			//cout << "a = "<<a<<" b= "<<b<<" c= "<<c<<" d = "<<d<<"\n";

			//cout << "w = "<<j<<" x= "<<l<<" y= "<<k<<" z= "<<m<<"\n";
			//cout << data->iloop22[a][b][c][d][j][l][k][m]<<"\n";

		    }
          //cin >> foo;
		}
	    }
	}
}

//Read the 2x1 internal loop data
for (i=1;i<=58;i++) in2 >> lineoftext; //get past text at top of file
for (i=1;i<=6;i++) { //read each row of tables
	for (e=1;e<=4;e++) {
		for (j=1;j<=66;j++) in2 >> lineoftext; //get past text in file
		in2 >> lineoftext;
		strcpy(base,lineoftext);
		strcpy(base+1,"\0");
		a = tonumi(base);
		for (j=1;j<=11;j++) in2 >> lineoftext; //get past text in file
		in2 >> lineoftext;
		strcpy(base,lineoftext);
		strcpy(base+1,"\0");
		b = tonumi(base);
		for (j=1;j<=35;j++) in2 >> lineoftext; //get past text in file
		for (c=1;c<=4;c++) {
			for (j=1;j<=6;j++) {
				switch (j) {
					case 1:
						f = 1;
						g = 4;
						break;
					case 2:
						f = 2;
						g = 3;
						break;
					case 3:
						f = 3;
						g = 2;
						break;
					case 4:
						f = 4;
						g = 1;
						break;
					case 5:
						f = 3;
						g = 4;
						break;
					case 6:
						f = 4;
						g = 3;
						break;
				}
				for (d=1;d<=4;d++) {
					in2 >> temp;
					data->iloop21[a][b][c][d][e][f][g]=short (floor(10*temp+0.5));
					//cout << a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f<<" "<<g<<"\n";
               //cout << temp<<"\n";
               //cout << "10*temp = "<<10*temp<<"\n";
					//cout << data->iloop21[a][b][c][d][e][f][g]<<"\n";
					//cin >> temp;
				}
            //cin >> temp;
			}
		}
	}

}
/*	Read info from triloops */
for (count=1;count<=3;count++)	tri >> lineoftext;//get past text in file
data->numoftriloops=0;
tri>>lineoftext;


for (count=1;count<=maxtloop&&!tri.eof();count++){

	//cout << lineoftext<<"\n";

	(data->numoftriloops)++;
	strcpy(base,lineoftext);
	strcpy(base+1,"\0");
	data->triloop[data->numoftriloops][0] = tonumi(base);
	strcpy(base,lineoftext+1);
	strcpy(base+1,"\0");
	data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
		5*tonumi(base);
	strcpy(base,lineoftext+2);
	strcpy(base+1,"\0");
	data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
		25*tonumi(base);
	strcpy(base,lineoftext+3);
	strcpy(base+1,"\0");
	data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
		125*tonumi(base);
	strcpy(base,lineoftext+4);
	strcpy(base+1,"\0");
	data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
		625*tonumi(base);
	tri >> temp;
	data->triloop[data->numoftriloops][1] = short (floor (10*temp+0.5));

   //cout << data->triloop[data->numoftriloops][1]<< "  "<<data->triloop[data->numoftriloops][0]<<"\n";

   tri >> lineoftext;
}


/* Read info from coax */
//add to the stack table the case where X (represented as 0) is looked up:

//this is the array that keeps track of flush coaxial stacking (no intervening nucs)

//data arrangement of coax: data->coax[a][b][c][d]
//5'b-c3'
//3'a d5'
//this means the helix backbone is continuous between nucs b and c


for (count=1;count<=42;count++) co1 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) co1 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {
				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->coax[j][i][k][l]=0;
				}
            else if ((i==5)||(j==5)||(k==5)||(l==5)) {
             	data->coax[j][i][k][l] = infinity;
            }
				else {
					co1 >> lineoftext;

               //cout << lineoftext <<"end\n";

					if (strcmp(lineoftext,".")){
						data->coax[j][i][k][l] =short (floor(10*(atof(lineoftext))+.5));
					}
					else data->coax[j][i][k][l] = infinity;

               //cin >> a;
				}
            //cout << j << " "<<i<<" "<<k<<" "<<l<<"  "<<data->coax[j][i][k][l]<<"\n";
			}
		}
      //cin>>foo;
	}
}

/* Read info from tstackcoax */
//add to the tstackh table the case where X (represented as 0) is looked up:

//data arrangement of tstackcoax:
//5'a-c -> strand continues into stack
//3'b-d -> strand does not continue to stack
//pair between a-b, c-d is a mismatch


for (count=1;count<=46;count++) co2 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if (!(i==0||i==5)) for (count=1;count<=60;count++) co2 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {
				if ((i==0)||(j==0)||(k==0)||(l==0)||(i==5)||(j==5)||(k==5)||(l==5)) {
					data->tstackcoax[i][j][k][l]=0;
				}
				else {
				    co2 >> lineoftext;
				    if (strcmp(lineoftext,".")){
					    data->tstackcoax[i][j][k][l] =short (floor(10*(atof(lineoftext))+.5));
				    }
				    else data->tstackcoax[i][j][k][l] = infinity;
				}
			}
		}
	}
}
/* Read info from coaxstack */
//add to the tstackh table the case where X (represented as 0) is looked up:

//data arrangement of coaxstack:
//5'a-c ->strand contnues into stack
//3'b d ->strand does not continue to stack
//pair between a-b, mismatch between c-d
//backbone is discontinuous between b and d


for (count=1;count<=46;count++) co3 >> lineoftext;//get past text in file
for (i=0;i<=4;i++) {
	if (!(i==0||i==5)) for (count=1;count<=60;count++) co3 >> lineoftext;
	for (k=0;k<=4;k++) {
		for (j=0;j<=4;j++) {
			for (l=0;l<=4;l++) {
				if ((i==0)||(j==0)||(k==0)||(l==0)||(i==5)||(j==5)||(k==5)||(l==5)) {
					data->coaxstack[i][j][k][l]=0;
				}
				else {
				    co3 >> lineoftext;
				    if (strcmp(lineoftext,".")){
					    data->coaxstack[i][j][k][l] =short (floor(10*(atof(lineoftext))+.5));
				    }
				    else data->coaxstack[i][j][k][l] = infinity;
				}
			}
		}
	}
}



/* Read info from tstack */
//this is the terminal mismatch data used in intermolecular folding
//add to the tstack table the case where X (represented as 0) is looked up.
//also add the case where 5 (the intermolecular linker) is looked up,
//this is actually a dangling end, not a terminal mismatch.

for (count=1;count<=46;count++) st2 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) st2 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {

				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->tstack[i][j][k][l]=0;
				}
            else if ((i==5)||(j==5)) {
             	data->tstack[i][j][k][l] = infinity;

            }
				else if ((k==5)||(l==5)) {
				//include "5", linker for intermolecular for case of flush ends
            	if ((k==5)&&(l==5)) {//flush end
						data->tstack[i][j][k][l]=0;
               }
               else if (k==5) {//5' dangling end
               	//look up number for dangling end
               	data->tstack[i][j][k][l] = data->dangle[i][j][l][2]+penalty2(i,j,data);
               }
               else if (l==5) {//3' dangling end
               	data->tstack[i][j][k][l] = data->dangle[i][j][k][1]+penalty2(i,j,data);
               }
				}
				else {
				    st2 >> lineoftext;
				    if (strcmp(lineoftext,".")){
					data->tstack[i][j][k][l] =short (floor (10*(atof(lineoftext))+.5));
				    }
				    else data->tstack[i][j][k][l] = infinity;
				}
            //cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->tstki[i][j][k][l]<<"\n";
			}
         //cin >> m;
		}
	}
}


/* Read info from tstackm */
//add to the tstackm table the case where X (represented as 0) is looked up:

for (count=1;count<=46;count++) tsm >> lineoftext;//get past text in file
for (i=0;i<=4;i++) {
	if (i!=0) for (count=1;count<=60;count++) tsm >> lineoftext;
	for (k=0;k<=4;k++) {
		for (j=0;j<=4;j++) {
			for (l=0;l<=4;l++) {
				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->tstkm[i][j][k][l]=0;
				}
				else {
				    tsm >> lineoftext;
				    if (strcmp(lineoftext,".")){
					    data->tstkm[i][j][k][l] =short (floor(10*(atof(lineoftext))+.5));
				    }
				    else data->tstkm[i][j][k][l] = infinity;
				}
            //cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->tstkh[i][j][k][l]<<"\n";
			}
         //cin >>m;
		}
	}
}

//data arrangement for 1x1 loop tables iloop11[a][b][c][d][e][f]:
//abc
//def


//Read the 1x1 internal loop data
//encode the data like:  abc
//                       def where b-e is a mismatch
for (i=1;i<=58;i++) i11 >> lineoftext; //get past text at top of file
for (i=1;i<=6;i++) { //read each row of table
	if (i==1) {
    	a = 1;
      d = 4;
   }
   else if (i==2) {
    	a = 2;
      d = 3;
   }
   else if (i==3) {
    	a = 3;
      d = 2;
   }
   else if (i==4) {
    	a = 4;
      d = 1;
   }
   else if (i==5) {
    	a = 3;
      d = 4;
   }
   else {
   	a = 4;
      d = 3;
   }
	for (j=1;j<=114;j++) i11 >> lineoftext;//get past text
   for (b=1;b<=4;b++) {
   	for (j=1;j<=6;j++) {
      	if (j==1) {
    			c = 1;
      		f = 4;
   		}
   		else if (j==2) {
    			c = 2;
      		f = 3;
   		}
   		else if (j==3) {
    			c = 3;
      		f = 2;
   		}
   		else if (j==4) {
    			c = 4;
      		f = 1;
   		}
   		else if (j==5) {
    			c = 3;
      		f = 4;
   		}
   		else {
   			c = 4;
      		f = 3;
   		}
   		for (e=1;e<=4;e++) {
         	i11 >> temp;
            data->iloop11[a][b][c][d][e][f]=short (floor(10*temp+0.5));

         }

      }
   }
}

return 1;


}



//read a ct file with sequence and structural information
#define linelength 20



void openct(structure *ct,char *ctfile) {
int count,i,j;
char base[2],header[ctheaderlength],line[linelength],temp[50];
ifstream in;
in.open(ctfile);
in >> count;
j = 0;


if (count == -100) { //this is a CCT formatted file:
	in >> ct->numofbases;
   ct->allocate(ct->numofbases);
   in >> ct->numofstructures;
   //in >> ct->ctlabel[1];
   for (count=1;count<=ct->numofstructures;count++) {
   	strcpy(ct->ctlabel[count],"\n");
   }
   for (i=1;i<=ct->numofbases;i++) {
   	in >> ct->numseq[i];
      ct->nucs[i]=*tobase(ct->numseq[i]);
      ct->hnumber[i] = i;
   }
   for (count=1;count<=ct->numofstructures;count++) {
    	for (i=1;i<=ct->numofbases;i++) {
       	in >> ct->basepr[count][i];
      }
   }
   return;

}
else {//this is a ct file:
//in >> ct->numofbases;
ct->allocate(count);
in.close();
in.open(ctfile);
for (ct->numofstructures = 1;((ct->numofstructures)<=(maxstructures))
		;(ct->numofstructures)++)	{

      strcpy (header,"");

	 if ((ct->numofstructures)==(maxstructures))
			errmsg (4,0);
	 in >> ct->numofbases;
	 strcpy(line,"");
	in.getline(header,ctheaderlength);

//	 do {
//	 	strcat(header,line);
//		if (in.eof()) {
//			ct->numofstructures--;
//			return;
//		}
//		in >> line;
//		strcat(header," ");
//	 }
//	 while(strcmp(line,"1"));
//	 in.putback(*line);*/

	if(in.eof()) {
		ct->numofstructures--;
		return;
	}

	strcpy((ct->ctlabel[ct->numofstructures]),header);
	 for (count=1;count<=((ct->numofbases));count++)	{

		//if(in.eof()) {
		//	ct->numofstructures--;
		//	return;
		//}
		in >> temp;//ignore base number in ctfile
		in >> base;//read the base
		strcpy(base+1,"\0");
      ct->nucs[count]=base[0];
		tonum(base,ct,count);//convert base to numeric
      if (ct->numseq[count]==5) {
      	ct->intermolecular = true;
       	ct->inter[j] = count;
         j++;
      }
		in >>	temp;//ignore numbering
		in >> temp;//ignore numbering
		in >> ct->basepr[ct->numofstructures][count];//read base pairing info
		in >> ct->hnumber[count];//read historical numbering
	 }
}
}
(ct->numofstructures)--;
return;
}

void tonum(char *base,structure *ct,int count)	{
if (!strcmp(base,"A")) (ct->numseq[count] = 1);
else if(!strcmp(base,"B")) {
	(ct->numseq[count] = 1);

}
else if(!strcmp(base,"a")) {
	ct->numseq[count]=1;
	ct->nnopair++;
	ct->nopair[ct->nnopair] = count;
}
else if(!strcmp(base,"C")) (ct->numseq[count] = 2);
else if(!strcmp(base,"Z")) {
	(ct->numseq[count] = 2);

}
else if(!strcmp(base,"c")) {
	ct->numseq[count] = 2;
	ct->nnopair++;
	ct->nopair[ct->nnopair] = count;
}
else if(!strcmp(base,"G")) (ct->numseq[count] = 3);
else if(!strcmp(base,"H")) {
	(ct->numseq[count] = 3);

}
else if(!strcmp(base,"g")) {
	ct->numseq[count] = 3;
	ct->nnopair++;
	ct->nopair[ct->nnopair] = count;
}

else if(!strcmp(base,"U")||!strcmp(base,"T")) (ct->numseq[count] = 4);
else if(!strcmp(base,"V")||!strcmp(base,"W")) {
	(ct->numseq[count] = 4);

}
else if(!strcmp(base,"u")||!strcmp(base,"t")) {
	ct->numseq[count] = 4;
	ct->nnopair++;
	ct->nopair[ct->nnopair] = count;
}

else if(!strcmp(base,"I")) {
	ct->numseq[count] = 5;
	ct->intermolecular= true;
}

else (ct->numseq[count]=0);  //this is for others, like X
return;
}

short int tonumi(char *base)	{
short int	a;
if (!strcmp(base,"A")||!strcmp(base,"B")) (a = 1);
else if(!strcmp(base,"C")||!strcmp(base,"Z")) (a = 2);
else if(!strcmp(base,"G")||!strcmp(base,"H")) (a = 3);
else if(!strcmp(base,"U")||!strcmp(base,"V")) (a = 4);
else if(!strcmp(base,"T")||!strcmp(base,"W")) (a = 4);
else (a=0);  //this is for others, like X
return a;
}

void push(stackstruct *stack,int a,int b,int c,int d)
{
(stack->sp)++;
stack->stk[stack->sp][0]= a;
stack->stk[stack->sp][1]= b;
stack->stk[stack->sp][2]= c;
stack->stk[stack->sp][3]= d;

}


void pull(stackstruct *stack,int *i,int *j,int *open,int *null,int *stz)
{
if (stack->sp==0) {
	*stz = 1;
	return;
}
else {
*stz = 0;
*i = stack->stk[stack->sp][0];
*j = stack->stk[stack->sp][1];
*open= stack->stk[stack->sp][2];
*null= stack->stk[stack->sp][3];
stack->sp--;

}
}








//this function calculates whether a terminal pair i,j requires the end penalty
inline int penalty(int i,int j,structure* ct, datatable *data) {
	
   if (ct->numseq[i]==4||ct->numseq[j]==4)
   	return data->auend;
	else return 0;//no end penalty
	


}

//this function calculates whether a terminal pair i,j requires the end penalty
int penalty2(int i,int j, datatable *data) {


   if (i==4||j==4)
   	return data->auend;
   else return 0;//no end penalty


}




#define seqline 300 //maximum line length in .seq file
//returns 0 on error
int openseq (structure *ct,char *seqfile) {
char temp[seqline],seq[seqline],base[seqline],test[seqline];
int i,j,length,nucs;

ct->nnopair = 0;

//nucs = 0;

FILE *se;


//read the sequence file to get the number of nucleotides
se=fopen(seqfile,"r");

do {
	fgets(temp,seqline,se);	
	strncpy(test,temp,1);
	strcpy(test+1,"\0");
} while (!strcmp(test,";"));

//fgets(ct->ctlabel[1],seqline,se);

strcpy(ct->ctlabel[1],temp);

nucs = 1;
while (1) {
	fgets(seq,seqline,se);
	length = strlen (seq);
	for (j=0;j<length;j++) {//up to length-1 bcs of /n character
		strncpy (base,seq+j,1);
		strcpy (base+1,"\0");
		if (!strcmp(base,"1")) break;
      if (!(!strcmp(base,"A")||!strcmp(base,"a")||!strcmp(base,"C")||!strcmp(base,"c")
      	||!strcmp(base,"G")||!strcmp(base,"g")||!strcmp(base,"T")
         ||!strcmp(base,"t")||!strcmp(base,"U")||!strcmp(base,"u")
         ||!strcmp(base,"X")||!strcmp(base,"x")||!strcmp(base," ")||
         !strcmp(base,"\n")||!strcmp(base,"N"))) return 0;

      if (strcmp(base," ")&&strcmp(base,"\n")) nucs++;
	}
	if (!strcmp(base,"1")) break;
}
ct->numofbases = nucs - 1;
nucs--;
fclose (se);

if (nucs==0) return 0;

ct->allocate(nucs);

//now read the file
se=fopen(seqfile,"r");


do {
	fgets(temp,seqline,se);
	strncpy(test,temp,1);
	strcpy(test+1,"\0");
} while (!strcmp(test,";"));

//fgets(ct->ctlabel[1],seqline,se);
strcpy(ct->ctlabel[1],temp);

i = 1;
while (1) {
	fgets(seq,seqline,se);
	length = strlen (seq);
	for (j=0;j<length;j++) {//up to length-1 bcs of /n character
		strncpy (base,seq+j,1);
		strcpy (base+1,"\0");
		if (!strcmp(base,"1")) break;
      if (!(!strcmp(base,"A")||!strcmp(base,"a")||!strcmp(base,"C")||!strcmp(base,"c")
      	||!strcmp(base,"G")||!strcmp(base,"g")||!strcmp(base,"T")
         ||!strcmp(base,"t")||!strcmp(base,"U")||!strcmp(base,"u")
         ||!strcmp(base,"X")||!strcmp(base,"x")||!strcmp(base," ")||
         !strcmp(base,"\n")||!strcmp(base,"N"))) return 0;
      if ((!strcmp(base,"A")||!strcmp(base,"a")||!strcmp(base,"C")||!strcmp(base,"c")
      	||!strcmp(base,"G")||!strcmp(base,"g")||!strcmp(base,"T")
         ||!strcmp(base,"t")||!strcmp(base,"U")||!strcmp(base,"u")
         ||!strcmp(base,"X")||!strcmp(base,"x")||
         !strcmp(base,"N"))) {
				tonum(base,ct,(i));
      		ct->nucs[i]=base[0];
      		ct->hnumber[i] = i;
            i++;
         }
      
	}
	if (!strcmp(base,"1")) break;
}
//ct->numofbases = i - 1;

fclose (se);




return 1;
}

//outputs a ct file
void ctout (structure *ct,char *ctoutfile, int totalscore) { 
int count,i;//length
char line[2*ctheaderlength],number[2*numlen];//base[2]

FILE *ctfile;
ctfile=fopen(ctoutfile,"w");

strcpy(line,"");
sprintf(line,"Score %5i",totalscore);
fputs (line,ctfile);

for (count=1;count<=(ct->numofstructures);count++) {

	strcpy(line,"");
   sprintf(line,"%5i",ct->numofbases);
   
	//itoa(ct->numofbases,number,10);
	//strcat(line,number);

   
	sgifix  //this corrects a difference on the sgi computers
	if (ct->energy[count]!=0) {
   		strcat(line,"  ENERGY = ");

		gcvt((float (ct->energy[count]))/10,6,number);
	
   		strcat(line,number);
   		strcat(line,"  ");
	}
   else strcat(line,"  ");
	strcat(line,ct->ctlabel[count]);
	fputs (line,ctfile);
	for (i=1;i<ct->numofbases;i++) {
		/*if (ct->numseq[i]==1)
			sprintf(line,"%5i A%8i%5i%5i%5i\n",
				i,(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i]);
		else if (ct->numseq[i]==2)
			sprintf(line,"%5i C%8i%5i%5i%5i\n",
				i,(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i]);
		else if (ct->numseq[i]==3)
			sprintf(line,"%5i G%8i%5i%5i%5i\n",
				i,(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i]);
		else if (ct->numseq[i]==4)
			sprintf(line,"%5i U%8i%5i%5i%5i\n",
				i,(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i]);
		else if (ct->numseq[i]==0)
			sprintf(line,"%5i X%8i%5i%5i%5i\n",
				i,(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i]);
      else if (ct->numseq[i]==5)
			sprintf(line,"%5i I%8i%5i%5i%5i\n",
				i,(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i]); */
      sprintf(line,"%5i%2c%8i%5i%5i%5i\n",
				i,ct->nucs[i],(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i]);
		fputs(line,ctfile);
	}
   i = ct->numofbases;
   /*if (ct->numseq[i]==1)
   	sprintf(line,"%5i A%8i%5i%5i%5i\n",
   		i,(i-1),0,ct->basepr[count][i],ct->hnumber[i]);
   else if (ct->numseq[i]==2)
			sprintf(line,"%5i C%8i%5i%5i%5i\n",
				i,(i-1),0,ct->basepr[count][i],ct->hnumber[i]);
   else if (ct->numseq[i]==3)
			sprintf(line,"%5i G%8i%5i%5i%5i\n",
				i,(i-1),0,ct->basepr[count][i],ct->hnumber[i]);
   else if (ct->numseq[i]==4)
			sprintf(line,"%5i U%8i%5i%5i%5i\n",
				i,(i-1),0,ct->basepr[count][i],ct->hnumber[i]);
   else if (ct->numseq[i]==0)
			sprintf(line,"%5i X%8i%5i%5i%5i\n",
				i,(i-1),0,ct->basepr[count][i],ct->hnumber[i]);
   else if (ct->numseq[i]==5)
			sprintf(line,"%5i I%8i%5i%5i%5i\n",
				i,(i-1),0,ct->basepr[count][i],ct->hnumber[i]); */
	sprintf(line,"%5i %1c%8i%5i%5i%5i\n",
				i,ct->nucs[i],(i-1),0,ct->basepr[count][i],ct->hnumber[i]);
   fputs(line,ctfile);

}

fclose (ctfile);
return;
}


inline void swap(int *a,int *b) {
int temp;

	temp = *a;
   *a = *b;
   *b = temp;

return;

}

inline void swap (short int *a,short int *b) {
short int temp;

	temp = *a;
   *a = *b;
   *b = temp;

return;

}







void digit(int row,int column,int pos,arraystruct* table)
{
//puts the number Column in row Row
char number[10],chr[100],ch [7];
int length,i;

strcpy(number,"");
itoa(column,number,10);//convert column from integer to character
length = strlen(number);

if (pos-length<0) {
	strcpy(table->array[row][pos],".");
	return;
}

for (i=pos;i>=(pos-length+1);i--) {
	strcpy(chr,table->array[row][i]);
	if (strcmp(chr," ")) {
		strcpy(table->array[row][pos],".");
	return;
	}
}

for (i=0;i<=(length-1);i++) {
	strcpy(ch,"");
   strcpy(ch,number+i);
   strcpy(ch+1,"\0");
	strcpy((table->array[row][pos-length+1+i]),ch);
}
return;
}
char *tobase (int i)

{  //function is a look up table to convert the base
	// 	integer represention to a familiar character

	if (i==1) return "A";
	 else if (i==2) return "C";
	 else if (i==3) return "G";
	 else if (i==4) return "U";
	 else if (i==0) return "X";
    else if (i==5) return "I";
	 else return "?";   //to track down mistakes

}

void sortstructures (structure *ct) {//this function reorders the structures by
													//the efn energies
register int c;
int cur,i;
char ctheader[ctheaderlength];

for (c = 2; c<=(ct->numofstructures);c++){

	cur = c;

	while (cur>1) {
		if ((ct->energy[cur])<(ct->energy[cur-1])) {
      	swap(&ct->energy[cur],&ct->energy[cur-1]);
         //also swap the ct labels:
         strcpy(ctheader, ct->ctlabel[cur]);
         strcpy(ct->ctlabel[cur],ct->ctlabel[cur-1]);
         strcpy(ct->ctlabel[cur-1],ctheader);
         for (i=1;i<=(ct->numofbases);i++) {
         	swap(&ct->basepr[cur][i],&ct->basepr[cur-1][i]);
         }
         cur--;
   	}
		else {
		break;
		}
	}
}

}

void de_allocate (int **v,int i) {//this function deallocates the memory used
												//in an array
int a;

	for (a=0;a<i;a++) {
   	delete[] v[a];
   }
   delete[] v;
}


void de_allocate (short int **v,int i) {//this function deallocates the memory used
												//in an array
int a;

	for (a=0;a<i;a++) {
   	delete[] v[a];
   }
   delete[] v;
}

void de_allocate (bool **v,int i) {//this function deallocates the memory used
												//in an array
int a;

	for (a=0;a<i;a++) {
   	delete[] v[a];
   }
   delete[] v;
}    



void cctout( structure *ct, char *filename) {
    int i, j;
    ofstream out(filename);
    
    
    out << "-100\n";
    out << ct->numofbases<<"\n";
    out << ct->numofstructures<<"\n";
    //out << ct->ctlabel[1];
    for (i=1;i<=ct->numofbases;i++) {
	out << ct->numseq[i]<<"\n";
    }
    for (i=1;i<=ct->numofstructures;i++) {
	for (j=1;j<=ct->numofbases;j++) {
	    out << ct->basepr[i][j]<<"\n";
	}
    }



}







