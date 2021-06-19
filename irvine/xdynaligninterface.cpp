#include "stdafx.h"
#include "xdynalign.cpp"



/*	Function getdat

	Function gets the names of data files to open

*/

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
		char *tloop, char *miscloop, char *danglef, char *int22,
      char *int21,char *coax, char *tstackcoax,
      char *coaxstack, char *tstack, char *tstackm,char *triloop, char *int11)

{
strcpy (loop,"loop.dat");
strcpy (stackf,"stack.dat");
strcpy (tstackh,"tstackh.dat");
strcpy (tstacki,"tstacki.dat");
strcpy (tloop,"tloop.dat");
strcpy (miscloop,"miscloop.dat");
strcpy (danglef,"dangle.dat");
strcpy (int22,"int22.dat");
strcpy (int21,"int21.dat");
strcpy (coax,"coaxial.dat");
strcpy (tstackcoax,"tstackcoax.dat");
strcpy (coaxstack,"coaxstack.dat");
strcpy (tstack,"tstack.dat");
strcpy (tstackm,"tstackm.dat");
strcpy (triloop,"triloop.dat");
strcpy (int11,"int11.dat");
}


int main(int argc, char* argv[])
{
	char inseq1[maxfil],inseq2[maxfil],inseq3[maxfil],outct[maxfil],outct2[maxfil],outct3[maxfil],aout[maxfil];
	structure ct1,ct2,ct3,ct_temp;
	char *alignment;
	int i,totalscore,imaxseparation,igapincrease;
	datatable data;
	char loop[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
		tloop[maxfil],miscloop[maxfil],danglef[maxfil],int22[maxfil],
      int21[maxfil],coax[maxfil],tstackcoax[maxfil],
      coaxstack[maxfil],tstack[maxfil],tstackm[maxfil],triloop[maxfil],int11[maxfil];
	short *align1,*align2;
	TProgressDialog progress;
	bool insert;
	float fgap;

	cout << "Enter the name of the first sequence:\n";
	cin >> inseq1;
	
	cout << "Enter the name of the second sequence:\n";
	cin >> inseq2;

	cout << "Enter the name of the third sequence:\n";
	cin >> inseq3;

	cout << "Enter the max separation:\n";
	cin >> imaxseparation;

	cout << "Enter the gap penalty:\n";
	cin >> fgap;

	igapincrease = (int) (fgap*10.0);

	cout << "Allow single BP inserts into one sequence? (1/0)";
	cin >> i;

	if (i) insert=true;
	else insert=false;
	
	openseq(&ct1,inseq1);
	openseq(&ct2,inseq2);
	openseq(&ct3,inseq3);

    i=1; // Assume sequence 1 is the smallest sequence
	if((ct2.numofbases<ct1.numofbases) && (ct2.numofbases<=ct3.numofbases)) i=2;
	else if((ct3.numofbases<ct1.numofbases) && (ct3.numofbases<ct2.numofbases)) i=3;

    if(i!=1)
	{
		if(i==2)
		{
		  cout << "Sawpping sequence 1 with sequence 2, because sequence 2 is the shortest sequence\n";
		  ct_temp=ct1; ct1=ct2; ct2=ct_temp;
        }
		else
		{
		  cout << "Sawpping sequence 1 with sequence 3, because sequence 3 is the shortest sequence\n";
		  ct_temp=ct1; ct1=ct3; ct3=ct_temp;
        }
    }

	cout << "Enter the name of the first output ct file:\n";
	cin >> outct;
	
	cout << "Enter the name of the second output ct file:\n";
	cin >> outct2;

	cout << "Enter the name of the third output ct file:\n";
	cin >> outct3;
	

	cout << "Enter name for the output of the alignment:\n";
	cin >> aout;
	

	getdat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,
   	int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11);

	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   	coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,&data)==0) {
      	cout << "A data file was lost";
         cin >> i;
   }



	//alignment = new char *[ct1.numofbases+ct2.numofbases];
	//for (i=1;i<ct1.numofbases+ct2.numofbases;i++) alignment[i] = new char [2];

	align1 = new short [ct1.numofbases+1];
	align2 = new short [ct2.numofbases+1];

	xdynalign(&ct1,&ct2,&ct3,align1,align2,imaxseparation,igapincrease,&data,insert,&totalscore,&progress);
	ctout(&ct1,outct,totalscore);   
	ctout(&ct2,outct2,totalscore);   
	ctout(&ct3,outct3,totalscore);  

	alignout(align1,align2,aout,&ct1,&ct2,&ct3,totalscore);



	//for (i=1;i<ct1.numofbases+ct2.numofbases;i++) delete[] alignment[i];
	//delete[] alignment;
   delete[] align1;
   delete[] align2;

	return 0;
}


void errmsg(int err,int erri) {

if (err==30) {
	cout << "End Reached at traceback #"<<erri<<"\n";
   
}
if (err==100) {
	cout << "error # "<<erri;
   
}
switch (err) {
	case 1:
   	cout << "Could not allocate enough memory";
      break;
   case 2:
   	cout << "Too many possible base pairs";
      break;
   case 3:
   	cout << "Too many helixes in multibranch loop";
   case 4:
   	cout << "Too many structures in CT file";
   default:
   	cout << "Unknown error";
}
cin >> err;

return;

}

void update(int i) {

	cout<< i<<"\n"<<flush;
}





