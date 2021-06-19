
#if !defined(xdynalignh)
#define xdynalignh

#include "tprogressdialog.h"


short int einternal(int i,int j,int ip,int jp,structure *ct, datatable *data);
inline short int ehairpin(short i, short j, structure *ct, datatable *data);
inline short int edangle5(int i,int j,int ip,structure* ct,datatable* data);
inline short int edangle3(int i,int j,int ip,structure* ct,datatable* data);
inline short int ebp(int i,int j,int ip,int jp,structure *ct, datatable *data);


void xdynalign(structure *ct1, structure *ct2, structure *ct3, short *alignment1,short *alignment2,
	short int maxseparation, short int gapincrease, datatable *data, bool singleinsert,
	int *totalscore=NULL, TProgressDialog *progress=NULL);

void alignout(short* align1,short* align2,char *aout, structure *ct1, structure *ct2, structure *ct3,int totalscore=500);


#endif

