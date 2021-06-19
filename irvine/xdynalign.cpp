/*                               -*- Mode: C -*- 
 * xdynalign.c --- simultaneous alignment and structure prediction of 3 RNA sequences
 * Author          : Beeta Masoumi
 * Created On      : July 2004
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Fri Sep 16 19:01:03 2005
 *
 * This copyrighted source code is freely distributed under the terms
 * of the GNU General Public License. 
 *
 * See the files COPYRIGHT and LICENSE for details.
 *
 * - Masoumi, B. and Turcotte, M. (2005) Simultaneous alignment and structure 
 *       prediction of three RNA sequences. Int. J. Bioinformatics Research and 
 *       Applications. Vol. 1, No. 2, pp. 230-245 
 *
 * - Masoumi, B. and Turcotte, M. (2005) Simultaneous alignment and structure 
 *       prediction of RNAs: Are three input sequences better than two? 
 *       In S. V. Sunderam et al., editor, 2005 International Conference on 
 *       Computational Science (ICCS 2005), Lecture Notes in Computer Science 3515,
 *       pages 936-943, Atlanta, USA, May 22-25 2005. 
 *
 * This work is based on the source code of Dynalign:
 *
 * RNA Secondary Structure Prediction, Using the Algorithm of Zuker
 * C++ version by David Mathews
 * Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004
 * Programmed for Isis Pharmaceuticals and 
 * the Turner Lab Chemistry Department, University of Rochester
 */

#include "stdafx.h"
#include "algorithm.cpp"
#include "xdynalign.h"
#include <math.h>
#include <string.h>
// #include <assert.h>
// #include <iostream.h>

// #define NDEBUG

#define maxloop 20 //maximum size of unpaired nucs on each side of an internal loop


// This function takes 3 sequences and creates an alignment. It requires:
// maxseparation (the largest difference in position that is permissible to align),
// maxgap (the maximum length of a gap),
// comp (a bonus for compensating changes),
// gapstart (an affine gap penalty for starting a gap), and
// gapincrease (an affine gap penalty for each gapped position)
// datatable is the thermodynamic data files.

void xdynalign(structure *ct1, structure *ct2, structure *ct3, short *alignment1, short *alignment2,
	       short int maxseparation, short int gapincrease, datatable *data,
	       bool singleinsert, int *totalscore, TProgressDialog *progress) {
    
  register short int i,j,k,l,m,n,a,b,aa,bb,size,c,d,e,f,g,h,en1,en2,en3,maxsep,en4,lowest,gap,iprime,kprime,mprime,u,aprime,asegond,p,q,nogaps;
  register short int ebp_i_j_ip1_jm1,ebp_ip1_jm1_ip2_jm2,einternal_i_j_c_d,ebp_k_l_kp1_lm1,ebp_kp1_lm1_kp2_lm2,einternal_k_l_e_f,ebp_m_n_mp1_nm1,
    ebp_mp1_nm1_mp2_nm2,ehairpin_i_j,ehairpin_k_l,ehairpin_m_n,
    l1,l2,l3,l4,l5,l6,max_l1_l3,tri_max_l1_l3_l5,max_l2_l4,
    l1pl2,l1pl2pl3,l1pl2pl3pl4,l1pl2pl3pl4pl5;

  short int ***mine;
  short int ******v,******w;
  short int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
		       {0,1,0,1,0,0},{0,0,0,0,0,0}};
  bool **pair;
  register bool found;
    
  //v[j][i][l][k][m][n] is the sum of the lowest free energies for fragment i to j and k to l and m to n
  //with i and j paired, k and l paired, m and n paired, and i aligned to k and j aligned to l and m aligned to n
    
  //w[j][i][l][k][m][n] is the sum of the lowest free energies for fragment i to j and k to l and m to n
  //with the fragments interior to multibranch loops and i aligned to k and m and j aligned to l and n
    
  //mine[i][k][m] is the minumum energy of the fragment of 1 to i with i aligned to k and m
  //with no constraints of the configuration
    
    
  //pair is a bool array to keep track of whether a nucleotide is forced single-stranded
  pair = new bool *[3];
  pair[0] = new bool [ct1->numofbases+1];
  pair[1] = new bool [ct2->numofbases+1];
  pair[2] = new bool [ct3->numofbases+1];
  for (i=0;i<=ct1->numofbases;i++) pair[0][i] = true;
  for (i=0;i<=ct2->numofbases;i++) pair[1][i] = true;
  for (i=0;i<=ct3->numofbases;i++) pair[2][i] = true;
    
  for(i=1;i<ct1->nnopair;i++) pair[0][ct1->nopair[i]]=false;
  for(i=1;i<ct2->nnopair;i++) pair[1][ct2->nopair[i]]=false;
  for(i=1;i<ct3->nnopair;i++) pair[2][ct3->nopair[i]]=false;	
    
  //allocate space in an array
  //refer to [j][i][l][k][n][m] because i<j so the i dimension can be smaller
  v = new short int *****[ct1->numofbases+1];
  for (i=0;i<=ct1->numofbases;i++) {
    v[i] = new short int ****[i];
    for (j=0;j<i;j++) {
      v[i][j] = new short int ***[2*maxseparation+2];
      for (k=0;k<2*maxseparation+2;k++) {
	v[i][j][k] = new short int **[2*maxseparation+2];
	for (l=0;l<2*maxseparation+2;l++) {
	  v[i][j][k][l] = new short int *[2*maxseparation+2];
	  for(m=0;m<2*maxseparation+2;m++) {
	    v[i][j][k][l][m] = new short int [2*maxseparation+2];
	    for(n=0;n<2*maxseparation+2;n++) v[i][j][k][l][m][n] = infinity;
	  }
	}
      }
    }
  }

  w = new short int *****[ct1->numofbases+1];
  for (i=0;i<=ct1->numofbases;i++) {
    w[i] = new short int ****[i];
    for (j=0;j<i;j++) {
      w[i][j] = new short int ***[2*maxseparation+2];
      for (k=0;k<2*maxseparation+2;k++) {
	w[i][j][k] = new short int **[2*maxseparation+2];
	for (l=0;l<2*maxseparation+2;l++) {
	  w[i][j][k][l] = new short int *[2*maxseparation+2];
	  for(m=0;m<2*maxseparation+2;m++) {
	    w[i][j][k][l][m] = new short int [2*maxseparation+2];
	    for(n=0;n<2*maxseparation+2;n++) w[i][j][k][l][m][n] = infinity;
	  }
	}
      }
    }
  }

  cout << "initialization done \n"<<flush;
  maxsep = maxseparation; //place maxseparation into a register int
  gap = gapincrease; //place the gap penalty into a register short
    
  mine = new short int **[ct1->numofbases+1];
  for (i=0;i<=ct1->numofbases;i++) {
    mine[i] = new short int *[2*maxsep+2];
    for (k=0; k<2*maxsep+2;k++) {
      mine[i][k] = new short int [2*maxsep+2];
    }
  }
  cout << "start 6 nested loops for filling v and w \n"<<flush;
    
  /*********************************start fill V and W ************************************/
  //look for alignment of every i paired to j in sequence 1 to k paired to l
  //in sequence 2 to m paired to n in sequence 3
  for (size=minloop;size<=ct1->numofbases;size++) {
    if (progress!=NULL&&size%10==0) progress->update(100*size/ct1->numofbases);
    cout <<size<<"\n"<<flush;
    for (i=1;i+size<=ct1->numofbases+1;i++) {
      j = i+size-1;	//size=j-i+1
      // Pre-calculate some energies for speeding-up
      ebp_i_j_ip1_jm1    =ebp(i,j,i+1,j-1,ct1, data);
      ebp_ip1_jm1_ip2_jm2=ebp( i+1,j-1,i+2,j-2,ct1,data );
      ehairpin_i_j=ehairpin(i,j,ct1,data);
      for (k=min(i+maxsep,ct2->numofbases);k>=i-maxsep&&k>=1;k--) {  
	a  = k-i+maxsep;
	for (l=max(max(j-maxsep,1),k+minloop);l<=j+maxsep&&l<=ct2->numofbases;l++) {  
	  b  = l-j+maxsep;
	  // Pre-calculate some energies for speeding-up
	  if(k< ct2->numofbases   ) ebp_k_l_kp1_lm1    =ebp(k,l,k+1,l-1,ct2, data);
	  if(k<(ct2->numofbases)-1) ebp_kp1_lm1_kp2_lm2=ebp( k+1,l-1,k+2,l-2,ct2,data );
	  ehairpin_k_l=ehairpin(k,l,ct2,data);
	  for (m=min(k+maxsep,ct3->numofbases);m>=k-maxsep&&m>=1;m--) {  
	    aa = m-k+maxsep;
	    for (n=max(max(l-maxsep,1),m+minloop);n<=l+maxsep&&n<=ct3->numofbases;n++) {

	      // Pre-calculate some energies for speeding-up
	      if(m< ct3->numofbases   ) ebp_m_n_mp1_nm1    =ebp(m,n,m+1,n-1,ct3, data); 
	      if(m<(ct3->numofbases)-1) ebp_mp1_nm1_mp2_nm2=ebp( m+1,n-1,m+2,n-2,ct3,data );
	      ehairpin_m_n=ehairpin(m,n,ct3,data);

	      //use a,aa,b and bb when refering to the energy arrays
	      bb = n-l+maxsep;
                            
                            
	      //filter out isolated base pairs	Note: e(f/g) can be 0 or 1 or 2
	      e = inc[ct1->numseq[i+1]][ct1->numseq[j-1]];
	      if (i-1>0&&j+1<=ct1->numofbases) e = e + inc[ct1->numseq[i-1]][ct1->numseq[j+1]];
                            
                            
	      f = inc[ct2->numseq[k+1]][ct2->numseq[l-1]];
	      if (k-1>0&&l+1<=ct2->numofbases) f = f + inc[ct2->numseq[k-1]][ct2->numseq[l+1]];
                            
                            
	      g = inc[ct3->numseq[m+1]][ct3->numseq[n-1]];
	      if (m-1>0&&n+1<=ct3->numofbases) g = g + inc[ct3->numseq[m-1]][ct3->numseq[n+1]];
                            
	      /********************************** start fill V *******************************************/
	      //Fill v, the best energy with i-j paired, k-l paired, m-n paired
	      //and their pairs aligned:
	      if (inc[ct1->numseq[i]][ct1->numseq[j]]&&
		  inc[ct2->numseq[k]][ct2->numseq[l]]&&
		  inc[ct3->numseq[m]][ct3->numseq[n]]&&
		  pair[0][i]&&pair[0][j]&&pair[1][k]&&pair[1][l]&&pair[2][m]&&pair[2][n]&&e&&f&&g) {
                                
                                
		/**************************************** V1 ************************************************/
		//Now consider hairpins for all
		nogaps=3*max(j-i,max(l-k,n-m))-j+i-l+k-n+m;	
		v[j][i][a][b][aa][bb] = ehairpin_i_j+ehairpin_k_l+ehairpin_m_n+gap*nogaps; 
                                
		/**************************************** V2 = en4 ******************************************/
		en4=infinity;
		//Now consider internal loops/one stacked and other not
		for (c=i+1;c<=i+maxloop&&c<j-minloop;c++) {
		  l1=c-i;
		  for (d=j-1;d>=j-maxloop&&d>c+minloop;d--) {
		    l2=j-d;
		    l1pl2=l1+l2;

		    // Cases for Seq1 
		    if ( c==i+1 && d==j-1 ) { // helical stacking 
		      en1= ebp_i_j_ip1_jm1;
		    } else if ( singleinsert && // base pair insertion 
				c==i+2 && d==j-2 && j-2>i+2 &&
				inc[ ct1->numseq[i+1] ][ ct1->numseq[j-1] ] &&
				inc[ ct1->numseq[i+2] ][ ct1->numseq[j-2] ] &&
				( ( e==k+1 && f==l-1 && inc[ ct2->numseq[e] ][ ct2->numseq[f] ] ) || 
				  ( g==m+1 && h==n-1 && inc[ ct3->numseq[g] ][ ct3->numseq[h] ] ) ) ) {
		      en1 = ebp_i_j_ip1_jm1 + ebp_ip1_jm1_ip2_jm2;
		    } else { // internal loop 
		      en1=einternal_i_j_c_d;
		    }
		    /************************************************************************************************/
		    /*										// Cases for Seq1                                         
												if ( c==i+1 && d==j-1 ) { // helical stacking 
												en1= ebp_i_j_ip1_jm1;
												} else { // internal loop 
												en1 = einternal(i,j,c,d,ct1,data);
												}
		    */
		    /************************************************************************************************/
		    for (e=max(k+1,c-maxsep);e<=k+maxloop&&e<l-minloop&&(e-c)<=maxsep;e++) {
		      l3=e-k;
		      max_l1_l3=max(l1,l3);
		      l1pl2pl3=l1pl2+l3;
		      for (f=min(l-1,d+maxsep);f>=l-maxloop&&f>e+minloop&&(d-f)<=maxsep;f--) {
			l4=l-f;
			l1pl2pl3pl4=l1pl2pl3+l4;
			max_l2_l4=max(l2,l4);
			// Cases for Seq2 
			if ( e==k+1 && f==l-1 ) { // helical stacking 
			  en2=ebp_k_l_kp1_lm1;
			} else if ( singleinsert && // base pair insertion 
				    e==k+2 && f==l-2 && l-2> k+2 && 
				    inc[ ct2->numseq[k+1] ][ ct2->numseq[l-1] ] &&
				    inc[ ct2->numseq[k+2] ][ ct2->numseq[l-2] ] && 
				    ( ( c==i+1 && d==j-1 && inc[ ct1->numseq[c] ][ ct1->numseq[d] ] ) || 
				      ( g==m+1 && h==n-1 && inc[ ct3->numseq[g] ][ ct3->numseq[h] ] ) ) ) {
			  en2 = ebp_k_l_kp1_lm1+ebp_kp1_lm1_kp2_lm2;
			} else { // internal loop 
			  en2 = einternal_k_l_e_f;
			}

			/************************************************************************************************/											    
			/*												// Cases for Seq2 
															if ( e==k+1 && f==l-1 ) { // helical stacking 
															en2=ebp_k_l_kp1_lm1;
															} else { // internal loop 
															en2 = einternal(k,l,e,f,ct2,data);
															}
			*/
			/************************************************************************************************/
			for (g=max(m+1,e-maxsep);g<=m+maxloop&&g<n-minloop&&(g-e)<=maxsep;g++) {
			  l5=g-m;
			  l1pl2pl3pl4pl5=l1pl2pl3pl4+l5;
			  tri_max_l1_l3_l5=3*max(max_l1_l3,l5);
			  for (h=min(n-1,f+maxsep);h>=n-maxloop&&h>g+minloop&&(f-h)<=maxsep;h--) {
			    l6=n-h;
			    // Cases for Seq3
			    if ( g==m+1 && h==n-1 ) { // helical stacking 
			      en3 = ebp_m_n_mp1_nm1;
			    } else if ( singleinsert && // base pair insertion 
					g==m+2 && h==n-2 && n-2>m+2 &&
					inc[ ct3->numseq[m+1] ][ ct3->numseq[n-1] ] &&
					inc[ ct3->numseq[m+2] ][ ct3->numseq[n-2] ] &&
					( ( c==i+1 && d==j-1 && inc[ ct1->numseq[c] ][ ct1->numseq[d] ] ) || 
					  ( e==k+1 && f==l-1 && inc[ ct2->numseq[e] ][ ct2->numseq[f] ] ) ) ) {
			      en3 = ebp_m_n_mp1_nm1 + ebp_mp1_nm1_mp2_nm2;
			    } else { // internal loop 
			      en3 = einternal(m,n,g,h,ct3,data);
			    }

			    /************************************************************************************************/
			    /*													    // Cases for Seq3
																    if ( g==m+1 && h==n-1 ) { // helical stacking 
																    en3 = ebp_m_n_mp1_nm1;
																    } else { // internal loop 
																    en3 = einternal(m,n,g,h,ct3,data);
																    }
			    */
			    /************************************************************************************************/
			    /*														
			    // Cases for Seq1 
			    if ( c==i+1 && d==j-1 ) { // helical stacking 


			    en1= ebp_i_j_ip1_jm1;

			    } else if ( singleinsert && // base pair insertion 
			    c==i+2 && d==j-2 && j-2>i+2 &&
			    inc[ ct1->numseq[i+1] ][ ct1->numseq[j-1] ] &&
			    inc[ ct1->numseq[i+2] ][ ct1->numseq[j-2] ] &&
			    ( ( e==k+1 && f==l-1 && inc[ ct2->numseq[e] ][ ct2->numseq[f] ] ) || 
			    ( g==m+1 && h==n-1 && inc[ ct3->numseq[g] ][ ct3->numseq[h] ] ) ) ) {
			    

			    en1 = ebp_i_j_ip1_jm1 + ebp_ip1_jm1_ip2_jm2;


			    } else { // internal loop 


			    en1=einternal_i_j_c_d;

			    }

			    // Cases for Seq2 

			    if ( e==k+1 && f==l-1 ) { // helical stacking 
	      

			    en2=ebp_k_l_kp1_lm1;

			    } else if ( singleinsert && // base pair insertion 
			    e==k+2 && f==l-2 && l-2> k+2 && 
			    inc[ ct2->numseq[k+1] ][ ct2->numseq[l-1] ] &&
			    inc[ ct2->numseq[k+2] ][ ct2->numseq[l-2] ] && 
			    ( ( c==i+1 && d==j-1 && inc[ ct1->numseq[c] ][ ct1->numseq[d] ] ) || 
			    ( g==m+1 && h==n-1 && inc[ ct3->numseq[g] ][ ct3->numseq[h] ] ) ) ) {


			    en2 = ebp_k_l_kp1_lm1+ebp_kp1_lm1_kp2_lm2;

			    } else { // internal loop 


			    en2 = einternal_k_l_e_f;

			    }

			    // Cases for Seq3
															
			    if ( g==m+1 && h==n-1 ) { // helical stacking 


			    en3 = ebp_m_n_mp1_nm1;

			    } else if ( singleinsert && // base pair insertion 
			    g==m+2 && h==n-2 && n-2>m+2 &&
			    inc[ ct3->numseq[m+1] ][ ct3->numseq[n-1] ] &&
			    inc[ ct3->numseq[m+2] ][ ct3->numseq[n-2] ] &&
			    ( ( c==i+1 && d==j-1 && inc[ ct1->numseq[c] ][ ct1->numseq[d] ] ) || 
			    ( e==k+1 && f==l-1 && inc[ ct2->numseq[e] ][ ct2->numseq[f] ] ) ) ) {


			    en3 = ebp_m_n_mp1_nm1 + ebp_mp1_nm1_mp2_nm2;

			    } else { // internal loop 

			    en3 = einternal(m,n,g,h,ct3,data);

			    }
			    */
			    nogaps=tri_max_l1_l3_l5+3*max(max_l2_l4,l6)-l1pl2pl3pl4pl5-l6;
			    en4 = min(en4,en1+en2+en3+v[d][c][e-c+maxsep][f-d+maxsep][g-e+maxsep][h-f+maxsep]+gap*nogaps); 

			  }
			}
		      }
		    }
		  }
		}
  							
		v[j][i][a][b][aa][bb] = min( en4, v[j][i][a][b][aa][bb] );                                                        

		en1=infinity;
		/********************************************* V3 = en1 ****************************************/
		//Now consider multiloop:
		//junction closed by i-j pair aligned with a and b and aa and bb
		//calculate the free energy of 3 fragments merged:
                                
		for (c=i+minloop+1;c<j-minloop;c++) {
		  for (d=max(k+minloop+1,c-maxsep);d<l-minloop&&d<=c+maxsep;d++) {
		    e = d-c+maxsep;
		    for(f=max(m+minloop+1,d-maxsep);f<n-minloop&&f<=d+maxsep;f++) {
		      g = f-d+maxsep;
                                            
                                            
		      //must consider whether an unpaired nucleotide is stacked
		      //is stacked onto each of the nucleotides in a helix
		      //There are 64 cases:
		      //consider them in this order (0 unstacked, 1 stacked)
		      //		i	j	k	l	m	n
		      //1		0	0	0	0	0	0
		      //2		0	0	0	0	0	1
		      //3		0	0	0	0	1	0
		      //4		0	0	0	0	1	1
		      //5		0	0	0	1	0	0
		      //6		0	0	0	1	0	1
		      //7		0	0	0	1	1	0
		      //8		0	0	0	1	1	1
		      //9		0	0	1	0	0	0
		      //10	0	0	1	0	0	1
		      //11	0	0	1	0	1	0
		      //12	0	0	1	0	1	1
		      //13	0	0	1	1	0	0
		      //14	0	0	1	1	0	1
		      //15	0	0	1	1	1	0
		      //16	0	0	1	1	1	1
		      //17	0	1	0	0	0	0
		      //18	0	1	0	0	0	1
		      //19	0	1	0	0	1	0
		      //20	0	1	0	0	1	1
		      //21	0	1	0	1	0	0
		      //22	0	1	0	1	0	1
		      //23	0	1	0	1	1	0
		      //24	0	1	0	1	1	1
		      //25	0	1	1	0	0	0
		      //26	0	1	1	0	0	1
		      //27	0	1	1	0	1	0
		      //28	0	1	1	0	1	1
		      //29	0	1	1	1	0	0
		      //30	0	1	1	1	0	1
		      //31	0	1	1	1	1	0
		      //32	0	1	1	1	1	1
		      //33	1	0	0	0	0	0
		      //34	1	0	0	0	0	1
		      //35	1	0	0	0	1	0
		      //36	1	0	0	0	1	1
		      //37	1	0	0	1	0	0
		      //38	1	0	0	1	0	1
		      //39	1	0	0	1	1	0
		      //40	1	0	0	1	1	1
		      //41	1	0	1	0	0	0
		      //42	1	0	1	0	0	1
		      //43	1	0	1	0	1	0
		      //44	1	0	1	0	1	1
		      //45	1	0	1	1	0	0
		      //46	1	0	1	1	0	1
		      //47	1	0	1	1	1	0
		      //48	1	0	1	1	1	1
		      //49	1	1	0	0	0	0
		      //50	1	1	0	0	0	1
		      //51	1	1	0	0	1	0
		      //52	1	1	0	0	1	1
		      //53	1	1	0	1	0	0
		      //54	1	1	0	1	0	1
		      //55	1	1	0	1	1	0
		      //56	1	1	0	1	1	1
		      //57	1	1	1	0	0	0
		      //58	1	1	1	0	0	1
		      //59	1	1	1	0	1	0
		      //60	1	1	1	0	1	1
		      //61	1	1	1	1	0	0
		      //62	1	1	1	1	0	1
		      //63	1	1	1	1	1	0
		      //64	1	1	1	1	1	1
                                            

                                            
                                            
		      // Case			a	b	aa	bb
		      //-----			---	---	---	---
		      // 39			-1	-1	1	1
		      // 37			-1	-1	0	1
		      // 40			-1	-1	1	0
		      // 38			-1	-1	0	0
                                            
		      // 52			-1	1	1	-1
		      // 50			-1	1	0	-1
		      // 51			-1	1	1	0
		      // 49			-1	1	0	0
                                            
		      // 35,56		-1	0	1	0
		      // 36			-1	0	1	-1
		      // 34			-1	0	0	-1
		      // 55			-1	0	1	1
		      // 53			-1	0	0	1
		      // 33,54		-1	0	0	0
                                            
		      // 14			1	-1	-1	0
		      // 16			1	-1	0	0
		      // 13			1	-1	-1	1
		      // 15			1	-1	0	1
                                            
		      // 26			1	1	-1	-1
		      // 28			1	1	0	-1
		      // 25			1	1	-1	0
		      // 27			1	1	0	0
                                            
		      // 10			1	0	-1	-1
		      // 12			1	0	0	-1
		      // 29			1	0	-1	1
		      // 31			1	0	0	1
		      // 9,30			1	0	-1	0
		      // 11,32		1	0	0	0
                                            
		      // 45			0	-1	-1	1
		      // 7			0	-1	1	1
		      // 5,47			0	-1	0	1
		      // 46			0	-1	-1	0
		      // 8			0	-1	1	0
		      // 6,48			0	-1	0	0
                                            
		      // 58			0	1	-1	-1
		      // 20			0	1	1	-1
		      // 18,60		0	1	0	-1
		      // 57			0	1	-1	0
		      // 19			0	1   1	0
		      // 17,59		0	1	0	0
                                            
		      // 42			0	0	-1	-1
		      // 4			0	0	1	-1
		      // 2,44			0	0	0	-1
		      // 61			0	0	-1	1
		      // 41,62		0	0	-1	0
		      // 3,24			0	0	1	0
		      // 23			0	0	1	1
		      // 21,63		0	0	0	1
		      // 1,22,43,64	0	0	0	0
                                            

		      if(b-1>=0) {
                                                
			if(a-1>0) {
                                                    
			  if(bb<2*maxsep+1)   
			    {
                                                        
			      if(aa<2*maxsep+1)		
				{
				  //case 39 - i, l and m stacked
				  en1 = min(en1,w[c][i+2][a-1][e][aa+1][g]
					    +w[j-1][c+1][e][b-1][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
					    +3*data->eparam[6]+edangle3(i,j,i+1,ct1,data)
					    +edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)+3*gap);
														
				}
                                                        
			      //case 37 - i and l stacked
			      en1 = min(en1,w[c][i+2][a-1][e][aa][g]
					+w[j-1][c+1][e][b-1][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
					+2*data->eparam[6]+edangle3(i,j,i+1,ct1,data)
					+edangle5(l,k,l-1,ct2,data)+4*gap);
														
			    }
                                                    
			  if(aa<2*maxsep+1)  // aa+1<2*maxsep+2
			    {
			      //case 40 - i, l, m, and n stacked
			      en1 = min(en1,w[c][i+2][a-1][e][aa+1][g]
					+w[j-1][c+1][e][b-1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
					+4*data->eparam[6]+edangle3(i,j,i+1,ct1,data)
					+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)
					+edangle5(n,m,n-1,ct3,data)+2*gap);
														
			    }
                                                    
			  //case 38 - i, l, and n stacked
			  en1 = min(en1,w[c][i+2][a-1][e][aa][g]
				    +w[j-1][c+1][e][b-1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				    +3*data->eparam[6]+edangle3(i,j,i+1,ct1,data)
				    +edangle5(l,k,l-1,ct2,data)+edangle5(n,m,n-1,ct3,data)+3*gap);
													
			}
                                                
			if(a<2*maxsep+1)  // a+1<2*maxsep+2
			  {
                                                    
			    if(bb+1<2*maxsep+2) {
			      //case 15 - k, l and m stacked
			      en1 = min(en1,w[c][i+1][a+1][e][aa][g]
					+w[j-1][c+1][e][b-1][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
					+3*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
					+edangle3(k,l,k+1,ct2,data)
					+edangle3(m,n,m+1,ct3,data)+3*gap);
														
			    }
			    //case 16 - k, l, m and n stacked
			    en1 = min(en1,w[c][i+1][a+1][e][aa][g]
				      +w[j-1][c+1][e][b-1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				      +4*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
				      +edangle5(n,m,n-1,ct3,data)
				      +edangle3(k,l,k+1,ct2,data)
				      +edangle3(m,n,m+1,ct3,data)+2*gap);
													
                                                    
			    if(aa-1>0) {
			      //case 14 - k, l and n stacked
			      en1 = min(en1,w[c][i+1][a+1][e][aa-1][g]
					+w[j-1][c+1][e][b-1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
					+3*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
					+edangle3(k,l,k+1,ct2,data)
					+edangle5(n,m,n-1,ct3,data)+3*gap);
														
			      if(bb+1<2*maxsep+2) {
				//case 13 - k and l stacked
				en1 = min(en1,w[c][i+1][a+1][e][aa-1][g]
					  +w[j-1][c+1][e][b-1][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
					  +2*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
					  +edangle3(k,l,k+1,ct2,data)+4*gap);
															
			      }
			    }
                                                    
			  }
                                                
			if(bb+1<2*maxsep+2) {
                                                    
			  if(aa-1>0) {
			    //case 45 - i, k and l stacked
			    en1 = min(en1,w[c][i+2][a][e][aa-1][g]
				      +w[j-1][c+1][e][b-1][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
				      +3*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
				      +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+3*gap);
														
			  }
                                                    
			  if(aa+1<2*maxsep+2) {
			    //case 7 - l and m stacked
			    en1 = min(en1,w[c][i+1][a][e][aa+1][g]
				      +w[j-1][c+1][e][b-1][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
				      +2*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
				      +edangle3(m,n,m+1,ct3,data)+4*gap);
														
			  }
                                                    
			  //case 5 - stack on l
			  en1 = min(en1,w[c][i+1][a][e][aa][g]
				    +w[j-1][c+1][e][b-1][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
				    +data->eparam[6]+edangle5(l,k,l-1,ct2,data)+2*gap);
													
                                                    
			  //case 47 - i, k, l and m stacked
			  en1 = min(en1,w[c][i+2][a][e][aa][g]
				    +w[j-1][c+1][e][b-1][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
				    +4*data->eparam[6]
				    +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)
				    +edangle5(l,k,l-1,ct2,data)+edangle3(m,b,m+1,ct3,data)+2*gap);
													
			}
                                                
			if(aa-1>0) {
			  //case 46 - i, k, l and n stacked
			  en1 = min(en1,w[c][i+2][a][e][aa-1][g]
				    +w[j-1][c+1][e][b-1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				    +4*data->eparam[6]
				    +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)
				    +edangle5(l,k,l-1,ct2,data)+edangle5(n,m,n-1,ct3,data)+2*gap);
													
			}
                                                
			if(aa+1<2*maxsep+2) {
			  //case 8 - l, m and n stacked
			  en1 = min(en1,w[c][i+1][a][e][aa+1][g]
				    +w[j-1][c+1][e][b-1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				    +3*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
				    +edangle3(m,n,m+1,ct3,data)
				    +edangle5(n,m,n-1,ct3,data)+3*gap);
													
			}
                                                
			//case 6 - l and n stacked
			en1 = min(en1,w[c][i+1][a][e][aa][g]
				  +w[j-1][c+1][e][b-1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				  +2*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
				  +edangle5(n,m,n-1,ct3,data)+gap);
												
                                                
			//case 48 - i, k, l, m and n stacked
			en1 = min(en1,w[c][i+2][a][e][aa][g]
				  +w[j-1][c+1][e][b-1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				  +5*data->eparam[6]+edangle3(m,n,m+1,ct3,data)
				  +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)
				  +edangle5(l,k,l-1,ct2,data)+edangle5(n,m,n-1,ct3,data)+gap);
												
		      }
                                            
		      if(b+1<2*maxsep+2) {
                                                
			if(a-1>0) {
                                                    
			  if(bb-1>=0) {
                                                        
			    if(aa+1<2*maxsep+2) {
			      //case 52 - i, j, m and n stacked
			      en1 = min(en1,w[c][i+2][a-1][e][aa+1][g]
					+w[j-2][c+1][e][b+1][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
					+4*data->eparam[6]
					+edangle5(j,i,j-1,ct1,data)+edangle5(n,m,n-1,ct3,data)
					+edangle3(i,j,i+1,ct1,data)+edangle3(m,n,m+1,ct3,data)+2*gap);
															
			    }
                                                        
			    //case 50 - i, j and n stacked
			    en1 = min(en1,w[c][i+2][a-1][e][aa][g]
				      +w[j-2][c+1][e][b+1][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
				      +3*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				      +edangle3(i,j,i+1,ct1,data)+edangle5(n,m,n-1,ct3,data)+3*gap);
														
			  }
                                                    
			  if(aa+1<2*maxsep+2) {
			    //case 51 - i, j and m stacked
			    en1 = min(en1,w[c][i+2][a-1][e][aa+1][g]
				      +w[j-2][c+1][e][b+1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				      +3*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				      +edangle3(i,j,i+1,ct1,data)+edangle3(m,n,m+1,ct3,data)+3*gap);
														
			  }
                                                    
			  //case 49 - i and j stacked
			  en1 = min(en1,w[c][i+2][a-1][e][aa][g]
				    +w[j-2][c+1][e][b+1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				    +2*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				    +edangle3(i,j,i+1,ct1,data)+4*gap);
													
			}
                                                
			if(a+1<2*maxsep+2) {
                                                    
			  if(bb-1>=0) {
                                                        
			    if(aa-1>0) {
			      //case 26 - j, k and n stacked
			      en1 = min(en1,w[c][i+1][a+1][e][aa-1][g]
					+w[j-2][c+1][e][b+1][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
					+3*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
					+edangle5(n,m,n-1,ct3,data)
					+edangle3(k,l,k+1,ct2,data)+3*gap);
															
			    }
                                                        
			    //case 28 - j, k, m and n stacked
			    en1 = min(en1,w[c][i+1][a+1][e][aa][g]
				      +w[j-2][c+1][e][b+1][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
				      +4*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				      +edangle5(n,m,n-1,ct3,data)
				      +edangle3(k,l,k+1,ct2,data)
				      +edangle5(m,n,m+1,ct3,data)+2*gap);
														
			  }
                                                    
			  if(aa-1>0) {
			    //case 25 - j and k stacked
			    en1 = min(en1,w[c][i+1][a+1][e][aa-1][g]
				      +w[j-2][c+1][e][b+1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				      +2*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				      +edangle3(k,l,k+1,ct2,data)+4*gap);
														
			  }
                                                    
			  //case 27 - j, k and m stacked
			  en1 = min(en1,w[c][i+1][a+1][e][aa][g]
				    +w[j-2][c+1][e][b+1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				    +3*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				    +edangle3(k,l,k+1,ct2,data)
				    +edangle3(m,n,m+1,ct3,data)+3*gap);
													
			}
                                                
			if(bb-1>=0) {
                                                    
			  if(aa-1>0) {
			    //case 58 - i, j, k and n stacked
			    en1 = min(en1,w[c][i+2][a][e][aa-1][g]
				      +w[j-2][c+1][e][b+1][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
				      +4*data->eparam[6]+edangle3(k,l,k+1,ct2,data)
				      +edangle5(j,i,j-1,ct1,data)+edangle5(n,m,n-1,ct3,data)
				      +edangle3(i,j,i+1,ct1,data)+2*gap);
														
			  }
                                                    
			  if(aa+1<2*maxsep+2) {
			    //case 20 - j, m and n stacked
			    en1 = min(en1,w[c][i+1][a][e][aa+1][g]
				      +w[j-2][c+1][e][b+1][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
				      +3*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				      +edangle5(n,m,n-1,ct3,data)
				      +edangle3(m,n,m+1,ct3,data)+3*gap);
														
			  }
                                                    
			  //case 18 - j and n stacked
			  en1 = min(en1,w[c][i+1][a][e][aa][g]
				    +w[j-2][c+1][e][b+1][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
				    +2*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				    +edangle5(n,m,n-1,ct3,data)+gap);
													
                                                    
			  //case 60 - i, j, k, m and n stacked
			  en1 = min(en1,w[c][i+2][a][e][aa][g]
				    +w[j-2][c+1][e][b+1][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
				    +5*data->eparam[6]+edangle3(k,l,k+1,ct2,data)
				    +edangle5(j,i,j-1,ct1,data)+edangle3(m,n,m+1,ct3,data)
				    +edangle3(i,j,i+1,ct1,data)+edangle5(n,m,n-1,ct3,data)+gap);
													
			}
                                                
			if(aa-1>0) {
			  //case 57 - i, j and k stacked
			  en1 = min(en1,w[c][i+2][a][e][aa-1][g]
				    +w[j-2][c+1][e][b+1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				    +3*data->eparam[6]+edangle3(k,l,k+1,ct2,data)
				    +edangle5(j,i,j-1,ct1,data)
				    +edangle3(i,j,i+1,ct1,data)+3*gap);
													
			}
                                                
			if(aa+1<2*maxsep+2) {
			  //case 19 - j and m stacked
			  en1 = min(en1,w[c][i+1][a][e][aa+1][g]
				    +w[j-2][c+1][e][b+1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				    +2*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				    +edangle3(m,n,m+1,ct3,data)+4*gap);
													
			}
                                                
			//case 17 - stack on j
			en1 = min(en1,w[c][i+1][a][e][aa][g]
				  +w[j-2][c+1][e][b+1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				  +data->eparam[6]+edangle5(j,i,j-1,ct1,data)+2*gap);
												
                                                
			//case 59 - i, j, k and m stacked
			en1 = min(en1,w[c][i+2][a][e][aa][g]
				  +w[j-2][c+1][e][b+1][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				  +4*data->eparam[6]+edangle3(k,l,k+1,ct2,data)
				  +edangle5(j,i,j-1,ct1,data)+edangle3(m,n,m+1,ct3,data)
				  +edangle3(i,j,i+1,ct1,data)+2*gap);
												
		      }
                                            
		      if(a-1>0) {
                                                
			if(bb-1>=0) {
                                                    
			  if(aa+1<2*maxsep+2) {
                                                        
			    //case 36 - i, m and n stacked
			    en1 = min(en1,w[c][i+2][a-1][e][aa+1][g]
				      +w[j-1][c+1][e][b][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
				      +3*data->eparam[6]+edangle3(i,j,i+1,ct1,data)
				      +edangle5(n,m,n-1,ct3,data)+edangle3(m,n,m+1,ct3,data)+3*gap);
														
			  }
                                                    
			  //case 34 - i and n stacked
			  en1 = min(en1,w[c][i+2][a-1][e][aa][g]
				    +w[j-1][c+1][e][b][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
				    +2*data->eparam[6]+edangle3(i,j,i+1,ct1,data)
				    +edangle5(n,m,n-1,ct3,data)+4*gap);
													
			}
                                                
			if(bb+1<2*maxsep+2) {
                                                    
			  if(aa+1<2*maxsep+2) {
			    //case 55 - i, j, l and m stacked
			    en1 = min(en1,w[c][i+2][a-1][e][aa+1][g]
				      +w[j-2][c+1][e][b][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
				      +4*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
				      +edangle5(j,i,j-1,ct1,data)+edangle5(m,n,m+1,ct3,data)
				      +edangle3(i,j,i+1,ct1,data)+2*gap);
														
			  }
                                                    
			  //case 53 - i, j and l stacked
			  en1 = min(en1,w[c][i+2][a-1][e][aa][g]
				    +w[j-2][c+1][e][b][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
				    +3*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
				    +edangle5(j,i,j-1,ct1,data)
				    +edangle3(i,j,i+1,ct1,data)+3*gap);
													
			}
                                                
			if(aa+1<2*maxsep+2) {
			  //case 56 - i, j, l, m and n stacked
			  en1 = min(en1,w[c][i+2][a-1][e][aa+1][g]
				    +w[j-2][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				    +5*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
				    +edangle5(j,i,j-1,ct1,data)+edangle5(m,n,m+1,ct3,data)
				    +edangle3(i,j,i+1,ct1,data)+edangle3(n,m,n-1,ct3,data)+gap);
													
                                                    
			  //case 35 - i and m stacked
			  en1 = min(en1,w[c][i+2][a-1][e][aa+1][g]
				    +w[j-1][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				    +2*data->eparam[6]+edangle3(i,j,i+1,ct1,data)
				    +edangle3(m,n,m+1,ct3,data)+gap);
													
			}
                                                
			//case 33 - stack on i
			en1 = min(en1,w[c][i+2][a-1][e][aa][g]
				  +w[j-1][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				  +data->eparam[6]+edangle3(i,j,i+1,ct1,data)+2*gap);
												
                                                
			//case 54 - i, j, l and n stacked
			en1 = min(en1,w[c][i+2][a-1][e][aa][g]
				  +w[j-2][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				  +4*data->eparam[6]+edangle5(l,k,l-1,ct2,data)
				  +edangle5(j,i,j-1,ct1,data)+edangle5(n,m,n-1,ct3,data)
				  +edangle3(i,j,i+1,ct1,data)+2*gap);
												
		      }
                                            
		      if(a+1<2*maxsep+2) {
                                                
			if(bb-1>=0) {
                                                    
			  if(aa-1>0) {
			    //case 10 - k and n stacked
			    en1 = min(en1,w[c][i+1][a+1][e][aa-1][g]+3*data->eparam[10]
				      +w[j-1][c+1][e][b][g][bb-1]+3*data->eparam[5]
				      +2*data->eparam[6]+edangle3(k,l,k+1,ct2,data)
				      +edangle5(n,m,n-1,ct3,data)+4*gap);
														
			  }
                                                    
			  //case 12 - k, m and n stacked
			  en1 = min(en1,w[c][i+1][a+1][e][aa][g]+3*data->eparam[10]
				    +w[j-1][c+1][e][b][g][bb-1]+3*data->eparam[5]
				    +3*data->eparam[6]+edangle3(k,l,k+1,ct2,data)
				    +edangle5(n,m,n-1,ct3,data)
				    +edangle3(m,n,m+1,ct3,data)+3*gap);
													
			}
                                                
			if(bb+1<2*maxsep+2) {
                                                    
			  if(aa-1>0) {
			    //case 29 - j, k and l stacked
			    en1 = min(en1,w[c][i+1][a+1][e][aa-1][g]
				      +w[j-2][c+1][e][b][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
				      +3*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				      +edangle5(l,k,l-1,ct2,data)+edangle3(k,l,k+1,ct2,data)+3*gap);
														
			  }
                                                    
			  //case 31 - j, k, l and m stacked
			  en1 = min(en1,w[c][i+1][a+1][e][aa][g]
				    +w[j-2][c+1][e][b][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
				    +edangle5(l,k,l-1,ct2,data)+edangle3(k,l,k+1,ct2,data)
				    +edangle3(m,n,m+1,ct3,data)+2*gap);
													
			}
                                                
			if(aa-1>0) {
			  //case 9 - stack on k
			  en1 = min(en1,w[c][i+1][a+1][e][aa-1][g]+3*data->eparam[10]
				    +w[j-1][c+1][e][b][g][bb]+3*data->eparam[5]
				    +data->eparam[6]+edangle3(k,l,k+1,ct2,data)+2*gap);
													
                                                    
			  //case 30 - j, k, l, and n stacked
			  en1 = min(en1,w[c][i+1][a+1][e][aa-1][g]
				    +w[j-2][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				    +4*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				    +edangle5(l,k,l-1,ct2,data)+edangle3(k,l,k+1,ct2,data)
				    +edangle5(n,m,n-1,ct3,data)+2*gap);
													
			}
                                                
			//case 11 - k and m stacked
			en1 = min(en1,w[c][i+1][a+1][e][aa][g]+3*data->eparam[10]
				  +w[j-1][c+1][e][b][g][bb]+3*data->eparam[5]
				  +2*data->eparam[6]+edangle3(k,l,k+1,ct2,data)
				  +edangle3(m,n,m+1,ct3,data)+gap);
												
                                                
			//case 32 - j, k, l, m, and n stacked
			en1 = min(en1,w[c][i+1][a+1][e][aa][g]
				  +w[j-2][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				  +edangle5(l,k,l-1,ct2,data)+edangle3(k,l,k+1,ct2,data)
				  +edangle5(n,m,n-1,ct3,data)+edangle3(m,n,m+1,ct3,data)+gap);
												
		      }
                                            
		      if(bb-1>=0) {
                                                
			if(aa-1>0) {
			  //case 42 - i, k and n stacked
			  en1 = min(en1,w[c][i+2][a][e][aa-1][g]
				    +w[j-1][c+1][e][b][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
				    +3*data->eparam[6]+edangle3(i,j,i+1,ct1,data)
				    +edangle5(n,m,n-1,ct3,data)+edangle3(k,l,k+1,ct2,data)+3*gap);
													
			}
                                                
			if(aa+1<2*maxsep+2) {
			  //case 4 - m and n stacked
			  en1 = min(en1,w[c][i+1][a][e][aa+1][g]
				    +w[j-1][c+1][e][b][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
				    +2*data->eparam[6]+edangle5(n,m,n-1,ct3,data)
				    +edangle3(m,n,m+1,ct3,data)+4*gap);
													
                                                    
			}
                                                
			//case 2 - stack on n
			en1 = min(en1,w[c][i+1][a][e][aa][g]
				  +w[j-1][c+1][e][b][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
				  +data->eparam[6]+edangle5(n,m,n-1,ct3,data)+2*gap);
												
                                                
			//case 44 - i, k, m and n stacked
			en1 = min(en1,w[c][i+2][a][e][aa][g]
				  +w[j-1][c+1][e][b][g][bb-1]+3*data->eparam[5]+3*data->eparam[10]
				  +4*data->eparam[6]+edangle3(i,j,i+1,ct1,data)
				  +edangle3(k,l,k+1,ct2,data)+edangle3(m,n,m+1,ct3,data)
				  +edangle5(n,m,n-1,ct3,data)+2*gap);
												
		      }
                                            
		      if(bb+1<2*maxsep+2) {
                                                
			if(aa-1>0) {
			  //case 61 - i, j, k and l stacked
			  en1 = min(en1,w[c][i+2][a][e][aa-1][g]
				    +w[j-2][c+1][e][b][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
				    +4*data->eparam[6]
				    +edangle5(j,i,j-1,ct1,data)
				    +edangle3(i,j,i+1,ct1,data)
				    +edangle5(l,k,l-1,ct2,data)
				    +edangle3(k,l,k+1,ct2,data)+2*gap);
													
			}
                                                
			//case 63 - i, j, k, l and m stacked
			en1 = min(en1,w[c][i+2][a][e][aa][g]
				  +w[j-2][c+1][e][b][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
				  +5*data->eparam[6]
				  +edangle5(j,i,j-1,ct1,data)
				  +edangle3(i,j,i+1,ct1,data)
				  +edangle5(l,k,l-1,ct2,data)
				  +edangle3(k,l,k+1,ct2,data)
				  +edangle3(m,n,m+1,ct3,data)+gap);
												
                                                
			//case 21 - j and l stacked
			en1 = min(en1,w[c][i+1][a][e][aa][g]
				  +w[j-2][c+1][e][b][g][bb+1 ]+3*data->eparam[5]+3*data->eparam[10]
				  +2*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				  +edangle5(l,k,l-1,ct2,data)+gap);
												
		      }
                                            
		      if(aa-1>0) {
			//case 41 - i and k stacked
			en1 = min(en1,w[c][i+2][a][e][aa-1][g]
				  +w[j-1][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				  +2*data->eparam[6]+edangle3(k,l,k+1,ct2,data)+
				  edangle3(i,j,i+1,ct1,data)+gap);
												
                                                
			//case 62 - i, j, k, l and n stacked
			en1 = min(en1,w[c][i+2][a][e][aa-1][g]
				  +w[j-2][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				  +5*data->eparam[6]
				  +edangle5(j,i,j-1,ct1,data)
				  +edangle3(i,j,i+1,ct1,data)
				  +edangle5(l,k,l-1,ct2,data)
				  +edangle3(k,l,k+1,ct2,data)
				  +edangle5(n,m,n-1,ct3,data)+gap);
												
		      }
                                            
		      if(aa+1<2*maxsep+2) {
			//case 3 - stack on m
			en1 = min(en1,w[c][i+1][a][e][aa+1][g]+3*data->eparam[10]
				  +w[j-1][c+1][e][b][g][bb]+3*data->eparam[5]
				  +data->eparam[6]+edangle3(m,n,m+1,ct3,data)+2*gap);
												
                                                
			//case 24 - j, l, m and n stacked
			en1 = min(en1,w[c][i+1][a][e][aa+1][g]
				  +w[j-2][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				  +4*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				  +edangle5(l,k,l-1,ct2,data)
				  +edangle5(n,m,n-1,ct3,data)
				  +edangle3(m,n,m+1,ct3,data)+2*gap);
												
			if(bb+1<2*maxsep+2) {
			  //case 23 - j and l and m stacked
			  en1 = min(en1,w[c][i+1][a][e][aa+1][g]
				    +w[j-2][c+1][e][b][g][bb+1]+3*data->eparam[5]+3*data->eparam[10]
				    +3*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				    +edangle5(l,k,l-1,ct2,data)
				    +edangle3(m,n,m+1,ct2,data)+3*gap);
													
			}
		      }
                                            
		      //case 1 - no stacks
		      en1 = min(en1,w[c][i+1][a][e][aa][g]
				+w[j-1][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]);
											
                                            
		      //case 22 - j, l and n stacked
		      en1 = min(en1,w[c][i+1][a][e][aa][g]
				+w[j-2][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				+3*data->eparam[6]+edangle5(j,i,j-1,ct1,data)
				+edangle5(l,k,l-1,ct2,data)
				+edangle5(n,m,n-1,ct2,data));
											
                                            
		      //case 43 - i, k and m stacked
		      en1 = min(en1,w[c][i+2][a][e][aa][g]
				+w[j-1][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				+3*data->eparam[6]+edangle3(i,j,i+1,ct1,data)
				+edangle3(k,l,k+1,ct2,data)+edangle3(m,n,m+1,ct3,data));
											
                                            
		      //case 64 - all six stacked
		      en1 = min(en1,w[c][i+2][a][e][aa][g]
				+w[j-2][c+1][e][b][g][bb]+3*data->eparam[5]+3*data->eparam[10]
				+6*data->eparam[6]
				+edangle5(n,m,n-1,ct3,data)
				+edangle3(m,n,m+1,ct3,data)
				+edangle5(l,k,l-1,ct2,data)
				+edangle3(k,l,k+1,ct2,data)
				+edangle5(j,i,j-1,ct1,data)
				+edangle3(i,j,i+1,ct1,data));
											
		    }
		  }
		}
                                
		en1 = en1 + penalty(i,j,ct1,data) + penalty(k,l,ct2,data) + penalty(m,n,ct3,data);
		v[j][i][a][b][aa][bb]=min(en1,v[j][i][a][b][aa][bb]);
                                
                                
                                
	      }
	      else v[j][i][a][b][aa][bb] = infinity;	
	      /************************************** end fill V *************************************/
                            
	      //Fill w, the best energy for a fragment in a multiloop:
                            
	      //save a temporaray value for w in en1, a register short int
                            
	      /******************************************** W1 & W2 at the same time *********************************************/
	      //consider the possibilities for adding nucleotides to existing fragments
                            
	      //		i	j	k	l	m	n
                            
	      //2		0	0	0	0	0	1
	      //3		0	0	0	0	1	0
	      //4		0	0	0	0	1	1
	      //5		0	0	0	1	0	0
	      //6		0	0	0	1	0	1
	      //7		0	0	0	1	1	0
	      //8		0	0	0	1	1	1
	      //9		0	0	1	0	0	0
	      //10	0	0	1	0	0	1
	      //11	0	0	1	0	1	0
	      //12	0	0	1	0	1	1
	      //13	0	0	1	1	0	0
	      //14	0	0	1	1	0	1
	      //15	0	0	1	1	1	0
	      //16	0	0	1	1	1	1
	      //17	0	1	0	0	0	0
	      //18	0	1	0	0	0	1
	      //19	0	1	0	0	1	0
	      //20	0	1	0	0	1	1
	      //21	0	1	0	1	0	0
	      //22	0	1	0	1	0	1
	      //23	0	1	0	1	1	0
	      //24	0	1	0	1	1	1
	      //25	0	1	1	0	0	0
	      //26	0	1	1	0	0	1
	      //27	0	1	1	0	1	0
	      //28	0	1	1	0	1	1
	      //29	0	1	1	1	0	0
	      //30	0	1	1	1	0	1
	      //31	0	1	1	1	1	0
	      //32	0	1	1	1	1	1
	      //33	1	0	0	0	0	0
	      //34	1	0	0	0	0	1
	      //35	1	0	0	0	1	0
	      //36	1	0	0	0	1	1
	      //37	1	0	0	1	0	0
	      //38	1	0	0	1	0	1
	      //39	1	0	0	1	1	0
	      //40	1	0	0	1	1	1
	      //41	1	0	1	0	0	0
	      //42	1	0	1	0	0	1
	      //43	1	0	1	0	1	0
	      //44	1	0	1	0	1	1
	      //45	1	0	1	1	0	0
	      //46	1	0	1	1	0	1
	      //47	1	0	1	1	1	0
	      //48	1	0	1	1	1	1
	      //49	1	1	0	0	0	0
	      //50	1	1	0	0	0	1
	      //51	1	1	0	0	1	0
	      //52	1	1	0	0	1	1
	      //53	1	1	0	1	0	0
	      //54	1	1	0	1	0	1
	      //55	1	1	0	1	1	0
	      //56	1	1	0	1	1	1
	      //57	1	1	1	0	0	0
	      //58	1	1	1	0	0	1
	      //59	1	1	1	0	1	0
	      //60	1	1	1	0	1	1
	      //61	1	1	1	1	0	0
	      //62	1	1	1	1	0	1
	      //63	1	1	1	1	1	0
	      //64	1	1	1	1	1	1
                            

	      // w1 Does not have case 1
	      //case 1 - nothing stacked
	      en1=v[j][i][a][b][aa][bb]+3*data->eparam[10]
		+penalty(i,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n,ct3,data);


	      //case 22 j,l,n
	      en1 = min(en1,w[j-1][i][a][b][aa][bb]+3*data->eparam[6]);

	      //case 22 - j, l and n stacked
	      en1 = min(en1,v[j-1][i][a][b][aa][bb]+3*data->eparam[10]
			+penalty(i,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n-1,ct3,data)
			+edangle3(l-1,k,l,ct2,data)+edangle3(j-1,i,j,ct1,data)+edangle3(n-1,m,n,ct3,data)
			+3*data->eparam[6]);

							                            
	      //case 43 i,k,m
	      en1 = min(en1,w[j][i+1][a][b][aa][bb]+3*data->eparam[6]);

	      //case 43 - i, k and m stacked
	      en1 = min(en1,v[j][i+1][a][b][aa][bb]+3*data->eparam[10]
			+penalty(i+1,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n,ct3,data)
			+edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle5(m+1,n,m,ct3,data)
			+3*data->eparam[6]);

                            
	      if(a>=1) {
                                
		if(j-1>i+1) {
                                    
		  if(b+1<2*maxseparation+2) {
                                        
		    if(aa+1<2*maxseparation+2) {
                                            
		      if(bb>=1) {
			//case 52 i,j,m,n	
			en1 = min(en1,w[j-1][i+1][a-1][b+1][aa+1][bb-1]+4*data->eparam[6]+2*gap);

			//case 52 - i, j, m and n stacked
			en1=min(en1,v[j-1][i+1][a-1][b+1][aa+1][bb-1]+3*data->eparam[10]
				+penalty(i+1,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n-1,ct3,data)
				+edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle5(m+1,n-1,m,ct3,data)
				+edangle3(n-1,m+1,n,ct3,data)
				+4*data->eparam[6]+2*gap);

		      }
                                            
		      //case 51 i,j,m
		      en1 = min(en1,w[j-1][i+1][a-1][b+1][aa+1][bb]+3*data->eparam[6]+3*gap);

		      //case 51 - i, j and m stacked
		      en1=min(en1,v[j-1][i+1][a-1][b+1][aa+1][bb]+3*data->eparam[10]
			      +penalty(i+1,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n,ct3,data)
			      +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle5(m+1,n,m,ct3,data)
			      +3*data->eparam[6]+3*gap);

		    }
                                        
		    if(bb>=1) {
		      //case 50 i,j,n
		      en1 = min(en1,w[j-1][i+1][a-1][b+1][aa][bb-1]+3*data->eparam[6]+3*gap);

		      //case 50 - i, j and n stacked
		      en1=min(en1,v[j-1][i+1][a-1][b+1][aa][bb-1]+3*data->eparam[10]
			      +penalty(i+1,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n-1,ct3,data)
			      +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle3(n-1,m,n,ct3,data)
			      +3*data->eparam[6]+3*gap);

		    }
                                        
		    //case 49 i,j
		    en1 = min(en1,w[j-1][i+1][a-1][b+1][aa][bb]+2*data->eparam[6]+4*gap);

		    //case 49 - i and j stacked
		    en1=min(en1,v[j-1][i+1][a-1][b+1][aa][bb]+3*data->eparam[10]
			    +penalty(i+1,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n,ct3,data)
			    +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)
			    +2*data->eparam[6]+4*gap);

		  }
                                    
		  if(aa+1<2*maxseparation+2) {
                                        
		    if(bb+1<2*maxseparation+2) {
		      //case 55 i,j,l,m
		      en1 = min(en1,w[j-1][i+1][a-1][b][aa+1][bb+1]+4*data->eparam[6]+2*gap);

		      //case 55 - i, j, l and m stacked
		      en1=min(en1,v[j-1][i+1][a-1][b][aa+1][bb+1]+3*data->eparam[10]
			      +penalty(i+1,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n,ct3,data)
			      +edangle3(j-1,i+1,j,ct1,data)+edangle5(i+1,j-1,i,ct1,data)+edangle3(l-1,k,l,ct2,data)
			      +edangle5(m+1,n,m,ct3,data)
			      +4*data->eparam[6]+2*gap);

		    }
                                        
		    //case 56 i,j,l,m,n
		    en1 = min(en1,w[j-1][i+1][a-1][b][aa+1][bb]+5*data->eparam[6]+gap);

		    //case 56 - i, j, l, m and n stacked
		    en1=min(en1,v[j-1][i+1][a-1][b][aa+1][bb]+3*data->eparam[10]
			    +penalty(i+1,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
			    +edangle3(j-1,i+1,j,ct1,data)+edangle5(i+1,j-1,i,ct1,data)+edangle3(l-1,k,l,ct2,data)
			    +edangle3(n-1,m+1,n,ct3,data)+edangle5(m+1,n-1,m,ct3,data)
			    +5*data->eparam[6]+gap);

		  }
                                    
		  if(bb+1<2*maxseparation+2) {
		    //case 53 i,j,l
		    en1 = min(en1,w[j-1][i+1][a-1][b][aa][bb+1]+3*data->eparam[6]+3*gap);

		    //case 53 - i, j and l stacked
		    en1=min(en1,v[j-1][i+1][a-1][b][aa][bb+1]+3*data->eparam[10]
			    +penalty(i+1,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n,ct3,data)
			    +edangle3(j-1,i+1,j,ct1,data)+edangle5(i+1,j-1,i,ct1,data)+edangle3(l-1,k,l,ct2,data)
			    +3*data->eparam[6]+3*gap);

		  }
                                    
		  //case 54 i,j,l,n
		  en1 = min(en1,w[j-1][i+1][a-1][b][aa][bb]+4*data->eparam[6]+2*gap);

		  //case 54 - i, j, l and n stacked
		  en1=min(en1,v[j-1][i+1][a-1][b][aa][bb]+3*data->eparam[10]
			  +penalty(i+1,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n-1,ct3,data)
			  +edangle3(j-1,i+1,j,ct1,data)+edangle5(i+1,j-1,i,ct1,data)+edangle3(l-1,k,l,ct2,data)
			  +edangle3(n-1,m,n,ct3,data)
			  +4*data->eparam[6]+2*gap);

		}
                                
		if(b>=1) {
                                    
		  if(aa+1<2*maxseparation+2) {
                                        
		    if(bb+1<2*maxseparation+2) {
		      //case 39 i,l,m
		      en1 = min(en1,w[j][i+1][a-1][b-1][aa+1][bb+1]+3*data->eparam[6]+3*gap);

		      //case 39 - i, l and m stacked
		      en1=min(en1,v[j][i+1][a-1][b-1][aa+1][bb+1]+3*data->eparam[10]
			      +penalty(i+1,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n,ct3,data)
			      +edangle5(i+1,j,i,ct1,data)+edangle3(l-1,k,l,ct2,data)+edangle5(m+1,n,m,ct3,data)
			      +3*data->eparam[6]+3*gap);

		    }
                                        
		    //case 40 i,l,m,n
		    en1 = min(en1,w[j][i+1][a-1][b-1][aa+1][bb]+4*data->eparam[6]+2*gap);

		    //case 40 - i, l, m and n stacked
		    en1=min(en1,v[j][i+1][a-1][b-1][aa+1][bb]+3*data->eparam[10]
			    +penalty(i+1,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
			    +edangle5(i+1,j,i,ct1,data)+edangle3(l-1,k,l,ct2,data)+edangle5(m+1,n-1,m,ct3,data)
			    +edangle3(n-1,m+1,n,ct3,data)
			    +4*data->eparam[6]+2*gap);

		  }
                                    
		  if(bb+1<2*maxseparation+2) {
		    //case 37 i,l
		    en1 = min(en1,w[j][i+1][a-1][b-1][aa][bb+1]+2*data->eparam[6]+4*gap);

		    //case 37 - i and l stacked
		    en1=min(en1,v[j][i+1][a-1][b-1][aa][bb+1]+3*data->eparam[10]
			    +penalty(i+1,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n,ct3,data)
			    +edangle5(i+1,j,i,ct1,data)+edangle3(l-1,k,l,ct2,data)
			    +2*data->eparam[6]+4*gap);

		  }
                                    
		  //case 38 i,l,n
		  en1 = min(en1,w[j][i+1][a-1][b-1][aa][bb]+3*data->eparam[6]+3*gap);

		  //case 38 - i, l and n stacked
		  en1=min(en1,v[j][i+1][a-1][b-1][aa][bb]+3*data->eparam[10]
			  +penalty(i+1,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n-1,ct3,data)
			  +edangle5(i+1,j,i,ct1,data)+edangle3(l-1,k,l,ct2,data)+edangle3(n-1,m,n,ct3,data)
			  +3*data->eparam[6]+3*gap);

		}
                                
		if(aa+1<2*maxseparation+2) {
                                    
		  if(bb>=1) {
		    //case 36 i,m,n
		    en1 = min(en1,w[j][i+1][a-1][b][aa+1][bb-1]+3*data->eparam[6]+3*gap);

		    //case 36 - i, m and n stacked
		    en1 = min(en1,v[j][i+1][a-1][b][aa+1][bb-1]+3*data->eparam[10]
			      +penalty(i+1,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n-1,ct3,data)
			      +edangle5(i+1,j,i,ct1,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
			      +3*data->eparam[6]+3*gap);

		  }
                                    
		  //case 35 i,m
		  en1 = min(en1,w[j][i+1][a-1][b][aa+1][bb]+2*data->eparam[6]+gap);

		  //case 35 - i and m stacked
		  en1 = min(en1,v[j][i+1][a-1][b][aa+1][bb]+3*data->eparam[10]
			    +penalty(i+1,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n,ct3,data)
			    +edangle5(i+1,j,i,ct1,data)+edangle5(m+1,n,m,ct3,data)
			    +2*data->eparam[6]+gap);

		}
                                
		if(bb>=1) {
		  //case 34 i,n
		  en1 = min(en1,w[j][i+1][a-1][b][aa][bb-1]+2*data->eparam[6]+4*gap);

		  //case 34 - i and n stacked
		  en1 = min(en1,v[j][i+1][a-1][b][aa][bb-1]+3*data->eparam[10]
			    +penalty(i+1,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n-1,ct3,data)
			    +edangle5(i+1,j,i,ct1,data)+edangle3(n-1,m,n,ct3,data)
			    +2*data->eparam[6]+4*gap);
									
		}
                                
		//case 33 i
		en1 = min(en1,w[j][i+1][a-1][b][aa][bb]+data->eparam[6]+2*gap);

		//case 33 - stack on i
		en1 = min(en1,v[j][i+1][a-1][b][aa][bb]+3*data->eparam[10]
			  +penalty(i+1,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n,ct3,data)
			  +edangle5(i+1,j,i,ct1,data)
			  +data->eparam[6]+2*gap);

	      }
                            
	      if(a+1<2*maxseparation+2) {
                                
		if(b>=1) {
                                    
		  if(aa>=1) {
                                        
		    if(bb+1<2*maxseparation+2) {
		      //case 13 k,l
		      en1 = min(en1,w[j][i][a+1][b-1][aa-1][bb+1]+2*data->eparam[6]+4*gap);

		      //case 13 - k and l stacked
		      en1=min(en1,v[j][i][a+1][b-1][aa-1][bb+1]+3*data->eparam[10]
			      +penalty(i,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n,ct3,data)
			      +edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)
			      +2*data->eparam[6]+4*gap);

		    }
                                        
		    //case 14 k,l,n
		    en1 = min(en1,w[j][i][a+1][b-1][aa-1][bb]+3*data->eparam[6]+3*gap);

		    //case 14 - k, l and n stacked
		    en1=min(en1,v[j][i][a+1][b-1][aa-1][bb]+3*data->eparam[10]
			    +penalty(i,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n-1,ct3,data)
			    +edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle3(n-1,m,n,ct3,data)
			    +3*data->eparam[6]+3*gap);

		  }
                                    
		  if(bb+1<2*maxseparation+2) {
		    //case 15 k,l,m
		    en1 = min(en1,w[j][i][a+1][b-1][aa][bb+1]+3*data->eparam[6]+3*gap);

		    //case 15 - k, l and m stacked
		    en1=min(en1,v[j][i][a+1][b-1][aa][bb+1]+3*data->eparam[10]
			    +penalty(i,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n,ct3,data)
			    +edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle5(m+1,n,m,ct3,data)
			    +3*data->eparam[6]+3*gap);

		  }
                                    
		  //case 16 k,l,m,n
		  en1 = min(en1,w[j][i][a+1][b-1][aa][bb]+4*data->eparam[6]+2*gap);

		  //case 16 - k, l, m and n stacked
		  en1=min(en1,v[j][i][a+1][b-1][aa][bb]+3*data->eparam[10]
			  +penalty(i,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
			  +edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle5(m+1,n-1,m,ct3,data)
			  +edangle3(n-1,m+1,n,ct3,data)
			  +4*data->eparam[6]+2*gap);

		}
                                
		if(b+1<2*maxseparation+2) {
                                    
		  if(aa>=1) {
                                        
		    if(bb>=1) {
		      //case 26 j,k,n
		      en1 = min(en1,w[j-1][i][a+1][b+1][aa-1][bb-1]+3*data->eparam[6]+3*gap);

		      //case 26 - j, k and n stacked
		      en1 = min(en1,v[j-1][i][a+1][b+1][aa-1][bb-1]+3*data->eparam[10]
				+penalty(i,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n-1,ct3,data)
				+edangle5(k+1,l,k,ct2,data)+edangle3(j-1,i,j,ct1,data)+edangle3(n-1,m,n,ct3,data)
				+3*data->eparam[6]+3*gap);

		    }
                                        
		    //case 25 j,k
		    en1 = min(en1,w[j-1][i][a+1][b+1][aa-1][bb]+2*data->eparam[6]+4*gap);

		    //case 25 - j and k stacked
		    en1 = min(en1,v[j-1][i][a+1][b+1][aa-1][bb]+3*data->eparam[10]
			      +penalty(i,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n,ct3,data)
			      +edangle5(k+1,l,k,ct2,data)+edangle3(j-1,i,j,ct1,data)
			      +2*data->eparam[6]+4*gap);

		  }
                                    
		  if(bb>=1) {
		    //case 28 j,k,m,n
		    en1 = min(en1,w[j-1][i][a+1][b+1][aa][bb-1]+4*data->eparam[6]+2*gap);

		    //case 28 - j, k, m and n stacked
		    en1 = min(en1,v[j-1][i][a+1][b+1][aa][bb-1]+3*data->eparam[10]
			      +penalty(i,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n-1,ct3,data)
			      +edangle5(k+1,l,k,ct2,data)+edangle3(j-1,i,j,ct1,data)+edangle5(m+1,n-1,m,ct3,data)
			      +edangle3(n-1,m+1,n,ct3,data)
			      +4*data->eparam[6]+2*gap);

		  }
                                    
		  //case 27 j,k,m
		  en1 = min(en1,w[j-1][i][a+1][b+1][aa][bb]+3*data->eparam[6]+3*gap);

		  //case 27 - j, k and m stacked
		  en1 = min(en1,v[j-1][i][a+1][b+1][aa][bb]+3*data->eparam[10]
			    +penalty(i,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n,ct3,data)
			    +edangle5(k+1,l,k,ct2,data)+edangle3(j-1,i,j,ct1,data)+edangle5(m+1,n,m,ct3,data)
			    +3*data->eparam[6]+3*gap);

		}
                                
		if(aa>=1) {
                                    
		  if(bb>=1) {
		    //case 10 k,n
		    en1 = min(en1,w[j][i][a+1][b][aa-1][bb-1]+2*data->eparam[6]+4*gap);

		    //case 10 - k and n stacked
		    en1=min(en1,v[j][i][a+1][b][aa-1][bb-1]+3*data->eparam[10]
			    +penalty(i,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n-1,ct3,data)
			    +edangle5(k+1,l,k,ct2,data)+edangle3(n-1,m,n,ct3,data)
			    +2*data->eparam[6]+4*gap);

		  }
                                    
		  if(bb+1<2*maxseparation+2) {
		    //case 29 j,k,l
		    en1 = min(en1,w[j-1][i][a+1][b][aa-1][bb+1]+3*data->eparam[6]+3*gap);

		    //case 29 - j, k and l stacked
		    en1 = min(en1,v[j-1][i][a+1][b][aa-1][bb+1]+3*data->eparam[10]
			      +penalty(i,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n,ct3,data)
			      +edangle3(l-1,k+1,l,ct2,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(j-1,i,j,ct1,data)
			      +3*data->eparam[6]+3*gap);

		  }
                                    
		  //case 9 k
		  en1 = min(en1,w[j][i][a+1][b][aa-1][bb]+data->eparam[6]+2*gap);

		  //case 9 - stack on k
		  en1=min(en1,v[j][i][a+1][b][aa-1][bb]+3*data->eparam[10]
			  +penalty(i,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n,ct3,data)
			  +edangle5(k+1,l,k,ct2,data)
			  +data->eparam[6]+2*gap);

                                    
		  //case 30 j,k,l,n
		  en1 = min(en1,w[j-1][i][a+1][b][aa-1][bb]+4*data->eparam[6]+2*gap);

		  //case 30 - j, k, l and n stacked
		  en1 = min(en1,v[j-1][i][a+1][b][aa-1][bb]+3*data->eparam[10]
			    +penalty(i,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n-1,ct3,data)
			    +edangle3(l-1,k+1,l,ct2,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(j-1,i,j,ct1,data)
			    +edangle3(n-1,m,n,ct3,data)
			    +4*data->eparam[6]+2*gap);

		}
                                
		if(bb>=1) {
		  //case 12 k,m,n
		  en1 = min(en1,w[j][i][a+1][b][aa][bb-1]+3*data->eparam[6]+3*gap);

		  //case 12 - k, m and n stacked
		  en1=min(en1,v[j][i][a+1][b][aa][bb-1]+3*data->eparam[10]
			  +penalty(i,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n-1,ct3,data)
			  +edangle5(k+1,l,k,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,m,ct3,data)
			  +3*data->eparam[6]+3*gap);

		}
                                
		if(bb+1<2*maxseparation+2) {
		  //case 31 j,k,l,m
		  en1 = min(en1,w[j-1][i][a+1][b][aa][bb+1]+4*data->eparam[6]+2*gap);

		  //case 31 - j, k, l and m stacked
		  en1 = min(en1,v[j-1][i][a+1][b][aa][bb+1]+3*data->eparam[10]
			    +penalty(i,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n,ct3,data)
			    +edangle3(l-1,k+1,l,ct2,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(j-1,i,j,ct1,data)
			    +edangle5(m+1,n,m,ct3,data)
			    +4*data->eparam[6]+2*gap);

		}
                                
		//case 11 k,m
		en1 = min(en1,w[j][i][a+1][b][aa][bb]+2*data->eparam[6]+gap);

		//case 11 - k and m stacked
		en1=min(en1,v[j][i][a+1][b][aa][bb]+3*data->eparam[10]
			+penalty(i,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n,ct3,data)
			+edangle5(k+1,l,k,ct2,data)+edangle5(m+1,n,m,ct3,data)
			+2*data->eparam[6]+gap);

                                
		//case 32 j,k,l,m,n
		en1 = min(en1,w[j-1][i][a+1][b][aa][bb]+5*data->eparam[6]+gap);

		//case 32 - j, k, l, m and n stacked
		en1 = min(en1,v[j-1][i][a+1][b][aa][bb]+3*data->eparam[10]
			  +penalty(i,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
			  +edangle3(l-1,k+1,l,ct2,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(j-1,i,j,ct1,data)
			  +edangle3(n-1,m+1,n,ct3,data)+edangle5(m+1,n-1,m,ct3,data)
			  +5*data->eparam[6]+gap);

	      }
                            
	      if(j-1>i+1) {
                                
		if(b+1<2*maxseparation+2) {
                                    
		  if(aa>=1) {
                                        
		    if(bb>=1) {
		      //case 58 i,j,k,n
		      en1 = min(en1,w[j-1][i+1][a][b+1][aa-1][bb-1]+4*data->eparam[6]+2*gap);

		      //case 58 - i, j, k and n stacked
		      en1 = min(en1,v[j-1][i+1][a][b+1][aa-1][bb-1]+3*data->eparam[10]
				+penalty(i+1,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n-1,ct3,data)
				+edangle3(j-1,i+1,j,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle5(i+1,j-1,i,ct1,data)
				+edangle3(n-1,m,n,ct3,data)
				+4*data->eparam[6]+2*gap);

		    }
                                        
		    //case 57 i,j,k
		    en1 = min(en1,w[j-1][i+1][a][b+1][aa-1][bb]+3*data->eparam[6]+3*gap);

		    //case 57 - i, j and k stacked
		    en1 = min(en1,v[j-1][i+1][a][b+1][aa-1][bb]+3*data->eparam[10]
			      +penalty(i+1,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n,ct3,data)
			      +edangle3(j-1,i+1,j,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle5(i+1,j-1,i,ct1,data)
			      +3*data->eparam[6]+3*gap);

		  }
                                    
		  if(bb>=1) {
		    //case 60 i,j,k,m,n
		    en1 = min(en1,w[j-1][i+1][a][b+1][aa][bb-1]+5*data->eparam[6]+gap);

		    //case 60 - i, j, k, m and n stacked
		    en1 = min(en1,v[j-1][i+1][a][b+1][aa][bb-1]+3*data->eparam[10]
			      +penalty(i+1,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n-1,ct3,data)
			      +edangle3(j-1,i+1,j,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle5(i+1,j-1,i,ct1,data)
			      +edangle3(n-1,m+1,n,ct3,data)+edangle5(m+1,n-1,m,ct3,data)
			      +5*data->eparam[6]+gap);
										 
		  }
                                    
		  //case 59 i,j,k,m
		  en1 = min(en1,w[j-1][i+1][a][b+1][aa][bb]+4*data->eparam[6]+2*gap);

		  //case 59 - i, j, k and m stacked
		  en1 = min(en1,v[j-1][i+1][a][b+1][aa][bb]+3*data->eparam[10]
			    +penalty(i+1,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n,ct3,data)
			    +edangle3(j-1,i+1,j,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle5(i+1,j-1,i,ct1,data)
			    +edangle5(m+1,n,m,ct3,data)
			    +4*data->eparam[6]+2*gap);

		}
                                
		if(aa>=1) {
                                    
		  if(bb+1<2*maxseparation+2) {
		    //case 61 i,j,k,l
		    en1 = min(en1,w[j-1][i+1][a][b][aa-1][bb+1]+4*data->eparam[6]+2*gap);

		    //case 61 - i, j, k and l stacked
		    en1=min(en1,v[j-1][i+1][a][b][aa-1][bb+1]+3*data->eparam[10]
			    +penalty(i+1,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n,ct3,data)
			    +edangle3(l-1,k+1,l,ct2,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(j-1,i+1,j,ct1,data)
			    +edangle5(i+1,j-1,i,ct1,data)
			    +4*data->eparam[6]+2*gap);

		  }
                                    
		  //case 62 i,j,k,l,n
		  en1 = min(en1,w[j-1][i+1][a][b][aa-1][bb]+5*data->eparam[6]+gap);

		  //case 62 - i, j, k, l and n stacked
		  en1=min(en1,v[j-1][i+1][a][b][aa-1][bb]+3*data->eparam[10]
			  +penalty(i+1,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n-1,ct3,data)
			  +edangle3(l-1,k+1,l,ct2,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(j-1,i+1,j,ct1,data)
			  +edangle5(i+1,j-1,i,ct1,data)+edangle3(n-1,m,n,ct3,data)
			  +5*data->eparam[6]+gap);

		}
                                
		if(bb+1<2*maxseparation+2) {
		  //case 63 i,j,k,l,m
		  en1 = min(en1,w[j-1][i+1][a][b][aa][bb+1]+5*data->eparam[6]+gap);

		  //case 63 - i, j, k, l and m stacked
		  en1=min(en1,v[j-1][i+1][a][b][aa][bb+1]+3*data->eparam[10]
			  +penalty(i+1,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n,ct3,data)
			  +edangle3(l-1,k+1,l,ct2,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(j-1,i+1,j,ct1,data)
			  +edangle5(i+1,j-1,i,ct1,data)+edangle5(m+1,n,m,ct3,data)
			  +5*data->eparam[6]+gap);

		}
                                
		//case 64 i,j,k,l,m,n
		en1 = min(en1,w[j-1][i+1][a][b][aa][bb]+6*data->eparam[6]);

		//case 64 - all six stacked
		en1=min(en1,v[j-1][i+1][a][b][aa][bb]+3*data->eparam[10]
			+penalty(i+1,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
			+edangle3(l-1,k+1,l,ct2,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(j-1,i+1,j,ct1,data)
			+edangle5(i+1,j-1,i,ct1,data)+edangle3(n-1,m+1,n,ct3,data)+edangle5(m+1,n-1,m,ct3,data)
			+6*data->eparam[6]);

	      }
                            
	      if(b>=1) {
                                
		if(aa>=1) {
                                    
		  if(bb+1<2*maxseparation+2) {
		    //case 45 i,k,l
		    en1 = min(en1,w[j][i+1][a][b-1][aa-1][bb+1]+3*data->eparam[6]+3*gap);

		    //case 45 - i, k and l stacked
		    en1=min(en1,v[j][i+1][a][b-1][aa-1][bb+1]+3*data->eparam[10]
			    +penalty(i+1,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n,ct3,data)
			    +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)
			    +3*data->eparam[6]+3*gap);

		  }
                                    
		  //case 46 i,k,l,n
		  en1 = min(en1,w[j][i+1][a][b-1][aa-1][bb]+4*data->eparam[6]+2*gap);

		  //case 46 - i, k, l and n stacked
		  en1=min(en1,v[j][i+1][a][b-1][aa-1][bb]+3*data->eparam[10]
			  +penalty(i+1,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n-1,ct3,data)
			  +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)
			  +edangle3(n-1,m,n,ct3,data)
			  +4*data->eparam[6]+2*gap);

		}
                                
		if(aa+1<2*maxseparation+2) {
                                    
		  if(bb+1<2*maxseparation+2) {
		    //case 7 l,m
		    en1 = min(en1,w[j][i][a][b-1][aa+1][bb+1]+2*data->eparam[6]+4*gap);

		    //case 7 - l and m stacked
		    en1=min(en1,v[j][i][a][b-1][aa+1][bb+1]+3*data->eparam[10]
			    +penalty(i,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n,ct3,data)
			    +edangle3(l-1,k,l,ct2,data)+edangle5(m+1,n,m,ct3,data)
			    +2*data->eparam[6]+4*gap);

		  }
                                    
		  //case 8 l,m,n
		  en1 = min(en1,w[j][i][a][b-1][aa+1][bb]+3*data->eparam[6]+3*gap);

		  //case 8 - l, m and n stacked
		  en1=min(en1,v[j][i][a][b-1][aa+1][bb]+3*data->eparam[10]
			  +penalty(i,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
			  +edangle3(l-1,k,l,ct2,data)+edangle3(n-1,m+1,n,ct3,data)+edangle5(m+1,n-1,m,ct3,data)
			  +3*data->eparam[6]+3*gap);

		}
                                
		if(bb+1<2*maxseparation+2) {
		  //case 5 l
		  en1 = min(en1,w[j][i][a][b-1][aa][bb+1]+data->eparam[6]+2*gap);

		  //case 5 - stack on l
		  en1=min(en1,v[j][i][a][b-1][aa][bb+1]+3*data->eparam[10]
			  +penalty(i,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n,ct3,data)
			  +edangle3(l-1,k,l,ct2,data)
			  +data->eparam[6]+2*gap);

                                    
		  //case 47 i,k,l,m
		  en1 = min(en1,w[j][i+1][a][b-1][aa][bb+1]+4*data->eparam[6]+2*gap);

		  //case 47 - i, k, l and m stacked
		  en1=min(en1,v[j][i+1][a][b-1][aa][bb+1]+3*data->eparam[10]
			  +penalty(i+1,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n,ct3,data)
			  +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)
			  +edangle5(m+1,n,m,ct3,data)
			  +4*data->eparam[6]+2*gap);

		}
                                
		//case 6 l,n
		en1 = min(en1,w[j][i][a][b-1][aa][bb]+2*data->eparam[6]+gap);

		//case 6 - l and n stacked
		en1=min(en1,v[j][i][a][b-1][aa][bb]+3*data->eparam[10]
			+penalty(i,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n-1,ct3,data)
			+edangle3(l-1,k,l,ct2,data)+edangle3(n-1,m,n,ct3,data)
			+2*data->eparam[6]+gap);

                                
		//case 48 i,k,l,m,n
		en1 = min(en1,w[j][i+1][a][b-1][aa][bb]+5*data->eparam[6]+gap);

		//case 48 - i, k, l, m and n stacked
		en1=min(en1,v[j][i+1][a][b-1][aa][bb]+3*data->eparam[10]
			+penalty(i+1,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
			+edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)
			+edangle3(n-1,m+1,n,ct3,data)+edangle5(m+1,n-1,m,ct3,data)
			+5*data->eparam[6]+gap);

	      }
                            
	      if(b+1<2*maxseparation+2) {
                                
		if(aa+1<2*maxseparation+2) {
                                    
		  if(bb>=1) {
		    //case 20 j,m,n
		    en1 = min(en1,w[j-1][i][a][b+1][aa+1][bb-1]+3*data->eparam[6]+3*gap);

		    //case 20 - j, m and n stacked
		    en1 = min(en1,v[j-1][i][a][b+1][aa+1][bb-1]+3*data->eparam[10]
			      +penalty(i,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n-1,ct3,data)
			      +edangle3(j-1,i,j,ct1,data)+edangle3(n-1,m+1,n,ct3,data)+edangle5(m+1,n-1,m,ct3,data)
			      +3*data->eparam[6]+3*gap);

		  }
                                    
		  //case 19 j,m
		  en1 = min(en1,w[j-1][i][a][b+1][aa+1][bb]+2*data->eparam[6]+4*gap);

		  //case 19 - j and m stacked
		  en1 = min(en1,v[j-1][i][a][b+1][aa+1][bb]+3*data->eparam[10]
			    +penalty(i,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n,ct3,data)
			    +edangle3(j-1,i,j,ct1,data)+edangle5(m+1,n,m,ct3,data)
			    +2*data->eparam[6]+4*gap);

		}
                                
		if(bb>=1) {
		  //case 18 j,n
		  en1 = min(en1,w[j-1][i][a][b+1][aa][bb-1]+2*data->eparam[6]+gap);

		  //case 18 - j and n stacked
		  en1 = min(en1,v[j-1][i][a][b+1][aa][bb-1]+3*data->eparam[10]
			    +penalty(i,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n-1,ct3,data)
			    +edangle3(j-1,i,j,ct1,data)+edangle3(n-1,m,n,ct3,data)
			    +2*data->eparam[6]+gap);

		}
                                
		//case 17 j
		en1 = min(en1,w[j-1][i][a][b+1][aa][bb]+data->eparam[6]+2*gap);

		//case 17 - stack on j
		en1 = min(en1,v[j-1][i][a][b+1][aa][bb]+3*data->eparam[10]
			  +penalty(i,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n,ct3,data)
			  +edangle3(j-1,i,j,ct1,data)
			  +data->eparam[6]+2*gap);

	      }
                            
	      if(aa>=1) {
                                
		if(bb>=1) {
		  //case 42 i,k,n
		  en1 = min(en1,w[j][i+1][a][b][aa-1][bb-1]+3*data->eparam[6]+3*gap);

		  //case 42 - i, k and n stacked
		  en1 = min(en1,v[j][i+1][a][b][aa-1][bb-1]+3*data->eparam[10]
			    +penalty(i+1,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n-1,ct3,data)
			    +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle3(n-1,m,n,ct3,data)
			    +3*data->eparam[6]+3*gap);

		}
                                
		//case 41 i,k
		en1 = min(en1,w[j][i+1][a][b][aa-1][bb]+2*data->eparam[6]+gap);

		//case 41 - i and k stacked
		en1 = min(en1,v[j][i+1][a][b][aa-1][bb]+3*data->eparam[10]
			  +penalty(i+1,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n,ct3,data)
			  +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l,k,ct2,data)
			  +2*data->eparam[6]+gap);

	      }
                            
	      if(aa+1<2*maxseparation+2) {
                                
		if(bb>=1) {
		  //case 4 m,n
		  en1 = min(en1,w[j][i][a][b][aa+1][bb-1]+2*data->eparam[6]+4*gap);

		  //case 4 - m and n stacked
		  en1=min(en1,v[j][i][a][b][aa+1][bb-1]+3*data->eparam[10]
			  +penalty(i,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n-1,ct3,data)
			  +edangle3(n-1,m+1,n,ct3,data)+edangle5(m+1,n-1,m,ct3,data)
			  +2*data->eparam[6]+4*gap);

		}
                                
		if(bb+1<2*maxseparation+2) {
		  //case 23 j,l,m
		  en1 = min(en1,w[j-1][i][a][b][aa+1][bb+1]+3*data->eparam[6]+3*gap);

		  //case 23 - j, l and m stacked
		  en1 = min(en1,v[j-1][i][a][b][aa+1][bb+1]+3*data->eparam[10]
			    +penalty(i,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n,ct3,data)
			    +edangle3(l-1,k,l,ct2,data)+edangle3(j-1,i,j,ct1,data)+edangle5(m+1,n,m,ct3,data)
			    +3*data->eparam[6]+3*gap);

		}
                                
		//case 3 m
		en1 = min(en1,w[j][i][a][b][aa+1][bb]+data->eparam[6]+2*gap);

		//case 3 - stack on m
		en1=min(en1,v[j][i][a][b][aa+1][bb]+3*data->eparam[10]
			+penalty(i,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n,ct3,data)
			+edangle5(m+1,n,m,ct3,data)
			+data->eparam[6]+2*gap);

                                
		//case 24 j,l,m,n
		en1 = min(en1,w[j-1][i][a][b][aa+1][bb]+4*data->eparam[6]+2*gap);

		//case 24 - j, l, m and n stacked
		en1 = min(en1,v[j-1][i][a][b][aa+1][bb]+3*data->eparam[10]
			  +penalty(i,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
			  +edangle3(l-1,k,l,ct2,data)+edangle3(j-1,i,j,ct1,data)+edangle3(n-1,m+1,n,ct3,data)
			  +edangle5(m+1,n-1,m,ct3,data)
			  +4*data->eparam[6]+2*gap);

	      }
                            
	      if(bb>=1) {
		//case 2 n
		en1 = min(en1,w[j][i][a][b][aa][bb-1]+data->eparam[6]+2*gap);

		//case 2 - stack on n
		en1=min(en1,v[j][i][a][b][aa][bb-1]+3*data->eparam[10]
			+penalty(i,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n-1,ct3,data)
			+edangle3(n-1,m,n,ct3,data)
			+data->eparam[6]+2*gap);

                                
		//case 44 i,k,m,n
		en1 = min(en1,w[j][i+1][a][b][aa][bb-1]+4*data->eparam[6]+2*gap);

		//case 44 - i, k, m and n stacked
		en1 = min(en1,v[j][i+1][a][b][aa][bb-1]+3*data->eparam[10]
			  +penalty(i+1,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n-1,ct3,data)
			  +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle5(m+1,n-1,m,ct3,data)
			  +edangle3(n-1,m+1,n,ct3,data)
			  +4*data->eparam[6]+2*gap);

	      }
                            
	      if(bb+1<2*maxseparation+2) {
		//case 21 j,l
		en1 = min(en1,w[j-1][i][a][b][aa][bb+1]+2*data->eparam[6]+gap);

		//case 21 - j and l stacked
		en1 = min(en1,v[j-1][i][a][b][aa][bb+1]+3*data->eparam[10]
			  +penalty(i,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n,ct3,data)
			  +edangle3(l-1,k,l,ct2,data)+edangle3(j-1,i,j,ct1,data)
			  +2*data->eparam[6]+gap);

	      }
                            
	      /************************************** W3 ****************************************/
                            
	      //calculate the free energy of 2 fragments merged:
	      for (c=i+minloop;c<j-minloop;c++) {
		for (d=max(k+minloop,c-maxsep);d<l-minloop&&d<=c+maxsep;d++) {
		  e = d-c+maxsep;
		  for (f=max(m+minloop,d-maxsep);f<n-minloop&&f<=d+maxsep;f++ ) {
		    g = f-d+maxsep;
		    en1 = min(en1,w[c][i][a][e][aa][g]+w[j][c+1][e][b][g][bb]);	

		  }
		}
	      }
                            
                            
	      w[j][i][a][b][aa][bb]=en1;
                            
                            
	    }
	  }
	}
      }
    }
  }
    
    
  //    cout << "end of 6 nested for loops \n"<<flush;
    
  /*************************************** end fill V and W ********************************/
    
  /****************************************start fill W9 **********************************/
  cout << "W9 initialization \n"<<flush;
  //initialize mine[i][k][m]:
  for (i=0;i<=ct1->numofbases;i++){
    for (k=0;k<2*maxsep+2;k++){
      for(m=0;m<2*maxsep+2;m++){
	mine[i][k][m] = gap*abs(m-maxsep);  
	//				cout << "init: i="<<i<<", k="<<k<<", m="<<m<<"\n"<<flush;	
      }
    }
  }
    
  cout << "start fill W9 \n"<<flush;
  lowest = infinity;
  //calculate exterior energy:
  //store this in mine[i][k][m];
  //mine[i][k][m] refers to the lowest free energy for fragment 1-i and 1-k and 1-m with i and k and m aligned
  for (i=1;i<=ct1->numofbases;i++) {
    c = min(ct2->numofbases,i+maxsep);
    for (k=max(1,i-maxsep);k<=c;k++) {
      a = k-i+maxsep;
      //            cout << "bif at "<<i<<","<<k<<"\n"<<flush;	
      f = min(ct3->numofbases,k+maxsep);
      for(m=max(1,k-maxsep);m<=f;m++) {
	aa = m-k+maxsep;
	//				cout << "fill entry: i="<<i<<", k="<<k<<", m="<<m<<"\n"<<flush;	
	//				cout << "fill entry: i="<<i<<", a="<<a<<", aa="<<aa<<"\n"<<flush;	
	e = mine[i-1][a][aa];	// W9(i-1, k-1, m-1)
                

	if(a>=1) {
	  if(aa+1<2*maxsep+2)
	    e = min(e,mine[i][a-1][aa+1]+2*gap);  //W9(i, k-1, m)
	  e = min(e,mine[i][a-1][aa]+gap);  //W9(i, k-1, m-1)
	}
	if(a+1<2*maxsep+2) {
	  if(aa>=1)
	    e = min(e,mine[i-1][a+1][aa-1]+gap);  //W9(i-1, k, m-1)
	  e = min(e,mine[i-1][a+1][aa]+2*gap);  //W9(i-1, k, m)
	}
	if(aa>=1)
	  e = min(e,mine[i][a][aa-1]+2*gap);  //W9(i, k, m-1)
	if(aa+1<2*maxsep+2)
	  e = min(e,mine[i-1][a][aa+1]+gap);  //W9(i-1, k-1, m)
                
	// check end
	for (iprime=0;iprime+minloop<i;iprime++) {
	  u = min(iprime+maxsep-1,ct2->numofbases);
	  u = min(u,k-minloop-1);
	  for (kprime=max(0,iprime-maxsep);kprime<=u;kprime++) {
	    aprime = kprime-iprime+maxsep;
	    g = min(kprime+maxsep-1,ct3->numofbases);
	    g = min(g,m-minloop-1);
	    for (mprime=max(0,kprime-maxsep);mprime<=g;mprime++) {
                            
	      asegond = mprime-kprime+maxsep;
                            
	      //check whether mine[i][a][aa] is split so that the lowest free energy is
	      //an exterior fragment from 1 to iprime and a helix from iprime+1 to i;
                            
	      //check all possible alignments (in index a(k) and aprime(kprime) and aa(m) and asegond(mprime))
                            
	      //must consider whether an unpaired nucleotide is stacked
	      //is stacked onto each of the nucleotides in a helix
	      //There are 64 cases:
	      //consider them in this order (0 unstacked, 1 stacked)
	      //	iprime+1	i	kprime+1	k	mprime+1	m
	      //1		0		0		0		0		0		1
	      //2		0		0		0		0		1		0
	      //3		0		0		0		0		1		0
	      //4		0		0		0		0		1		1
	      //5		0		0		0		1		0		0
	      //6		0		0		0		1		0		1
	      //7		0		0		0		1		1		0
	      //8		0		0		0		1		1		1
	      //9		0		0		1		0		0		0
	      //10	0		0		1		0		0		1
	      //11	0		0		1		0		1		0
	      //12	0		0		1		0		1		1
	      //13	0		0  		1		1		0		0
	      //14	0		0		1		1		0		1
	      //15	0		0		1		1		1		0
	      //16	0		0		1		1		1		1
	      //17	0		1		0		0		0		0
	      //18	0		1		0		0		0		1
	      //19	0		1		0		0		1		0
	      //20	0		1		0		0		1		1
	      //21	0		1		0		1		0		0
	      //22	0		1		0		1		0		1
	      //23	0		1		0		1		1		0
	      //24	0		1		0		1		1		1
	      //25	0		1		1		0		0		0
	      //26	0		1		1		0		0		1
	      //27	0		1		1		0		1		0
	      //28	0		1		1		0		1		1
	      //29	0		1		1		1		0		0
	      //30	0		1		1		1		0		1
	      //31	0		1		1		1		1		0
	      //32	0		1		1		1		1		1
	      //33	1		0		0		0		0		0
	      //34	1		0		0		0		0		1
	      //35	1		0		0		0		1		0
	      //36	1		0		0		0		1		1
	      //37	1		0		0		1		0		0
	      //38	1		0		0		1		0		1
	      //39	1		0		0		1		1		0
	      //40	1		0		0		1		1		1
	      //41	1		0		1		0		0		0
	      //42	1		0		1		0		0		1
	      //43	1		0		1		0		1		0
	      //44	1		0		1		0		1		1
	      //45	1		0		1		1		0		0
	      //46	1		0		1		1		0		1
	      //47	1		0		1		1		1		0
	      //48	1		0		1		1		1		1
	      //49	1		1		0		0		0		0
	      //50	1		1		0		0		0		1
	      //51	1		1		0		0		1		0
	      //52	1		1		0		0		1		1
	      //53	1		1		0		1		0		0
	      //54	1		1		0		1		0		1
	      //55	1		1		0		1		1		0
	      //56	1		1		0		1		1		1
	      //57	1		1		1		0		0		0
	      //58	1		1		1		0		0		1
	      //59	1		1		1		0		1		0
	      //60	1		1		1		0		1		1
	      //61	1		1		1		1		0		0
	      //62	1		1		1		1		0		1
	      //63	1		1		1		1		1		0
	      //64	1		1		1		1		1		1
                            
                            
	      //note that for these exterior loops:
	      //iprime<i
	      //kprime<k
	      //mprime<m
                            
	      if(aprime>=1) {
                                
		if(a>=1) {
                                    
		  if(asegond+1<2*maxsep+2) {
                                        
		    if(aa+1<2*maxsep+2) {
		      //case 39 - iprime,k,mprime
		      e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime-1][a-1][asegond+1][aa+1]
			      +penalty(i,iprime+2,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m,mprime+2,ct3,data)
			      +edangle5(iprime+2,i,iprime+1,ct1,data)+edangle3(k-1,kprime+1,k,ct2,data)+edangle5(mprime+2,m,mprime+1,ct3,data)
			      +3*gap);
		    }
                                        
		    //case 40 - iprime,k,mprime,m
		    e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime-1][a-1][asegond+1][aa]
			    +penalty(i,iprime+2,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			    +edangle5(iprime+2,i,iprime+1,ct1,data)+edangle3(k-1,kprime+1,k,ct2,data)+edangle5(mprime+2,m-1,mprime+2,ct3,data)
			    +edangle3(m-1,mprime+2,m,ct3,data)
			    +2*gap);
		  }
                                    
		  if(aa+1<2*maxsep+2) {
		    //case 37 - iprime,k
		    e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime-1][a-1][asegond][aa+1]
			    +penalty(i,iprime+2,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m,mprime+1,ct3,data)
			    +edangle5(iprime+2,i,iprime+1,ct1,data)+edangle3(k-1,kprime+1,k,ct2,data)
			    +4*gap);
		  }
                                    
		  //case 38 - iprime,k,m
		  e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime-1][a-1][asegond][aa]
			  +penalty(i,iprime+2,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			  +edangle5(iprime+2,i,iprime+1,ct1,data)+edangle3(k-1,kprime+1,k,ct2,data)+edangle3(m-1,mprime+1,m,ct3,data)
			  +3*gap);
		}
                                
		if(a+1<2*maxsep+2) {
                                    
		  if(asegond+1<2*maxsep+2) {
                                        
		    if(aa>=1) {
		      //case 52 - iprime,i,mprime,m
		      e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime-1][a+1][asegond+1][aa-1]
			      +penalty(i-1,iprime+2,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			      +edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle5(mprime+2,m-1,mprime+1,ct3,data)
			      +edangle3(m-1,mprime+2,m,ct3,data)
			      +2*gap);
		    }
                                        
		    //case 51 - iprime,i,mprime
		    e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime-1][a+1][asegond+1][aa]
			    +penalty(i-1,iprime+2,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m,mprime+2,ct3,data)
			    +edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle5(mprime+2,m,mprime+1,ct3,data)
			    +3*gap);
		  }
                                    
		  if(aa>=1) {
		    //case 50 - iprime,i,m
		    e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime-1][a+1][asegond][aa-1]
			    +penalty(i-1,iprime+2,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			    +edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle3(m-1,mprime+1,m,ct3,data)
			    +3*gap);
		  }
                                    
		  //case 49 - iprime,i
		  e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime-1][a+1][asegond][aa]
			  +penalty(i-1,iprime+2,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m,mprime+1,ct3,data)
			  +edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)
			  +4*gap);
		}
                                
		if(asegond+1<2*maxsep+2) {
                                    
		  if(aa>=1) {
		    //case 36 - iprime,mprime,m
		    e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime-1][a][asegond+1][aa-1]
			    +penalty(i,iprime+2,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			    +edangle5(iprime+2,i,iprime+1,ct1,data)+edangle5(mprime+2,m-1,mprime+2,ct3,data)+edangle3(m-1,mprime+2,m,ct3,data)
			    +3*gap);
		  }
                                    
		  if(aa+1<2*maxsep+2) {
		    //case 55 - iprime,i,k,mprime
		    e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime-1][a][asegond+1][aa+1]
			    +penalty(i-1,iprime+2,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m,mprime+2,ct3,data)
			    +edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle3(k-1,kprime+1,k,ct2,data)
			    +edangle5(mprime+2,m-1,mprime+1,ct3,data)
			    +2*gap);
		  }
                                    
		  //case 35 - iprime,mprime
		  e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime-1][a][asegond+1][aa]
			  +penalty(i,iprime+2,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m,mprime+2,ct3,data)
			  +edangle5(iprime+2,i,iprime+1,ct1,data)+edangle5(mprime+2,m,mprime+1,ct3,data)
			  +gap);
                                    
		  //case 56 - iprime,i,k,mprime,m
		  e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime-1][a][asegond+1][aa]
			  +penalty(i-1,iprime+2,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			  +edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle3(k-1,kprime+1,k,ct2,data)
			  +edangle5(mprime+2,m-1,mprime+2,ct3,data)+edangle3(m-1,mprime+2,m,ct3,data)
			  +gap);
		}
                                
		if(aa>=1) {
		  //case 34 - iprime,m
		  e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime-1][a][asegond][aa-1]
			  +penalty(i,iprime+2,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			  +edangle5(iprime+2,i,iprime+1,ct1,data)+edangle3(m-1,mprime,m,ct3,data)
			  +4*gap);
		}
                                
		if(aa+1<2*maxsep+2) {
		  //case 53 - iprime,i,k
		  e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime-1][a][asegond][aa+1]
			  +penalty(i-1,iprime+2,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m,mprime+1,ct3,data)
			  +edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle3(k-1,kprime+1,k,ct2,data)
			  +3*gap);
		}
                                
		//case 33 - iprime
		e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime-1][a][asegond][aa]
			+penalty(i,iprime+2,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m,mprime+1,ct3,data)
			+edangle5(iprime+2,i,iprime+1,ct1,data)
			+2*gap);
                                
		//case 54 - iprime,i,k,m
		e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime-1][a][asegond][aa]
			+penalty(i-1,iprime+2,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			+edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle3(k-1,kprime+1,k,ct2,data)
			+edangle3(m-1,mprime+1,m,ct3,data)
			+2*gap);
	      }
                            
	      if(aprime+1<2*maxsep+2) {
                                
		if(a>=1) {
                                    
		  if(asegond>=1) {
                                        
		    if(aa+1<2*maxsep+2) {
		      //case 13 - kprime,k
		      e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime+1][a-1][asegond-1][aa+1]
			      +penalty(i,iprime+1,ct1,data)+penalty(k-1,kprime+2,ct2,data)+penalty(m,mprime+1,ct3,data)
			      +edangle3(k-1,kprime+2,k,ct2,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)
			      +4*gap);
		    }
                                        
		    //case 14 - kprime,k,m
		    e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime+1][a-1][asegond-1][aa]
			    +penalty(i,iprime+1,ct1,data)+penalty(k-1,kprime+2,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			    +edangle3(k-1,kprime+2,k,ct2,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)+edangle3(m-1,mprime+1,m,ct3,data)
			    +3*gap);
		  }
                                    
		  if(aa+1<2*maxsep+2) {
		    //case 15 - kprime,k,mprime
		    e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime+1][a-1][asegond][aa+1]
			    +penalty(i,iprime+1,ct1,data)+penalty(k-1,kprime+2,ct2,data)+penalty(m,mprime+2,ct3,data)
			    +edangle3(k-1,kprime+2,k,ct2,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)+edangle5(mprime+2,m,mprime+1,ct3,data)
			    +3*gap);
		  }
                                    
		  //case 16 - kprime,k,mprime,m
		  e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime+1][a-1][asegond][aa]
			  +penalty(i,iprime+1,ct1,data)+penalty(k-1,kprime+2,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			  +edangle3(k-1,kprime+2,k,ct2,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)+edangle5(mprime+2,m-1,mprime+1,ct3,data)
			  +edangle3(m-1,mprime+2,m,ct3,data)
			  +2*gap);
		}
                                
		if(a+1<2*maxsep+2) {
                                    
		  if(asegond>=1) {
                                        
		    if(aa>=1) {
		      //case 26 - i,kprime,m
		      e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime+1][a+1][asegond-1][aa-1]
			      +penalty(i-1,iprime+1,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			      +edangle5(kprime+2,k,kprime+1,ct2,data)+edangle3(i-1,iprime+1,i,ct1,data)+edangle3(m-1,mprime+1,m,ct3,data)
			      +3*gap);
		    }
                                        
		    //case 25 - i,kprime
		    e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime+1][a+1][asegond-1][aa]
			    +penalty(i-1,iprime+1,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m,mprime+1,ct3,data)
			    +edangle5(kprime+2,k,kprime+1,ct2,data)+edangle3(i-1,iprime+1,i,ct1,data)
			    +4*gap);
		  }
                                    
		  if(aa>=1) {
		    //case 28 - i,kprime,mprime,m
		    e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime+1][a+1][asegond][aa-1]
			    +penalty(i-1,iprime+1,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			    +edangle5(kprime+2,k,kprime+1,ct2,data)+edangle3(i-1,iprime+1,i,ct1,data)+edangle3(m-1,mprime+2,m,ct3,data)
			    +edangle5(mprime+2,m-1,mprime+1,ct3,data)
			    +2*gap);
		  }
                                    
		  //case 27 - i,kprime,mprime
		  e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime+1][a+1][asegond][aa]
			  +penalty(i-1,iprime+1,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m,mprime+2,ct3,data)
			  +edangle5(kprime+2,k,kprime+1,ct2,data)+edangle3(i-1,iprime+1,i,ct1,data)+edangle5(mprime+2,m,mprime+1,ct3,data)
			  +3*gap);
		}
                                
		if(asegond>=1) {
                                    
		  if(aa>=1) {
		    //case 10 - kprime,m
		    e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime+1][a][asegond-1][aa-1]
			    +penalty(i,iprime+1,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			    +edangle5(kprime+2,k,kprime+1,ct2,data)+edangle3(m-1,mprime+1,m,ct3,data)
			    +4*gap);
		  }
                                    
		  if(aa+1<2*maxsep+2) {
		    //case 29 - i,kprime,k
		    e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime+1][a][asegond-1][aa+1]
			    +penalty(i-1,iprime+1,ct1,data)+penalty(kprime+2,k-1,ct2,data)+penalty(m,mprime+1,ct3,data)
			    +edangle3(i-1,iprime+1,i,ct1,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)+edangle3(k-1,kprime+2,k,ct2,data)
			    +3*gap);
		  }
                                    
		  //case 9 - kprime
		  e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime+1][a][asegond-1][aa]
			  +penalty(i,iprime+1,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m,mprime+1,ct3,data)
			  +edangle5(kprime+2,k,kprime+1,ct2,data)
			  +2*gap);
                                    
		  //case 30 - i,kprime,k,m
		  e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime+1][a][asegond-1][aa]
			  +penalty(i-1,iprime+1,ct1,data)+penalty(kprime+2,k-1,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			  +edangle3(i-1,iprime+1,i,ct1,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)+edangle3(k-1,kprime+2,k,ct2,data)
			  +edangle3(m-1,mprime+2,m,ct3,data)
			  +2*gap);
		}
                                
		if(aa>=1) {
		  //case 12 - kprime,mprime,m
		  e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime+1][a][asegond][aa-1]
			  +penalty(i,iprime+1,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			  +edangle5(kprime+2,k,kprime+1,ct2,data)+edangle5(mprime+2,m-1,mprime+1,ct3,data)+edangle3(m-1,mprime+2,m,ct3,data)
			  +3*gap);
		}
                                
		if(aa+1<2*maxsep+2) {
		  //case 31 - i,kprime,k,mprime
		  e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime+1][a][asegond][aa+1]
			  +penalty(i-1,iprime+1,ct1,data)+penalty(kprime+2,k-1,ct2,data)+penalty(m,mprime+2,ct3,data)
			  +edangle3(i-1,iprime+1,i,ct1,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)+edangle3(k-1,kprime+2,k,ct2,data)
			  +edangle5(mprime+2,m,mprime+1,ct3,data)
			  +2*gap);
		}
                                
		//case 11 - kprime,mprime
		e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime+1][a][asegond][aa]
			+penalty(i,iprime+1,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m,mprime+2,ct3,data)
			+edangle5(kprime+2,k,kprime+1,ct2,data)+edangle5(mprime+2,m,mprime+1,ct3,data)
			+gap);
                                
		//case 32 - i,kprime,k,mprime,m
		e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime+1][a][asegond][aa]
			+penalty(i-1,iprime+1,ct1,data)+penalty(kprime+2,k-1,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			+edangle3(i-1,iprime+1,i,ct1,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)+edangle3(k-1,kprime+2,k,ct2,data)
			+edangle3(m-1,mprime+2,m,ct3,data)+edangle5(mprime+2,m-1,mprime+1,ct3,data)
			+gap);
	      }
                            
	      if(a>=1) {
                                
		if(asegond>=1) {
                                    
		  if(aa+1<2*maxsep+2) {
		    //case 45 - iprime,kprime,k
		    e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime][a-1][asegond-1][aa+1]
			    +penalty(i,iprime+2,ct1,data)+penalty(k-1,kprime+2,ct2,data)+penalty(m,mprime+1,ct3,data)
			    +edangle5(iprime+2,i,iprime+1,ct1,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)+edangle3(k-1,kprime+2,k,ct2,data)
			    +3*gap);
		  }
                                    
		  //case 46 - iprime,kprime,k,m
		  e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime][a-1][asegond-1][aa]
			  +penalty(i,iprime+2,ct1,data)+penalty(k-1,kprime+2,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			  +edangle5(iprime+2,i,iprime+1,ct1,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)+edangle3(k-1,kprime+2,k,ct2,data)
			  +edangle3(m-1,mprime+1,m,ct3,data)
			  +2*gap);
		}
                                
		if(asegond+1<2*maxsep+2) {
                                    
		  if(aa+1<2*maxsep+2) {
		    //case 7 - k,mprime
		    e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime][a-1][asegond+1][aa+1]
			    +penalty(i,iprime+1,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m,mprime+2,ct3,data)
			    +edangle3(k-1,kprime+1,k,ct2,data)+edangle5(mprime+2,m,mprime+1,ct3,data)
			    +4*gap);
		  }
                                    
		  //case 8 - k,mprime,m
		  e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime][a-1][asegond+1][aa]
			  +penalty(i,iprime+1,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			  +edangle3(k-1,kprime+1,k,ct2,data)+edangle3(m-1,mprime+2,m,ct3,data)+edangle5(mprime+2,m-1,mprime+1,ct3,data)
			  +3*gap);
		}
                                
		if(aa+1<2*maxsep+2) {
		  //case 5 - k
		  e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime][a-1][asegond][aa+1]
			  +penalty(i,iprime+1,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m,mprime+1,ct3,data)
			  +edangle3(k-1,kprime+1,k,ct2,data)
			  +2*gap);
                                    
		  //case 47 - iprime,kprime,k,mprime
		  e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime][a-1][asegond][aa+1]
			  +penalty(i,iprime+2,ct1,data)+penalty(k-1,kprime+2,ct2,data)+penalty(m,mprime+2,ct3,data)
			  +edangle5(iprime+2,i,iprime+1,ct1,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)+edangle3(k-1,kprime+2,k,ct2,data)
			  +edangle5(mprime+2,m,mprime+1,ct3,data)
			  +2*gap);
		}
                                
		//case 6 - k,m
		e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime][a-1][asegond][aa]
			+penalty(i,iprime+1,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			+edangle3(k-1,kprime+1,k,ct2,data)+edangle3(m-1,mprime+1,m,ct3,data)
			+gap);
                                
		//case 48 - iprime,kprime,k,mprime,m
		e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime][a-1][asegond][aa]
			+penalty(i,iprime+2,ct1,data)+penalty(k-1,kprime+2,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			+edangle5(iprime+2,i,iprime+1,ct1,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)+edangle3(k-1,kprime+2,k,ct2,data)
			+edangle5(mprime+2,m-1,mprime+1,ct3,data)+edangle3(m-1,mprime+2,m,ct3,data)
			+gap);
	      }
                            
	      if(a+1<2*maxsep+2) {
                                
		if(asegond>=1) {
                                    
		  if(aa>=1) {
		    //case 58 - iprime,i,kprime,m
		    e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime][a+1][asegond-1][aa-1]
			    +penalty(i-1,iprime+2,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			    +edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle5(kprime+2,k,kprime+1,ct2,data)
			    +edangle3(m-1,mprime+1,m,ct3,data)
			    +2*gap);
		  }
                                    
		  //case 57 - iprime,i,kprime
		  e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime][a+1][asegond-1][aa]
			  +penalty(i-1,iprime+2,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m,mprime+1,ct3,data)
			  +edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle5(kprime+2,k,kprime+1,ct2,data)
			  +3*gap);
		}
                                
		if(asegond+1<2*maxsep+2) {
                                    
		  if(aa>=1) {
		    //case 20 - i,mprime,m
		    e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime][a+1][asegond+1][aa-1]
			    +penalty(i-1,iprime+1,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			    +edangle3(i-1,iprime+1,i,ct1,data)+edangle3(m-1,mprime+2,m,ct3,data)+edangle5(mprime+2,m-1,mprime+1,ct3,data)
			    +3*gap);
		  }
                                    
		  //case 19 - i,mprime
		  e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime][a+1][asegond+1][aa]
			  +penalty(i-1,iprime+1,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m,mprime+2,ct3,data)
			  +edangle3(i-1,iprime+1,i,ct1,data)+edangle5(mprime+2,m,mprime+1,ct3,data)
			  +4*gap);
		}
                                
		if(aa>=1) {
		  //case 18 - i,m
		  e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime][a+1][asegond][aa-1]
			  +penalty(i-1,iprime+1,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			  +edangle3(i-1,iprime+1,i,ct1,data)+edangle3(m-1,mprime+1,m,ct3,data)
			  +gap);
                                    
		  //case 60 - iprime,i,kprime,mprime,m
		  e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime][a+1][asegond][aa-1]
			  +penalty(i-1,iprime+2,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			  +edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle5(kprime+2,k,kprime+1,ct2,data)
			  +edangle5(mprime+2,m-1,mprime+1,ct3,data)+edangle3(m-1,mprime+2,m,ct3,data)
			  +gap);
		}
                                
		//case 17 - i
		e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime][a+1][asegond][aa]
			+penalty(i-1,iprime+1,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m,mprime+1,ct3,data)
			+edangle3(i-1,iprime+1,i,ct1,data)
			+2*gap);
                                
		//case 59 - iprime,i,kprime,mprime
		e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime][a+1][asegond][aa]
			+penalty(i-1,iprime+2,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m,mprime+2,ct3,data)
			+edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle5(kprime+2,k,kprime+1,ct2,data)
			+edangle5(mprime+2,m-1,mprime+1,ct3,data)
			+2*gap);
	      }
                            
	      if(asegond>=1) {
                                
		if(aa>=1) {
		  //case 42 - iprime,kprime,m
		  e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime][a][asegond-1][aa-1]
			  +penalty(i,iprime+2,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			  +edangle5(iprime+2,i,iprime+1,ct1,data)+edangle5(kprime+2,k,kprime+1,ct2,data)+edangle3(m-1,mprime+1,m,ct3,data)
			  +3*gap);
		}
                                
		if(aa+1<2*maxsep+2) {
		  //case 61 - iprime,i,kprime,k
		  e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime][a][asegond-1][aa+1]
			  +penalty(i-1,iprime+2,ct1,data)+penalty(kprime+2,k-1,ct2,data)+penalty(m,mprime+1,ct3,data)
			  +edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)
			  +edangle3(k-1,kprime+2,k,ct2,data)
			  +2*gap);
		}
                                
		//case 41 - iprime,kprime
		e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime][a][asegond-1][aa]
			+penalty(i,iprime+2,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m,mprime+1,ct3,data)
			+edangle5(iprime+2,i,iprime+1,ct1,data)+edangle5(kprime+2,k,kprime+1,ct2,data)
			+gap);
                                
		//case 62 - iprime,i,kprime,k,m
		e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime][a][asegond-1][aa]
			+penalty(i-1,iprime+2,ct1,data)+penalty(kprime+2,k-1,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			+edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)
			+edangle3(k-1,kprime+2,k,ct2,data)+edangle3(m-1,mprime+2,m,ct3,data)
			+gap);
	      }
                            
	      if(asegond+1<2*maxsep+2) {
                                
		if(aa>=1) {
		  //case 4 - m,mprime
		  e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime][a][asegond+1][aa-1]
			  +penalty(i,iprime+1,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			  +edangle3(m-1,mprime+2,m,ct3,data)+edangle5(mprime+2,m-1,mprime+1,ct3,data)
			  +4*gap);
		}
                                
		if(aa+1<2*maxsep+2) {
		  //case 23 - i,k,mprime
		  e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime][a][asegond+1][aa+1]
			  +penalty(i-1,iprime+1,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m,mprime+2,ct3,data)
			  +edangle3(k-1,kprime+1,k,ct2,data)+edangle3(i-1,iprime+1,i,ct1,data)+edangle5(mprime+2,m,mprime+1,ct3,data)
			  +3*gap);
		}
                                
		//case 3 - mprime
		e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime][a][asegond+1][aa]
			+penalty(i,iprime+1,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m,mprime+2,ct3,data)
			+edangle5(mprime+2,m,mprime+1,ct3,data)
			+2*gap);
                                
		//case 24 - i,k,mprime,m
		e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime][a][asegond+1][aa]
			+penalty(i-1,iprime+1,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			+edangle3(k-1,kprime+1,k,ct2,data)+edangle3(i-1,iprime+1,i,ct1,data)+edangle3(m-1,mprime+2,m,ct3,data)
			+edangle5(mprime+2,m-1,mprime-1,ct3,data)
			+2*gap);
	      }
                            
	      if(aa>=1) {
		//case 2 - m
		e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime][a][asegond][aa-1]
			+penalty(i,iprime+1,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m-1,mprime+1,ct3,data)
			+edangle3(m-1,mprime+1,m,ct3,data)
			+2*gap);
                                
		//case 44 - iprime,kprime,mprime,m
		e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime][a][asegond][aa-1]
			+penalty(i,iprime+2,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m-1,mprime+2,ct3,data)
			+edangle5(iprime+2,i,iprime+1,ct1,data)+edangle5(kprime+2,k,kprime+1,ct2,data)+edangle5(mprime+2,m-1,mprime+1,ct3,data)
			+edangle3(m-1,mprime+2,m,ct3,data)
			+2*gap);
	      }
                            
	      if(aa+1<2*maxsep+2) {
		//case 21 - i,k
		e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime][a][asegond][aa+1]
			+penalty(i-1,iprime+1,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m,mprime+1,ct3,data)
			+edangle3(k-1,kprime+1,k,ct2,data)+edangle3(i-1,iprime+1,i,ct1,data)
			+gap);
                                
		//case 63 - iprime,i,kprime,k,mprime
		e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime][a][asegond][aa+1]
			+penalty(i-1,iprime+2,ct1,data)+penalty(kprime+2,k-1,ct2,data)+penalty(m,mprime+2,ct3,data)
			+edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)
			+edangle3(k-1,kprime+2,k,ct2,data)+edangle5(mprime+2,m-1,mprime+1,ct3,data)
			+gap);
	      }
                            
	      //case 1 -
	      e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+1][aprime][a][asegond][aa]
		      +penalty(i,iprime+1,ct1,data)+penalty(k,kprime+1,ct2,data)+penalty(m,mprime+1,ct3,data));
                            
	      //case 22 - i,k,m
	      e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+1][aprime][a][asegond][aa]
		      +penalty(i-1,iprime+1,ct1,data)+penalty(k-1,kprime+1,ct2,data)+penalty(m-1,mprime+1,ct3,data)
		      +edangle3(k-1,kprime+1,k,ct2,data)+edangle3(i-1,iprime+1,i,ct1,data)+edangle3(m-1,mprime+1,m,ct3,data));
                            
	      //case 43 - iprime,kprime,mprime
	      e = min(e,mine[iprime][aprime][asegond]+v[i][iprime+2][aprime][a][asegond][aa]
		      +penalty(i,iprime+2,ct1,data)+penalty(kprime+2,k,ct2,data)+penalty(m,mprime+2,ct3,data)
		      +edangle5(iprime+2,i,iprime+1,ct1,data)+edangle5(kprime+2,k,kprime+1,ct2,data)+edangle5(mprime+2,m,mprime+1,ct3,data));
                            
	      //case 64 - iprime,i,kprime,k,mprime,m
	      e = min(e,mine[iprime][aprime][asegond]+v[i-1][iprime+2][aprime][a][asegond][aa]
		      +penalty(i-1,iprime+2,ct1,data)+penalty(kprime+2,k-1,ct2,data)+penalty(m-1,mprime+2,ct3,data)
		      +edangle5(iprime+2,i-1,iprime+1,ct1,data)+edangle3(i-1,iprime+2,i,ct1,data)+edangle5(kprime+2,k-1,kprime+1,ct2,data)
		      +edangle3(k-1,kprime+2,k,ct2,data)+edangle3(m-1,mprime+2,m,ct3,data)+edangle5(mprime+2,m-1,mprime+1,ct3,data));
                            
                            
                            
	    }
	  }
	}
	mine[i][a][aa] = e;
	nogaps=3*max((ct1->numofbases)-i,max((ct2->numofbases)-k,(ct3->numofbases)-m))-(ct1->numofbases)+i-(ct2->numofbases)+k-(ct3->numofbases)+m;
	lowest = min(lowest,e+(gap*nogaps));
	//				cout << "fill exit: i="<<i<<", k="<<k<<", m="<<m<<"\n"<<flush;
	//              cout << "fill exit: i="<<i<<", a="<<a<<", aa="<<aa<<"\n"<<flush;
	//				cout << "fill exit; lowest="<<lowest<<"\n"<<flush;
      }
    }
  }
  cout << "end of fill W9 \n"<<flush;
  /************************************* end fill W9 ************************************/
    
  /************************************* Insert Traceback ************************************/
    
  //declarations for traceback portion:
  stackclass stack;
  register int constantclosure;
    
  if (totalscore!=NULL)
    *totalscore = lowest;
    
    
  if (lowest==0) {
    cout << "Lowest is zero!!! This is too bad!!!\n"<<flush;
    //return with empty CT files
    return;
        
  }
    
  //now traceback:
  //find all exterior fragments and place them on a stack
    
  /*********************find i_lowest, k_lowest, m_lowest***************************************/
  found = false; 
  for (i=1;i<=ct1->numofbases&&!found;i++) {
    for (k=0;k<=2*maxsep+1&&!found;k++) {
      p = k+i-maxsep;  // old k
      if (p>0&&p<=ct2->numofbases) {
	for(m=0;m<=2*maxsep+1&&!found;m++) {
	  //check
	  q=m+p-maxsep; // old m 
	  //check end
	  if(q>0&&q<=ct3->numofbases) {
                        

	    nogaps=3*max((ct1->numofbases)-i,max((ct2->numofbases)-p,(ct3->numofbases)-q))-(ct1->numofbases)+i-(ct2->numofbases)+p-(ct3->numofbases)+q;
						
	    if (mine[i][k][m]+gap*nogaps==lowest) { 
	      en2=i;
	      a=k;
	      aa=m;
	      found = true;
	      en4=mine[i][k][m];

	    }
                        
                        
	  }
	}
      }
    }
  }
  if (!found) {
    cout << "Traceback error at start!\n"<<flush;
    return;
  }
  cout << "end of find lowest \n"<<flush;

  /******************************** end find i_lowest, k_lowest, m_lowest ******************************************************/
  //initialize the structure portion of the ct files:
  for (i=1;i<=ct1->numofbases;i++) ct1->basepr[1][i]=0;
  for (i=1;i<=ct2->numofbases;i++) ct2->basepr[1][i]=0;
  for (i=1;i<=ct3->numofbases;i++) ct3->basepr[1][i]=0;
  ct1->numofstructures=1;
  ct2->numofstructures=1;
  ct3->numofstructures=1;
    
  //initialize the alignment
  for (i=1;i<=ct1->numofbases;i++) alignment1[i] = 0;
  for (k=1;k<=ct2->numofbases;k++) alignment2[k] = 0;
    
  //	cout << "end of find init align \n"<<flush;
    
  /******************************** find exterior ******************************************************/
       
       
  i = en2;
   
  while(i>0) {
       
    found = false;
    k = a+i-maxsep;
    m = aa+k-maxsep;
		
    if (en4==mine[i-1][a][aa]) {
      found = true;
      en4 = mine[i-1][a][aa];		//W9(i-1, k-1, m-1)
      i--;			
    }
        
    if (aa+1<2*maxsep+2 && !found)
      {
	if(a>=1 && !found) {

	  if (en4==mine[i][a-1][aa+1]+2*gap){
	    found = true;
	    en4 = mine[i][a-1][aa+1];		//W9(i, k-1, m)
	    a--;
	    aa++;				
	  }
	}


	if(en4==mine[i-1][a][aa+1]+gap && !found) {
	  found = true;
	  en4 = mine[i-1][a][aa+1];		//W9(i-1, k-1, m)
	  i--;
	  aa++;				                
	}
      }
		
		
    if(a>=1 && !found) {

      if (en4==mine[i][a-1][aa]+gap && !found){  
	found = true;
	en4 = mine[i][a-1][aa];		//W9(i, k-1, m-1)
	a--;				
      }
    }

    if (a+1<2*maxsep+2 && !found) {
      if(aa>=1 && !found) {

	if (en4==mine[i-1][a+1][aa-1]+gap){
	  found = true;
	  en4 = mine[i-1][a+1][aa-1];		//W9(i-1, k, m-1)
	  i--;
	  a++;
	  aa--;				
	}
      }


      if (en4==mine[i-1][a+1][aa]+2*gap && !found){
	found = true;
	en4 = mine[i-1][a+1][aa];		//W9(i-1, k, m)
	i--;
	a++;				
      }
    }

    if(aa>=1 && !found) {

      if (en4==mine[i][a][aa-1]+2*gap){
	found=true;
	en4 = mine[i][a][aa-1];		//W9(i, k, m-1)
	aa--;				
      }
    }
		
    for (j=0;j+minloop<i&&!found;j++) {
      for (l=max(0,j-maxsep);l<j+maxsep&&l<=ct2->numofbases&&l+minloop<k&&!found;l++) {
	b = l-j+maxsep;
	for(n=max(0,l-maxsep);n<l+maxsep&&n<=ct3->numofbases&&n+minloop<m&&!found;n++) {
	  bb= n-l+maxsep;
                    
	  //check whether mine[i][a][aa] is split so that the lowest free energy is
	  //an exterior fragment from 1 to j and a helix from j+1 to i;
                    
	  //check all possible alignments (in index a(k), b(l), aa(m) and bb(n))
                    
	  //must consider whether an unpaired nucleotide is stacked
	  //is stacked onto each of the nucleotides in a helix
	  //There are 64 cases:
	  //consider them in this order (0 unstacked, 1 stacked)
	  //		j+1	i	l+1	k	n+1	m
	  //1		0	0	0	0	0	1
	  //2		0	0	0	0	1	0
	  //3		0	0	0	0	1	0
	  //4		0	0	0	0	1	1
	  //5		0	0	0	1	0	0
	  //6		0	0	0	1	0	1
	  //7		0	0	0	1	1	0
	  //8		0	0	0	1	1	1
	  //9		0	0	1	0	0	0
	  //10	0	0	1	0	0	1
	  //11	0	0	1	0	1	0
	  //12	0	0	1	0	1	1
	  //13	0	0	1	1	0	0
	  //14	0	0	1	1	0	1
	  //15	0	0	1	1	1	0
	  //16	0	0	1	1	1	1
	  //17	0	1	0	0	0	0
	  //18	0	1	0	0	0	1
	  //19	0	1	0	0	1	0
	  //20	0	1	0	0	1	1
	  //21	0	1	0	1	0	0
	  //22	0	1	0	1	0	1
	  //23	0	1	0	1	1	0
	  //24	0	1	0	1	1	1
	  //25	0	1	1	0	0	0
	  //26	0	1	1	0	0	1
	  //27	0	1	1	0	1	0
	  //28	0	1	1	0	1	1
	  //29	0	1	1	1	0	0
	  //30	0	1	1	1	0	1
	  //31	0	1	1	1	1	0
	  //32	0	1	1	1	1	1
	  //33	1	0	0	0	0	0
	  //34	1	0	0	0	0	1
	  //35	1	0	0	0	1	0
	  //36	1	0	0	0	1	1
	  //37	1	0	0	1	0	0
	  //38	1	0	0	1	0	1
	  //39	1	0	0	1	1	0
	  //40	1	0	0	1	1	1
	  //41	1	0	1	0	0	0
	  //42	1	0	1	0	0	1
	  //43	1	0	1	0	1	0
	  //44	1	0	1	0	1	1
	  //45	1	0	1	1	0	0
	  //46	1	0	1	1	0	1
	  //47	1	0	1	1	1	0
	  //48	1	0	1	1	1	1
	  //49	1	1	0	0	0	0
	  //50	1	1	0	0	0	1
	  //51	1	1	0	0	1	0
	  //52	1	1	0	0	1	1
	  //53	1	1	0	1	0	0
	  //54	1	1	0	1	0	1
	  //55	1	1	0	1	1	0
	  //56	1	1	0	1	1	1
	  //57	1	1	1	0	0	0
	  //58	1	1	1	0	0	1
	  //59	1	1	1	0	1	0
	  //60	1	1	1	0	1	1
	  //61	1	1	1	1	0	0
	  //62	1	1	1	1	0	1
	  //63	1	1	1	1	1	0
	  //64	1	1	1	1	1	1
                    
	  //note that for these eterior loops:
	  //j<i
	  //l<k
	  //n<m
                    
                    
	  if(b>=1) {
                        
	    if(a>=1) {
                            
	      if(bb+1<2*maxsep+2) {
                                
		if(aa+1<2*maxsep+2) {
		  //case 39 - j,k,n

		  if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b-1][a-1][bb+1][aa+1]
				 +penalty(i,j+2,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m,n+2,ct3,data)
				 +edangle5(j+2,i,j+1,ct1,data)+edangle3(k-1,l+1,k,ct2,data)+edangle5(n+2,m,n+1,ct3,data)
				 +3*gap){
		    stack.push(j+2,i,b-1,a-1,bb+1,aa+1,v[i][j+2][b-1][a-1][bb+1][aa+1]);
		    found = true;
		    en4=mine[j][b][bb];
		    i = j;
		    a = b;
		    aa = bb;										
		  }
		}
                                
		//case 40 - j,k,n,m

		if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b-1][a-1][bb+1][aa]
			       +penalty(i,j+2,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m-1,n+2,ct3,data)
			       +edangle5(j+2,i,j+1,ct1,data)+edangle3(k-1,l+1,k,ct2,data)+edangle5(n+2,m-1,n+2,ct3,data)
			       +edangle3(m-1,n+2,m,ct3,data)
			       +2*gap){
		  stack.push(j+2,i,b-1,a-1,bb+1,aa,v[i][j+2][b-1][a-1][bb+1][aa]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      if(aa+1<2*maxsep+2) {
		//case 37 - j,k

		if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b-1][a-1][bb][aa+1]
			       +penalty(i,j+2,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m,n+1,ct3,data)
			       +edangle5(j+2,i,j+1,ct1,data)+edangle3(k-1,l+1,k,ct2,data)
			       +4*gap){
		  stack.push(j+2,i,b-1,a-1,bb,aa+1,v[i][j+2][b-1][a-1][bb][aa+1]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      //case 38 - j,k,m

	      if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b-1][a-1][bb][aa]
			     +penalty(i,j+2,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m-1,n+1,ct3,data)
			     +edangle5(j+2,i,j+1,ct1,data)+edangle3(k-1,l+1,k,ct2,data)+edangle3(m-1,n+1,m,ct3,data)
			     +3*gap){
		stack.push(j+2,i,b-1,a-1,bb,aa,v[i][j+2][b-1][a-1][bb][aa]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(a+1<2*maxsep+2) {
                            
	      if(bb+1<2*maxsep+2) {
                                
		if(aa>=1) {
		  //case 52 - j,i,n,m

		  if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b-1][a+1][bb+1][aa-1]
				 +penalty(i-1,j+2,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m-1,n+2,ct3,data)
				 +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle5(n+2,m-1,n+1,ct3,data)
				 +edangle3(m-1,n+2,m,ct3,data)
				 +2*gap){
		    stack.push(j+2,i-1,b-1,a+1,bb+1,aa-1,v[i-1][j+2][b-1][a+1][bb+1][aa-1]);
		    found = true;
		    en4=mine[j][b][bb];
		    i = j;
		    a = b;
		    aa = bb;										
		  }
		}
                                
		//case 51 - j,i,n

		if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b-1][a+1][bb+1][aa]
			       +penalty(i-1,j+2,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m,n+2,ct3,data)
			       +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle5(n+2,m,n+1,ct3,data)
			       +3*gap){
		  stack.push(j+2,i-1,b-1,a+1,bb+1,aa,v[i-1][j+2][b-1][a+1][bb+1][aa]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      if(aa>=1) {
		//case 50 - j,i,m

		if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b-1][a+1][bb][aa-1]
			       +penalty(i-1,j+2,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m-1,n+1,ct3,data)
			       +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle3(m-1,n+1,m,ct3,data)
			       +3*gap){
		  stack.push(j+2,i-1,b-1,a+1,bb,aa-1,v[i-1][j+2][b-1][a+1][bb][aa-1]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      //case 49 - j,i

	      if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b-1][a+1][bb][aa]
			     +penalty(i-1,j+2,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m,n+1,ct3,data)
			     +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)
			     +4*gap){
		stack.push(j+2,i-1,b-1,a+1,bb,aa,v[i-1][j+2][b-1][a+1][bb][aa]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(bb+1<2*maxsep+2) {
                            
	      if(aa>=1) {
		//case 36 - j,n,m

		if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b-1][a][bb+1][aa-1]
			       +penalty(i,j+2,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m-1,n+2,ct3,data)
			       +edangle5(j+2,i,j+1,ct1,data)+edangle5(n+2,m-1,n+2,ct3,data)+edangle3(m-1,n+2,m,ct3,data)
			       +3*gap){
		  stack.push(j+2,i,b-1,a,bb+1,aa-1,v[i][j+2][b-1][a][bb+1][aa-1]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      if(aa+1<2*maxsep+2) {
		//case 55 - j,i,k,n

		if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b-1][a][bb+1][aa+1]
			       +penalty(i-1,j+2,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m,n+2,ct3,data)
			       +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle3(k-1,l+1,k,ct2,data)
			       +edangle5(n+2,m-1,n+1,ct3,data)
			       +2*gap){
		  stack.push(j+2,i-1,b-1,a,bb+1,aa+1,v[i-1][j+2][b-1][a][bb+1][aa+1]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      //case 35 - j,n

	      if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b-1][a][bb+1][aa]
			     +penalty(i,j+2,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m,n+2,ct3,data)
			     +edangle5(j+2,i,j+1,ct1,data)+edangle5(n+2,m,n+1,ct3,data)
			     +gap){
		stack.push(j+2,i,b-1,a,bb+1,aa,v[i][j+2][b-1][a][bb+1][aa]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
                            
	      //case 56 - j,i,k,n,m

	      if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b-1][a][bb+1][aa]
			     +penalty(i-1,j+2,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m-1,n+2,ct3,data)
			     +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle3(k-1,l+1,k,ct2,data)
			     +edangle5(n+2,m-1,n+2,ct3,data)+edangle3(m-1,n+2,m,ct3,data)
			     +gap){
		stack.push(j+2,i-1,b-1,a,bb+1,aa,v[i-1][j+2][b-1][a][bb+1][aa]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(aa>=1) {
	      //case 34 - j,m

	      if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b-1][a][bb][aa-1]
			     +penalty(i,j+2,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m-1,n+1,ct3,data)
			     +edangle5(j+2,i,j+1,ct1,data)+edangle3(m-1,n,m,ct3,data)
			     +4*gap){
		stack.push(j+2,i,b-1,a,bb,aa-1,v[i][j+2][b-1][a][bb][aa-1]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(aa+1<2*maxsep+2) {
	      //case 53 - j,i,k

	      if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b-1][a][bb][aa+1]
			     +penalty(i-1,j+2,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m,n+1,ct3,data)
			     +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle3(k-1,l+1,k,ct2,data)
			     +3*gap){
		stack.push(j+2,i-1,b-1,a,bb,aa+1,v[i-1][j+2][b-1][a][bb][aa+1]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    //case 33 - j

	    if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b-1][a][bb][aa]
			   +penalty(i,j+2,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m,n+1,ct3,data)
			   +edangle5(j+2,i,j+1,ct1,data)
			   +2*gap){
	      stack.push(j+2,i,b-1,a,bb,aa,v[i][j+2][b-1][a][bb][aa]);
	      found = true;
	      en4=mine[j][b][bb];
	      i = j;
	      a = b;
	      aa = bb;							
	    }
                        
	    //case 54 - j,i,k,m

	    if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b-1][a][bb][aa]
			   +penalty(i-1,j+2,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m-1,n+1,ct3,data)
			   +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle3(k-1,l+1,k,ct2,data)
			   +edangle3(m-1,n+1,m,ct3,data)
			   +2*gap){
	      stack.push(j+2,i-1,b-1,a,bb,aa,v[i-1][j+2][b-1][a][bb][aa]);
	      found = true;
	      en4=mine[j][b][bb];
	      i = j;
	      a = b;
	      aa = bb;							
	    }
	  }
                    
	  if(b+1<2*maxsep+2) {
                        
	    if(a>=1) {
                            
	      if(bb>=1) {
                                
		if(aa+1<2*maxsep+2) {
		  //case 13 - l,k

		  if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b+1][a-1][bb-1][aa+1]
				 +penalty(i,j+1,ct1,data)+penalty(k-1,l+2,ct2,data)+penalty(m,n+1,ct3,data)
				 +edangle3(k-1,l+2,k,ct2,data)+edangle5(l+2,k-1,l+1,ct2,data)
				 +4*gap){
		    stack.push(j+1,i,b+1,a-1,bb-1,aa+1,v[i][j+1][b+1][a-1][bb-1][aa+1]);
		    found = true;
		    en4=mine[j][b][bb];
		    i = j;
		    a = b;
		    aa = bb;										
		  }
		}
                                
		//case 14 - l,k,m

		if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b+1][a-1][bb-1][aa]
			       +penalty(i,j+1,ct1,data)+penalty(k-1,l+2,ct2,data)+penalty(m-1,n+1,ct3,data)
			       +edangle3(k-1,l+2,k,ct2,data)+edangle5(l+2,k-1,l+1,ct2,data)+edangle3(m-1,n+1,m,ct3,data)
			       +3*gap){
		  stack.push(j+1,i,b+1,a-1,bb-1,aa,v[i][j+1][b+1][a-1][bb-1][aa]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      if(aa+1<2*maxsep+2) {
		//case 15 - l,k,n

		if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b+1][a-1][bb][aa+1]
			       +penalty(i,j+1,ct1,data)+penalty(k-1,l+2,ct2,data)+penalty(m,n+2,ct3,data)
			       +edangle3(k-1,l+2,k,ct2,data)+edangle5(l+2,k-1,l+1,ct2,data)+edangle5(n+2,m,n+1,ct3,data)
			       +3*gap){
		  stack.push(j+1,i,b+1,a-1,bb,aa+1,v[i][j+1][b+1][a-1][bb][aa+1]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      //case 16 - l,k,n,m

	      if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b+1][a-1][bb][aa]
			     +penalty(i,j+1,ct1,data)+penalty(k-1,l+2,ct2,data)+penalty(m-1,n+2,ct3,data)
			     +edangle3(k-1,l+2,k,ct2,data)+edangle5(l+2,k-1,l+1,ct2,data)+edangle5(n+2,m-1,n+1,ct3,data)
			     +edangle3(m-1,n+2,m,ct3,data)
			     +2*gap){
		stack.push(j+1,i,b+1,a-1,bb,aa,v[i][j+1][b+1][a-1][bb][aa]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(a+1<2*maxsep+2) {
                            
	      if(bb>=1) {
                                
		if(a>=1) {
		  //case 26 - i,l,m

		  if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b+1][a+1][bb-1][aa-1]
				 +penalty(i-1,j+1,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m-1,n+1,ct3,data)
				 +edangle5(l+2,k,l+1,ct2,data)+edangle3(i-1,j+1,i,ct1,data)+edangle3(m-1,n+1,m,ct3,data)
				 +3*gap){
		    stack.push(j+1,i-1,b+1,a+1,bb-1,aa-1,v[i-1][j+1][b+1][a+1][bb-1][aa-1]);
		    found = true;
		    en4=mine[j][b][bb];
		    i = j;
		    a = b;
		    aa = bb;										
		  }
		}
                                
		//case 25 - i,l

		if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b+1][a+1][bb-1][aa]
			       +penalty(i-1,j+1,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m,n+1,ct3,data)
			       +edangle5(l+2,k,l+1,ct2,data)+edangle3(i-1,j+1,i,ct1,data)
			       +4*gap){
		  stack.push(j+1,i-1,b+1,a+1,bb-1,aa,v[i-1][j+1][b+1][a+1][bb-1][aa]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      if(aa>=1) {
		//case 28 - i,l,m,n

		if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b+1][a+1][bb][aa-1]
			       +penalty(i-1,j+1,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m-1,n+2,ct3,data)
			       +edangle5(l+2,k,l+1,ct2,data)+edangle3(i-1,j+1,i,ct1,data)+edangle3(m-1,n+2,m,ct3,data)
			       +edangle5(n+2,m-1,n+1,ct3,data)
			       +2*gap){
		  stack.push(j+1,i-1,b+1,a+1,bb,aa-1,v[i-1][j+1][b+1][a+1][bb][aa-1]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      //case 27 - i,l,n

	      if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b+1][a+1][bb][aa]
			     +penalty(i-1,j+1,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m,n+2,ct3,data)
			     +edangle5(l+2,k,l+1,ct2,data)+edangle3(i-1,j+1,i,ct1,data)+edangle5(n+2,m,n+1,ct3,data)
			     +3*gap){
		stack.push(j+1,i-1,b+1,a+1,bb,aa,v[i-1][j+1][b+1][a+1][bb][aa]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(bb>=1) {
                            
	      if(aa>=1) {
		//case 10 - l,m

		if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b+1][a][bb-1][aa-1]
			       +penalty(i,j+1,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m-1,n+1,ct3,data)
			       +edangle5(l+2,k,l+1,ct2,data)+edangle3(m-1,n+1,m,ct3,data)
			       +4*gap){
		  stack.push(j+1,i,b+1,a,bb-1,aa-1,v[i][j+1][b+1][a][bb-1][aa-1]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      if(aa+1<2*maxsep+2) {
		//case 29 - i,l,k

		if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b+1][a][bb-1][aa+1]
			       +penalty(i-1,j+1,ct1,data)+penalty(l+2,k-1,ct2,data)+penalty(m,n+1,ct3,data)
			       +edangle3(i-1,j+1,i,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)+edangle3(k-1,l+2,k,ct2,data)
			       +3*gap){
		  stack.push(j+1,i-1,b+1,a,bb-1,aa+1,v[i-1][j+1][b+1][a][bb-1][aa+1]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      //case 9 - l

	      if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b+1][a][bb-1][aa]
			     +penalty(i,j+1,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m,n+1,ct3,data)
			     +edangle5(l+2,k,l+1,ct2,data)
			     +2*gap){
		stack.push(j+1,i,b+1,a,bb-1,aa,v[i][j+1][b+1][a][bb-1][aa]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
                            
	      //case 30 - i,l,k,m

	      if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b+1][a][bb-1][aa]
			     +penalty(i-1,j+1,ct1,data)+penalty(l+2,k-1,ct2,data)+penalty(m-1,n+1,ct3,data)
			     +edangle3(i-1,j+1,i,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)+edangle3(k-1,l+2,k,ct2,data)
			     +edangle3(m-1,n+2,m,ct3,data)
			     +2*gap){
		stack.push(j+1,i-1,b+1,a,bb-1,aa,v[i-1][j+1][b+1][a][bb-1][aa]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(aa>=1) {
	      //case 12 - l,n,m

	      if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b+1][a][bb][aa-1]
			     +penalty(i,j+1,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m-1,n+2,ct3,data)
			     +edangle5(l+2,k,l+1,ct2,data)+edangle5(n+2,m-1,n+1,ct3,data)+edangle3(m-1,n+2,m,ct3,data)
			     +3*gap){
		stack.push(j+1,i,b+1,a,bb,aa-1,v[i][j+1][b+1][a][bb][aa-1]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(aa+1<2*maxsep+2) {
	      //case 31 - i,l,k,n

	      if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b+1][a][bb][aa+1]
			     +penalty(i-1,j+1,ct1,data)+penalty(l+2,k-1,ct2,data)+penalty(m,n+2,ct3,data)
			     +edangle3(i-1,j+1,i,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)+edangle3(k-1,l+2,k,ct2,data)
			     +edangle5(n+2,m,n+1,ct3,data)
			     +2*gap){
		stack.push(j+1,i-1,b+1,a,bb,aa+1,v[i-1][j+1][b+1][a][bb][aa+1]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    //case 11 - l,n

	    if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b+1][a][bb][aa]
			   +penalty(i,j+1,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m,n+2,ct3,data)
			   +edangle5(l+2,k,l+1,ct2,data)+edangle5(n+2,m,n+1,ct3,data)
			   +gap){
	      stack.push(j+1,i,b+1,a,bb,aa,v[i][j+1][b+1][a][bb][aa]);
	      found = true;
	      en4=mine[j][b][bb];
	      i = j;
	      a = b;
	      aa = bb;							
	    }
                        
	    //case 32 - i,l,k,n,m

	    if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b+1][a][bb][aa]
			   +penalty(i-1,j+1,ct1,data)+penalty(l+2,k-1,ct2,data)+penalty(m-1,n+2,ct3,data)
			   +edangle3(i-1,j+1,i,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)+edangle3(k-1,l+2,k,ct2,data)
			   +edangle3(m-1,n+2,m,ct3,data)+edangle5(n+2,m-1,n+1,ct3,data)
			   +gap){
	      stack.push(j+1,i-1,b+1,a,bb,aa,v[i-1][j+1][b+1][a][bb][aa]);
	      found = true;
	      en4=mine[j][b][bb];
	      i = j;
	      a = b;
	      aa = bb;							
	    }
	  }
                    
                    
	  if(a>=1) {
                        
	    if(bb>=1) {
                            
	      if(aa+1<2*maxsep+2) {
		//case 45 - j,l,k

		if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b][a-1][bb-1][aa+1]
			       +penalty(i,j+2,ct1,data)+penalty(k-1,l+2,ct2,data)+penalty(m,n+1,ct3,data)
			       +edangle5(j+2,i,j+1,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)+edangle3(k-1,l+2,k,ct2,data)
			       +3*gap){
		  stack.push(j+2,i,b,a-1,bb-1,aa+1,v[i][j+2][b][a-1][bb-1][aa+1]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      //case 46 - j,l,k,m

	      if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b][a-1][bb-1][aa]
			     +penalty(i,j+2,ct1,data)+penalty(k-1,l+2,ct2,data)+penalty(m-1,n+1,ct3,data)
			     +edangle5(j+2,i,j+1,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)+edangle3(k-1,l+2,k,ct2,data)
			     +edangle3(m-1,n+1,m,ct3,data)
			     +2*gap){
		stack.push(j+2,i,b,a-1,bb-1,aa,v[i][j+2][b][a-1][bb-1][aa]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(bb+1<2*maxsep+2) {
                            
	      if(aa+1<2*maxsep+2) {
		//case 7 - k,n

		if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b][a-1][bb+1][aa+1]
			       +penalty(i,j+1,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m,n+2,ct3,data)
			       +edangle3(k-1,l+1,k,ct2,data)+edangle5(n+2,m,n+1,ct3,data)
			       +4*gap){
		  stack.push(j+1,i,b,a-1,bb+1,aa+1,v[i][j+1][b][a-1][bb+1][aa+1]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;																		
		}
	      }
                            
	      //case 8 - k,n,m

	      if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b][a-1][bb+1][aa]
			     +penalty(i,j+1,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m-1,n+2,ct3,data)
			     +edangle3(k-1,l+1,k,ct2,data)+edangle3(m-1,n+2,m,ct3,data)+edangle5(n+2,m-1,n+1,ct3,data)
			     +3*gap){
		stack.push(j+1,i,b,a-1,bb+1,aa,v[i][j+1][b][a-1][bb+1][aa]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(aa+1<2*maxsep+2) {
	      //case 5 - k

	      if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b][a-1][bb][aa+1]
			     +penalty(i,j+1,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m,n+1,ct3,data)
			     +edangle3(k-1,l+1,k,ct2,data)
			     +2*gap){
		stack.push(j+1,i,b,a-1,bb,aa+1,v[i][j+1][b][a-1][bb][aa+1]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
                            
	      //case 47 - j,l,k,n

	      if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b][a-1][bb][aa+1]
			     +penalty(i,j+2,ct1,data)+penalty(k-1,l+2,ct2,data)+penalty(m,n+2,ct3,data)
			     +edangle5(j+2,i,j+1,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)+edangle3(k-1,l+2,k,ct2,data)
			     +edangle5(n+2,m,n+1,ct3,data)
			     +2*gap){
		stack.push(j+2,i,b,a-1,bb,aa+1,v[i][j+2][b][a-1][bb][aa+1]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    //case 6 - k,m

	    if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b][a-1][bb][aa]
			   +penalty(i,j+1,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m-1,n+1,ct3,data)
			   +edangle3(k-1,l+1,k,ct2,data)+edangle3(m-1,n+1,m,ct3,data)
			   +gap){
	      stack.push(j+1,i,b,a-1,bb,aa,v[i][j+1][b][a-1][bb][aa]);
	      found = true;
	      en4=mine[j][b][bb];
	      i = j;
	      a = b;
	      aa = bb;							
	    }
                        
	    //case 48 - j,l,k,n,m

	    if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b][a-1][bb][aa]
			   +penalty(i,j+2,ct1,data)+penalty(k-1,l+2,ct2,data)+penalty(m-1,n+2,ct3,data)
			   +edangle5(j+2,i,j+1,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)+edangle3(k-1,l+2,k,ct2,data)
			   +edangle5(n+2,m-1,n+1,ct3,data)+edangle3(m-1,n+2,m,ct3,data)
			   +gap){
	      stack.push(j+2,i,b,a-1,bb,aa,v[i][j+2][b][a-1][bb][aa]);
	      found = true;
	      en4=mine[j][b][bb];
	      i = j;
	      a = b;
	      aa = bb;							
	    }
	  }
                    
	  if(a+1<2*maxsep+2) {
                        
	    if(bb>=1) {
                            
	      if(aa>=1) {
		//case 58 - j,i,l,m

		if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b][a+1][bb-1][aa-1]
			       +penalty(i-1,j+2,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m-1,n+1,ct3,data)
			       +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle5(l+2,k,l+1,ct2,data)
			       +edangle3(m-1,n+1,m,ct3,data)
			       +2*gap){
		  stack.push(j+2,i-1,b,a+1,bb-1,aa-1,v[i-1][j+2][b][a+1][bb-1][aa-1]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      //case 57 - j,i,l

	      if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b][a+1][bb-1][aa]
			     +penalty(i-1,j+2,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m,n+1,ct3,data)
			     +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle5(l+2,k,l+1,ct2,data)
			     +3*gap){
		stack.push(j+2,i-1,b,a+1,bb-1,aa,v[i-1][j+2][b][a+1][bb-1][aa]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(bb+1<2*maxsep+2) {
                            
	      if(aa>=1) {
		//case 20 - i,n,m

		if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b][a+1][bb+1][aa-1]
			       +penalty(i-1,j+1,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m-1,n+2,ct3,data)
			       +edangle3(i-1,j+1,i,ct1,data)+edangle3(m-1,n+2,m,ct3,data)+edangle5(n+2,m-1,n+1,ct3,data)
			       +3*gap){
		  stack.push(j+1,i-1,b,a+1,bb+1,aa-1,v[i-1][j+1][b][a+1][bb+1][aa-1]);
		  found = true;
		  en4=mine[j][b][bb];
		  i = j;
		  a = b;
		  aa = bb;									
		}
	      }
                            
	      //case 19 - i,n

	      if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b][a+1][bb+1][aa]
			     +penalty(i-1,j+1,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m,n+2,ct3,data)
			     +edangle3(i-1,j+1,i,ct1,data)+edangle5(n+2,m,n+1,ct3,data)
			     +4*gap){
		stack.push(j+1,i-1,b,a+1,bb+1,aa,v[i-1][j+1][b][a+1][bb+1][aa]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(aa>=1) {
	      //case 18 - i,m

	      if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b][a+1][bb][aa-1]
			     +penalty(i-1,j+1,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m-1,n+1,ct3,data)
			     +edangle3(i-1,j+1,i,ct1,data)+edangle3(m-1,n+1,m,ct3,data)
			     +gap){
		stack.push(j+1,i-1,b,a+1,bb,aa-1,v[i-1][j+1][b][a+1][bb][aa-1]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
                            
	      //case 60 - j,i,l,n,m

	      if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b][a+1][bb][aa-1]
			     +penalty(i-1,j+2,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m-1,n+2,ct3,data)
			     +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle5(l+2,k,l+1,ct2,data)
			     +edangle5(n+2,m-1,n+1,ct3,data)+edangle3(m-1,n+2,m,ct3,data)
			     +gap){
		stack.push(j+2,i-1,b,a+1,bb,aa-1,v[i-1][j+2][b][a+1][bb][aa-1]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    //case 17 - i

	    if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b][a+1][bb][aa]
			   +penalty(i-1,j+1,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m,n+1,ct3,data)
			   +edangle3(i-1,j+1,i,ct1,data)
			   +2*gap){
	      stack.push(j+1,i-1,b,a+1,bb,aa,v[i-1][j+1][b][a+1][bb][aa]);
	      found = true;
	      en4=mine[j][b][bb];
	      i = j;
	      a = b;
	      aa = bb;							
	    }
                        
	    //case 59 - j,i,l,n

	    if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b][a+1][bb][aa]
			   +penalty(i-1,j+2,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m,n+2,ct3,data)
			   +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle5(l+2,k,l+1,ct2,data)
			   +edangle5(n+2,m-1,n+1,ct3,data)
			   +2*gap){
	      stack.push(j+2,i-1,b,a+1,bb,aa,v[i-1][j+2][b][a+1][bb][aa]);
	      found = true;
	      en4=mine[j][b][bb];
	      i = j;
	      a = b;
	      aa = bb;							
	    }
	  }
                    
	  if(bb>=1) {
                        
	    if(aa>=1) {
	      //case 42 - j,l,m

	      if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b][a][bb-1][aa-1]
			     +penalty(i,j+2,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m-1,n+1,ct3,data)
			     +edangle5(j+2,i,j+1,ct1,data)+edangle5(l+2,k,l+1,ct2,data)+edangle3(m-1,n+1,m,ct3,data)
			     +3*gap){
		stack.push(j+2,i,b,a,bb-1,aa-1,v[i][j+2][b][a][bb-1][aa-1]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    if(aa+1<2*maxsep+2) {
	      //case 61 - j,i,l,k

	      if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b][a][bb-1][aa+1]
			     +penalty(i-1,j+2,ct1,data)+penalty(l+2,k-1,ct2,data)+penalty(m,n+1,ct3,data)
			     +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)
			     +edangle3(k-1,l+2,k,ct2,data)
			     +2*gap){
		stack.push(j+2,i-1,b,a,bb-1,aa+1,v[i-1][j+2][b][a][bb-1][aa+1]);
		found = true;
		en4=mine[j][b][bb];
		i = j;
		a = b;
		aa = bb;								
	      }
	    }
                        
	    //case 41 - j,l

	    if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b][a][bb-1][aa]
			   +penalty(i,j+2,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m,n+1,ct3,data)
			   +edangle5(j+2,i,j+1,ct1,data)+edangle5(l+2,k,l+1,ct2,data)
			   +gap){
	      stack.push(j+2,i,b,a,bb-1,aa,v[i][j+2][b][a][bb-1][aa]);
	      found = true;
	      en4=mine[j][b][bb];
	      i = j;
	      a = b;
	      aa = bb;							
	    }
                        
	    //case 62 - j,i,l,k,m

	    if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b][a][bb-1][aa]
			   +penalty(i-1,j+2,ct1,data)+penalty(l+2,k-1,ct2,data)+penalty(m-1,n+1,ct3,data)
			   +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)
			   +edangle3(k-1,l+2,k,ct2,data)+edangle3(m-1,n+2,m,ct3,data)
			   +gap){
	      stack.push(j+2,i-1,b,a,bb-1,aa,v[i-1][j+2][b][a][bb-1][aa]);
	      found = true;
	      en4=mine[j][b][bb];
	      i = j;
	      a = b;
	      aa = bb;							
	    }
	  }
                    
	  if(bb+1<2*maxsep+2) {
                        
	    if(aa>=1) {
	      //case 4 - n,m

							if(!found) if (en4==mine[j][b][bb]+v[i][j+1][b][a][bb+1][aa-1]
                            +penalty(i,j+1,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m-1,n+2,ct3,data)
                            +edangle3(m-1,n+2,m,ct3,data)+edangle5(n+2,m-1,n+1,ct3,data)
                            +4*gap){
                                stack.push(j+1,i,b,a,bb+1,aa-1,v[i][j+1][b][a][bb+1][aa-1]);
								found = true;
                                en4 = mine[j][b][bb];
                                i=j;
                                a=b;
                                aa=bb;								
                            }
                        }
                        
                        if(aa+1<2*maxsep+2) {
                            //case 23 - i,k,n

							if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b][a][bb+1][aa+1]
                            +penalty(i-1,j+1,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m,n+2,ct3,data)
                            +edangle3(k-1,l+1,k,ct2,data)+edangle3(i-1,j+1,i,ct1,data)+edangle5(n+2,m,n+1,ct3,data)
                            +3*gap){
                                stack.push(j+1,i-1,b,a,bb+1,aa+1,v[i-1][j+1][b][a][bb+1][aa+1]);
                                found = true;
                                en4=mine[j][b][bb];
                                i = j;
                                a = b;
                                aa = bb;								
                            }
                        }
                        
                        //case 3 - n

						if(!found) if (en4 == mine[j][b][bb]+v[i][j+1][b][a][bb+1][aa]
                        +penalty(i,j+1,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m,n+2,ct3,data)
                        +edangle5(n+2,m,n+1,ct3,data)
                        +2*gap) {                            
                            stack.push(j+1,i,b,a,bb+1,aa,v[i][j+1][b][a][bb+1][aa]);
							found = true;
                            en4 = mine[j][b][bb];
                            i=j;
                            a=b;
                            aa=bb;							
                        }
                        
                        //case 24 - i,k,n,m

						if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b][a][bb+1][aa]
                        +penalty(i-1,j+1,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m-1,n+2,ct3,data)
                        +edangle3(k-1,l+1,k,ct2,data)+edangle3(i-1,j+1,i,ct1,data)+edangle3(m-1,n+2,m,ct3,data)
                        +edangle5(n+2,m-1,n-1,ct3,data)
                        +2*gap){
                            stack.push(j+1,i-1,b,a,bb+1,aa,v[i-1][j+1][b][a][bb+1][aa]);
                            found = true;
                            en4=mine[j][b][bb];
                            i = j;
                            a = b;
                            aa = bb;							
                        }
                    }
                    
                    if(aa>=1) {
                        //case 2 - m

						if(!found) if (en4 == mine[j][b][bb]+v[i][j+1][b][a][bb][aa-1]
                        +penalty(i,j+1,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m-1,n+1,ct3,data)
                        +edangle3(m-1,n+1,m,ct3,data)
                        +2*gap) {                            
                            stack.push(j+1,i,b,a,bb,aa-1,v[i][j+1][b][a][bb][aa-1]);
							found = true;
                            en4 = mine[j][b][bb];
                            i=j;
                            a=b;
                            aa=bb;							
                        }
                        
                        //case 44 - j,l,n,m

						if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b][a][bb][aa-1]
                        +penalty(i,j+2,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m-1,n+2,ct3,data)
                        +edangle5(j+2,i,j+1,ct1,data)+edangle5(l+2,k,l+1,ct2,data)+edangle5(n+2,m-1,n+1,ct3,data)
                        +edangle3(m-1,n+2,m,ct3,data)
                        +2*gap){
                            stack.push(j+2,i,b,a,bb,aa-1,v[i][j+2][b][a][bb][aa-1]);
                            found = true;
                            en4=mine[j][b][bb];
                            i = j;
                            a = b;
                            aa = bb;							
                        }
                    }
                    
                    if(aa+1<2*maxsep+2) {
                        //case 21 - i,k

						if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b][a][bb][aa+1]
                        +penalty(i-1,j+1,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m,n+1,ct3,data)
                        +edangle3(k-1,l+1,k,ct2,data)+edangle3(i-1,j+1,i,ct1,data)
                        +gap){
                            stack.push(j+1,i-1,b,a,bb,aa+1,v[i-1][j+1][b][a][bb][aa+1]);
                            found = true;
                            en4=mine[j][b][bb];
                            i = j;
                            a = b;
                            aa = bb;							
                        }
                        
                        //case 63 - j,i,l,k,n

						if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b][a][bb][aa+1]
                        +penalty(i-1,j+2,ct1,data)+penalty(l+2,k-1,ct2,data)+penalty(m,n+2,ct3,data)
                        +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)
                        +edangle3(k-1,l+2,k,ct2,data)+edangle5(n+2,m-1,n+1,ct3,data)
                        +gap){
                            stack.push(j+2,i-1,b,a,bb,aa+1,v[i-1][j+2][b][a][bb][aa+1]);
                            found = true;
                            en4=mine[j][b][bb];
                            i = j;
                            a = b;
                            aa = bb;							
                        }
                    }
                    
                    //case 1

					if(!found) if (en4 == mine[j][b][bb]+v[i][j+1][b][a][bb][aa]
                    +penalty(i,j+1,ct1,data)+penalty(k,l+1,ct2,data)+penalty(m,n+1,ct3,data)) {
                        stack.push(j+1,i,b,a,bb,aa,v[i][j+1][b][a][bb][aa]);
                        found = true;
                        en4 = mine[j][b][bb];
                        i=j;
                        a=b;
                        aa=bb;						
                    }
                    
                    //case 22 - i,k,m

					if(!found) if (en4==mine[j][b][bb]+v[i-1][j+1][b][a][bb][aa]
                    +penalty(i-1,j+1,ct1,data)+penalty(k-1,l+1,ct2,data)+penalty(m-1,n+1,ct3,data)
                    +edangle3(k-1,l+1,k,ct2,data)+edangle3(i-1,j+1,i,ct1,data)+edangle3(m-1,n+1,m,ct3,data)){
                        stack.push(j+1,i-1,b,a,bb,aa,v[i-1][j+1][b][a][bb][aa]);
                        found = true;
                        en4=mine[j][b][bb];						
                        i = j;
                        a = b;
                        aa = bb;						
                    }
                    
                    //case 43 - j,l,n

					if(!found) if (en4==mine[j][b][bb]+v[i][j+2][b][a][bb][aa]
                    +penalty(i,j+2,ct1,data)+penalty(l+2,k,ct2,data)+penalty(m,n+2,ct3,data)
                    +edangle5(j+2,i,j+1,ct1,data)+edangle5(l+2,k,l+1,ct2,data)+edangle5(n+2,m,n+1,ct3,data)){
                        stack.push(j+2,i,b,a,bb,aa,v[i][j+2][b][a][bb][aa]);
                        found = true;
                        en4=mine[j][b][bb];
                        i = j;
                        a = b;
                        aa = bb;						
                    }
                    
                    //case 64 - j,i,l,k,n,m

					if(!found) if (en4==mine[j][b][bb]+v[i-1][j+2][b][a][bb][aa]
                    +penalty(i-1,j+2,ct1,data)+penalty(l+2,k-1,ct2,data)+penalty(m-1,n+2,ct3,data)
                    +edangle5(j+2,i-1,j+1,ct1,data)+edangle3(i-1,j+2,i,ct1,data)+edangle5(l+2,k-1,l+1,ct2,data)
                    +edangle3(k-1,l+2,k,ct2,data)+edangle3(m-1,n+2,m,ct3,data)+edangle5(n+2,m-1,n+1,ct3,data)){
                        stack.push(j+2,i-1,b,a,bb,aa,v[i-1][j+2][b][a][bb][aa]);
                        found = true;
                        en4=mine[j][b][bb];
                        i = j;
                        a = b;
                        aa = bb;						
                    }
                    
                    
                    
                }
                
            }
        }
        
        if(!found) {//bifuraction was never found--issue error
            cout << "Bifurcation not found!\n"<<flush;
            return;
        }
    }
    
    cout << "end of find exterior \n"<<flush;

    /*************************************** find interior *****************************************/
    // while stack not empty
	while (stack.pull(&i,&j, &a, &b, &aa, &bb, &en4)) {
        cout << "pulled "<<i<< " "<<j<<" "<<a<<" "<<b<<" "<<aa<<" "<<bb<<" "<<en4<<"\n"<<flush;	
        
        
        
        
        k = a+i-maxsep;
        l = b+j-maxsep;
        m = aa+k-maxsep;
        n = bb+l-maxsep;
        
        if (en4==v[j][i][a][b][aa][bb]) {
            //i and j are paired and aligned to a and b and also aa and bb
            ct1->basepr[1][i]=j;
            ct1->basepr[1][j]=i;
            ct2->basepr[1][k]=l;
            ct2->basepr[1][l]=k;
            ct3->basepr[1][m]=n;
            ct3->basepr[1][n]=m;
            alignment1[i]=k;
            alignment1[j]=l;
            alignment2[k]=m;
            alignment2[l]=n;
			cout << "Aligned "<<i<< " "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<n<<" "<<en4<<"\n"<<flush;	
            
            //now find the next pair:
            /************************************ traceback V1 & V2 ****************************************/
                nogaps=3*max(j-i,max(l-k,n-m))-j+i-l+k-n+m;	
				if (en4!=ehairpin(i,j,ct1,data)+ehairpin(k,l,ct2,data)+ehairpin(m,n,ct3,data)+gap*nogaps){ 
                //the fragment does not close hairpins, therefore, the internal
                // fragment does need to be characterized
                
                //first check for internal loop:
                found= false;
                
                //Now consider internal loops/one stacked and other not
                for (c=i+1;c<=i+maxloop&&c<j-minloop&&!found;c++) {
                    for (d=j-1;d>=j-maxloop&&d>c+minloop;d--) {
                        for (e=max(k+1,c-maxsep);e<=k+maxloop&&e<l-minloop&&(e-c)<=maxsep&&!found;e++) {
                            for (f=min(l-1,d+maxsep);f>=l-maxloop&&f>e+minloop&&(d-f)<=maxsep&&!found;f--) {
                                for (g=max(m+1,e-maxsep);g<=m+maxloop&&g<n-minloop&&(g-e)<=maxsep&&!found;g++) {
                                    for (h=min(n-1,f+maxsep);h>=n-maxloop&&h>g+minloop&&(f-h)<=maxsep&&!found;h--) {  
                                        
                                        /* Cases for Seq1 */
										if ( c==i+1 && d==j-1 ) { /* helical stacking */
											en1 = ebp(i,j,i+1,j-1,ct1, data);
										} 
										else { /* internal loop */
										    en1 = einternal(i,j,c,d,ct1,data);
									    }


										/* Cases for Seq2 */
										if ( e==k+1 && f==l-1 ) { /* helical stacking */	      
										    en2 = ebp(k,l,k+1,l-1,ct2, data);
										} 
										else { /* internal loop */
											en2 = einternal(k,l,e,f,ct2,data);
									    }


										/* Cases for Seq3 */															
										if ( g==m+1 && h==n-1 ) { /* helical stacking */
										    en3 = ebp(m,n,m+1,n-1,ct3, data); 
										} 
										else { /* internal loop */
											en3 = einternal(m,n,g,h,ct3,data);
										}
										nogaps=3*max(c-i,max(e-k,g-m))-c+i-e+k-g+m+
											   3*max(j-d,max(l-f,n-h))-j+d-l+f-n+h;
										if(en4== en1+en2+en3+v[d][c][e-c+maxsep][f-d+maxsep][g-e+maxsep][h-f+maxsep]+gap*nogaps){
											found = true;
											stack.push(c,d,e-c+maxsep,f-d+maxsep,g-e+maxsep,h-f+maxsep,v[d][c][e-c+maxsep][f-d+maxsep][g-e+maxsep][h-f+maxsep]);											
										}


										/* base pair insertion */
										else if (c==i+2 && d==j-2 && j-2>i+2 &&
											inc[ ct1->numseq[i+1] ][ ct1->numseq[j-1] ] &&
											inc[ ct1->numseq[i+2] ][ ct1->numseq[j-2] ] &&
											( ( e==k+1 && f==l-1 && inc[ ct2->numseq[e] ][ ct2->numseq[f] ] ) || 
											( g==m+1 && h==n-1 && inc[ ct3->numseq[g] ][ ct3->numseq[h] ] ) ) ) {
			    
											en1 = ebp( i,j,i+1,j-1,ct1,data ) + ebp( i+1,j-1,i+2,j-2,ct1,data );
											nogaps=3*max(c-i,max(e-k,g-m))-c+i-e+k-g+m+
												   3*max(j-d,max(l-f,n-h))-j+d-l+f-n+h;
											if(en4== en1+en2+en3+v[d][c][e-c+maxsep][f-d+maxsep][g-e+maxsep][h-f+maxsep]+gap*nogaps){
												//base pairs with single base pair insertion in ct1
												ct1->basepr[1][i+1]=j-1;
												ct1->basepr[1][j-1]=i+1;
												stack.push(c,d,e-c+maxsep,f-d+maxsep,g-e+maxsep,h-f+maxsep,v[d][c][e-c+maxsep][f-d+maxsep][g-e+maxsep][h-f+maxsep]);
												found = true;
											}
									    }

										/* base pair insertion */
										else if (e==k+2 && f==l-2 && l-2> k+2 && 
											inc[ ct2->numseq[k+1] ][ ct2->numseq[l-1] ] &&
											inc[ ct2->numseq[k+2] ][ ct2->numseq[l-2] ] && 
											( ( c==i+1 && d==j-1 && inc[ ct1->numseq[c] ][ ct1->numseq[d] ] ) || 
											( g==m+1 && h==n-1 && inc[ ct3->numseq[g] ][ ct3->numseq[h] ] ) ) ) {

											en2 = ebp( k,l,k+1,l-1,ct2,data ) + ebp( k+1,l-1,k+2,l-2,ct2,data );
											nogaps=3*max(c-i,max(e-k,g-m))-c+i-e+k-g+m+
												   3*max(j-d,max(l-f,n-h))-j+d-l+f-n+h;
											if(en4== en1+en2+en3+v[d][c][e-c+maxsep][f-d+maxsep][g-e+maxsep][h-f+maxsep]+gap*nogaps){
												//base pairs with single base pair insertion in ct2
												ct2->basepr[1][k+1]=l-1;
												ct2->basepr[1][l-1]=k+1;
												stack.push(c,d,e-c+maxsep,f-d+maxsep,g-e+maxsep,h-f+maxsep,v[d][c][e-c+maxsep][f-d+maxsep][g-e+maxsep][h-f+maxsep]);
												found = true;
											}

										}

										/* base pair insertion */
										else if (g==m+2 && h==n-2 && n-2>m+2 &&
											inc[ ct3->numseq[m+1] ][ ct3->numseq[n-1] ] &&
											inc[ ct3->numseq[m+2] ][ ct3->numseq[n-2] ] &&
											( ( c==i+1 && d==j-1 && inc[ ct1->numseq[c] ][ ct1->numseq[d] ] ) || 
											( e==k+1 && f==l-1 && inc[ ct2->numseq[e] ][ ct2->numseq[f] ] ) ) ) {

											en3 = ebp( m,n,m+1,n-1,ct3,data) + ebp( m+1,n-1,m+2,n-2,ct3,data );
											nogaps=3*max(c-i,max(e-k,g-m))-c+i-e+k-g+m+
												   3*max(j-d,max(l-f,n-h))-j+d-l+f-n+h;
											if(en4== en1+en2+en3+v[d][c][e-c+maxsep][f-d+maxsep][g-e+maxsep][h-f+maxsep]+gap*nogaps){
												//base pairs with single base pair insertion in ct3
												ct3->basepr[1][m+1]=n-1;
												ct3->basepr[1][n-1]=m+1;
												stack.push(c,d,e-c+maxsep,f-d+maxsep,g-e+maxsep,h-f+maxsep,v[d][c][e-c+maxsep][f-d+maxsep][g-e+maxsep][h-f+maxsep]);
												found = true;
											}

										} 



									}
								}
							}
						}
					}
				}
  							
                /*********************************** traceback V3 ***************************************/
                //Now consider multiloop:
                //junction closed by i-j pair aligned with a and b and also aa and bb
                //calculate the free energy of 3 fragments merged:
                for (c=i+minloop+1;c<j-minloop&&!found;c++) {
                    for (d=max(k+minloop+1,c-maxsep);d<l-minloop&&d<=c+maxsep&&!found;d++) {
                        e = d-c+maxsep;
                        for (f=max(m+minloop+1,d-maxsep);f<n-minloop&&f<=d+maxsep&&!found;f++) {
                            
                            constantclosure = penalty(i,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n,ct3,data);
                            
                            g = f-d+maxsep;
                            
                            //must consider whether an unpaired nucleotide is stacked
                            //is stacked onto each of the nucleotides in a helix
                            //There are 64 cases:
                            //consider rthem in this order (0 unstacked, 1 stacked)
                            //		i	j	k	l	n	n
                            //1		0	0	0	0	0	1
                            //2		0	0	0	0	1	0
                            //3		0	0	0	0	1	0
                            //4		0	0	0	0	1	1
                            //5		0	0	0	1	0	0
                            //6		0	0	0	1	0	1
                            //7		0	0	0	1	1	0
                            //8		0	0	0	1	1	1
                            //9		0	0	1	0	0	0
                            //10	0	0	1	0	0	1
                            //11	0	0	1	0	1	0
                            //12	0	0	1	0	1	1
                            //13	0	0	1	1	0	0
                            //14	0	0	1	1	0	1
                            //15	0	0	1	1	1	0
                            //16	0	0	1	1	1	1
                            //17	0	1	0	0	0	0
                            //18	0	1	0	0	0	1
                            //19	0	1	0	0	1	0
                            //20	0	1	0	0	1	1
                            //21	0	1	0	1	0	0
                            //22	0	1	0	1	0	1
                            //23	0	1	0	1	1	0
                            //24	0	1	0	1	1	1
                            //25	0	1	1	0	0	0
                            //26	0	1	1	0	0	1
                            //27	0	1	1	0	1	0
                            //28	0	1	1	0	1	1
                            //29	0	1	1	1	0	0
                            //30	0	1	1	1	0	1
                            //31	0	1	1	1	1	0
                            //32	0	1	1	1	1	1
                            //33	1	0	0	0	0	0
                            //34	1	0	0	0	0	1
                            //35	1	0	0	0	1	0
                            //36	1	0	0	0	1	1
                            //37	1	0	0	1	0	0
                            //38	1	0	0	1	0	1
                            //39	1	0	0	1	1	0
                            //40	1	0	0	1	1	1
                            //41	1	0	1	0	0	0
                            //42	1	0	1	0	0	1
                            //43	1	0	1	0	1	0
                            //44	1	0	1	0	1	1
                            //45	1	0	1	1	0	0
                            //46	1	0	1	1	0	1
                            //47	1	0	1	1	1	0
                            //48	1	0	1	1	1	1
                            //49	1	1	0	0	0	0
                            //50	1	1	0	0	0	1
                            //51	1	1	0	0	1	0
                            //52	1	1	0	0	1	1
                            //53	1	1	0	1	0	0
                            //54	1	1	0	1	0	1
                            //55	1	1	0	1	1	0
                            //56	1	1	0	1	1	1
                            //57	1	1	1	0	0	0
                            //58	1	1	1	0	0	1
                            //59	1	1	1	0	1	0
                            //60	1	1	1	0	1	1
                            //61	1	1	1	1	0	0
                            //62	1	1	1	1	0	1
                            //63	1	1	1	1	1	0
                            //64	1	1	1	1	1	1


                            
                            if(a>=1) {
                                
                                if(b>=1) {
                                    
                                    if(aa+1<2*maxsep+2) {
                                        
                                        if(bb+1<2*maxsep+2) {
                                            //case 39 - i,l,m
                                            if(en4==w[c][i+2][a-1][e][aa+1][g]+w[j-1][c+1][e][b-1][g][bb+1]
                                            +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                            +edangle3(i,j,i+1,ct1,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)
                                            +3*data->eparam[6]+3*gap
                                            &&!found) {
                                                stack.push(i+2,c,a-1,e,aa+1,g,w[c][i+2][a-1][e][aa+1][g]);

												stack.push(c+1,j-1,e,b-1,g,bb+1,w[j-1][c+1][e][b-1][g][bb+1]);

                                                found=true;
                                            }
                                        }
                                        
                                        //case 40 - i,l,m,n
                                        if(en4==w[c][i+2][a-1][e][aa+1][g]+w[j-1][c+1][e][b-1][g][bb]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle3(i,j,i+1,ct1,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                        +4*data->eparam[6]+2*gap
                                        &&!found) {
                                            stack.push(i+2,c,a-1,e,aa+1,g,w[c][i+2][a-1][e][aa+1][g]);

                                            stack.push(c+1,j-1,e,b-1,g,bb,w[j-1][c+1][e][b-1][g][bb]);
											 
                                            found=true;
                                        }
                                    }
                                    
                                    if(bb+1<2*maxsep+2) {
                                        //case 37 - i,l
                                        if(en4==w[c][i+2][a-1][e][aa][g]+w[j-1][c+1][e][b-1][g][bb+1]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle3(i,j,i+1,ct1,data)+edangle5(l,k,l-1,ct2,data)
                                        +2*data->eparam[6]+4*gap
                                        &&!found) {
                                            stack.push(i+2,c,a-1,e,aa,g,w[c][i+2][a-1][e][aa][g]);
											 
                                            stack.push(c+1,j-1,e,b-1,g,bb+1,w[j-1][c+1][e][b-1][g][bb+1]);
											 
                                            found=true;
                                        }
                                    }
                                    
                                    //case 38 - i,l,n
                                    if(en4==w[c][i+2][a-1][e][aa][g]+w[j-1][c+1][e][b-1][g][bb]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(i,j,i+1,ct1,data)+edangle5(l,k,l-1,ct2,data)+edangle5(n,m,n-1,ct3,data)
                                    +3*data->eparam[6]+3*gap
                                    &&!found) {
                                        stack.push(i+2,c,a-1,e,aa,g,w[c][i+2][a-1][e][aa][g]);
										 
                                        stack.push(c+1,j-1,e,b-1,g,bb,w[j-1][c+1][e][b-1][g][bb]);
										 
                                        found=true;
                                    }
                                }
                                
                                if(b+1<2*maxsep+2) {
                                    
                                    if(aa+1<2*maxsep+2) {
                                        
                                        if(bb>=1) {
                                            //case 52 - i,j,m,n
                                            if(en4==w[c][i+2][a-1][e][aa+1][g]+w[j-2][c+1][e][b+1][g][bb-1]
                                            +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                            +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                            +4*data->eparam[6]+2*gap
                                            &&!found) {
                                                stack.push(i+2,c,a-1,e,aa+1,g,w[c][i+2][a-1][e][aa+1][g]);
										
                                                stack.push(c+1,j-2,e,b+1,g,bb-1,w[j-2][c+1][e][b+1][g][bb-1]);
										
                                                found=true;
                                            }
                                        }
                                        
                                        //case 51 - i,j,m
                                        if(en4==w[c][i+2][a-1][e][aa+1][g]+w[j-2][c+1][e][b+1][g][bb]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle3(m,n,m+1,ct3,data)
                                        +3*data->eparam[6]+3*gap
                                        &&!found) {
                                            stack.push(i+2,c,a-1,e,aa+1,g,w[c][i+2][a-1][e][aa+1][g]);
										
                                            stack.push(c+1,j-2,e,b+1,g,bb,w[j-2][c+1][e][b+1][g][bb]);
										
                                            found=true;
                                        }
                                    }
                                    
                                    if(bb>=1) {
                                        //case 50 - i,j,n
                                        if(en4==w[c][i+2][a-1][e][aa][g]+w[j-2][c+1][e][b+1][g][bb-1]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle5(n,m,n-1,ct3,data)
                                        +3*data->eparam[6]+3*gap
                                        &&!found) {
                                            stack.push(i+2,c,a-1,e,aa,g,w[c][i+2][a-1][e][aa][g]);
											 
                                            stack.push(c+1,j-2,e,b+1,g,bb-1,w[j-2][c+1][e][b+1][g][bb-1]);
											 
                                            found=true;
                                        }
                                    }
                                    
                                    //case 49 - i,j
                                    if(en4==w[c][i+2][a-1][e][aa][g]+w[j-2][c+1][e][b+1][g][bb]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)
                                    +2*data->eparam[6]+4*gap
                                    &&!found) {
                                        stack.push(i+2,c,a-1,e,aa,g,w[c][i+2][a-1][e][aa][g]);
										 
                                        stack.push(c+1,j-2,e,b+1,g,bb,w[j-2][c+1][e][b+1][g][bb]);
										 
                                        found=true;
                                    }
                                }
                                
                                if(aa+1<2*maxsep+2) {
                                    
                                    if(bb>=1) {
                                        //case 36 - i,m,n
                                        if(en4==w[c][i+2][a-1][e][aa+1][g]+w[j-1][c+1][e][b][g][bb-1]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle3(i,j,i+1,ct1,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                        +3*data->eparam[6]+3*gap
                                        &&!found) {
                                            stack.push(i+2,c,a-1,e,aa+1,g,w[c][i+2][a-1][e][aa+1][g]);
										
                                            stack.push(c+1,j-1,e,b,g,bb-1,w[j-1][c+1][e][b][g][bb-1]);
										
                                            found=true;
                                        }
                                    }
                                    
                                    if(bb+1<2*maxsep+2) {
                                        //case 55 - i,j,l,m
                                        if(en4==w[c][i+2][a-1][e][aa+1][g]+w[j-2][c+1][e][b][g][bb+1]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)
                                        +4*data->eparam[6]+2*gap
                                        &&!found) {
                                            stack.push(i+2,c,a-1,e,aa+1,g,w[c][i+2][a-1][e][aa+1][g]);
										
                                            stack.push(c+1,j-2,e,b,g,bb+1,w[j-2][c+1][e][b][g][bb+1]);
										
                                            found=true;
                                        }
                                    }
                                    
                                    //case 35 - i,m
                                    if(en4==w[c][i+2][a-1][e][aa+1][g]+w[j-1][c+1][e][b][g][bb]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(i,j,i+1,ct1,data)+edangle3(m,n,m+1,ct3,data)
                                    +2*data->eparam[6]+gap
                                    &&!found) {
                                        stack.push(i+2,c,a-1,e,aa+1,g,w[c][i+2][a-1][e][aa+1][g]);
										 
                                        stack.push(c+1,j-1,e,b,g,bb,w[j-1][c+1][e][b][g][bb]);
										 
                                        found=true;
                                    }
                                    
                                    //case 56 - i,j,l,m,n
                                    if(en4==w[c][i+2][a-1][e][aa+1][g]+w[j-2][c+1][e][b][g][bb]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)+edangle3(n,m,n-1,ct3,data)
                                    +5*data->eparam[6]+gap
                                    &&!found) {
                                        stack.push(i+2,c,a-1,e,aa+1,g,w[c][i+2][a-1][e][aa+1][g]);
										 
                                        stack.push(c+1,j-2,e,b,g,bb,w[j-2][c+1][e][b][g][bb]);
										 
                                        found=true;
                                    }
                                }
                                
                                if(bb>=1) {
                                    //case 34 - i,n
                                    if(en4==w[c][i+2][a-1][e][aa][g]+w[j-1][c+1][e][b][g][bb-1]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(i,j,i+1,ct1,data)+edangle5(n,m,n-1,ct3,data)
                                    +2*data->eparam[6]+4*gap
                                    &&!found) {
                                        stack.push(i+2,c,a-1,e,aa,g,w[c][i+2][a-1][e][aa][g]);
										 
                                        stack.push(c+1,j-1,e,b,g,bb-1,w[j-1][c+1][e][b][g][bb-1]);
										 
                                        found=true;
                                    }
                                }
                                
                                if(bb+1<2*maxsep+2) {
                                    //case 53 - i,j,l
                                    if(en4==w[c][i+2][a-1][e][aa][g]+w[j-2][c+1][e][b][g][bb+1]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)
                                    +3*data->eparam[6]+3*gap
                                    &&!found) {
                                        stack.push(i+2,c,a-1,e,aa,g,w[c][i+2][a-1][e][aa][g]);
										 
                                        stack.push(c+1,j-2,e,b,g,bb+1,w[j-2][c+1][e][b][g][bb+1]);
										 
                                        found=true;
                                    }
                                }
                                
                                //case 33 - i
                                if(en4==w[c][i+2][a-1][e][aa][g]+w[j-1][c+1][e][b][g][bb]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle3(i,j,i+1,ct1,data)
                                +1*data->eparam[6]+2*gap
                                &&!found) {
                                    stack.push(i+2,c,a-1,e,aa,g,w[c][i+2][a-1][e][aa][g]);
									 
                                    stack.push(c+1,j-1,e,b,g,bb,w[j-1][c+1][e][b][g][bb]);
									 
                                    found=true;
                                }
                                
                                //case 54 - i,j,l,n
                                if(en4==w[c][i+2][a-1][e][aa][g]+w[j-2][c+1][e][b][g][bb]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+edangle5(n,m,n-1,ct3,data)
                                +4*data->eparam[6]+2*gap
                                &&!found) {
                                    stack.push(i+2,c,a-1,e,aa,g,w[c][i+2][a-1][e][aa][g]);
									 
                                    stack.push(c+1,j-2,e,b,g,bb,w[j-2][c+1][e][b][g][bb]);
									 
                                    found=true;
                                }
                            }
                            
                            if(a+1<2*maxsep+2) {
                                
                                if(b>=1) {
                                    
                                    if(aa>=1) {
                                        
                                        if(bb+1<2*maxsep+2) {
                                            //case 13 - k,l
                                            if(en4==w[c][i+1][a+1][e][aa-1][g]+w[j-1][c+1][e][b-1][g][bb+1]
                                            +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                            +edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)
                                            +2*data->eparam[6]+4*gap
                                            &&!found) {
                                                stack.push(i+1,c,a+1,e,aa-1,g,w[c][i+1][a+1][e][aa-1][g]);
									
                                                stack.push(c+1,j-1,e,b-1,g,bb+1,w[j-1][c+1][e][b-1][g][bb+1]);
									
                                                found=true;
                                            }
                                        }
                                        
                                        //case 14 - k,l,n
                                        if(en4==w[c][i+1][a+1][e][aa-1][g]+w[j-1][c+1][e][b-1][g][bb]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)+edangle5(n,m,n-1,ct3,data)
                                        +3*data->eparam[6]+3*gap
                                        &&!found) {
                                            stack.push(i+1,c,a+1,e,aa-1,g,w[c][i+1][a+1][e][aa-1][g]);
									
                                            stack.push(c+1,j-1,e,b-1,g,bb,w[j-1][c+1][e][b-1][g][bb]);
									
                                            found=true;
                                        }
                                    }
                                    
                                    if(bb+1<2*maxsep+2) {
                                        //case 15 - k,l,m
                                        if(en4==w[c][i+1][a+1][e][aa][g]+w[j-1][c+1][e][b-1][g][bb+1]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)
                                        +3*data->eparam[6]+3*gap
                                        &&!found) {
                                            stack.push(i+1,c,a+1,e,aa,g,w[c][i+1][a+1][e][aa][g]);
											
                                            stack.push(c+1,j-1,e,b-1,g,bb+1,w[j-1][c+1][e][b-1][g][bb+1]);
											
                                            found=true;
                                        }
                                    }
                                    
                                    //case 16 - k,l,m,n
                                    if(en4==w[c][i+1][a+1][e][aa][g]+w[j-1][c+1][e][b-1][g][bb]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                    +4*data->eparam[6]+2*gap
                                    &&!found) {
                                        stack.push(i+1,c,a+1,e,aa,g,w[c][i+1][a+1][e][aa][g]);
										 
                                        stack.push(c+1,j-1,e,b-1,g,bb,w[j-1][c+1][e][b-1][g][bb]);
										 
                                        found=true;
                                    }
                                }
                                
                                if(b+1<2*maxsep+2) {
                                    
                                    if(aa>=1) {
                                        
                                        if(bb>=1) {
                                            //case 26 - j,k,n
                                            if(en4==w[c][i+1][a+1][e][aa-1][g]+w[j-2][c+1][e][b+1][g][bb-1]
                                            +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                            +edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(n,m,n-1,ct3,data)
                                            +3*data->eparam[6]+3*gap
                                            &&!found) {
                                                stack.push(i+1,c,a+1,e,aa-1,g,w[c][i+1][a+1][e][aa-1][g]);
										
                                                stack.push(c+1,j-2,e,b+1,g,bb-1,w[j-2][c+1][e][b+1][g][bb-1]);
										
                                                found=true;
                                            }
                                        }
                                        
                                        //case 25 - j,k
                                        if(en4==w[c][i+1][a+1][e][aa-1][g]+w[j-2][c+1][e][b+1][g][bb]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)
                                        +2*data->eparam[6]+4*gap
                                        &&!found) {
                                            stack.push(i+1,c,a+1,e,aa-1,g,w[c][i+1][a+1][e][aa-1][g]);
										
                                            stack.push(c+1,j-2,e,b+1,g,bb,w[j-2][c+1][e][b+1][g][bb]);
										
                                            found=true;
                                        }
                                    }
                                    
                                    if(bb>=1) {
                                        //case 28 - j,k,m,n
                                        if(en4==w[c][i+1][a+1][e][aa][g]+w[j-2][c+1][e][b+1][g][bb-1]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                        +4*data->eparam[6]+2*gap
                                        &&!found) {
                                            stack.push(i+1,c,a+1,e,aa,g,w[c][i+1][a+1][e][aa][g]);
											 
                                            stack.push(c+1,j-2,e,b+1,g,bb-1,w[j-2][c+1][e][b+1][g][bb-1]);
											 
                                            found=true;
                                        }
                                    }
                                    
                                    //case 27 - j,k,m
                                    if(en4==w[c][i+1][a+1][e][aa][g]+w[j-2][c+1][e][b+1][g][bb]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle3(m,n,m+1,ct3,data)
                                    +3*data->eparam[6]+3*gap
                                    &&!found) {
                                        stack.push(i+1,c,a+1,e,aa,g,w[c][i+1][a+1][e][aa][g]);
										 
                                        stack.push(c+1,j-2,e,b+1,g,bb,w[j-2][c+1][e][b+1][g][bb]);
										 
                                        found=true;
                                    }
                                }
                                
                                if(aa>=1) {
                                    
                                    if(bb>=1) {
                                        //case 10 - k,n
                                        if(en4==w[c][i+1][a+1][e][aa-1][g]+w[j-1][c+1][e][b][g][bb-1]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle3(k,l,k+1,ct2,data)+edangle5(n,m,n-1,ct3,data)
                                        +2*data->eparam[6]+4*gap
                                        &&!found) {
                                            stack.push(i+1,c,a+1,e,aa-1,g,w[c][i+1][a+1][e][aa-1][g]);
										
                                            stack.push(c+1,j-1,e,b,g,bb-1,w[j-1][c+1][e][b][g][bb-1]);
										
                                            found=true;
                                        }
                                    }
                                    
                                    if(bb+1<2*maxsep+2) {
                                        //case 29 - j,k,l
                                        if(en4==w[c][i+1][a+1][e][aa-1][g]+w[j-2][c+1][e][b][g][bb+1]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)
                                        +3*data->eparam[6]+3*gap
                                        &&!found) {
                                            stack.push(i+1,c,a+1,e,aa-1,g,w[c][i+1][a+1][e][aa-1][g]);
										
                                            stack.push(c+1,j-2,e,b,g,bb+1,w[j-2][c+1][e][b][g][bb+1]);
										
                                            found=true;
                                        }
                                    }
                                    
                                    //case 9 - k
                                    if(en4==w[c][i+1][a+1][e][aa-1][g]+w[j-1][c+1][e][b][g][bb]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(k,l,k+1,ct2,data)
                                    +data->eparam[6]+2*gap
                                    &&!found) {
                                        stack.push(i+1,c,a+1,e,aa-1,g,w[c][i+1][a+1][e][aa-1][g]);
										
                                        stack.push(c+1,j-1,e,b,g,bb,w[j-1][c+1][e][b][g][bb]);
										
                                        found=true;
                                    }
                                    
                                    //case 30 - j,k,l,n
                                    if(en4==w[c][i+1][a+1][e][aa-1][g]+w[j-2][c+1][e][b][g][bb]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)+edangle5(n,m,n-1,ct3,data)
                                    +4*data->eparam[6]+2*gap
                                    &&!found) {
                                        stack.push(i+1,c,a+1,e,aa-1,g,w[c][i+1][a+1][e][aa-1][g]);
										
                                        stack.push(c+1,j-2,e,b,g,bb,w[j-2][c+1][e][b][g][bb]);
										
                                        found=true;
                                    }
                                }
                                
                                if(bb>=1) {
                                    //case 12 - k,m,n
                                    if(en4==w[c][i+1][a+1][e][aa][g]+w[j-1][c+1][e][b][g][bb-1]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(k,l,k+1,ct2,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                    +3*data->eparam[6]+3*gap
                                    &&!found) {
                                        stack.push(i+1,c,a+1,e,aa,g,w[c][i+1][a+1][e][aa][g]);
										
                                        stack.push(c+1,j-1,e,b,g,bb-1,w[j-1][c+1][e][b][g][bb-1]);
										
                                        found=true;
                                    }
                                }
                                
                                if(bb+1<2*maxsep+2) {
                                    //case 31 - j,k,l,m
                                    if(en4==w[c][i+1][a+1][e][aa][g]+w[j-2][c+1][e][b][g][bb+1]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)
                                    +4*data->eparam[6]+2*gap
                                    &&!found) {
                                        stack.push(i+1,c,a+1,e,aa,g,w[c][i+1][a+1][e][aa][g]);
										
                                        stack.push(c+1,j-2,e,b,g,bb+1,w[j-2][c+1][e][b][g][bb+1]);
										
                                        found=true;
                                    }
                                }
                                
                                //case 11 - k,m
                                if(en4==w[c][i+1][a+1][e][aa][g]+w[j-1][c+1][e][b][g][bb]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle3(k,l,k+1,ct2,data)+edangle3(m,n,m+1,ct3,data)
                                +2*data->eparam[6]+gap
                                &&!found) {
                                    stack.push(i+1,c,a+1,e,aa,g,w[c][i+1][a+1][e][aa][g]);
									 
                                    stack.push(c+1,j-1,e,b,g,bb,w[j-1][c+1][e][b][g][bb]);
									 
                                    found=true;
                                }
                                
                                //case 32 - j,k,l,m,n
                                if(en4==w[c][i+1][a+1][e][aa][g]+w[j-2][c+1][e][b][g][bb]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                +5*data->eparam[6]+gap
                                &&!found) {
                                    stack.push(i+1,c,a+1,e,aa,g,w[c][i+1][a+1][e][aa][g]);
									 
                                    stack.push(c+1,j-2,e,b,g,bb,w[j-2][c+1][e][b][g][bb]);
									 
                                    found=true;
                                }
                            }
                            
                            if(b>=1) {
                                
                                if(aa>=1) {
                                    
                                    if(bb+1<2*maxsep+2) {
                                        //case 45 - i,k,l
                                        if(en4==w[c][i+2][a][e][aa-1][g]+w[j-1][c+1][e][b-1][g][bb+1]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)
                                        +3*data->eparam[6]+3*gap
                                        &&!found) {
                                            stack.push(i+2,c,a,e,aa-1,g,w[c][i+2][a][e][aa-1][g]);
									
                                            stack.push(c+1,j-1,e,b-1,g,bb+1,w[j-1][c+1][e][b-1][g][bb+1]);
									
                                            found=true;
                                        }
                                    }
                                    
                                    //case 46 - i,k,l,n
                                    if(en4==w[c][i+2][a][e][aa-1][g]+w[j-1][c+1][e][b-1][g][bb]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)+edangle5(n,m,n-1,ct3,data)
                                    +4*data->eparam[6]+2*gap
                                    &&!found) {
                                        stack.push(i+2,c,a,e,aa-1,g,w[c][i+2][a][e][aa-1][g]);
									
                                        stack.push(c+1,j-1,e,b-1,g,bb,w[j-1][c+1][e][b-1][g][bb]);
									
                                        found=true;
                                    }
                                }
                                
                                if(aa+1<2*maxsep+2) {
                                    
                                    if(bb+1<2*maxsep+2) {
                                        //case 7 - l,m
                                        if(en4==w[c][i+1][a][e][aa+1][g]+w[j-1][c+1][e][b-1][g][bb+1]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)
                                        +2*data->eparam[6]+4*gap
                                        &&!found) {
                                            stack.push(i+1,c,a,e,aa+1,g,w[c][i+1][a][e][aa+1][g]);
											 
                                            stack.push(c+1,j-1,e,b-1,g,bb+1,w[j-1][c+1][e][b-1][g][bb+1]);
											 
                                            found=true;
                                        }
                                    }
                                    
                                    //case 8 - l,m,n
                                    if(en4==w[c][i+1][a][e][aa+1][g]+w[j-1][c+1][e][b-1][g][bb]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                    +3*data->eparam[6]+3*gap
                                    &&!found) {
                                        stack.push(i+1,c,a,e,aa+1,g,w[c][i+1][a][e][aa+1][g]);
										 
                                        stack.push(c+1,j-1,e,b-1,g,bb,w[j-1][c+1][e][b-1][g][bb]);
										 
                                        found=true;
                                    }
                                }
                                
                                if(bb+1<2*maxsep+2) {
                                    //case 5 - l
                                    if(en4==w[c][i+1][a][e][aa][g]+w[j-1][c+1][e][b-1][g][bb+1]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle5(l,k,l-1,ct2,data)
                                    +data->eparam[6]+2*gap
                                    &&!found) {
                                        stack.push(i+1,c,a,e,aa,g,w[c][i+1][a][e][aa][g]);
										 
                                        stack.push(c+1,j-1,e,b-1,g,bb+1,w[j-1][c+1][e][b-1][g][bb+1]);
										 
                                        found=true;
                                    }
                                    
                                    //case 47 - i,k,l,m
                                    if(en4==w[c][i+2][a][e][aa][g]+w[j-1][c+1][e][b-1][g][bb+1]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)
                                    +4*data->eparam[6]+2*gap
                                    &&!found) {
                                        stack.push(i+2,c,a,e,aa,g,w[c][i+2][a][e][aa][g]);
										 
                                        stack.push(c+1,j-1,e,b-1,g,bb+1,w[j-1][c+1][e][b-1][g][bb+1]);
										 
                                        found=true;
                                    }
                                }
                                
                                //case 6 - l,n
                                if(en4==w[c][i+1][a][e][aa][g]+w[j-1][c+1][e][b-1][g][bb]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle5(l,k,l-1,ct2,data)+edangle5(n,m,n-1,ct3,data)
                                +2*data->eparam[6]+gap
                                &&!found) {
                                    stack.push(i+1,c,a,e,aa,g,w[c][i+1][a][e][aa][g]);
									 
                                    stack.push(c+1,j-1,e,b-1,g,bb,w[j-1][c+1][e][b-1][g][bb]);
									 
                                    found=true;
                                }
                                
                                //case 48 - i,k,l,m,n
                                if(en4==w[c][i+2][a][e][aa][g]+w[j-1][c+1][e][b-1][g][bb]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                +5*data->eparam[6]+gap
                                &&!found) {
                                    stack.push(i+2,c,a,e,aa,g,w[c][i+2][a][e][aa][g]);
									 
                                    stack.push(c+1,j-1,e,b-1,g,bb,w[j-1][c+1][e][b-1][g][bb]);
									 
                                    found=true;
                                }
                            }
                            
                            if(b+1<2*maxsep+2) {
                                
                                if(aa>=1) {
                                    
                                    if(bb>=1) {
                                        //case 58 - i,j,k,n
                                        if(en4==w[c][i+2][a][e][aa-1][g]+w[j-2][c+1][e][b+1][g][bb-1]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(n,m,n-1,ct3,data)
                                        +4*data->eparam[6]+2*gap
                                        &&!found) {
                                            stack.push(i+2,c,a,e,aa-1,g,w[c][i+2][a][e][aa-1][g]);
									
                                            stack.push(c+1,j-2,e,b+1,g,bb-1,w[j-2][c+1][e][b+1][g][bb-1]);
									
                                            found=true;
                                        }
                                    }
                                    
                                    //case 57 - i,j,k
                                    if(en4==w[c][i+2][a][e][aa-1][g]+w[j-2][c+1][e][b+1][g][bb]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)
                                    +3*data->eparam[6]+3*gap
                                    &&!found) {
                                        stack.push(i+2,c,a,e,aa-1,g,w[c][i+2][a][e][aa-1][g]);
									
                                        stack.push(c+1,j-2,e,b+1,g,bb,w[j-2][c+1][e][b+1][g][bb]);
									
                                        found=true;
                                    }
                                }
                                
                                if(aa+1<2*maxsep+2) {
                                    
                                    if(bb>=1) {
                                        //case 20 - j,m,n
                                        if(en4==w[c][i+1][a][e][aa+1][g]+w[j-2][c+1][e][b+1][g][bb-1]
                                        +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                        +edangle5(j,i,j-1,ct1,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                        +3*data->eparam[6]+3*gap
                                        &&!found) {
                                            stack.push(i+1,c,a,e,aa+1,g,w[c][i+1][a][e][aa+1][g]);
											 
                                            stack.push(c+1,j-2,e,b+1,g,bb-1,w[j-2][c+1][e][b+1][g][bb-1]);
											 
                                            found=true;
                                        }
                                    }
                                    
                                    //case 19 - j,m
                                    if(en4==w[c][i+1][a][e][aa+1][g]+w[j-2][c+1][e][b+1][g][bb]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle5(j,i,j-1,ct1,data)+edangle3(m,n,m+1,ct3,data)
                                    +2*data->eparam[6]+4*gap
                                    &&!found) {
                                        stack.push(i+1,c,a,e,aa+1,g,w[c][i+1][a][e][aa+1][g]);
										 
                                        stack.push(c+1,j-2,e,b+1,g,bb,w[j-2][c+1][e][b+1][g][bb]);
										 
                                        found=true;
                                    }
                                }
                                
                                if(bb>=1) {
                                    //case 18 - j,n
                                    if(en4==w[c][i+1][a][e][aa][g]+w[j-2][c+1][e][b+1][g][bb-1]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle5(j,i,j-1,ct1,data)+edangle5(n,m,n-1,ct3,data)
                                    +2*data->eparam[6]+gap
                                    &&!found) {
                                        stack.push(i+1,c,a,e,aa,g,w[c][i+1][a][e][aa][g]);
										 
                                        stack.push(c+1,j-2,e,b+1,g,bb-1,w[j-2][c+1][e][b+1][g][bb-1]);
										 
                                        found=true;
                                    }
                                    
                                    //case 60 - i,j,k,m,n
                                    if(en4==w[c][i+2][a][e][aa][g]+w[j-2][c+1][e][b+1][g][bb-1]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                    +5*data->eparam[6]+gap
                                    &&!found) {
                                        stack.push(i+2,c,a,e,aa,g,w[c][i+2][a][e][aa][g]);
										 
                                        stack.push(c+1,j-2,e,b+1,g,bb-1,w[j-2][c+1][e][b+1][g][bb-1]);
										 
                                        found=true;
                                    }
                                }
                                
                                //case 17 - j
                                if(en4==w[c][i+1][a][e][aa][g]+w[j-2][c+1][e][b+1][g][bb]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle5(j,i,j-1,ct1,data)
                                +1*data->eparam[6]+2*gap
                                &&!found) {
                                    stack.push(i+1,c,a,e,aa,g,w[c][i+1][a][e][aa][g]);
									 
                                    stack.push(c+1,j-2,e,b+1,g,bb,w[j-2][c+1][e][b+1][g][bb]);
									 
                                    found=true;
                                }
                                
                                //case 59 - i,j,k,m
                                if(en4==w[c][i+2][a][e][aa][g]+w[j-2][c+1][e][b+1][g][bb]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle3(m,n,m+1,ct3,data)
                                +4*data->eparam[6]+2*gap
                                &&!found) {
                                    stack.push(i+2,c,a,e,aa,g,w[c][i+2][a][e][aa][g]);
									 
                                    stack.push(c+1,j-2,e,b+1,g,bb,w[j-2][c+1][e][b+1][g][bb]);
									 
                                    found=true;
                                }
                            }
                            
                            if(aa>=1) {
                                
                                if(bb>=1) {
                                    //case 42 - i,k,n
                                    if(en4==w[c][i+2][a][e][aa-1][g]+w[j-1][c+1][e][b][g][bb-1]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(n,m,n-1,ct3,data)
                                    +3*data->eparam[6]+3*gap
                                    &&!found) {
                                        stack.push(i+2,c,a,e,aa-1,g,w[c][i+2][a][e][aa-1][g]);
									
                                        stack.push(c+1,j-1,e,b,g,bb-1,w[j-1][c+1][e][b][g][bb-1]);
									
                                        found=true;
                                    }
                                }
                                
                                if(bb+1<2*maxsep+2) {
                                    //case 61 - i,j,k,l
                                    if(en4==w[c][i+2][a][e][aa-1][g]+w[j-2][c+1][e][b][g][bb+1]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)
                                    +4*data->eparam[6]+2*gap
                                    &&!found) {
                                        stack.push(i+2,c,a,e,aa-1,g,w[c][i+2][a][e][aa-1][g]);
									
                                        stack.push(c+1,j-2,e,b,g,bb+1,w[j-2][c+1][e][b][g][bb+1]);
									
                                        found=true;
                                    }
                                }
                                
                                //case 41 - i,k
                                if(en4==w[c][i+2][a][e][aa-1][g]+w[j-1][c+1][e][b][g][bb]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)
                                +2*data->eparam[6]+gap
                                &&!found) {
                                    stack.push(i+2,c,a,e,aa-1,g,w[c][i+2][a][e][aa-1][g]);
									
                                    stack.push(c+1,j-1,e,b,g,bb,w[j-1][c+1][e][b][g][bb]);
									
                                    found=true;
                                }
                                
                                //case 62 - i,j,k,l,n
                                if(en4==w[c][i+2][a][e][aa-1][g]+w[j-2][c+1][e][b][g][bb]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)+edangle5(n,m,n-1,ct3,data)
                                +5*data->eparam[6]+gap
                                &&!found) {
                                    stack.push(i+2,c,a,e,aa-1,g,w[c][i+2][a][e][aa-1][g]);
									
                                    stack.push(c+1,j-2,e,b,g,bb,w[j-2][c+1][e][b][g][bb]);
									
                                    found=true;
                                }
                            }
                            
                            if(aa+1<2*maxsep+2) {
                                
                                if(bb>=1) {
                                    //case 4 - m,n
                                    if(en4==w[c][i+1][a][e][aa+1][g]+w[j-1][c+1][e][b][g][bb-1]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                    +2*data->eparam[6]+4*gap
                                    &&!found) {
                                        stack.push(i+1,c,a,e,aa+1,g,w[c][i+1][a][e][aa+1][g]);
									
                                        stack.push(c+1,j-1,e,b,g,bb-1,w[j-1][c+1][e][b][g][bb-1]);
									
                                        found=true;
                                    }
                                }
                                
                                if(bb+1<2*maxsep+2) {
                                    //case 23 - j,l,m
                                    if(en4==w[c][i+1][a][e][aa+1][g]+w[j-2][c+1][e][b][g][bb+1]
                                    +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                    +edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)
                                    +3*data->eparam[6]+3*gap
                                    &&!found) {
                                        stack.push(i+1,c,a,e,aa+1,g,w[c][i+1][a][e][aa+1][g]);
									
                                        stack.push(c+1,j-2,e,b,g,bb+1,w[j-2][c+1][e][b][g][bb+1]);
									
                                        found=true;
                                    }
                                }
                                
                                //case 3 - m
                                if(en4==w[c][i+1][a][e][aa+1][g]+w[j-1][c+1][e][b][g][bb]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle3(m,n,m+1,ct3,data)
                                +1*data->eparam[6]+2*gap
                                &&!found) {
                                    stack.push(i+1,c,a,e,aa+1,g,w[c][i+1][a][e][aa+1][g]);
									 
                                    stack.push(c+1,j-1,e,b,g,bb,w[j-1][c+1][e][b][g][bb]);
									 
                                    found=true;
                                }
                                
                                //case 24 - j,l,m,n
                                if(en4==w[c][i+1][a][e][aa+1][g]+w[j-2][c+1][e][b][g][bb]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                +4*data->eparam[6]+2*gap
                                &&!found) {
                                    stack.push(i+1,c,a,e,aa+1,g,w[c][i+1][a][e][aa+1][g]);
									 
                                    stack.push(c+1,j-2,e,b,g,bb,w[j-2][c+1][e][b][g][bb]);
									 
                                    found=true;
                                }
                            }
                            
                            if(bb>=1) {
                                //case 2 - n
                                if(en4==w[c][i+1][a][e][aa][g]+w[j-1][c+1][e][b][g][bb-1]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle5(n,m,n-1,ct3,data)
                                +1*data->eparam[6]+2*gap
                                &&!found) {
                                    stack.push(i+1,c,a,e,aa,g,w[c][i+1][a][e][aa][g]);
									 
                                    stack.push(c+1,j-1,e,b,g,bb-1,w[j-1][c+1][e][b][g][bb-1]);
									 
                                    found=true;
                                }
                                
                                //case 44 - i,k,m,n
                                if(en4==w[c][i+2][a][e][aa][g]+w[j-1][c+1][e][b][g][bb-1]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                                +4*data->eparam[6]+2*gap
                                &&!found) {
                                    stack.push(i+2,c,a,e,aa,g,w[c][i+2][a][e][aa][g]);
									 
                                    stack.push(c+1,j-1,e,b,g,bb-1,w[j-1][c+1][e][b][g][bb-1]);
									 
                                    found=true;
                                }
                            }
                            
                            if(bb+1<2*maxsep+2) {
                                //case 21 - j,l
                                if(en4==w[c][i+1][a][e][aa][g]+w[j-2][c+1][e][b][g][bb+1]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)
                                +2*data->eparam[6]+gap
                                &&!found) {
                                    stack.push(i+1,c,a,e,aa,g,w[c][i+1][a][e][aa][g]);
									 
                                    stack.push(c+1,j-2,e,b,g,bb+1,w[j-2][c+1][e][b][g][bb+1]);
									 
                                    found=true;
                                }
                                
                                //case 63 - i,j,k,l,m
                                if(en4==w[c][i+2][a][e][aa][g]+w[j-2][c+1][e][b][g][bb+1]
                                +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                                +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)
                                +5*data->eparam[6]+gap
                                &&!found) {
                                    stack.push(i+2,c,a,e,aa,g,w[c][i+2][a][e][aa][g]);
									 
                                    stack.push(c+1,j-2,e,b,g,bb+1,w[j-2][c+1][e][b][g][bb+1]);
									 
                                    found=true;
                                }
                            }
                            
                            //case 1 -
                            if(en4==w[c][i+1][a][e][aa][g]+w[j-1][c+1][e][b][g][bb]
                            +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                            &&!found) {
                                stack.push(i+1,c,a,e,aa,g,w[c][i+1][a][e][aa][g]);
								 
                                stack.push(c+1,j-1,e,b,g,bb,w[j-1][c+1][e][b][g][bb]);
								 
                                found=true;
                            }                     
											
							
							//case 22 - j,l,n
                            if(en4==w[c][i+1][a][e][aa][g]+w[j-2][c+1][e][b][g][bb]
                            +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                            +edangle5(j,i,j-1,ct1,data)+edangle5(l,k,l-1,ct2,data)+edangle5(n,m,n-1,ct3,data)
                            +3*data->eparam[6]
                            &&!found) {
                                stack.push(i+1,c,a,e,aa,g,w[c][i+1][a][e][aa][g]);
								 
                                stack.push(c+1,j-2,e,b,g,bb,w[j-2][c+1][e][b][g][bb]);
								 
                                found=true;
                            }
                            
                            //case 43 - i,k,m
                            if(en4==w[c][i+2][a][e][aa][g]+w[j-1][c+1][e][b][g][bb]
                            +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                            +edangle3(i,j,i+1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle3(m,n,m+1,ct3,data)
                            +3*data->eparam[6]
                            &&!found) {
                                stack.push(i+2,c,a,e,aa,g,w[c][i+2][a][e][aa][g]);
								 
                                stack.push(c+1,j-1,e,b,g,bb,w[j-1][c+1][e][b][g][bb]);
								 
                                found=true;
                            }
                            
                            //case 64 - i,j,k,l,m,n
                            if(en4==w[c][i+2][a][e][aa][g]+w[j-2][c+1][e][b][g][bb]
                            +3*data->eparam[5]+3*data->eparam[10]+constantclosure
                            +edangle3(i,j,i+1,ct1,data)+edangle5(j,i,j-1,ct1,data)+edangle3(k,l,k+1,ct2,data)+edangle5(l,k,l-1,ct2,data)+edangle3(m,n,m+1,ct3,data)+edangle5(n,m,n-1,ct3,data)
                            +6*data->eparam[6]
                            &&!found) {
                                stack.push(i+2,c,a,e,aa,g,w[c][i+2][a][e][aa][g]);
								 
                                stack.push(c+1,j-2,e,b,g,bb,w[j-2][c+1][e][b][g][bb]);
								 
                                found=true;
                            }
                            
                        }
                    }
                }
                
                
                
            }
            
            if (!found) {
                cout << "Traceback error for V!\n"<<flush;
                return;
            }
            
            
        }
        // ******************start from W1->W2->W3
        else {//v[i][j][a][b][aa][bb]!=en4
            found = false;
            /************************************ traceback W1 ********************************/
            if(j-1>i+1) {
                
                if(a>=1) {
                    
                    if(b+1<2*maxsep+2) {
                        
                        if(aa+1<2*maxsep+2) {
                            
                            if(bb>=1) {
                                //case 52 - i,j,m,n
                                if(en4==w[j-1][i+1][a-1][b+1][aa+1][bb-1]
                                +4*data->eparam[6]+2*gap
                                &&!found) {
                                    found = true;
                                    stack.push(i+1,j-1,a-1,b+1,aa+1,bb-1,w[j-1][i+1][a-1][b+1][aa+1][bb-1]);
									
                                }
                            }
                            
                            //case 51 - i,j,m
                            if(en4==w[j-1][i+1][a-1][b+1][aa+1][bb]
                            +3*data->eparam[6]+3*gap
                            &&!found) {
                                found = true;
                                stack.push(i+1,j-1,a-1,b+1,aa+1,bb,w[j-1][i+1][a-1][b+1][aa+1][bb]);
							
                            }
                        }
                        
                        if(bb>=1) {
                            //case 50 - i,j,n
                            if(en4==w[j-1][i+1][a-1][b+1][aa][bb-1]
                            +3*data->eparam[6]+3*gap
                            &&!found) {
                                found = true;
                                stack.push(i+1,j-1,a-1,b+1,aa,bb-1,w[j-1][i+1][a-1][b+1][aa][bb-1]);
							
                            }
                        }
                        
                        //case 49 - i,j
                        if(en4==w[j-1][i+1][a-1][b+1][aa][bb]
                        +2*data->eparam[6]+4*gap
                        &&!found) {
                            found = true;
                            stack.push(i+1,j-1,a-1,b+1,aa,bb,w[j-1][i+1][a-1][b+1][aa][bb]);
							
                        }
                    }
                    
                    if(aa+1<2*maxsep+2) {
                        
                        if(bb+1<2*maxsep+2) {
                            //case 55 - i,j,l,m
                            if(en4==w[j-1][i+1][a-1][b][aa+1][bb+1]
                            +4*data->eparam[6]+2*gap
                            &&!found) {
                                found = true;
                                stack.push(i+1,j-1,a-1,b,aa+1,bb+1,w[j-1][i+1][a-1][b][aa+1][bb+1]);
							
                            }
                        }
                        
                        //case 56 - i,j,l,m,n
                        if(en4==w[j-1][i+1][a-1][b][aa+1][bb]
                        +5*data->eparam[6]+gap
                        &&!found) {
                            found = true;
                            stack.push(i+1,j-1,a-1,b,aa+1,bb,w[j-1][i+1][a-1][b][aa+1][bb]);
							 // cout << "Case 134:En=" << en4<<"\n"<<flush;
                        }
                    }
                    
                    if(bb+1<2*maxsep+2) {
                        //case 53 - i,j,l
                        if(en4==w[j-1][i+1][a-1][b][aa][bb+1]
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                            stack.push(i+1,j-1,a-1,b,aa,bb+1,w[j-1][i+1][a-1][b][aa][bb+1]);
							 
                        }
                    }
                    
                    //case 54 - i,j,l,n
                    if(en4==w[j-1][i+1][a-1][b][aa][bb]
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        stack.push(i+1,j-1,a-1,b,aa,bb,w[j-1][i+1][a-1][b][aa][bb]);
						 
                    }
                }
                
                if(b+1<2*maxsep+2) {
                    
                    if(aa>=1) {
                        
                        if(bb>=1) {
                            //case 58 - i,j,k,n
                            if(en4==w[j-1][i+1][a][b+1][aa-1][bb-1]
                            +4*data->eparam[6]+2*gap
                            &&!found) {
                                found = true;
                                stack.push(i+1,j-1,a,b+1,aa-1,bb-1,w[j-1][i+1][a][b+1][aa-1][bb-1]);
						
                            }
                        }
                        
                        //case 57 - i,j,k
                        if(en4==w[j-1][i+1][a][b+1][aa-1][bb]
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                            stack.push(i+1,j-1,a,b+1,aa-1,bb,w[j-1][i+1][a][b+1][aa-1][bb]);
						
                        }
                    }
                    
                    if(bb>=1) {
                        //case 60 - i,j,k,m,n
                        if(en4==w[j-1][i+1][a][b+1][aa][bb-1]
                        +5*data->eparam[6]+gap
                        &&!found) {
                            found = true;
                            stack.push(i+1,j-1,a,b+1,aa,bb-1,w[j-1][i+1][a][b+1][aa][bb-1]);
							 
                        }
                    }
                    
                    //case 59 - i,j,k,m
                    if(en4==w[j-1][i+1][a][b+1][aa][bb]
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        stack.push(i+1,j-1,a,b+1,aa,bb,w[j-1][i+1][a][b+1][aa][bb]);
						 
                    }
                }
                
                if(aa>=1) {
                    
                    if(bb+1<2*maxsep+2) {
                        //case 61 - i,j,k,l
                        if(en4==w[j-1][i+1][a][b][aa-1][bb+1]
                        +4*data->eparam[6]+2*gap
                        &&!found) {
                            found = true;
                            stack.push(i+1,j-1,a,b,aa-1,bb+1,w[j-1][i+1][a][b][aa-1][bb+1]);
						
                        }
                    }
                    
                    //case 62 - i,j,k,l,n
                    if(en4==w[j-1][i+1][a][b][aa-1][bb]
                    +5*data->eparam[6]+gap
                    &&!found) {
                        found = true;
                        stack.push(i+1,j-1,a,b,aa-1,bb,w[j-1][i+1][a][b][aa-1][bb]);
						
                    }
                }
                
                if(bb+1<2*maxsep+2) {
                    //case 63 - i,j,k,l,m
                    if(en4==w[j-1][i+1][a][b][aa][bb+1]
                    +5*data->eparam[6]+gap
                    &&!found) {
                        found = true;
                        stack.push(i+1,j-1,a,b,aa,bb+1,w[j-1][i+1][a][b][aa][bb+1]);
						
                    }
                }
                
                //case 64 - i,j,k,l,m,n
                if(en4==w[j-1][i+1][a][b][aa][bb]
                +6*data->eparam[6]
                &&!found) {
                    found = true;
                    stack.push(i+1,j-1,a,b,aa,bb,w[j-1][i+1][a][b][aa][bb]);
					 
                }
            }
            
            if(a>=1) {
                
                if(b>=1) {
                    
                    if(aa+1<2*maxsep+2) {
                        
                        if(bb+1<2*maxsep+2) {
                            //case 39 - i,l,m
                            if(en4==w[j][i+1][a-1][b-1][aa+1][bb+1]
                            +3*data->eparam[6]+3*gap
                            &&!found) {
                                found = true;
                                stack.push(i+1,j,a-1,b-1,aa+1,bb+1,w[j][i+1][a-1][b-1][aa+1][bb+1]);
								 
                            }
                        }
                        
                        //case 40 - i,l,m,n
                        if(en4==w[j][i+1][a-1][b-1][aa+1][bb]
                        +4*data->eparam[6]+2*gap
                        &&!found) {
                            found = true;
                            stack.push(i+1,j,a-1,b-1,aa+1,bb,w[j][i+1][a-1][b-1][aa+1][bb]);
							 
                        }
                    }
                    
                    if(bb+1<2*maxsep+2) {
                        //case 37 - i,l
                        if(en4==w[j][i+1][a-1][b-1][aa][bb+1]
                        +2*data->eparam[6]+4*gap
                        &&!found) {
                            found = true;
                            stack.push(i+1,j,a-1,b-1,aa,bb+1,w[j][i+1][a-1][b-1][aa][bb+1]);
							 
                        }
                    }
                    
                    //case 38 - i,l,n
                    if(en4==w[j][i+1][a-1][b-1][aa][bb]
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                        stack.push(i+1,j,a-1,b-1,aa,bb,w[j][i+1][a-1][b-1][aa][bb]);
						 
                    }
                }
                
                if(aa+1<2*maxsep+2) {
                    
                    if(bb>=1) {
                        //case 36 - i,m,n
                        if(en4==w[j][i+1][a-1][b][aa+1][bb-1]
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                            stack.push(i+1,j,a-1,b,aa+1,bb-1,w[j][i+1][a-1][b][aa+1][bb-1]);
						
                        }
                    }
                    
                    //case 35 - i,m
                    if(en4==w[j][i+1][a-1][b][aa+1][bb]
                    +2*data->eparam[6]+gap
                    &&!found) {
                        found = true;
                        stack.push(i+1,j,a-1,b,aa+1,bb,w[j][i+1][a-1][b][aa+1][bb]);
						
                    }
                }
                
                if(bb>=1) {
                    //case 34 - i,n
                    if(en4==w[j][i+1][a-1][b][aa][bb-1]
                    +2*data->eparam[6]+4*gap
                    &&!found) {
                        found = true;
                        stack.push(i+1,j,a-1,b,aa,bb-1,w[j][i+1][a-1][b][aa][bb-1]);
						
                    }
                }
                
                //case 33 - i
                if(en4==w[j][i+1][a-1][b][aa][bb]
                +data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    stack.push(i+1,j,a-1,b,aa,bb,w[j][i+1][a-1][b][aa][bb]);
					 
                }
            }
            
            if(a+1<2*maxsep+2) {
                
                if(b>=1) {
                    
                    if(aa>=1) {
                        
                        if(bb+1<2*maxsep+2) {
                            //case 13 - k,l
                            if(en4==w[j][i][a+1][b-1][aa-1][bb+1]
                            +2*data->eparam[6]+4*gap
                            &&!found) {
                                found = true;
                                stack.push(i,j,a+1,b-1,aa-1,bb+1,w[j][i][a+1][b-1][aa-1][bb+1]);
					
                            }
                        }
                        
                        //case 14 - k,l,n
                        if(en4==w[j][i][a+1][b-1][aa-1][bb]
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                            stack.push(i,j,a+1,b-1,aa-1,bb,w[j][i][a+1][b-1][aa-1][bb]);
					
                        }
                    }
                    
                    if(bb+1<2*maxsep+2) {
                        //case 15 - k,l,m
                        if(en4==w[j][i][a+1][b-1][aa][bb+1]
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                            stack.push(i,j,a+1,b-1,aa,bb+1,w[j][i][a+1][b-1][aa][bb+1]);
					
                        }
                    }
                    
                    //case 16 - k,l,m,n
                    if(en4==w[j][i][a+1][b-1][aa][bb]
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        stack.push(i,j,a+1,b-1,aa,bb,w[j][i][a+1][b-1][aa][bb]);
					
                    }
                }
                
                if(b+1<2*maxsep+2) {
                    
                    if(aa>=1) {
                        
                        if(bb>=1) {
                            //case 26 - j,k,n
                            if(en4==w[j-1][i][a+1][b+1][aa-1][bb-1]
                            +3*data->eparam[6]+3*gap
                            &&!found) {
                                found = true;
                                stack.push(i,j-1,a+1,b+1,aa-1,bb-1,w[j-1][i][a+1][b+1][aa-1][bb-1]);
					
                            }
                        }
                        
                        //case 25 - j,k
                        if(en4==w[j-1][i][a+1][b+1][aa-1][bb]
                        +2*data->eparam[6]+4*gap
                        &&!found) {
                            found = true;
                            stack.push(i,j-1,a+1,b+1,aa-1,bb,w[j-1][i][a+1][b+1][aa-1][bb]);
					
                        }
                    }
                    
                    if(bb>=1) {
                        //case 28 - j,k,m,n
                        if(en4==w[j-1][i][a+1][b+1][aa][bb-1]
                        +4*data->eparam[6]+2*gap
                        &&!found) {
                            found = true;
                            stack.push(i,j-1,a+1,b+1,aa,bb-1,w[j-1][i][a+1][b+1][aa][bb-1]);
							 
                        }
                    }
                    
                    //case 27 - j,k,m
                    if(en4==w[j-1][i][a+1][b+1][aa][bb]
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                        stack.push(i,j-1,a+1,b+1,aa,bb,w[j-1][i][a+1][b+1][aa][bb]);
						 
                    }
                }
                
                if(aa>=1) {
                    
                    if(bb>=1) {
                        //case 10 - k,n
                        if(en4==w[j][i][a+1][b][aa-1][bb-1]
                        +2*data->eparam[6]+4*gap
                        &&!found) {
                            found = true;
                            stack.push(i,j,a+1,b,aa-1,bb-1,w[j][i][a+1][b][aa-1][bb-1]);
						
                        }
                    }
                    
                    if(bb+1<2*maxsep+2) {
                        //case 29 - j,k,l
                        if(en4==w[j-1][i][a+1][b][aa-1][bb+1]
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                            stack.push(i,j-1,a+1,b,aa-1,bb+1,w[j-1][i][a+1][b][aa-1][bb+1]);
						
                        }
                    }
                    
                    //case 9 - k
                    if(en4==w[j][i][a+1][b][aa-1][bb]
                    +data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        stack.push(i,j,a+1,b,aa-1,bb,w[j][i][a+1][b][aa-1][bb]);
						
                    }
                    
                    //case 30 - j,k,l,n
                    if(en4==w[j-1][i][a+1][b][aa-1][bb]
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        stack.push(i,j-1,a+1,b,aa-1,bb,w[j-1][i][a+1][b][aa-1][bb]);
						
                    }
                }
                
                if(bb>=1) {
                    //case 12 - k,m,n
                    if(en4==w[j][i][a+1][b][aa][bb-1]
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                        stack.push(i,j,a+1,b,aa,bb-1,w[j][i][a+1][b][aa][bb-1]);
						 
                    }
                }
                
                if(bb+1<2*maxsep+2) {
                    //case 31 - j,k,l,m
                    if(en4==w[j-1][i][a+1][b][aa][bb+1]
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        stack.push(i,j-1,a+1,b,aa,bb+1,w[j-1][i][a+1][b][aa][bb+1]);
						 
                    }
                }
                
                //case 11 - k,m
                if(en4==w[j][i][a+1][b][aa][bb]
                +2*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    stack.push(i,j,a+1,b,aa,bb,w[j][i][a+1][b][aa][bb]);
					 
                }
                
                //case 32 - j,k,l,m,n
                if(en4==w[j-1][i][a+1][b][aa][bb]
                +5*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    stack.push(i,j-1,a+1,b,aa,bb,w[j-1][i][a+1][b][aa][bb]);
					 
                }
            }
            
            if(b>=1) {
                
                if(aa>=1) {
                    
                    if(bb+1<2*maxsep+2) {
                        //case 45 - i,k,l
                        if(en4==w[j][i+1][a][b-1][aa-1][bb+1]
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                            stack.push(i+1,j,a,b-1,aa-1,bb+1,w[j][i+1][a][b-1][aa-1][bb+1]);
					
                        }
                    }
                    
                    //case 46 - i,k,l,n
                    if(en4==w[j][i+1][a][b-1][aa-1][bb]
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        stack.push(i+1,j,a,b-1,aa-1,bb,w[j][i+1][a][b-1][aa-1][bb]);
						
                    }
                }
                
                if(aa+1<2*maxsep+2) {
                    
                    if(bb+1<2*maxsep+2) {
                        //case 7 - l,m
                        if(en4==w[j][i][a][b-1][aa+1][bb+1]
                        +2*data->eparam[6]+4*gap
                        &&!found) {
                            found = true;
                            stack.push(i,j,a,b-1,aa+1,bb+1,w[j][i][a][b-1][aa+1][bb+1]);
						
                        }
                    }
                    
                    //case 8 - l,m,n
                    if(en4==w[j][i][a][b-1][aa+1][bb]
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                        stack.push(i,j,a,b-1,aa+1,bb,w[j][i][a][b-1][aa+1][bb]);
						
                    }
                }
                
                if(bb+1<2*maxsep+2) {
                    //case 5 - l
                    if(en4==w[j][i][a][b-1][aa][bb+1]
                    +data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        stack.push(i,j,a,b-1,aa,bb+1,w[j][i][a][b-1][aa][bb+1]);
						
                    }
                    
                    //case 47 - i,k,l,m
                    if(en4==w[j][i+1][a][b-1][aa][bb+1]
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        stack.push(i+1,j,a,b-1,aa,bb+1,w[j][i+1][a][b-1][aa][bb+1]);
						
                    }
                }
                
                //case 6 - l,n
                if(en4==w[j][i][a][b-1][aa][bb]
                +2*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    stack.push(i,j,a,b-1,aa,bb,w[j][i][a][b-1][aa][bb]);
					 
                }
                
                //case 48 - i,k,l,m,n
                if(en4==w[j][i+1][a][b-1][aa][bb]
                +5*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    stack.push(i+1,j,a,b-1,aa,bb,w[j][i+1][a][b-1][aa][bb]);
					 
                }
            }
            
            if(b+1<2*maxsep+2) {
                
                if(aa+1<2*maxsep+2) {
                    
                    if(bb>=1) {
                        //case 20 - j,m,n
                        if(en4==w[j-1][i][a][b+1][aa+1][bb-1]
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                            stack.push(i,j-1,a,b+1,aa+1,bb-1,w[j-1][i][a][b+1][aa+1][bb-1]);
					
                        }
                    }
                    
                    //case 19 - j,m
                    if(en4==w[j-1][i][a][b+1][aa+1][bb]
                    +2*data->eparam[6]+4*gap
                    &&!found) {
                        found = true;
                        stack.push(i,j-1,a,b+1,aa+1,bb,w[j-1][i][a][b+1][aa+1][bb]);
					
                    }
                }
                
                if(bb>=1) {
                    //case 18 - j,n
                    if(en4==w[j-1][i][a][b+1][aa][bb-1]
                    +2*data->eparam[6]+gap
                    &&!found) {
                        found = true;
                        stack.push(i,j-1,a,b+1,aa,bb-1,w[j-1][i][a][b+1][aa][bb-1]);
					
                    }
                }
                
                //case 17 - j
                if(en4==w[j-1][i][a][b+1][aa][bb]
                +data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    stack.push(i,j-1,a,b+1,aa,bb,w[j-1][i][a][b+1][aa][bb]);
					
                }
            }
            
            if(aa>=1) {
                
                if(bb>=1) {
                    //case 42 - i,k,n
                    if(en4==w[j][i+1][a][b][aa-1][bb-1]
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                        stack.push(i+1,j,a,b,aa-1,bb-1,w[j][i+1][a][b][aa-1][bb-1]);
					
                    }
                }
                
                //case 41 - i,k
                if(en4==w[j][i+1][a][b][aa-1][bb]
                +2*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    stack.push(i+1,j,a,b,aa-1,bb,w[j][i+1][a][b][aa-1][bb]);
					
                }
            }
            
            if(aa+1<2*maxsep+2) {
                
                if(bb>=1) {
                    //case 4 - m,n
                    if(en4==w[j][i][a][b][aa+1][bb-1]
                    +2*data->eparam[6]+4*gap
                    &&!found) {
                        found = true;
                        stack.push(i,j,a,b,aa+1,bb-1,w[j][i][a][b][aa+1][bb-1]);
					
                    }
                }
                
                if(bb+1<2*maxsep+2) {
                    //case 23 - j,l,m
                    if(en4==w[j-1][i][a][b][aa+1][bb+1]
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                        stack.push(i,j-1,a,b,aa+1,bb+1,w[j-1][i][a][b][aa+1][bb+1]);
					
                    }
                }
                
                //case 3 - m
                if(en4==w[j][i][a][b][aa+1][bb]
                +data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    stack.push(i,j,a,b,aa+1,bb,w[j][i][a][b][aa+1][bb]);
					
                }
                
                //case 24 - j,l,m,n
                if(en4==w[j-1][i][a][b][aa+1][bb]
                +4*data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    stack.push(i,j-1,a,b,aa+1,bb,w[j-1][i][a][b][aa+1][bb]);
					
                }
            }
            
            if(bb>=1) {
                //case 2 - n
                if(en4==w[j][i][a][b][aa][bb-1]
                +data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    stack.push(i,j,a,b,aa,bb-1,w[j][i][a][b][aa][bb-1]);
					
                }
                
                //case 44 - i,k,m,n
                if(en4==w[j][i+1][a][b][aa][bb-1]
                +4*data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    stack.push(i+1,j,a,b,aa,bb-1,w[j][i+1][a][b][aa][bb-1]);
					
                }
            }
            
            if(bb+1<2*maxsep+2) {
                //case 21 - j,l
                if(en4==w[j-1][i][a][b][aa][bb+1]
                +2*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    stack.push(i,j-1,a,b,aa,bb+1,w[j-1][i][a][b][aa][bb+1]);
					
                }
            }
            
            //case 22 - j,l,n
            if(en4==w[j-1][i][a][b][aa][bb]
            +3*data->eparam[6]
            &&!found) {
                found = true;
                stack.push(i,j-1,a,b,aa,bb,w[j-1][i][a][b][aa][bb]);
				
            }
            
            //case 43 - i,k,m
            if(en4==w[j][i+1][a][b][aa][bb]
            +3*data->eparam[6]
            &&!found) {
                found = true;
                stack.push(i+1,j,a,b,aa,bb,w[j][i+1][a][b][aa][bb]);
				
            }
            
            
            /************************** traceback W2 *************************************/
            
            
            
            
            if(a>=1) {
                
                if(b>=1) {
                    
                    if(aa+1<2*maxsep+2) {
                        
                        if(bb+1<2*maxsep+2) {
                            //case 39 - i,l,m
                            if(en4==v[j][i+1][a-1][b-1][aa+1][bb+1]
                            +3*data->eparam[10]
                            +penalty(i+1,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n,ct3,data)
                            +edangle5(i+1,j,i,ct1,data)+edangle3(l-1,k,l,ct2,data)+edangle5(m+1,n,m,ct3,data)
                            +3*data->eparam[6]+3*gap
                            &&!found) {
                                found = true;
                
								stack.push(i+1,j,a-1,b-1,aa+1,bb+1,v[j][i+1][a-1][b-1][aa+1][bb+1]);
                            }
                        }
                        
                        //case 40 - i,l,m,n
                        if(en4==v[j][i+1][a-1][b-1][aa+1][bb]
                        +3*data->eparam[10]
                        +penalty(i+1,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
                        +edangle5(i+1,j,i,ct1,data)+edangle3(l-1,k,l,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                        +4*data->eparam[6]+2*gap
                        &&!found) {
                            found = true;
                
							stack.push(i+1,j,a-1,b-1,aa+1,bb,v[j][i+1][a-1][b-1][aa+1][bb]);
                        }
                    }
                    
                    if(bb+1<2*maxsep+2) {
                        //case 37 - i,l
                        if(en4==v[j][i+1][a-1][b-1][aa][bb+1]
                        +3*data->eparam[10]
                        +penalty(i+1,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n,ct3,data)
                        +edangle5(i+1,j,i,ct1,data)+edangle3(l-1,k,l,ct2,data)
                        +2*data->eparam[6]+4*gap
                        &&!found) {
                            found = true;
                            
							stack.push(i+1,j,a-1,b-1,aa,bb+1,v[j][i+1][a-1][b-1][aa][bb+1]);
                        }
                    }
                    
                    //case 38 - i,l,n
                    if(en4==v[j][i+1][a-1][b-1][aa][bb]
                    +3*data->eparam[10]
                    +penalty(i+1,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n-1,ct3,data)
                    +edangle5(i+1,j,i,ct1,data)+edangle3(l-1,k,l,ct2,data)+edangle3(n-1,m,n,ct3,data)
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                        
						stack.push(i+1,j,a-1,b-1,aa,bb,v[j][i+1][a-1][b-1][aa][bb]);
                    }
                }
                
                if(b+1<2*maxsep+2) {
                    
                    if(aa+1<2*maxsep+2) {
                        
                        if(bb>=1) {
                            //case 52 - i,j,m,n
                            if(en4==v[j-1][i+1][a-1][b+1][aa+1][bb-1]
                            +3*data->eparam[10]
                            +penalty(i+1,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n-1,ct3,data)
                            +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                            +4*data->eparam[6]+2*gap
                            &&!found) {
                                found = true;
                        
								stack.push(i+1,j-1,a-1,b+1,aa+1,bb-1,v[j-1][i+1][a-1][b+1][aa+1][bb-1]);
                            }
                        }
                        
                        //case 51 - i,j,m
                        if(en4==v[j-1][i+1][a-1][b+1][aa+1][bb]
                        +3*data->eparam[10]
                        +penalty(i+1,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n,ct3,data)
                        +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle5(m+1,n,m,ct3,data)
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                        
							stack.push(i+1,j-1,a-1,b+1,aa+1,bb,v[j-1][i+1][a-1][b+1][aa+1][bb]);
                        }
                    }
                    
                    if(bb>=1) {
                        //case 50 - i,j,n
                        if(en4==v[j-1][i+1][a-1][b+1][aa][bb-1]
                        +3*data->eparam[10]
                        +penalty(i+1,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n-1,ct3,data)
                        +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle3(n-1,m,n,ct3,data)
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                         
							stack.push(i+1,j-1,a-1,b+1,aa,bb-1,v[j-1][i+1][a-1][b+1][aa][bb-1]);
                        }
                    }
                    
                    //case 49 - i,j
                    if(en4==v[j-1][i+1][a-1][b+1][aa][bb]
                    +3*data->eparam[10]
                    +penalty(i+1,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n,ct3,data)
                    +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)
                    +2*data->eparam[6]+4*gap
                    &&!found) {
                        found = true;
                        
						stack.push(i+1,j-1,a-1,b+1,aa,bb,v[j-1][i+1][a-1][b+1][aa][bb]);
                    }
                }
                
                if(aa+1<2*maxsep+2) {
                    
                    if(bb>=1) {
                        //case 36 - i,m,n
                        if(en4==v[j][i+1][a-1][b][aa+1][bb-1]
                        +3*data->eparam[10]
                        +penalty(i+1,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n-1,ct3,data)
                        +edangle5(i+1,j,i,ct1,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                        
							stack.push(i+1,j,a-1,b,aa+1,bb-1,v[j][i+1][a-1][b][aa+1][bb-1]);
                        }
                    }
                    
                    if(bb+1<2*maxsep+2) {
                        //case 55 - i,j,l,m
                        if(en4==v[j-1][i+1][a-1][b][aa+1][bb+1]
                        +3*data->eparam[10]
                        +penalty(i+1,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n,ct3,data)
                        +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle3(l-1,k,l,ct2,data)+edangle5(m+1,n,m,ct3,data)
                        +4*data->eparam[6]+2*gap
                        &&!found) {
                            found = true;
                        
							stack.push(i+1,j-1,a-1,b,aa+1,bb+1,v[j-1][i+1][a-1][b][aa+1][bb+1]);
                        }
                    }
                    
                    //case 35 - i,m
                    if(en4==v[j][i+1][a-1][b][aa+1][bb]
                    +3*data->eparam[10]
                    +penalty(i+1,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n,ct3,data)
                    +edangle5(i+1,j,i,ct1,data)+edangle5(m+1,n,m,ct3,data)
                    +2*data->eparam[6]+gap
                    &&!found) {
                        found = true;
                       
						stack.push(i+1,j,a-1,b,aa+1,bb,v[j][i+1][a-1][b][aa+1][bb]);
                    }
                    
                    //case 56 - i,j,l,m,n
                    if(en4==v[j-1][i+1][a-1][b][aa+1][bb]
                    +3*data->eparam[10]
                    +penalty(i+1,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
                    +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle3(l-1,k,l,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                    +5*data->eparam[6]+gap
                    &&!found) {
                        found = true;
                       
						stack.push(i+1,j-1,a-1,b,aa+1,bb,v[j-1][i+1][a-1][b][aa+1][bb]);
                    }
                }
                
                if(bb>=1) {
                    //case 34 - i,n
                    if(en4==v[j][i+1][a-1][b][aa][bb-1]
                    +3*data->eparam[10]
                    +penalty(i+1,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n-1,ct3,data)
                    +edangle5(i+1,j,i,ct1,data)+edangle3(n-1,m,n,ct3,data)
                    +2*data->eparam[6]+4*gap
                    &&!found) {
                        found = true;
                       
						stack.push(i+1,j,a-1,b,aa,bb-1,v[j][i+1][a-1][b][aa][bb-1]);
                    }
                }
                
                if(bb+1<2*maxsep+2) {
                    //case 53 - i,j,l
                    if(en4==v[j-1][i+1][a-1][b][aa][bb+1]
                    +3*data->eparam[10]
                    +penalty(i+1,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n,ct3,data)
                    +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle3(l-1,k,l,ct2,data)
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                       
						stack.push(i+1,j-1,a-1,b,aa,bb+1,v[j-1][i+1][a-1][b][aa][bb+1]);
                    }
                }
                
                //case 33 - i
                if(en4==v[j][i+1][a-1][b][aa][bb]
                +3*data->eparam[10]
                +penalty(i+1,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n,ct3,data)
                +edangle5(i+1,j,i,ct1,data)
                +data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    
					stack.push(i+1,j,a-1,b,aa,bb,v[j][i+1][a-1][b][aa][bb]);
                }
                
                //case 54 - i,j,l,n
                if(en4==v[j-1][i+1][a-1][b][aa][bb]
                +3*data->eparam[10]
                +penalty(i+1,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n-1,ct3,data)
                +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle3(l-1,k,l,ct2,data)+edangle3(n-1,m,n,ct3,data)
                +4*data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    
					stack.push(i+1,j-1,a-1,b,aa,bb,v[j-1][i+1][a-1][b][aa][bb]);
                }
            }
            
            if(a+1<2*maxsep+2) {
                
                if(b>=1) {
                    
                    if(aa>=1) {
                        
                        if(bb+1<2*maxsep+2) {
                            //case 13 - k,l
                            if(en4==v[j][i][a+1][b-1][aa-1][bb+1]
                            +3*data->eparam[10]
                            +penalty(i,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n,ct3,data)
                            +edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)
                            +2*data->eparam[6]+4*gap
                            &&!found) {
                                found = true;
                    
								stack.push(i,j,a+1,b-1,aa-1,bb+1,v[j][i][a+1][b-1][aa-1][bb+1]);
                            }
                        }
                        
                        //case 14 - k,l,n
                        if(en4==v[j][i][a+1][b-1][aa-1][bb]
                        +3*data->eparam[10]
                        +penalty(i,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n-1,ct3,data)
                        +edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle3(n-1,m,n,ct3,data)
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                    
							stack.push(i,j,a+1,b-1,aa-1,bb,v[j][i][a+1][b-1][aa-1][bb]);
                        }
                    }
                    
                    if(bb+1<2*maxsep+2) {
                        //case 15 - k,l,m
                        if(en4==v[j][i][a+1][b-1][aa][bb+1]
                        +3*data->eparam[10]
                        +penalty(i,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n,ct3,data)
                        +edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle5(m+1,n,m,ct3,data)
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                          
							stack.push(i,j,a+1,b-1,aa,bb+1,v[j][i][a+1][b-1][aa][bb+1]);
                        }
                    }
                    
                    //case 16 - k,l,m,n
                    if(en4==v[j][i][a+1][b-1][aa][bb]
                    +3*data->eparam[10]
                    +penalty(i,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
                    +edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        
						stack.push(i,j,a+1,b-1,aa,bb,v[j][i][a+1][b-1][aa][bb]);
                    }
                }
                
                if(b+1<2*maxsep+2) {
                    
                    if(aa>=1) {
                        
                        if(bb>=1) {
                            //case 26 - j,k,n
                            if(en4==v[j-1][i][a+1][b+1][aa-1][bb-1]
                            +3*data->eparam[10]
                            +penalty(i,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n-1,ct3,data)
                            +edangle3(j-1,i,j,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle3(n-1,m,n,ct3,data)
                            +3*data->eparam[6]+3*gap
                            &&!found) {
                                found = true;
                        
								stack.push(i,j-1,a+1,b+1,aa-1,bb-1,v[j-1][i][a+1][b+1][aa-1][bb-1]);
                            }
                        }
                        
                        //case 25 - j,k
                        if(en4==v[j-1][i][a+1][b+1][aa-1][bb]
                        +3*data->eparam[10]
                        +penalty(i,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n,ct3,data)
                        +edangle3(j-1,i,j,ct1,data)+edangle5(k+1,l,k,ct2,data)
                        +2*data->eparam[6]+4*gap
                        &&!found) {
                            found = true;
                        
							stack.push(i,j-1,a+1,b+1,aa-1,bb,v[j-1][i][a+1][b+1][aa-1][bb]);
                        }
                    }
                    
                    if(bb>=1) {
                        //case 28 - j,k,m,n
                        if(en4==v[j-1][i][a+1][b+1][aa][bb-1]
                        +3*data->eparam[10]
                        +penalty(i,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n-1,ct3,data)
                        +edangle3(j-1,i,j,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                        +4*data->eparam[6]+2*gap
                        &&!found) {
                            found = true;
                            
							stack.push(i,j-1,a+1,b+1,aa,bb-1,v[j-1][i][a+1][b+1][aa][bb-1]);
                        }
                    }
                    
                    //case 27 - j,k,m
                    if(en4==v[j-1][i][a+1][b+1][aa][bb]
                    +3*data->eparam[10]
                    +penalty(i,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n,ct3,data)
                    +edangle3(j-1,i,j,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle5(m+1,n,m,ct3,data)
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                        
                    }
                }
                
                if(aa>=1) {
                    
                    if(bb>=1) {
                        //case 10 - k,n
                        if(en4==v[j][i][a+1][b][aa-1][bb-1]
                        +3*data->eparam[10]
                        +penalty(i,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n-1,ct3,data)
                        +edangle5(k+1,l,k,ct2,data)+edangle3(n-1,m,n,ct3,data)
                        +2*data->eparam[6]+4*gap
                        &&!found) {
                            found = true;
                        
							stack.push(i,j,a+1,b,aa-1,bb-1,v[j][i][a+1][b][aa-1][bb-1]);
                        }
                    }
                    
                    if(bb+1<2*maxsep+2) {
                        //case 29 - j,k,l
                        if(en4==v[j-1][i][a+1][b][aa-1][bb+1]
                        +3*data->eparam[10]
                        +penalty(i,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n,ct3,data)
                        +edangle3(j-1,i,j,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                        
							stack.push(i,j-1,a+1,b,aa-1,bb+1,v[j-1][i][a+1][b][aa-1][bb+1]);
                        }
                    }
                    
                    //case 9 - k
                    if(en4==v[j][i][a+1][b][aa-1][bb]
                    +3*data->eparam[10]
                    +penalty(i,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n,ct3,data)
                    +edangle5(k+1,l,k,ct2,data)
                    +data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        
						stack.push(i,j,a+1,b,aa-1,bb,v[j][i][a+1][b][aa-1][bb]);
                    }
                    
                    //case 30 - j,k,l,n
                    if(en4==v[j-1][i][a+1][b][aa-1][bb]
                    +3*data->eparam[10]
                    +penalty(i,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n-1,ct3,data)
                    +edangle3(j-1,i,j,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle3(n-1,m,n,ct3,data)
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        
						stack.push(i,j-1,a+1,b,aa-1,bb,v[j-1][i][a+1][b][aa-1][bb]);
                    }
                }
                
                if(bb>=1) {
                    //case 12 - k,m,n
                    if(en4==v[j][i][a+1][b][aa][bb-1]
                    +3*data->eparam[10]
                    +penalty(i,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n-1,ct3,data)
                    +edangle5(k+1,l,k,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                        
						stack.push(i,j,a+1,b,aa,bb-1,v[j][i][a+1][b][aa][bb-1]);
                    }
                }
                
                if(bb+1<2*maxsep+2) {
                    //case 31 - j,k,l,m
                    if(en4==v[j-1][i][a+1][b][aa][bb+1]
                    +3*data->eparam[10]
                    +penalty(i,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n,ct3,data)
                    +edangle3(j-1,i,j,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle5(m+1,n,m,ct3,data)
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        
						stack.push(i,j-1,a+1,b,aa,bb+1,v[j-1][i][a+1][b][aa][bb+1]);
                    }
                }
                
                //case 11 - k,m
                if(en4==v[j][i][a+1][b][aa][bb]
                +3*data->eparam[10]
                +penalty(i,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n,ct3,data)
                +edangle5(k+1,l,k,ct2,data)+edangle5(m+1,n,m,ct3,data)
                +2*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    
					stack.push(i,j,a+1,b,aa,bb,v[j][i][a+1][b][aa][bb]);
                }
                
                //case 32 - j,k,l,m,n
                if(en4==v[j-1][i][a+1][b][aa][bb]
                +3*data->eparam[10]
                +penalty(i,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
                +edangle3(j-1,i,j,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                +5*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    
					stack.push(i,j-1,a+1,b,aa,bb,v[j-1][i][a+1][b][aa][bb]);
                }
            }
            
            if(b>=1) {
                
                if(aa>=1) {
                    
                    if(bb+1<2*maxsep+2) {
                        //case 45 - i,k,l
                        if(en4==v[j][i+1][a][b-1][aa-1][bb+1]
                        +3*data->eparam[10]
                        +penalty(i+1,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n,ct3,data)
                        +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                    
							stack.push(i+1,j,a,b-1,aa-1,bb+1,v[j][i+1][a][b-1][aa-1][bb+1]);
                        }
                    }
                    
                    //case 46 - i,k,l,n
                    if(en4==v[j][i+1][a][b-1][aa-1][bb]
                    +3*data->eparam[10]
                    +penalty(i+1,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n-1,ct3,data)
                    +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle3(n-1,m,n,ct3,data)
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                    
						stack.push(i+1,j,a,b-1,aa-1,bb,v[j][i+1][a][b-1][aa-1][bb]);
                    }
                }
                
                if(aa+1<2*maxsep+2) {
                    
                    if(bb+1<2*maxsep+2) {
                        //case 7 - l,m
                        if(en4==v[j][i][a][b-1][aa+1][bb+1]
                        +3*data->eparam[10]
                        +penalty(i,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n,ct3,data)
                        +edangle3(l-1,k,l,ct2,data)+edangle5(m+1,n,m,ct3,data)
                        +2*data->eparam[6]+4*gap
                        &&!found) {
                            found = true;
                        
							stack.push(i,j,a,b-1,aa+1,bb+1,v[j][i][a][b-1][aa+1][bb+1]);
                        }
                    }
                    
                    //case 8 - l,m,n
                    if(en4==v[j][i][a][b-1][aa+1][bb]
                    +3*data->eparam[10]
                    +penalty(i,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
                    +edangle3(l-1,k,l,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                        
						stack.push(i,j,a,b-1,aa+1,bb,v[j][i][a][b-1][aa+1][bb]);
                    }
                }
                
                if(bb+1<2*maxsep+2) {
                    //case 5 - l
                    if(en4==v[j][i][a][b-1][aa][bb+1]
                    +3*data->eparam[10]
                    +penalty(i,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n,ct3,data)
                    +edangle3(l-1,k,l,ct2,data)
                    +data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        
						stack.push(i,j,a,b-1,aa,bb+1,v[j][i][a][b-1][aa][bb+1]);
                    }
                    
                    //case 47 - i,k,l,m
                    if(en4==v[j][i+1][a][b-1][aa][bb+1]
                    +3*data->eparam[10]
                    +penalty(i+1,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n,ct3,data)
                    +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle5(m+1,n,m,ct3,data)
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        
						stack.push(i+1,j,a,b-1,aa,bb+1,v[j][i+1][a][b-1][aa][bb+1]);
                    }
                }
                
                //case 6 - l,n
                if(en4==v[j][i][a][b-1][aa][bb]
                +3*data->eparam[10]
                +penalty(i,j,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n-1,ct3,data)
                +edangle3(l-1,k,l,ct2,data)+edangle3(n-1,m,n,ct3,data)
                +2*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    
					stack.push(i,j,a,b-1,aa,bb,v[j][i][a][b-1][aa][bb]);
                }
                
                //case 48 - i,k,l,m,n
                if(en4==v[j][i+1][a][b-1][aa][bb]
                +3*data->eparam[10]
                +penalty(i+1,j,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
                +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                +5*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    
					stack.push(i+1,j,a,b-1,aa,bb,v[j][i+1][a][b-1][aa][bb]);
                }
            }
            
            if(b+1<2*maxsep+2) {
                
                if(aa>=1) {
                    
                    if(bb>=1) {
                        //case 58 - i,j,k,n
                        if(en4==v[j-1][i+1][a][b+1][aa-1][bb-1]
                        +3*data->eparam[10]
                        +penalty(i+1,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n-1,ct3,data)
                        +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle3(n-1,m,n,ct3,data)
                        +4*data->eparam[6]+2*gap
                        &&!found) {
                            found = true;
                    
							stack.push(i+1,j-1,a,b+1,aa-1,bb-1,v[j-1][i+1][a][b+1][aa-1][bb-1]);
                        }
                    }
                    
                    //case 57 - i,j,k
                    if(en4==v[j-1][i+1][a][b+1][aa-1][bb]
                    +3*data->eparam[10]
                    +penalty(i+1,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n,ct3,data)
                    +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle5(k+1,l,k,ct2,data)
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                    
						stack.push(i+1,j-1,a,b+1,aa-1,bb,v[j-1][i+1][a][b+1][aa-1][bb]);
                    }
                }
                
                if(aa+1<2*maxsep+2) {
                    
                    if(bb>=1) {
                        //case 20 - j,m,n
                        if(en4==v[j-1][i][a][b+1][aa+1][bb-1]
                        +3*data->eparam[10]
                        +penalty(i,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n-1,ct3,data)
                        +edangle3(j-1,i,j,ct1,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                        +3*data->eparam[6]+3*gap
                        &&!found) {
                            found = true;
                    
							stack.push(i,j-1,a,b+1,aa+1,bb-1,v[j-1][i][a][b+1][aa+1][bb-1]);
                        }
                    }
                    
                    //case 19 - j,m
                    if(en4==v[j-1][i][a][b+1][aa+1][bb]
                    +3*data->eparam[10]
                    +penalty(i,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n,ct3,data)
                    +edangle3(j-1,i,j,ct1,data)+edangle5(m+1,n,m,ct3,data)
                    +2*data->eparam[6]+4*gap
                    &&!found) {
                        found = true;
                    
						stack.push(i,j-1,a,b+1,aa+1,bb,v[j-1][i][a][b+1][aa+1][bb]);
                    }
                }
                
                if(bb>=1) {
                    //case 18 - j,n
                    if(en4==v[j-1][i][a][b+1][aa][bb-1]
                    +3*data->eparam[10]
                    +penalty(i,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n-1,ct3,data)
                    +edangle3(j-1,i,j,ct1,data)+edangle3(n-1,m,n,ct3,data)
                    +2*data->eparam[6]+gap
                    &&!found) {
                        found = true;
                       
						stack.push(i,j-1,a,b+1,aa,bb-1,v[j-1][i][a][b+1][aa][bb-1]);
                    }
                    
                    //case 60 - i,j,k,m,n
                    if(en4==v[j-1][i+1][a][b+1][aa][bb-1]
                    +3*data->eparam[10]
                    +penalty(i+1,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n-1,ct3,data)
                    +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                    +5*data->eparam[6]+gap
                    &&!found) {
                        found = true;
                       
						stack.push(i+1,j-1,a,b+1,aa,bb-1,v[j-1][i+1][a][b+1][aa][bb-1]);
                    }
                }
                
                //case 17 - j
                if(en4==v[j-1][i][a][b+1][aa][bb]
                +3*data->eparam[10]
                +penalty(i,j-1,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n,ct3,data)
                +edangle3(j-1,i,j,ct1,data)
                +data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    
					stack.push(i,j-1,a,b+1,aa,bb,v[j-1][i][a][b+1][aa][bb]);
                }
                
                //case 59 - i,j,k,m
                if(en4==v[j-1][i+1][a][b+1][aa][bb]
                +3*data->eparam[10]
                +penalty(i+1,j-1,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n,ct3,data)
                +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle5(m+1,n,m,ct3,data)
                +4*data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    
					stack.push(i+1,j-1,a,b+1,aa,bb,v[j-1][i+1][a][b+1][aa][bb]);
                }
            }
            
            if(aa>=1) {
                
                if(bb>=1) {
                    //case 42 - i,k,n
                    if(en4==v[j][i+1][a][b][aa-1][bb-1]
                    +3*data->eparam[10]
                    +penalty(i+1,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n-1,ct3,data)
                    +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle3(n-1,m,n,ct3,data)
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                    
						stack.push(i+1,j,a,b,aa-1,bb-1,v[j][i+1][a][b][aa-1][bb-1]);
                    }
                }
                
                if(bb+1<2*maxsep+2) {
                    //case 61 - i,j,k,l
                    if(en4==v[j-1][i+1][a][b][aa-1][bb+1]
                    +3*data->eparam[10]
                    +penalty(i+1,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n,ct3,data)
                    +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)
                    +4*data->eparam[6]+2*gap
                    &&!found) {
                        found = true;
                        
						stack.push(i+1,j-1,a,b,aa-1,bb+1,v[j-1][i+1][a][b][aa-1][bb+1]);
                    }
                }
                
                //case 41 - i,k
                if(en4==v[j][i+1][a][b][aa-1][bb]
                +3*data->eparam[10]
                +penalty(i+1,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m,n,ct3,data)
                +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l,k,ct2,data)
                +2*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    
					stack.push(i+1,j,a,b,aa-1,bb,v[j][i+1][a][b][aa-1][bb]);
                }
                
                //case 62 - i,j,k,l,n
                if(en4==v[j-1][i+1][a][b][aa-1][bb]
                +3*data->eparam[10]
                +penalty(i+1,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m,n-1,ct3,data)
                +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle3(n-1,m,n,ct3,data)
                +5*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    
					stack.push(i+1,j-1,a,b,aa-1,bb,v[j-1][i+1][a][b][aa-1][bb]);
                }
            }
            
            if(aa+1<2*maxsep+2) {
                
                if(bb>=1) {
                    //case 4 - m,n
                    if(en4==v[j][i][a][b][aa+1][bb-1]
                    +3*data->eparam[10]
                    +penalty(i,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n-1,ct3,data)
                    +edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                    +2*data->eparam[6]+4*gap
                    &&!found) {
                        found = true;
                    
						stack.push(i,j,a,b,aa+1,bb-1,v[j][i][a][b][aa+1][bb-1]);
                    }
                }
                
                if(bb+1<2*maxsep+2) {
                    //case 23 - j,l,m
                    if(en4==v[j-1][i][a][b][aa+1][bb+1]
                    +3*data->eparam[10]
                    +penalty(i,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n,ct3,data)
                    +edangle3(j-1,i,j,ct1,data)+edangle3(l-1,k,l,ct2,data)+edangle5(m+1,n,m,ct3,data)
                    +3*data->eparam[6]+3*gap
                    &&!found) {
                        found = true;
                    
						stack.push(i,j-1,a,b,aa+1,bb+1,v[j-1][i][a][b][aa+1][bb+1]);
                    }
                }
                
                //case 3 - m
                if(en4==v[j][i][a][b][aa+1][bb]
                +3*data->eparam[10]
                +penalty(i,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m+1,n,ct3,data)
                +edangle5(m+1,n,m,ct3,data)
                +data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    
					stack.push(i,j,a,b,aa+1,bb,v[j][i][a][b][aa+1][bb]);
                }
                
                //case 24 - j,l,m,n
                if(en4==v[j-1][i][a][b][aa+1][bb]
                +3*data->eparam[10]
                +penalty(i,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
                +edangle3(j-1,i,j,ct1,data)+edangle3(l-1,k,l,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                +4*data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    
					stack.push(i,j-1,a,b,aa+1,bb,v[j-1][i][a][b][aa+1][bb]);
                }
            }
            
            if(bb>=1) {
                //case 2 - n
                if(en4==v[j][i][a][b][aa][bb-1]
                +3*data->eparam[10]
                +penalty(i,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n-1,ct3,data)
                +edangle3(n-1,m,n,ct3,data)
                +data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    
					stack.push(i,j,a,b,aa,bb-1,v[j][i][a][b][aa][bb-1]);
                }
                
                //case 44 - i,k,m,n
                if(en4==v[j][i+1][a][b][aa][bb-1]
                +3*data->eparam[10]
                +penalty(i+1,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n-1,ct3,data)
                +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
                +4*data->eparam[6]+2*gap
                &&!found) {
                    found = true;
                    
					stack.push(i+1,j,a,b,aa,bb-1,v[j][i+1][a][b][aa][bb-1]);
                }
            }
            
            if(bb+1<2*maxsep+2) {
                //case 21 - j,l
                if(en4==v[j-1][i][a][b][aa][bb+1]
                +3*data->eparam[10]
                +penalty(i,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n,ct3,data)
                +edangle3(j-1,i,j,ct1,data)+edangle3(l-1,k,l,ct2,data)
                +2*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    
					stack.push(i,j-1,a,b,aa,bb+1,v[j-1][i][a][b][aa][bb+1]);
                }
                
                //case 63 - i,j,k,l,m
                if(en4==v[j-1][i+1][a][b][aa][bb+1]
                +3*data->eparam[10]
                +penalty(i+1,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n,ct3,data)
                +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle5(m+1,n,m,ct3,data)
                +5*data->eparam[6]+gap
                &&!found) {
                    found = true;
                    
					stack.push(i+1,j-1,a,b,aa,bb+1,v[j-1][i+1][a][b][aa][bb+1]);
                }
            }
            
            //case 1 -
            if(en4==v[j][i][a][b][aa][bb]
            +3*data->eparam[10]
            +penalty(i,j,ct1,data)+penalty(k,l,ct2,data)+penalty(m,n,ct3,data)
            &&!found) {
                found = true;
                
				stack.push(i,j,a,b,aa,bb,v[j][i][a][b][aa][bb]);
            }
            
            //case 22 - j,l,n
            if(en4==v[j-1][i][a][b][aa][bb]
            +3*data->eparam[10]
            +penalty(i,j-1,ct1,data)+penalty(k,l-1,ct2,data)+penalty(m,n-1,ct3,data)
            +edangle3(j-1,i,j,ct1,data)+edangle3(l-1,k,l,ct2,data)+edangle3(n-1,m,n,ct3,data)
            +3*data->eparam[6]
            &&!found) {
                found = true;
                
				stack.push(i,j-1,a,b,aa,bb,v[j-1][i][a][b][aa][bb]);
            }
            
            //case 43 - i,k,m
            if(en4==v[j][i+1][a][b][aa][bb]
            +3*data->eparam[10]
            +penalty(i+1,j,ct1,data)+penalty(k+1,l,ct2,data)+penalty(m+1,n,ct3,data)
            +edangle5(i+1,j,i,ct1,data)+edangle5(k+1,l,k,ct2,data)+edangle5(m+1,n,m,ct3,data)
            +3*data->eparam[6]
            &&!found) {
                found = true;
                
				stack.push(i+1,j,a,b,aa,bb,v[j][i+1][a][b][aa][bb]);
            }

		
            //case 64 - i,j,k,l,m,n
            if(en4==v[j-1][i+1][a][b][aa][bb]
            +3*data->eparam[10]
            +penalty(i+1,j-1,ct1,data)+penalty(k+1,l-1,ct2,data)+penalty(m+1,n-1,ct3,data)
            +edangle5(i+1,j-1,i,ct1,data)+edangle3(j-1,i+1,j,ct1,data)+edangle5(k+1,l-1,k,ct2,data)+edangle3(l-1,k+1,l,ct2,data)+edangle5(m+1,n-1,m,ct3,data)+edangle3(n-1,m+1,n,ct3,data)
            +6*data->eparam[6]
            &&!found) {
                found = true;
        
				stack.push(i+1,j-1,a,b,aa,bb,v[j-1][i+1][a][b][aa][bb]);
            }
            
            
            
            /***************************** traceback W3 *********************************/
            //calculate the free energy of 2 fragments merged:
            for (c=i+minloop;c<j-minloop&&!found;c++) {
                for (d=max(k+minloop,c-maxsep);d<l-minloop&&d<=c+maxsep&&!found;d++) {
                    e = d-c+maxsep;
                    for (f=max(m+minloop,d-maxsep);f<n-minloop&&f<=d+maxsep&&!found;f++) {
                        g=f-d+maxsep;
                        
                        if(en4 ==w[c][i][a][e][aa][g]+w[j][c+1][e][b][g][bb]) {
                            found = true;
                            stack.push(i,c,a,e,aa,g,w[c][i][a][e][aa][g]);
		
                            stack.push(c+1,j,e,b,g,bb,w[j][c+1][e][b][g][bb]);
		
                            
                        }
                        
                    }
                }
            }
            
            if (!found) {
                cout << "Traceback error for W!\n"<<flush;
                return;
            }
            
        }
        
        
    }	//end of while stack not empty
    
    cout << "end of find interior \n"<<flush;
    
    
    /***************************** end of traceback ***********************/
    
    for (i=0;i<=ct1->numofbases;i++) {
        
        for (j=0;j<i;j++) {
            
            for (k=0;k<2*maxseparation+2;k++) {
                
                for(l=0;l<2*maxseparation+2;l++) {
                    
                    for(m=0;m<2*maxseparation+2;m++) {
                        
                        delete[] v[i][j][k][l][m];
                    }
                    
                    delete[] v[i][j][k][l];
                }
                
                delete[] v[i][j][k];
            }
            delete[] v[i][j];
            
        }
        delete[] v[i];
    }
    delete[] v;
    
    for (i=0;i<=ct1->numofbases;i++) {
        
        for (j=0;j<i;j++) {
            
            for (k=0;k<2*maxseparation+2;k++) {
                
                for(l=0;l<2*maxseparation+2;l++) {
                    
                    for(m=0;m<2*maxseparation+2;m++) {
                        
                        delete[] w[i][j][k][l][m];
                        
                    }
                    delete[] w[i][j][k][l];
                    
                }
                
                delete[] w[i][j][k];
            }
            delete[] w[i][j];
            
        }
        delete[] w[i];
    }
    delete[] w;
    
    
    for (i=0;i<=ct1->numofbases;i++) {
        for (k=0; k<2*maxseparation+2;k++) {
            delete[] mine[i][k];
        }
        delete[] mine[i];
        
    }
    delete[] mine;

    delete[] pair[0];
    delete[] pair[1];
    delete[] pair[2];
    delete[] pair;
    
    cout << "end of find xdynalign function \n"<<flush;
}


void alignout(short* align1,short* align2,char *aout, structure *ct1, structure *ct2, structure *ct3 ,int totalscore) {
  ofstream out;

  char *line1,*line2,*line3,*line4;
  short i,k,m;
  short j,l,n,t,len_max;

  line1 = new char [ct1->numofbases+ct2->numofbases+ct3->numofbases+100];
  line2 = new char [ct1->numofbases+ct2->numofbases+ct3->numofbases+100];
  line3 = new char [ct1->numofbases+ct2->numofbases+ct3->numofbases+100];
  line4 = new char [ct1->numofbases+ct2->numofbases+ct3->numofbases+100];

  out.open(aout);

  strcpy(line1,"\n");
  strcpy(line2,"\n");
  strcpy(line3,"\n");
  strcpy(line4,"\n");

  i=1; k=1; m=1; 
  for(j=1;j<=ct1->numofbases;j++)
    {
      if(align1[j]>0)
	{
	  l=align1[j];
	  n=align2[l];
	  len_max=max(j-i+1,max(l-k+1,n-m+1));


	  // Fill in Seq1

	  for(t=0;t<len_max-j+i-1;t++) strcat(line1,"-");
	  for(t=i;t<=j;t++)
	    {
	      line1[strlen(line1)+1]='\0';
	      line1[strlen(line1)]=ct1->nucs[t];
	    }
	  i=j+1;

	  // Fill in Seq2
	  for(t=0;t<len_max-l+k-1;t++) strcat(line2,"-");
	  for(t=k;t<=l;t++)
	    {
	      line2[strlen(line2)+1]='\0';
	      line2[strlen(line2)]=ct2->nucs[t];
	    }
	  k=l+1;

	  // Fill in Seq3
	  for(t=0;t<len_max-n+m-1;t++) strcat(line3,"-");
	  for(t=m;t<=n;t++)
	    {
	      line3[strlen(line3)+1]='\0';
	      line3[strlen(line3)]=ct3->nucs[t];
	    }
	  m=n+1;

	  // Mark the alignment
	  for(t=0;t<len_max-1;t++) strcat(line4," ");
	  strcat(line4,"^");
	}
    }

  // Print all left overs
  j=ct1->numofbases;
  l=ct2->numofbases;
  n=ct3->numofbases;
  len_max=max(j-i+1,max(l-k+1,n-m+1));

  // Fill in Seq1
  for(t=i;t<=j;t++)
    {
      line1[strlen(line1)+1]='\0';
      line1[strlen(line1)]=ct1->nucs[t];
    }
  for(t=0;t<len_max-j+i-1;t++) strcat(line1,"-");

  // Fill in Seq2
  for(t=k;t<=l;t++)
    {
      line2[strlen(line2)+1]='\0';
      line2[strlen(line2)]=ct2->nucs[t];
    }
  for(t=0;t<len_max-l+k-1;t++) strcat(line2,"-");

  // Fill in Seq3
  for(t=m;t<=n;t++)
    {
      line3[strlen(line3)+1]='\0';
      line3[strlen(line3)]=ct3->nucs[t];
    }
  for(t=0;t<len_max-n+m-1;t++) strcat(line3,"-");

  // Mark the alignment
  for(t=0;t<len_max;t++) strcat(line4," ");

  if (totalscore!=500) out<<"Score= "<<totalscore<<"\n";	
  out <<line1<<"\n"<<line2<<"\n"<<line3<<"\n"<<line4<<"\n";
 
 
  out.close();
  delete[] line1;
  delete[] line2;
  delete[] line3;
  delete[] line4;
}






//calculate the energy of a bulge/internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
short int einternal(int i,int j,int ip,int jp,structure *ct, datatable *data) {
  int energy,size,size1,size2,loginc, lopsid;//tlink,count,key,e[4]
  /* size,size1,size2 = size of a loop
     energy = energy calculated
     loginc = the value of a log used in large hairpin loops
  */
    
    
  size1 = ip-i-1;
  size2 = j - jp - 1;
    
    
    
  //a typical internal or bulge loop:
  //size1 = ip-i-1;
  //size2 = j - jp - 1;
  if (size1==0||size2==0) {//bulge loop
        
        
    size = size1+size2;
    if (size==1) {
      energy = data->stack[ct->numseq[i]][ct->numseq[j]]
	[ct->numseq[ip]][ct->numseq[jp]]
	+ data->bulge[size] + data->eparam[2];
    }
    else if (size>30) {
            
      loginc = int((data->prelog)*log(double ((size)/30.0)));
      energy = data->bulge[30] + loginc + data->eparam[2];
      energy = energy + penalty(i,j,ct,data) + penalty(jp,ip,ct,data);
            
    }
    else {
      energy = data->bulge[size] + data->eparam[2];
      energy = energy + penalty(i,j,ct,data) + penalty(jp,ip,ct,data);
    }
  }
  else {//internal loop
    size = size1 + size2;
    lopsid = abs(size1-size2);
        
    if (size>30) {
            
      loginc = int((data->prelog)*log((double ((size))/30.0)));
      if ((size1==1||size2==1)&&data->gail) {
	energy = data->tstki[ct->numseq[i]][ct->numseq[j]][1][1] +
	  data->tstki[ct->numseq[jp]][ct->numseq[ip]][1][1] +
	  data->inter[30] + loginc + data->eparam[3] +
	  min(data->maxpen,(lopsid*
			    data->poppen[min(2,min(size1,size2))]));
                
      }
            
      else {
	energy = data->tstki[ct->numseq[i]][ct->numseq[j]]
	  [ct->numseq[i+1]][ct->numseq[j-1]] +
	  data->tstki[ct->numseq[jp]][ct->numseq[ip]]
	  [ct->numseq[jp+1]][ct->numseq[ip-1]] +
	  data->inter[30] + loginc + data->eparam[3] +
	  min(data->maxpen,(lopsid*
			    data->poppen[min(2,min(size1,size2))]));
      }
    }
    else if ((size1==2)&&(size2==2))//2x2 internal loop
      energy = data->iloop22[ct->numseq[i]][ct->numseq[ip]]
	[ct->numseq[j]][ct->numseq[jp]]
	[ct->numseq[i+1]][ct->numseq[i+2]]
	[ct->numseq[j-1]][ct->numseq[j-2]];
        
        
    else if ((size1==1)&&(size2==2)) {//2x1 internal loop
      energy = data->iloop21[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]]
	[ct->numseq[j-1]][ct->numseq[jp+1]][ct->numseq[ip]][ct->numseq[jp]];
            
            
    }
    else if ((size1==2)&&(size2==1)) {//1x2 internal loop
      energy = data->iloop21[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]]
	[ct->numseq[ip-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]];

    }
        
    else if (size==2) //a single mismatch
      energy = data->iloop11[ct->numseq[i]][ct->numseq[i+1]][ct->numseq[ip]]
	[ct->numseq[j]][ct->numseq[j-1]][ct->numseq[jp]];

    else if ((size1==1||size2==1)&&data->gail) { //this loop is lopsided
      //note, we treat this case as if we had a loop composed of all As
      //if and only if the gail rule is set to 1 in miscloop.dat
      energy = data->tstki[ct->numseq[i]][ct->numseq[j]][1][1] +
	data->tstki[ct->numseq[jp]][ct->numseq[ip]][1][1] +
	data->inter[size] + data->eparam[3] +
	min(data->maxpen,(lopsid*
			  data->poppen[min(2,min(size1,size2))]));
    }
        
    else {
                    
      energy = data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
	data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
	data->inter[size] + data->eparam[3] +
	min(data->maxpen,(lopsid*data->poppen[min(2,min(size1,size2))]));
    }
  }
    
  return energy;
}



inline short int ehairpin(short i, short j, structure *ct, datatable *data) {
    
  int energy,size,loginc,tlink,count,key,k;
  size = j-i-1;
    
    
    
  if (size>30) {
        
    loginc = int((data->prelog)*log((double ((size))/30.0)));
        
    energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
      [ct->numseq[i+1]][ct->numseq[j-1]]
      + data->hairpin[30]+loginc+data->eparam[4];
  }
  else if (size<3) {
    energy = data->hairpin[size] + data->eparam[4];
    if (ct->numseq[i]==4||ct->numseq[j]==4) energy = energy+6;
  }
  else if (size==4) {
    tlink = 0;
    key = (ct->numseq[j])*3125 + (ct->numseq[i+4])*625 +
      (ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
    for (count=1;count<=data->numoftloops&&tlink==0;count++) {
      if (key==data->tloop[count][0]) tlink = data->tloop[count][1];
    }
    energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
      [ct->numseq[i+1]][ct->numseq[j-1]]
      + data->hairpin[size] + data->eparam[4] + tlink;
  }
  else if (size==3) {
    tlink = 0;
    key = (ct->numseq[j])*625 +
      (ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
    for (count=1;count<=data->numoftriloops&&tlink==0;count++) {
      if (key==data->triloop[count][0]) tlink = data->triloop[count][1];
    }
        
    energy =	data->hairpin[size] + data->eparam[4] + tlink
      +penalty(i,j,ct,data);
  }
    
  else {
    energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
      [ct->numseq[i+1]][ct->numseq[j-1]]
      + data->hairpin[size] + data->eparam[4];
  }
    
  //check for GU closeure preceded by GG
  if (ct->numseq[i]==3&&ct->numseq[j]==4) {
    if ((i>2&&i<ct->numofbases)||(i>ct->numofbases+2))
      if (ct->numseq[i-1]==3&&ct->numseq[i-2]==3) {
                
	energy = energy + data->gubonus;
                
                
                
      }
  }
    
  //check for a poly-c loop
  tlink = 1;
  for (k=1;(k<=size)&&(tlink==1);k++) {
    if (ct->numseq[i+k] != 2) tlink = 0;
  }
  if (tlink==1) {  //this is a poly c loop so penalize
    if (size==3) energy = energy + data->c3;
    else energy = energy + data->cint + size*data->cslope;
  }
    
  return energy;
}

inline short int edangle5(int i,int j,int ip,structure* ct,datatable* data) {
  return data->dangle[ct->numseq[j]][ct->numseq[i]][ct->numseq[ip]][2];
}

inline short int edangle3(int i,int j,int ip,structure* ct,datatable* data) {
  return data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][1];
}

inline short int ebp(int i,int j,int ip,int jp,structure *ct, datatable *data) {
  return data->stack[(ct->numseq[i])][(ct->numseq[j])][(ct->numseq[ip])][(ct->numseq[jp])];
}
