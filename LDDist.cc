#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "LDDexcept.h"
#define TINY 1.0e-20	// A small number.
#define BIGINT 30000	// A biggish number.


const char VERSTR[] =  "LogDet-distances PERL module in C++, version 1.3.2 2005-06-28";
const float LIMITFRACTION=0.5; // The fraction to subsitute when the distance matrix has a zero on the diagonal
const int TaxNumCutOff = 50;
const int RepNums = TaxNumCutOff*(TaxNumCutOff-1)*(TaxNumCutOff-2)*(TaxNumCutOff-3)/8;

const char* GetVersion () {
	return VERSTR;
} 

class TCM {
	private:
		int FXYSIZE;				//Size of comparison matrix (4 or 20)
		int AMBIG;
		int GAP;
		int GXYSIZE;
		
		int TaxNum;					//Number of taxa (rows)
		int SiteNum;				//Number of sites (columns)
		
		int *theTCM;				//The actual taxon-character matrix
		float *PairWise;			//Matrix of pair-wise distances
		
		int ClassNum;				//Number of rate classes
		int *RateClasses;
		int *Constants;				//Mask for constant sites
		float InvariantFraction;	//Fraction of constant sites treated as invariant and excluded
		
		float *Fxy;
		int *Gxy;

		int FullFour (int i, int j, int k, int l, int site);
		float QEval (int i, int j, int k, int l);
		int AAtoNum(char residue);
		int DNAtoNum(char nucleotide);
		

		void NormFxy();
		float LogDet();

	public:

		TCM();
		// Do not use default constructor
		
		//TCM(char** newmatrix);
		
		~TCM(void);
		// Do not use default destructor
		
		//TCM operator = (const TCM &old_TCM)
		// Use default assignment operator

		//TCM (const TCM &oldTCM)
		// Use default copy constructor
		
		float SidowCRCInvariants (void); 
		void SWgrouping(int numberofclasses);
		void SetInvariants (float theFraction);
		float SteelCRCInvariants (void);
		int LogDetDistances(int bootstrap);

		void setTCM(char **newmatrix, int DNA);
		float getDistance (int thePair);
		
		int GetRateClass (int Site);
		int SetRateClass (int Site, int theClass);
		int NewRateClasses (int NumOfClasses);

};


 TCM::TCM () {

 
 }
 
 
 TCM::~TCM (void) {

	delete [] Fxy;
	delete [] Gxy;
	delete [] theTCM;
	delete [] RateClasses;
	delete [] Constants;
	delete [] PairWise;
}

float TCM::getDistance (int thePair) {
	if (thePair>((TaxNum*(TaxNum-1)/2)-1)) {
		throw RangeError();
	}
	return PairWise[thePair];
}

void TCM::setTCM(char **newmatrix, int DNA) {


// Initialize some constants
	if (DNA==0) {
		FXYSIZE=20; //I.e., AA
	} else {
		FXYSIZE=4;  //I.e., DNA 
	};
	AMBIG=FXYSIZE;
	GAP=AMBIG+1;
	GXYSIZE=GAP+1;
//Init done



	Fxy = (float*)(new float[FXYSIZE*FXYSIZE]);
	Gxy = (int*)(new int[GXYSIZE*GXYSIZE]);
	
	srand(time(NULL));
	int sites=0, taxa=0;
	while (newmatrix[0][sites]) sites++;
	while (newmatrix[taxa]) taxa++;


	theTCM= (int*) (new int[taxa*sites]);

	RateClasses = new int[sites];
	ClassNum=1;
	Constants= new int[sites];
	InvariantFraction = 0;
	
	for (int i=0;i<sites;i++) {
		RateClasses[i]=0;
		Constants[i]=0;
	}
	int i=0;

	while (newmatrix[i]) {
		int j=0;
		while (newmatrix[i][j]) {

			theTCM[i*sites+j]=(DNA) ? TCM::DNAtoNum(newmatrix[i][j]) : TCM::AAtoNum(newmatrix[i][j]);//Convert DNA or AA to numerals
			j++;
		}

		i++;
	}
	TaxNum=taxa;
	SiteNum=sites;
	PairWise = new float[TaxNum*(TaxNum-1)/2];

}


inline int TCM::FullFour (int i, int j, int k, int l, int site) {

if (theTCM[SiteNum*i+site]==GAP || 
	theTCM[SiteNum*i+site]==AMBIG || 
	theTCM[SiteNum*j+site]==GAP || 
	theTCM[SiteNum*j+site]==AMBIG || 
	theTCM[SiteNum*k+site]==GAP || 
	theTCM[SiteNum*k+site]==AMBIG || 
	theTCM[SiteNum*l+site]==GAP || 
	theTCM[SiteNum*l+site]==AMBIG) {
		return 0;
	} else {
		return 1;
	}
}


inline float TCM::QEval (int i, int j, int k, int l) {

	int lmij,lmkl,lmik,lmjl,lmil,lmjk,lmijkl,lmikjl,lmiljk;
	int total;
	int unvaried;
	

	int mij=0;
	int mik=0;
	int mil=0;
	int mjk=0;
	int mjl=0;
	int mkl=0;
	
	int mijkl=0;
	int mikjl=0;
	int miljk=0;
	unvaried=0;
	total=0;
	
	for (int site=0;site<SiteNum;site++) {
		if (FullFour(i,j,k,l,site)) {
			total++;
			lmij = int (!(theTCM[SiteNum*i+site]==theTCM[SiteNum*j+site]));
			lmik = int (!(theTCM[SiteNum*i+site]==theTCM[SiteNum*k+site]));
			lmil = int (!(theTCM[SiteNum*i+site]==theTCM[SiteNum*l+site]));
			lmjk = int (!(theTCM[SiteNum*j+site]==theTCM[SiteNum*k+site]));
			lmjl = int (!(theTCM[SiteNum*j+site]==theTCM[SiteNum*l+site]));
			lmkl = int (!(theTCM[SiteNum*k+site]==theTCM[SiteNum*l+site]));
			
			mij+=lmij;
			mik+=lmik;
			mil+=lmil;
			mjk+=lmjk;
			mjl+=lmjl;
			mkl+=lmkl;
			
			mijkl+=int(lmij && lmkl);
			mikjl+=int(lmik && lmjl);
			miljk+=int(lmil && lmjk);
			unvaried+= (!(lmij||lmik||lmil||lmjk||lmjl||lmkl));
		}
	}

	int Term1 = mij*mkl/mijkl;
	int Term2 = mik*mjl/mikjl;
	int Term3 = mil*mjk/miljk;


	
	int MaxTerm = (Term1>Term2 ? (Term1>Term3 ? Term1 : (Term3>Term2 ? Term3 : Term2)) : (Term2>Term3 ? Term2 : Term3));

	float TempSum = (MaxTerm>total ? 1.0 : ((MaxTerm>(total-unvaried)) ? (float) ((float) MaxTerm/total) : (float) (1.0-(float)unvaried/total)));
	return ((1.0-TempSum)*total/unvaried);

}

float TCM::SteelCRCInvariants (void) {
// Estimates the fraction of constant sites that are invariant by Capture-Recapture according to Steel et al.
// Steel, M., Huson, D. and Lockhart, P.J. (2000) Invariable sites models and their use in phylogeny reconstruction. Syst. Biol., 49, 225-232.
// Returns the fraction of the constant sites that are also estimated to be invariant 



const int pick = 4;
const int scale = 2;
int RepsToDo = RepNums/(scale*(TaxNum/pick));

int replicates = 0;
float pSt = 0;
int theTaxa[pick];
int *theArr;
theArr = new int[TaxNum];

if (TaxNum>TaxNumCutOff) {
	for (int Rep=0;Rep<RepsToDo;Rep++) {
		// Initialize the list
		for (int i=0;i<TaxNum;i++) {
				theArr[i] = i;			
		}
		int CurrentSize=TaxNum-1;
	
		while (CurrentSize>pick) {
			// Choose four
			for (int tx=0;tx<pick;tx++) {
					int choice = rand()%(CurrentSize);
					theTaxa[tx] = theArr[choice];
					theArr[choice]=theArr[CurrentSize];
					CurrentSize--;
			}
			pSt+=TCM::QEval(theTaxa[0],theTaxa[1],theTaxa[2],theTaxa[3]);
			replicates++;
			
		}
	}
} else {
	for (int l=3;l<TaxNum;l++) {
		for (int k=2;k<l;k++) {
			for (int j=1;j<k;j++) {
				for (int i=0;i<j;i++) {
					pSt+=TCM::QEval(i,j,k,l);
					replicates++;
				}
			}
		}
	}
}
//return vSum/replicates;
return pSt/replicates;
}

float TCM::SidowCRCInvariants (void) {
// Estimates the fraction of constant sites that are invariant by Capture-Recapture according to Sidow et al.
// Sidow, A., Nguyen, T. and Speed, T.P. (1992) Estimating the fraction of invariable codons with a capture-recapture method. J. Mol. Evol., 35, 253-260.
// Returns the fraction of the constant sites that are also estimated to be invariant 

// Standard genetic code
const char FirstNuc[] = "GMAGTGCGCAYAATCTATTG";
const char SecNuc  [] = "CGAAGAAGATTATTCCCGAT";

int Nnn=0; // Positions with a change in both first and second base
int Nnp=0; // Positions with a change in first base
int Npn=0; // Positions with a change in second base
int JustThirdChanged = 0; //Positions with an AA substitution due to a change in third base

for (int site=0;site<SiteNum;site++) {
	int FirstChanged = 0;
	int SecondChanged = 0;
	int NotSame = 0;
	int CompAA = BIGINT;
	
	if (!(theTCM[site]==21 || theTCM[site]==20)) {
		CompAA = theTCM[site];
	}
	 
	for (int Seq1=1;Seq1<TaxNum; Seq1++) {
		if (theTCM[SiteNum*Seq1+site]==21 || theTCM[SiteNum*Seq1+site]==20) {
			// If the first sequence contains a gap or a wildcard, do nothing
		}  else {
			if (CompAA > GXYSIZE) {
				CompAA = theTCM[SiteNum*Seq1+site];
			} else {
				if (CompAA != theTCM[SiteNum*Seq1+site]) {
					NotSame = 1;
				}
			}
			for (int Seq2=0;Seq2<Seq1; Seq2++) {			
				if (theTCM[SiteNum*Seq2+site]==21 || theTCM[SiteNum*Seq2+site]==20) {
					// If the second sequence contains a gap or a wildcard, do nothing
				} else {			
					if (FirstNuc[theTCM[SiteNum*Seq1+site]]!=FirstNuc[theTCM[SiteNum*Seq2+site]]) {
						FirstChanged=1;
					}
					if (SecNuc[theTCM[SiteNum*Seq1+site]]!=SecNuc[theTCM[SiteNum*Seq2+site]]) {
						SecondChanged=1;
					}
				}
			}
		}
	}

	if (FirstChanged) {Nnp++;};
	if (SecondChanged) {Npn++;};
	if (FirstChanged && SecondChanged) {Nnn++;};
	if (!(FirstChanged || SecondChanged) && NotSame) {
		JustThirdChanged++;
	};

}

float Npp = float(Nnp*Npn)/float (Nnn);
int ConstSites = SiteNum - JustThirdChanged - (Nnp+Npn-Nnn); // Constant sites are total sites - sites varying in 3rd positions - sites varying in 1st, 2nd or both
float InvSites = float (SiteNum) - float (JustThirdChanged) - Npp; // Invariant sites are total sites minus sites not examined minus sites estimated to be variable

return InvSites/ConstSites;

//return (1- float(Nnp*Npn)/(float (Nnn)*SiteNum));
}

void TCM::SetInvariants (float theFraction)
{
//int ConstantSites=0;
	if (theFraction>1) {
		throw RangeError();
	} else {
		InvariantFraction = theFraction;

		for (int site=0;site<SiteNum;site++) {
			Constants[site]=1;
			for (int taxon=0;taxon<TaxNum;taxon++) {
				if (theTCM[site]!=theTCM[SiteNum*taxon+site] || theTCM[SiteNum*taxon+site]==GAP || theTCM[SiteNum*taxon+site]==AMBIG) {
					Constants[site]=0;
					break;
				}
			}
		}
	}
}

void TCM::SWgrouping(int numberofclasses)
{
//
//	Calculates the Shannon-Wiener information index
//	for an array as -sum((pi)(log2 (pi))
//	The actual values used are -sum((pi)(log (pi))
//	

float normaliser1 = log(FXYSIZE);
//const float normaliser = 2;

float *SWs;
SWs = new float[SiteNum];

float max = 0;
float *classes;

	for (int site=0;site<SiteNum;site++) {
		float SW=0;
		int *tmpSW = new int[GXYSIZE];
		for (int AA=0;AA<GXYSIZE;AA++) {
			tmpSW[AA]=0; // init tmpSW
		}
		for (int taxon=0;taxon<TaxNum;taxon++) {
			tmpSW[theTCM[taxon*SiteNum+site]]++;
		}
		for (int AA=0;AA<FXYSIZE;AA++) {			
			SW-=(tmpSW[AA])?(tmpSW[AA]/(float)(TaxNum-tmpSW[AMBIG]-tmpSW[GAP]))*log(tmpSW[AA]/(float)(TaxNum-tmpSW[AMBIG]-tmpSW[GAP])):0;
		}
		SWs[site]=SW;

if (max < SW) {
	max=SW;
}

} //

	float normaliser2 = max+0.001;
	for (int site=0;site<SiteNum;site++) {
		RateClasses[site]=int(SWs[site]*(float)numberofclasses/normaliser2);

	}
	
	ClassNum=numberofclasses;


}

int TCM::GetRateClass (int Site)
{
	if (Site<0 || Site>=SiteNum) {
		throw ArrayOutOfBound();
	} else {
		return RateClasses[Site];
	}
}

int TCM::SetRateClass (int Site, int theClass)
{
	if (Site<0 || Site>=SiteNum || theClass>=ClassNum || theClass<0) {
		throw ArrayOutOfBound();
	} else {
		int OldClass=RateClasses[Site];
		RateClasses[Site] = theClass;
		return OldClass;
	}
}

int TCM::NewRateClasses (int NumOfClasses)
{
	ClassNum=NumOfClasses;
	for (int i=0;i<SiteNum;i++) {
		RateClasses[i]=0;
	}

}
double Det(float *a, int n)
//Given a matrix a[0..n-1][0..n-1], this routine replaces it by the LU decomposition of a rowwise
//permutation of itself. It returns the determinant of the matrix;
//The routine originates from the TNT/JAMA template library, http://math.nist.gov/tnt/, and it is
//implemented here directly in the LDDist module due to some compatibility issues
//

{
// Use a "left-looking", dot-product, Crout/Doolittle algorithm.

int pivsign = 1;

int *piv = new int[n];
for (int i = 0; i < n; i++) {
	piv[i] = i;
}

double *LUcolj = new double[n];


// Outer loop.

for (int j = 0; j < n; j++) {

// Make a copy of the j-th column to localize references.

	for (int i = 0; i < n; i++) {
		LUcolj[i] = a[i*n+j];
	}

	// Apply previous transformations.

	for (int i = 0; i < n; i++) {

		// Most of the time is spent in the following dot product.

		int kmax = i>j ? j : i;
		double s = 0.0;
		for (int k = 0; k < kmax; k++) {
			s += a[i*n+k]*LUcolj[k];
		}

		a[i*n+j] = LUcolj[i] -= s;
	}

	// Find pivot and exchange if necessary.

	int p = j;
	for (int i = j+1; i < n; i++) {
		if (fabs(LUcolj[i]) > fabs(LUcolj[p])) {
			p = i;
		}
	}
	if (p != j) {
		for (int k = 0; k < n; k++) {
			double t = a[p*n+k]; 
			a[p*n+k] = a[j*n+k]; 
			a[j*n+k] = t;
		}
		int k = piv[p]; 
		piv[p] = piv[j]; 
		piv[j] = k;
		pivsign = -pivsign;
	}

	// Compute multipliers.

	if ((j < n) && (a[j*n+j] != 0.0)) {
		for (int i = j+1; i < n; i++) {
			a[i*n+j] /= a[j*n+j];
		}
	}
}
double d = double(pivsign);
for (int row = 0; row < n; row++) {
	if (a[row*n+row] == 0.0) throw SingularMatrix();
		d *= a[row*n+row];
	}
return d;
}



void TCM::NormFxy ()
{
	float mxsum=0;
	for (int Row=0;Row<FXYSIZE;Row++) {
		for (int Col=0;Col<FXYSIZE;Col++) {
			mxsum+= Fxy[FXYSIZE*Row+Col];
		}
	}
	
	for (int Row=0;Row<FXYSIZE;Row++) {
		for (int Col=0;Col<FXYSIZE;Col++) {
			Fxy[FXYSIZE*Row+Col]/=mxsum;
		}
	}

}

float TCM::LogDet () 
{
// Debug output
//cout<<"Before MxNorm\nFxy=[[";
//for (int S1=0;S1<FXYSIZE;S1++) {
//	for (int S2=0;S2<FXYSIZE;S2++) {
//		cout << Fxy[FXYSIZE*S1+S2]<< ",";
//	}
//	cout << "];[";
//}
//cout<<"\n";
// end debug output

double *PiX = new double[FXYSIZE];
double *PiY = new double[FXYSIZE];
double normaliser = 1;


NormFxy(); //Normalise Fxy to 1

for (int pos =0;pos<FXYSIZE;pos++) {
	PiX[pos]=0;
	PiY[pos]=0;
}

for (int Row=0;Row<FXYSIZE;Row++) {
	for (int Col=0;Col<FXYSIZE;Col++) {
		PiX[Row] += Fxy[FXYSIZE*Row+Col];
		PiY[Row] += Fxy[FXYSIZE*Col+Row];
	}
}

for (int Row=0;Row<FXYSIZE;Row++) {
	normaliser *= PiX[Row]*PiY[Row];
}

normaliser = log(normaliser)/2;

//float rCorrect = FXYSIZE; // Original correction
//New correction according to Tamura & Kumar, 2002

float SquareSum = 0;
for (int i = 0;i < FXYSIZE; i++) {
	SquareSum += (PiX[i]+PiY[i])*(PiX[i]+PiY[i]);
}
float rCorrect = (FXYSIZE-1)/(1-SquareSum/4); // 4 is 2 squared, ie 2*100%

double theDet = (double) fabs(Det((float*)Fxy,FXYSIZE));
float LD = (float) (normaliser-log(theDet))/rCorrect;

delete [] PiY;
delete [] PiX;


return LD;
}


int TCM::LogDetDistances (int bootstrap) {


int *bootweight;
bootweight= new int[SiteNum];

//
// If a bootstrap replicate, generate a random bootstrap weightvector
// Else, fill with weight=1 for all sites
//

if (bootstrap) {
	for (int site=0;site<SiteNum;site++){
		bootweight[site]=0;
	}
	for (int site=0;site<SiteNum;site++) {
		int thesite = rand()%(SiteNum-1);
		bootweight[thesite]++;
	}	
	int slask = 0;
} else {
	for (int site=0;site<SiteNum;site++){
		bootweight[site]=1;
	}

}

//
// The actual distances
//

int pair=0;	// Init pairwise array index.

for (int taxon2=1;taxon2<TaxNum;taxon2++) {
	for (int taxon1=0;taxon1<taxon2;taxon1++) {
		float thedistance=0;
		float currentsize = 0;
		
		int *AAConstant=new int[GXYSIZE];
		int *rowsum= new int [FXYSIZE];
		int *colsum= new int [FXYSIZE];
		
		//Initialize AAConstant, rowsum,colsum
		for (int i=0;i<FXYSIZE;i++) {
			AAConstant[i]= 0;
			rowsum[i]= 0;
			colsum[i]= 0;
		};
		AAConstant[FXYSIZE]=0;
		AAConstant[FXYSIZE+1]=0;
		//Finished init
		
		// initilize local Gxy matrix
		for (int AA1=0;AA1<GXYSIZE;AA1++) {
			for (int AA2=0;AA2<GXYSIZE;AA2++) {
				Gxy[GXYSIZE*AA1+AA2]=0;
			}
		}
		
		for (int rateclass=0;rateclass<ClassNum;rateclass++) {
			int classize = 0;
			int currentconstants = 0;
			for (int site=0;site<SiteNum;site++) {
				if (rateclass==RateClasses[site]) {
					Gxy[GXYSIZE*theTCM[SiteNum*taxon1+site]+theTCM[SiteNum*taxon2+site]] += bootweight[site]; //count the different matches, including - and X
					classize+=bootweight[site];
				}
				if (Constants[site]) {
					AAConstant[theTCM[SiteNum*taxon1+site]]+=bootweight[site];	// Only used for rateclass=0, no reinit necessary - for now!!!!
				}
			}
	
			for (int AA1=0;AA1<FXYSIZE;AA1++) {
				for (int AA2=0;AA2<FXYSIZE;AA2++) {
					rowsum[AA1]+=Gxy[GXYSIZE*AA1+AA2];
					colsum[AA1]+=Gxy[GXYSIZE*AA2+AA1];
				}
			}
			for (int AA1=0;AA1<FXYSIZE;AA1++) {
				for (int AA2=0;AA2<FXYSIZE;AA2++) {
					float singlemiss = colsum[AA2] ? (float) (Gxy[GXYSIZE*GAP+AA2]+Gxy[GXYSIZE*AMBIG+AA2])/ (float) colsum[AA2] : 0;
					singlemiss +=  rowsum[AA1] ? (float) (Gxy[GXYSIZE*AA1+GAP]+Gxy[GXYSIZE*AA1+AMBIG])/ (float) rowsum[AA1] : 0;
					if (rateclass==0 && AA1==AA2) {
						Fxy[FXYSIZE*AA1+AA2]=(float(Gxy[GXYSIZE*AA1+AA2])-float(AAConstant[AA1])*InvariantFraction)*(1.0+singlemiss);
						currentconstants+=AAConstant[AA1];
					} else {
						Fxy[FXYSIZE*AA1+AA2]=float(Gxy[GXYSIZE*AA1+AA2])*(1.0+singlemiss);
					}
				}
				if (Fxy[FXYSIZE*AA1+AA1] == 0) {
					Fxy[FXYSIZE*AA1+AA1]=LIMITFRACTION;
				}
			}
			thedistance += (float(classize)-float(currentconstants)*InvariantFraction) * LogDet();
			currentsize += (float(classize)-float(currentconstants)*InvariantFraction);
		}
		thedistance/=currentsize;
		PairWise[pair]=thedistance;
		pair++;
	}
}
delete [] bootweight;
return pair;
}

inline int TCM::AAtoNum (char residue) {
int code;
	switch (residue) {
		case 'A' :
			code = 0;
			break;
		case 'R' :
			code = 1;
			break;
		case 'N' :
			code = 2;
			break;
		case 'D' :
			code = 3;
			break;
		case 'C' :
			code = 4;
			break;
		case 'E' :
			code = 5;
			break;
		case 'Q' :
			code = 6;
			break;
		case 'G' :
			code = 7;
			break;
		case 'H' :
			code = 8;
			break;
		case 'I' :
			code = 9;
			break;
		case 'L' :
			code = 10;
			break;
		case 'K' :
			code = 11;
			break;
		case 'M' :
			code = 12;
			break;
		case 'F' :
			code = 13;
			break;
		case 'P' :
			code = 14;
			break;
		case 'S' :
			code = 15;
			break;
		case 'T' :
			code = 16;
			break;
		case 'W' :
			code = 17;
			break;
		case 'Y' :
			code = 18;
			break;
		case 'V' :
			code = 19;
			break;
		case 'X' :
			code = 20;
			break;
		case '-' :
			code = 21;
			break;
		default :

			throw IllegalSymbol();
	}
	return code;
}

inline int TCM::DNAtoNum (char nucleotide) {
int code;
	switch (nucleotide) {
		case 'A' :
			code = 0;
			break;
		case 'C' :
			code = 1;
			break;
		case 'G' :
			code = 2;
			break;
		case 'T' :
			code = 3;
			break;
		case 'S' :
			code = 4;
			break;
		case 'W' :
			code = 4;
			break;
		case 'K' :
			code = 4;
			break;
		case 'M' :
			code = 4;
			break;
		case 'R' :
			code = 4;
			break;
		case 'Y' :
			code = 4;
			break;
		case 'B' :
			code = 4;
			break;
		case 'D' :
			code = 4;
			break;
		case 'H' :
			code = 4;
			break;
		case 'V' :
			code = 4;
			break;
		case 'N' :
			code = 4;
			break;
		case 'X' :
			code = 4;
			break;
		case '-' :
			code = 5;
			break;
		default :

			throw IllegalSymbol();
	}
	return code;
}
