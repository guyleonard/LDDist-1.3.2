const char* GetVersion ();

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
