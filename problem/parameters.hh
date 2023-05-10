#ifndef PROBLEM_PARAMETERS_HH_
#define PROBLEM_PARAMETERS_HH_

/* IDs:
 * m -> no. of dissolved species
 * n -> no. of mineral species
 * Ns -> total no. of species (Ns=m+n)
 * Ne -> total number of equilibrium reactions
 * Nk -> total number of kinetic reactions
 * nPDE -> no. of PDEs in the model (nPDE = Ns-Ne)
 * nEQC -> no. of equality constraints
 * nIQC -> no. of inequalita constraints (arising from eqb precipitation/dissolution of mineral species)
 * nEQC+nIQC=Ne
 */
class Tags{

public:

	/*************************************/
	// REACTION NETWORK AND VARIABLES TAGS
	/*************************************/
	// SPECIES
	const static int m 	= 7; // aqueous species
	const static int n 	= 1; // solid species
	constexpr static int Ns	= m+n; // total no. of species
	//EQUILIBRIUM REACTIONS (includes aq. as well as mineral reactions, i.e., Ne=nEQC+nIQC)
	const static int Ne	= 2;
	const static int Re_AB1 = 0; // CO2 + H2O <-> HCO3_1n + H_1p
	const static int Re_AB2 = 1; // HCO3_1n <-> CO3_2n + H_1p
	//KINETIC REACTIONS
	const static int Nk = 3;
	const static int Rk_OMdegSO4	= 0; // OM degradation with sulfate reduction:  	OM + (a/2) SO4_2n --> (a/2+b-2c) HCO3_1n + (a/2) HS_1n + b NH4_1p
	const static int Rk_OMdegCH4	= 1; // OM degradation with methane production: 	OM + (b-2c) H2O --> (b-2c) HCO3_1n + (a/2) CH4 + (a/2-b+2c) CO2 + b NH4_1p + c HPO4_2n
	const static int Rk_AOM			= 2; // Anaerobic Oxidation of Methane:				CH4 + SO4_2n --> HCO3_1n + HS_1- + H2O
	//
	constexpr static int nIQC = 0;  // no. of inequality constraints (i.e., eqb. mineral precipitation-dissolution)
	constexpr static int nEQC = Ne; // no. of equality constraints (i.e., eqb. aqueous reactions. e.g., acid-base equilibria)
	constexpr static int nPDE = Ns-nEQC;
	constexpr static int nPrC = nPDE+nIQC; 	// num. of primary species
	constexpr static int nEqns = nPrC+4; 	// +4 for [Pw, Sh, Sg, T]
	// VARIABLES TAGS
	/* !!ORDERING IS IMPORTANT!!
	 * order:
	 * 1. aq primary species,
	 * 2. solid primary species,
	 * 3. non-eqb mineral species,
	 * 4. eqb. mineral species,
	 * 5. Pw, Sh, Sg, T
	 **/
	//primary species -> aqueous
	const static int SO4_2n		= 0;
	const static int NH4_1p		= 1;
	const static int CH4		= 2;
	const static int CO2		= 3;
	const static int Cl_1n		= 4;
	//primary species -> solid
	const static int OM			= 5;
	//non-eqb. mineral (solid) species
	//none
	//eqb. mineral (solid) species -> (IQC) - pseudo-primary
	//none
	//secondary (aq.) species -> (EQC)
	const static int CO3_2n		= 7;
	const static int HCO3_1n	= 6;
	// phase primary variables
	const static int Pw = nPDE+nIQC+0;
	const static int Sh = nPDE+nIQC+1;
	const static int Sg = nPDE+nIQC+2;
	const static int T  = nPDE+nIQC+3;

	//post-processing variable order
	const static int NPP = Ns+12;
	// 0-to-Ns-1 --> all species
	// Ns onwards --> other vars, as below
	const static int PP_pw = Ns+0;
	const static int PP_pg = Ns+1;
	const static int PP_pc = Ns+2;
	const static int PP_pe = Ns+3;
	const static int PP_sh = Ns+4;
	const static int PP_sg = Ns+5;
	const static int PP_sw = Ns+6;
	const static int PP_T  = Ns+7;
	const static int PP_K  = Ns+8;
	const static int PP_sal= Ns+9;
	const static int PP_zf = Ns+10;
	const static int PP_CeqCH4= Ns+11;
	const static int PP_CeqHYD= Ns+12;

	//VERIFY THE TAGS
	void check(){

		if ( Ne != nEQC+nIQC ){
			std::cout<< "ERROR: Ne != nEQC+nIQC" << std::endl;
			exit(0);
		}

		std::cout<< "All TAGS are verified!" << std::endl;
	}

	//GET PRINT NAME FOR THE SPECIES
	std::string getSpeciesName( int tag )const{

		std::string name;

		switch(tag) {
			case Cl_1n:
				name = "Cl^1-";
				break;
			case NH4_1p:
				name = "NH4^1+";
				break;
			case SO4_2n:
				name = "SO4^2-";
				break;
			case CH4:
				name = "CH4";
				break;
			case CO2:
				name = "CO2";
				break;
			case CO3_2n:
				name = "CO3^2-";
				break;
			case HCO3_1n:
				name = "HCO3^1-";
				break;
			case OM:
				name = "organic_matter";
				break;
		}
	    return name;
	}

	//GET PRINT NAME FOR THE KINETIC REACTIONS
	std::string getKineticReactionName( int tag )const{

		std::string name;

		switch(tag) {
    		case Rk_OMdegSO4:
				name = "OM-degradation-with-SO4-reduction";
				break;
    		case Rk_OMdegCH4:
				name = "OM-degradation-with-CH4-production";
				break;
    		case Rk_AOM:
    			name = "Anaerobic-Oxidation-of-Methane";
		}
	    return name;
	}

	//GET PRINT NAME (VTK) FOR THE PRIMARY VARIABLES
	std::string getPrimaryVariablesName( int tag )const{

		std::string name;

		switch(tag) {
			case NH4_1p:
				name = "NH4^1+";
				break;
			case Cl_1n:
				name = "Cl^1-";
				break;
			case SO4_2n:
				name = "SO4^2-";
				break;
			case CH4:
				name = "CH4";
				break;
			case CO2:
				name = "CO2";
				break;
			case OM:
				name = "organic_carbon";
				break;
			case Pw:
				name = "pw";
				break;
			case Sh:
				name = "sh";
				break;
			case Sg:
				name = "sg";
				break;
			case T:
				name = "T";
				break;
		}
	    return name;
	}

	// BOUNDARY TYPE TAGS
	const static int neumann	= 0;
	const static int dirichlet	= 1;
};

class ConversionFactors{

public:

	constexpr static double F_Mkg 	= 1.01/1000.;
	constexpr static double F_Corg 	= (1./100.) * (2.5/12.0) * 1.e6;

	// atomin C:N:P ratio
	constexpr static double CNP_a = 106.0;
	constexpr static double CNP_b = 16.0;
	constexpr static double CNP_c = 1.0;

	double Fpw( double porosity ) const {

		double F_pw = (1.-porosity)/porosity;
		return F_pw;
	}

	double Fds( double porosity ) const {

		double F_ds = porosity/(1.-porosity);
		return F_ds;
	}

};

class Matrices{

public:

	// STOICHIOMETRIC MATRIX FOR EQUILIBRIUM REACTIONS (size: NexNs)
	std::vector< std::vector<double> > S_eqb() const {

		int nrow = Tags::Ne;
		int ncol = Tags::Ns;
		if( Tags::Ne == 0 ){ //dummy values
			nrow = 1;
			ncol = 1;
		}
		std::vector< std::vector<double> > S(nrow);
		for( int i=0;i<nrow;i++ ){
		  S[i] = std::vector<double>(ncol,0.);
		}

		int r = Tags::Re_AB1;
		S[r][Tags::CO2] 	= -1.;
		S[r][Tags::HCO3_1n]	= 1.;

		r = Tags::Re_AB2;
		S[r][Tags::HCO3_1n] = -1.;
		S[r][Tags::CO3_2n]	= 1.;

		return S;
	}

	// STOICHIOMETRIC MATRIX FOR KINETIC REACTIONS (size: NkxNs)
	std::vector< std::vector<double> > S_kin() const {

		std::vector< std::vector<double> > S(Tags::Nk);
		for( int i=0;i<Tags::Nk;i++ ){
		  S[i] = std::vector<double>(Tags::Ns,0.);
		}

		// atomin C:N:P ratio
		ConversionFactors factor;
		double a = factor.CNP_a;
		double b = factor.CNP_b;
		double c = factor.CNP_c;

		// OM degradation with sulfate reduction
		// OM + (a/2) SO4_2n --> (a/2+b-2c) HCO3_1n + (a/2) HS_1n + b NH4_1p + c HPO4_2n + (a/2-b+2c) CO2 + (a/2-b+2c) H2O
		int r = Tags::Rk_OMdegSO4;
		S[r][Tags::OM] 		= -a;
		S[r][Tags::SO4_2n] 	= -((1./2.)*a);
		S[r][Tags::HCO3_1n]	= ((1./2.)*a+b-2.*c);
		S[r][Tags::NH4_1p] 	= (b);
		S[r][Tags::CO2] 	= ((1./2.)*a-b+2.*c);

		// OM degradation with methane production
		// OM + (b-2c) H2O --> (b-2c) HCO3_1n + (a/2) CH4 + (a/2-b+2c) CO2 + b NH4_1p + c HPO4_2n
		r = Tags::Rk_OMdegCH4;
		S[r][Tags::OM] 		= -a;
		S[r][Tags::HCO3_1n]	= (b-2.*c);
		S[r][Tags::CH4] 	= ((1./2.)*a);
		S[r][Tags::CO2] 	= ((1./2.)*a-b+2.*c);
		S[r][Tags::NH4_1p] 	= (b);

		// Anaerobic oxidation of methane
		// CH4 + SO4_2n --> HCO3_1n + HS_1- + H2O
		r = Tags::Rk_AOM;
		S[r][Tags::CH4] 	= -1.;
		S[r][Tags::SO4_2n] 	= -1.;
		S[r][Tags::HCO3_1n]	= 1.;

		return S;

	}

	// MOBILITY MATRIX
	/* diagonal matrix of size NsxNs
	 * M[i][i] = 1 for i corresponding to dissolved species
	 * M[i][i] = 0 for i corresponding to solid species
	 * */
	std::vector< std::vector<double> > M() const {

		std::vector< std::vector<double> > m(Tags::Ns); //diagonal matrix, Ns rows, Ns columns
		for( int i=0;i<Tags::Ns;i++ ){
		  m[i] = std::vector<double>(Tags::Ns,1.); //default: all species are dissolved
		}

		for(int i=0; i<Tags::Ns; i++){
			if( i==Tags::OM ) m[i][i] = 0.;
		}

		return m;
	}

};

template<typename GV,typename PTree>
class Domain{
private:

	const GV& gv;
	const PTree& ptree;
	const static int dim = GV::dimension;
	constexpr static double eps = 1.0e-6;

	double Xmax;
	double Zmax;
#if DIMENSION==3
	double Ymax;
#endif

public:

	Domain (const GV& gv_ ,
			const PTree& ptree_)
	: gv( gv_ ),
	  ptree(ptree_)
	{
		Xmax = ptree.get("domain.yasp.Xmax",(double)1.0);//m
		Xmax *= 100.0;	//convert to cm
#if DIMENSION==2
		Zmax = ptree.get("domain.yasp.Zmax",(double)400.0);//m
		Zmax *= 100.0;	//convert to cm
#elif DIMENSION==3
		Ymax = ptree.get("domain.yasp.Ymax",(double)1.0);//m
		Ymax *= 100.0;	//convert to cm
		Zmax = ptree.get("domain.yasp.Zmax",(double)400.0);//m
		Zmax *= 100.0;	//convert to cm
#endif

	}

	double SeaFloorHeight( double xpos )const{
		return Zmax; //cm
	}

	bool isTopBoundary( Dune::FieldVector< double, dim > globalPos/*cm*/ ) const {
		if( globalPos[dim-1] > Zmax/*cm*/ - eps ){
			return true;
		}
		else return false;
	}

	bool isBottomBoundary( Dune::FieldVector< double, dim > globalPos/*cm*/ ) const {
		if( globalPos[dim-1] < 0. + eps ){
			return true;
		}
		else return false;
	}

};

template<typename GV,typename PTree>
class Properties
{
private:
	const GV& gv;
	const PTree& ptree;
	double *time;
	double *dt;

	const static int dim = GV::dimension;

	ConversionFactors factor;

	double Xmax;
	double Zmax;
#if DIMENSION==3
	double Ymax;
#endif

	double v_0;
	double v_inf;
	double phi_0;
	double phi_inf;
	double exponent_beta;

	double permeability;
	double swr, sgr;
	double lambda_BC, pe_BC, pcmax;
	double exponent_m;

	bool gravity_flag;
	double gravity_magnitude;

	double Ru;

	double ch4_accentricity;
	double ch4_Tcr;
	double ch4_Pcr;
	double h2o_Tcr;
	double h2o_Pcr;

	double Mch4;
	double Mh2o;
	double Mhyd;
	double Nhyd;

	double kr_hyddissociation;
	double kr_hydformation;
	double kr_hyddissolution;
	double kr_hydprecipitation;

	double reference_salinity;

	double fcorg;
	double bwc_Cl_1n;
	double bwc_NH4_1p;
	double bwc_SO4_2n;
	double bwc_CH4;
	double bwc_CO2;
	double bwc_OM;

	double D_Cl_1n;
	double D_NH4_1p;
	double D_SO4_2n;
	double D_CH4;
	double D_CO2;
	double D_CO3_2n;
	double D_HCO3_1n;

	double Db0;
	double zb;

	double Keq_AB1;
	double Keq_AB2;

	double k_AOM;
	double K_SO42n;
	double K_C;
	double a0;
	double kx0;
	double a1;

	double Pw_SF;
	double rhow_ref;
	double T_SF;
	double gradT;

public:

	Properties (const GV& gv_ ,
				const PTree& ptree_,
				double *time_,
				double *dt_)
	: gv( gv_ ),
	  ptree(ptree_),
	  time(time_),
	  dt(dt_)
	{
		Xmax = ptree.get("domain.yasp.Xmax",(double)1.0);//m
		Xmax *= 100.0;	//convert to cm
#if DIMENSION==2
		Zmax = ptree.get("domain.yasp.Zmax",(double)400.0);//m
		Zmax *= 100.0;	//convert to cm
#elif DIMENSION==3
		Ymax = ptree.get("domain.yasp.Ymax",(double)1.0);//m
		Ymax *= 100.0;	//convert to cm
		Zmax = ptree.get("domain.yasp.Zmax",(double)400.0);//m
		Zmax *= 100.0;	//convert to cm
#endif

		v_0 = 0.; 															// read in cm/yr
		v_inf = ptree.get("parameters.burial_velocity",(double)0.0);		// read in cm/yr
		phi_0 = ptree.get("parameters.porosity.phi_0",(double)0.78); 		//--
		phi_inf = ptree.get("parameters.porosity.phi_inf",(double)0.65);	//--
		exponent_beta = ptree.get("parameters.porosity.beta",(double)0.2);	//--

		permeability  = ptree.get("parameters.permeability.K0",(double)1.e-16);			// read in m^2
		permeability *= 100.*100.; //convert to cm^2
		swr = ptree.get("parameters.capillary_pressure.swr",(double)0.); 				//--
		sgr = ptree.get("parameters.capillary_pressure.sgr",(double)0.);				//--
		lambda_BC = ptree.get("parameters.capillary_pressure.lambda",(double)1.2);		//--
		pe_BC  = ptree.get("parameters.capillary_pressure.entry_pressure",(double)0.0); // read in Pa
		pe_BC *= 1.e-6; 																// converted to MPa
		pcmax  = ptree.get("parameters.capillary_pressure.pc_max",(double)1.e8) ;	 	// read in Pa
		pcmax *= 1.e-6; 																// converted to MPa
		exponent_m = ptree.get("parameters.permeability.m",(double)3.0);				//--

		gravity_flag = true;
		gravity_magnitude  = 9.81; // Pa.m^2/kg
		gravity_magnitude *= 1.e-8;// converted to MPa.L/g.cm

		/* http://en.wikipedia.org/wiki/Gas_constant */
		Ru = 8.314462175; /* [J*mol^-1*K^-1] */
		ch4_accentricity = 0.011;
		ch4_Tcr = -82.7 + 273.15 ; 	/* [K]  */
		ch4_Pcr = 45.96 * 1.0e5 ; 	/* [Pa] */
		h2o_Tcr = 647.096 ; 		/* [K]  */
		h2o_Pcr = 22.064 * 1.0e6 ; 	/* [Pa] */

		Mch4 = 16.04 * 1.0e-3; 	/* [g/mmol] */
		Mh2o = 18.0  * 1.0e-3; 	/* [g/mmol] */
		Mhyd = 119.5 * 1.0e-3; 	/* [g/mmol] */
		Nhyd = 5.9;				//--

		kr_hyddissociation  = ptree.get("parameters.hydrate_phase_change.dissociation_rate",(double)0.); 	/*mol/m².Pa.s*/
		kr_hyddissociation *= 364.25*86400.0*1.e8;															// convert to (cm^3/L)*(mmol/yr.cm^2.MPa)
		kr_hydformation 	= ptree.get("parameters.hydrate_phase_change.formation_rate",(double)0.); 		/*mol/m².Pa.s*/
		kr_hydformation	   *= 364.25*86400.0*1.e8;															// convert to (cm^3/L)*(mmol/yr.cm^2.MPa)
		kr_hyddissolution   = ptree.get("parameters.hydrate_phase_change.dissolution_rate",(double)0.02); 	/*TODO*/
		kr_hyddissolution  *= 1.e-7;																		// convert to TODO
		kr_hydprecipitation = ptree.get("parameters.hydrate_phase_change.precipitation_rate",(double)0.02); 	/*TODO*/
		kr_hydprecipitation*= 1.e-7;																		// convert to TODO

		reference_salinity = ptree.get("parameters.reference_state.salinity",(double)35.0); /*g/kg*/

		fcorg = factor.F_Corg;
		bwc_Cl_1n 	= ptree.get( "parameters.bottom_water_concentration.Cl_1n"	, (double)555.0	);  //read in mmol/Lpw
		bwc_NH4_1p 	= ptree.get( "parameters.bottom_water_concentration.NH4_1p"	, (double)0.003	);  //read in mmol/Lpw
		bwc_SO4_2n 	= ptree.get( "parameters.bottom_water_concentration.SO4_2n"	, (double)30.0	);  //read in mmol/Lpw
		bwc_CH4 	= ptree.get( "parameters.bottom_water_concentration.CH4"	, (double)0.0	);  //read in mmol/Lpw
		bwc_CO2 	= ptree.get( "parameters.bottom_water_concentration.CO2"	, (double)5.28e-2); //read in mmol/Lpw
		bwc_OM 		= ptree.get( "parameters.bottom_water_concentration.OM"		, (double)1.0	);  //read in wt %
		bwc_OM 	   *= fcorg; 																		// convert to mmol/Lds

		D_Cl_1n 	= ptree.get( "parameters.aqueous_diffusion.Cl_1n"	, (double)188.11); //read in cm^2/yr
		D_NH4_1p 	= ptree.get( "parameters.aqueous_diffusion.NH4_1p"	, (double)188.11); //read in cm^2/yr
		D_SO4_2n 	= ptree.get( "parameters.aqueous_diffusion.SO4_2n"	, (double)188.11); //read in cm^2/yr
		D_CH4 		= ptree.get( "parameters.aqueous_diffusion.CH4"		, (double)188.11); //read in cm^2/yr
		D_CO2 		= ptree.get( "parameters.aqueous_diffusion.CO2"		, (double)188.11); //read in cm^2/yr
		D_CO3_2n 	= ptree.get( "parameters.aqueous_diffusion.CO3_2n"	, (double)188.11); //read in cm^2/yr
		D_HCO3_1n 	= ptree.get( "parameters.aqueous_diffusion.HCO3_1n"	, (double)188.11); //read in cm^2/yr

		Db0 		= ptree.get( "parameters.bioturbation.Db0", (double)300.0); //read in cm^2/yr
		zb			= ptree.get( "parameters.bioturbation.zb" , (double)2.0); 	//read in cm

		Keq_AB1 = ptree.get( "parameters.equilibrium_constant.AB1", (double) 1.3139e-6) ; //TODO
		Keq_AB2 = ptree.get( "parameters.equilibrium_constant.AB2", (double) 6.0940e-10); //TODO

		k_AOM 		= ptree.get( "parameters.kinetic_rate.AOM"					, (double)0.001	); //L/(mmol.yr)
		K_SO42n 	= ptree.get( "parameters.kinetic_rate.K_SO42n"				, (double)1.0	); //mmol/Lpw
		K_C 		= ptree.get( "parameters.kinetic_rate.K_C"					, (double)40.0	); //mmol/Lpw
		a0 			= ptree.get( "parameters.kinetic_rate.POC_age_coefficient"	, (double)1000.0); //yr
		kx0 		= ptree.get( "parameters.kinetic_rate.POC_intrinsic_rate"	, (double)0.16	); //TODO
		a1 			= ptree.get( "parameters.kinetic_rate.POC_decreace_exponent", (double)0.95); //--

		rhow_ref = 1027.0; //g/L-->kg/m^3
		double hwc = ptree.get("parameters.water_column_height",(double)800.0);//in m
		Pw_SF = rhow_ref*9.81*hwc; //Pa
		T_SF = ptree.get("parameters.bottom_water_temperature",(double)4.0);//in degC
		T_SF += 273.15; //convert to K
		gradT = ptree.get("parameters.regional_thermal_gradient",(double)24.0);//in degC/km
		gradT *= 1./1000.0; // convert to degC/m
		gradT *= 1.e-2; //convert to degC/cm
	}

	double Salinity( double C_Cl/*mmol/Lpw*/ )const{
		return 0.00180665 * C_Cl/*mmol/Lpw*/ * 35.453/*g/mol*/; // in g/kg or ppt
	}

	/**********************************************************************
	 * PHASE PROPERTIES
	 * GAS:
	 * 		molar mass
	 * 		density,
	 * 		compressibility factor,
	 * 		dynamic viscosity
	 * 		thermal conductivity
	 * 		heat capacity Cp
	 * 		heat capacity Cv
	 * 		maximum solubility (i.e., equilibrium concentration: C_CH4_eqb)
	 * WATER:
	 * 		molar mass
	 * 		density,
	 * 		dynamic viscosity
	 * 		thermal conductivity
	 * 		heat capacity Cp
	 * 		heat capacity Cv(=Cp)
	 * HYDRATE:
	 * 		molar mass
	 * 		density
	 * 		thermal conductivity
	 * 		heat capacity Cp
	 * 		heat capacity Cv(=Cp)
	 * 		hydration number Nh
	 * SEDIMENT:
	 * 		density
	 * 		thermal conductivity
	 * 		heat capacity Cp
	 * 		heat capacity Cv(=Cp)
	 */

	// GAS PHASE
	double MolarMassCH4()const{
		return Mch4; 	/* [g/mmol] */
	}

	double CompressibilityFactorCH4( double Pg/*MPa*/, double T/*K*/ )const{
		// FrauenhoferCOMSOL model
		double C0 = 1.0003;
		double C1 = 2.1173e-8;
		double C2 = 5.0207e-16;
		double C3 = 7.5669e-23;
		double C4 = 1.3259e-30;

		Pg *= 1.e6; //convert to Pa
		double compressibilityFactor = 0.8;//std::min( C0 - C1*Pg - C2*Pg*Pg + C3*Pg*Pg*Pg - C4*Pg*Pg*Pg*Pg , 1.0 );

		return compressibilityFactor;//std::max( compressibilityFactor , 0.5 ) ; //--
	}

	double GasPhaseDensity( double Pg/*MPa*/, double zf/*--*/, double T/*K*/ )const{
		Pg *= 1.e6; //convert to Pa
		double rho/*kg/m^3*/ = Pg/( zf * (Ru/Mch4) * T); // Ru/Mch4--> dimensionality is consistent --> Ru: Pa.m^3/mol/K ; Mch4: g/mmol => kg/mol
		return rho; //g/L
	}

	double GasPhaseViscosity( double Pg/*MPa*/, double zf/*--*/, double T/*K*/ )const{
		Pg *= 1.e6; //convert to Pa.
		// Sutherland Correlation:
		// ref: http://portal.tpu.ru/SHARED/n/NATASHA/Material/Tab3/Glava_1.pdf
		double C = 162; // empirical constant
		double mu_0 = 1.0707e-5; //Pa.s
		// ref for mu_0 :http://www.pipeflowcalculations.com/tables/gas.php
//		mu_0 *= (  1. - (1./(1.0707e-5)) * (  4.8134e-14 * Pg
//											- 4.1719e-20 * Pg * Pg
//											+ 7.3232e-28 * Pg * Pg * Pg )
//				); // Pa.s -> ref: Frauenhofer Comsol Model
		double mu/*Pa.s*/ = mu_0;// * (273.15 + C ) * ( pow( (T/273.15), 1.5) / ( T + C ) ) ;
		mu *= 1.e-6 * (1./(364.25*86400.0)); //convert to MPa.yr
		return mu ; //MPa.yr
	}

	double GasPhaseThermalConductivity( double Pg/*MPa*/, double zf/*--*/, double T/*K*/ )const{
		Pg *= 1.e6; //convert to Pa.
		// REFERENCE: " Thermal Conductivity of Methane for temperatures between 110 K and 310 K with Pressures upto 70 MPa
		//			  author : H.M. Roder
		//			  Journal: Journal of Thermophysics, volume 6, No 2, pages 119-142
		// assumption: dilute gas , therefore, density effects are neglected
		double A0 = - 0.8863333440 * 1.e-2 ;
		double A1 =   0.2419639784 * 1.e-3 ;
		double A2 = - 0.6997019196 * 1.e-6 ;
		double A3 =   0.1224609018 * 1.e-8 ;
		double kth/*[W/(m*K)]*/ = A0 + A1 * T + A2 * T*T + A3 * T*T*T ;
		kth *= (364.25*86400.0*1000.0)/(100.0); //convert
		return kth; // (J/(a*cm*K)) * (cm^3/L)
	}

	double GasPhaseCp( double Pg/*MPa*/, double zf/*--*/, double T/*K*/ )const{
		Pg *= 1.e6; //convert to Pa.

		// IDEAL Cp
		/* REF: 1D Modelling of Hydrate Decomposition in Porous Media, by F. Esmailzadeh, M.E. Zeighami, J. Fathi */
		double Ai = 1.238 ;
		double Bi = 0.00313 ;
		double Ci = 7.905*1.0e-7 ;
		double Di = -6.858*1.0e-10 ;
		double Cp_ideal/*[J/(kg*K)]*/ = ( Ai + Bi*T + Ci*T*T +Di*T*T*T ) * 1000.0 ;
		Cp_ideal *= 1.e-3; //convert to J/(g*K)

//		// RESIDUAL Cp
//		// Based on Peng Robinson's EoS
//		double omega/*-*/= ch4_accentricity;
//		double Tc/*K*/	 = ch4_Tcr;
//		double Pc/*Pa*/	 = ch4_Pcr;
//		double kappa = 0.;
//		if( omega <= 0.49 ){
//			kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega ;
//		}
//		else{
//			kappa = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega ;
//		}
//		double ac = pow( (1 + kappa * ( 1 - sqrt(T/Tc) ) ) , 2 );
//		double b = 0.07780 * ( Ru * Tc / Pc );
//		double a_T = 0.45724 * ( ( Ru * Ru * Tc * Tc ) / Pc ) * ac ;
//		double da_T = kappa * ac * ( ( kappa / Tc ) - (( 1 + kappa )/( sqrt( T*Tc ) ) ) );
//		double dda_T = ( kappa * ac * ( 1. + kappa ) ) / ( 2. * sqrt( T*Tc ));
//		double A = ( a_T * Pg ) / ( pow( ( Ru * T ) , 2 ) ) ;
//		double B = ( b * Pg ) / ( Ru * T ) ;
//		double M = ( zf*zf + 2*B*zf - B*B ) / ( zf - B ) ;
//		double N = ( da_T * B ) / ( b * Ru );
//		double Cp_residual/*[J/(kg*K)]*/ =   dda_T * ( T / (2*sqrt(2) * b ) ) * log((zf+(sqrt(2)+1)*B)/(zf-(sqrt(2)-1)*B))
//				+ ( Ru * pow( M-N , 2 ) ) / ( M*M - 2.*A * (zf+B) )
//				- Ru ;
//		Cp_residual *= 1.e-3; //convert to J/(g*K)

		return Cp_ideal;// + Cp_residual; //J/(g*K)
	}

	double GasPhaseCv( double Pg/*MPa*/, double zf/*--*/, double T/*K*/ )const{
		double Cv/*[J/(g.K)]*/ =   GasPhaseCp( Pg/*MPa*/,zf,T ); // NOTE: Cases are checked in Cp_ideal function!

//		// Based on Peng Robinson's EoS
//		Pg *= 1.e6; //convert to Pa.
//		double omega = ch4_accentricity;
//		double Tc/*K*/	 = ch4_Tcr;
//		double Pc/*Pa*/	 = ch4_Pcr;
//		double kappa = 0.;
//		if( omega <= 0.49 ){
//			kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega ;
//		}
//		else{
//			kappa = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega ;
//		}
//		double ac = pow( (1 + kappa * ( 1 - sqrt(T/Tc) ) ) , 2 );
//		double b = 0.07780 * ( Ru * Tc / Pc );
//		double a_T = 0.45724 * ( ( Ru * Ru * Tc * Tc ) / Pc ) * ac ;
//		double da_T = kappa * ac * ( ( kappa / Tc ) - (( 1 + kappa )/( sqrt( T*Tc ) ) ) );
//		double dda_T = ( kappa * ac * ( 1. + kappa ) ) / ( 2. * sqrt( T*Tc ));
//		double A = ( a_T * Pg ) / ( pow( ( Ru * T ) , 2 ) ) ;
//		double B = ( b * Pg ) / ( Ru * T ) ;
//		double M = ( zf*zf + 2*B*zf - B*B ) / ( zf - B ) ;
//		double N = ( da_T * B ) / ( b * Ru );
//
//		double nonidealfactor = ( pow( M-N , 2 ) ) / ( M*M - 2.*A * (zf+B) );
//		Cv += (-1.) * Ru * nonidealfactor * (1.e-3/*conv. factor*/) ;

		return Cv;  /* [J/(g*K)] */
	}

	double EquilibriumConcentrationCH4( const typename GV::Traits::template Codim<0>::Entity& element,
			 	 	 	 	 	 	    const Dune::FieldVector<double,dim>& xlocal,
										double Pg/*MPa*/, double zf/*--*/, double T/*K*/, double S/*g/kg*/ )const{
//		Pg *= 1.e6; //convert to Pa.
//		double Tr = T/100.;
//		double Ceq = 1.e-3*std::exp( -417.5053 + 599.8626*(1./Tr) + 380.3636*std::log(Tr) - 62.0764*Tr
//									 + S * (-0.06423 + 0.03498*Tr - 0.0052732*Tr*Tr )
//									);
//		Ceq *= 20.*1./GasPhaseDensity(Pg/1.e6,zf,T);

		Pg *= 1.e6; //convert to Pa.

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);//cm
		double Pg0/*Pa*/ = Pw_SF/*Pa*/
						 + rhow_ref/*g/L*/
						 	 * abs(g()[dim-1])/*(MPa.L)/(g.cm)*/
		  	  	  	   	 	 * (Zmax/*cm*/-x[dim-1]/*cm*/)
							 * 1.e6 /*MPa to Pa*/;
		double T0/*K*/ = T_SF/*K*/ + gradT/*degC/cm*/*(Zmax/*cm*/-x[dim-1]/*cm*/);
		double Cl0/*mmol/Lpw*/ = BottomWaterConcentration(element,xlocal,0.,0.,Tags::Cl_1n);
		double S0 = Salinity(Cl0);


		double Tr = T/190.6;
		double Pr = Pg/(46.41e5);
		double Vr = zf*Tr/Pr;
		double a1=0.0872553928,
			   a2=-0.752599476,
			   a3=0.375419887,
			   a4=0.0107291342,
			   a5=0.00549626360,
			   a6=-0.0184772802,
			   a7=3.18993183e-4,
			   a8=2.11079375e-4,
			   a9=2.01682801e-5,
			   a10=-1.65606189e-5,
			   a11=1.19614546e-4,
			   a12=-1.08087289e-4,
			   Alpha=0.0448262295,
			   Beta=0.75397,
			   Gamma=0.077167;
		double B = a1 + a2/(Tr*Tr) + a3/(Tr*Tr*Tr);
		double C = a4 + a5/(Tr*Tr) + a6/(Tr*Tr*Tr);
		double D = a7 + a8/(Tr*Tr) + a9/(Tr*Tr*Tr);
		double E = a10 + a11/(Tr*Tr) + a12/(Tr*Tr*Tr);
		double F = Alpha/(Tr*Tr*Tr);
		double G = F/(2.*Gamma)*( Beta + 1.0 - (Beta + 1.0 + (Gamma/Vr*Vr) ) * std::exp(-Gamma/(Vr*Vr)) );
		double lnPHI = zf - 1.0 - std::log(zf) + B/Vr + C/(2.0*Vr*Vr) + D/(4.0*Vr*Vr*Vr*Vr) + E/(5.0*Vr*Vr*Vr*Vr*Vr) + G;
		double PHI = std::exp(lnPHI);

		double K = std::exp( -417.5053 + 599.8626*100./T + 380.3636*std::log(T/100.) - 62.0764*T/100. + S*( -0.06423 + 0.03498*T/100. - 0.0052732*(T/100.)*(T/100.)) );

		double Ceq = ((1000.0-S)/1000.) * WaterPhaseDensity(Pg,T,S) * std::exp(  std::log( PHI*Pg/1.e5)
																					- std::log(0.996*1.01325)
																   	   	 	   	  	+ std::log(1.e-9*K)
																   	   	   	   	   	- ( 37.0 + 0.04*(T0-298.15)*(Pg0*1.e-5-1.) ) / (T*83.14467)
																					- std::log(1.0-0.00100511*S)
																  	  	      	  );

//		std::cout<< lnPHI << '\t'
//				 << PHI << '\t'
//				 << K << '\t'
//				 << Ceq << std::endl;

		return 0.715*Ceq; //mmol/Lpw
	}

	// WATER PHASE
	double MolarMassH2O()const{
		return Mh2o; //g/mmol
	}

	double WaterPhaseDensity( double Pw/*MPa*/, double T/*K*/, double S/*g/kg*/ )const{
		Pw*= 1.e6; //convert to Pa.
		/* averages values & expansion coefficients: ρ0=1027 kg/m^3,  T0=10°C,  S_0=35 g/kg
		 * Thermal expansion: \alpha_T=0.15 kg/(m^3 °C)
		 * Salinity contraction: \alpha_S=0.78 kg/(m^3 kg/g)
		 * Pressure compressibility: \alpha_P=0.0045 kg/(m^3 dbar)
		 * UNESCO EOS-80 : Equation of state for seawater
		 * We use a linear EOS (web ref:http://www.ccpo.odu.edu/~atkinson/OEAS405/Chapter2_Stratified_Ocean/Lec_04_DensityEOS.pdf)
		 */
		double rho_0/*kg/m^3*/ = 1027.0;
		double T_0/*degC*/ 	   = 10.;
		double S_0/*g/kg*/	   = 35.0;
		double alpha_T = -0.15;
		double alpha_S = 0.78;
		double alpha_P = 0.0045;

		double rho/*kg/m^3*/ = rho_0;// + ( alpha_P*(Pw*1.e-4) + alpha_T*((T-273.15)-T_0) + alpha_S*(S-S_0) );
		return rho; //g/L
	}

	double WaterPhaseViscosity( double Pw/*MPa*/, double T/*K*/, double S/*g/kg*/ )const{
		Pw*= 1.e6; //convert to Pa.
		double mu_0 = 0.001792 ; // kg/m/s
		double a = - 1.94 ;
		double b = - 4.80 ;
		double c =  6.74 ;
		double T0 = 273.15 ; // K
		double Tr = T0/T;
		double mu/*Pa.s*/ = mu_0 * exp( a + b * Tr + c * Tr*Tr );
		return mu * ( 1.e-6 * (1./(364.25*86400.0)) ); //MPa.a
	}

	double WaterPhaseThermalConductivity( double Pw/*MPa*/, double T/*K*/, double S/*g/kg*/ )const{
		Pw*= 1.e6; //convert to Pa.
		double kth/*W.m^-1 K^-1*/ = 0.57153;//*( 1 + 0.003*(T-273.15) - 1.025e-5*(T-273.15)*(T-273.15) + 6.53e-10*Pw - 0.29*S );
		return kth * (( 364.25*86400.0*1000.0)/(100.0)); // (J/(a*cm*K))*(cm^3/L)
	}

	double WaterPhaseCp( double Pw/*MPa*/, double T/*K*/, double S/*g/kg*/ )const{
		Pw*= 1.e6; //convert to Pa.
		double Cp/*J*kg^-1*K^-1*/ = 3945.0 ;
		return Cp * 1.e-3; //J/(g*K)
	}

	double WaterPhaseCv( double Pw/*MPa*/, double T/*K*/, double S/*g/kg*/ )const{
		return WaterPhaseCp(Pw/*MPa*/,T,S); //J/(g*K)
	}

	// HYDRATE PHASE
	double MolarMassHYD()const{
		return Mhyd; //g/mmol
	}

	double HdrationNumber()const{
		return Nhyd; //--
	}

	double HydratePhaseDensity( double Peff/*MPa*/, double T/*K*/ )const{
		return 920.0; //g/L
	}

	double HydatePhaseThermalConductivity( double Peff/*MPa*/, double T/*K*/ )const{
		return 0.5 * ((364.25*86400.0*1000.0)/(100.0)); // (J/(a*cm*K))*(cm^3/L)
	}

	double HydratePhaseCv( double Peff/*MPa*/, double T/*K*/ )const{
		double Cp/*J/kg.K*/ = ( 1.9370547e-05*T*T*T - 1.5151760e-02*T*T + 3.9553876*T - 342.70565 )*1.0e3;
		return Cp * 1.e-3; //J/(g*K)
	}

	double EquilibriumConcentrationHydrate( double pw/*MPa*/,
											double T/*K*/,
											double S/*g/kg*/)const{
		pw *= 1.e6; //convert to Pa

		double r=0.;

		double C_eq_h = std::exp( -2.5640213e5 - 1.6448053e2*T + 9.1089042e-2*T*T
								  + 4.90352929e6/T + 4.93009113e4*std::log(T)
								  + S*( - 614.627279
										- 0.40619264*T
										+ 2.305936e-4*T*T
										+ 11.5242524e3/T
										+ 118.748089*std::log(T)
								  	  )
								 );
		C_eq_h *= 0.715*WaterPhaseDensity(pw*1.e-6,T,S);

		return C_eq_h; //mmol/Lpw
	}

	// SEDIMENT PHASE
	double SedimentDensity( double Peff/*MPa*/, double T/*K*/ )const{
		return 2600.0;//g/L
	}

	double SedimentThermalConductivity( double Peff/*MPa*/, double T/*K*/ )const{
		return 3.0 * ((364.25*86400.0*1000.0)/(100.0)); // (J/(a*cm*K))*(cm^3/L)
	}

	double SedimentCv( double Peff/*MPa*/, double T/*K*/ )const{
		return 1000.0 * 1.e-3; //J/(g*K)
	}

	/**********************************************************************
	 * HYDRATE REACTION KINETICS
	 * source terms: q_g, q_w, q_h
	 */

	double HPDReactionRate( double pw/*MPa*/,double T/*K*/, double S/*g/kg*/, double C_CH4/*mmol/Lpw*/, double Sh /*-*/, double por)const{
		pw *= 1.e6; //convert to Pa

		double r=0.;

		double k_HP = kr_hydprecipitation * 10*SedimentDensity(pw*1.e-6,T)*(1.-por)/(por*MolarMassHYD()); //rate of precipitation
		double k_HD = kr_hyddissolution * 10.*SedimentDensity(pw*1.e-6,T)*(1.-por)/(por*MolarMassHYD()); //rate of dissolution
		double C_eq_h = std::exp( -2.5640213e5 - 1.6448053e2*T + 9.1089042e-2*T*T
								  + 4.90352929e6/T + 4.93009113e4*std::log(T)
								  + S*( - 614.627279
										- 0.40619264*T
										+ 2.305936e-4*T*T
										+ 11.5242524e3/T
										+ 118.748089*std::log(T)
								  	  )
								 );
		C_eq_h *= 0.715*WaterPhaseDensity(pw*1.e-6,T,S);

		if( (C_CH4/C_eq_h - 1.0) > 0. ){
//			std::cout<< C_eq_h << '\t' << (C_CH4/C_eq_h - 1.0) << std::endl;
			r = k_HP*(1.-Sh)*(C_CH4/C_eq_h - 1.0);
		}
		else r = k_HP * Sh * (C_CH4/C_eq_h - 1.0);

		return r; //mmol/L/a
	}

	double HeatOfHydrateDissolution( double pw/*MPa*/,double T/*K*/, double S/*g/kg*/, double C_CH4/*mmol/Lpw*/, double Sh /*-*/, double por)const{
		double rhyd = HPDReactionRate(pw/*MPa*/,T,S,C_CH4,Sh,por); //mmol/L/a
		double Q/*J/L/a*/ = rhyd * ( 56.599 - 0.016744*( T ) );
		return Q; //J/a/L
	}

	double HydrateEquilibriumPressure( double T/*K*/, double S/*g/kg*/ )const{
		double p_eq/*MPa*/ = exp(- 1.6892692e3 - 0.15162984*T + 5.60912482e4/T + 2.72067506E2*log(T)
								 + S*( 2.06621298e4 + 14.18127*T - 8.3417066e-3*T*T - 3.78757519E5/T - 4.01544417e3*log(T) )
								 + S*S*( 6.04018045e2 + 0.420253434*T - 2.501548e-4*T*T - 1.09745335E4/T - 1.17640966e2*log(T) ) );
		return p_eq; //MPa
	}

	double HFDReactionRate( double pg/*MPa*/, double zf/*--*/, double T/*K*/, double S/*g/kg*/, double sw, double sh, double porosity )const{
		double p_eq/*MPa*/ = exp(- 1.6892692e3 - 0.15162984*T + 5.60912482e4/T + 2.72067506E2*log(T)
								 + S*( 2.06621298e4 + 14.18127*T - 8.3417066e-3*T*T - 3.78757519E5/T - 4.01544417e3*log(T) )
								 + S*S*( 6.04018045e2 + 0.420253434*T - 2.501548e-4*T*T - 1.09745335E4/T - 1.17640966e2*log(T) ) );
		double A_s/*cm^2/cm^3*/ = 1.e6/*A_0*/
								* 1.e-2/*conversion factor*/
								* pow( porosity*(1.-sh), 3./2. );
		double potential_P/*-*/ = p_eq/pg - 1.  ;

		double r = 0.;
		if(potential_P > 0. ) r/*mmol/L/a*/ = kr_hyddissociation /*(mmol/(cm^2.MPa.a))*(cm^3/L)*/
											* A_s 		 /*cm^2/cm^3*/
											* sh 		 /*--*/
											* pg	 	 /*MPa*/
											* potential_P/*--*/ ;
		else r/*mmol/L/a*/	= kr_hydformation 	/*(mmol/(cm^2.MPa.a))*(cm^3/L)*/
							* A_s 		 	/*cm^2/cm^3*/
							* sw*(1.-sw-sh) /*--*/
							* pg 			/*MPa*/
							* potential_P	/*--*/;

		return r; //mmol/L/a
	}

	double HeatOfHydrateDissociation( double pg/*MPa*/, double zf/*--*/, double T/*K*/, double S/*g/kg*/, double sw, double sh, double porosity )const{
		double rhyd = HFDReactionRate(pg/*MPa*/,zf,T,S,sw,sh,porosity); //mmol/L/a
		double Q/*J/L/a*/ = - rhyd * ( 56.599 - 0.016744*( T ) );
		return Q; //J/a/L
	}

	/**********************************************************************
	 * SEDIMENT AND HYDRAULIC PROPERTIES
	 * porosity
	 * permeability ( K_abs * KSF(Sh) )
	 * Swe
	 * Pc (Brooks Corey, PcSF(Sh) )
	 * krw
	 * krg
	 * tortuosity
	 */

	double SedimentPorosity
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal,
	 double time/*a*/, double dt/*a*/) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal); //cm

		double phi = phi_inf + (phi_0 - phi_inf) * std::exp (- exponent_beta/*1/cm*/ * (Zmax/*cm*/-x[dim-1]/*cm*/) );
		return phi; //-- [L_p/L_T]
	}

	// SEDIMENT PERMEABILITY (unit: cm^2)
	Dune::FieldVector<double,dim>
	SedimentPermeability
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal,
	 double time/*a*/, double dt/*a*/,
	 double Sh) const {
		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);//cm

		Dune::FieldVector<double,dim> K(permeability); //cm^2

		double SF = std::pow( (1.0 - Sh) , (5*exponent_m + 4)/(2*exponent_m) );
		K *= SF;

		return K; //cm^2
	}

	// EFFECTIVE WETTING SATURATION
	double EffectiveSw( const typename GV::Traits::template Codim<0>::Entity& element,
 	 	 	 			const Dune::FieldVector<double,dim>& xlocal,
						double time/*a*/, double dt/*a*/,
						double sw, double sh )const{
		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);//cm
		double se = std::max( ( sw - swr ) / (1. - sh - swr - sgr) , 0.0 ) ; //--
		return std::min( se , 1.0 );
	}

	// CAPILLARY PRESSURE
	double CapillaryPressure(const typename GV::Traits::template Codim<0>::Entity& element,
			 	 	 	 	 const Dune::FieldVector<double,dim>& xlocal,
							 double time/*a*/, double dt/*a*/,
							 double sw, double sh )const{
		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);//cm
		auto swe = EffectiveSw( element,xlocal,time,dt,sw,sh );
		double pc/*MPa*/ = std::min( pe_BC/*MPa*/ * std::pow(swe,-(1./lambda_BC)) , pcmax/*MPa*/ );

		double eta/*-*/ = (lambda_BC*exponent_m - 1.)/(lambda_BC*exponent_m);
		double SF/*-*/  = std::pow( (1.-sh) , -eta );

		return pc*SF; //MPa
	}

	// TORTUOSITY
	double Tortuosity(const typename GV::Traits::template Codim<0>::Entity& element,
			 	 	  const Dune::FieldVector<double,dim>& xlocal,
					  double time/*a*/, double dt/*a*/) const {

		double porosity = SedimentPorosity(element,xlocal,time,dt);
		double tau = 1. - 2.*log(porosity);
		return 1./tau;//--
	}

	// RELATIVE PERMEABILITIES
	double krw( const typename GV::Traits::template Codim<0>::Entity& element,
	   	 	 	const Dune::FieldVector<double,dim>& xlocal,
				double time/*a*/, double dt/*a*/,
				double sw, double sh ) const {

		double swe = EffectiveSw(element,xlocal,time,dt,sw,sh);

		double kr = std::pow(swe, (2.0/lambda_BC + 3.0) );
		return kr ;//--
	}

	double krg( const typename GV::Traits::template Codim<0>::Entity& element,
   	 	 		const Dune::FieldVector<double,dim>& xlocal,
				double time/*a*/, double dt/*a*/,
				double sw, double sh ) const {

		double swe = EffectiveSw(element,xlocal,time,dt,sw,sh);

		double kr = std::pow(1.0-swe, 2.0) * ( 1.0 - std::pow(swe, (2.0/lambda_BC + 1.0) ) );
		return kr;//--
	}

	/**********************************************************************/
	// BURIAL VELOCITY (vs) (unit: cm/annum)
	Dune::FieldVector<double,dim>
	AqueousBurialVelocity
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal,
	 double time/*a*/, double dt/*a*/) const {

		double porosity = SedimentPorosity(element,xlocal,time,dt);
		Dune::FieldVector<double,dim> v_s(0.);
		v_s[dim-1] = (-1.) * v_inf/*cm/a*/ * phi_inf;

		return v_s;//cm/a
	}

	Dune::FieldVector<double,dim>
	SolidsBurialVelocity
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal,
	 double time/*a*/, double dt/*a*/) const {

		double porosity = SedimentPorosity(element,xlocal,time,dt);
		Dune::FieldVector<double,dim> v_s(0.);
		v_s[dim-1] = (-1.) * v_inf/*cm/a*/ * (1.-phi_inf);

		return v_s; //cm/a
	}

	// GRAVITY VECTOR (units: (MPa.L)/(g.cm))
	Dune::FieldVector<double,dim>
	g( ) const {
		Dune::FieldVector<double,dim> gravity( 0. );
		double g = 0.;
		if(gravity_flag) g/*(MPa.L)/(g.cm)*/ = gravity_magnitude;
		gravity[ dim - 1 ] = g;
		gravity[0] = 0.;
		return gravity; //(MPa.L)/(g.cm)
	}

	/**********************************************************************/
	// DIFFUSION SOCEFFICIENTS FOR DISSOLVED SPECIES
	Dune::FieldVector<double,dim>
	AqueousDiffusionCoefficient(const typename GV::Traits::template Codim<0>::Entity& element,
								const Dune::FieldVector<double,dim>& xlocal,
								double time/*a*/, double dt/*a*/, int tag ) const {

		Dune::FieldVector<double,dim> D(0.); // cm^2/annum

		switch(tag) {
			case Tags::NH4_1p:
				D[0] = D_NH4_1p;
				break;
			case Tags::Cl_1n:
				D[0] = D_Cl_1n;
				break;
			case Tags::SO4_2n:
				D[0] = D_SO4_2n;
				break;
			case Tags::CH4:
				D[0] = D_CH4;
				break;
			case Tags::CO2:
				D[0] = D_CO2;
				break;
			case Tags::HCO3_1n:
				D[0] = D_HCO3_1n;
				break;
			case Tags::CO3_2n:
				D[0] = D_CO3_2n;
				break;
			default : D[0] = 0.; break;
		}
		for( int i=1; i<dim; i++ ) D[i] = D[0];

	    return D; // cm^2/annum
	}

	// DIFFUSION COEFFICIENTS FOR SOLID SPECIES
	Dune::FieldVector<double,dim>
	SolidsDiffusionCoefficient( const typename GV::Traits::template Codim<0>::Entity& element,
								const Dune::FieldVector<double,dim>& xlocal,
								double time/*a*/, double dt/*a*/, int tag ) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);
		double porosity = SedimentPorosity(element,xlocal,time,dt);

		Dune::FieldVector<double,dim> D(0.); // cm^2/annum
//		D[0] = porosity*porosity*D_OM/*cm^2/a*/;
//		D[1] = D[0];

		D[0] = Db0 * std::exp( - 0.5 * std::pow( (Zmax-x[1])/zb , 2.0 ) );
		D[1] = D[0];

		return D; // cm^2/annum

	}

	/**********************************************************************/
	double
	BottomWaterConcentration(const typename GV::Traits::template Codim<0>::Entity& element,
			 	 	 	 	 const Dune::FieldVector<double,dim>& xlocal,
							 double time/*a*/, double dt/*a*/, int tag ) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);//cm

		  double bwc = 0.;
			  switch(tag) {
			  case Tags::NH4_1p:
				  bwc = bwc_NH4_1p; //mmol/Lpw
				  break;
			  case Tags::Cl_1n:
				  bwc = bwc_Cl_1n; //mmol/Lpw
				  break;
			  case Tags::SO4_2n:
				  bwc = bwc_SO4_2n; //mmol/Lpw
				  break;
			  case Tags::CH4:
				  bwc = bwc_CH4; 	//mmol/Lpw
				  break;
			  case Tags::CO2:
				  bwc = bwc_CO2; 	//mmol/Lpw
				  break;
			  case Tags::OM:
				  bwc = bwc_OM; 	//mmol/Lds
				  break;
		  }

		  return bwc; //mmol/Lpw or mmol/Lds
	}

	/**********************************************************************/
	// KINETIC REACTION RATES (SOURCE TERMS)
	double KineticReactionRate(const typename GV::Traits::template Codim<0>::Entity& element,
	 	 	 	 	 	 	   const Dune::FieldVector<double,dim>& xlocal,
							   double time/*a*/, double dt/*a*/,
							   int tag, std::vector<double> C ) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);//cm

		double porosity = SedimentPorosity(element,xlocal,time,dt);

		double DIC/*mmol/Lpw*/ 	= C[Tags::CO2] + C[Tags::HCO3_1n] + C[Tags::CO3_2n];
		double w 				= v_inf * (1.-phi_inf)/(1.-porosity);
		double kx 				= kx0 * std::pow( ( a0/*a*/ + (Zmax/*cm*/-x[dim-1]/*cm*/)/w/*cm/a*/ ) , -a1 );
		double RPOC				= (1./ConversionFactors::CNP_a) * kx * C[Tags::OM] * K_C / (DIC + C[Tags::CH4] + K_C);

		// Monod constants
		double R_SO4 	= C[Tags::SO4_2n]/(K_SO42n+C[Tags::SO4_2n]);
		double R_CH4 	= K_SO42n/(K_SO42n+C[Tags::SO4_2n]);

		double R=0.;
		switch(tag) {
			case Tags::Rk_OMdegSO4 :
				R = 0.5 * RPOC * R_SO4 * (1.-porosity);
				break;
			case Tags::Rk_OMdegCH4 :
				R = 0.5 * RPOC * R_CH4 * (1.-porosity);
				break;
			case Tags::Rk_AOM :
				R = k_AOM /*L/(mmol.a*/ * C[Tags::CH4] * C[Tags::SO4_2n] * porosity;
				break;
		}
		return R; //mmol/L_T/a
	}

	// ACTIVITY
	double Log10Activity( int tag, std::vector<double> C ) const {

		double activity = std::max( C[tag], 1.e-20 ) * ConversionFactors::F_Mkg ; //C*F_Mkg -> molal conc. [mol/kg_H2O]
		double logA = std::log10(activity);
		return logA;
	}

	// EQUILIBRIUM CONSTANTS FOR THE MASS ACTION LAWS
	double EquilibriumConstant( int tag ) const {

		double K = 0.;

		switch(tag) {
		    case Tags::Re_AB1:
		    	K = Keq_AB1;
		        break;
		    case Tags::Re_AB2:
		    	K = Keq_AB2;
		        break;
		}
	    return K;
	}

	double Log10EquilibriumConstant( int tag ) const {

		double logK = std::log10(EquilibriumConstant(tag));
	    return logK;
	}
	/**********************************************************************/

};

template<typename GV,typename PTree>
class PropertiesAndParameters
{
private:
	  const GV& gv;
	  const PTree& ptree;
	  double *time;
	  double *dt;

	  const static int dim = GV::dimension;

public:

	//PARAMETERS AND PROPERTIES
  	Tags tag;
  	Matrices matrix;
  	ConversionFactors factor;
  	Domain<GV,PTree> domain;
  	Properties<GV,PTree> property;

  	//! construct from grid view
  	PropertiesAndParameters ( const GV& gv_ ,
			 	 	 	 	  const PTree& ptree_,
							  double *time_,
							  double *dt_  )
	: gv( gv_ ) ,
	  ptree(ptree_),
	  time(time_),
	  dt(dt_),
	  domain(gv_,ptree_),
	  property(gv_,ptree_,time_,dt_)
  	{}

	/******************************************************************************/

  	template <typename T>
  	std::string to_string_with_precision_scientific(const T a_value, const int n = 6) const {
  	    std::ostringstream out;
  	    out.precision(n);
  	    out << std::scientific << a_value;
  	    return out.str();
  	}

  	template <typename T>
  	std::string to_string_with_precision_fixed(const T a_value, const int n = 6) const {
  	    std::ostringstream out;
  	    out.precision(n);
  	    out << std::fixed << a_value;
  	    return out.str();
  	}

	/******************************************************************************/

  	void ReportStatistics( std::string file_name,
  						   double time,
						   double dt,
						   int total_newton_iterations,
						   double clock_time_elapsed ) {

  		std::fstream result;

  		if(time == 0. ){
  			result.open(file_name, std::fstream::out | std::fstream::trunc);
  			result	<< "time [s]" << '\t'
  					<< "dt [s]"	<< '\t'
					<< "total no. of newton iterations" << '\t'
					<< "clock time [s]"
  					<< std::endl;
  			result.close();
  		}

  		result.open(file_name, std::fstream::app);
  		double t_new = time+dt;

		result	<< time	<< '\t'
				<< dt	<< '\t'
				<< total_newton_iterations << '\t'
				<< clock_time_elapsed
				<< std::endl;
		result.close();
  	}

	/******************************************************************************/

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};

#endif /* PROBLEM_PARAMETERS_HH_ */
