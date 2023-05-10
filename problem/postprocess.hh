#ifndef PROBLEM_POSTPROCESS_HH_
#define PROBLEM_POSTPROCESS_HH_


template< class GV, class PARAMS, class GFS, class U >
class PostProcess{

public:
	typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
    typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
	typedef typename U::template LocalView<LFSCache> VectorView;

	typedef typename GV::Traits::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::IndexSet IndexSet;
    static const int dim = GV::dimension;

    PostProcess(std::string file_name_,
    			const GV& gv_,
				const PARAMS& param_,
				GFS 	gfs_,
				U		*u_,
				std::vector<std::vector<double>> *Sk_,
				std::vector<std::vector<double>> *Ses_inv_,
				std::vector<std::vector<double>> *Sep_star_,
				std::vector<std::vector<double>> *matU_,
				double *time_	,
				double *dt_ )
    :file_name(file_name_),
	 gv(gv_),
	 param( param_ ),
	 gfs(gfs_),
	 u(u_),
	 Sk(Sk_),
	 Ses_inv(Ses_inv_),
	 Sep_star(Sep_star_),
	 matU(matU_),
	 time( time_ ),
	 dt( dt_ )
    {}

	virtual ~PostProcess(){}

	void print_all_qoi_per_timeoutput(){

		std::string name = file_name;
		name += "_qoi_";
		auto t_new = (*time) + (*dt);
		t_new *= (1./1000.);
//		name += to_string_with_precision_fixed( t_new , 2 ); //TODO: Check why this fnc isn't working
  	    std::ostringstream t;
  	    t.precision(0);
  	    t << std::fixed << t_new;
  	    name += t.str();
		name += ".txt";
		std::fstream result;
		result.open(name, std::fstream::out | std::fstream::trunc );
		if (!result.is_open()){
			std::cout << "Could not open file:" << name << std::endl;
			exit(0);
		}

		result	<< "z[m]" << " ";
		result << "pw[MPa]" << " ";
		result << "pg[MPa]" << " ";
		result << "pc[MPa]" << " ";
		result << "pe[MPa]" << " ";
		result << "sh[-]" << " ";
		result << "sg[-]" << " ";
		result << "sw[-]" << " ";
		result << "T[K]" << " ";
		result << "K[m^2]" << " ";
		result << "salinity[g/kg]" << " ";
		result << "zf[-]" << " ";
		result << "Ceq_CH4[mmol/Lpw]" << " ";
		result << "Ceq_HYD[mmol/Lpw]" << " ";
		for(int i=0; i<param.tag.Ns; i++){
			result <<  param.tag.getSpeciesName(i) << " ";
		}
		result << std::endl;
		result.flush();

		LFS lfs(gfs);
		LFSCache lfs_cache(lfs);
		VectorView u_view( *u );

		// Iterate over each element
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			lfs.bind(*self);
			lfs_cache.update();
			u_view.bind(lfs_cache);
		    std::vector<double> ul(lfs.size());
		    u_view.read( ul );

			// Get Element geometry
			typedef typename LeafIterator::Entity::Geometry ElementGeometry;
			const auto& geo = (*self).geometry();
	        auto cell_center_local 	= geo.local(geo.center());
	        auto cell_center_global = geo.center();
//
			// compute PVs at local center
	        std::vector<double> C(param.tag.Ns,0.);
	        // assign primary C values
	        for(int i=0; i<param.tag.nPrC; i++ ){
	        	C[i] = ul[lfs.child(i).localIndex(0)];
	        }

	        if( param.tag.Ne>0 ){
				// evaluate secondary C values as functions of primary C's using the equilibrium constraints
				// log C_s = Ses_inv * log K + Sep_star * log C_p
				std::vector<double> Log10K(param.tag.Ne,0.);
				for(int i=0; i<param.tag.Ne; i++){
					Log10K[i] = param.property.Log10EquilibriumConstant(i);
				}
				auto Log10K_star  = operation.multiplyMatrixVector( (*Ses_inv ), Log10K  );

				std::vector<double> Log10Cp(param.tag.nPDE,0.);
				for(int i=0; i<param.tag.nPDE; i++){
					Log10Cp[i] = param.property.Log10Activity(i,C);
				}
				auto Log10Cp_star = operation.multiplyMatrixVector( (*Sep_star), Log10Cp );

				auto Fmkg = param.factor.F_Mkg;
				for(int i=param.tag.nIQC; i<param.tag.Ne; i++ ){
					C[param.tag.nPDE+i] = std::pow(10.,Log10K_star[i] + Log10Cp_star[i] )/Fmkg;
				}
	        }

	        double pw = ul[lfs.child(Tags::Pw).localIndex(0)];//in MPa
	        double sh = ul[lfs.child(Tags::Sh).localIndex(0)];
	        double sg = ul[lfs.child(Tags::Sg).localIndex(0)];
	        double sw = 1.-sg-sh;
	        double T  = ul[lfs.child(Tags::T ).localIndex(0)];// in K
	        double pc = param.property.CapillaryPressure((*self),cell_center_local,(*time),(*dt),sw,sh);
	        double pg = pw + pc;
	        double S = param.property.Salinity(C[param.tag.Cl_1n]); // in g/kg
	        double pe = param.property.HydrateEquilibriumPressure(T,S);
	        double K = param.property.SedimentPermeability((*self),cell_center_local,(*time),(*dt),sh)[0];
	        double Ceq_HYD = param.property.EquilibriumConcentrationHydrate(pw/*MPa*/,T/*K*/,S/*g/kg*/);
	        double zf = param.property.CompressibilityFactorCH4(pg,T);
	        double Ceq_CH4 = param.property.EquilibriumConcentrationCH4((*self),cell_center_local,pg/*MPa*/,zf/*--*/,T/*K*/,S/*g/kg*/);

			result	<< (param.domain.SeaFloorHeight(cell_center_global[0])-cell_center_global[dim-1])/100.0 << " ";
			result << pw << " ";
			result << pg << " ";
			result << pc << " ";
			result << pe << " ";
			result << sh << " ";
			result << sg << " ";
			result << sw << " ";
			result << T << " ";
			result << K << " ";
			result << S << " ";
			result << zf << " ";
			result << Ceq_CH4 << " ";
			result << Ceq_HYD << " ";
			for(int i=0; i<param.tag.Ns; i++){
				result << C[i] << " ";
			}
			result << std::endl;
			result.flush();

		}//END: loop over volumes

		result.close();

	}//END: print_all_qoi_per_timeoutput()


	void print_qoi( int tag ){

		auto t_new = (*time) + (*dt);

		std::fstream result;
		std::string name = file_name;
		name +="_";
		for(int i=0; i<param.tag.Ns; i++){
			if(tag==i) name += param.tag.getSpeciesName(i);
		}
		if(tag==Tags::PP_pw ) 	name += "pw";
		if(tag==Tags::PP_pg ) 	name += "pg";
		if(tag==Tags::PP_pc ) 	name += "pc";
		if(tag==Tags::PP_pe ) 	name += "pe";
		if(tag==Tags::PP_sh ) 	name += "sh";
		if(tag==Tags::PP_sg ) 	name += "sg";
		if(tag==Tags::PP_sw ) 	name += "sw";
		if(tag==Tags::PP_T  ) 	name += "T";
		if(tag==Tags::PP_K  ) 	name += "K";
		if(tag==Tags::PP_sal) 	name += "S";
		if(tag==Tags::PP_zf ) 	name += "zf";
		if(tag==Tags::PP_CeqCH4)name += "Ceq_CH4";
		if(tag==Tags::PP_CeqHYD)name += "Ceq_HYD";
		name += "_profile";
		name += ".txt";
		if( (*time)<1.e-6 ){
			result.open(name, std::fstream::out | std::fstream::trunc );
			if (!result.is_open()){
				std::cout << "Could not open file:" << name << std::endl;
				exit(0);
			}
			result	<< "time[yr]" << " ";
		}else{
			result.open(name, std::fstream::out | std::fstream::app );
			result << t_new << " ";
		}

		LFS lfs(gfs);
		LFSCache lfs_cache(lfs);
		VectorView u_view( *u );

		// Iterate over each element
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			lfs.bind(*self);
			lfs_cache.update();
			u_view.bind(lfs_cache);
		    std::vector<double> ul(lfs.size());
		    u_view.read( ul );

			// Get Element geometry
			typedef typename LeafIterator::Entity::Geometry ElementGeometry;
			const auto& geo = (*self).geometry();
	        auto cell_center_local 	= geo.local(geo.center());
	        auto cell_center_global = geo.center();

			// compute PVs at local center
	        std::vector<double> C(param.tag.Ns,0.);
	        // assign primary C values
	        for(int i=0; i<param.tag.nPrC; i++ ){
	        	C[i] = ul[lfs.child(i).localIndex(0)];
	        }

	        if( param.tag.Ne>0 ){
				// evaluate secondary C values as functions of primary C's using the equilibrium constraints
				// log C_s = Ses_inv * log K + Sep_star * log C_p
				std::vector<double> Log10K(param.tag.Ne,0.);
				for(int i=0; i<param.tag.Ne; i++){
					Log10K[i] = param.property.Log10EquilibriumConstant(i);
				}
				auto Log10K_star  = operation.multiplyMatrixVector( (*Ses_inv ), Log10K  );

				std::vector<double> Log10Cp(param.tag.nPDE,0.);
				for(int i=0; i<param.tag.nPDE; i++){
					Log10Cp[i] = param.property.Log10Activity(i,C);
				}
				auto Log10Cp_star = operation.multiplyMatrixVector( (*Sep_star), Log10Cp );

				auto Fmkg = param.factor.F_Mkg;
				for(int i=param.tag.nIQC; i<param.tag.Ne; i++ ){
					C[param.tag.nPDE+i] = std::pow(10.,Log10K_star[i] + Log10Cp_star[i] )/Fmkg;
				}
	        }

	        double pw = ul[lfs.child(Tags::Pw).localIndex(0)];//in MPa
	        double sh = ul[lfs.child(Tags::Sh).localIndex(0)];
	        double sg = ul[lfs.child(Tags::Sg).localIndex(0)];
	        double sw = 1.-sg-sh;
	        double T  = ul[lfs.child(Tags::T ).localIndex(0)];// in K
	        double pc = param.property.CapillaryPressure((*self),cell_center_local,(*time),(*dt),sw,sh);
	        double pg = pw + pc;
	        double S = param.property.Salinity(C[param.tag.Cl_1n]); // in g/kg
	        double pe = param.property.HydrateEquilibriumPressure(T,S);
	        double K = param.property.SedimentPermeability((*self),cell_center_local,(*time),(*dt),sh)[0];
	        double Ceq_HYD = param.property.EquilibriumConcentrationHydrate(pw/*MPa*/,T/*K*/,S/*g/kg*/);
	        double zf = param.property.CompressibilityFactorCH4(pg,T);
	        double Ceq_CH4 = param.property.EquilibriumConcentrationCH4((*self),cell_center_local,pg/*MPa*/,zf/*--*/,T/*K*/,S/*g/kg*/);

			if( (*time)<1.e-6 ){
				result	<< (param.domain.SeaFloorHeight(cell_center_global[0])-cell_center_global[dim-1])/100.0 << " ";
			}else{
				for(int i=0; i<param.tag.Ns; i++){
					if(tag==i) result << C[i] << " ";
				}
				if(tag==Tags::PP_pw ) result << pw << " ";
				if(tag==Tags::PP_pg ) result << pg << " ";
				if(tag==Tags::PP_pc ) result << pc << " ";
				if(tag==Tags::PP_pe ) result << pe << " ";
				if(tag==Tags::PP_sh ) result << sh << " ";
				if(tag==Tags::PP_sg ) result << sg << " ";
				if(tag==Tags::PP_sw ) result << sw << " ";
				if(tag==Tags::PP_T  ) result << T  << " ";
				if(tag==Tags::PP_K  ) result << K  << " ";
				if(tag==Tags::PP_sal) result << S  << " ";
				if(tag==Tags::PP_zf ) result << zf << " ";
				if(tag==Tags::PP_CeqCH4) result << Ceq_CH4 << " ";
				if(tag==Tags::PP_CeqHYD) result << Ceq_HYD << " ";
			}

		}//END: loop over volumes

		result << std::endl;
		result.flush();
		result.close();

	} //END: print_qoi

private:
	std::string file_name;
    const GV& gv;
	const PARAMS& param;
	GFS gfs;
	U *u;
	std::vector<std::vector<double>> *matU;
	std::vector<std::vector<double>> *Ses_inv;
	std::vector<std::vector<double>> *Sep_star;
	std::vector<std::vector<double>> *Sk;
	double *time;
	double *dt;
	Operations operation;
};


#endif /* PROBLEM_POSTPROCESS_HH_ */
