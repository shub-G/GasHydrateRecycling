/*
 * localoperator.hh
 *
 *  Created on: Aug 30, 2021
 *      Author: sgupta
 */

#ifndef OPERATORS_LOCALOPERATOR_HH_
#define OPERATORS_LOCALOPERATOR_HH_

template <	class GV, class Params, class BC >
class LocalOperator :
  public Dune::PDELab::NumericalJacobianApplyVolume		< LocalOperator<GV,Params,BC> >,
  public Dune::PDELab::NumericalJacobianVolume			< LocalOperator<GV,Params,BC> >,
  public Dune::PDELab::NumericalJacobianApplySkeleton	< LocalOperator<GV,Params,BC> >,
  public Dune::PDELab::NumericalJacobianSkeleton		< LocalOperator<GV,Params,BC> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary	< LocalOperator<GV,Params,BC> >,
  public Dune::PDELab::NumericalJacobianBoundary		< LocalOperator<GV,Params,BC> >,
  public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
    const GV& gv;
	const Params&	  param;
	const BC&	 	  bc;
	std::vector<std::vector<double>> *matU;
	std::vector<std::vector<double>> *Ses_inv;
	std::vector<std::vector<double>> *Sep_star;
	std::vector<std::vector<double>> *Sk;
	double 			 *time;
	double 			 *dt;
	constexpr static double eps 	= 1.0e-6;
	constexpr static double pi 		= atan(1.)*4 ;
	Operations operation;

public:
	  // pattern assembly flags
	  enum { doPatternVolume	= true };
	  enum { doPatternSkeleton	= true };

	  // residual assembly flags
	  enum { doAlphaVolume  	= true };
	  enum { doAlphaSkeleton	= true };
	  enum { doAlphaBoundary	= true };

	  typedef typename GV::IndexSet IndexSet;

	  // constructor stores parameters
	  LocalOperator(  const GV& 	 gv_,
			  	  	  const Params&	 param_,
					  const BC& 	 bc_	,
					  std::vector<std::vector<double>> *Sk_,
					  std::vector<std::vector<double>> *Ses_inv_,
					  std::vector<std::vector<double>> *Sep_star_,
					  std::vector<std::vector<double>> *matU_,
					  double *time_	,
					  double *dt_	)
	  :  gv(gv_),
		 param( param_ ),
		 bc( bc_ ),
		 Sk(Sk_),
		 Ses_inv(Ses_inv_),
		 Sep_star(Sep_star_),
		 matU(matU_),
		 time( time_ ),
		 dt( dt_ )
	  {}

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
			// Reference to cell
	        const auto& cell = eg.entity();
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);

	        // get geometry
	        auto geo = eg.geometry();

			// dimension
//			const auto dim = geo.mydimension;
	        const auto dim = gv.dimension;

	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0);
	        auto cell_center_global = cell.geometry().global(cell_center_local);
	        auto cell_volume = geo.volume();


	        // compute phase variables at local center
	        auto pw = x(lfsu.child(param.tag.Pw),0);
	        auto sh = x(lfsu.child(param.tag.Sh),0);
	        auto sg = x(lfsu.child(param.tag.Sg),0);
	        auto T  = x(lfsu.child(param.tag.T ),0);

	        // secondary phase variables
	        auto sw = 1.-sg-sh;
	        auto pc = param.property.CapillaryPressure(cell,cell_center_local,(*time),(*dt),sw,sh);
	        auto pg = pw + pc;
	        auto peff = (sw*pw+sg*pg)/(1.-sh);
	        auto zg = param.property.CompressibilityFactorCH4(pg,T);
	        auto por = param.property.SedimentPorosity(cell,cell_center_local,(*time),(*dt));

			// compute species concentrations at local center
	        std::vector<double> C(param.tag.Ns,0.);
	        // assign primary C values
	        for(int i=0; i<param.tag.nPrC; i++ ){
	        	C[i] = x(lfsu.child(i),0);
	        }

	        if(param.tag.Ne>0){ // only when equilibrium reactions are considered
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

				std::vector<double> Log10C(param.tag.Ns,0.);
				for(int i=0; i<param.tag.Ns; i++){
					Log10C[i] = param.property.Log10Activity(i,C);
				}

				auto Se = param.matrix.S_eqb();
				auto SeLogC = operation.multiplyMatrixVector(Se,Log10C);
				auto constraint = Log10K;
				for( int i=0; i<constraint.size(); i++ ){
					constraint[i] -= SeLogC[i];
				}

				// RESIDUALS: INEQUALITY CONSTRAINTS FOR MINERAL PRECIPITATION
				double alpha = 1.;
				for( int i=0; i<param.tag.nIQC; i++ ){
					if( C[param.tag.nPDE+i] - alpha * constraint[i] > 0. ){
						r.accumulate(lfsu.child(param.tag.nPDE+i) , 0., alpha * constraint[i]*cell_volume );
					}else{
						r.accumulate(lfsu.child(param.tag.nPDE+i) , 0., C[param.tag.nPDE+i]*cell_volume );
					}
				}
	        }

	        //evaluate salinity
	        double salinity = param.property.Salinity(C[param.tag.Cl_1n]);

	        // Kinetic reaction source terms
	        std::vector<double> Rk(param.tag.Nk,0.);
	        for(int i=0; i<param.tag.Nk; i++){
	        	Rk[i] = param.property.KineticReactionRate(cell,cell_center_local,(*time),(*dt),i,C);
	        }
	        auto Sk_t = operation.transposeMatrix( (*Sk) );
	        auto Sk_t_corr = Sk_t;
	        auto SkRk = operation.multiplyMatrixVector( Sk_t_corr,Rk );
	        auto USkRk = operation.multiplyMatrixVector( (*matU),SkRk );

	        // RESIDUALS: SPECIES MASS BALANCE
	        for( int i=0; i<param.tag.nPDE; i++ ){
		        	r.accumulate(lfsu.child(i) , 0., -USkRk[i]*cell_volume );
		        	if( i==param.tag.CH4 ){
		        		auto source_g = param.property.HFDReactionRate(pg,zg,T,salinity,sw,sh,por);
		        		r.accumulate(lfsu.child(i) , 0., -source_g*cell_volume );
		        	}
	        }

	        // REDIDUAL: MASS BALANCE HYDRATE PHASE
	        auto source_HFD_h = (-1.) * param.property.MolarMassHYD() * param.property.HFDReactionRate(pg,zg,T,salinity,sw,sh,por);
	        auto source_HPD_h = (+1.) * param.property.MolarMassHYD() * param.property.HPDReactionRate(pw,T,salinity,C[param.tag.CH4],sh,por);
	        r.accumulate(lfsu.child(param.tag.Sh) , 0., -(source_HFD_h+source_HPD_h)*cell_volume );

	        // RESIDUAL: DISSOLVED METHANE --> FREE GAS TRANSITION
	        auto C_ch4_eq = param.property.EquilibriumConcentrationCH4(cell,cell_center_local,pg,zg,T,salinity);
	        double alpha = 1.0;
			if( sg - alpha * (C_ch4_eq - C[param.tag.CH4]) > 0. ){
//				std::cout<< cell_center_global[1]/100.0 << '\t' << C_ch4_eq-C[param.tag.CH4] << std::endl;
				r.accumulate(lfsu.child(param.tag.Sg) , 0., alpha * (C_ch4_eq - C[param.tag.CH4])*cell_volume );
			}else{
				r.accumulate(lfsu.child(param.tag.Sg) , 0., sg*cell_volume );
			}

			// RESIDUAL: ENERGY BALANCE
	        auto source_heat = param.property.HeatOfHydrateDissociation(pg,zg,T,salinity,sw,sh,por)
	        				 + param.property.HeatOfHydrateDissolution(pw,T,salinity,C[param.tag.CH4],sh,por);
	        r.accumulate(lfsu.child(param.tag.T) , 0., -source_heat*cell_volume );

	  }

	  // skeleton integral depending on test and ansatz functions
	  // each face is only visited ONCE!
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_skeleton (const IG& ig,
			  	  	  	   const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
						   const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
						   R& r_s, R& r_n) const
	  {
	        // References to inside and outside cells
	        const auto& cell_inside  = ig.inside();
	        const auto& cell_outside = ig.outside();
			const IndexSet &indexSet = gv.indexSet();
			int inside_cell_number  = indexSet.index(cell_inside);
			int outside_cell_number = indexSet.index(cell_outside);

	        // get geometries
	        auto geo = ig.geometry();
//	        const auto dim = geo.mydimension;
	        const auto dim = gv.dimension;
	        auto geo_inside  = cell_inside.geometry();
	        auto geo_outside = cell_outside.geometry();

	        // cell geometries
	        auto ref_el_inside 	= referenceElement(geo_inside);
	        auto ref_el_outside = referenceElement(geo_outside);
	        auto inside_cell_center_local 	= ref_el_inside.position(0,0);
	        auto outside_cell_center_local 	= ref_el_outside.position(0,0);
	        auto inside_cell_center_global 	= geo_inside.center();
	        auto outside_cell_center_global = geo_outside.center();

	        // distance of cell centers
	        auto d = outside_cell_center_global;
	        d -= inside_cell_center_global;
	        auto distance = d.two_norm();

	        // face geometry
	        auto ref_el = referenceElement(geo);
	        auto face_center_local = ref_el.position(0,0);
	        auto face_center_global = geo.center();
	        auto face_volume = geo.volume();
	        auto normal = ig.unitOuterNormal(face_center_local);

	        // compute phase variables at local self and neighbour centers
	        auto pw_s = x_s(lfsu_s.child(param.tag.Pw),0);
	        auto sh_s = x_s(lfsu_s.child(param.tag.Sh),0);
	        auto sg_s = x_s(lfsu_s.child(param.tag.Sg),0);
	        auto T_s  = x_s(lfsu_s.child(param.tag.T ),0);

	        auto pw_n = x_n(lfsu_n.child(param.tag.Pw),0);
	        auto sh_n = x_n(lfsu_n.child(param.tag.Sh),0);
	        auto sg_n = x_n(lfsu_n.child(param.tag.Sg),0);
	        auto T_n  = x_n(lfsu_n.child(param.tag.T ),0);

			// compute species concentrations at local self and neighbour centers
	        std::vector<double> C_s(param.tag.Ns,0.);
	        for(int i=0; i<param.tag.nPrC; i++ ){
	        	C_s[i] = x_s(lfsu_s.child(i),0);
	        }
	        std::vector<double> C_n(param.tag.Ns,0.);
	        for(int i=0; i<param.tag.nPrC; i++ ){
	        	C_n[i] = x_n(lfsu_n.child(i),0);
	        }

	        if( param.tag.Ne>0 ){ // only if the equilibrium reactions are considered
				// evaluate secondary C values as functions of primary C's using the equilibrium constraints
				// log C_s = Ses_inv * log K + Sep_star * log C_p
				std::vector<double> Log10K(param.tag.Ne,0.);
				for(int i=0; i<param.tag.Ne; i++){
					Log10K[i] = param.property.Log10EquilibriumConstant(i);
				}
				auto Log10K_star  = operation.multiplyMatrixVector( (*Ses_inv ), Log10K  );

				std::vector<double> Log10Cp_s(param.tag.nPDE,0.);
				std::vector<double> Log10Cp_n(param.tag.nPDE,0.);
				for(int i=0; i<param.tag.nPDE; i++){
					Log10Cp_s[i] = param.property.Log10Activity(i,C_s);
					Log10Cp_n[i] = param.property.Log10Activity(i,C_n);
				}
				auto Log10Cp_star_s = operation.multiplyMatrixVector( (*Sep_star), Log10Cp_s );
				auto Log10Cp_star_n = operation.multiplyMatrixVector( (*Sep_star), Log10Cp_n );

				auto Fmkg = param.factor.F_Mkg;
				for(int i=param.tag.nIQC; i<param.tag.Ne; i++ ){
					C_s[param.tag.nPDE+i] = std::pow(10.,Log10K_star[i] + Log10Cp_star_s[i])/Fmkg;
					C_n[param.tag.nPDE+i] = std::pow(10.,Log10K_star[i] + Log10Cp_star_n[i])/Fmkg;
				}
	        }

	        // salinity
	        auto salinity_s = param.property.Salinity(C_s[param.tag.Cl_1n]);
	        auto salinity_n = param.property.Salinity(C_n[param.tag.Cl_1n]);

	        // secondary phase variables and properties
	        auto sw_s = 1.-sg_s-sh_s;
	        auto pc_s = param.property.CapillaryPressure(cell_inside,inside_cell_center_local,(*time),(*dt),sw_s,sh_s);
	        auto pg_s = pw_s + pc_s;
	        auto peff_s = (sw_s*pw_s+sg_s*pg_s)/(1.-sh_s);
	        auto por_s = param.property.SedimentPorosity(cell_inside,inside_cell_center_local,(*time),(*dt));
	        auto krg_s = param.property.krg(cell_inside,inside_cell_center_local,(*time),(*dt),sw_s,sh_s);
	        auto krw_s = param.property.krw(cell_inside,inside_cell_center_local,(*time),(*dt),sw_s,sh_s);;
	        auto K_s = std::abs(param.property.SedimentPermeability(cell_inside,inside_cell_center_local,(*time),(*dt),sh_s)
	        		 * ig.unitOuterNormal(face_center_local));
	        auto zg_s = param.property.CompressibilityFactorCH4(pg_s,T_s);
	        auto mug_s = param.property.GasPhaseViscosity(pg_s,zg_s,T_s);
	        auto muw_s = param.property.WaterPhaseViscosity(pw_s,T_s,salinity_s);
	        auto rhog_s = param.property.GasPhaseDensity(pg_s,zg_s,T_s);
	        auto rhow_s = param.property.WaterPhaseDensity(pw_s,T_s,salinity_s);
	        auto rhoh_s = param.property.HydratePhaseDensity(peff_s,T_s);
	        auto rhos_s = param.property.SedimentDensity(peff_s,T_s);
	        auto cpg_s = param.property.GasPhaseCp(pg_s,zg_s,T_s);
	        auto cpw_s = param.property.WaterPhaseCp(pw_s,T_s,salinity_s);
	        auto ktg_s = param.property.GasPhaseThermalConductivity(pg_s,zg_s,T_s);
	        auto ktw_s = param.property.WaterPhaseThermalConductivity(pw_s,T_s,salinity_s);
	        auto kth_s = param.property.HydatePhaseThermalConductivity(peff_s,T_s);
	        auto kts_s = param.property.SedimentThermalConductivity(peff_s,T_s);
	        auto kteff_s = (1.-por_s) * kts_s
	        			 + por_s * sw_s * ktw_s
						 + por_s * sg_s * ktg_s
						 + por_s * sh_s * kth_s;

	        auto sw_n = 1.-sg_n-sh_n;
	        auto pc_n = param.property.CapillaryPressure(cell_outside,outside_cell_center_local,(*time),(*dt),sw_n,sh_n);
	        auto pg_n = pw_n + pc_n;
	        auto peff_n = (sw_n*pw_n+sg_n*pg_n)/(1.-sh_n);
	        auto por_n = param.property.SedimentPorosity(cell_outside,outside_cell_center_local,(*time),(*dt));
	        auto krg_n = param.property.krg(cell_outside,outside_cell_center_local,(*time),(*dt),sw_n,sh_n);
	        auto krw_n = param.property.krw(cell_outside,outside_cell_center_local,(*time),(*dt),sw_n,sh_n);;
	        auto K_n = std::abs(param.property.SedimentPermeability(cell_outside,outside_cell_center_local,(*time),(*dt),sh_n)
	        		 * ig.unitOuterNormal(face_center_local));
	        auto zg_n = param.property.CompressibilityFactorCH4(pg_n,T_n);
	        auto mug_n = param.property.GasPhaseViscosity(pg_n,zg_n,T_n);
	        auto muw_n = param.property.WaterPhaseViscosity(pw_n,T_n,salinity_n);
	        auto rhog_n = param.property.GasPhaseDensity(pg_n,zg_n,T_n);
	        auto rhow_n = param.property.WaterPhaseDensity(pw_n,T_n,salinity_n);
	        auto rhoh_n = param.property.HydratePhaseDensity(peff_n,T_n);
	        auto rhos_n = param.property.SedimentDensity(peff_n,T_n);
	        auto cpg_n = param.property.GasPhaseCp(pg_n,zg_n,T_n);
	        auto cpw_n = param.property.WaterPhaseCp(pw_n,T_n,salinity_n);
	        auto ktg_n = param.property.GasPhaseThermalConductivity(pg_n,zg_n,T_n);
	        auto ktw_n = param.property.WaterPhaseThermalConductivity(pw_n,T_n,salinity_n);
	        auto kth_n = param.property.HydatePhaseThermalConductivity(peff_n,T_n);
	        auto kts_n = param.property.SedimentThermalConductivity(peff_n,T_n);
	        auto kteff_n = (1.-por_n) * kts_n
	        			 + por_n * sw_n * ktw_n
						 + por_n * sg_n * ktg_n
						 + por_n * sh_n * kth_n;

	        // properties at interface
	        auto K_int  = 2.*K_s*K_n;
	        if(K_s+K_n>0) K_int *= 1./( K_s+K_n );
	        auto kteff_int = 2.*kteff_s*kteff_n;
	        if(kteff_s+kteff_n>0.) kteff_int *= 1./(kteff_s+kteff_n);
	        auto rhog_int = 0.5*(rhog_s + rhog_n);
	        auto rhow_int = 0.5*(rhow_s + rhow_n);

	        // phase velocities
			auto normalvelocity_sed_s = param.property.SolidsBurialVelocity(cell_inside,inside_cell_center_local,(*time),(*dt))
					* ig.unitOuterNormal(face_center_local) ;
			auto normalvelocity_sed_aq = param.property.AqueousBurialVelocity(cell_inside,inside_cell_center_local,(*time),(*dt))
					* ig.unitOuterNormal(face_center_local) ;
			double omegaup_sed_s_s = 0., omegaup_sed_s_n = 0.;
			double omegaup_sed_aq_s = 0., omegaup_sed_aq_n = 0.;
			if( (-1.)*normalvelocity_sed_s>0.){
				omegaup_sed_s_s = 0.;
				omegaup_sed_s_n = 1.;
			}else{
				omegaup_sed_s_s = 1.;
				omegaup_sed_s_n = 0.;
			}
			if( (-1.)*normalvelocity_sed_aq>0.){
				omegaup_sed_aq_s = 0.;
				omegaup_sed_aq_n = 1.;
			}else{
				omegaup_sed_aq_s = 1.;
				omegaup_sed_aq_n = 0.;
			}
			auto gravity = param.property.g() * ig.unitOuterNormal(face_center_local) ;
			//GAS PHASE
			auto normalpotential_g = (pg_n - pg_s)/distance
								   + rhog_int * gravity ;
			auto normalvelocity_g_tmp = - K_int
										* (krg_s/mug_s)
										* normalpotential_g;
			double omegaup_g_s = 0., omegaup_g_n = 0.;
//			if( normalvelocity_g_tmp+sg_s*normalvelocity_sed_aq<0.){
			if( normalpotential_g+normalvelocity_sed_aq>0.){
//			if( normalpotential_g>0.){
				omegaup_g_s = 0.;
				omegaup_g_n = 1.;
			}else{
				omegaup_g_s = 1.;
				omegaup_g_n = 0.;
			}
			auto normalvelocity_g = - K_int
									* ( omegaup_g_s * (krg_s/mug_s) + omegaup_g_n * (krg_n/mug_n) )
									* normalpotential_g;
			double normalvelocity_g_total = normalvelocity_g + (omegaup_sed_aq_s*sg_s+omegaup_sed_aq_n*sg_n) * normalvelocity_sed_aq;
//			double normalvelocity_g_total = normalvelocity_g + (omegaup_g_s*sg_s+omegaup_g_n*sg_n) * normalvelocity_sed_aq;
			//WATER PHASE
			auto normalpotential_w = (pw_n - pw_s)/distance
								   + rhow_int * gravity ;
			auto normalvelocity_w_tmp = - K_int
										* (krw_s/muw_s)
										* normalpotential_w;
			double omegaup_w_s = 0., omegaup_w_n = 0.;
//			if( normalvelocity_w_tmp+sw_s*normalvelocity_sed_aq<0.){
			if( normalpotential_w+normalvelocity_sed_aq>0.){
//			if( normalpotential_w>0.){
				omegaup_w_s = 0.;
				omegaup_w_n = 1.;
			}else{
				omegaup_w_s = 1.;
				omegaup_w_n = 0.;
			}
			auto normalvelocity_w = - K_int
									* ( omegaup_w_s * (krw_s/muw_s) + omegaup_w_n * (krw_n/muw_n) )
									* normalpotential_w;
			double normalvelocity_w_total = normalvelocity_w + (omegaup_sed_aq_s*sw_s+omegaup_sed_aq_n*sw_n) * normalvelocity_sed_aq;
//			double normalvelocity_w_total = normalvelocity_w + (omegaup_w_s*sw_s+omegaup_w_n*sw_n) * normalvelocity_sed_aq;
			//HYDRATE PHASE
			auto normalvelocity_h = 0.;
			double omegaup_h_s = 0., omegaup_h_n = 0.;
			if( normalvelocity_h<0.){
				omegaup_h_s = 0.;
				omegaup_h_n = 1.;
			}else{
				omegaup_h_s = 1.;
				omegaup_h_n = 0.;
			}
			double normalvelocity_h_total = normalvelocity_h + (omegaup_sed_aq_s*sh_s+omegaup_sed_aq_n*sh_n) * normalvelocity_sed_aq;

			/****************************************************/
			// RESIDUAL: PW EQN
			auto term_pw = (omegaup_w_s * rhow_s + omegaup_w_n * rhow_n) * normalvelocity_w_total
						 + (omegaup_g_s * rhog_s + omegaup_g_n * rhog_n) * normalvelocity_g_total
						 + (omegaup_h_s * rhoh_s + omegaup_h_n * rhoh_n) * normalvelocity_h_total;
			r_s.accumulate(lfsu_s.child(param.tag.Pw) , 0, +term_pw*face_volume);
			r_n.accumulate(lfsu_n.child(param.tag.Pw) , 0, -term_pw*face_volume);

			// RESIDUAL: SH EQN
			auto term_sh = (omegaup_h_s * rhoh_s + omegaup_h_n * rhoh_n) * normalvelocity_h_total;
			r_s.accumulate(lfsu_s.child(param.tag.Sh) , 0, +term_sh*face_volume);
			r_n.accumulate(lfsu_n.child(param.tag.Sh) , 0, -term_sh*face_volume);

			// RESIDUAL: ENERGY BALANCE
			auto T_ref = 273.15; //K
			auto convectiveflux_HEAT = (  omegaup_g_s * rhog_s * cpg_s * (T_s-T_ref) + omegaup_g_n * rhog_n * cpg_n * (T_n-T_ref)  ) * normalvelocity_g_total
									 + (  omegaup_w_s * rhow_s * cpw_s * (T_s-T_ref) + omegaup_w_n * rhow_n * cpw_n * (T_n-T_ref)  ) * normalvelocity_w_total;
			auto diffusiveflux_HEAT = - kteff_int * ( T_n - T_s )/distance;
			auto term_T = convectiveflux_HEAT + diffusiveflux_HEAT ;
			r_s.accumulate(lfsu_s.child(param.tag.T) , 0, +term_T*face_volume);
			r_n.accumulate(lfsu_n.child(param.tag.T) , 0, -term_T*face_volume);
			/*****************************************************/

			// mobility matrix
			auto M = param.matrix.M();
			// mass diffusion coefficients
	        std::vector<double> D_s(param.tag.Ns,0.);
	        std::vector<double> D_n(param.tag.Ns,0.);
	        std::vector<double> D_int(param.tag.Ns,0.);
	        std::vector<double> Ds_s(param.tag.Ns,0.);
	        std::vector<double> Ds_n(param.tag.Ns,0.);
	        std::vector<double> Ds_int(param.tag.Ns,0.);
	        double tau_s = param.property.Tortuosity(cell_inside,inside_cell_center_local,(*time),(*dt));
	        double tau_n = param.property.Tortuosity(cell_outside,outside_cell_center_local,(*time),(*dt));
	        for(int i=0; i<param.tag.Ns; i++ ){
	        	D_s[i] = std::abs(param.property.AqueousDiffusionCoefficient(cell_inside,inside_cell_center_local,(*time),(*dt),i)
	        					* normal );
	        	D_n[i] = std::abs(param.property.AqueousDiffusionCoefficient(cell_outside,outside_cell_center_local,(*time),(*dt),i)
	        					* normal );
	        	D_int[i] = 2.*(por_s*sw_s*tau_s*D_s[i])*(por_n*sw_n*tau_n*D_n[i]);
	        	if(por_s*sw_s*tau_s*D_s[i]+por_n*sw_n*tau_n*D_n[i] > 0.) D_int[i] *= 1./(por_s*sw_s*tau_s*D_s[i]+por_n*sw_n*tau_n*D_n[i]);

		        Ds_s[i] = std::abs(param.property.SolidsDiffusionCoefficient(cell_inside,inside_cell_center_local,(*time),(*dt),i)
								* normal ); //cm^2 . a^-1
		        Ds_n[i] = std::abs(param.property.SolidsDiffusionCoefficient(cell_outside,outside_cell_center_local,(*time),(*dt),i)
								* normal ); //cm^2 . a^-1
		        Ds_int[i] = 2.*(1.-por_s)*Ds_s[i]*(1.-por_n)*Ds_n[i];
		        if( (1.-por_s)*Ds_s[i]+(1.-por_n)*Ds_n[i] > 0. ) Ds_int[i] *= 1./((1.-por_s)*Ds_s[i]+(1.-por_n)*Ds_n[i]);

	        }
			// upwinding
	        std::vector<double> C_up(param.tag.Ns,0.);
	        for(int i=0; i<param.tag.Ns; i++ ){
	        	C_up[i] = C_s[i];
				if( M[i][i]*normalvelocity_w_total + (1.-M[i][i])*normalvelocity_sed_s < 0.){
					C_up[i] = C_n[i];
				}
	        }
	        // transport operators
	        std::vector<double> L(param.tag.Ns,0.);
	        for(int i=0; i<param.tag.Ns; i++ ){
	        	L[i] = M[i][i] 		* ( C_up[i] * normalvelocity_w_total	- D_int[i]  * ( C_n[i] - C_s[i] )/distance )
	        		 + (1.-M[i][i]) * ( C_up[i] * normalvelocity_sed_s		- Ds_int[i] * ( C_n[i] - C_s[i] )/distance );
	        	if( i==param.tag.CH4 ){
	        		L[i] += (omegaup_g_s*rhog_s + omegaup_g_n*rhog_n)*(1./param.property.MolarMassCH4()) * normalvelocity_g_total;
	        	}
	        }
	        auto UL = operation.multiplyMatrixVector( (*matU),L );

	        /****************************************************/
			//RESIDUALS: SPECIES MASS BALANCE
	        for(int i=0; i<param.tag.nPDE; i++){
				r_s.accumulate(lfsu_s.child(i) , 0, +UL[i]*face_volume);
				r_n.accumulate(lfsu_n.child(i) , 0, -UL[i]*face_volume);
	        }
	        /****************************************************/

	  }

	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_boundary ( const IG& ig,
			  	  	  	  	const LFSU& lfsu, const X& x, const LFSV& lfsv,
							R& r ) const
	  {
	        // References to inside and outside cells
	        const auto& cell_inside  = ig.inside();
			const IndexSet &indexSet = gv.indexSet();
			int inside_cell_number  = indexSet.index(cell_inside);

	        // get geometries
	        auto geo = ig.geometry();
//	        const auto dim = geo.mydimension;
	        const auto dim = gv.dimension;
	        auto geo_inside = cell_inside.geometry();

	        // cell geometries
	        auto ref_el_inside 	= referenceElement(geo_inside);
	        auto inside_cell_center_local 	= ref_el_inside.position(0,0);
	        auto inside_cell_center_global 	= geo_inside.center();

	        // face geometry
	        auto ref_el = referenceElement(geo);
	        auto face_center_local = ref_el.position(0,0);
	        auto face_center_global = geo.center();
	        auto face_volume = geo.volume();
	        auto normal = ig.unitOuterNormal(face_center_local);

	        // distance of cell centers
	        auto d = geo.global(face_center_local);
	        d -= inside_cell_center_global;
	        auto distance = d.two_norm();

			// compute primary species at local self centers
	        std::vector<double> C_s(param.tag.Ns,0.);
	        for(int i=0; i<param.tag.nPrC; i++ ){
	        	C_s[i] = x(lfsu.child(i),0);
	        }

			// evaluate primary species at neighbour center
	        std::vector<double> C_n(param.tag.Ns,0.);
	        C_n = C_s;
	        // BCs
	        std::vector<int> 	bct( param.tag.Ns,ConvectionDiffusionBoundaryConditions::Neumann);
	        std::vector<double>	bcv(param.tag.Ns,0.);
	        for(int i=0; i<param.tag.nPrC; i++){
	        	bct[i] = bc.type( ig , face_center_local, (*time), (*dt), i );
	        	bcv[i] = bc.value(ig , face_center_local, (*time), (*dt), i );
	        	if(bct[i]==ConvectionDiffusionBoundaryConditions::Dirichlet ){
	        		C_n[i] = bcv[i];
//	        		if( i==param.tag.CH4 and ( sg_s - (param.property.EquilibriumConcentrationCH4(pg_s,zg_s,T_s)-C_s[i]) > 0. ) )
//	        			C_n[i] = param.property.EquilibriumConcentrationCH4(pg_s,zg_s,T_s);
	        	}
	        }
			// compute secondary species concentrations
	        if( param.tag.Ne>0 ){ // only if equilibrium reactions are considered
				// evaluate secondary C values as functions of primary C's using the equilibrium constraints
				// log C_s = Ses_inv * log K + Sep_star * log C_p
				std::vector<double> Log10K(param.tag.Ne,0.);
				for(int i=0; i<param.tag.Ne; i++){
					Log10K[i] = param.property.Log10EquilibriumConstant(i);
				}
				auto Log10K_star  = operation.multiplyMatrixVector( (*Ses_inv ), Log10K  );

				std::vector<double> Log10Cp_s(param.tag.nPDE,0.);
				std::vector<double> Log10Cp_n(param.tag.nPDE,0.);
				for(int i=0; i<param.tag.nPDE; i++){
					Log10Cp_s[i] = param.property.Log10Activity(i,C_s);
					Log10Cp_n[i] = param.property.Log10Activity(i,C_n);
				}
				auto Log10Cp_star_s = operation.multiplyMatrixVector( (*Sep_star), Log10Cp_s );
				auto Log10Cp_star_n = operation.multiplyMatrixVector( (*Sep_star), Log10Cp_n );

				auto Fmkg = param.factor.F_Mkg;
				for(int i=param.tag.nIQC; i<param.tag.Ne; i++ ){
					C_s[param.tag.nPDE+i] = std::pow(10.,Log10K_star[i] + Log10Cp_star_s[i])/Fmkg;
					C_n[param.tag.nPDE+i] = std::pow(10.,Log10K_star[i] + Log10Cp_star_n[i])/Fmkg;
				}
	        }

	        // salinity
	        auto salinity_s = param.property.Salinity(C_s[param.tag.Cl_1n]);
	        auto salinity_n = param.property.Salinity(C_n[param.tag.Cl_1n]);

	        // compute phase variables at local self center
	        auto pw_s = x(lfsu.child(param.tag.Pw),0);
	        auto sh_s = x(lfsu.child(param.tag.Sh),0);
	        auto sg_s = x(lfsu.child(param.tag.Sg),0);
	        auto T_s  = x(lfsu.child(param.tag.T ),0);
	        // secondary phase variables and properties at local self center
	        auto sw_s = 1.-sg_s-sh_s;
	        auto pc_s = param.property.CapillaryPressure(cell_inside,inside_cell_center_local,(*time),(*dt),sw_s,sh_s);
	        auto pg_s = pw_s + pc_s;
	        auto peff_s = (sw_s*pw_s+sg_s*pg_s)/(1.-sh_s);
	        auto por_s = param.property.SedimentPorosity(cell_inside,inside_cell_center_local,(*time),(*dt));
	        auto krg_s = param.property.krg(cell_inside,inside_cell_center_local,(*time),(*dt),sw_s,sh_s);
	        auto krw_s = param.property.krw(cell_inside,inside_cell_center_local,(*time),(*dt),sw_s,sh_s);;
	        auto K_s = std::abs(param.property.SedimentPermeability(cell_inside,inside_cell_center_local,(*time),(*dt),sh_s)
	        		 * ig.unitOuterNormal(face_center_local));
	        auto zg_s = param.property.CompressibilityFactorCH4(pg_s,T_s);
	        auto mug_s = param.property.GasPhaseViscosity(pg_s,zg_s,T_s);
	        auto muw_s = param.property.WaterPhaseViscosity(pw_s,T_s,salinity_s);
	        auto rhog_s = param.property.GasPhaseDensity(pg_s,zg_s,T_s);
	        auto rhow_s = param.property.WaterPhaseDensity(pw_s,T_s,salinity_s);
	        auto rhoh_s = param.property.HydratePhaseDensity(peff_s,T_s);
	        auto rhos_s = param.property.SedimentDensity(peff_s,T_s);
	        auto cpg_s = param.property.GasPhaseCp(pg_s,zg_s,T_s);
	        auto cpw_s = param.property.WaterPhaseCp(pw_s,T_s,salinity_s);
	        auto ktg_s = param.property.GasPhaseThermalConductivity(pg_s,zg_s,T_s);
	        auto ktw_s = param.property.WaterPhaseThermalConductivity(pw_s,T_s,salinity_s);
	        auto kth_s = param.property.HydatePhaseThermalConductivity(peff_s,T_s);
	        auto kts_s = param.property.SedimentThermalConductivity(peff_s,T_s);
	        auto kteff_s = (1.-por_s) * kts_s
	        			 + por_s * sw_s * ktw_s
						 + por_s * sg_s * ktg_s
						 + por_s * sh_s * kth_s;

	        // compute phase variables at local boundary face center
	        auto bct_pw = bc.type( ig , face_center_local, (*time), (*dt), param.tag.Pw );
	        auto bct_sg = bc.type( ig , face_center_local, (*time), (*dt), param.tag.Sg );
	        auto bct_T  = bc.type( ig , face_center_local, (*time), (*dt), param.tag.T  );
	        auto bcv_pw = bc.value(ig , face_center_local, (*time), (*dt), param.tag.Pw );
	        auto bcv_sg = bc.value(ig , face_center_local, (*time), (*dt), param.tag.Sg );
	        auto bcv_T  = bc.value(ig , face_center_local, (*time), (*dt), param.tag.T  );
	        auto pw_n = pw_s;
	        if( bct_pw == ConvectionDiffusionBoundaryConditions::Dirichlet ) pw_n = bcv_pw;
	        auto sg_n = sg_s;
//	        if( bct_sg == ConvectionDiffusionBoundaryConditions::Dirichlet and ( sg_s - (param.property.EquilibriumConcentrationCH4(pg_s,zg_s,T_s)-C_s[param.tag.CH4]) > 0. ) )
//	        	sg_n = bcv_sg;
	        auto T_n  = T_s;
	        if( bct_T  == ConvectionDiffusionBoundaryConditions::Dirichlet ) T_n  = bcv_T ;
	        auto sh_n = sh_s;
	        // secondary phase variables and properties at local neighbour center
	        auto sw_n = 1.-sg_n-sh_n;
	        auto pc_n = param.property.CapillaryPressure(cell_inside,inside_cell_center_local,(*time),(*dt),sw_n,sh_n);
	        auto pg_n = pw_n + pc_n;
	        auto peff_n = (sw_n*pw_n+sg_n*pg_n)/(1.-sh_n);
	        auto por_n = param.property.SedimentPorosity(cell_inside,inside_cell_center_local,(*time),(*dt));
	        auto krg_n = param.property.krg(cell_inside,inside_cell_center_local,(*time),(*dt),sw_n,sh_n);
	        auto krw_n = param.property.krw(cell_inside,inside_cell_center_local,(*time),(*dt),sw_n,sh_n);;
	        auto K_n = std::abs(param.property.SedimentPermeability(cell_inside,inside_cell_center_local,(*time),(*dt),sh_n)
	        		 * ig.unitOuterNormal(face_center_local));
	        auto zg_n = param.property.CompressibilityFactorCH4(pg_n,T_n);
	        auto mug_n = param.property.GasPhaseViscosity(pg_n,zg_n,T_n);
	        auto muw_n = param.property.WaterPhaseViscosity(pw_n,T_n,salinity_n);
	        auto rhog_n = param.property.GasPhaseDensity(pg_n,zg_n,T_n);
	        auto rhow_n = param.property.WaterPhaseDensity(pw_n,T_n,salinity_n);
	        auto rhoh_n = param.property.HydratePhaseDensity(peff_n,T_n);
	        auto rhos_n = param.property.SedimentDensity(peff_n,T_n);
	        auto cpg_n = param.property.GasPhaseCp(pg_n,zg_n,T_n);
	        auto cpw_n = param.property.WaterPhaseCp(pw_n,T_n,salinity_n);
	        auto ktg_n = param.property.GasPhaseThermalConductivity(pg_n,zg_n,T_n);
	        auto ktw_n = param.property.WaterPhaseThermalConductivity(pw_n,T_n,salinity_n);
	        auto kth_n = param.property.HydatePhaseThermalConductivity(peff_n,T_n);
	        auto kts_n = param.property.SedimentThermalConductivity(peff_n,T_n);
	        auto kteff_n = (1.-por_n) * kts_n
	        			 + por_n * sw_n * ktw_n
						 + por_n * sg_n * ktg_n
						 + por_n * sh_n * kth_n;

	        // properties at interface
	        auto K_int  = 2.*K_s*K_n;
	        if(K_s+K_n>0) K_int *= 1./( K_s+K_n );
	        auto kteff_int = 2.*kteff_s*kteff_n;
	        if(kteff_s+kteff_n>0.) kteff_int *= 1./(kteff_s+kteff_n);
	        auto rhog_int = 0.5*(rhog_s + rhog_n);
	        auto rhow_int = 0.5*(rhow_s + rhow_n);

	        // phase velocities
			auto normalvelocity_sed_s = param.property.SolidsBurialVelocity(cell_inside,inside_cell_center_local,(*time),(*dt))
					* ig.unitOuterNormal(face_center_local) ;
			auto normalvelocity_sed_aq = param.property.AqueousBurialVelocity(cell_inside,inside_cell_center_local,(*time),(*dt))
					* ig.unitOuterNormal(face_center_local) ;
			double omegaup_sed_s_s = 0., omegaup_sed_s_n = 0.;
			if( (-1.)*normalvelocity_sed_s>0.){
				omegaup_sed_s_s = 0.;
				omegaup_sed_s_n = 1.;
			}else{
				omegaup_sed_s_s = 1.;
				omegaup_sed_s_n = 0.;
			}
			double omegaup_sed_aq_s = 0., omegaup_sed_aq_n = 0.;
			if( (-1.)*normalvelocity_sed_aq>0.){
				omegaup_sed_aq_s = 0.;
				omegaup_sed_aq_n = 1.;
			}else{
				omegaup_sed_aq_s = 1.;
				omegaup_sed_aq_n = 0.;
			}
			auto gravity = param.property.g() * ig.unitOuterNormal(face_center_local) ;
			//GAS PHASE
			auto normalpotential_g = (pg_n - pg_s)/distance
								   + rhog_int * gravity ;
			if( bct_sg == ConvectionDiffusionBoundaryConditions::Neumann ) normalpotential_g=0.;
			auto normalvelocity_g_tmp = - K_int
										* (krg_s/mug_s)
										* normalpotential_g;
			double omegaup_g_s = 0., omegaup_g_n = 0.;
//			if( normalvelocity_g_tmp+sg_s*normalvelocity_sed_aq<0.){
			if( normalpotential_g+normalvelocity_sed_aq>0.){
//			if( normalpotential_g>0.){
				omegaup_g_s = 0.;
				omegaup_g_n = 1.;
			}else{
				omegaup_g_s = 1.;
				omegaup_g_n = 0.;
			}
			auto normalvelocity_g = - K_int
									* ( omegaup_g_s * (krg_s/mug_s) + omegaup_g_n * (krg_n/mug_n) )
									* normalpotential_g;
			double normalvelocity_g_total = normalvelocity_g + (omegaup_sed_aq_s*sg_s+omegaup_sed_aq_n*sg_n) * normalvelocity_sed_aq;
//			double normalvelocity_g_total = normalvelocity_g + (omegaup_g_s*sg_s+omegaup_g_n*sg_n) * normalvelocity_sed_aq;
			//WATER PHASE
			auto normalpotential_w = (pw_n - pw_s)/distance
								   + rhow_int * gravity ;
			if( bct_pw == ConvectionDiffusionBoundaryConditions::Neumann ) normalpotential_w=0.;
			auto normalvelocity_w_tmp = - K_int
									* (krw_s/muw_s)
									* normalpotential_w;
			double omegaup_w_s = 0., omegaup_w_n = 0.;
//			if( normalvelocity_w_tmp+sw_s*normalvelocity_sed_aq<0.){
			if( normalpotential_w+normalvelocity_sed_aq>0.){
//			if( normalpotential_w>0.){
				omegaup_w_s = 0.;
				omegaup_w_n = 1.;
			}else{
				omegaup_w_s = 1.;
				omegaup_w_n = 0.;
			}
			auto normalvelocity_w = - K_int
									* ( omegaup_w_s * (krw_s/muw_s) + omegaup_w_n * (krw_n/muw_n) )
									* normalpotential_w;
			double normalvelocity_w_total = normalvelocity_w + (omegaup_sed_aq_s*sw_s+omegaup_sed_aq_n*sw_n) * normalvelocity_sed_aq;
//			double normalvelocity_w_total = normalvelocity_w + (omegaup_w_s*por_s*sw_s+omegaup_w_n*por_n*sw_n) * normalvelocity_sed_aq;
//			if( bct_pw == ConvectionDiffusionBoundaryConditions::Neumann ) normalvelocity_w_total=0.;
			//HYDRATE PHASE
			auto normalvelocity_h = 0.;
			double omegaup_h_s = 0., omegaup_h_n = 0.;
			if( normalvelocity_h<0.){
				omegaup_h_s = 0.;
				omegaup_h_n = 1.;
			}else{
				omegaup_h_s = 1.;
				omegaup_h_n = 0.;
			}
			double normalvelocity_h_total = normalvelocity_h + (omegaup_sed_aq_s*sh_s+omegaup_sed_aq_n*sh_n) * normalvelocity_sed_aq;

			/****************************************************/
			// RESIDUAL: PW EQN
			auto term_pw = (omegaup_w_s * rhow_s + omegaup_w_n * rhow_n) * normalvelocity_w_total
						 + (omegaup_g_s * rhog_s + omegaup_g_n * rhog_n) * normalvelocity_g_total
						 + (omegaup_h_s * rhoh_s + omegaup_h_n * rhoh_n) * normalvelocity_h_total;
			r.accumulate(lfsu.child(param.tag.Pw) , 0, +term_pw*face_volume);

			// RESIDUAL: SH EQN
			auto term_sh = (omegaup_h_s * rhoh_s + omegaup_h_n * rhoh_n) * normalvelocity_h_total;
			r.accumulate(lfsu.child(param.tag.Sh) , 0, +term_sh*face_volume);

			// RESIDUAL: ENERGY BALANCE
			auto T_ref = 273.15; //K
			auto convectiveflux_HEAT = (  omegaup_g_s * rhog_s * cpg_s * (T_s-T_ref) + omegaup_g_n * rhog_n * cpg_n * (T_n-T_ref)  ) * normalvelocity_g_total
									 + (  omegaup_w_s * rhow_s * cpw_s * (T_s-T_ref) + omegaup_w_n * rhow_n * cpw_n * (T_n-T_ref)  ) * normalvelocity_w_total;
			auto diffusiveflux_HEAT = - kteff_int * ( T_n - T_s )/distance;
			if( bct_T == ConvectionDiffusionBoundaryConditions::Neumann ) {
				diffusiveflux_HEAT = - kteff_int * bcv_T;
//				convectiveflux_HEAT = 0.;
			}
			auto term_T = convectiveflux_HEAT + diffusiveflux_HEAT ;
			r.accumulate(lfsu.child(param.tag.T) , 0, +term_T*face_volume);
			/*****************************************************/

			// mobility matrix
			auto M = param.matrix.M();
			// mass diffusion coefficients
	        std::vector<double> D_s(param.tag.Ns,0.);
	        std::vector<double> D_n(param.tag.Ns,0.);
	        std::vector<double> D_int(param.tag.Ns,0.);
	        std::vector<double> Ds_s(param.tag.Ns,0.);
	        std::vector<double> Ds_n(param.tag.Ns,0.);
	        std::vector<double> Ds_int(param.tag.Ns,0.);
	        double tau_s = param.property.Tortuosity(cell_inside,inside_cell_center_local,(*time),(*dt));
	        double tau_n = param.property.Tortuosity(cell_inside,inside_cell_center_local,(*time),(*dt));
	        for(int i=0; i<param.tag.Ns; i++ ){
	        	D_s[i] = std::abs(param.property.AqueousDiffusionCoefficient(cell_inside,inside_cell_center_local,(*time),(*dt),i)
	        					* normal );
	        	D_n[i] = std::abs(param.property.AqueousDiffusionCoefficient(cell_inside,inside_cell_center_local,(*time),(*dt),i)
	        					* normal );
	        	D_int[i] = 2.*(por_s*sw_s*tau_s*D_s[i])*(por_n*sw_n*tau_n*D_n[i]);
	        	if(por_s*sw_s*tau_s*D_s[i]+por_n*sw_n*tau_n*D_n[i] > 0.) D_int[i] *= 1./(por_s*sw_s*tau_s*D_s[i]+por_n*sw_n*tau_n*D_n[i]);

		        Ds_s[i] = std::abs(param.property.SolidsDiffusionCoefficient(cell_inside,inside_cell_center_local,(*time),(*dt),i)
								* normal ); //cm^2 . a^-1
		        Ds_n[i] = std::abs(param.property.SolidsDiffusionCoefficient(cell_inside,inside_cell_center_local,(*time),(*dt),i)
								* normal ); //cm^2 . a^-1
		        Ds_int[i] = 2.*(1.-por_s)*Ds_s[i]*(1.-por_n)*Ds_n[i];
		        if( (1.-por_s)*Ds_s[i]+(1.-por_n)*Ds_n[i] > 0. ) Ds_int[i] *= 1./((1.-por_s)*Ds_s[i]+(1.-por_n)*Ds_n[i]);
	        }
			// upwinding
	        std::vector<double> C_up(param.tag.Ns,0.);
	        for(int i=0; i<param.tag.Ns; i++ ){
	        	C_up[i] = C_s[i];
				if( M[i][i]*normalvelocity_w_total + (1.-M[i][i])*normalvelocity_sed_s < 0.){
					C_up[i] = C_n[i];
				}
	        }
	        // transport operators
	        std::vector<double> L(param.tag.Ns,0.);
	        for(int i=0; i<param.tag.Ns; i++ ){
	        	if( bct[i]==ConvectionDiffusionBoundaryConditions::Neumann ){
					L[i] = M[i][i] 		* ( C_up[i] * normalvelocity_w_total	- D_int[i]  * bcv[i] )
						 + (1.-M[i][i]) * ( C_up[i] * normalvelocity_sed_s		- Ds_int[i] * bcv[i] );
	        	}else{
					L[i] = M[i][i] 		* ( C_up[i] * normalvelocity_w_total	- D_int[i]  * ( C_n[i] - C_s[i] )/distance )
						 + (1.-M[i][i]) * ( C_up[i] * normalvelocity_sed_s		- Ds_int[i] * ( C_n[i] - C_s[i] )/distance );
	        	}
				if( i==param.tag.CH4 ){
					L[i] += (omegaup_g_s*rhog_s + omegaup_g_n*rhog_n)*(1./param.property.MolarMassCH4()) * normalvelocity_g_total;
				}
	        }
	        auto UL = operation.multiplyMatrixVector( (*matU),L );

	        /****************************************************/
			//RESIDUALS: SPECIES MASS BALANCE
	        for(int i=0; i<param.tag.nPDE; i++){
				r.accumulate(lfsu.child(i) , 0, +UL[i]*face_volume);
	        }
	        /****************************************************/

	  }

};

#endif /* OPERATORS_LOCALOPERATOR_HH_ */
