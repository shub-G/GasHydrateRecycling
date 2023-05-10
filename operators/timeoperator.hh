/*
 * timeoperator.hh
 *
 *  Created on: Aug 30, 2021
 *      Author: sgupta
 */

#ifndef OPERATORS_TIMEOPERATOR_HH_
#define OPERATORS_TIMEOPERATOR_HH_

template < class GV, class Params >
class TimeOperator
  : public Dune::PDELab::NumericalJacobianApplyVolume< TimeOperator<GV,Params> >,
    public Dune::PDELab::NumericalJacobianVolume	 < TimeOperator<GV,Params> >,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
	const GV& gv;
	const Params&	  param;
	std::vector<std::vector<double>> *matU;
	std::vector<std::vector<double>> *Ses_inv;
	std::vector<std::vector<double>> *Sep_star;
	double 			 *time;
	double 			 *dt;
	constexpr static double eps 	= 1.0e-6;
	constexpr static double pi 		= atan(1.)*4 ;

	Operations operation;

public:
	  // pattern assembly flags
	  enum { doPatternVolume = true };

	  // residual assembly flags
	  enum { doAlphaVolume = true };

	  typedef typename GV::IndexSet IndexSet;

	  // constructor remembers parameters
	  TimeOperator( const GV& gv_, const Params& param_,
			  	  	std::vector<std::vector<double>> *Ses_inv_,
					std::vector<std::vector<double>> *Sep_star_,
			  	  	std::vector<std::vector<double>> *matU_,
					double			*time_	,
					double 			*dt_	 )
	  :  gv(gv_),
		 param(param_),
		 Ses_inv(Ses_inv_),
		 Sep_star(Sep_star_),
		 matU(matU_),
		 time( time_ ),
		 dt( dt_ )
	  {}

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)const{

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
	        auto cell_center_global = geo.global(cell_center_local);
	        auto cell_volume = geo.volume();

	        // compute phase variables at local center
	        auto pw = x(lfsu.child(param.tag.Pw),0);
	        auto sh = x(lfsu.child(param.tag.Sh),0);
	        auto sg = x(lfsu.child(param.tag.Sg),0);
	        auto T  = x(lfsu.child(param.tag.T ),0);

			// compute species concentrations at local center
	        std::vector<double> C(param.tag.Ns,0.);
	        for(int i=0; i<param.tag.nPrC; i++ ){
	        	C[i] = x(lfsu.child(i),0);
	        }
	        if( param.tag.Ne>0 ){ // only if equilibrium reactions are considered
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

				for(int i=param.tag.nIQC; i<param.tag.Ne; i++ ){
					C[param.tag.nPDE+i] = std::pow(10.,Log10K_star[i] + Log10Cp_star[i] )/param.factor.F_Mkg;
				}
	        }

	        // secondary phase variables
	        auto salinity = param.property.Salinity(C[param.tag.Cl_1n]);
	        auto por = param.property.SedimentPorosity(cell,cell_center_local,(*time),(*dt));
	        auto sw = 1.-sg-sh;
	        auto pc = param.property.CapillaryPressure(cell,cell_center_local,(*time),(*dt),sw,sh);
	        auto pg = pw + pc;
	        auto peff = (sw*pw+sg*pg)/(1.-sh);
	        auto zg = param.property.CompressibilityFactorCH4(pg,T);
	        auto rhog = param.property.GasPhaseDensity(pg,zg,T);
	        auto rhow = param.property.WaterPhaseDensity(pw,T,salinity);
	        auto rhoh = param.property.HydratePhaseDensity(peff,T);
	        auto rhos = param.property.SedimentDensity(peff,T);
	        auto cvg = param.property.GasPhaseCv(pg,zg,T);
	        auto cvw = param.property.WaterPhaseCv(pw,T,salinity);
	        auto cvh = param.property.HydratePhaseCv(peff,T);
	        auto cvs = param.property.SedimentCv(peff,T);


	        auto M = param.matrix.M();
	        for(int i=0; i<C.size(); i++){
	        	C[i] *= ( por * sw * M[i][i] + (1.-por) * ( 1. - M[i][i]) );
	        	if( i==param.tag.CH4 ) C[i] += por * sg * (rhog/param.property.MolarMassCH4());
	        }
	        auto UC = operation.multiplyMatrixVector( (*matU),C );

	        /********************
			 * RESIDUALS
			 *******************/

	        // C EQNS: SPECIES MASS BALANCES
	        for(int i=0; i<param.tag.nPDE; i++){
				r.accumulate(lfsu.child(i) , 0, +UC[i]*cell_volume);
	        }

	        //PW EQN: TOTAL MASS BALANCE
	        auto term_pw = por * rhow * sw
	        			 + por * rhog * sg
						 + por * rhoh * sh;
	        r.accumulate(lfsu.child(param.tag.Pw) , 0 , term_pw*cell_volume);

	        //SH EQN: HYDRATE-PHASE MASS BALANCE
	        auto term_sh = por * rhoh * sh;
	        r.accumulate(lfsu.child(param.tag.Sh) , 0 , term_sh*cell_volume);

	        //T EQN: ENERGY BALANCE
	        auto term_T = por * rhow * sw * cvw
	        			+ por * rhog * sg * cvg
						+ por * rhoh * sh * cvh
						+ (1.-por) * rhos * cvs;
	        term_T *= T;
	        r.accumulate(lfsu.child(param.tag.T ) , 0 , term_T*cell_volume);

	  }

};

#endif /* OPERATORS_TIMEOPERATOR_HH_ */
