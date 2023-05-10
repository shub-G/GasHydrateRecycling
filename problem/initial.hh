#ifndef PROBLEM_INITIAL_HH_
#define PROBLEM_INITIAL_HH_

template<typename GV,typename PTree,typename Parameters>
class InitialConditions
{
private:
	  const GV& gv;
	  const PTree& ptree;
	  const Parameters& parameter;
	  const static int dim = GV::dimension;

	  double Pw_SF;
	  double T_SF;
	  double gradT;
	  double rhow_ref;

public:

	  //! construct from grid view
	  InitialConditions ( const GV& gv_, const PTree& ptree_, const Parameters& parameter_ )
	  : gv( gv_ ) , ptree(ptree_), parameter(parameter_)
	  {
		  rhow_ref = 1027.0; //g/L-->kg/m^3
		  double hwc = ptree.get("parameters.water_column_height",(double)800.0);//in m
		  Pw_SF = rhow_ref*9.81*hwc; //Pa
		  Pw_SF *= 1.e-6; // convert to MPa
		  T_SF = ptree.get("parameters.bottom_water_temperature",(double)4.0);//in degC
		  T_SF += 273.15; //convert to K
		  gradT = ptree.get("parameters.regional_thermal_gradient",(double)24.0);//in degC/km
		  gradT *= 1./1000.0; // convert to degC/m
		  gradT *= 1.e-2; //convert to degC/cm
	  }

	  double evaluate (	const typename GV::Traits::template Codim<0>::Entity& element,
			  			const Dune::FieldVector<double,dim>& xlocal,
						int tag ) const {

		  Dune::FieldVector<double,dim> xglobal/*cm*/ = element.geometry().global(xlocal);

		  double icv=0.;
		  if( tag < parameter.tag.nPrC ){
			  icv/*mmol/Lpw or mmol/Lds*/ = 0.;//parameter.property.BottomWaterConcentration(element,xlocal,0.,0.,tag);
		  }

		  if( tag == parameter.tag.Cl_1n ){
			  icv/*mmol/Lpw*/ = parameter.property.BottomWaterConcentration(element,xlocal,0.,0.,tag);
		  }

		  if( tag == parameter.tag.OM ){
			  icv/*mmol/Lpw or mmol/Lds*/ = 0.;//parameter.property.BottomWaterConcentration(element,xlocal,0.,0.,tag);
		  }

		  if( tag == parameter.tag.Pw ){
			  icv/*MPa*/ = Pw_SF/*MPa*/
					  	 + rhow_ref/*g/L*/
						   * abs(parameter.property.g()[dim-1])/*(MPa.L)/(g.cm)*/
			  	  	  	   * (parameter.domain.SeaFloorHeight(xglobal[0])/*cm*/-xglobal[dim-1]/*cm*/)
//						 + parameter.property.SedimentPorosity(element,xlocal,0.,0.)
//						   * 1.0/*Sw0*/
//						   * parameter.property.BurialVelocity(element,xlocal,0.,0.)[dim-1]
//						   * parameter.property.WaterPhaseViscosity(Pw_SF,T_SF)/parameter.property.SedimentPermeability(element,xlocal,0.,0.,0.)[dim-1]
//						   *(parameter.domain.SeaFloorHeight(xglobal[0])/*cm*/-xglobal[dim-1]/*cm*/)
						 ;
		  }

		  if( tag == parameter.tag.Sg ){
			  icv/*-*/ = 0.0;
		  }

		  if( tag == parameter.tag.Sh ){
			  icv/*-*/ = 0.0;
		  }

		  if( tag == parameter.tag.T  ){
			  icv/*K*/ = T_SF/*K*/ + gradT/*degC/cm*/*(parameter.domain.SeaFloorHeight(xglobal[0])/*cm*/-xglobal[dim-1]/*cm*/);
		  }

		  return icv;
	  }

	  //! get a reference to the grid view
	  inline const GV& getGridView () {return gv;}
};

#endif /* PROBLEM_INITIAL_HH_ */
