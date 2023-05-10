#ifndef PROBLEM_BOUNDARY_HH_
#define PROBLEM_BOUNDARY_HH_

struct ConvectionDiffusionBoundaryConditions
{
  enum Type { Dirichlet=1, Neumann=-1 };
};

template<typename GV,typename PTree,typename Parameters>
class BoundaryConditions
{
private :
	const GV& gv ;
	const PTree& ptree;
	const Parameters& parameter;
	InitialConditions<GV,PTree,Parameters> ic;
	const static int dim = GV::dimension;

	double gradT;

	typedef ConvectionDiffusionBoundaryConditions::Type BCType;

public :

	// ! construct from gridview
	BoundaryConditions (const GV& gv_,const PTree& ptree_,const Parameters& parameter_ )
	: gv ( gv_ ),
	  ptree(ptree_),
	  parameter(parameter_),
	  ic( gv_,ptree_,parameter_)
	{
		  gradT = ptree.get("parameters.regional_thermal_gradient",(double)24.0);//in degC/km
		  gradT *= 1./1000.0; // convert to degC/m
		  gradT *= 1.e-2; //convert to degC/cm
	}

	/* boundary types */
	template<typename I>
	int type ( I& intersection,
			   const Dune::FieldVector<double,dim-1>& xlocal,
			   double time/*a*/,
			   double dt/*a*/,
			   int tag) const {

		auto xglobal = intersection.geometry().global(xlocal);//cm

		int bct = BCType::Neumann; //DEFAULT is NEUMANN

		if ( parameter.domain.isTopBoundary(xglobal) ){
			bct = BCType::Dirichlet;
		}
		return bct;
	}


	/* boundary values */
	template<typename I>
	double value (	I& intersection,
					const Dune::FieldVector<double,dim-1>& xlocal,
					double time/*a*/,
					double dt/*a*/,
					int tag) const {

		double bcvalue= 0.; // DEFAULT is NEUMANN 0

		auto xglobal = intersection.geometry().global(xlocal);//cm
		const auto& element = intersection.inside();
        auto element_xlocal = referenceElement(element.geometry()).position(0,0);

        auto t_new =time+dt;

		if ( parameter.domain.isTopBoundary(xglobal) ){
			if(tag<parameter.tag.nPrC) bcvalue = parameter.property.BottomWaterConcentration(element,element_xlocal,time,dt,tag);
			else bcvalue = ic.evaluate(element,element_xlocal,tag);
		}

		if( parameter.domain.isBottomBoundary(xglobal) and tag==parameter.tag.T ){
			bcvalue = gradT;
		}

		return bcvalue;
	}

	// ! get a reference to the gridview
	inline const GV& getGridView () { return gv ; }
};

#endif /* PROBLEM_BOUNDARY_HH_ */
