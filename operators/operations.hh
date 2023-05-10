/*
 * operations.hh
 *
 *  Created on: Aug 30, 2021
 *      Author: sgupta
 */

#ifndef OPERATORS_OPERATIONS_HH_
#define OPERATORS_OPERATIONS_HH_

class Operations{
private:
  	Tags tag;
  	Matrices matrix;
public:

	// PRIMARY STOICHIOMETRIC MATRIX (size: NexnPDE)
	std::vector< std::vector<double> >
	S_eqb_primary() const {

		int nrow = Tags::Ne;
		int ncol = Tags::nPDE;
		if( Tags::Ne == 0 ){ //dummy values
			nrow = 1;
			ncol = 1;
		}
		std::vector< std::vector<double> > Sep(nrow);
		for( int i=0;i<nrow;i++ ){
		  Sep[i] = std::vector<double>(ncol,0.);
		}

		if(Tags::Ne>0){
			auto Se = matrix.S_eqb();

			for(int i=0; i<Sep.size(); i++ ){
				for(int j=0; j<Sep.at(0).size(); j++ ){
					Sep[i][j] = Se[i][j];
				}
			}
		}

		return Sep;
	}

	// SECONDARY STOICHIOMETRIC MATRIX (size: NexNe)
	std::vector< std::vector<double> >
	S_eqb_secondary() const {

		int nrow = Tags::Ne;
		int ncol = Tags::Ne;
		if( Tags::Ne == 0 ){ //dummy values
			nrow = 1;
			ncol = 1;
		}

		std::vector< std::vector<double> > Ses(nrow);
		for( int i=0;i<nrow;i++ ){
		  Ses[i] = std::vector<double>(ncol,0.);
		}

		if( Tags::Ne >0 ){
			auto Se = matrix.S_eqb();

			for(int i=0; i<Ses.size(); i++ ){
				for(int j=0; j<Ses.size(); j++ ){
					Ses[i][j] = Se[i][Tags::nPDE+j];
				}
			}
		}

		return Ses;
	}

	// SECONDARY STOICHIOMETRIC MATRIX INVERSE (size: NexNe)
	std::vector< std::vector<double> >
	S_eqb_secondary_inverse() const {
		auto Ses = S_eqb_secondary();
		auto Ses_inv = Ses; // initialize
		if(Tags::Ne>0){ //invert only if equilibrium reactions are considered
			Ses_inv = inverseMatrix( Ses );
		}
		return Ses_inv;
	}

	// Sep^star = - Ses_inv * Sep -> (size: NexnPDE)
	std::vector< std::vector<double> >
	S_eqb_primary_star() const {

		if(Tags::Ne>0){
			auto Ses_inv = S_eqb_secondary_inverse();
			auto Sep = S_eqb_primary();
			auto Sep_star = multiplyMatrixMatrix( Ses_inv,Sep );

			for(int i=0; i<Sep_star.size(); i++){
				for(int j=0; j<Sep_star.at(0).size(); j++){
					Sep_star[i][j] *= -1.;
				}
			}

			return Sep_star;

		}else{ //dummy
			std::vector< std::vector<double> > Sep_star(1);
			for( int i=0;i<1;i++ ){
				Sep_star[i] = std::vector<double>(1,0.);
			}
			return Sep_star;
		}
	}

	// COMPONENT MATRIX (size: nPDExNs)
	// U = [ I | Sep_star^t ]
	std::vector< std::vector<double> >
	ComponentMatrix() const {

		std::vector< std::vector<double> > u(Tags::nPDE); //nPDE rows, Ns columns
		for( int i=0;i<Tags::nPDE;i++ ){
		  u[i] = std::vector<double>(Tags::Ns,0.);
		}

		for( int i=0; i<Tags::nPDE; i++ ){
			u[i][i] = 1.;
		}

		if(Tags::Ne>0){ //if Ne=0, nPDE==Ns
			auto Sep_star = S_eqb_primary_star();
			auto Sep_star_t = transposeMatrix(Sep_star);

			for( int i=0; i<Sep_star_t.size(); i++ ){
				for( int j=0; j<Sep_star_t.at(0).size(); j++ ){
					u[i][Tags::nPDE+j] = Sep_star_t[i][j];
				}
			}
		}

		return u;
	}

	void checkComponentMatrix() const {
		if( Tags::Ne>0 ){
			// component matrix
			auto U = ComponentMatrix();
			// Equilibrium stoichiometric matrix
			auto Se = matrix.S_eqb();
			auto Se_t = transposeMatrix(Se);
			// check U*Se_T = 0
			auto u = multiplyMatrixMatrix(U,Se_t);
			std::string name = "U x Se^t";
			printMatrix(name,u);

			for(int i=0; i<u.size(); i++){
				for(int j=0; j<u.at(0).size(); j++){
					if( u[i][j] != 0 ){
						std::cout<< "component matrix check failed! " << std::endl;
						exit(0);
					}
				}
			}
		}else{
			std::cout<< "component matrix check cannot be done because Ne=0." << std::endl;
		}
	}

	/********************************************************************************/

	std::vector<double>
	multiplyMatrixVector(std::vector<std::vector<double> > &Matrix,
						 std::vector<double> &Vector) const {
		int nrows = Matrix.size();
		int ncols = Matrix.at(0).size();

		std::vector<double> Target(nrows,0.);

		for(int i=0;i<nrows;i++){
			Target.at(i) = 0.;
			for(int j=0;j<ncols;j++){
				Target.at(i) += Matrix.at(i).at(j) * Vector.at(j);
			}
		}

		return Target;
	}

	double multiplyVectorVector( std::vector<double> &A,
						  	  	 std::vector<double> &B) const {
		int dim = A.size();
		double target = 0.;
		for(int i=0;i<dim;i++){
			target += A.at(i)*B.at(i);
		}
		return target;
	}

	std::vector<std::vector<double> >
	multiplyMatrixMatrix( std::vector<std::vector<double> > &MatrixA,
					  	  std::vector<std::vector<double> > &MatrixB) const {

		std::vector<std::vector<double>> MatrixAB(MatrixA.size());
		for( int i=0;i<MatrixA.size();i++ ){
			MatrixAB[i] = std::vector<double>(MatrixB.at(0).size(),0.);
		}
		for(int i=0;i<MatrixA.size();i++){
			for(int j=0;j<MatrixB.at(0).size();j++){
				MatrixAB.at(i).at(j) = 0.;
				for(int k=0;k<MatrixA.at(0).size();k++){
					MatrixAB.at(i).at(j) += MatrixA.at(i).at(k) * MatrixB.at(k).at(j);
				}
			}
		}
		return MatrixAB;
	}

	std::vector<std::vector<double> >
	transposeMatrix( std::vector<std::vector<double> > &Matrix ) const {
		int nrows = Matrix.size();
		int ncols = Matrix.at(0).size();

		std::vector<std::vector<double>> Matrix_t(ncols);
		for( int i=0;i<ncols;i++ ){
			Matrix_t[i] = std::vector<double>(nrows,0.);
		}

		for(int i=0;i<nrows;i++){
			for(int j=0;j<ncols;j++){
				Matrix_t.at(j).at(i) = Matrix.at(i).at(j);
			}
		}

		return Matrix_t;
	}

	std::vector<std::vector<double> >
	inverseMatrix(std::vector<std::vector<double> > &Matrix ) const {
		int nrows = Matrix.size();
		int ncols = Matrix.at(0).size();

		std::vector<std::vector<double>> mat(nrows);
		for( int i=0;i<nrows;i++ ){
		  mat[i] = std::vector<double>(2*ncols,0.);
		}

		for( int i=0; i<nrows; i++ ){
			for(int j=0;j<ncols;j++ ){
				mat[i][j] = Matrix[i][j];
			}
			for(int j=ncols;j<2*ncols;j++ ){
				mat[i][ncols+i] = 1.;
			}
		}

		auto mat_tmp = mat;
		for( int i = 0; i<nrows; i++ ){
			double pivot = mat[i][i];
//			if (pivot == 0){
//				for( int i = 0; i<nrows; i++ ){
//					for( int j =0; j<2*ncols; j++ ){
//						mat[i][j] =
//					}
//				}
//			}
//			pivot = mat[i][i];
			for( int j =0; j<2*ncols; j++ ){
				mat[i][j] *= 1./pivot;
			}
			for( int k = 0; k<nrows; k++ ){
				for( int j =0; j<2*ncols; j++ ){
					if( k!=i){
						mat[k][j] = mat[k][j] - mat_tmp[k][i]*mat[i][j];
					}
				}
			}
			mat_tmp = mat;
		}

		std::vector<std::vector<double>> Matrix_inv(nrows);
		for( int i=0;i<nrows;i++ ){
			Matrix_inv[i] = std::vector<double>(ncols,0.);
		}
		for( int i = 0; i<nrows; i++ ){
			for( int j = 0; j<ncols; j++ ){
				Matrix_inv[i][j] = mat[i][j+ncols];
			}
		}

//		name = "check";
//		std::vector<std::vector<double>> I(nrows);
//		for( int i = 0; i<nrows; i++ ){
//			I[i] = std::vector<double>(ncols,0.);
//		}
//		multiplyMatrixMatrix(Matrix,Matrix_inv,I);
//		printMatrix(name,I);

		return Matrix_inv;
	}

	void printMatrix( std::string &name,
			          std::vector<std::vector<double>> &Matrix ) const {
		std::cout<< "matrix : " << name << '\n';
		for( int i=0;i<Matrix.size();i++ ){
			for( int j=0;j<Matrix.at(0).size();j++ ){
				std::cout<< Matrix[i][j]  << '\t';
			}
			std::cout<<std::endl;
		}
	}

	void printVector( std::string &name,
			          std::vector<double> &Vector ) const {
		std::cout<< "vector : " << name << '\n';
		for( int i=0;i<Vector.size();i++ ){
				std::cout<< Vector[i]  << '\t';
		}
		std::cout<<std::endl;
	}

	void check() const {
		std::vector<std::vector<double>> matA(3);
		for( int i=0;i<3;i++ ){
		  matA[i] = std::vector<double>(5,0.);
		}
		matA[0][0] = 2.0; matA[0][1] = 1.0; matA[0][2] = 1.0; matA[0][3] = 0.0; matA[0][4] = 2.5;
		matA[1][0] = 0.0; matA[1][1] = 0.5; matA[1][2] = 3.0; matA[1][3] = 1.5; matA[1][4] = 1.0;
		matA[2][0] = 2.0; matA[2][1] = 2.5; matA[2][2] = 1.0; matA[2][3] = 1.0; matA[2][4] = 5.0;

		std::string name = "A";
		printMatrix( name,matA );

		std::vector<std::vector<double>> matB(5);
		for( int i=0;i<5;i++ ){
			matB[i] = std::vector<double>(3,0.);
		}
		matB[0][0] = 1.5; matB[0][1] = 3.0; matB[0][2] = 3.0;
		matB[1][0] = 0.0; matB[1][1] = 1.5; matB[1][2] = 0.0;
		matB[2][0] = 1.0; matB[2][1] = 2.5; matB[2][2] = 2.5;
		matB[3][0] = 3.0; matB[3][1] = 1.5; matB[3][2] = 1.0;
		matB[4][0] = 1.5; matB[4][1] = 0.5; matB[4][2] = 1.5;

		name = "B";
		printMatrix( name,matB );

		std::vector<double> vecA(5,0.);
		vecA[0] = 1.5; vecA[1] = 0.0; vecA[2] = 3.5; vecA[3] = 0.0; vecA[4] = 1.0;

		name = "A";
		printVector( name,vecA );

		std::vector<double> vecB(3,0.);
		vecB[0] = 5.0; vecB[1] = 0.5; vecB[1] = 5.5;

		name = "B";
		printVector( name,vecB );

		auto matAB = multiplyMatrixMatrix(matA,matB);

		name = "matA x matB";
		printMatrix( name,matAB );

		auto matAB_t = transposeMatrix(matAB);

		name = "matAB^t";
		printMatrix( name,matAB_t );

		auto vecC = multiplyMatrixVector(matA,vecA);

		name = "matA x vecA";
		printVector( name,vecC );

		auto vecD = multiplyMatrixVector(matB,vecB);

		name = "matB x vecB";
		printVector( name,vecD );

		auto AxAt = multiplyVectorVector( vecA,vecA );
		std::cout << "vecA * vecAt = " << AxAt << std::endl;

		auto BxBt = multiplyVectorVector( vecB,vecB );
		std::cout << "vecB * vecBt = " << BxBt << std::endl;

		auto matAB_inv = inverseMatrix(matAB);
		name = "matAB^-1";
		printMatrix( name,matAB_inv );
	}

};

#endif /* OPERATORS_OPERATIONS_HH_ */
