//#include "NLevel.h"
//#include "Constants.h"
//
//using namespace std;
//
//NLevel::NLevel(const Array& frequencyv) :
//		_frequencyv(frequencyv)
//{
//	// Using NIST data (http://physics.nist.gov/cgi-bin/ASD/lines1.pl?unit=1&line_out=0&bibrefs=1&show_obs_wl=1&show_calc_wl=1&A_out=0&intens_out=1&allowed_out=1&forbid_out=1&conf_out=1&term_out=1&enrg_out=1&J_out=1&g_out=0&spectra=H%20I)
//	// 1s, 2p, 2s
//	// 3, 4, 5
//	_Ev << 0., 82258.9191133, 82258.9543992821,
//			97492.304, 102823.904, 105291.657;
//	// Energy in cm^-1, so multiply with hc
//	_Ev *= Constant::PLANCKLIGHT;
//
//	_gv << 1, 1, 3, 9, 16, 25;
//
//	_nv = Eigen::VectorXd(_Ev.size());
//	_nv(0) = _n;
//
//	double A2p1 = 6.2649e8;
//
//	double A2s1 = 2.495e-6;
//
//	double A31 = 5.5751e7;
//	double A32 = 4.4101e+07;
//	double A32p = 3 * A32 / 4.;
//	double A32s = A32 / 4.;
//
//	double A41 = 1.2785e+07;
//	double A42 = 8.4193e+06;
//	double A42p = 3 * A42 / 4.;
//	double A42s = A42 / 4.;
//	double A43 = 8.9860e+06;
//
//	double A51 = 4.1250e+06;
//	double A52 = 2.5304e+06;
//	double A52p = 3 * A52 / 4.;
//	double A52s = A52 / 4.;
//	double A53 = 2.2008e+06;
//	double A54 = 2.6993e+06;
//
//	_Avv << 0, 0, 0, 0, 0, 0,
//			A2p1, 0, 0, 0, 0, 0,
//			A2s1, 0, 0, 0, 0, 0,
//			A31, A32p, A32s, 0, 0, 0,
//			A41, A42p, A42s, A43, 0, 0,
//			A51, A52p, A52s, A53, A54, 0;
//
//	// approximation used:
//	// A_n,2p = 3 * A_n,2s = 3/4 * A_n,2 from NIST
//
//}
//
//void NLevel::solveBalance(double atomDensity, double electronDensity, double protonDensity, double T,
//		const Array& specificIntensityv, const Array& sourcev, const Array& sinkv)
//{
//	if (specificIntensityv.size() != _frequencyv.size())
//		throw range_error("Given ISRF and wavelength vectors do not have the same size");
//	if (sourcev.size() != _Ev.size() || sinkv.size() != _Ev.size())
//		throw range_error("Source and/or sink term vector(s) of wrong size");
//
//	_n = atomDensity;
//	_ne = electronDensity;
//	_np = protonDensity;
//	_T = T;
//
//	if (_n > 0)
//	{
//		// Calculate Cij (needs to happen bfore the first call to the line profile. Maybe the dependencies need to be made
//		// clearer by adopting a more functional programming-like style)
//		prepareCollisionMatrix();
//
//#ifdef REPORT_LINE_QUALITY
//		double norm = TemplatedUtils::integrate<double>(_frequencyv, lineProfile(1, 0));
//		DEBUG("line profile norm = " << norm << endl);
//#endif
//
//		// Calculate BijPij (needs to be redone at each temperature because the line profile can change)
//		// Also needs the Cij to calculate collisional broadening
//		prepareAbsorptionMatrix(specificIntensityv);
//
//#ifdef PRINT_MATRICES
//		DEBUG("Aij" << endl << _Avv << endl << endl);
//		DEBUG("BPij" << endl << _BPvv << endl << endl);
//		DEBUG("Cij" << endl << _Cvv << endl << endl);
//#endif
//		// Calculate Fij and bi and solve F.n = b
//		solveRateEquations(Eigen::Map<const Eigen::VectorXd>(&sourcev[0], sourcev.size()),
//				Eigen::Map<const Eigen::VectorXd>(&sinkv[0], sinkv.size()), 0);
//	}
//	else
//	{
//		_nv = Eigen::VectorXd::Zero(_Ev.size());
//	}
//}
//
//
//
//
