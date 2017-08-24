#include "GrainProperties.h"

using Type = GasModule::GrainType;

GrainProperties::SfcInteractionPar::SfcInteractionPar(double EH2, double Es, double EHp, double EHc,
                                                      double aSqrt, double nuH2, double nuHc,
                                                      double F)
                : _valid{true}, _eH2{EH2}, _es{Es}, _eHp{EHp}, _eHc{EHc}, _aSqrt{aSqrt},
                  _nuH2{nuH2}, _nuHc{nuHc}, _f{F}
{
}

GrainProperties::SfcInteractionPar GrainProperties::sfcInteractionPar(Type t)
{
	if (t == Type::CAR)
		return {520, 260, 800, 30000, 14, 3e12, 1.3e13, 1e-10};
	else if (t == Type::SIL)
		return {320, 110, 450, 30000, 14.4, 3e12, 1.3e13, 1e-10};
	else
		return {};
}

bool GrainProperties::heatingAvailable(Type t)
{
	return t == Type::CAR || t == Type::SIL;
}
