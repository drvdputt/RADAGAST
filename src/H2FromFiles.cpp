#include "H2FromFiles.h"

// TODO: make model

H2FromFiles::H2FromFiles()
{
}

H2FromFiles::~H2FromFiles() = default;

int H2FromFiles::numLv() const
{
	return 1;
}

EVector H2FromFiles::ev() const
{
	return EVector::Zero(1);
}

EVector H2FromFiles::gv() const
{
	return EVector::Constant(1, 1);
}

EMatrix H2FromFiles::avv() const
{
	return EMatrix::Zero(1, 1);
}

EMatrix H2FromFiles::extraAvv() const
{
	return EMatrix::Zero(1, 1);
}

EMatrix H2FromFiles::cvv(double T, double ne, double np) const
{
	return EMatrix::Zero(1, 1);
}

EVector H2FromFiles::sourcev(double T, double ne, double np) const
{
	return EVector::Zero(1);
}

EVector H2FromFiles::sinkv(double T, double ne, double np) const
{
	return EVector::Zero(1);
}
