// #define REPORT_LINE_QUALITY
#ifdef REPORT_LINE_QUALITY
		// Full integral of the line profile, to check the discretization of the output
		// (emission) grid.
		double maxNorm = 0, minNorm = 1e9;
		forActiveLinesDo([&](size_t upper, size_t lower) {
			auto lp = lineProfile(upper, lower, s.T, s.cvv);

			Array lpv(specificIntensity.frequencyv().size());
			for (size_t i = 0; i < lpv.size(); i++)
				lpv[i] = lp(specificIntensity.frequencyv()[i]);

			double norm = TemplatedUtils::integrate<double>(
			                specificIntensity.frequencyv(), lpv);
			DEBUG("on the input intensity grid, line  "
			      << upper << " --> " << lower << " has norm " << norm << endl);
			maxNorm = max(norm, maxNorm);
			minNorm = min(norm, minNorm);
		});
		DEBUG("Max profile norm = " << maxNorm << endl);
		DEBUG("Min profile norm = " << minNorm << endl);
		maxNorm = 0;
		minNorm = 1e9; // any large number will do
		// Integral over contant spectrum, using the LineProfile class. An integration
		// grid is chosen internally, and we check the quality of it here (at least for
		// now.)
		double minFreq = specificIntensity.freqMin();
		double maxFreq = specificIntensity.freqMax();
		Array someFreqs = {minFreq, (minFreq + maxFreq) / 2, maxFreq};
		Spectrum flat(someFreqs, Array(1, someFreqs.size()));
		forActiveLinesDo([&](size_t upper, size_t lower) {
			auto lp = lineProfile(upper, lower, s);
			double norm = lp.integrateSpectrum(flat);
			DEBUG("LineProfile " << upper << " --> " << lower << " has norm "
			                     << norm << endl);
			maxNorm = max(norm, maxNorm);
			minNorm = min(norm, minNorm);
		});
		DEBUG("Max profile norm = " << maxNorm << endl);
		DEBUG("Min profile norm = " << minNorm << endl);
#endif /* REPORT_LINE_QUALITY */
