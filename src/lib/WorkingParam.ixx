export module WorkingParam;

import std;

import Options;

export struct WorkingParam
{
	double aa1_cutoff = 0;
	double aas_cutoff = 0;
	double aan_cutoff = 0;
	size_t len_upper_bound = 0;
	size_t len_lower_bound = 0;

	size_t len_eff = 0;
	size_t aln_cover_flag = 0;
	size_t min_aln_lenS = 0;
	size_t min_aln_lenL = 0;
	size_t required_aa1 = 0;
	size_t required_aas = 0; /* or aa2 */
	size_t required_aan = 0;

	WorkingParam() = default; // default constructor
	WorkingParam(double a1, double a2 = 0, double an = 0) : aa1_cutoff{ a1 }, aas_cutoff{ a2 }, aan_cutoff{ an } {}

	constexpr void ControlShortCoverage(const int len)
	{
		len_eff = len;
		aln_cover_flag = 0;
		// has alignment coverage control
		if (options.short_coverage > 0.0 || options.min_control > 0) {
			aln_cover_flag = 1;
			min_aln_lenS = std::max<size_t>({
				static_cast<size_t>(len * options.short_coverage),
				(len - options.short_control),
				options.min_control
			});
		}
		if (options.global_identity == 0)
			len_eff = min_aln_lenS;
	}


	constexpr void ControlLongCoverage(const int len2)
	{
		if (aln_cover_flag == 0) return;
		min_aln_lenL = std::max<size_t>({
			static_cast<size_t>(len2 * options.long_coverage),
			(len2 - options.long_control),
			options.min_control
		});
	}


	constexpr void ComputeRequiredBases(const int NAA, const int ss)
	{
		// d: distance, fraction of errors;
		// e: number of errors;
		// g: length of the maximum gap;
		// m: word length;
		// n: sequence length;
		// alignment length = n - g + 1;
		// d = e / (n - g + 1);
		// e >= 1, so that, g <= n + 1 - 1/d
		// word count = (n - g - m + 1) - (e - 1)*m;
		//            = (n - g - m + 1) - (d*(n - g + 1) - 1)*m
		//            = (n - g + 1) - d*m*(n - g + 1)
		//            = (n - g + 1)*(1 - d*m)
		// minimum word count is reached when g == n + 1 - 1/d
		// so, minimum word count = 1/d - m.
		// if g == band_width: word count = (n - band + 1)*(1 - d*m);
		if (options.useDistance()) {
			// int band = options.band_width + 1;
			int invd = int(1.0 / (options.distance_thd + 1E-9));
			// int k = len_eff < invd ? len_eff : invd;
			int ks = len_eff - ss + 1;
			int kn = len_eff - NAA + 1;
			int ks2 = invd - ss;
			int kn2 = invd - NAA;
			// int ks3 = int((len_eff - band + 1.0)*(1.0 - options.distance_thd * ss));
			// int kn3 = int((len_eff - band + 1.0)*(1.0 - options.distance_thd * NAA));
			//if( ks3 > ks2 ) ks2 = ks3;
			//if( kn3 > kn2 ) kn2 = kn3;
			required_aa1 = required_aas = (ks2 < ks ? ks2 : ks);
			required_aan = kn2 < kn ? kn2 : kn;
			if (required_aa1 <= 0) required_aa1 = required_aas = 1;
			if (required_aan <= 0) required_aan = 1;
			//required_aa1 = required_aas = required_aan = 0;
			return;
		}

		double significance = (1.0 - aa1_cutoff);
		size_t i = ss * std::llround(significance * len_eff);
		// (N-K)-K*(1-C)*N = C*K*N-(K-1)*N-K = (C*K-K+1)*N-K
		required_aa1 = std::llround(len_eff - ss - i);
		required_aas = required_aa1;

		size_t j = NAA * std::llround(significance * len_eff);
		required_aa1 = std::llround(len_eff - NAA - j);

		int aa1_old = int(aa1_cutoff * len_eff) - ss + 1;
		int aas_old = int(aas_cutoff * len_eff);
		int aan_old = int(aan_cutoff * len_eff);

		double thd = options.cluster_thd;
		//double rest = (len_eff - ss) / double(len_eff * ss);
		double rest = (len_eff - NAA) / double(len_eff * NAA);
		double thd0 = 1.0 - rest;
		double fnew = 0;
		double fold = 1;
		if (thd > thd0) {
			fnew = (thd - thd0) / rest;
			fold = 1.0 - fnew;
		}
		//printf( "%g %g %g\n", thd, thd0, fnew );

		required_aa1 = (int)(fnew * required_aa1 + fold * aa1_old);
		required_aas = (int)(fnew * required_aas + fold * aas_old);
		required_aan = (int)(fnew * required_aan + fold * aan_old);
	}
};
