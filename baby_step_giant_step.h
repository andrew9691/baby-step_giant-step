#ifndef COMPATALG_H
#define COMPATALG_H

#include <NTL\ZZ.h>
#include <NTL\ZZ_p.h>
#include <NTL\tools.h>
#include <math.h>
#include "Bin-Method.h"
#include <vector>
#include <algorithm>
#include <pair.h>
#include <utility>

std::vector<std::pair<std::size_t, NTL::ZZ>> giant_steps_table(const NTL::ZZ& r, const NTL::ZZ& p, const NTL::ZZ& n)
{
	std::vector<std::pair<std::size_t, NTL::ZZ>> v;
	for (std::size_t i = 1; i <= n; i++)
		v.emplace_back(std::make_pair(i, bin_method_mod(r, i, p)));
	std::sort(v.begin(), v.end(), [](std::pair<std::size_t, NTL::ZZ> x, std::pair<std::size_t, NTL::ZZ> y) {return x.second < y.second; });
	return v;
}

NTL::ZZ_p bs_gs_alg(const NTL::ZZ& a, const NTL::ZZ& b, const NTL::ZZ& p)
{
	NTL::ZZ n = NTL::SqrRoot(p) + 1;
	NTL::ZZ r = NTL::PowerMod(a, n, p);
	std::vector<std::pair<std::size_t, NTL::ZZ>> giant_steps = giant_steps_table(r, p, n);

	std::size_t i, j, average_ind, fst_ind, last_ind;
	for (j = 1; j <= n; j++)
	{
		NTL::ZZ new_baby_step = (b * NTL::PowerMod(a, j, p)) % p;
		average_ind = 0;
		fst_ind = 0;
		last_ind = giant_steps.size() - 1;

		while (fst_ind < last_ind)
		{
			average_ind = fst_ind + (last_ind - fst_ind) / 2;
			if (new_baby_step <= giant_steps[average_ind].second)
				last_ind = average_ind;
			else
				fst_ind = average_ind + 1;
		}

		if (giant_steps[last_ind].second == new_baby_step)
		{
			i = giant_steps[last_ind].first;
			break;
		}
	}

	NTL::ZZ_p::init(p - 1);
	return NTL::conv<NTL::ZZ_p>(n * i - j);
}

#endif
