/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef TARGET_ITERATOR_H_
#define TARGET_ITERATOR_H_

#define DP_STAT

#include "../dp.h"

template<int _n>
struct TargetIterator
{

	TargetIterator(vector<DpTarget>::const_iterator subject_begin, vector<DpTarget>::const_iterator subject_end) :
		next(0),
		n_targets(int(subject_end - subject_begin)),
		subject_begin(subject_begin)
	{
		for (; next < std::min(_n, n_targets); ++next) {
			pos[next] = 0;
			target[next] = next;
			active.push_back(next);
		}
	}

	TargetIterator(vector<DpTarget>::const_iterator subject_begin, vector<DpTarget>::const_iterator subject_end, int i1, int qlen) :
		next(0),
		n_targets(int(subject_end - subject_begin)),
		cols(0),
		subject_begin(subject_begin)
	{
		for (; next < std::min(_n, n_targets); ++next) {
			const DpTarget &t = subject_begin[next];
			pos[next] = i1 - (t.d_end - 1);
			const int j1 = std::min(qlen - 1 - t.d_begin, (int)(t.seq.length() - 1)) + 1;
			cols = std::max(cols, j1 - pos[next]);
			target[next] = next;
			active.push_back(next);
		}
	}

	char operator[](int channel)
	{
		if (pos[channel] >= 0) {
#ifdef DP_STAT
			++live;
#endif
			return subject_begin[target[channel]].seq[pos[channel]];
		} else
			return value_traits.mask_char;
	}

#ifdef __SSSE3__
	__m128i get() 
	{
		int16_t s[8];
#ifdef DP_STAT
		live = 0;
#endif
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			s[channel] = (*this)[channel];
		}
		return _mm_loadu_si128((const __m128i*)s);
	}
#else
	uint64_t get() 
	{
		uint64_t dst = 0;
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			dst |= uint64_t((*this)[channel]) << (8 * channel);
		}
		return dst;
	}
#endif

	bool init_target(int i, int channel)
	{
		if (next < n_targets) {
			pos[channel] = 0;
			target[channel] = next++;
			return true;
		}
		active.erase(i);
		return false;
	}

	bool inc(int channel)
	{
		++pos[channel];
		if (pos[channel] >= (int)subject_begin[target[channel]].seq.length())
			return false;
		return true;
	}

	int pos[_n], target[_n], next, n_targets, cols;
#ifdef DP_STAT
	int live;
#endif
	Static_vector<int, _n> active;
	const vector<DpTarget>::const_iterator subject_begin;
};

#endif
