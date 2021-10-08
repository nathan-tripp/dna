#include "sequence.hpp"

long double 
sequence_t::at_ratio(const basepair_t* start, const basepair_t* end)
{
	std::size_t			adenine(0), thymine(0), guanine(0), cytosine(0);
	const std::size_t	total(end - start);
	long double			ratio(0);

	if (nullptr == start || nullptr == end || end < start)
		throw std::runtime_error("dna_t::at_ratio(): invalid parameter(s)");

	for (std::size_t idx = 0; idx < total; idx++) {
		basepair_t  bp(start[idx]);
		uint8_t     left(NUCLEO_LEFT(bp));
		uint8_t     right(NUCLEO_RIGHT(bp));

		count_nucleotides(left, adenine, cytosine, guanine, thymine);
		count_nucleotides(right, adenine, cytosine, guanine, thymine);
	}

//	if (adenine != thymine || guanine != cytosine)
//		throw std::runtime_error("invalid molecule (mismatched pairs");

	ratio = double(adenine + thymine) / double(adenine + thymine + guanine + cytosine);
	ratio *= 100.0;
	return ratio;
}

long double 
sequence_t::gc_ratio(const basepair_t* start, const basepair_t* end)
{
	std::size_t 		adenine(0), thymine(0), guanine(0), cytosine(0);
	const std::size_t	total(end - start);
	long double	 		ratio(0);

	if (nullptr == start || nullptr == end || end < start)
		throw std::runtime_error("dna_t::gc_ratio(): invalid parameter(s)");

	for (std::size_t idx = 0; idx < total; idx++) {
		basepair_t 	bp(start[idx]);
		uint8_t		left(NUCLEO_LEFT(bp));
		uint8_t		right(NUCLEO_RIGHT(bp));

		count_nucleotides(left, adenine, cytosine, guanine, thymine);
		count_nucleotides(right, adenine, cytosine, guanine, thymine);
	}

//	if (adenine != thymine || guanine != cytosine)
//		throw std::runtime_error("invalid molecule (mismatched pairs");

	ratio = double(guanine + cytosine) / double(adenine + thymine + guanine + cytosine);
	ratio *= 100.0;
	return ratio;
}

void 
sequence_t::count_nucleotides(const uint8_t val, std::size_t& adenine, std::size_t& cytosine, std::size_t& guanine, std::size_t& thymine)
{
	switch (val) {
	case ADENINE:
		adenine++;
		break;
	case CYTOSINE:
		cytosine++;
		break;
	case GUANINE:
		guanine++;
		break;
	case THYMINE:
		thymine++;
		break;
	default: // should be impossible
		throw std::runtime_error("dna_t::count_nucleotides(): invalid nucleotide encountered");
		break;
	}

	return;
}

bool 
sequence_t::is_at_rich(const basepair_t* start, const basepair_t* end)
{
	long double at(at_ratio(start, end));
	long double gc(gc_ratio(start, end));

	if (at > gc)
		return true;

	return false;
}

bool 
sequence_t::is_gc_rich(const basepair_t* start, const basepair_t* end)
{
	return ! is_at_rich(start, end);
}

inline basepair_t* 
sequence_t::create_sequence(seq_element_t& seq)
{
	basepair_t* ret(nullptr);
	const std::size_t	length(seq.mer_length * sizeof(basepair_t));
	const std::size_t	seq_length(length * seq.repeat_cnt);
	std::size_t			at(0);
	std::size_t			gc(0);
	long double			calc_at(0.0);
	long double			calc_gc(0.0);

	if (seq.mer_length > std::numeric_limits< std::size_t >::max() / sizeof(basepair_t))
		throw std::runtime_error("dna_t::create_sequence(): multiplicative overflow");
	if (length > std::numeric_limits< std::size_t >::max() / seq.repeat_cnt)
		throw std::runtime_error("dna_t::create_sequence(): multiplicative overflow");

	ret = new basepair_t[seq.mer_length * seq.repeat_cnt];
	at = static_cast< std::size_t >((seq.at_ratio / 100) * seq.mer_length);
	gc = static_cast< std::size_t >((seq.gc_ratio / 100) * seq.mer_length);

	std::memset(ret, 0, seq_length);

	for (std::size_t rep = 0; rep < seq.repeat_cnt; rep += seq.mer_length) {
		std::size_t tmp_at(at);
		std::size_t tmp_gc(gc);

		for (std::size_t idx = 0; idx < seq.mer_length; idx++) {
			if (0 != tmp_at) {
				if (m_rnd.byte() & 1) {
					SET_NUCLEO_LEFT(ret[rep + idx], ADENINE);
					SET_NUCLEO_RIGHT(ret[rep + idx], THYMINE);
				}
				else {
					SET_NUCLEO_LEFT(ret[rep + idx], THYMINE);
					SET_NUCLEO_RIGHT(ret[rep + idx], ADENINE);
				}
				tmp_at--;

			}
			else {
				if (m_rnd.byte() & 1) {
					SET_NUCLEO_LEFT(ret[rep + idx], GUANINE);
					SET_NUCLEO_RIGHT(ret[rep + idx], CYTOSINE);
				}
				else {
					SET_NUCLEO_LEFT(ret[rep + idx], CYTOSINE);
					SET_NUCLEO_RIGHT(ret[rep + idx], GUANINE);
				}

				tmp_gc--;
			}
		}
	}

	for (std::size_t idx = 0; idx < seq.repeat_cnt; idx++) {
		calc_at = at_ratio(ret, ret + length);
		calc_gc = gc_ratio(ret, ret + length);

		if (calc_at != seq.at_ratio)
			seq.at_ratio = calc_at;
		if (calc_gc != seq.gc_ratio)
			seq.gc_ratio = calc_gc;
	}

	/* XXX something is getting corrupted somewhere-- this should be the same result as the above loop.

	calc_at = at_ratio(ret, ret+seq_length);
	calc_gc = gc_ratio(ret, ret+seq_length);

	if (calc_at != seq.at_ratio)
		throw std::runtime_error("dna_t::create_sequence(): fragment has a different AT-ratio than the whole");
	if (calc_gc != seq.gc_ratio)
		throw std::runtime_error("dna_t::create_sequence(): fragment has a different GC-ratio than the whole");
	*/

	return ret;
}

sequence_t::sequence_t(const std::size_t max_length)
{
	const std::size_t	osize(get_random(3, 8));
	std::size_t			offset(0);
	std::size_t			remain(max_length);

	m_origin.count	= osize;
	m_origin.origin = new seq_element_t[osize];
	
	std::memset(m_origin.origin, 0, osize * sizeof(seq_element_t));

	for (std::size_t idx = 0; idx < osize; idx++) {
		long double at((m_rnd.byte() % 100) * 1.0);
		long double gc(100.0 - at);

		if (0 != max_length) {
			std::size_t length(0);
			std::size_t cnt(0);

			do {
				length	= get_random(5, 10);
				cnt		= get_random(5, 10);

				if (length > std::numeric_limits< std::size_t >::max() / cnt)
					throw std::invalid_argument("sequence_t::sequence_t(): impossible exception (multiplicative overflow)");

			} while (length * cnt > remain);

			remain -= length * cnt;

			m_origin.origin[idx].mer_length = length;
			m_origin.origin[idx].repeat_cnt = cnt;
		
		} else {
			m_origin.origin[idx].mer_length = get_random(5, 10);
			m_origin.origin[idx].repeat_cnt = get_random(5, 10);
		}

		m_origin.origin[idx].at_ratio	= at;
		m_origin.origin[idx].gc_ratio	= gc;
		m_origin.origin[idx].ptr		= create_sequence(m_origin.origin[idx]);
	}

	return;
}

sequence_t::~sequence_t(void)
{

	for (std::size_t idx = 0; idx < m_origin.count; idx++) {
		delete[] m_origin.origin[idx].ptr;
		m_origin.origin[idx].ptr = nullptr;
	}

	delete[] m_origin.origin;
	m_origin.origin = nullptr;
	return;
}

origin_t& 
sequence_t::origin(void)
{
	return m_origin;
}

std::size_t 
sequence_t::origin_length(void) const
{
	std::size_t total(0);

	for (std::size_t idx = 0; idx < m_origin.count; idx++) {
		std::size_t length(m_origin.origin[idx].mer_length);
		std::size_t count(m_origin.origin[idx].repeat_cnt);

		if (length > std::numeric_limits< std::size_t >::max() / count)
			throw std::runtime_error("sequence_t::origin_length(): multiplicative overflow");
		if (total > std::numeric_limits< std::size_t >::max() - (length * count))
			throw std::runtime_error("sequence_t::origin_length(): additive overflow");

		total += length * count;
	}

	return total;
}

bool 
sequence_t::origin_ptr(basepair_t** dst, std::size_t& len) const
{
	basepair_t* tmp(nullptr);
	std::size_t total(origin_length());
	std::size_t offset(0);

	if (nullptr == dst || nullptr != *dst)
		return false;

	tmp = new basepair_t[total];
	len = total;

	/* we know there cannot be an overflow here due to the call to origin_length() */
	for (std::size_t idx = 0; idx < m_origin.count; idx++) {
		for (std::size_t cnt = 0; cnt < m_origin.origin[idx].repeat_cnt; cnt++) {
			std::memcpy(&tmp[offset], m_origin.origin[idx].ptr, m_origin.origin[idx].mer_length);
			offset += m_origin.origin[idx].mer_length;
		}
	}

	*dst = tmp;
	return true;
}

bool 
sequence_t::find_n_consecutive_repeats(basepair_t* start, basepair_t* end, const seq_element_t& seq, basepair_t** dst, basepair_t** cur) const
{
	std::size_t			repeats(seq.repeat_cnt);
	std::size_t			length(seq.mer_length);
	const long double	aratio(seq.at_ratio);
	const long double	gratio(seq.gc_ratio);
	basepair_t*			ptr(nullptr);
	basepair_t*			curr(nullptr == cur || nullptr == *cur ? start : *cur);

	if (nullptr != dst)
		*dst = nullptr;

	while (curr < (end - length)) {
		for (std::size_t idx = 0; idx < repeats; idx++) {
			const long double at(at_ratio(curr, curr + length));
			const long double gc(gc_ratio(curr, curr + length));

			if (at != aratio || gc != gratio) {
				curr++;
				ptr = nullptr;
				break;
			} else {
				if (0 == idx)
					ptr = curr;

				curr += length;
			}
		}

		if (nullptr != ptr) {
			if (nullptr != dst)
				*dst = ptr;
			if (nullptr != cur)
				*cur = curr;

			return true;
		}
	}

	return false;
}

bool
sequence_t::find_data(basepair_t* start, basepair_t* end, basepair_t** origin) const
{
	std::size_t length(origin_length());

	if (false == find_origin(start, end, origin))
		return false;

	if (nullptr != origin)
		*origin = (*origin + length); 

	return true;
}

bool 
sequence_t::find_origin(basepair_t* start, basepair_t* end, basepair_t** origin) const
{
	basepair_t* ptr(nullptr);
	basepair_t* curr(start);
	basepair_t* tmp(nullptr);

	if (nullptr != origin)
		*origin = nullptr;

	for (std::size_t idx = 0; idx < m_origin.count; ) {
		if (0 == idx) {
			if (true == find_n_consecutive_repeats(curr, end, m_origin.origin[idx], &ptr, &tmp)) {
				idx++;
				curr = tmp;
				continue;
			}

			curr++;
		} else {
			if (true == find_n_consecutive_repeats(curr, end, m_origin.origin[idx], nullptr, &tmp)) {
				idx++;
				curr = tmp;
				continue;
			} else {
				ptr = nullptr;
				idx = 0;
			}
		}

	}

	if (nullptr != origin)
		*origin = ptr;

	return ptr != nullptr;
}
