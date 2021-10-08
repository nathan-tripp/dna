#include "dna.hpp"

void
dna_t::set_byte(basepair_t* ptr, uint8_t val)
{
	for (std::size_t idx = 0; idx < 4; idx++) {
		uint8_t v((val & 0xC0) >> 6);

		SET_NUCLEO_LEFT(ptr[idx], v & 3);
		SET_NUCLEO_RIGHT(ptr[idx], ~v & 3);
		val <<= 2;
	}

	return;
}

uint8_t
dna_t::get_byte(basepair_t* ptr, const std::size_t mlen, bool left)
{
	uint8_t v(0);

	for (std::size_t idx = 0; idx < 4; idx++) {
		if (true == left) {
			v |= NUCLEO_LEFT_SPLIT(ptr[idx]); //NUCLEO_LEFT(ptr[idx]);
			NUCLEO_CLEAR_LEFT_SPLIT(ptr[idx]);
		} else
			v |= ~(NUCLEO_RIGHT(ptr[idx])) & 0x03;

		if (3 != idx)
			v <<= 2;
	}

	return v;
}

const std::size_t 
dna_t::calculate_length(const std::size_t len, const std::size_t rep_cnt, const std::size_t max_origin)
{
	std::size_t cnt(0), ret(0);
	rand_t		rnd;

	cnt = len;

	if (0 != rep_cnt) {
		if (rep_cnt % dna_t::g_telomere_len)
			cnt += 1;

		if (rep_cnt > (std::numeric_limits< std::size_t >::max() - cnt) / dna_t::g_telomere_len)
			throw std::invalid_argument("dna_t::dna_t(): invalid telomere count (multiplicative overflow)");

		cnt += rep_cnt * dna_t::g_telomere_len;
	}
	else
		cnt += rnd.word() * dna_t::g_telomere_len;

	if (0 != max_origin) {
		if (cnt > std::numeric_limits< std::size_t >::max() - max_origin)
			throw std::invalid_argument("dna_t::dna_t(): invalid origin length (additive overflow)");

		cnt += max_origin;
	}
	else {
		std::size_t olen(rnd.word());

		if (cnt > std::numeric_limits< std::size_t >::max() - olen)
			throw std::invalid_argument("dna_t::dna_t(): invalid origin length (additive overflow)");

		cnt += olen;
	}

	do {
		ret = rnd.qword() % cnt;
	} while (ret < (cnt / 2));

	if (ret > std::numeric_limits< std::size_t >::max() - cnt)
		throw std::invalid_argument("dna_t::dna_t(): invalid padding length (additive overflow)");

	ret += cnt;
	return ret;
}

dna_t::dna_t(const uint8_t* ptr, const std::size_t len, const std::size_t rep_cnt, const std::size_t max_origin)
	: m_three(nullptr), m_five(nullptr), m_length(calculate_length(len, rep_cnt, max_origin)), m_repeats(rep_cnt)
{
	rand_t		rnd;
	std::size_t offset(rnd.qword() % (m_length / 3));
	basepair_t* origin(nullptr);
	std::size_t origin_len(0);

	m_vec.resize(m_length);
	m_three = &m_vec[0]; // new basepair_t[m_length];
	m_five	= m_three + m_length - 1;

	if (m_length > std::numeric_limits< std::size_t >::max() / sizeof(basepair_t))
		throw std::invalid_argument("dna_t::dna_t(): invalid parameter length");

	std::memset(m_three, 0, m_length * sizeof(basepair_t));

	if (false == m_sequence.origin_ptr(&origin, origin_len))
		throw std::invalid_argument("dna_t::dna_t(): failure during retrieval of origin data");

	std::memcpy(&m_three[offset], origin, origin_len);
	set(ptr, len);

	return;
}

dna_t::~dna_t(void)
{
	delete[] m_three;
	m_three = m_five = nullptr;
	return;
}

// TTAGGG
const basepair_t* dna_t::g_telomere = reinterpret_cast< const basepair_t* >("\x30\x30\x0c\x18\x18\x18"); 
const std::size_t dna_t::g_telomere_len = 6;

void 
dna_t::helix_buffer(basepair_t* ptr, const std::size_t len)
{
	std::size_t end(0);

	for (std::size_t idx = 0; idx < len - (dna_t::g_telomere_len*2); idx += dna_t::g_telomere_len) {
		SPLIT_NUCLEO_RIGHT(ptr[idx + dna_t::g_telomere_len + 0]);
		SPLIT_NUCLEO_RIGHT(ptr[idx + 0]);
		SET_NUCLEO_RIGHT(ptr[idx + dna_t::g_telomere_len + 0], NUCLEO_RIGHT_SPLIT(ptr[idx + 0]));
		NUCLEO_CLEAR_RIGHT_SPLIT(ptr[idx + 0]);

		SPLIT_NUCLEO_RIGHT(ptr[idx + dna_t::g_telomere_len + 1]);
		SPLIT_NUCLEO_RIGHT(ptr[idx + 1]);
		SET_NUCLEO_RIGHT(ptr[idx + dna_t::g_telomere_len + 1], NUCLEO_RIGHT_SPLIT(ptr[idx + 1]));
		NUCLEO_CLEAR_RIGHT_SPLIT(ptr[idx + 1]);

		SPLIT_NUCLEO_RIGHT(ptr[idx + dna_t::g_telomere_len + 2]);
		SPLIT_NUCLEO_RIGHT(ptr[idx + 2]);
		SET_NUCLEO_RIGHT(ptr[idx + dna_t::g_telomere_len + 2], NUCLEO_RIGHT_SPLIT(ptr[idx + 2]));
		NUCLEO_CLEAR_RIGHT_SPLIT(ptr[idx + 2]);

		SPLIT_NUCLEO_RIGHT(ptr[idx + dna_t::g_telomere_len + 3]);
		SPLIT_NUCLEO_RIGHT(ptr[idx + 3]);
		SET_NUCLEO_RIGHT(ptr[idx + dna_t::g_telomere_len + 3], NUCLEO_RIGHT_SPLIT(ptr[idx + 3]));
		NUCLEO_CLEAR_RIGHT_SPLIT(ptr[idx + 3]);

		SPLIT_NUCLEO_RIGHT(ptr[idx + dna_t::g_telomere_len + 4]);
		SPLIT_NUCLEO_RIGHT(ptr[idx + 4]);
		SET_NUCLEO_RIGHT(ptr[idx + dna_t::g_telomere_len + 4], NUCLEO_RIGHT_SPLIT(ptr[idx + 4]));
		NUCLEO_CLEAR_RIGHT_SPLIT(ptr[idx + 4]);

		SPLIT_NUCLEO_RIGHT(ptr[idx + dna_t::g_telomere_len + 5]);
		SPLIT_NUCLEO_RIGHT(ptr[idx + 5]);
		SET_NUCLEO_RIGHT(ptr[idx + dna_t::g_telomere_len + 5], NUCLEO_RIGHT_SPLIT(ptr[idx + 5]));
		NUCLEO_CLEAR_RIGHT_SPLIT(ptr[idx + 5]);
	}

	for (std::size_t idx = 0; idx < len; idx++)
		NUCLEO_CLEAR_RIGHT_SPLIT(ptr[idx]);

	return;
}

void 
dna_t::telomere_fill(const std::size_t offset)
{
	std::size_t len(m_five - m_three);
	std::size_t rem((len - offset) % dna_t::g_telomere_len);
	std::size_t tidx(0);
	std::size_t cnt(m_repeats);

	if (cnt) {
		for (std::size_t idx = offset; idx < len - rem; idx += dna_t::g_telomere_len) {
			std::memcpy(&m_three[idx], dna_t::g_telomere, dna_t::g_telomere_len);

			if (cnt)
				cnt--;
			else
				break;
		}

		cnt--;

		if (!cnt)
			return;

		for (std::size_t idx = len - rem; idx < len; idx++)
			m_three[idx] = dna_t::g_telomere[tidx++];
	} else {
		for (std::size_t idx = offset; idx < len - rem; idx += dna_t::g_telomere_len) 
			std::memcpy(&m_three[idx], dna_t::g_telomere, dna_t::g_telomere_len);

		for (std::size_t idx = len - rem; idx < len; idx++)
			m_three[idx] = dna_t::g_telomere[tidx++];
	}

	return;
}

void
dna_t::set(const uint8_t* ptr, const std::size_t len)
{
	std::size_t			length(len);
	std::size_t			lidx(0); //dna_t::g_telomere_len);
	std::size_t			oidx(0); //dna_t::g_telomere_len + (length * 4));
	basepair_t*			origin(find_start_of_data());

	//if (false == m_sequence.find_origin(m_three, m_five, &origin))
	if (nullptr == origin)
		throw std::invalid_argument("...");

	lidx = origin - m_three;
	oidx = (origin - m_three) + (length * 4);

	if (len > this->length())
		length = this->length();

	if (length > std::numeric_limits< std::size_t >::max() - (dna_t::g_telomere_len*2) || this->length() < dna_t::g_telomere_len * 2 + length)
			throw std::invalid_argument("dna_t::set(): invalid parameters (length)");

	//std::memcpy(m_three, dna_t::g_telomere, dna_t::g_telomere_len);

	for (std::size_t idx = 0; idx < length; idx++, lidx += 4)
		set_byte(&m_three[lidx], ptr[idx]);

//	for (std::size_t idx = 0; idx < length; idx++, lidx += 4)
//		set_byte(&m_three[lidx], ptr[idx]);

	telomere_fill(oidx);
	dna_t::helix_buffer(m_three, this->length());
	return;
}

bool
dna_t::is_telomere_right(basepair_t* ptr)
{
	for (std::size_t idx = 0; idx < dna_t::g_telomere_len; idx++)
		if (NUCLEO_RIGHT(ptr[idx]) != NUCLEO_RIGHT(dna_t::g_telomere[idx]))
			return false;

	return true;
}

bool
dna_t::is_telomere_left(basepair_t* ptr)
{
	for (std::size_t idx = 0; idx < dna_t::g_telomere_len; idx++)
		if (NUCLEO_LEFT(ptr[idx]) != NUCLEO_LEFT(dna_t::g_telomere[idx]))
			return false;

	return true;
}

basepair_t* 
dna_t::find_end_telomere(void)
{
	const uint8_t telo_right(NUCLEO_RIGHT(dna_t::g_telomere[dna_t::g_telomere_len - 1]));

	while (m_five > m_three) {
		while (*m_five != dna_t::g_telomere[dna_t::g_telomere_len - 1])
			m_five--;

		if (m_five - (dna_t::g_telomere_len*2) <= m_three)
			break;

		if (false == is_telomere_right(m_five - dna_t::g_telomere_len + 1)) {
			m_five--;
			continue;
		}

		if (false == is_telomere_left(m_five - (dna_t::g_telomere_len * 2) + 1)) {
			m_five--;
			continue;
		}

		m_five--;
		if (0 != m_repeats)
			m_repeats--;

		return m_five - (dna_t::g_telomere_len * 2) + 1;
	}

	return nullptr;

	for (basepair_t* end_ptr = m_five; end_ptr > m_three + dna_t::g_telomere_len; end_ptr--) {
		if (NUCLEO_RIGHT(*end_ptr) == telo_right) {
			if (end_ptr - m_three <= dna_t::g_telomere_len)
				continue;
			else if (true == is_telomere_right(end_ptr - dna_t::g_telomere_len + 1)) {
				if (end_ptr - m_three <= (dna_t::g_telomere_len * 2))
					continue;

				if (false == is_telomere_left(end_ptr - (dna_t::g_telomere_len * 2) + 1))
					continue;

				m_five -= 1;
				return end_ptr - (dna_t::g_telomere_len * 2) + 1;
			}
		}
	}

	return nullptr;
}

basepair_t* 
dna_t::find_start_of_data(void) const
{
	basepair_t* origin(nullptr);
	basepair_t* start(m_three);
	basepair_t* end(m_five);

	if (false == m_sequence.find_data(m_three, m_five, &origin))
		return nullptr;

	return origin;
}

basepair_t* 
dna_t::find_first_telomere(basepair_t* ptr, const std::size_t len) const
{
	for (basepair_t* start_ptr = ptr; start_ptr < ptr + len - dna_t::g_telomere_len; start_ptr++)
		if (0 == std::memcmp(start_ptr, dna_t::g_telomere, dna_t::g_telomere_len))
			return start_ptr;

	return nullptr;
}

void 
dna_t::set_output(uint8_t* dst, basepair_t* src, const std::size_t len, const bool is_odd, const std::size_t rem)
{
	if (true == is_odd) {
		for (std::size_t idx = 0; idx < len - 1; idx++)
			dst[idx] = get_byte(src + (idx * 4), 4);

		dst[len - 1] = get_byte(src + ((len - 1) * 4), rem);

	}
	else
		for (std::size_t idx = 0; idx < len; idx++)
			dst[idx] = get_byte(src + (idx * 4), 4);

	return;
}

uint8_t*
dna_t::get(std::size_t& out_len)
{
	basepair_t*		end_ptr(find_end_telomere());
	basepair_t*		start_ptr(nullptr); // find_start_of_data());
	basepair_t*		eod_ptr(nullptr); // find_first_telomere(start_ptr, (nullptr == start_ptr ? 0 : this->length() - (start_ptr - m_three) + dna_t::g_telomere_len)));
	basepair_t*		curr_ptr(start_ptr);
	std::size_t		len(0);
	std::size_t		rem(0);
	uint8_t*		ret(nullptr);

	out_len = 0;

	if (nullptr == end_ptr) // || nullptr == start_ptr || nullptr == eod_ptr)
		return nullptr;

	for (std::size_t off = 0; off < dna_t::g_telomere_len; off++) 
		SET_NUCLEO_RIGHT(end_ptr[off], NUCLEO_RIGHT(end_ptr[dna_t::g_telomere_len + off]));

	for (curr_ptr = m_three; curr_ptr < end_ptr; curr_ptr++) {
		if (curr_ptr + dna_t::g_telomere_len < m_five) {
//			SPLIT_NUCLEO_LEFT(*curr_ptr);
			SET_NUCLEO_RIGHT(*curr_ptr, NUCLEO_RIGHT(curr_ptr[dna_t::g_telomere_len]));
//			SET_NUCLEO_LEFT(*curr_ptr, (~(NUCLEO_RIGHT(*curr_ptr) & 0x03)));
		}
	}

	start_ptr	= find_start_of_data();
	eod_ptr = find_first_telomere(start_ptr, (nullptr == start_ptr ? 0 : this->length() - (start_ptr - m_three))); // +dna_t::g_telomere_len));
	len = eod_ptr - start_ptr; // -dna_t::g_telomere_len;

	if (nullptr == start_ptr || nullptr == eod_ptr)
		return nullptr;

	if (0 == len)
		return nullptr;
	else if (len % 4) {
		rem = len % 4;
		len = (len / 4) + 1;
	}
	else
		len = (len / 4);

	ret = new uint8_t[len + 1];
	std::memset(ret, 0, len + 1);

	for (curr_ptr = m_three; curr_ptr < end_ptr; curr_ptr++) {
		SPLIT_NUCLEO_LEFT(*curr_ptr);
		SET_NUCLEO_LEFT(*curr_ptr, (~(NUCLEO_RIGHT(*curr_ptr) & 0x03)));
	}

	set_output(ret, start_ptr, len, (end_ptr - start_ptr) % 4, rem);
	out_len = len;
	
	for (curr_ptr = m_three; curr_ptr < m_five; curr_ptr++)
		NUCLEO_CLEAR_LEFT_SPLIT(*curr_ptr);

	telomere_fill(eod_ptr - m_three);
	helix_buffer(m_three, this->length()); 
	return ret;
}
