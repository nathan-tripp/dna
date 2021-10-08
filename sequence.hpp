#pragma once
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include "global.hpp"
#include "rand.hpp"

typedef struct {
	std::size_t mer_length;
	std::size_t repeat_cnt;
	long double	at_ratio;
	long double	gc_ratio;
	basepair_t* ptr;
} seq_element_t;

typedef struct {
	std::size_t		count;
	seq_element_t*	origin;
} origin_t;

class sequence_t
{
	private:
	protected:
		origin_t	m_origin;
		rand_t		m_rnd;

		static long double at_ratio(const basepair_t* start, const basepair_t* end);
		static long double gc_ratio(const basepair_t* start, const basepair_t* end);

		static void count_nucleotides(const uint8_t, std::size_t&, std::size_t&, std::size_t&, std::size_t&);

		static bool is_at_rich(const basepair_t* start, const basepair_t* end);
		static bool is_gc_rich(const basepair_t* start, const basepair_t* end);

		inline basepair_t* create_sequence(seq_element_t&);

		inline std::size_t
		get_random(const std::size_t min, const std::size_t max)
		{
			std::size_t val(0);

			do {
				val = m_rnd.qword() % max;
			} while (min > val);

			return val;
		}

		virtual bool find_n_consecutive_repeats(basepair_t* start, basepair_t* end, const seq_element_t& seq, basepair_t** dst = nullptr, basepair_t** cur = nullptr) const;

	public:
		sequence_t(const std::size_t max_length = 0);
		virtual ~sequence_t(void);
		virtual origin_t& origin(void);
		virtual bool origin_ptr(basepair_t** dst, std::size_t& len) const;
		virtual std::size_t origin_length(void) const;
		virtual bool find_origin(basepair_t* start, basepair_t* end, basepair_t** origin) const;
		virtual bool find_data(basepair_t* start, basepair_t* end, basepair_t** origin) const;
};
