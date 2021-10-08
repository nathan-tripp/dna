#pragma once
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <limits>
#include <vector>
#include "global.hpp"
#include "sequence.hpp"
#include "rand.hpp"

class dna_t
{
	private:
		static const basepair_t* g_telomere;
		static const std::size_t g_telomere_len;

	protected:
		std::vector< basepair_t >	m_vec;
		basepair_t*					m_three;
		mutable basepair_t*			m_five;
		const std::size_t			m_length;
		std::size_t					m_repeats;
		sequence_t					m_sequence;

		static void set_byte(basepair_t* ptr, uint8_t val);
		static uint8_t get_byte(basepair_t* ptr, const std::size_t mlen = 4, bool left = true);
		static void set_output(uint8_t* dst, basepair_t* src, const std::size_t len, const bool is_odd = false, const std::size_t rem = 0);
		static void helix_buffer(basepair_t* ptr, const std::size_t len);

		static bool is_telomere_right(basepair_t* ptr);
		static bool is_telomere_left(basepair_t* ptr);

		virtual basepair_t* find_end_telomere(void);
		virtual basepair_t* find_start_of_data(void) const;
		virtual void telomere_fill(const std::size_t offset);
		virtual basepair_t* find_first_telomere(basepair_t* ptr, const std::size_t len) const;

		inline std::size_t length(void) { if (m_five < m_three) return 0; return m_five - m_three; }


		const std::size_t calculate_length(const std::size_t len, const std::size_t rep_cnt, const std::size_t max_origin);
		dna_t(const uint8_t* ptr, const std::size_t len, const std::size_t rep_cnt = 0, const std::size_t max_origin = 0);
		
	public:
		static dna_t* 
		new_dna(const uint8_t* ptr, const std::size_t len, const std::size_t rep_cnt = 0, const std::size_t max_origin = 0)
		{
			if (nullptr == ptr || 0 == len)
				return nullptr;

			return new dna_t(ptr, len, rep_cnt, max_origin);
		}

		virtual ~dna_t(void);

		virtual void set(const uint8_t* ptr, const std::size_t len);
		virtual uint8_t* get(std::size_t& out_len);
};
