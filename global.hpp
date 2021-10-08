#pragma once
#include <cstdint>

#define ADENINE 0
#define GUANINE 1
#define CYTOSINE 2
#define THYMINE 3

typedef uint8_t basepair_t;

#define NUCLEO_COMPLEMENT(x) ((~x) & 0x03)

#define NUCLEO_LEFT(x) ((x & 0x30) >> 4)
#define NUCLEO_LEFT_SPLIT(x) ((x & 0xC0) >> 6)
#define NUCLEO_CLEAR_LEFT(X) (x &= 0xCF)
#define NUCLEO_CLEAR_LEFT_SPLIT(x) (x &= 0x3E)
#define SET_NUCLEO_LEFT(x,y) do {			\
	x &= 0xCF;								\
	(x |= ((y & 0x03) << 4));				\
} while(0)

#define SPLIT_NUCLEO_LEFT(x) do {			\
	if (! (x & 0x01))						\
		(x |= (0x01 | ((x << 2) & 0xC0)));	\
} while(0)

#define NUCLEO_RIGHT(x) ((x & 0x0C) >> 2)
#define NUCLEO_RIGHT_SPLIT(x) (x & 0x03)
#define NUCLEO_CLEAR_RIGHT(x) (x &= 0xF3)
#define NUCLEO_CLEAR_RIGHT_SPLIT(x) (x &= 0x7C)

#define SET_NUCLEO_RIGHT(x,y) do {			\
	x &= 0xF3;								\
	(x |= ((y & 0x03) << 2));				\
} while(0) 

#define SPLIT_NUCLEO_RIGHT(x) do {			\
	if (! (x & 0x80))						\
		(x |= (0x80|(x >> 0x02) & 0x03));	\
} while(0)

#define SPLIT_SET_NUCLEO_RIGHT(x, y) do {	\
	if (! x & 0x80)							\
		x |= ((x & 0x0C) >> 2) & 0x03;		\
											\
	x &= 0xF3;								\
	 (x |= ((y & 0x03) << 2));				\
} while (0)
