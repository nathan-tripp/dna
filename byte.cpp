
/*
Actually I think this makes more sense
    pair[0] = (val & 0xC0)|(((~(val & 0xC0) >> 6) & 0x03) << 4)|(val & 0x0C)|((~(val & 0x0C) >> 2) & 0x03);
    pair[1] = (((~(val & 0x30) >> 4) & 3) << 6)|(val & 0x30)|(((~(val & 0x03) & 0x03) << 2))|(val & 0x03);


*/
uint8_t*
encode_byte(uint8_t val)
{
    uint8_t* ret(new uint8_t[4]);
//  uint8_t tmp(0); 

//  ret[0] = 0xF0;
//  ret[1] = 0x80;
//  ret[2] = 0x80;
//  ret[3] = 0x80;

    ret[0] = 0xF0 | (((val & 0xC0) >> 6) << 1);
    ret[1] = 0x80 | (val & 0x30) | (val & 0x0c);
    ret[2] = 0x80 | ((~((val & 0x0c) >> 2) & 0x03) << 4) | ((val & 0x03) << 2) | ((~val & 0x03));
    ret[3] = 0x80 | ((((~(val & 0xC0) >> 6)) & 0x03) << 4) | ((~((val & 0x30) >> 4) & 0x03) << 2);

/*  for (std::size_t idx = 0; idx < 4; idx++) {
        tmp = ((val & 0xC0) >> 6);
        val <<= 2;

        switch (idx) {
            case 0:
                ret[0] |= (tmp << 1);
                ret[3] |= ((~tmp & 0x03) << 4);
                break;
            case 1:
                ret[1] |= (tmp << 4);
                ret[3] |= ((~tmp & 0x03) << 2);
                break;
            case 2:
                ret[1] |= (tmp << 2);
                ret[2] |= (~tmp & 0x03) << 4;
                break;
            default:
                ret[2] |= (tmp << 2);
                ret[2] |= (~tmp & 0x03);
                break;
        }
    }
*/
    return ret;
}

uint8_t
decode_byte(uint8_t* ptr)
{
    uint8_t left(0), right(0);

    left    |= (((ptr[0] & 0x06) >> 1) << 6) | (((ptr[1] & 0x30) >> 4) << 4) | (((ptr[1] & 0x0C) >> 2) << 2) | (((ptr[2] & 0x0C) >> 2));
    right   |= (((ptr[3] & 0x30) >> 4) << 6) | (((ptr[3] & 0x0C) >> 2) << 4) | (((ptr[2] & 0x30) >> 4) << 2) | (ptr[2] & 0x03);

    if (0xFF != (left^right))
        throw std::runtime_error("decode_byte(): AT/GC-ratio mismatch");

    return left;
}
