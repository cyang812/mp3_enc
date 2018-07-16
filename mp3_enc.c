/***************************************************************************
 *             __________               __   ___.
 *   Open      \______   \ ____   ____ |  | _\_ |__   _______  ___
 *   Source     |       _//  _ \_/ ___\|  |/ /| __ \ /  _ \  \/  /
 *   Jukebox    |    |   (  <_> )  \___|    < | \_\ (  <_> > <  <
 *   Firmware   |____|_  /\____/ \___  >__|_ \|___  /\____/__/\_ \
 *                     \/            \/     \/    \/            \/
 * $Id$
 *
 * Copyright (C) 2006 Antonius Hellmann
 * Copyright (C) 2006-2013 Michael Sevakis
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed on an "AS IS" basis, WITHOUT WARRANTY OF ANY
 * KIND, either express or implied.
 *
 ****************************************************************************/

//    Shine is an MP3 encoder
//    Copyright (C) 1999-2000  Gabriel Bouvigne
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Library General Public
//    License as published by the Free Software Foundation; either
//    version 2 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Library General Public License for more details.

#include "mp3_enc.h"
#include "enc_rom.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

int enc_stream_write(uint8_t *data, uint32_t size);
int enc_stream_lseek(uint32_t x, uint32_t y);
uint32_t enc_pcmbuf_read(uint8_t *buf, uint32_t size);
int enc_pcmbuf_advance(void);

static config_t  cfg                        ;
static short     mfbuf       [2*(1152+512)]
    /* for memcpy and 32-bit access */ ; /* 6656 Bytes */
static int       sb_data     [2][2][18][32] ; /* 9216 Bytes */
static int       mdct_freq   [576]          ; /* 2304 Bytes */
static char      mdct_sign   [576]          ; /*  576 Bytes */
static short     enc_data    [576]          ; /* 1152 Bytes */
static uint32_t  scalefac    [22]           ; /*   88 Bytes */
static bf_data   coded_data                 ; /* 1448 Bytes */
static uint8_t   band_scale_f[22]           ; /*   22 Bytes */

static void putbits(uint32_t val, uint32_t nbit)
{
    int new_bitpos = coded_data.bitpos + nbit;
    int ptrpos     = coded_data.bitpos >> 5;

    val = val & (0xffffffff >> (32 - nbit));

    /* data fit in one uint32_t */
    if (((new_bitpos - 1) >> 5) == ptrpos)
    {
        coded_data.bbuf[ptrpos] |= val << ((32 - new_bitpos) & 31);
    }
    else
    {
        coded_data.bbuf[ptrpos  ] |= val >> ((new_bitpos - 32) & 31);
        coded_data.bbuf[ptrpos+1] |= val << ((32 - new_bitpos) & 31);
    }

    coded_data.bitpos = new_bitpos;
}

#define PUTLONG_INIT() \
    uint32_t _cc = 0, _sz = 0

#define putlong_flush() \
    putbits(_cc, _sz)

#if 0
#define putlong(c, s)  \
    ({ uint32_t _c = (c), _s = (s);                                   \
       if (_s + _sz <= 32) { _cc = (_cc << _s) | _c; _sz += _s;     } \
       else                { putlong_flush(); _cc = _c; _sz = _s; } })
#else 1
#define putlong(c, s)  \
    { uint32_t _c = (c), _s = (s);                                   \
       if (_s + _sz <= 32) { _cc = (_cc << _s) | _c; _sz += _s;     } \
       else                { putlong_flush(); _cc = _c; _sz = _s; } }
#endif

static inline uint32_t encode_header(int padding, long bitr_id)
{
    /*
     * MPEG header layout:
     * AAAAAAAA AAABBCCD EEEEFFGH IIJJKLMM
     * A (31-21) = frame sync
     * B (20-19) = MPEG type
     * C (18-17) = MPEG layer
     * D (16)    = protection bit
     * E (15-12) = bitrate index
     * F (11-10) = samplerate index
     * G (9)     = padding bit
     * H (8)     = private bit
     * I (7-6)   = channel mode
     * J (5-4)   = mode extension (jstereo only)
     * K (3)     = copyright bit
     * L (2)     = original
     * M (1-0)   = emphasis
     */
    return (0xffe00000           )  /* frame sync (AAAAAAAAA AAA) */
         | (0x2             << 19)  /* mp3 type (upper):  1 (BB)  */
         | (cfg.mpg.type    << 19)
         | (0x1             << 17)  /* mp3 layer:        01 (CC)  */
         | (0x1             << 16)  /* mp3 crc:           1 (D)   */
         | (bitr_id         << 12)
         | (cfg.mpg.smpl_id << 10)
         | (padding         <<  9)
         | (cfg.mpg.mode    <<  6)
         | (0x1             <<  2); /* mp3 org:           1 (L)   */
    /* no emphasis (bits 0-1) */
}

static long calc_frame_size(int bitr_id, long *frac)
{
    unsigned long v = bitr_index[cfg.mpg.type][bitr_id];
    v = 576 * 16000 * v / (2 - cfg.mpg.type);
    v /= cfg.samplerate;

    if (frac)
        *frac = v % 64;

    return v / 64;
}

static void encode_side_info(side_info_t si[2][2])
{
    PUTLONG_INIT();

    putbits(encode_header(cfg.mpg.padding, cfg.mpg.bitr_id), 32);

    if (cfg.mpg.type == 1)
    {
        /* MPEG1 */
        #if 0
        if (cfg.channels == 2) putlong(0, 20);
        else                   putlong(0, 18);
        #else
        if (cfg.channels == 2) putlong(0, 20)
        else                   putlong(0, 18)
        #endif
        for (int gr = 0; gr < cfg.granules; gr++)
        {
            for (int ch = 0; ch < cfg.channels; ch++)
            {
                side_info_t *gi = &si[gr][ch];

                putlong(gi->part2_3_length + 42, 12); /* add scale_facs array size */
                putlong(gi->address3 >> 1,        9);
                putlong(gi->global_gain,          8);
                putlong(9,                        4); /* set scale_facs compr type */
                putlong(gi->table_select[0],      6);
                putlong(gi->table_select[1],      5);
                putlong(gi->table_select[2],      5);
                putlong(gi->region_0_1,           7);
                putlong(1                  ,      2); /* set scale_facs to 1bit */
                putlong(gi->table_select[3],      1);
            }
        }
    }
    else
    {
        /* MPEG2 */
        #if 0
        if (cfg.channels == 2) putlong(0, 10);
        else                   putlong(0,  9);
        #else
        if (cfg.channels == 2) putlong(0, 10)
        else                   putlong(0,  9)
        #endif

        for (int ch = 0; ch < cfg.channels; ch++)
        {
            side_info_t *gi = &si[0][ch];

            putlong(gi->part2_3_length + 42, 12); /* add scale_facs array size */
            putlong(gi->address3 >> 1,        9);
            putlong(gi->global_gain,          8);
            putlong(0xca,                     9); /* set scale_facs compr type */
            putlong(gi->table_select[0],      6);
            putlong(gi->table_select[1],      5);
            putlong(gi->table_select[2],      5);
            putlong(gi->region_0_1     ,      7);
            putlong(1                  ,      1); /* set scale_facs to 1bit */
            putlong(gi->table_select[3],      1);
        }
    }

    putlong_flush();
}

/* Implements the pseudocode of page 98 of the IS */
static int huffman_code(short *ix, char *xr_sign, uint32_t begin, uint32_t end,
                        int table)
{
    if (table == 0)
        return 0;

    PUTLONG_INIT();
    int sumbit = 0;

    #define sign_x xr_sign[i+0]
    #define sign_y xr_sign[i+1]

    if (table > 15)
    {
        /* ESC-table is used */
        uint32_t linbits = ht_big[table - 16].linbits;
        const uint16_t *hffcode = table < 24 ? t16HB : t24HB;
        const uint8_t  *hlen    = table < 24 ? t16l  : t24l;
        uint32_t xl = 0, yl = 0;

        for (uint32_t i = begin; i < end; i += 2)
        {
            int x = ix[i+0];
            int y = ix[i+1];

            if (x > 14) { xl = x - 15;  x = 15; }
            if (y > 14) { yl = y - 15;  y = 15; }

            uint32_t idx  = x * 16 + y;
            uint32_t code = hffcode[idx];
            int      bit  = hlen   [idx];

            if (x)
            {
                if (x > 14)
                {
                    code = (code << linbits) | xl;
                    bit += linbits;
                }

                code = (code << 1) | sign_x;
                bit += 1;
            }

            if (y)
            {
                if (y > 14)
                {
                    if (bit + linbits + 1 > 32)
                    {
                        putlong(code, bit);
                        sumbit += bit;
                        code = bit = 0;
                    }

                    code = (code << linbits) | yl;
                    bit += linbits;
                }

                code = (code << 1) | sign_y;
                bit += 1;
            }

            putlong(code, bit);
            sumbit += bit;
        }
    }
    else
    {
        /* No ESC-words */
        const struct huffcodetab *h = &ht[table];

        for (uint32_t i = begin; i < end; i += 2)
        {
            int x = ix[i+0];
            int y = ix[i+1];

            uint32_t idx  = x * h->len + y;
            uint32_t code = h->table[idx];
            int      bit  = h->hlen [idx];

            if (x)
            {
                code = (code << 1) | sign_x;
                bit += 1;
            }

            if (y)
            {
                code = (code << 1) | sign_y;
                bit += 1;
            }

            putlong(code, bit);
            sumbit += bit;
        }
    }

    putlong_flush();

    return sumbit;
}

static int huffman_cod1(short *ix, char *xr_sign, uint32_t begin,
                        uint32_t end, int tbl)
{
    PUTLONG_INIT();

    #define sgnv xr_sign[i+0]
    #define sgnw xr_sign[i+1]
    #define sgnx xr_sign[i+2]
    #define sgny xr_sign[i+3]

    int sumbit = 0, s = 0, l = 0;
    for (uint32_t i = begin; i < end; i += 4)
    {
        int v = ix[i+0];
        int w = ix[i+1];
        int x = ix[i+2];
        int y = ix[i+3];

        uint32_t p = (v << 3) + (w << 2) + (x << 1) + y;

        switch (p)
        {
        case  0: l=0; s = 0; break;
        case  1: l=1; s =                                           sgny; break;
        case  2: l=1; s =                              sgnx;              break;
        case  3: l=2; s =                             (sgnx << 1) + sgny; break;
        case  4: l=1; s =                sgnw;                            break;
        case  5: l=2; s =               (sgnw << 1)               + sgny; break;
        case  6: l=2; s =               (sgnw << 1) +  sgnx;              break;
        case  7: l=3; s =               (sgnw << 2) + (sgnx << 1) + sgny; break;
        case  8: l=1; s =  sgnv;                                          break;
        case  9: l=2; s = (sgnv << 1)                             + sgny; break;
        case 10: l=2; s = (sgnv << 1)               +  sgnx;              break;
        case 11: l=3; s = (sgnv << 2)               + (sgnx << 1) + sgny; break;
        case 12: l=2; s = (sgnv << 1) +  sgnw;                            break;
        case 13: l=3; s = (sgnv << 2) + (sgnw << 1)               + sgny; break;
        case 14: l=3; s = (sgnv << 2) + (sgnw << 1) +  sgnx;              break;
        case 15: l=4; s = (sgnv << 3) + (sgnw << 2) + (sgnx << 1) + sgny; break;
        }

        uint32_t d = (ht_count[tbl][0][p] << l) + s;
        l = ht_count[tbl][1][p];
        putlong(d, l);
        sumbit += l;
    }

    putlong_flush();

    return sumbit;
}

/* Note the discussion of huffmancodebits() on pages 28 and 29 of the IS,
   as well as the definitions of the side information on pages 26 and 27. */
static void huffman_code_bits(short *ix, char *xr_sign, side_info_t *gi)
{
    int region1 = gi->address1;
    int region2 = gi->address2;
    int bigvals = gi->address3;
    int count1  = bigvals + (gi->count1 << 2);
    int i, v;

    for (i = v = 0; i < 32; i += 2)
        v |= band_scale_f[i >> 1] << (30 - i);

    putbits(v, 32); // store scale_facs (part1)

    for (v = 0; i < 42; i += 2)
        v |= band_scale_f[i >> 1] << (40 - i);

    putbits(v, 10); // store scale_facs (part2)

    int bits = 0;

    if (region1 > 0)
        bits += huffman_code(ix, xr_sign,       0, region1, gi->table_select[0]);

    if (region2 > region1)
        bits += huffman_code(ix, xr_sign, region1, region2, gi->table_select[1]);

    if (bigvals > region2)
        bits += huffman_code(ix, xr_sign, region2, bigvals, gi->table_select[2]);

    if (count1 > bigvals)
        bits += huffman_cod1(ix, xr_sign, bigvals,  count1, gi->table_select[3]);

    int stuff_bits = gi->part2_3_length - bits;

    if (stuff_bits > 0)
    {
        int stuff_words = stuff_bits >> 5;
        int remain_bits = stuff_bits & 31;

        if (remain_bits)
            putbits(~0, remain_bits);

        while (stuff_words--)
            putbits(~0, 32); /* Huffman code tables leed to padding ones */
    }
}

/*************************************************************************/
/* Function: Count the number of bits necessary to code the subregion.   */
/*************************************************************************/
static int count_bit1(short *ix, uint32_t start, uint32_t end, int *bits)
{
    int sum = 0;

    for (uint32_t i = start; i < end; i += 2)
        sum += t1l[4 + ix[i] * 2 + ix[i + 1]];

    *bits = sum;

    return 1; /* this is table1 */
}

static int count_bigv(short *ix, uint32_t start, uint32_t end, int table0,
                      int table1, int *bits)
{
    uint32_t sum = 0, bigv = 0;

    /* ESC-table is used */
    for (uint32_t i = start; i < end; i += 2)
    {
        uint32_t x = ix[i+0];
        uint32_t y = ix[i+1];

        if (x > 14) { x = 15; bigv++; }
        if (y > 14) { y = 15; bigv++; }

        sum += tab1624[x * 16 + y];
    }

    int sum0 = (sum  >>  16)  + bigv * ht_big[table0].linbits;
    int sum1 = (sum & 0xffff) + bigv * ht_big[table1].linbits;

    if (sum0 <= sum1)
    {
        *bits = sum0;
        return table0;
    }
    else
    {
        *bits = sum1;
        return table1;
    }
}

static int find_best_2(short *ix, uint32_t start, uint32_t end,
                       const uint32_t *table, uint32_t len, int *bits)
{
    uint32_t sum = 0;

    for (uint32_t i = start; i < end; i += 2)
        sum += table[ix[i] * len + ix[i + 1]];

    if ((sum & 0xffff) <= (sum >> 16))
    {
        *bits = sum & 0xffff;
        return 1;
    }
    else
    {
        *bits = sum >> 16;
        return 0;
    }
}

static int find_best_3(short *ix, uint32_t start, uint32_t end,
                       const uint32_t *table, uint32_t len, int *bits)
{
    int sum1 = 0;
    int sum2 = 0;
    int sum3 = 0;

    /* avoid overflow in packed additions: 78*13 < 1024 */
    for (uint32_t i = start; i < end;)
    {
        uint32_t j = i + 2*78 > end ? end : i + 2*78;

        uint32_t sum = 0;
        for (; i < j; i += 2)
            sum += table[ix[i] * len + ix[i + 1]];

        sum1 += (sum >> 20);
        sum2 += (sum >> 10) & 0x3ff;
        sum3 += (sum >>  0) & 0x3ff;
    }

    int r = 0;
    if (sum1 > sum2) { sum1 = sum2; r = 1; }
    if (sum1 > sum3) { sum1 = sum3; r = 2; }

    *bits = sum1;

    return r;
}

/***************************************************************************/
/*  Choose the Huffman table that will encode ix[begin..end] with          */
/*  the fewest bits.                                                       */
/*  Note: This code contains knowledge about the sizes and characteristic  */
/*  of the Huffman tables as defined in the IS (Table B.7), and will not   */
/*  work with any arbitrary tables.                                        */
/***************************************************************************/
static int choose_table(short *ix, uint32_t begin, uint32_t end, int *bits)
{
    int max = 0;
    for (uint32_t i = begin; i < end; i++)
    {
        if (ix[i] > max)
            max = ix[i];
    }

    if (max < 16)
    {
        /* tables without linbits */
        /* indx: 0  1  2  3  4  5  6  7  8  9 10 11 12  13 14  15 */
        /*  len: 0, 2, 3, 3, 0, 4, 4, 6, 6, 6, 8, 8, 8, 16, 0, 16 */
        switch (max)
        {
        case 0:  return  0;
        case 1:  return       count_bit1(ix, begin, end,              bits);
        case 2:  return  2 + find_best_2(ix, begin, end,   tab23,  3, bits);
        case 3:  return  5 + find_best_2(ix, begin, end,   tab56,  4, bits);
        case 4:
        case 5:  return  7 + find_best_3(ix, begin, end,  tab789,  6, bits);
        case 6:
        case 7:  return 10 + find_best_3(ix, begin, end,  tabABC,  8, bits);
        default: return 13 + find_best_2(ix, begin, end, tab1315, 16, bits)*2;
        }
    }
    else
    {
        /* tables with linbits */
        max -= 15;

        int table0, table1;

        for (table0 = 0; table0 < 8; table0++)
        {
            if (ht_big[table0].linmax >= max)
                break;
        }

        for (table1 = 8; table1 < 16; table1++)
        {
            if (ht_big[table1].linmax >= max)
                break;
        }

        return 16 + count_bigv(ix, begin, end, table0, table1, bits);
    }
}

/*************************************************************************/
/* Function: Calculation of rzero, count1, address3                      */
/* (Partitions ix into big values, quadruples and zeros).                */
/*************************************************************************/
static int calc_runlen(short *ix, side_info_t *si)
{
    int i = 576;

    while (i -= 2)
    {
        if (*(uint32_t *)&ix[i - 2]) /* !!!! short *ix; !!!!! */
            break;
    }

    si->count1 = 0;

    int sum = 0;
    for ( ; i > 3; i -= 4)
    {
        int v = ix[i-1];
        int w = ix[i-2];
        int x = ix[i-3];
        int y = ix[i-4];

        if ((v | w | x | y) <= 1)
        {
            int p = (y << 3) + (x << 2) + (w << 1) + (v);
            sum += tab01[p];
            si->count1++;
        }
        else
        {
            break;
        }
    }

    si->address3 = i;

    if ((sum >> 16) < (sum & 0xffff))
    {
        si->table_select[3] = 0;
        return sum >> 16;
    }
    else
    {
        si->table_select[3] = 1;
        return sum & 0xffff;
    }
}

/*************************************************************************/
/*   Function: Quantization of the vector xr ( -> ix)                    */
/*************************************************************************/
static int quantize_int(int *xr, short *ix, side_info_t *si)
{
    static const unsigned int frac_pow[] =
        { 0x10000, 0xd745, 0xb505, 0x9838 };

    unsigned s = frac_pow[si->quantStep & 3] >> si->quantStep / 4;

    /* check for possible 'out of range' values */
    if (((si->max_val + 256) >> 8) * s >= (65536 << 8))
        return 0;

    if (((si->max_val + 256) >> 8) * s < (4096 << 8))
    {
        /* all values fit the table size */
        for (int i = 576; i--; )
            ix[i] = int2idx[(xr[i] * s + 0x8000) >> 16];
    }
    else
    {
        /* check each index wether it fits the table */
        for (int i = 576; i--; )
        {
            unsigned idx = (xr[i] * s + 0x08000) >> 16;

            if (idx > 4095) ix[i] = int2idx[(idx + 8) >> 4] << 3;
            else            ix[i] = int2idx[idx];
        }
    }

    return 1;
}

/*************************************************************************/
/* subdivides the bigvalue region which will use separate Huffman tables */
/*************************************************************************/
static void subdivide(side_info_t *si)
{
    if (!si->address3)
    {
        /* no bigvalue region */
        si->region_0_1 = 0;
        si->address1   = 0;
        si->address2   = 0;
    }
    else
    {
        /* Calculate scale factor band index */
        int scfb = 1;
        while (scalefac[scfb] < si->address3)
            scfb++;

        int count0 = subdv_table[scfb].region0_cnt;
        int count1 = subdv_table[scfb].region1_cnt;

        si->region_0_1 = (count0 << 3) | count1;
        si->address1   = scalefac[count0];
        si->address2   = scalefac[count0 + count1 + 1];
    }
}

/*******************************************************************/
/* Count the number of bits necessary to code the bigvalues region */
/*******************************************************************/
static int bigv_bitcount(short *ix, side_info_t *gi)
{
    int b1 = 0, b2 = 0, b3 = 0;

    /* Select huffman code tables for bigvalues regions */
    gi->table_select[0] = 0;
    gi->table_select[1] = 0;
    gi->table_select[2] = 0;

    if (gi->address1 > 0)            /* region0 */
        gi->table_select[0] = choose_table(ix,            0, gi->address1, &b1);

    if (gi->address2 > gi->address1) /* region1 */
        gi->table_select[1] = choose_table(ix, gi->address1, gi->address2, &b2);

    if (gi->address3 > gi->address2) /* region2 */
        gi->table_select[2] = choose_table(ix, gi->address2, gi->address3, &b3);

    return b1 + b2 + b3;
}

static int quantize_and_count_bits(int *xr, short *ix, side_info_t *si)
{
    int bits = 10000;

    if (quantize_int(xr, ix, si))
    {
        bits = calc_runlen(ix, si);      /* rzero,count1,address3  */
        subdivide(si);                   /* bigvalues sfb division */
        bits += bigv_bitcount(ix,si);    /* bit count */
    }

    return bits;
}

/************************************************************************/
/* The code selects the best quantStep for a particular set of scalefacs*/
/************************************************************************/
static int inner_loop(int *xr, int max_bits, side_info_t *si)
{
    int bits;

    while ((bits = quantize_and_count_bits(xr, enc_data, si)) < max_bits - 64)
    {
        if (si->quantStep == 0)
            break;

        if (si->quantStep <= 2)
            si->quantStep  = 0;
        else
            si->quantStep -= 2;
    }

    while (bits > max_bits)
    {
        si->quantStep++;
        bits = quantize_and_count_bits(xr, enc_data, si);
    }

    return bits;
}

static void iteration_loop(int *xr, side_info_t *si, int gr_cnt)
{
    int max_bits = cfg.mean_bits;

    /* distribute reserved bits to remaining granules */
    int tar_bits = max_bits + (cfg.resv_size / gr_cnt & ~7);

    if (tar_bits > max_bits + max_bits / 2)
        tar_bits = max_bits + max_bits / 2;

    si->part2_3_length = inner_loop(xr, tar_bits, si);
    si->global_gain    = si->quantStep + 142 - si->additStep;

    /* unused bits of the reservoir can be used for remaining granules */
    cfg.resv_size += max_bits - si->part2_3_length;

    /* end: distribute the reserved bits to one or two granules */
    if (gr_cnt == 1)
    {
        si->part2_3_length += cfg.resv_size;

        /* mp3 format allows max 12bits for granule length */
        if (si->part2_3_length > 4092)
        {
            int remain = (si->part2_3_length - 4092 + 31) >> 5;
            si->part2_3_length    -= remain << 5;
            si[-1].part2_3_length += remain << 5;

            while (remain--)
                putbits(~0, 32);
        }
    }
}

/* returns sum_j=0^31 a[j]*cos(PI*j*(k+1/2)/32), 0<=k<32 */
/* Generic CPU */
static void window_subband1_s_(const short *wk, int a[32])
{
    for (int k = 0; k < 18; k++, wk += 64, a += 32)
    {
        const short *wp = enwindow;
        const short *x1 = wk;
        const short *x2 = x1 - 124;
        int s, t;

        /* x1[-572] .... x1[448] = 1022 */
        /* 18*4*16*32 */
        for (int i = -15; i < 0; i++)
        {
            s  = (int)x2[-224*2] * wp[ 0];  t  = (int)x1[ 224*2] * wp[ 0];
            s += (int)x2[-160*2] * wp[ 1];  t += (int)x1[ 160*2] * wp[ 1];
            s += (int)x2[- 96*2] * wp[ 2];  t += (int)x1[  96*2] * wp[ 2];
            s += (int)x2[- 32*2] * wp[ 3];  t += (int)x1[  32*2] * wp[ 3];
            s += (int)x2[  32*2] * wp[ 4];  t += (int)x1[- 32*2] * wp[ 4];
            s += (int)x2[  96*2] * wp[ 5];  t += (int)x1[- 96*2] * wp[ 5];
            s += (int)x2[ 160*2] * wp[ 6];  t += (int)x1[-160*2] * wp[ 6];
            s += (int)x2[ 224*2] * wp[ 7];  t += (int)x1[-224*2] * wp[ 7];
            s += (int)x1[-256*2] * wp[ 8];  t -= (int)x2[ 256*2] * wp[ 8];
            s += (int)x1[-192*2] * wp[ 9];  t -= (int)x2[ 192*2] * wp[ 9];
            s += (int)x1[-128*2] * wp[10];  t -= (int)x2[ 128*2] * wp[10];
            s += (int)x1[- 64*2] * wp[11];  t -= (int)x2[  64*2] * wp[11];
            s += (int)x1[   0*2] * wp[12];  t -= (int)x2[   0*2] * wp[12];
            s += (int)x1[  64*2] * wp[13];  t -= (int)x2[- 64*2] * wp[13];
            s += (int)x1[ 128*2] * wp[14];  t -= (int)x2[-128*2] * wp[14];
            s += (int)x1[ 192*2] * wp[15];  t -= (int)x2[-192*2] * wp[15];

            a[30+i*2] =  shft4(t)          + shft13(s) * wp[16];
            a[31+i*2] = shft13(t) * wp[17] - shft13(s) * wp[18];
            wp += 20;
            x1 -=  2;
            x2 +=  2;
        }

        t  =  (int)x1[- 16*2]              * wp[ 8];  s  = (int)x1[ -32*2] * wp[0];
        t += ((int)x1[- 48*2] - x1[ 16*2]) * wp[ 9];  s += (int)x1[ -96*2] * wp[1];
        t += ((int)x1[- 80*2] + x1[ 48*2]) * wp[10];  s += (int)x1[-160*2] * wp[2];
        t += ((int)x1[-112*2] - x1[ 80*2]) * wp[11];  s += (int)x1[-224*2] * wp[3];
        t += ((int)x1[-144*2] + x1[112*2]) * wp[12];  s += (int)x1[  32*2] * wp[4];
        t += ((int)x1[-176*2] - x1[144*2]) * wp[13];  s += (int)x1[  96*2] * wp[5];
        t += ((int)x1[-208*2] + x1[176*2]) * wp[14];  s += (int)x1[ 160*2] * wp[6];
        t += ((int)x1[-240*2] - x1[208*2]) * wp[15];  s += (int)x1[ 224*2] * wp[7];

        int u = shft4(s - t);
        int v = shft4(s + t);
        t = a[14];
        s = a[15] - t;

        a[31] = v + t;   /* A0 */
        a[30] = u + s;   /* A1 */
        a[15] = u - s;   /* A2 */
        a[14] = v - t;   /* A3 */
    }
}

static inline void window_subband1_s(const short *wk, int a0[32], int a1[32])
{
    window_subband1_s_(wk    , a0);
    window_subband1_s_(wk + 1, a1);
}

static void window_subband1_m(const short *wk, int a[32])
{
    for (int k = 0; k < 18; k++, wk += 32, a += 32)
    {
        const short *wp = enwindow;
        const short *x1 = wk;
        const short *x2 = x1 - 62;
        int s, t;

        /* x1[-286] .... x1[224] = 511 */
        /* 18*2*16*32 */
        for (int i = -15; i < 0; i++)
        {
            s  = (int)x2[-224] * wp[ 0];  t  = (int)x1[ 224] * wp[ 0];
            s += (int)x2[-160] * wp[ 1];  t += (int)x1[ 160] * wp[ 1];
            s += (int)x2[- 96] * wp[ 2];  t += (int)x1[  96] * wp[ 2];
            s += (int)x2[- 32] * wp[ 3];  t += (int)x1[  32] * wp[ 3];
            s += (int)x2[  32] * wp[ 4];  t += (int)x1[- 32] * wp[ 4];
            s += (int)x2[  96] * wp[ 5];  t += (int)x1[- 96] * wp[ 5];
            s += (int)x2[ 160] * wp[ 6];  t += (int)x1[-160] * wp[ 6];
            s += (int)x2[ 224] * wp[ 7];  t += (int)x1[-224] * wp[ 7];
            s += (int)x1[-256] * wp[ 8];  t -= (int)x2[ 256] * wp[ 8];
            s += (int)x1[-192] * wp[ 9];  t -= (int)x2[ 192] * wp[ 9];
            s += (int)x1[-128] * wp[10];  t -= (int)x2[ 128] * wp[10];
            s += (int)x1[- 64] * wp[11];  t -= (int)x2[  64] * wp[11];
            s += (int)x1[   0] * wp[12];  t -= (int)x2[   0] * wp[12];
            s += (int)x1[  64] * wp[13];  t -= (int)x2[- 64] * wp[13];
            s += (int)x1[ 128] * wp[14];  t -= (int)x2[-128] * wp[14];
            s += (int)x1[ 192] * wp[15];  t -= (int)x2[-192] * wp[15];

            a[30+i*2] =  shft4(t)          + shft13(s) * wp[16];
            a[31+i*2] = shft13(t) * wp[17] - shft13(s) * wp[18];
            wp += 20;
            x1--;
            x2++;
        }

        t  =  (int)x1[- 16]            * wp[ 8];  s  = (int)x1[ -32] * wp[0];
        t += ((int)x1[- 48] - x1[ 16]) * wp[ 9];  s += (int)x1[ -96] * wp[1];
        t += ((int)x1[- 80] + x1[ 48]) * wp[10];  s += (int)x1[-160] * wp[2];
        t += ((int)x1[-112] - x1[ 80]) * wp[11];  s += (int)x1[-224] * wp[3];
        t += ((int)x1[-144] + x1[112]) * wp[12];  s += (int)x1[  32] * wp[4];
        t += ((int)x1[-176] - x1[144]) * wp[13];  s += (int)x1[  96] * wp[5];
        t += ((int)x1[-208] + x1[176]) * wp[14];  s += (int)x1[ 160] * wp[6];
        t += ((int)x1[-240] - x1[208]) * wp[15];  s += (int)x1[ 224] * wp[7];

        int u = shft4(s - t);
        int v = shft4(s + t);
        t = a[14];
        s = a[15] - t;

        a[31] = v + t;   /* A0 */
        a[30] = u + s;   /* A1 */
        a[15] = u - s;   /* A2 */
        a[14] = v - t;   /* A3 */
    }
}

static void window_subband2_(int a[32])
{
    /* 36864=4*18*16*32 */
    const short * const wp = enwindow + 20 * 15;
    for (int k = 0; k < 18; k++, a += 32)
    {
        int xr;

        xr = a[28] - a[0];  a[0] += a[28];  a[28] = shft9(xr) * wp[-2*20+17];
        xr = a[29] - a[1];  a[1] += a[29];  a[29] = shft9(xr) * wp[-2*20+17];
        xr = a[26] - a[2];  a[2] += a[26];  a[26] = shft9(xr) * wp[-4*20+17];
        xr = a[27] - a[3];  a[3] += a[27];  a[27] = shft9(xr) * wp[-4*20+17];
        xr = a[24] - a[4];  a[4] += a[24];  a[24] = shft9(xr) * wp[-6*20+17];
        xr = a[25] - a[5];  a[5] += a[25];  a[25] = shft9(xr) * wp[-6*20+17];
        xr = a[22] - a[6];  a[6] += a[22];  a[22] = shft9(xr) * SQRT        ;
        xr = a[23] - a[7];  a[7] += a[23];  a[23] = shft9(xr) * SQRT  - a[7];
        a[ 7] -= a[ 6];
        a[22] -= a[ 7];
        a[23] -= a[22];

        xr = a[ 6];  a[ 6] = a[31] - xr;  a[31] = a[31] + xr;
        xr = a[ 7];  a[ 7] = a[30] - xr;  a[30] = a[30] + xr;
        xr = a[22];  a[22] = a[15] - xr;  a[15] = a[15] + xr;
        xr = a[23];  a[23] = a[14] - xr;  a[14] = a[14] + xr;

        xr = a[20] - a[ 8];  a[ 8] += a[20];  a[20] = shft9(xr) * wp[-10*20+17];
        xr = a[21] - a[ 9];  a[ 9] += a[21];  a[21] = shft9(xr) * wp[-10*20+17];
        xr = a[18] - a[10];  a[10] += a[18];  a[18] = shft9(xr) * wp[-12*20+17];
        xr = a[19] - a[11];  a[11] += a[19];  a[19] = shft9(xr) * wp[-12*20+17];
        xr = a[16] - a[12];  a[12] += a[16];  a[16] = shft9(xr) * wp[-14*20+17];
        xr = a[17] - a[13];  a[13] += a[17];  a[17] = shft9(xr) * wp[-14*20+17];
        xr =-a[20] + a[24];  a[20] += a[24];  a[24] = shft9(xr) * wp[-12*20+17];
        xr =-a[21] + a[25];  a[21] += a[25];  a[25] = shft9(xr) * wp[-12*20+17];
        xr = a[ 4] - a[ 8];  a[ 4] += a[ 8];  a[ 8] = shft9(xr) * wp[-12*20+17];
        xr = a[ 5] - a[ 9];  a[ 5] += a[ 9];  a[ 9] = shft9(xr) * wp[-12*20+17];
        xr = a[ 0] - a[12];  a[ 0] += a[12];  a[12] = shft9(xr) * wp[ -4*20+17];
        xr = a[ 1] - a[13];  a[ 1] += a[13];  a[13] = shft9(xr) * wp[ -4*20+17];
        xr = a[16] - a[28];  a[16] += a[28];  a[28] = shft9(xr) * wp[ -4*20+17];
        xr =-a[17] + a[29];  a[17] += a[29];  a[29] = shft9(xr) * wp[ -4*20+17];

        xr = SQRT * shft9(a[ 2] - a[10]);  a[ 2] += a[10];  a[10] = xr;
        xr = SQRT * shft9(a[ 3] - a[11]);  a[ 3] += a[11];  a[11] = xr;
        xr = SQRT * shft9(a[26] - a[18]);  a[18] += a[26];  a[26] = xr - a[18];
        xr = SQRT * shft9(a[27] - a[19]);  a[19] += a[27];  a[27] = xr - a[19];

        xr = a[ 2];  a[19] -= a[ 3];  a[ 3] -= xr;  a[ 2] = a[31] - xr;  a[31] += xr;
        xr = a[ 3];  a[11] -= a[19];  a[18] -= xr;  a[ 3] = a[30] - xr;  a[30] += xr;
        xr = a[18];  a[27] -= a[11];  a[19] -= xr;  a[18] = a[15] - xr;  a[15] += xr;

        xr = a[19];  a[10] -= xr;  a[19] = a[14] - xr;  a[14] += xr;
        xr = a[10];  a[11] -= xr;  a[10] = a[23] - xr;  a[23] += xr;
        xr = a[11];  a[26] -= xr;  a[11] = a[22] - xr;  a[22] += xr;
        xr = a[26];  a[27] -= xr;  a[26] = a[ 7] - xr;  a[ 7] += xr;

        xr = a[27];  a[27] = a[6] - xr;  a[6] += xr;

        xr = SQRT * shft9(a[ 0] - a[ 4]);  a[ 0] += a[ 4];  a[ 4] = xr;
        xr = SQRT * shft9(a[ 1] - a[ 5]);  a[ 1] += a[ 5];  a[ 5] = xr;
        xr = SQRT * shft9(a[16] - a[20]);  a[16] += a[20];  a[20] = xr;
        xr = SQRT * shft9(a[17] - a[21]);  a[17] += a[21];  a[21] = xr;
        xr =-SQRT * shft9(a[ 8] - a[12]);  a[ 8] += a[12];  a[12] = xr - a[ 8];
        xr =-SQRT * shft9(a[ 9] - a[13]);  a[ 9] += a[13];  a[13] = xr - a[ 9];
        xr =-SQRT * shft9(a[25] - a[29]);  a[25] += a[29];  a[29] = xr - a[25];
        xr =-SQRT * shft9(a[24] + a[28]);  a[24] -= a[28];  a[28] = xr - a[24];

        xr = a[24] - a[16]; a[24] = xr;
        xr = a[20] - xr;    a[20] = xr;
        xr = a[28] - xr;    a[28] = xr;

        xr = a[25] - a[17]; a[25] = xr;
        xr = a[21] - xr;    a[21] = xr;
        xr = a[29] - xr;    a[29] = xr;

        xr = a[17] - a[1];  a[17] = xr;
        xr = a[ 9] - xr;    a[ 9] = xr;
        xr = a[25] - xr;    a[25] = xr;
        xr = a[ 5] - xr;    a[ 5] = xr;
        xr = a[21] - xr;    a[21] = xr;
        xr = a[13] - xr;    a[13] = xr;
        xr = a[29] - xr;    a[29] = xr;

        xr = a[ 1] - a[0];  a[ 1] = xr;
        xr = a[16] - xr;    a[16] = xr;
        xr = a[17] - xr;    a[17] = xr;
        xr = a[ 8] - xr;    a[ 8] = xr;
        xr = a[ 9] - xr;    a[ 9] = xr;
        xr = a[24] - xr;    a[24] = xr;
        xr = a[25] - xr;    a[25] = xr;
        xr = a[ 4] - xr;    a[ 4] = xr;
        xr = a[ 5] - xr;    a[ 5] = xr;
        xr = a[20] - xr;    a[20] = xr;
        xr = a[21] - xr;    a[21] = xr;
        xr = a[12] - xr;    a[12] = xr;
        xr = a[13] - xr;    a[13] = xr;
        xr = a[28] - xr;    a[28] = xr;
        xr = a[29] - xr;    a[29] = xr;

        xr = a[ 0];  a[ 0] += a[31];  a[31] -= xr;
        xr = a[ 1];  a[ 1] += a[30];  a[30] -= xr;
        xr = a[16];  a[16] += a[15];  a[15] -= xr;
        xr = a[17];  a[17] += a[14];  a[14] -= xr;
        xr = a[ 8];  a[ 8] += a[23];  a[23] -= xr;
        xr = a[ 9];  a[ 9] += a[22];  a[22] -= xr;
        xr = a[24];  a[24] += a[ 7];  a[ 7] -= xr;
        xr = a[25];  a[25] += a[ 6];  a[ 6] -= xr;
        xr = a[ 4];  a[ 4] += a[27];  a[27] -= xr;
        xr = a[ 5];  a[ 5] += a[26];  a[26] -= xr;
        xr = a[20];  a[20] += a[11];  a[11] -= xr;
        xr = a[21];  a[21] += a[10];  a[10] -= xr;
        xr = a[12];  a[12] += a[19];  a[19] -= xr;
        xr = a[13];  a[13] += a[18];  a[18] -= xr;
        xr = a[28];  a[28] += a[ 3];  a[ 3] -= xr;
        xr = a[29];  a[29] += a[ 2];  a[ 2] -= xr;

        /* Compensate for inversion in the analysis filter */
        if (k & 1)
        {
            for (int band = 1; band < 32; band += 2)
                a[band] *= -1;
        }
    }
}

static inline void window_subband2_m(int a0[32])
{
    window_subband2_(a0);
}

static inline void window_subband2_s(int a0[32], int a1[32])
{
    window_subband2_(a0);
    window_subband2_(a1);
}

static inline void mdct_long(int *out, int *in)
{
    int ct, st;
    int tc1, tc2, tc3, tc4, ts5, ts6, ts7, ts8;
    int ts1, ts2, ts3, ts4, tc5, tc6, tc7, tc8;

    /* 1,2, 5,6, 9,10, 13,14, 17 */
    tc1 = in[17] - in[ 9];
    tc3 = in[15] - in[11];
    tc4 = in[14] - in[12];
    ts5 = in[ 0] + in[ 8];
    ts6 = in[ 1] + in[ 7];
    ts7 = in[ 2] + in[ 6];
    ts8 = in[ 3] + in[ 5];

    out[17] = (ts5 + ts7 - ts8) * cx[8] - (ts6 - in[4]) * cx[8];
    st      = (ts5 + ts7 - ts8) * cx[7] + (ts6 - in[4]) * cx[8];
    ct      = (tc1 - tc3 - tc4) * cx[6];
    out[5]  = ct + st;
    out[6]  = ct - st;

    tc2     = (in[16] - in[10]) * cx[6];
    ts6     =  ts6 * cx[7] + in[4] * cx[8];

    ct      =  tc1 * cx[0] + tc2 + tc3 * cx[1] + tc4 * cx[2];
    st      = -ts5 * cx[4] + ts6 - ts7 * cx[5] + ts8 * cx[3];
    out[1]  = ct + st;
    out[2]  = ct - st;

    ct      =  tc1 * cx[1] - tc2 - tc3 * cx[2] + tc4 * cx[0];
    st      = -ts5 * cx[5] + ts6 - ts7 * cx[3] + ts8 * cx[4];
    out[ 9] = ct + st;
    out[10] = ct - st;

    ct      = tc1 * cx[2] - tc2 + tc3 * cx[0] - tc4 * cx[1];
    st      = ts5 * cx[3] - ts6 + ts7 * cx[4] - ts8 * cx[5];
    out[13] = ct + st;
    out[14] = ct - st;

    ts1 = in[ 8] - in[ 0];
    ts3 = in[ 6] - in[ 2];
    ts4 = in[ 5] - in[ 3];
    tc5 = in[17] + in[ 9];
    tc6 = in[16] + in[10];
    tc7 = in[15] + in[11];
    tc8 = in[14] + in[12];

    out[0]  = (tc5 + tc7 + tc8) * cx[8] + (tc6 + in[13]) * cx[8];
    ct      = (tc5 + tc7 + tc8) * cx[7] - (tc6 + in[13]) * cx[8];
    st      = (ts1 - ts3 + ts4) * cx[6];
    out[11] = ct + st;
    out[12] = ct - st;

    ts2     = (in[7] - in[1]) * cx[6];
    tc6     = in[13] * cx[8] - tc6 * cx[7];

    ct      = tc5 * cx[3] - tc6 + tc7 * cx[4] + tc8 * cx[5];
    st      = ts1 * cx[2] + ts2 + ts3 * cx[0] + ts4 * cx[1];
    out[3]  = ct + st;
    out[4]  = ct - st;

    ct      =-tc5 * cx[5] + tc6 - tc7 * cx[3] - tc8 * cx[4];
    st      = ts1 * cx[1] + ts2 - ts3 * cx[2] - ts4 * cx[0];
    out[7]  = ct + st;
    out[8]  = ct - st;

    ct      =-tc5 * cx[4] + tc6 - tc7 * cx[5] - tc8 * cx[3];
    st      = ts1 * cx[0] - ts2 + ts3 * cx[1] - ts4 * cx[2];
    out[15] = ct + st;
    out[16] = ct - st;
}

static int find_bitrate_index(int type, int bitrate, bool stereo)
{
    if (type == 1 && !stereo && bitrate > 160)
        bitrate = 160;

//    return ci->round_value_to_list32(bitrate, &bitr_index[type][1], 14, true) + 1;

	int idx = 0;
	for(int i=0; i < 15; i++)
	{
		if(bitrate == bitr_index[type][i])
		{
			idx = i + 1;
			break;
		}
	}
}

static int find_samplerate_index(long freq, int *mp3_type)
{
    int mpeg = freq >= (32000 + 24000) / 2 ? 1 : 0;
//    int i = ci->round_value_to_list32(freq, sampr_index[mpeg], 3, true);
	int i;
    if(mpeg) //MPEG 1
    {
		switch (freq)
		{
			case 44100:
				i = 0;
				break;
			case 48000:
				i = 1;
				break;
			case 32000:
				i = 2;
				break;
			default:
				break;
		}
    }
    else	//MPEG 2
    {
		switch (freq)
		{
			case 22050:
				i = 0;
				break;
			case 24000:
				i = 1;
				break;
			case 16000:
				i = 2;
				break;
			default:
				break;
		}	
    }
    
    *mp3_type = mpeg;
    return i;
}

static void mp3_encoder_reset(void)
{
    memset(&cfg.cod_info, 0, sizeof (cfg.cod_info));
    memset(mfbuf        , 0, sizeof (mfbuf       ));
    memset(sb_data      , 0, sizeof (sb_data     ));
    memset(mdct_freq    , 0, sizeof (mdct_freq   ));
    memset(mdct_sign    , 0, sizeof (mdct_sign   ));
    memset(enc_data     , 0, sizeof (enc_data    ));
    memset(&coded_data  , 0, sizeof (coded_data  ));
    memset(band_scale_f , 0, sizeof (band_scale_f));
    cfg.slot_lag = 0;
}

static void mp3_encoder_init(unsigned long sample_rate, int num_channels,
                             unsigned long bitrate)
{
    mp3_encoder_reset();

    const bool stereo  = num_channels > 1;
    cfg.channels       = stereo ? 2 : 1;
    cfg.mpg.mode       = stereo ? 0 : 3; /* 0=stereo, 3=mono */
    cfg.mpg.smpl_id    = find_samplerate_index(sample_rate, &cfg.mpg.type);
    cfg.samplerate     = sampr_index[cfg.mpg.type][cfg.mpg.smpl_id];
    cfg.src_samplerate = sample_rate;
    cfg.mpg.bitr_id    = find_bitrate_index(cfg.mpg.type, bitrate, stereo);
    cfg.mpg.bitrate    = bitr_index[cfg.mpg.type][cfg.mpg.bitr_id];//128
    cfg.mpg.num_bands  = num_bands[stereo ? cfg.mpg.type : 2][cfg.mpg.bitr_id];

    if (cfg.mpg.type == 1)
    {
        cfg.granules       = 2;
        cfg.samp_per_frame = 1152;
        cfg.flush_frames   = 2;
    }
    else
    {
        cfg.granules       = 1;
        cfg.samp_per_frame = 576;
        cfg.flush_frames   = 3;
    }

    cfg.delay   = 576-16;
    cfg.padding = 3*576+16;

    cfg.samp_buffer = mfbuf + cfg.channels*512;

    memcpy(scalefac, sf_band[cfg.mpg.smpl_id + 3*cfg.mpg.type], sizeof (scalefac));

    /* Figure average number of 'bytes' per frame */
    cfg.byte_per_frame = calc_frame_size(cfg.mpg.bitr_id, &cfg.frac_per_frame);
    cfg.sideinfo_len   = 32 + (cfg.mpg.type ? (cfg.channels == 1 ? 136 : 256)
                                            : (cfg.channels == 1 ?  72 : 136));

    #if 0
    cfg.req_byte_per_frame = ALIGN_UP(cfg.byte_per_frame + 1,
                                      sizeof (uint32_t));
	#else
	printf("cfg.byte_per_frame = %d\n", cfg.byte_per_frame);

	cfg.req_byte_per_frame = cfg.byte_per_frame + 3;
	#endif
}

static void set_scale_facs(int *mdct_freq)
{
    int avrg_freq_val = 0;

    /* calc average of first 256 frequency values */
    for (int i = 0; i < 256; i++)
        avrg_freq_val += mdct_freq[i];

    avrg_freq_val >>= 8;

    /* if max of current band is smaller than average, increase precision */
    /* last band keeps untouched (not scaled) */
    for (unsigned int k = 0, is = 0; is < scalefac[20]; k++)
    {
        int max_freq_val = 0;
        unsigned ie = scalefac[k];

        for (unsigned int i = is; i < ie; i++)
        {
            if (max_freq_val < mdct_freq[i])
                max_freq_val = mdct_freq[i];
        }

        unsigned int s = 0;
        for (; s < 3; s++)
        {
            if ((max_freq_val << s) > avrg_freq_val)
                break;
        }

        band_scale_f[k] = (unsigned char)s;

        for (unsigned int i = is; s && i < ie; i++)
            mdct_freq[i] <<= s;

        is = ie;
    }
}

static inline void window_subband(int gr, int sbd[2][2][18][32])
{
    int *a0 = sbd[0][1-gr][0];

    if (cfg.channels == 1)
    {
        window_subband1_m(mfbuf + 286 + gr*576, a0);
        window_subband2_m(a0);
    }
    else
    {
        int *a1 = sbd[1][1-gr][0];
        window_subband1_s(mfbuf + 572 + gr*1152, a0, a1);
        window_subband2_s(a0, a1);
    }
}

static inline void get_sb_data(int gr)
{
    window_subband(gr, sb_data);
}

/* Encode one mp3 frame */
static void compress_frame(void)
{
    coded_data.bitpos = cfg.sideinfo_len; /* leave space for mp3 header */
    memset(coded_data.bbuf, 0, sizeof (coded_data.bbuf));

    if ((cfg.slot_lag += cfg.frac_per_frame) >= 64)
    {
        /* Padding for this frame */
        cfg.slot_lag   -= 64;
        cfg.mpg.padding = 1;
    }
    else
    {
        cfg.mpg.padding = 0;
    }

    cfg.mean_bits = (8 * cfg.byte_per_frame + 8 * cfg.mpg.padding
                       - cfg.sideinfo_len) / cfg.granules / cfg.channels
                       - 42; // reserved for scale_facs

    cfg.resv_size = 0;

    int gr_cnt = cfg.granules * cfg.channels;
    for (int gr = 0; gr < cfg.granules; gr++)
    {
        get_sb_data(gr);

        for (int ch = 0; ch < cfg.channels; ch++, gr_cnt--)
        {
            /* Perform imdct of 18 previous + 18 current subband samples */
            /* for integer precision do this loop again (if neccessary)  */
            int shift = 14 - (cfg.cod_info[gr][ch].additStep >> 2);

            for (int ii = 0; ii < 3; ii++)
            {
                int *mdct = mdct_freq;

                cfg.cod_info[gr][ch].additStep = 4 * (14 - shift);

                for (int band = 0; band < cfg.mpg.num_bands; band++, mdct += 18)
                {
                    int *band0 = sb_data[ch][  gr][0] + order[band];
                    int *band1 = sb_data[ch][1-gr][0] + order[band];
                    int work[18];

                    /* 9216=4*32*9*8 */
                    for (int k = -9; k < 0; k++)
                    {
                        int a = shft_n(band1[(k+9)*32], shift);
                        int b = shft_n(band1[(8-k)*32], shift);
                        int c = shft_n(band0[(k+9)*32], shift);
                        int d = shft_n(band0[(8-k)*32], shift);

                        work[k+ 9] = shft16(a * win[k+ 9][0] +
                                            b * win[k+ 9][1] +
                                            c * win[k+ 9][2] +
                                            d * win[k+ 9][3]);

                        work[k+18] = shft16(c * win[k+18][0] +
                                            d * win[k+18][1] +
                                            a * win[k+18][2] +
                                            b * win[k+18][3]);
                    }

                    /* 7200=4*18*100 */
                    mdct_long(mdct, work);

                    /* Perform aliasing reduction butterfly */
                    if (band != 0)
                    {
                        for (int k = 7; k >= 0; k--)
                        {
                            int bu = shft15(mdct[k]) * ca[k] +
                                     shft15(mdct[-1-k]) * cs[k];
                            int bd = shft15(mdct[k]) * cs[k] -
                                     shft15(mdct[-1-k]) * ca[k];
                            mdct[-1-k] = bu;
                            mdct[ k  ] = bd;
                        }
                    }
                }

                uint32_t max = 0;

                for (int k = 0; k < 576; k++)
                {
                    if (mdct_freq[k] < 0)
                    {
                        mdct_sign[k] = 1; /* negative */
                        mdct_freq[k] = shft13(-mdct_freq[k]);
                    }
                    else
                    {
                        mdct_sign[k] = 0; /* positive */
                        mdct_freq[k] = shft13(mdct_freq[k]);
                    }

                    if (max < (uint32_t)mdct_freq[k])
                        max = (uint32_t)mdct_freq[k];
                }

                cfg.cod_info[gr][ch].max_val = max;

                /* calc new shift for higher integer precision */
                int i = 0;
                while (max < (0x7800u >> i))   i++, shift--;
                while ((max >> i) >= 0x10000u) i++, shift++;
                if (i == 0) break;
                if (shift < 0) shift = 0;
            }

            cfg.cod_info[gr][ch].quantStep +=
                                cfg.cod_info[gr][ch].additStep;

            set_scale_facs(mdct_freq);

            /* bit and noise allocation */
            iteration_loop(mdct_freq, &cfg.cod_info[gr][ch], gr_cnt);

            /* write the frame to the bitstream */
            huffman_code_bits(enc_data, mdct_sign, &cfg.cod_info[gr][ch]);

            cfg.cod_info[gr][ch].quantStep -=
                                cfg.cod_info[gr][ch].additStep;

            if (cfg.granules == 1)
            {
                memcpy(sb_data[ch][0], sb_data[ch][1],sizeof (sb_data[ch][0]));
            }
        }
    }

    /* finish this chunk by adding sideinfo header data */
    coded_data.bitpos = 0;
    encode_side_info(cfg.cod_info);

}

/* Process the PCM data in the encoder's input buffer */
static void mp3_enc_encode_frame(void)
{
    compress_frame();

    /* shift out old samples */
    memmove(mfbuf, mfbuf + cfg.channels*cfg.granules*576,cfg.channels*2*512);
}

static inline unsigned short swap16(unsigned short value)
	   /*
		 result[15..8] = value[ 7..0];
		 result[ 7..0] = value[15..8];
	   */
{
   return (value >> 8) | (value << 8);
}

static inline unsigned long swap32(unsigned long value)
	  /*
		result[31..24] = value[ 7.. 0];
		result[23..16] = value[15.. 8];
		result[15.. 8] = value[23..16];
		result[ 7.. 0] = value[31..24];
	  */
{
  unsigned long hi = swap16(value >> 16);
  unsigned long lo = swap16(value & 0xffff);
  return (lo << 16) | hi;
}


/* Get the last encoded frame */
static size_t mp3_enc_get_frame(uint8_t *outbuf)
{
    long size = cfg.byte_per_frame + cfg.mpg.padding;

#define ROCKBOX_LITTLE_ENDIAN

#ifdef ROCKBOX_LITTLE_ENDIAN
    /* convert frame to big endian */
    const uint32_t *src = coded_data.bbuf;
    uint32_t *dst = (uint32_t *)outbuf;

    for (long i = 0; i < size; i += sizeof (uint32_t))
        *dst++ = swap32(*src++);
#else
    memcpy(outbuf, coded_data.bbuf, size);
#endif /* ROCKBOX_LITTLE_ENDIAN */

    return size;
}


/*======== Codec section ========*/

/* CRC code lovingly ripped from:
 * github.com/CFR-maniac/lame/blob/master/libmp3lame/VbrTag.c */

/* Lookup table for fast CRC computation
 * See 'crc_update_lookup'
 * Uses the polynomial x^16+x^15+x^2+1 */
static const uint16_t crc16_lookup[256] =
{
    0x0000, 0xc0c1, 0xc181, 0x0140, 0xc301, 0x03c0, 0x0280, 0xc241,
    0xc601, 0x06c0, 0x0780, 0xc741, 0x0500, 0xc5c1, 0xc481, 0x0440,
    0xcc01, 0x0cc0, 0x0d80, 0xcd41, 0x0f00, 0xcfc1, 0xce81, 0x0e40,
    0x0a00, 0xcac1, 0xcb81, 0x0b40, 0xc901, 0x09c0, 0x0880, 0xc841,
    0xd801, 0x18c0, 0x1980, 0xd941, 0x1b00, 0xdbc1, 0xda81, 0x1a40,
    0x1e00, 0xdec1, 0xdf81, 0x1f40, 0xdd01, 0x1dc0, 0x1c80, 0xdc41,
    0x1400, 0xd4c1, 0xd581, 0x1540, 0xd701, 0x17c0, 0x1680, 0xd641,
    0xd201, 0x12c0, 0x1380, 0xd341, 0x1100, 0xd1c1, 0xd081, 0x1040,
    0xf001, 0x30c0, 0x3180, 0xf141, 0x3300, 0xf3c1, 0xf281, 0x3240,
    0x3600, 0xf6c1, 0xf781, 0x3740, 0xf501, 0x35c0, 0x3480, 0xf441,
    0x3c00, 0xfcc1, 0xfd81, 0x3d40, 0xff01, 0x3fc0, 0x3e80, 0xfe41,
    0xfa01, 0x3ac0, 0x3b80, 0xfb41, 0x3900, 0xf9c1, 0xf881, 0x3840,
    0x2800, 0xe8c1, 0xe981, 0x2940, 0xeb01, 0x2bc0, 0x2a80, 0xea41,
    0xee01, 0x2ec0, 0x2f80, 0xef41, 0x2d00, 0xedc1, 0xec81, 0x2c40,
    0xe401, 0x24c0, 0x2580, 0xe541, 0x2700, 0xe7c1, 0xe681, 0x2640,
    0x2200, 0xe2c1, 0xe381, 0x2340, 0xe101, 0x21c0, 0x2080, 0xe041,
    0xa001, 0x60c0, 0x6180, 0xa141, 0x6300, 0xa3c1, 0xa281, 0x6240,
    0x6600, 0xa6c1, 0xa781, 0x6740, 0xa501, 0x65c0, 0x6480, 0xa441,
    0x6c00, 0xacc1, 0xad81, 0x6d40, 0xaf01, 0x6fc0, 0x6e80, 0xae41,
    0xaa01, 0x6ac0, 0x6b80, 0xab41, 0x6900, 0xa9c1, 0xa881, 0x6840,
    0x7800, 0xb8c1, 0xb981, 0x7940, 0xbb01, 0x7bc0, 0x7a80, 0xba41,
    0xbe01, 0x7ec0, 0x7f80, 0xbf41, 0x7d00, 0xbdc1, 0xbc81, 0x7c40,
    0xb401, 0x74c0, 0x7580, 0xb541, 0x7700, 0xb7c1, 0xb681, 0x7640,
    0x7200, 0xb2c1, 0xb381, 0x7340, 0xb101, 0x71c0, 0x7080, 0xb041,
    0x5000, 0x90c1, 0x9181, 0x5140, 0x9301, 0x53c0, 0x5280, 0x9241,
    0x9601, 0x56c0, 0x5780, 0x9741, 0x5500, 0x95c1, 0x9481, 0x5440,
    0x9c01, 0x5cc0, 0x5d80, 0x9d41, 0x5f00, 0x9fc1, 0x9e81, 0x5e40,
    0x5a00, 0x9ac1, 0x9b81, 0x5b40, 0x9901, 0x59c0, 0x5880, 0x9841,
    0x8801, 0x48c0, 0x4980, 0x8941, 0x4b00, 0x8bc1, 0x8a81, 0x4a40,
    0x4e00, 0x8ec1, 0x8f81, 0x4f40, 0x8d01, 0x4dc0, 0x4c80, 0x8c41,
    0x4400, 0x84c1, 0x8581, 0x4540, 0x8701, 0x47c0, 0x4680, 0x8641,
    0x8201, 0x42c0, 0x4380, 0x8341, 0x4100, 0x81c1, 0x8081, 0x4040
};

static uint32_t header_size ;
static unsigned int mp3_crc16 ;

/* fast CRC-16 computation - uses table crc16_lookup 8*/
static inline unsigned int crc_update_lookup(unsigned int value,
                                             unsigned int crc)
{
    unsigned int tmp = crc ^ value;
    crc = (crc >> 8) ^ crc16_lookup[tmp & 0xff];
    return crc & 0xffff;
}

/* Calculate position of 'Info' header */
static int get_info_offset(uint32_t header)
{
    uint32_t type = (header & (0x3 << 19)) >> 19;
    uint32_t mode = (header & (0x3 << 6)) >> 6;

    return type == 3 ? (mode == 3 ? 21 : 36) : (mode == 3 ? 13 : 21);
}

/* Write very basic 'Info' header with delay, padding and a bit of
 * miscellaneous info. */
static bool write_info_header(bool first_encode)
{
    uint32_t size = cfg.byte_per_frame;

    /* By default the MP3 frame header for the info frame is the same as
       unpadded audio frames */
    uint32_t header = encode_header(0, cfg.mpg.bitr_id);

    int i = get_info_offset(header);

    if (i + 8 + 36 > size)
    {
        /* The default frame size too small so find the smallest one that
           may accomodate it by increasing the bit rate for this empty
           MP3 frame */
        int j;
        for (j = cfg.mpg.bitr_id + 1; j < 15; j++)
        {
            size = calc_frame_size(j, NULL);

            if (size >= i + 8 + 36)
                break;
        }

        if (j >= 15)
        {
            /* Shouldn't really happen but... */
            header_size = -1;
            return true;
        }

        header = encode_header(0, j);
        /* Info offset won't change */
    }

    uint8_t frame[size];
    memset(frame, 0, size);

    frame[0] = header >> 24;
    frame[1] = header >> 16;
    frame[2] = header >>  8;
    frame[3] = header >>  0;

    /* 'Info' header (CBR 'Xing') */
    memcpy(&frame[i], "Info", 4);

    /* flags = 0; Info contains no other sections and is 8 bytes */

    /* Just mark the LAMEness to indicate header presence; we're not
       actually _the_ LAME so 'rbshn' is the version we give */
    memcpy(&frame[i + 8], "LAMErbshn", 9);

    /* Fill-in some info about us
     * reference: http://gabriel.mp3-tech.org/mp3infotag.html
     */

    /* Revision + VBR method:
     * [7:4] = Revision (0 ??)
     * [3:0] = VBR method (CBR)
     */
    frame[i + 17] = (0 << 4) | (1 << 0);

    /* If first frame since encoder reset is long gone (not unlikely in
       prerecording), then the delay is long passed and no trimming done
       at the start */
    unsigned int delay = first_encode ? cfg.delay : 0;
    unsigned int padding = cfg.padding;

    /* Delay and padding:
     * [23:12] = delay
     * [11: 0] = padding
     */
    frame[i + 29] = delay >> 4;
    frame[i + 30] = (delay << 4) | (padding >> 8);
    frame[i + 31] = padding;

    /* Misc:
     * [7:6] = source frequency
     * [  5] = unwise settings (of course not :)
     * [4:2] = stereo mode (mono or stereo)
     * [1:0] = noise shaping (who knows, 0)
     */
    uint8_t misc;

    if (cfg.src_samplerate <= 32000)
        misc = (0 << 6);
    else if (cfg.src_samplerate <= 44100)
        misc = (1 << 6);
    else if (cfg.src_samplerate <= 48000)
        misc = (2 << 6);
    else /* > 48000 */
        misc = (3 << 6);

    if (cfg.channels > 1)
        misc |= (1 << 2); /* Stereo */

    frame[i + 32] = misc;

    if (enc_stream_write(frame, size) != size)
    {
        enc_stream_lseek(0, SEEK_SET);
        header_size = -1;
        return false;
    }

    header_size = size;
    return true;
}

static inline void encode_frame(void)
{
    mp3_enc_encode_frame();
}

static inline size_t get_frame(uint8_t *buffer)
{
    return mp3_enc_get_frame(buffer);
}

#if 1  //for test
/* this is called for each file to process */
FILE *in_file, *out_file;
uint8_t out_buf[1024];  //usually is 418 byte
enum codec_status codec_run(void)
{
    mp3_encoder_reset();

	int loop_cnt = 0;
	while(enc_pcmbuf_read(cfg.samp_buffer, cfg.samp_per_frame * 2))
    {
    	loop_cnt++;
    	printf("loop%d, ", loop_cnt);
    	encode_frame();

    	int size = get_frame(out_buf);
    	printf("size = %d\n", size);
    	enc_stream_write(out_buf, size);
    } /* while */

    return CODEC_OK;
}

#define IN_FILE   "F:\\c\\mp3_enc\\test_file\\in.wav"
#define OUT_FILE  "F:\\c\\mp3_enc\\test_file\\out.mp3"

int main()
{
	printf("mp3 enc!\n");

    in_file = fopen(IN_FILE, "rb");
    out_file = fopen(OUT_FILE, "wb");

    if(!in_file || !out_file)
    {
        printf("open file err!\n");
    }

	mp3_encoder_init(44100, 2, 128);  //default config, do not change
	codec_run();

	fclose(in_file);
	fclose(out_file);

	return 0;
}
#endif

/*============================================================================
* mp3_enc port
*===========================================================================*/
uint8_t enc_out_buf[1024];  //usually is 418 byte
int enc_stream_write(uint8_t *data, uint32_t size)
{
	return fwrite(data, sizeof(char), size, out_file);
}

int enc_stream_lseek(uint32_t x, uint32_t y)
{

}

uint32_t enc_pcmbuf_read(uint8_t *buf, uint32_t size)
{
	uint32_t read_cnt = 0;
	read_cnt = fread(buf, sizeof(short), size, in_file);

	printf("read_cnt = %d ", read_cnt);
	return read_cnt;
}

int enc_pcmbuf_advance()
{

}

/*============================================================================
* mp3_enc api
*===========================================================================*/
void mp3_enc_init()
{
	mp3_encoder_init(48000, 2, 128);  

	printf("mp3 enc init succ!\n");
}

int mp3_enc_run()
{
	if(enc_pcmbuf_read(cfg.samp_buffer, cfg.samp_per_frame * 2))
	{
		encode_frame();

		int size = get_frame(enc_out_buf);
		enc_stream_write(enc_out_buf, size);
		printf("size = %d\n", size);
		return 1;
	}
	else
		return 0;
}
