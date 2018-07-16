#ifndef MP3_ENC_H
#define MP3_ENC_H

#include <stdint.h>

/* codec return codes */
enum codec_status {
    CODEC_OK = 0,
    CODEC_ERROR = -1,
};

typedef struct
{
    int   type; /* 0=(MPEG2 - 22.05,24,16kHz) 1=(MPEG1 - 44.1,48,32kHz) */
    int   mode; /* 0=stereo, 1=jstereo, 2=dual, 3=mono  */
    int   bitrate;
    int   padding;
    int   num_bands;
    long  bitr_id;
    int   smpl_id;
} mpeg_t;

/* Side information */
typedef struct
{
    uint32_t part2_3_length;
    int      count1;          /* number of 0-1-quadruples  */
    uint32_t global_gain;
    uint32_t table_select[4];
    uint32_t region_0_1;
    uint32_t address1;
    uint32_t address2;
    uint32_t address3;
    long     quantStep;
    long     additStep;
    uint32_t max_val;
} side_info_t;

typedef struct
{
    side_info_t cod_info[2][2];
    mpeg_t   mpg;
    long     frac_per_frame;
    long     byte_per_frame;
    long     req_byte_per_frame;
    long     slot_lag;
    int      sideinfo_len;
    int      mean_bits;
    int      resv_size;
    int      channels;
    int      granules;
    long     src_samplerate;
    long     samplerate;
    short   *samp_buffer;
    unsigned samp_per_frame;
    int      flush_frames;
    int      delay;
    int      padding;
} config_t;

typedef struct
{
    int      bitpos;   /* current bitpos for writing */
    uint32_t bbuf[362];
} bf_data;

struct huffcodetab
{
    int            len;   /* max. index                  */
    const uint8_t *table; /* pointer to array[len][len]  */
    const uint8_t *hlen;  /* pointer to array[len][len]  */
};

struct huffcodebig
{
    int len;     /* max. index                   */
    int linbits; /* number of linbits            */
    int linmax;  /* max number stored in linbits */
};

#define shft4(x)     (((x) +     8) >>  4)
#define shft9(x)     (((x) +   256) >>  9)
#define shft13(x)    (((x) +  4096) >> 13)
#define shft15(x)    (((x) + 16384) >> 15)
#define shft16(x)    (((x) + 32768) >> 16)
#define shft_n(x, n) ( (x) >> (n))
#define SQRT         724 /* sqrt(2) * 512 */


#endif
