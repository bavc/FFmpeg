/*
 ** <h3>values graph</h3>
 * copyright (c) 2010 Mark Heath mjpeg0 @ silicontrip dot org
 * http://silicontrip.net/~mark/lavtools/
 *
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include <stdio.h>

#include "libavutil/opt.h"
#include "libavutil/pixdesc.h"
#include "internal.h"

// dB = 10 * log10 (r1/r2)

enum FilterMode {
    FILTER_NONE = -1,
    FILTER_TOUT,
    FILTER_VREP,
    FILTER_RANGE,
    FILTER_HEADSWITCHING,
    FILT_NUMB
};


typedef struct
{
    const AVClass *class;

    FILE  *fh;
    char *filename;

    int chromah;
    int chromaw;
     int hsub;
     int vsub;

    int fc;

    int fs;
    int cfs;

    enum FilterMode outfilter;
    int filters;

    AVFrame *frame_prev;

    int filter_vrep_prev;

    char* vrep_line;
    unsigned int * filter_head_border;

} valuesContext;


/* Prototypes for filter functions */

static int filter_tout(valuesContext *values, const AVFrame *p, int y, int w, int h);
static int filter_vrep(valuesContext *values, const AVFrame *p, int y, int w, int h);
static int filter_range(valuesContext *values, const AVFrame *p,int y, int w, int h);
static int filter_head(valuesContext *values, const AVFrame *p, int y, int w, int h);


static void filter_init_tout(valuesContext *values, const AVFrame *p, int w, int h);
static void filter_init_vrep(valuesContext *values, const AVFrame *p, int w, int h);
static void filter_init_range(valuesContext *values, const AVFrame *p, int w, int h);
static void filter_init_head(valuesContext *values, const AVFrame *p, int w, int h);

static void filter_uninit_tout(valuesContext *values);
static void filter_uninit_vrep(valuesContext *values);
static void filter_uninit_range(valuesContext *values);
static void filter_uninit_head(valuesContext *values);



static int (*filter_call[FILT_NUMB])(valuesContext *values, const AVFrame *p, int y, int w, int h) = {
    filter_tout,
    filter_vrep,
    filter_range,
    filter_head
};

static void (*filter_init[FILT_NUMB])(valuesContext *values, const AVFrame *p,  int w, int h) = {
    filter_init_tout,
    filter_init_vrep,
    filter_init_range,
    filter_init_head
};

static void (*filter_uninit[FILT_NUMB])(valuesContext *values) = {
    filter_uninit_tout,
    filter_uninit_vrep,
    filter_uninit_range,
    filter_uninit_head
};


static const char *const filter_metanames[] = { "TOUT", "VREP", "RANG", "HEAD", NULL };
static const char *const filter_names[]     = { "tout", "vrep", "rang", "head", NULL };

/* end of filter definitions */


#define OFFSET(x) offsetof(valuesContext, x)
#define FLAGS AV_OPT_FLAG_FILTERING_PARAM|AV_OPT_FLAG_VIDEO_PARAM

static const AVOption values_options[] = {
    {"filename", "set output file", OFFSET(filename), AV_OPT_TYPE_STRING, {.str=NULL}, .flags=FLAGS},
    {"f",        "set output file", OFFSET(filename), AV_OPT_TYPE_STRING, {.str=NULL}, .flags=FLAGS},
    {"out", "set video filter", OFFSET(outfilter), AV_OPT_TYPE_INT, {.i64=FILTER_NONE}, -1, FILT_NUMB-1, FLAGS, "out"},
        {"tout", "", 0, AV_OPT_TYPE_CONST, {.i64=FILTER_TOUT},  0, 0, FLAGS, "out"},
        {"vrep", "", 0, AV_OPT_TYPE_CONST, {.i64=FILTER_VREP},  0, 0, FLAGS, "out"},
        {"rang", "", 0, AV_OPT_TYPE_CONST, {.i64=FILTER_RANGE}, 0, 0, FLAGS, "out"},
        {"head", "", 0, AV_OPT_TYPE_CONST, {.i64=FILTER_HEADSWITCHING}, 0, 0, FLAGS, "out"},
    {"stat", "statistics filters", OFFSET(filters), AV_OPT_TYPE_FLAGS, {.i64=0}, 0, INT_MAX, FLAGS, "filters"},
        {"tout", "", 0, AV_OPT_TYPE_CONST, {.i64=1<<FILTER_TOUT},          0, 0, FLAGS, "filters"},
        {"vrep", "", 0, AV_OPT_TYPE_CONST, {.i64=1<<FILTER_VREP},          0, 0, FLAGS, "filters"},
        {"rang", "", 0, AV_OPT_TYPE_CONST, {.i64=1<<FILTER_RANGE},         0, 0, FLAGS, "filters"},
        {"head", "", 0, AV_OPT_TYPE_CONST, {.i64=1<<FILTER_HEADSWITCHING}, 0, 0, FLAGS, "filters"},
    {NULL}
};

AVFILTER_DEFINE_CLASS(values);

static av_cold int init(AVFilterContext *ctx)
{
    valuesContext *values = ctx->priv;

    if (values->filename)
        values->fh = fopen(values->filename, "w");

    av_log(ctx, AV_LOG_DEBUG, "<<< init().\n");

    return 0;
}

static av_cold void uninit(AVFilterContext *ctx)
{
    valuesContext *values = ctx->priv;
    if (values->fh)
        fclose(values->fh);
    av_frame_free(&values->frame_prev);
}

static int query_formats(AVFilterContext *ctx)
{
    // will want more
    enum PixelFormat pix_fmts[] = {
        PIX_FMT_YUV444P,  PIX_FMT_YUV422P,  PIX_FMT_YUV420P, PIX_FMT_YUV411P,
        PIX_FMT_NONE
    };

    ff_set_common_formats(ctx, ff_make_format_list(pix_fmts));
    return 0;
}

static int config_props(AVFilterLink *outlink)
{
    AVFilterContext *ctx = outlink->src;
    valuesContext *values = ctx->priv;
    AVFilterLink *inlink = outlink->src->inputs[0];
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(outlink->format);
    values->hsub = desc->log2_chroma_w;
    values->vsub = desc->log2_chroma_h;

    av_log(ctx, AV_LOG_DEBUG, ">>> config_props().\n");

    outlink->w = inlink->w;
    outlink->h = inlink->h;

    values->chromaw = FF_CEIL_RSHIFT(inlink->w, values->hsub);
    values->chromah = FF_CEIL_RSHIFT(inlink->h, values->vsub);

    values->fs = inlink->w * inlink->h;
    values->cfs = values->chromaw * values->chromah;

    av_log(ctx, AV_LOG_DEBUG, "<<< config_props().\n");

    return 0;
}

static void filter_init_head(valuesContext *values, const AVFrame *p, int w, int h)
{

    int y;
    int tol = 16; // this needs to be configurable.
    int lw = p->linesize[0];
    unsigned int *order;
    int median;
    //int switching =1;

    // should check if this fails
    values->filter_head_border = (unsigned int*)malloc (h * sizeof(unsigned int));

    order = (unsigned int*)malloc (h * sizeof(unsigned int));


    for (y=0; y< h; y++)
    {

        int yw = y*lw;
        int x;
        // try not to get fooled by non matted video.
        int col = 16;

        values->filter_head_border[y] = 0;


        if (p->data[0][yw] <= 16)
            col = p->data[0][yw];

        for (x=0; x< w; x++)
        {
            order[y] = x;

            if (abs(col - p->data[0][yw + x]) > tol)
            {
                values->filter_head_border[y] = x;
                break;
            }

        }
    }

    // there is probably a better and easier (and faster) sorting method
    for (y=0; y< h; y++ )
    {
        int min = order[y];
        int select = y;
        int x;

        for (x=y; x<h; x++)
        {
            if (order[x] < min)
            {
                min = order[x];
                select = x;
            }

            if (select != y)
            {
                min = order[y];
                order[y] = order[select];
                order[select] = min;
            }
        }

    }

    // especially if all I want now is the 50th percentile.
    // or 98th depending on what mood the head switching pattern is in.

    median = order[98 * h / 100] ;

    // remove possible matting
    for (y=0;y<h;y++)
    {
        if (values->filter_head_border[y] < median)
            values->filter_head_border[y] = 0;
    }

    /*
    for (y=h-1; y>0; y--)
    {

       if (values->filter_head_border[y] <= median )
           switching = 0;

        if (!switching)
            values->filter_head_border[y] = 0;

    }
    */
    free (order);

}


static int filter_head(valuesContext *values, const AVFrame *p, int y, int w, int h)
{
    return values->filter_head_border[y];
}

static void filter_uninit_head(valuesContext *values)
{
    free(values->filter_head_border);

}

static void filter_init_range(valuesContext *values, const AVFrame *p, int w, int h)
{
    ;
}

static void filter_uninit_range(valuesContext *values)
{
    ;
}


static int filter_range(valuesContext *values, const AVFrame *p, int y, int w, int h)
{
    int x, score = 0;
    const int yc = FF_CEIL_RSHIFT(y, values->vsub);
    const uint8_t *pluma    = &p->data[0][y  * p->linesize[0]];
    const uint8_t *pchromau = &p->data[1][yc * p->linesize[1]];
    const uint8_t *pchromav = &p->data[2][yc * p->linesize[2]];

    for (x = 0; x < w; x++) {
        const int xc = FF_CEIL_RSHIFT(x, values->hsub);
        const int luma    = pluma[x];
        const int chromau = pchromau[xc];
        const int chromav = pchromav[xc];
        score += luma    < 16 || luma    > 235 ||
                 chromau < 16 || chromau > 240 ||
                 chromav < 16 || chromav > 240;
    }
    return score;
}

static void filter_init_tout(valuesContext *values, const AVFrame *p, int w, int h)
{
    ;
}

static void filter_uninit_tout(valuesContext *values)
{
    ;
}

static int filter_tout_outlier(uint8_t x, uint8_t y, uint8_t z)
{
    return ((abs(x - y) + abs (z - y)) / 2) - abs(z - x) > 4; // make 4 configurable?
}

static int filter_tout(valuesContext *values, const AVFrame *f, int y, int w, int h)
{
    const uint8_t *p = f->data[0];
    int lw = f->linesize[0];
    int x, score = 0;

    if (y - 1 < 0 || y + 1 >= h)
        return 0;

    for (x = 1; x < w - 1; x++) {

            // detect two pixels above and below (to eliminate interlace artefacts)
            // should check that video format is infact interlace.

#define FILTER(i, j) \
filter_tout_outlier(p[(y-j) * lw + x + i], \
                    p[    y * lw + x + i], \
                    p[(y+j) * lw + x + i])

#define FILTER3(j) (FILTER(-1, j) && FILTER(0, j) && FILTER(1, j))

        if (y-2 >= 0 && y+2 < h)
            score += FILTER3(2) && FILTER3(1);
        else
            score += FILTER3(1);
    }
    return score;
}

static void filter_init_vrep(valuesContext *values, const AVFrame *p, int w, int h)
{

    int i,y;
    int lw = p->linesize[0];

    values->vrep_line = (char *) malloc (h);

    for (y=4;y<h;y++)
    {
        int totdiff = 0;

        int y2lw = (y-4) * lw;
        int ylw = y * lw;

        for (i = 0; i < w; i++)
            totdiff += abs(p->data[0][y2lw + i] - p->data[0][ylw + i]);

        /* this value should be definable */
        if (totdiff < w ) {
            values->vrep_line[y] = 1;
        }

    }

}

static int filter_vrep(valuesContext *values, const AVFrame *p, int y, int w, int h)
{
    int x, score = 0;

    for (x = 0; x < w; x++)
        score += values->vrep_line[y];
    return score;
}

static void filter_uninit_vrep(valuesContext *values)
{
    free (values->vrep_line);
}




#define DEPTH 256

static int filter_frame(AVFilterLink *link, AVFrame *in)
{

    AVFilterContext *ctx = link->src;
    valuesContext *values = link->dst->priv;
    AVFilterLink *outlink = link->dst->outputs[0];
    AVFrame *out;
    int i, j, direct = 0;
    int  w = 0,  cw = 0, // in
        ow = 0, cow = 0, // out
        pw = 0, cpw = 0; // prev
    int yuv;
    int fil;
    char metabuf[128];
    unsigned int histy[DEPTH] = {0},
                 histu[DEPTH] = {0},
                 histv[DEPTH] = {0}; // limited to 8 bit data.
    int miny, minu, minv;
    int maxy, maxu, maxv;
    int lowy  = -1, lowu  = -1, lowv  = -1;
    int highy = -1, highu = -1, highv = -1;
    int lowp, highp, clowp, chighp;
    int accy, accu, accv;
    int toty = 0, totu = 0, totv = 0;
    int dify = 0, difu = 0, difv = 0;
    int dify1 = 0, dify2 = 0;

    int filtot[FILT_NUMB] = {0};
    AVFrame *prev;

    av_log(ctx, AV_LOG_DEBUG, ">>> filter_frame().\n");

    if (!values->frame_prev)
        values->frame_prev = av_frame_clone(in);
    prev = values->frame_prev;

    if (av_frame_is_writable(in) || values->outfilter == FILTER_NONE) {
        out = in;
        direct = 1;
    } else {
        out = ff_get_video_buffer(outlink, outlink->w, outlink->h);
        if (!out) {
            av_frame_free(&in);
            return AVERROR(ENOMEM);
        }
        av_frame_copy_props(out, in);
    }

    miny = in->data[0][0];
    maxy = in->data[0][0];

    minu = in->data[1][0];
    maxu = in->data[1][0];

    minv = in->data[2][0];
    maxv = in->data[2][0];

    av_log(ctx, AV_LOG_DEBUG, "    filter_frame() for(j=0...)\n");

    for (fil = 0; fil < FILT_NUMB; fil ++) {
        if ((values->filters & 1<<fil) || values->outfilter == fil) {
            filter_init[fil](values, in,link->w,link->h);
        }
    }


    for (j = 0; j < link->h; j++) {
        for (i = 0; i < link->w; i++) {

            yuv = in->data[0][w+i];

            if (yuv > maxy) maxy=yuv;
            if (yuv < miny) miny=yuv;

            toty += yuv;

            histy[yuv]++;

            if (!direct)
                out->data[0][ow+i] = in->data[0][w+i]; // or 16;

            dify += abs(in->data[0][w+i] - prev->data[0][pw+i]);


            //if (in->interlaced_frame && (j % 2 == 0)) // every second line
            if (j % 2 == 0) // every second line
            {

               // get bottom field
               // should check that we are not currently at the bottom line.
                // but who has heard of an interlaced file with odd vertical dimentions?
                int yuvi = in->data[0][w+in->linesize[0]+i];

                // dif2 = diff bottom field with top field

                dify2 += abs(yuv - yuvi);

                if (in->top_field_first)
                {
                    // dif1 = diff top field with prev bottom field
                    dify1 += abs(yuv - prev->data[0][pw+prev->linesize[0]+i]);
                } else {
                    // dif1 = diff bottom field with prev top field
                    dify1 += abs(yuvi - prev->data[0][pw+i]);
                }

            }

            if (i < values->chromaw && j < values->chromah) {

                yuv = in->data[1][cw+i];
                if (yuv > maxu) maxu=yuv;
                if (yuv < minu) minu=yuv;
                totu += yuv;
                histu[yuv]++;

                yuv = in->data[2][cw+i];
                if (yuv > maxv) maxv=yuv;
                if (yuv < minv) minv=yuv;
                totv += yuv;
                histv[yuv]++;

                difu += abs(in->data[1][cw+i] - prev->data[1][cpw+i]);
                difv += abs(in->data[2][cw+i] - prev->data[2][cpw+i]);

                if (!direct) {
                    out->data[1][cow+i] = in->data[1][cw+i]; // or 128
                    out->data[2][cow+i] = in->data[2][cw+i]; // or 128
                }
            }

            // magic filter array
        }

            for (fil = 0; fil < FILT_NUMB; fil ++) {
                if (values->filters & 1<<fil) {
                    int ret = filter_call[fil](values, in, j, link->w, link->h);
                    if (ret) {
                        filtot[fil] += ret;
                        if (values->outfilter == fil)
                            out->data[0][ow+i] = 235;
                    }
                }
            }

        ow  += out->linesize[0];
        cow += out->linesize[1];

        w  += in->linesize[0];
        cw += in->linesize[1];

        pw  += prev->linesize[0];
        cpw += prev->linesize[1];
    }


    // find low / high based on histogram percentile
    // these only need to be calculated once.

    lowp   = values->fs  * 10 / 100;
    highp  = values->fs  * 95 / 100;
    clowp  = values->cfs * 10 / 100;
    chighp = values->cfs * 95 / 100;

    accy = 0; accu=0; accv=0;
    for (fil = 0; fil < DEPTH; fil++) {
        accy += histy[fil];
        accu += histu[fil];
        accv += histv[fil];

        if (lowy == -1 && accy >=  lowp) lowy = fil;
        if (lowu == -1 && accu >= clowp) lowu = fil;
        if (lowv == -1 && accv >= clowp) lowv = fil;

        if (highy == -1 && accy >=  highp) highy = fil;
        if (highu == -1 && accu >= chighp) highu = fil;
        if (highv == -1 && accv >= chighp) highv = fil;
    }

    av_log(ctx, AV_LOG_DEBUG, "    filter_frame() av_frame_free()\n");

    av_frame_free(&values->frame_prev);

    av_log(ctx, AV_LOG_DEBUG, "    filter_frame() prev->data = in->data\n");

    values->frame_prev = av_frame_clone(in);

#define SET_META(key, fmt, val) do {                                \
    snprintf(metabuf, sizeof(metabuf), fmt, val);                   \
    av_dict_set(&out->metadata, "lavfi.values." key, metabuf, 0);   \
} while (0)

    SET_META("YMIN",  "%d", miny);
    SET_META("YLOW",  "%d", lowy);
    SET_META("YAVG",  "%g", 1.0 * toty / values->fs);
    SET_META("YHIGH", "%d", highy);
    SET_META("YMAX",  "%d", maxy);
    SET_META("UMIN",  "%d", minu);
    SET_META("ULOW",  "%d", lowu);
    SET_META("UAVG",  "%g", 1.0 * totu / values->cfs);
    SET_META("UHIGH", "%d", highu);
    SET_META("UMAX",  "%d", maxu);
    SET_META("VMIN",  "%d", minv);
    SET_META("VLOW",  "%d", lowv);
    SET_META("VAVG",  "%g", 1.0 * totv / values->cfs);
    SET_META("VHIGH", "%d", highv);
    SET_META("VMAX",  "%d", maxv);
    SET_META("YDIF",  "%g", 1.0 * dify / values->fs);
    SET_META("UDIF",  "%g", 1.0 * difu / values->cfs);
    SET_META("VDIF",  "%g", 1.0 * difv / values->cfs);
    SET_META("YDIF1",  "%g", 1.0 * dify1 / (values->fs/2));
    SET_META("YDIF2",  "%g", 1.0 * dify2 / (values->fs/2));

    av_log(ctx, AV_LOG_DEBUG, "    filter_frame() for (fil = 0; fil < FILT_NUMB; fil ++).\n");

    for (fil = 0; fil < FILT_NUMB; fil ++) {
        if (values->filters & 1<<fil) {
       //     SET_META(filter_metanames[fil],"%g",1.0 * filtot[fil]/values->fs);

            char metaname[128];
            snprintf(metabuf,sizeof(metabuf),"%g",1.0 * filtot[fil]/values->fs);
            snprintf(metaname,sizeof(metaname),"lavfi.values.%s",filter_metanames[fil]);
            av_dict_set(&out->metadata,metaname,metabuf,0);

        }
    }

    for (fil = 0; fil < FILT_NUMB; fil ++) {
        if ((values->filters & 1<<fil) || values->outfilter == fil) {
            filter_uninit[fil](values);
        }
    }

    values->fc++;
    if (!direct)
        av_frame_free(&in);
    av_log(ctx, AV_LOG_DEBUG, "<<< filter_frame().\n");
    return ff_filter_frame(outlink, out);
}

static const AVFilterPad values_inputs[] = {
    {
        .name            = "default",
        .type            = AVMEDIA_TYPE_VIDEO,
        .min_perms       = AV_PERM_READ,
        .filter_frame = filter_frame,
    },
    { NULL }
};

static const AVFilterPad values_outputs[] = {
    {
        .name            = "default",
        .config_props     = config_props,
        .type            = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};

AVFilter ff_vf_values = {
    .name          = "values",
    .description   = ".",
    .init          = init,
    .uninit        = uninit,
    .query_formats = query_formats,
    .priv_size     = sizeof(valuesContext),
    .inputs        = values_inputs,
    .outputs       = values_outputs,
    .priv_class    = &values_class,
};
