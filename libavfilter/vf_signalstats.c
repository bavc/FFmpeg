/*
 * Copyright (c) 2010 Mark Heath mjpeg0 @ silicontrip dot org
 * Copyright (c) 2014 Clément Bœsch
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

#include "libavutil/opt.h"
#include "libavutil/pixdesc.h"
#include "libavcodec/mathops.h"
#include "libavutil/libm.h"
#include "libavutil/mathematics.h"
#include "internal.h"

enum FilterMode {
    FILTER_NONE = -1,
    FILTER_TOUT,
    FILTER_VREP,
    FILTER_RANGE,
    FILTER_HEADSWITCHING,
    FILT_NUMB
};

typedef struct {
    const AVClass *class;
    FILE *fh;
    char *filename;
    int chromah;
    int chromaw;
    int hsub;
    int vsub;
    int fs;
    int cfs;
    enum FilterMode outfilter;
    int filters;
    AVFrame *frame_prev;
    char *vrep_line;
    unsigned int *filter_head_border;
    unsigned int *filter_head_order;
} signalstatsContext;

#define OFFSET(x) offsetof(signalstatsContext, x)
#define FLAGS AV_OPT_FLAG_FILTERING_PARAM|AV_OPT_FLAG_VIDEO_PARAM

static const AVOption signalstats_options[] = {
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

AVFILTER_DEFINE_CLASS(signalstats);

static av_cold int init(AVFilterContext *ctx)
{
    signalstatsContext *signalstats = ctx->priv;

    if (signalstats->filename)
        signalstats->fh = fopen(signalstats->filename, "w");

    if (signalstats->outfilter != FILTER_NONE)
        signalstats->filters |= 1 << signalstats->outfilter;

    return 0;
}

static av_cold void uninit(AVFilterContext *ctx)
{
    signalstatsContext *signalstats = ctx->priv;
    if (signalstats->fh)
        fclose(signalstats->fh);
    av_frame_free(&signalstats->frame_prev);

    av_freep(&signalstats->filter_head_border);
    av_freep(&signalstats->filter_head_order);

    av_freep(&signalstats->vrep_line);
}

static int query_formats(AVFilterContext *ctx)
{
    // TODO: add more
    enum AVPixelFormat pix_fmts[] = {
        AV_PIX_FMT_YUV444P, AV_PIX_FMT_YUV422P, AV_PIX_FMT_YUV420P, AV_PIX_FMT_YUV411P,
        AV_PIX_FMT_NONE
    };

    ff_set_common_formats(ctx, ff_make_format_list(pix_fmts));
    return 0;
}

static int config_props(AVFilterLink *outlink)
{
    AVFilterContext *ctx = outlink->src;
    signalstatsContext *signalstats = ctx->priv;
    AVFilterLink *inlink = outlink->src->inputs[0];
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(outlink->format);
    signalstats->hsub = desc->log2_chroma_w;
    signalstats->vsub = desc->log2_chroma_h;

    outlink->w = inlink->w;
    outlink->h = inlink->h;

    signalstats->chromaw = FF_CEIL_RSHIFT(inlink->w, signalstats->hsub);
    signalstats->chromah = FF_CEIL_RSHIFT(inlink->h, signalstats->vsub);

    signalstats->fs = inlink->w * inlink->h;
    signalstats->cfs = signalstats->chromaw * signalstats->chromah;

    if (signalstats->filters & 1<<FILTER_HEADSWITCHING) {
        signalstats->filter_head_border = av_malloc(inlink->h * sizeof(*signalstats->filter_head_border));
        signalstats->filter_head_order  = av_malloc(inlink->h * sizeof(*signalstats->filter_head_order));
        if (!signalstats->filter_head_border || !signalstats->filter_head_order)
            return AVERROR(ENOMEM);
    }

    if (signalstats->filters & 1<<FILTER_VREP) {
        signalstats->vrep_line = av_malloc(inlink->h * sizeof(*signalstats->vrep_line));
        if (!signalstats->vrep_line)
            return AVERROR(ENOMEM);
    }

    return 0;
}

static void burn_frame(AVFrame *f, int x, int y)
{
    f->data[0][y * f->linesize[0] + x] = 235;
}

static void filter_init_head(signalstatsContext *signalstats, const AVFrame *p, int w, int h)
{
    int y;
    int tol = 16; // this needs to be configurable.
    int lw = p->linesize[0];
    unsigned int *order = signalstats->filter_head_order;
    int median;
    //int switching =1;

    for (y = 0; y < h; y++) {
        int yw = y*lw;
        int x;
        // try not to get fooled by non matted video.
        int col = 16;

        signalstats->filter_head_border[y] = 0;

        if (p->data[0][yw] <= 16)
            col = p->data[0][yw];

        for (x = 0; x < w; x++) {
            order[y] = x;
            if (abs(col - p->data[0][yw + x]) > tol) {
                signalstats->filter_head_border[y] = x;
                break;
            }
        }
    }

    // there is probably a better and easier (and faster) sorting method
    for (y = 0; y < h; y++ ) {
        int min = order[y];
        int select = y;
        int x;

        for (x = y; x < h; x++) {
            if (order[x] < min) {
                min = order[x];
                select = x;
            }

            if (select != y) {
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
    for (y = 0; y < h; y++)
        if (signalstats->filter_head_border[y] < median)
            signalstats->filter_head_border[y] = 0;

    /*
    for (y = h - 1; y > 0; y--) {
       if (signalstats->filter_head_border[y] <= median )
           switching = 0;
        if (!switching)
            signalstats->filter_head_border[y] = 0;
    }
    */
}

static int filter_head(signalstatsContext *signalstats, const AVFrame *in, AVFrame *out, int y, int w, int h)
{
    const int filt = signalstats->filter_head_border[y];
    if (out && filt) {
        int x;
        for (x = 0; x < w; x++)
            burn_frame(out, x, y);
    }
    return filt;
}

static void filter_init_range(signalstatsContext *signalstats, const AVFrame *p, int w, int h)
{
}

static int filter_range(signalstatsContext *signalstats, const AVFrame *in, AVFrame *out, int y, int w, int h)
{
    int x, score = 0;
    const int yc = FF_CEIL_RSHIFT(y, signalstats->vsub);
    const uint8_t *pluma    = &in->data[0][y  * in->linesize[0]];
    const uint8_t *pchromau = &in->data[1][yc * in->linesize[1]];
    const uint8_t *pchromav = &in->data[2][yc * in->linesize[2]];

    for (x = 0; x < w; x++) {
        const int xc = FF_CEIL_RSHIFT(x, signalstats->hsub);
        const int luma    = pluma[x];
        const int chromau = pchromau[xc];
        const int chromav = pchromav[xc];
        const int filt = luma    < 16 || luma    > 235 ||
                         chromau < 16 || chromau > 240 ||
                         chromav < 16 || chromav > 240;
        score += filt;
        if (out && filt)
            burn_frame(out, x, y);
    }
    return score;
}

static void filter_init_tout(signalstatsContext *signalstats, const AVFrame *p, int w, int h)
{
}

static int filter_tout_outlier(uint8_t x, uint8_t y, uint8_t z)
{
    return ((abs(x - y) + abs (z - y)) / 2) - abs(z - x) > 4; // make 4 configurable?
}

static int filter_tout(signalstatsContext *signalstats, const AVFrame *in, AVFrame *out, int y, int w, int h)
{
    const uint8_t *p = in->data[0];
    int lw = in->linesize[0];
    int x, score = 0, filt;

    if (y - 1 < 0 || y + 1 >= h)
        return 0;

    // detect two pixels above and below (to eliminate interlace artefacts)
    // should check that video format is infact interlace.

#define FILTER(i, j) \
filter_tout_outlier(p[(y-j) * lw + x + i], \
                    p[    y * lw + x + i], \
                    p[(y+j) * lw + x + i])

#define FILTER3(j) (FILTER(-1, j) && FILTER(0, j) && FILTER(1, j))

    if (y-2 >= 0 && y+2 < h) {
        for (x = 1; x < w - 1; x++) {
            filt = FILTER3(2) && FILTER3(1);
            score += filt;
            if (filt && out)
                burn_frame(out, x, y);
        }
    } else {
        for (x = 1; x < w - 1; x++) {
            filt = FILTER3(1);
            score += filt;
            if (filt && out)
                burn_frame(out, x, y);
        }
    }
    return score;
}

#define VREP_START 4

static void filter_init_vrep(signalstatsContext *signalstats, const AVFrame *p, int w, int h)
{
    int i, y;
    int lw = p->linesize[0];

    for (y = VREP_START; y < h; y++) {
        int totdiff = 0;
        int y2lw = (y - VREP_START) * lw;
        int ylw = y * lw;

        for (i = 0; i < w; i++)
            totdiff += abs(p->data[0][y2lw + i] - p->data[0][ylw + i]);

        /* this value should be definable */
        signalstats->vrep_line[y] = totdiff < w;
    }
}

static int filter_vrep(signalstatsContext *signalstats, const AVFrame *in, AVFrame *out, int y, int w, int h)
{
    int x, score = 0;

    if (y < VREP_START)
        return 0;

    for (x = 0; x < w; x++) {
        if (signalstats->vrep_line[y]) {
            score++;
            if (out)
                burn_frame(out, x, y);
        }
    }
    return score;
}

static const char *const filter_metanames[] = { "TOUT", "VREP", "RANG", "HEAD", NULL };

static int (*filter_call[FILT_NUMB])(signalstatsContext *signalstats, const AVFrame *in, AVFrame *out, int y, int w, int h) = {
    filter_tout,
    filter_vrep,
    filter_range,
    filter_head,
};

static void (*filter_init[FILT_NUMB])(signalstatsContext *signalstats, const AVFrame *p,  int w, int h) = {
    filter_init_tout,
    filter_init_vrep,
    filter_init_range,
    filter_init_head,
};

#define DEPTH 256

static int filter_frame(AVFilterLink *link, AVFrame *in)
{
    signalstatsContext *signalstats = link->dst->priv;
    AVFilterLink *outlink = link->dst->outputs[0];
    AVFrame *out = in;
    int i, j;
    int  w = 0,  cw = 0, // in
        pw = 0, cpw = 0; // prev
    int yuv,yuvu,yuvv;
    int fil;
    char metabuf[128];
    unsigned int histy[DEPTH] = {0},
                 histu[DEPTH] = {0},
                 histv[DEPTH] = {0},
    histhue[360] = {0},
                histsat[DEPTH] = {0}; // limited to 8 bit data.
    int miny = -1, minu = -1, minv = -1;
    int maxy = -1, maxu = -1, maxv = -1;
    int lowy  = -1, lowu  = -1, lowv  = -1;
    int highy = -1, highu = -1, highv = -1;
    int minsat=-1,maxsat=-1,lowsat=-1,highsat=-1;
    int lowp, highp, clowp, chighp;
    int accy, accu, accv;
    int accsat,acchue=0;
    int medhue,modhue,maxhue;
    int toty = 0, totu = 0, totv = 0,totsat=0;
    int tothue = 0;
    int dify = 0, difu = 0, difv = 0;
    int dify1 = 0, dify2 = 0;

    int filtot[FILT_NUMB] = {0};
    AVFrame *prev;

    if (!signalstats->frame_prev)
        signalstats->frame_prev = av_frame_clone(in);
    
    prev = signalstats->frame_prev;

    if (signalstats->outfilter != FILTER_NONE)
        out = av_frame_clone(in);

    for (fil = 0; fil < FILT_NUMB; fil ++) {
        if (signalstats->filters & 1<<fil)
            filter_init[fil](signalstats, in, link->w, link->h);
    }

    // Calculate luma histogram and difference with previous frame or field.
    for (j = 0; j < link->h; j++) {
        for (i = 0; i < link->w; i++) {

            yuv = in->data[0][w + i];
            histy[yuv]++;

            dify += abs(in->data[0][w + i] - prev->data[0][pw + i]);

            //if (in->interlaced_frame && (j % 2 == 0)) // every second line
            if (!(j & 1)) { // every second line

               // get bottom field
               // should check that we are not currently at the bottom line.
                // but who has heard of an interlaced file with odd vertical dimentions?
                int yuvi = in->data[0][w + in->linesize[0] + i];

                // dif2 = diff bottom field with top field

                dify2 += abs(yuv - yuvi);

                if (in->top_field_first)
                {
                    // dif1 = diff top field with prev bottom field
                    dify1 += abs(yuv - prev->data[0][pw + prev->linesize[0] + i]);
                } else {
                    // dif1 = diff bottom field with prev top field
                    dify1 += abs(yuvi - prev->data[0][pw + i]);
                }

            }
        }
        w  += in->linesize[0];
        pw  += prev->linesize[0];
    }

    // Calculate chroma histogram and difference with previous frame or field.
    for (j = 0; j < signalstats->chromah; j++) {
        for (i = 0; i < signalstats->chromaw; i++) {
            yuvu = in->data[1][cw+i];
            yuvv = in->data[2][cw+i];
            histu[yuvu]++;
            difu += abs(in->data[1][cw+i] - prev->data[1][cpw+i]);
            histv[yuvv]++;
            difv += abs(in->data[2][cw+i] - prev->data[2][cpw+i]);
            
            // int or round?
           int sat = ff_sqrt((yuvu-128) * (yuvu-128) + (yuvv-128)* (yuvv-128));
            histsat[sat]++;
            int hue = floor( (180 / M_PI ) *  atan2f (yuvu-128, yuvv-128) + 180);
            histhue[hue] ++;
        }
        cw += in->linesize[1];
        cpw += prev->linesize[1];
    }

    for (j = 0; j < link->h; j++) {
        for (fil = 0; fil < FILT_NUMB; fil ++) {
            if (signalstats->filters & 1<<fil) {
                AVFrame *dbg = out != in && signalstats->outfilter == fil ? out : NULL;
                filtot[fil] += filter_call[fil](signalstats, in, dbg, j, link->w, link->h);
            }
        }
    }

    // find low / high based on histogram percentile
    // these only need to be calculated once.

    lowp   = signalstats->fs  * 10 / 100;
    highp  = signalstats->fs  * 95 / 100;
    clowp  = signalstats->cfs * 10 / 100;
    chighp = signalstats->cfs * 95 / 100;

    accy = 0; accu=0; accv=0; accsat =0;
    for (fil = 0; fil < DEPTH; fil++) {
        if (miny < 0 && histy[fil]) miny = fil;
        if (minu < 0 && histu[fil]) minu = fil;
        if (minv < 0 && histv[fil]) minv = fil;
        if (minsat <0 && histsat[fil]) minsat = fil;

        if (histy[fil]) maxy = fil;
        if (histu[fil]) maxu = fil;
        if (histv[fil]) maxv = fil;
        if (histsat[fil]) maxsat = fil;

        
        toty += histy[fil] * fil;
        totu += histu[fil] * fil;
        totv += histv[fil] * fil;
        totsat += histsat[fil] * fil;

        accy += histy[fil];
        accu += histu[fil];
        accv += histv[fil];
        accsat += histsat[fil];

        if (lowy == -1 && accy >=  lowp) lowy = fil;
        if (lowu == -1 && accu >= clowp) lowu = fil;
        if (lowv == -1 && accv >= clowp) lowv = fil;
        if (lowsat == -1 && accsat >= clowp) lowsat = fil;

        
        if (highy == -1 && accy >=  highp) highy = fil;
        if (highu == -1 && accu >= chighp) highu = fil;
        if (highv == -1 && accv >= chighp) highv = fil;
        if (highsat == -1 && accsat >= chighp) highsat = fil;

    }
    
    maxhue = histhue[0]; modhue = 0;
    medhue = -1;
    for (fil = 0; fil < 360; fil++) {
        tothue += histhue[fil] * fil;
        acchue += histhue[fil];
        
        if (medhue == -1 && acchue > signalstats->cfs / 2) medhue = fil;
        if (histhue[fil] > maxhue ) {
            maxhue = histhue[fil];
            modhue = fil;
        }
    }

    av_frame_free(&signalstats->frame_prev);
    signalstats->frame_prev = av_frame_clone(in);

#define SET_META(key, fmt, val) do {                                \
    snprintf(metabuf, sizeof(metabuf), fmt, val);                   \
    av_dict_set(&out->metadata, "lavfi.values." key, metabuf, 0);   \
} while (0)

    SET_META("YMIN",  "%d", miny);
    SET_META("YLOW",  "%d", lowy);
    SET_META("YAVG",  "%g", 1.0 * toty / signalstats->fs);
    SET_META("YHIGH", "%d", highy);
    SET_META("YMAX",  "%d", maxy);
    
    SET_META("UMIN",  "%d", minu);
    SET_META("ULOW",  "%d", lowu);
    SET_META("UAVG",  "%g", 1.0 * totu / signalstats->cfs);
    SET_META("UHIGH", "%d", highu);
    SET_META("UMAX",  "%d", maxu);
    
    SET_META("VMIN",  "%d", minv);
    SET_META("VLOW",  "%d", lowv);
    SET_META("VAVG",  "%g", 1.0 * totv / signalstats->cfs);
    SET_META("VHIGH", "%d", highv);
    SET_META("VMAX",  "%d", maxv);
    
    SET_META("SATMIN", "%d", minsat);
    SET_META("SATLOW", "%d", lowsat);
    SET_META("SATAVG", "%g", 1.0 * totsat / signalstats->cfs);
    SET_META("SATHIGH", "%d", highsat);
    SET_META("SATMAX", "%d", maxsat);
    
    SET_META("HUEMOD","%d",modhue);
    SET_META("HUEMED","%d",medhue);
    SET_META("HUEAVG","%g",1.0 * tothue / signalstats->cfs);

    SET_META("YDIF",  "%g", 1.0 * dify / signalstats->fs);
    SET_META("UDIF",  "%g", 1.0 * difu / signalstats->cfs);
    SET_META("VDIF",  "%g", 1.0 * difv / signalstats->cfs);
    SET_META("YDIF1", "%g", 1.0 * dify1 / (signalstats->fs/2));
    SET_META("YDIF2", "%g", 1.0 * dify2 / (signalstats->fs/2));

    for (fil = 0; fil < FILT_NUMB; fil ++) {
        if (signalstats->filters & 1<<fil) {
            char metaname[128];
            snprintf(metabuf,  sizeof(metabuf),  "%g", 1.0 * filtot[fil] / signalstats->fs);
            snprintf(metaname, sizeof(metaname), "lavfi.values.%s", filter_metanames[fil]);
            av_dict_set(&out->metadata, metaname, metabuf, 0);
        }
    }

    if (in != out)
        av_frame_free(&in);
    return ff_filter_frame(outlink, out);
}

static const AVFilterPad signalstats_inputs[] = {
    {
        .name           = "default",
        .type           = AVMEDIA_TYPE_VIDEO,
      //  .min_perms      = AV_PERM_READ,
        .filter_frame   = filter_frame,
    },
    { NULL }
};

static const AVFilterPad signalstats_outputs[] = {
    {
        .name           = "default",
        .config_props   = config_props,
        .type           = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};

AVFilter ff_vf_signalstats = {
    .name          = "signalstats",
    .description   = "Extract various metrics.",
    .init          = init,
    .uninit        = uninit,
    .query_formats = query_formats,
    .priv_size     = sizeof(signalstatsContext),
    .inputs        = signalstats_inputs,
    .outputs       = signalstats_outputs,
    .priv_class    = &signalstats_class,
};
