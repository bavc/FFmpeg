/*
 ** <h3>values graph</h3>
 * copyright (c) 2010 Mark Heath mjpeg0 @ silicontrip dot net 
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

#include "libavcodec/avcodec.h"
#include "libavformat/avformat.h"
#include "libavutil/pixfmt.h"
#include "libavutil/pixdesc.h"
#include "libswscale/swscale.h"
#include "avfilter.h"
#include "formats.h"
#include "internal.h"
#include "video.h"

#include "libavutil/opt.h"


#include <strings.h>
#include <fcntl.h>

/* time
 10 April 06:58 - 07:38
 15 April 16:15 - 16:50
 16 April 07:10 - 7:40
 17 April 16:20 - 17:00
 18 April 7:10 - 8:10, 12:50 - 13:05
 */


/* Prototypes for filter functions */

static int filter_tout(AVFrame *p, int x, int y, int w, int h);

enum FilterMode {
    FILTER_NONE = -1,
    FILTER_TOUT,
    FILT_NUMB
};

static int (*filter_call[FILT_NUMB])(AVFrame *p, int x, int y, int w, int h) = {
    filter_tout
};

static const char *const filter_metanames[] = { "TOUT", NULL };
static const char *const filter_names[] = { "tout", NULL };




typedef struct
{
    const AVClass *class;

    
	FILE  *fh;
    char *filename;
    
    int chromah;
    int chromaw;
    int fc;
    
    int fs;
    int cfs;
    
    enum FilterMode outfilter;
    int filter[FILT_NUMB];
    char *statistics_str;
    
    AVFrame *frame_prev;
    
} valuesContext;


#define OFFSET(x) offsetof(valuesContext, x)
#define FLAGS AV_OPT_FLAG_FILTERING_PARAM|AV_OPT_FLAG_VIDEO_PARAM

int filter_tout_outlier(uint8_t x, uint8_t y, uint8_t z);

static const AVOption values_options[]= {
    {"filename", "set output file", OFFSET(filename), AV_OPT_TYPE_STRING, {.str=NULL},  CHAR_MIN, CHAR_MAX},
    {"f", "set output file", OFFSET(filename), AV_OPT_TYPE_STRING, {.str=NULL},  CHAR_MIN, CHAR_MAX},
    {"out", "set video filter", OFFSET(outfilter), AV_OPT_TYPE_INT, {.i64=FILTER_NONE}, -1, FILT_NUMB-1,FLAGS,"out"},
    {"tout", "", 0, AV_OPT_TYPE_CONST, {.i64=FILTER_TOUT}, 0,0,FLAGS,"out"},
    {"stat","set | seperated statistics filter", OFFSET(statistics_str),AV_OPT_TYPE_STRING, {.str=NULL},  CHAR_MIN, CHAR_MAX},
    {NULL}
};

AVFILTER_DEFINE_CLASS(values);


// static av_cold int init(AVFilterContext *ctx, const char *args, void *opaque)
static av_cold int init(AVFilterContext *ctx, const char *args)
{
    valuesContext *values = ctx->priv;
    int ret;

    av_log(ctx, AV_LOG_DEBUG, ">>> init().\n");

    values->class = &values_class;

    av_opt_set_defaults(values);
    
    av_log(ctx, AV_LOG_DEBUG, "    init() av_set_options_string.\n");

    
    if ((ret = av_set_options_string(values, args, "=", ":")) < 0)
        return ret;

    // parse statistics filter string

    do {
        char *next,*cur=values->statistics_str;
        int fil,ok;
        while (cur) {
            next = strchr(cur,'|');
            if (next)
                *next++=0;
            
            ok=0;
            for (fil = 0; fil < FILT_NUMB; fil ++) {
                
                
                if (strcmp(filter_names[fil],cur)==0)
                {
                    av_log(ctx,AV_LOG_DEBUG, "Found filter: %s\n",filter_names[fil]);

                    ok = 1;
                    values->filter[fil] = 1;
                }
            }
            if (!ok)
            {
                av_log(ctx, AV_LOG_ERROR, "Error parsing: %s.\n", cur);
                av_opt_free(values);
                return AVERROR(EINVAL);
            }
            cur = next;

        }
    } while(0);
    //
    
    values->fc = 0;
    values->fh = NULL;
    
    if (values->filename != NULL)
        values->fh = fopen (values->filename,"w");
    
    
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

    int hsub, vsub;
    
    av_log(ctx, AV_LOG_DEBUG, ">>> config_props().\n");

    
	avcodec_get_chroma_sub_sample(outlink->format, &hsub, &vsub);
	
	outlink->w = inlink->w;
    outlink->h = inlink->h;
	
    values->chromaw = inlink->w >> hsub;
    values->chromah = inlink->h >> vsub;

    values->fs = inlink->w * inlink->h;
    values->cfs = values->chromaw * values->chromah;
    
    
    values->frame_prev = av_frame_alloc();
    values->frame_prev->width  = inlink->w;
    values->frame_prev->height = inlink->h;
    values->frame_prev->format = inlink->format;
    
    if (av_frame_get_buffer(values->frame_prev, 32) < 0)
    {
        av_frame_free(&values->frame_prev);
        return AVERROR(EINVAL);
    }
    
    memset(values->frame_prev->data[0],16,values->frame_prev->linesize[0]*values->frame_prev->height);
    memset(values->frame_prev->data[1],128,values->frame_prev->linesize[1]*values->frame_prev->height >> vsub);
    memset(values->frame_prev->data[2],128,values->frame_prev->linesize[2]*values->frame_prev->height >> vsub);
    
    av_log(ctx, AV_LOG_DEBUG, "<<< config_props().\n");

	return 0;
	
}


int filter_tout_outlier(uint8_t x, uint8_t y, uint8_t z)
{
	
	int dif;
	
	dif =  ((abs(x - y) + abs (z - y) ) / 2) - abs(z-x);
	
	//fprintf(stderr,"dif: %d %d\n",dif2,dif1);
	
	// Will make this configurable by command line option.
	return dif>4?1:0;
	
}


static int filter_tout(AVFrame *p, int x, int y, int w, int h) {
	
	int lw = p->linesize[0];
    
	if ((x-1 < 0) || (x+1 > w) || (y-1 < 0) || (y+1 >= h)) {
		return 0;
	} else {
		int i;
        
        
		for (i=-1; i<2; i++)
		{
			// detect two pixels above and below (to eliminate interlace artefacts)
			if ((y-2 >=0) && (y+2 < h)) {
				if (!filter_tout_outlier(p->data[0][(y-2) * lw + i+x], p->data[0][y * lw + i+x], p->data[0][(y+2) * lw + i+x]))
					return 0;
			}
			
			if (!filter_tout_outlier(p->data[0][(y-1) * lw + i+x], p->data[0][y * lw + i+x], p->data[0][(y+1) * lw + i+x]))
				return 0;
			
		}
	}
	return 1;
	
}

static int filter_frame(AVFilterLink *link, AVFrame *in)
{
   
    valuesContext *values = link->dst->priv;
    AVFilterLink *outlink = link->dst->outputs[0];
    AVFrame *out; // = link->dst->outputs[0]->outpic;
    
	int i,j;
	int cw =0 ,w=0,ow=0,cow=0;
    int yuv;
    int fil;
    char metabuf[128];

    
	int miny,minu,minv;
	int toty=0,totu=0,totv=0;
    int dify=0,difu=0,difv=0;
	int maxy,maxu,maxv;

    int filtot[FILT_NUMB];
    
    for (i=0; i<FILT_NUMB; i++)
        filtot[i]=0;
    
    out = ff_get_video_buffer(outlink, outlink->w, outlink->h);
    av_frame_copy_props(out, in);

    miny = in->data[0][0];
	maxy = in->data[0][0];
	
	minu = in->data[1][0];
	maxu = in->data[1][0];
	
	minv = in->data[2][0];
	maxv = in->data[2][0];

    for (j=0; j<link->h; j++) {
        for (i=0;i<link->w;i++) {
           
            yuv = in->data[0][w+i];

            if (yuv > maxy) maxy=yuv;
			if (yuv < miny) miny=yuv;
			toty += yuv;

            out->data[0][ow+i] = 16;

            dify  += abs(in->data[0][w+i] - values->frame_prev->data[0][w+i]);

            
            if (i<values->chromaw && j<values->chromah) {
				yuv = in->data[1][cw+i];
				if (yuv > maxu) maxu=yuv;
				if (yuv < minu) minu=yuv;
				totu += yuv;
				
				yuv = in->data[2][cw+i];
				if (yuv > maxv) maxv=yuv;
				if (yuv < minv) minv=yuv;
				totv += yuv;

                difu  += abs(in->data[1][cw+i] - values->frame_prev->data[1][cw+i]);
                difv  += abs(in->data[2][cw+i] - values->frame_prev->data[2][cw+i]);

                
                out->data[1][cow+i] = 128;
                out->data[2][cow+i] = 128;
                
			}
            
            // magic filter
            // options to disable and enable them
            // option to pick one for video out
            
            for (fil = 0; fil < FILT_NUMB; fil ++) {
                if (filter_call[fil](in,i,j,link->w,link->h)) {
                    if (values->filter[fil] || values->outfilter == fil)
                    {
                        values->filter[fil]=1;
                        filtot[fil] ++;
                        if (values->outfilter == fil)
                        {
                            out->data[0][ow+i] = 235;
                        }
                    }
                }
            }            
        }
        
        cow += out->linesize[1];
        ow += out->linesize[0];
        w += in->linesize[0];
        cw += in->linesize[1];
	}
	
    for (j=0;j<link->h;j++) {
        memcpy(values->frame_prev->data[0] + j * values->frame_prev->linesize[0], in->data[0]+ j * in->linesize[0], link->w);
        memcpy(values->frame_prev->data[1] + j * values->frame_prev->linesize[1], in->data[1]+ j * in->linesize[1], values->chromaw);
        memcpy(values->frame_prev->data[2] + j * values->frame_prev->linesize[2], in->data[2]+ j * in->linesize[2], values->chromaw);

    }
    
    snprintf(metabuf,sizeof(metabuf),"%d",miny);
    av_dict_set(&out->metadata,"lavfi.values.YMIN",metabuf,0);

    snprintf(metabuf,sizeof(metabuf),"%g",1.0 * toty / values->fs);
    av_dict_set(&out->metadata,"lavfi.values.YAVG",metabuf,0);

    snprintf(metabuf,sizeof(metabuf),"%d",maxy);
    av_dict_set(&out->metadata,"lavfi.values.YMAX",metabuf,0);

    snprintf(metabuf,sizeof(metabuf),"%d",minu);
    av_dict_set(&out->metadata,"lavfi.values.UMIN",metabuf,0);
    
    snprintf(metabuf,sizeof(metabuf),"%g",1.0 * totu / values->cfs);
    av_dict_set(&out->metadata,"lavfi.values.UAVG",metabuf,0);
    
    snprintf(metabuf,sizeof(metabuf),"%d",maxu);
    av_dict_set(&out->metadata,"lavfi.values.UMAX",metabuf,0);

    snprintf(metabuf,sizeof(metabuf),"%d",minv);
    av_dict_set(&out->metadata,"lavfi.values.VMIN",metabuf,0);
    
    snprintf(metabuf,sizeof(metabuf),"%g",1.0 * totv / values->cfs);
    av_dict_set(&out->metadata,"lavfi.values.VAVG",metabuf,0);
    
    snprintf(metabuf,sizeof(metabuf),"%d",maxv);
    av_dict_set(&out->metadata,"lavfi.values.VMAX",metabuf,0);

    
    snprintf(metabuf,sizeof(metabuf),"%g",1.0 * difv / values->fs);
    av_dict_set(&out->metadata,"lavfi.values.YDIF",metabuf,0);

    snprintf(metabuf,sizeof(metabuf),"%g",1.0 * difu / values->cfs);
    av_dict_set(&out->metadata,"lavfi.values.UDIF",metabuf,0);

    snprintf(metabuf,sizeof(metabuf),"%g",1.0 * difv / values->cfs);
    av_dict_set(&out->metadata,"lavfi.values.VDIF",metabuf,0);

    
    
    for (fil = 0; fil < FILT_NUMB; fil ++) {
        if (values->filter[fil]) {
            char metaname[128];
            snprintf(metabuf,sizeof(metabuf),"%g",1.0 * filtot[fil]/values->fs);
            snprintf(metaname,sizeof(metaname),"lavfi.values.%s",filter_metanames[fil]);
            av_dict_set(&out->metadata,metaname,metabuf,0);
        }
    }
    
    
    if (values->fh != NULL) {
        fprintf(values->fh,"%d %d %g %d %d %g %d %d %g %d",values->fc,
                miny,1.0 * toty / values->fs, maxy,
                minu,1.0 * totu / values->cfs, maxu,
                minv,1.0 * totv / values->cfs, maxv);
    
        for (fil = 0; fil < FILT_NUMB; fil ++) {
            if (values->filter[fil]) {
                fprintf (values->fh," %g",1.0 * filtot[fil] / values->fs);
            }
        }
    
        fprintf(values->fh,"\n");
    }
    
    values->fc++;
	av_frame_free(&in);
    return ff_filter_frame(outlink, out);
}

static const AVFilterPad avfilter_vf_values_inputs[] = {
    {
        .name            = "default",
        .type            = AVMEDIA_TYPE_VIDEO,
        .min_perms       = AV_PERM_READ,
        .filter_frame = filter_frame,
        
    },
    { NULL }
};

static const AVFilterPad avfilter_vf_values_outputs[] = {
    {
        .name            = "default",
        .config_props     = config_props,
        .type            = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};


AVFilter avfilter_vf_values = {
    .name      = "values",
	.description = ".",
	
	.init = init,
	.uninit    = uninit,
    .query_formats = query_formats,
	
	.priv_size = sizeof(valuesContext),
	
    .inputs    = avfilter_vf_values_inputs,
    .outputs   = avfilter_vf_values_outputs,
    .priv_class = &values_class,

};
