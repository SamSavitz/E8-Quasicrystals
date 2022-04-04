/*
 * Based on Fabrice Bellard's encode_video.c:
 * video encoding with libavcodec API example
 * Fabrice Bellard's encode_video.c
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>

extern "C" {
    #include <libavcodec/avcodec.h>
    #include <libavutil/opt.h>
    #include <libavutil/imgutils.h>
}

using namespace std;

struct RGB {
    unsigned char r, g, b;

    RGB(  ) {  }

    RGB(unsigned char r, unsigned char g, unsigned char b) :
	r(r),
	g(g),
	b(b)
	{  }
};

struct YUV {
    unsigned char y, u, v;

    YUV(const RGB c) {
        int yt =  77*c.r + 150*c.g +  29*c.b;
        int ut = -43*c.r -  84*c.g + 127*c.b;
        int vt = 127*c.r - 106*c.g -  21*c.b;

        y =  (yt + 128) >> 8       ;
        u = ((ut + 128) >> 8) + 128;
        v = ((vt + 128) >> 8) + 128;
    }
};

struct YUV420 {
    unsigned char y[4];
    unsigned char u, v;

    YUV420() = default;

    YUV420(const YUV a, const YUV b, const YUV c, const YUV d) {
        y[0] = a.y;
        y[1] = b.y;
        y[2] = c.y;
        y[3] = d.y;

        u = ((int) a.u + (int) b.u + (int) c.u + (int) d.u + 2) >> 2;
        v = ((int) a.v + (int) b.v + (int) c.v + (int) d.v + 2) >> 2;
    }
};

void rgb2yuv420(const int W, const int H, const RGB* const rgb, YUV420* const yuv) {
    for (unsigned int y = 0; y < H/2; ++y)
        for (unsigned int x = 0; x < W/2; ++x)
            yuv[W/2*y + x] = YUV420(
                YUV(rgb[W*y + x        ] ),
                YUV(rgb[W*y + x +     1] ),
                YUV(rgb[W*y + x + W    ] ),
                YUV(rgb[W*y + x + W + 1] )
            );
}

class RGBEncoder {
public:
    RGBEncoder(const char* filename, const int W, const int H, const int framerate, const double bpp, const char* codec_name="libx264rgb") :
        codec(avcodec_find_encoder_by_name(codec_name)) {

        if (!codec)
            throw "Codec not found";

        c = avcodec_alloc_context3(codec);

        if (!c)
            throw "Unable to allocate context";

        pkt = av_packet_alloc();

        if (!pkt)
            throw "Unable to allocate packet";

        c->bit_rate = round(bpp*W*H*framerate);
        c->width = W;
        c->height = H;
        c->time_base = (AVRational) {1, framerate};
        c->framerate = (AVRational) {framerate, 1};
        c->gop_size = 10;
        c->max_b_frames = 1;
        c->pix_fmt = AV_PIX_FMT_RGB24;

        if (codec->id == AV_CODEC_ID_H264)
            av_opt_set(c->priv_data, "preset", "slow", 0);

        if (avcodec_open2(c, codec, NULL) < 0)
            throw "Codec error";

        f.open(filename);

        frame = av_frame_alloc();

        if (!frame)
            throw "Unable to allocate frame";

        frame->format = c->pix_fmt;
        frame->width  = W;
        frame->height = H;

        if (av_frame_get_buffer(frame, 0) < 0)
            throw "Frame error";
    }

    void writeFrame(const RGB* const data) {
        if (av_frame_make_writable(frame) < 0)
            throw "Error making frame writable";

        for (int y = 0; y < c->height; ++y)
            for (int x = 0; x < c->width; ++x) {
                const RGB& p = data[frame->width*y + x];

                frame->data[0] [frame->linesize[0]*y + 3*x    ] = p.r;
                frame->data[0] [frame->linesize[0]*y + 3*x + 1] = p.g;
                frame->data[0] [frame->linesize[0]*y + 3*x + 2] = p.b;
            }

        frame->pts = frameNumber;

        /* encode the image */
        encode();

        frameNumber++;
    }

    void close() {
        static const char endcode[] = "\0\0\x01\xb7";

        frame = nullptr;
        encode();

        f.flush();

        if (codec->id == AV_CODEC_ID_MPEG1VIDEO || codec->id == AV_CODEC_ID_MPEG2VIDEO)
            f.write(endcode, sizeof(endcode));

        f.close();

        avcodec_free_context(&c);
        av_frame_free(&frame);
        av_packet_free(&pkt);
    }

    ~RGBEncoder() {
        close();
    }

protected:
    const AVCodec* codec;
    AVCodecContext* c = nullptr;
    ofstream f;
    AVFrame* frame;
    AVPacket* pkt;
    int frameNumber = 0;

    void encode() {
        auto ret = avcodec_send_frame(c, frame);

        if (ret < 0)
            throw "Error sending a frame for encoding";

        while (ret >= 0) {
            ret = avcodec_receive_packet(c, pkt);

            if (ret == AVERROR(EAGAIN) || ret == AVERROR_EOF)
                return;
            else if (ret < 0)
                throw "Error during encoding";
            else {
                //printf("Write packet %3d (size=%5d)\n", pkt->pts, pkt->size);
                f.write((const char*) pkt->data, pkt->size);
                av_packet_unref(pkt);
            }
        }
    }
};

class YUV420Encoder {
public:
    YUV420Encoder(const char* filename, const int W, const int H, const int framerate, const double bpp, const char* codec_name="libx264") :
        codec(avcodec_find_encoder_by_name(codec_name)) {

        if (!codec)
            throw "Codec not found";

        c = avcodec_alloc_context3(codec);

        if (!c)
            throw "Unable to allocate context";

        pkt = av_packet_alloc();

        if (!pkt)
            throw "Unable to allocate packet";

        c->bit_rate = round(bpp*W*H*framerate);
        c->width = W;
        c->height = H;
        c->time_base = (AVRational) {1, framerate};
        c->framerate = (AVRational) {framerate, 1};
        c->gop_size = 10;
        c->max_b_frames = 1;
        c->pix_fmt = AV_PIX_FMT_YUV420P;

        if (codec->id == AV_CODEC_ID_H264)
            av_opt_set(c->priv_data, "preset", "slow", 0);

        if (avcodec_open2(c, codec, NULL) < 0)
            throw "Codec error";

        f.open(filename);

        frame = av_frame_alloc();

        if (!frame)
            throw "Unable to allocate frame";

        frame->format = c->pix_fmt;
        frame->width  = W;
        frame->height = H;

        if (av_frame_get_buffer(frame, 0) < 0)
            throw "Frame error";
    }

    void writeFrame(const YUV420* const data) {
        if (av_frame_make_writable(frame) < 0)
            throw "Error making frame writable";

        for (int y = 0; y < c->height/2; ++y)
            for (int x = 0; x < c->width/2; ++x) {
                const auto& p = data[frame->width*y + x];
                const auto LS = frame->linesize[0];
                const auto i = 2*(LS*y + x);

                frame->data[0] [i         ] = p.y[0];
                frame->data[0] [i +      1] = p.y[1];
                frame->data[0] [i + LS    ] = p.y[2];
                frame->data[0] [i + LS + 1] = p.y[3];

                frame->data[1] [frame->linesize[1]*y + x] = p.u;
                frame->data[2] [frame->linesize[2]*y + x] = p.v;
            }

        frame->pts = frameNumber;

        /* encode the image */
        encode();

        frameNumber++;
    }

    void close() {
        static const char endcode[] = "\0\0\x01\xb7";

        frame = nullptr;
        encode();

        f.flush();

        if (codec->id == AV_CODEC_ID_MPEG1VIDEO || codec->id == AV_CODEC_ID_MPEG2VIDEO)
            f.write(endcode, sizeof(endcode));

        f.close();

        avcodec_free_context(&c);
        av_frame_free(&frame);
        av_packet_free(&pkt);
    }

    ~YUV420Encoder() {
        close();
    }

protected:
    const AVCodec* codec;
    AVCodecContext* c = nullptr;
    ofstream f;
    AVFrame* frame;
    AVPacket* pkt;
    int frameNumber = 0;

    void encode() {
        auto ret = avcodec_send_frame(c, frame);

        if (ret < 0)
            throw "Error sending a frame for encoding";

        while (ret >= 0) {
            ret = avcodec_receive_packet(c, pkt);

            if (ret == AVERROR(EAGAIN) || ret == AVERROR_EOF)
                return;
            else if (ret < 0)
                throw "Error during encoding";
            else {
                f.write((const char*) pkt->data, pkt->size);
                av_packet_unref(pkt);
            }
        }
    }
};
