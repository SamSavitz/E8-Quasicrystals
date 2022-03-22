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
};

class Encoder {
public:
    Encoder(const char* filename, const int W, const int H, const int framerate, const double bpp, const char* codec_name="libx264rgb") :
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
        /* resolution must be a multiple of two */
        c->width = W;
        c->height = H;
        /* frames per second */
        c->time_base = (AVRational) {1, framerate};
        c->framerate = (AVRational) {framerate, 1};

        /* emit one intra frame every ten frames
        * check frame pict_type before passing frame
        * to encoder, if frame->pict_type is AV_PICTURE_TYPE_I
        * then gop_size is ignored and the output of encoder
        * will always be I frame irrespective to gop_size
        */
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
        /* Make sure the frame data is writable.
        On the first round, the frame is fresh from av_frame_get_buffer()
        and therefore we know it is writable.
        But on the next rounds, encode() will have called
        avcodec_send_frame(), and the codec may have kept a reference to
        the frame in its internal structures, that makes the frame
        unwritable.
        av_frame_make_writable() checks that and allocates a new buffer
        for the frame only if necessary.
        */
        if (av_frame_make_writable(frame) < 0)
            throw "Error making frame writable";

        /* Prepare a dummy image.
        In real code, this is where you would have your own logic for
        filling the frame. FFmpeg does not care what you put in the
        frame.
        */
        /* Y */
        /* Cb and Cr */
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

        /* Add sequence end code to have a real MPEG file.
        It makes only sense because this tiny examples writes packets
        directly. This is called "elementary stream" and only works for some
        codecs. To create a valid file, you usually need to write packets
        into a proper file format or protocol; see muxing.c.
        */
        if (codec->id == AV_CODEC_ID_MPEG1VIDEO || codec->id == AV_CODEC_ID_MPEG2VIDEO)
            f.write(endcode, sizeof(endcode));

        f.close();

        avcodec_free_context(&c);
        av_frame_free(&frame);
        av_packet_free(&pkt);
    }

    ~Encoder() {
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
        //if (frame)
        //    printf("Send frame %3d\n", frame->pts);

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
