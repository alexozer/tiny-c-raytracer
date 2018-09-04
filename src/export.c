//
// Created by alex on 9/4/18.
//

#include "export.h"

#include <png.h>
#include <stdlib.h>

void write_png_file(const char *filename, unsigned width, unsigned height, uint8_t *img_buffer) {
    FILE *fp = fopen(filename, "wb");
    if(!fp) abort();

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) abort();

    png_infop info = png_create_info_struct(png);
    if (!info) abort();

    if (setjmp(png_jmpbuf(png))) abort();

    png_init_io(png, fp);

    // Output is 8bit depth, RGBA format.
    png_set_IHDR(
            png,
            info,
            width, height,
            8,
            PNG_COLOR_TYPE_RGBA,
            PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT,
            PNG_FILTER_TYPE_DEFAULT
    );
    png_write_info(png, info);

    // To remove the alpha channel for PNG_COLOR_TYPE_RGB format,
    // Use png_set_filler().
    // png_set_filler(png, 0, PNG_FILLER_AFTER);

    // Construct array of row pointers
    uint8_t **rows = malloc(height * sizeof(uintptr_t));
    for (unsigned i = 0; i < height; i++) {
        rows[i] = &img_buffer[i * width * 4];
    }
    png_write_image(png, rows);
    free(rows);

    png_write_end(png, NULL);

    free(info);
    free(png);
    fclose(fp);
}

