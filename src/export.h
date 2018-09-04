#pragma once

#include <stdint.h>

void write_png_file(const char *filename, unsigned width, unsigned height, uint8_t *img_buffer);
