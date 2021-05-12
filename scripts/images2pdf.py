#!/usr/bin/env python3

import argparse
from PIL import Image

def init_argparser():

    p = argparse.ArgumentParser()
    p.add_argument('-n', type=int, default=5)
    p.add_argument('-o', '--outfile', default='outplot.pdf')
    p.add_argument('infiles', nargs='+')

    return p

def batch(seq, size):
    return [
        seq[i:i + size]
        for i in range(0, len(seq), size)
    ]


def images2pdf( args ):

    pages = []
    counter = 0
    for file_list in batch(args.infiles, args.n):

        imgs = [ Image.open(infile) for infile in file_list ]
        imgs = [ img.resize((img.width//2, img.height//2)) for img in imgs ]

        img_height = sum( [ img.height for img in imgs ] )
        img_width = max( [ img.width for img in imgs ] )

        dest = Image.new('RGB', (img_width, img_height) )

        y_offset = 0
        for img in imgs:
            dest.paste(img, (0, y_offset))
            y_offset += img.height

        pages.append( dest )
        counter += len(imgs)

    pages[0].save(args.outfile, format='PDF', save_all=True, append_images = pages[1:])
    print('Processing depthplot for %d samples' % counter)


def main():
    images2pdf( init_argparser().parse_args() )


if __name__ == '__main__':
    main()






