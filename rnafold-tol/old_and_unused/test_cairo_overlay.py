import cairo

fo = file('test_cairo_overlay.png', 'w')

h=2000
w=2000

surface = cairo.ImageSurface(cairo.Format.RGB24, w, h)
ctx = cairo.Context( surface )


def addPanel(filename):
    ctx.save()

    im1 = cairo.ImageSurface.create_from_png( file( filename, 'r') )

    # Set origin for this layer
    ctx.translate( 0, 0 )

    imgpat = cairo.SurfacePattern( im1 )

    imh = im1.get_height()
    imw = im1.get_width()

    scale_w = imw/w
    scale_h = imh/h
    compromise_scale = max(scale_w, scale_h)
    

    # Scale source image
    scaler = cairo.Matrix()
    scaler.scale(compromise_scale, compromise_scale)
    imgpat.set_matrix(scaler)
    imgpat.set_filter(cairo.FILTER_BEST)



    ctx.set_source(imgpat)

    ctx.rectangle( 0, 0, w, h )

    ctx.fill()
    ctx.restore()
    


# ------------------------- Background -------------------------
ctx.set_source_rgb( 1.0, 1.0, 1.0 )
ctx.paint()

addPanel("/home/michael/rnafold/20180508/abs(dLFE)/log/test2/pca_profiles_density.png")
addPanel("/home/michael/rnafold/20180508/abs(dLFE)/log/test2/pca_profiles.png")

surface.write_to_png(fo)

surface.finish()
