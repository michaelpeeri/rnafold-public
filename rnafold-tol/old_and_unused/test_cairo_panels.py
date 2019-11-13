import json
#import argparse
import sys
import cairo

#from math import pi

#(width, height) = (1280, 768)
#fileOut = 'test_cairo_panels.py.out.png'
#fileImage1 = '/home/michael/OneDrive/Figures/20171206/tree_phenotypes_regression.out.by.taxgroups.10.Terrabacteria.pdf.png'
#fileImage2 = '/home/michael/OneDrive/Figures/20171206/tree_phenotypes_regression.out.by.taxgroups.66.all.dLFE.pdf.png'
#fileImage3 = '/home/michael/OneDrive/Figures/20171206/tree_phenotypes_regression.out.by.taxgroups.10.all.shuffled_only.pdf.png'


#config_file = "test_cairo_panels_test01.json"
config_file = sys.argv[1]
annotationsFont = "Arial"


def addPanel(ctx, label, filename, x, y, w, h, config, panelConfig):
    ctx.save()
    #ctx.new_path()

    im1 = cairo.ImageSurface.create_from_png( file( filename, 'r') )

    # Set origin for this layer
    ctx.translate( x, y )

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



    #ctx.set_source_surface(im1, 0, 0)  # draw image at specific position (in context-relative units)
    ctx.set_source(imgpat)
    #ctx.move_to(245, 245)


    ctx.rectangle( 0, 0, w, h )
    #ctx.scale(0.2, 0.2)

    ctx.fill()
    #ctx.paint()
    #ctx.save()
    ctx.restore()


    ctx.save()
    ctx.translate( x, y-0.02 )
    ctx.set_source_rgb( 0.0, 0.0, 0.0 )
    ctx.select_font_face( annotationsFont )
    ctx.set_font_size( 50 )
    ctx.show_text( label )
    ctx.restore()


def addPanelUsingConfig( ctx, config, panelConfig ):
    outputWidth  = config["plotOutput"]["width"]
    outputHeight = config["plotOutput"]["height"]
    
    #print(panelConfig)
    return addPanel(
        ctx,
        panelConfig["label"],
        panelConfig["fileName"],
        panelConfig["x"] * outputWidth,
        panelConfig["y"] * outputHeight,
        panelConfig["w"] * outputWidth,
        panelConfig["h"] * outputHeight,
        config,
        panelConfig
                   )


def addAnnotation( ctx, x, y, size, config, annotationConfig ):
    ctx.save()
    ctx.translate( x, y )
    ctx.set_source_rgb( 0.0, 0.0, 0.0 )
    ctx.select_font_face( annotationsFont ) 
    ctx.set_font_size( size )
    ctx.show_text( annotationConfig["text"] )
    ctx.restore()


def addAnnotationUsingConfig( ctx, config, annotationConfig ):
    outputWidth  = config["plotOutput"]["width"]
    outputHeight = config["plotOutput"]["height"]

    addAnnotation( ctx,
                   annotationConfig["x"] * outputWidth ,
                   annotationConfig["y"] * outputHeight,
                   annotationConfig["size"],
                   config,
                   annotationConfig
    )
    
def createPlot( config ):
    fo = file(config["plotOutput"]["fileName"], 'w')


    surface = cairo.ImageSurface(cairo.Format.RGB24, config["plotOutput"]["width"], config["plotOutput"]["height"] )

    ctx = cairo.Context( surface )

    # ------------------------- Background -------------------------
    ctx.set_source_rgb( 1.0, 1.0, 1.0 )
    ctx.paint()


    #ctx.set_source_rgb( 0.3, 0.2, 0.5 )
    #ctx.arc( width/2, height/2, height/2, 0, 2*pi )
    #ctx.stroke()


    for panelConfig in config["panels"]:
        addPanelUsingConfig( ctx, config, panelConfig )

    if "annotations" in config:
        for annotationConfig in config["annotations"]:
            addAnnotationUsingConfig( ctx, config, annotationConfig )



    # ------------------------- Image 1 -------------------------

    #addPanel( "A", fileImage2, width*0.01, height*0.10, width*0.8, height*0.9, config )

    # ------------------------- Image 2 -------------------------

    #addPanel( "B", fileImage1, width*0.72, height*0.08, width*0.3, height*0.27, config )

    #addPanel( "C", fileImage1, width*0.72, height*0.45, width*0.5, height*0.27, config )


    #addPanel( "A", fileImage2, width*0.05, height*0.05, width*0.29, height*0.9, config )
    #addPanel( "B", fileImage2, width*(0.05 + 0.33), height*0.05, width*0.29, height*0.9, config )
    #addPanel( "C", fileImage2, width*(0.05 + 0.66), height*0.05, width*0.29, height*0.9, config )

    


    # ------------------------- Finalize -------------------------

    #ctx.set_source_rgb( 0.3, 0.2, 0.5 )
    #ctx.arc( 30, 30, 30, 0, 2*pi )
    #ctx.stroke()


    surface.write_to_png(fo)

    surface.finish()
    



config = None
with open(config_file, "r") as config_file:
    config = json.loads(config_file.read())
#print(config)
createPlot(config)



