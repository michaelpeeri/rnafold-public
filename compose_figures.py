import json
import sys
import cairo


config_file = sys.argv[1]
annotationsFont = "Arial"



def addAnnotation( ctx, x, y, size, _, annotationConfig ):
    ctx.save()
    ctx.translate( x, y )
    ctx.set_source_rgb( 0.0, 0.0, 0.0 )
    ctx.select_font_face( annotationsFont ) 
    ctx.set_font_size( size )
    ctx.show_text( annotationConfig["text"] )
    ctx.restore()



def addPanel(ctx, label, filename, x, y, w, h, config, panelConfig):
    ctx.save()

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



    ctx.set_source(imgpat)

    ctx.rectangle( 0, 0, w, h )

    ctx.fill()
    ctx.restore()


    # Add the annotation
    # ctx.save()
    # ctx.translate( x, y-0.02 )
    # ctx.set_source_rgb( 0.0, 0.0, 0.0 )
    # ctx.select_font_face( annotationsFont )
    # ctx.set_font_size( 50 )
    # ctx.show_text( label )
    # ctx.restore()
    addAnnotation( ctx, x, y-0.02, 50, config, dict(text=label) )


def addPanelUsingConfig( ctx, config, panelConfig ):
    outputWidth  = config["plotOutput"]["width"]
    outputHeight = config["plotOutput"]["height"]
    
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


    for panelConfig in config["panels"]:
        addPanelUsingConfig( ctx, config, panelConfig )

    # TODO - fix bug causing annotations after the first to appear in the wrong place
    if "annotations" in config:
        for annotationConfig in config["annotations"]:
            addAnnotationUsingConfig( ctx, config, annotationConfig )


    surface.write_to_png(fo)

    surface.finish()
    



config = None
with open(config_file, "r") as config_file:
    config = json.loads(config_file.read())
#print(config)
createPlot(config)



