def addMargins(x0, x1, margin=0.1):
    if( x0 > x1 ):
        (x0,x1) = (x1,x0)
    assert(x1 >= x0)

    _margin = abs(x1-x0) * 0.1
    y0, y1 = 0, 0
    y0 = x0 - _margin

    y1 = x1 + _margin

    return (y0,y1)
