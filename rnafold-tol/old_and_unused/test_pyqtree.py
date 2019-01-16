from pyqtree import Index

spindex = Index(bbox=(0,0,100,100))

for item in items:
    speindex.insert( item, item.bbox )
    
